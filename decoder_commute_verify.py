"""
Check whether a learned decoder (decoder_C_{i}.txt) fixes the data error for a given
LUT syndrome pattern, using the same commute constraints as proof_protocol.py:
  fixed_j_x = data_j_x XOR dec_j_x, fixed_j_z = data_j_z XOR dec_j_z
and require fixed to commute with every line in log_txt then stab_txt (symplectic pairs).
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Tuple

from z3 import Bool, BoolVal, is_true, parse_smt2_string, simplify, substitute

DEC_LINE_START = re.compile(r"^(dec\d+_[xz])\t(.*)$")


def decoder_file_index(first_stabilizer_index: int) -> int:
    """decoder_C_{i}.txt with i = 4 - first_stab (first_stabilizer_index in 0..3)."""
    return 4 - first_stabilizer_index


def parse_decoder_c_file(path: Path) -> Tuple[List[str], Dict[str, str]]:
    """Return (meas_var_names in order, dec_name -> sexpr string)."""
    text = path.read_text(encoding="utf-8")
    meas_names: List[str] = []
    formulas: Dict[str, str] = {}
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("meas_var_names:"):
            rest = line.split(":", 1)[1].strip()
            meas_names = [x.strip() for x in rest.split(",") if x.strip()]
            i += 1
            continue
        if line.startswith("#") or not line.strip():
            i += 1
            continue
        m = DEC_LINE_START.match(line)
        if m:
            key = m.group(1).strip()
            parts: List[str] = [m.group(2).strip()]
            i += 1
            while i < len(lines):
                nxt = lines[i]
                if DEC_LINE_START.match(nxt):
                    break
                if nxt.strip():
                    parts.append(nxt.strip())
                i += 1
            formulas[key] = " ".join(parts)
            continue
        i += 1
    if not meas_names:
        raise ValueError(f"No meas_var_names in {path}")
    if not formulas:
        raise ValueError(f"No decoder formulas in {path}")
    return meas_names, formulas


def measurement_assignment_from_bitstring6(
    first_stabilizer_index: int, bitstring6: str
) -> Dict[str, bool]:
    """
    Rounds r_0..r_{t-1}: successful first subround [0,0] -> ancX0=0, flagZ0=0.
    Round r_t: triggering (s,f) = bitstring6[0], bitstring6[1].
    Round r_{t+1}: raw four syndromes on ancX0..3 = bitstring6[2:6].
    """
    t = first_stabilizer_index
    if len(bitstring6) != 6 or any(c not in "01" for c in bitstring6):
        raise ValueError(f"bitstring6 must be 6 bits of 0/1, got {bitstring6!r}")
    bits = [c == "1" for c in bitstring6]
    assign: Dict[str, bool] = {}
    for k in range(t):
        assign[f"r_{k}_ancX0"] = False
        assign[f"r_{k}_flagZ0"] = False
    assign[f"r_{t}_ancX0"] = bits[0]
    assign[f"r_{t}_flagZ0"] = bits[1]
    for j in range(4):
        assign[f"r_{t + 1}_ancX{j}"] = bits[2 + j]
    return assign


def _eval_dec_formula(
    dec_name: str, sexpr: str, meas_names: List[str], assign: Dict[str, bool]
) -> bool:
    decls = {m: Bool(m) for m in meas_names}
    decls[dec_name] = Bool(dec_name)
    smt = f"(assert (= {dec_name} {sexpr}))"
    asts = parse_smt2_string(smt, decls=decls)
    if not asts:
        raise ValueError(f"parse_smt2_string returned empty for {dec_name}")
    eq = asts[0]
    rhs = eq.arg(1)
    subs = [(Bool(m), BoolVal(assign[m])) for m in meas_names]
    val = simplify(substitute(rhs, subs))
    return bool(is_true(val))


def load_log_then_stab_pairs(log_path: Path, stab_path: Path) -> List[Tuple[str, str]]:
    """Same order as proof_protocol proof: log lines then stab lines."""
    pairs: List[Tuple[str, str]] = []
    for p in (log_path, stab_path):
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            xs, zs = line.split()
            pairs.append((xs, zs))
    return pairs


def fixed_commutes_with_all_generators(
    fixed_x: List[bool], fixed_z: List[bool], pairs: List[Tuple[str, str]]
) -> bool:
    n = len(fixed_x)
    for xs, zs in pairs:
        if len(xs) != n or len(zs) != n:
            raise ValueError(f"Generator length mismatch: n_data={n}, xs={len(xs)}")
        acc = False
        for j in range(n):
            if zs[j] == "1":
                acc ^= fixed_x[j]
            if xs[j] == "1":
                acc ^= fixed_z[j]
        if acc:
            return False
    return True


def verify_decoder_commute(
    *,
    first_stabilizer_index: int,
    bitstring6: str,
    data_x: List[bool],
    data_z: List[bool],
    decoder_c_path: Path,
    log_path: Path,
    stab_path: Path,
) -> Tuple[bool, Dict[str, bool]]:
    """
    Returns (all_commute_holds, dec_values_by_name).
    """
    meas_names, dec_formulas = parse_decoder_c_file(decoder_c_path)
    assign = measurement_assignment_from_bitstring6(first_stabilizer_index, bitstring6)
    for m in meas_names:
        if m not in assign:
            raise KeyError(
                f"Decoder expects meas var {m!r} but assignment builder has no key "
                f"(first_stab={first_stabilizer_index})"
            )
    dec_vals: Dict[str, bool] = {}
    for dec_name, sexpr in sorted(dec_formulas.items()):
        dec_vals[dec_name] = _eval_dec_formula(dec_name, sexpr, meas_names, assign)

    n_data = len(data_x)
    if len(data_z) != n_data:
        raise ValueError("data_x and data_z length mismatch")
    fixed_x: List[bool] = []
    fixed_z: List[bool] = []
    for j in range(n_data):
        decx = dec_vals.get(f"dec{j}_x", False)
        decz = dec_vals.get(f"dec{j}_z", False)
        fixed_x.append(data_x[j] ^ decx)
        fixed_z.append(data_z[j] ^ decz)

    pairs = load_log_then_stab_pairs(log_path, stab_path)
    ok = fixed_commutes_with_all_generators(fixed_x, fixed_z, pairs)
    return ok, dec_vals
