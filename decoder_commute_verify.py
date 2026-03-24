"""
Check whether a learned decoder (decoder_C_{i}.txt) fixes the data error for a given
LUT syndrome pattern, using the same commute constraints as proof_protocol.py:
  fixed_j_x = data_j_x XOR dec_j_x, fixed_j_z = data_j_z XOR dec_j_z
and require fixed to commute with every line in log_txt then stab_txt (symplectic pairs).
"""
from __future__ import annotations

import builtins
import re
from pathlib import Path
from typing import Dict, List, Tuple

# Prefer bundled z3 shared library path (important for frozen executables).
_HERE = Path(__file__).resolve().parent
_BUNDLED_Z3_LIB_DIR = _HERE / "z3" / "lib"
if _BUNDLED_Z3_LIB_DIR.is_dir():
    builtins.Z3_LIB_DIRS = [str(_BUNDLED_Z3_LIB_DIR)]

from z3 import Bool, BoolVal, is_true, parse_smt2_string, simplify, substitute

DEC_LINE_START = re.compile(r"^(dec\d+_[xz])\t(.*)$")
_DECODER_CACHE: Dict[str, Tuple[List[str], Dict[str, str]]] = {}
_PAIR_CACHE: Dict[Tuple[str, str], List[Tuple[str, str]]] = {}


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


def _load_pairs_from_file(path: Path) -> List[Tuple[str, str]]:
    pairs: List[Tuple[str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        xs, zs = line.split()
        pairs.append((xs, zs))
    return pairs


def load_log_then_stab_pairs(log_path: Path, stab_path: Path) -> List[Tuple[str, str]]:
    """Same order as proof_protocol proof: log lines then stab lines."""
    return _load_pairs_from_file(log_path) + _load_pairs_from_file(stab_path)


def preload_decoder_assets(
    *,
    decoder_paths: List[Path],
    log_path: Path,
    stab_path: Path,
) -> None:
    """Preload decoder formulas and log/stab pairs into process-local cache."""
    for p in decoder_paths:
        key = str(p.resolve())
        if key not in _DECODER_CACHE:
            _DECODER_CACHE[key] = parse_decoder_c_file(p)
    pair_key = (str(log_path.resolve()), str(stab_path.resolve()))
    if pair_key not in _PAIR_CACHE:
        _PAIR_CACHE[pair_key] = load_log_then_stab_pairs(log_path, stab_path)


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
    dec_key = str(decoder_c_path.resolve())
    cached_dec = _DECODER_CACHE.get(dec_key)
    if cached_dec is None:
        cached_dec = parse_decoder_c_file(decoder_c_path)
        _DECODER_CACHE[dec_key] = cached_dec
    meas_names, dec_formulas = cached_dec
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

    pair_key = (str(log_path.resolve()), str(stab_path.resolve()))
    pairs = _PAIR_CACHE.get(pair_key)
    if pairs is None:
        pairs = load_log_then_stab_pairs(log_path, stab_path)
        _PAIR_CACHE[pair_key] = pairs
    ok = fixed_commutes_with_all_generators(fixed_x, fixed_z, pairs)
    return ok, dec_vals


def symplectic_inner_product_parity(
    fixed_x: List[bool], fixed_z: List[bool], xs: str, zs: str
) -> bool:
    """
    Symplectic inner product <fixed, (xs,zs)> mod 2.
    True => anticommute with that generator; False => commute.
    """
    n = len(fixed_x)
    acc = False
    for j in range(n):
        if zs[j] == "1":
            acc ^= fixed_x[j]
        if xs[j] == "1":
            acc ^= fixed_z[j]
    return acc


def format_decoder_commute_fail_report(
    *,
    first_stabilizer_index: int,
    bitstring6: str,
    data_x: List[bool],
    data_z: List[bool],
    decoder_c_path: Path,
    log_path: Path,
    stab_path: Path,
    fault_detail_lines: Optional[List[str]] = None,
    smt_fault_block: Optional[str] = None,
    stopped_after_first_fail: bool = False,
) -> str:
    """
    Human-readable report when decoder commute FAILs (same checks as verify_decoder_commute).
    fault_detail_lines: optional lines (no leading indent) under "fault (2q gate injection...)".
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

    log_pairs = _load_pairs_from_file(log_path)
    stab_pairs = _load_pairs_from_file(stab_path)

    dec_x_list = [int(dec_vals.get(f"dec{j}_x", False)) for j in range(n_data)]
    dec_z_list = [int(dec_vals.get(f"dec{j}_z", False)) for j in range(n_data)]
    d_x_list = [int(x) for x in data_x]
    d_z_list = [int(z) for z in data_z]
    fx_list = [int(x) for x in fixed_x]
    fz_list = [int(z) for z in fixed_z]

    lines: List[str] = []
    lines.append(f"decoder FAIL  file={decoder_c_path.name}")
    lines.append(f"  first_stabilizer_index={first_stabilizer_index}  bitstring6={bitstring6}")

    if fault_detail_lines:
        lines.append("  fault (2q gate injection for this DFS path):")
        for ln in fault_detail_lines:
            lines.append(f"    {ln}")

    lines.append("  meas vars (order as in decoder_C_*.txt), one per line:")
    for m in meas_names:
        lines.append(f"    {m} = {int(assign[m])}")
    lines.append(f"  dec_x = {dec_x_list}")
    lines.append(f"  dec_z = {dec_z_list}")
    lines.append(f"  data_x = {d_x_list}")
    lines.append(f"  data_z = {d_z_list}")
    lines.append(f"  fixed_x (list) = {fx_list}  fixed_z (list) = {fz_list}")

    lines.append(
        "  commute vs generators (log_txt lines then stab_txt lines; "
        "commute_ok True => Pauli commutes):"
    )
    for i, (xs, zs) in enumerate(log_pairs):
        ac = symplectic_inner_product_parity(fixed_x, fixed_z, xs, zs)
        ip = int(ac)
        ok_g = not ac
        kind = "anticommute" if ac else "commute"
        lines.append(
            f"    log[{i}] ({xs} {zs}): inner_parity={ip}  {kind}  ok={int(ok_g)}"
        )
    for i, (xs, zs) in enumerate(stab_pairs):
        ac = symplectic_inner_product_parity(fixed_x, fixed_z, xs, zs)
        ip = int(ac)
        ok_g = not ac
        kind = "anticommute" if ac else "commute"
        lines.append(
            f"    stab[{i}] ({xs} {zs}): inner_parity={ip}  {kind}  ok={int(ok_g)}"
        )

    if smt_fault_block:
        lines.append("")
        lines.append("  SMT fault asserts (paste after declare-fun block in output_C.txt):")
        for sl in smt_fault_block.splitlines():
            lines.append(f"    {sl}")

    if stopped_after_first_fail:
        lines.append("Stopped after first decoder commute FAIL (partial enumeration).")

    return "\n".join(lines)
