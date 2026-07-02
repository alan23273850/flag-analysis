"""
Projected weighted model counting (PWMC) pipeline for [[5,1,3]] physical error
rate. Computes the same quantity as monte_carlo_table_513.py:

    ratio_fail = decoder_fail / total_runs

but via exact summation over weighted Boolean models using `gpmc -mode=3`.

For each "first triggering stabilizer" t in {0, 1, 2, 3} we:
  1. Symbolically simulate sim_513.run_protocol_once but with fault indicator
     Boolean variables (one event per fault site, mutually exclusive).
  2. Build the failure event:
         path_condition(t) AND NOT all_commute(data XOR decoder)
     using the previously-learned decoder formulas in
     `decoder/<backend>/path_{4 - t}.txt` and the log/stab generators.
  3. Tseitin-convert to CNF, emit DIMACS+weights in MCC2021 format.
  4. Invoke `gpmc -mode=3` and parse the weight result.
The final estimate is the sum across t in 0..3.

This module is self-contained: it does not depend on
proof_protocol_boolean. It only reuses the QASM gate evolution helpers from
flag_analysis and the decoder-file parser from decoder_commute_verify.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from qiskit import QuantumCircuit
from z3 import (
    And,
    Bool,
    BoolRef,
    BoolVal,
    Goal,
    Not,
    Or,
    Then,
    Xor,
    is_const,
    is_false,
    is_not,
    is_true,
    parse_smt2_string,
    simplify,
    substitute,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from flag_analysis import (  # noqa: E402
    CircuitXZ,
    QubitXZ,
    apply_qasm_gate_into_state,
    detect_qubit_groups,
    new_clean_circuit_state,
)
from decoder_commute_verify import decoder_dir, parse_decoder_c_file  # noqa: E402

ORIGIN_DIR = PROJECT_ROOT / "[[5,1,3]]_origin"
FLAG_QASM = ORIGIN_DIR / "[[5,1,3]]_origin_flag_syndrome.qasm "
RAW_QASM = ORIGIN_DIR / "raw_[[5,1,3]]_origin_flag_syndrome.qasm"
STAB_TXT = ORIGIN_DIR / "[[5,1,3]]_origin.txt"
LOG_TXT = ORIGIN_DIR / "[[5,1,3]]_log_op.txt"

GPMC_BIN = os.environ.get("GPMC_BIN", "gpmc")


# ---------------------------------------------------------------------------
# Fault-site data structures and helpers
# ---------------------------------------------------------------------------

@dataclass
class FaultSite:
    """One fault site with a compact 'trigger + free flip bits' encoding.

    Encoding (`kind` -> structure):
      - "2q":  show vars = [t, fx0, fz0, fx1, fz1] (5 vars), weighted lits =
        {t: (p/15, 1-p)} and unit-weight free bits. The constraint
            t <-> (fx0 OR fz0 OR fx1 OR fz1)
        in `constraint` couples the bits to t so the joint distribution over
        the 5 vars reproduces sim_513's 2q fault model exactly:
            P(t=F, all bits=0) = 1 - p
            P(t=T, bits = any of 15 non-zero patterns) = p/15 each
      - "1q_depol":  show vars = [t, fx, fz] (3 vars), weighted lits =
        {t: (tp/3, 1-tp)}, constraint `t <-> (fx OR fz)`. Joint:
            P(t=F, fx=fz=F) = 1 - tp
            P(t=T, (fx,fz) in {(T,F),(T,T),(F,T)}) = tp/3 each
      - "prep":  show vars = [f] (1 var), weighted lit = (2p/3, 1-2p/3),
        no constraint (single Bernoulli).

    `weighted_lits` maps each show var (BoolRef) to (w_true, w_false).
    `free_bits` are the show vars whose literal weights are both 1 (omitted
    from the explicit weight comments).
    `flip_bits` are the bits that drive Pauli flips in the symbolic state
    (subset of `free_bits` for 2q/1q, or `[trigger]` for prep).
    """

    site_id: str
    kind: str
    trigger: Optional[BoolRef]
    free_bits: List[BoolRef]
    flip_bits: List[BoolRef]
    weighted_lits: Dict[str, Tuple[float, float]]
    constraint: BoolRef

    def all_show_vars(self) -> List[BoolRef]:
        out: List[BoolRef] = []
        if self.trigger is not None:
            out.append(self.trigger)
        out.extend(self.free_bits)
        return out


def _build_reg_index_map(qc: QuantumCircuit) -> Dict[Tuple[str, int], int]:
    out: Dict[Tuple[str, int], int] = {}
    g = 0
    for qreg in qc.qregs:
        for j in range(qreg.size):
            out[(qreg.name, j)] = g
            g += 1
    return out


def _make_bool(name: str) -> BoolRef:
    return Bool(name)


def _add_2q_fault(
    state: CircuitXZ,
    q0: int,
    q1: int,
    *,
    p: float,
    site_id: str,
) -> FaultSite:
    """Compact 5-var encoding for 2q gate fault with sim_513 distribution."""
    trig = _make_bool(f"t_{site_id}")
    fx0 = _make_bool(f"b_{site_id}_fx0")
    fz0 = _make_bool(f"b_{site_id}_fz0")
    fx1 = _make_bool(f"b_{site_id}_fx1")
    fz1 = _make_bool(f"b_{site_id}_fz1")
    state.qubits[q0].x = Xor(state.qubits[q0].x, fx0)
    state.qubits[q0].z = Xor(state.qubits[q0].z, fz0)
    state.qubits[q1].x = Xor(state.qubits[q1].x, fx1)
    state.qubits[q1].z = Xor(state.qubits[q1].z, fz1)
    any_bit = Or(fx0, fz0, fx1, fz1)
    constraint = And(Or(Not(trig), any_bit), Or(trig, Not(any_bit)))
    weights = {str(trig): (p / 15.0, 1.0 - p)}
    return FaultSite(
        site_id=site_id,
        kind="2q",
        trigger=trig,
        free_bits=[fx0, fz0, fx1, fz1],
        flip_bits=[fx0, fz0, fx1, fz1],
        weighted_lits=weights,
        constraint=constraint,
    )


def _add_1q_depol_fault(
    state: CircuitXZ,
    q: int,
    *,
    trigger_p: float,
    site_id: str,
) -> FaultSite:
    """Compact 3-var encoding for 1q depolarizing fault."""
    tp = max(0.0, min(1.0, trigger_p))
    trig = _make_bool(f"t_{site_id}")
    fx = _make_bool(f"b_{site_id}_fx")
    fz = _make_bool(f"b_{site_id}_fz")
    state.qubits[q].x = Xor(state.qubits[q].x, fx)
    state.qubits[q].z = Xor(state.qubits[q].z, fz)
    any_bit = Or(fx, fz)
    constraint = And(Or(Not(trig), any_bit), Or(trig, Not(any_bit)))
    weights = {str(trig): (tp / 3.0, 1.0 - tp)}
    return FaultSite(
        site_id=site_id,
        kind="1q_depol",
        trigger=trig,
        free_bits=[fx, fz],
        flip_bits=[fx, fz],
        weighted_lits=weights,
        constraint=constraint,
    )


def _add_prep_fault(
    state: CircuitXZ,
    q: int,
    *,
    axis: str,
    p: float,
    site_id: str,
) -> FaultSite:
    """Single-Bernoulli encoding for prep fault on a fixed axis."""
    tp = min(1.0, (2.0 * p) / 3.0)
    flip = _make_bool(f"f_{site_id}")
    if axis == "x":
        state.qubits[q].x = Xor(state.qubits[q].x, flip)
    elif axis == "z":
        state.qubits[q].z = Xor(state.qubits[q].z, flip)
    else:
        raise ValueError(f"axis must be x or z, got {axis}")
    weights = {str(flip): (tp, 1.0 - tp)}
    return FaultSite(
        site_id=site_id,
        kind="prep",
        trigger=None,
        free_bits=[flip],
        flip_bits=[flip],
        weighted_lits=weights,
        constraint=BoolVal(True),
    )


# ---------------------------------------------------------------------------
# Symbolic protocol simulation per path
# ---------------------------------------------------------------------------

@dataclass
class PathSymbolic:
    t: int
    fault_sites: List[FaultSite]
    meas_exprs: Dict[str, BoolRef]
    data_x: List[BoolRef]
    data_z: List[BoolRef]
    path_condition: BoolRef


def _sym_run_gate_slice(
    qc: QuantumCircuit,
    state: CircuitXZ,
    start: int,
    end: int,
    *,
    base_p: float,
    gamma: float,
    phase: str,
    fault_sites: List[FaultSite],
) -> None:
    participating = set()
    for instr, qargs, _ in qc.data[start:end]:
        _ = instr
        participating.update(qc.find_bit(q).index for q in qargs)

    for gidx in range(start, end):
        instr, qargs, _ = qc.data[gidx]
        name = instr.name
        qidxs = [qc.find_bit(q).index for q in qargs]
        apply_qasm_gate_into_state(state, name, qidxs)
        if name in ("cx", "cnot", "cz"):
            site = _add_2q_fault(
                state,
                qidxs[0],
                qidxs[1],
                p=base_p,
                site_id=f"{phase}_g{gidx}_2q",
            )
            fault_sites.append(site)
        active = set(qidxs)
        for qi in participating:
            if qi in active:
                continue
            site = _add_1q_depol_fault(
                state,
                qi,
                trigger_p=gamma * base_p,
                site_id=f"{phase}_g{gidx}_idle_q{qi}",
            )
            fault_sites.append(site)


def _copy_data_error_only_sym(
    src_state: CircuitXZ,
    src_qc: QuantumCircuit,
    dst_qc: QuantumCircuit,
) -> CircuitXZ:
    src_groups = detect_qubit_groups(src_qc)
    dst_groups = detect_qubit_groups(dst_qc)
    src_data = sorted(src_groups["data"])
    dst_data = sorted(dst_groups["data"])
    if len(src_data) != len(dst_data):
        raise ValueError("data qubit count mismatch")
    dst_state = new_clean_circuit_state(dst_qc.num_qubits)
    for si, di in zip(src_data, dst_data):
        dst_state.qubits[di].x = src_state.qubits[si].x
        dst_state.qubits[di].z = src_state.qubits[si].z
    return dst_state


def build_path_symbolic(
    t: int,
    *,
    p: float,
    gamma: float,
    beta: float,
) -> PathSymbolic:
    """Build symbolic state for the case where stab `t` is the first to trigger."""
    if not (0 <= t <= 3):
        raise ValueError("t must be in 0..3")

    flag_qc = QuantumCircuit.from_qasm_file(str(FLAG_QASM))
    raw_qc = QuantumCircuit.from_qasm_file(str(RAW_QASM))

    first_state = new_clean_circuit_state(flag_qc.num_qubits)
    fault_sites: List[FaultSite] = []
    meas_exprs: Dict[str, BoolRef] = {}

    idx_flag = _build_reg_index_map(flag_qc)
    anc_idx = idx_flag[("ancX", 0)]
    flag_idx = idx_flag[("flagZ", 0)]

    fault_sites.append(
        _add_prep_fault(
            first_state, anc_idx, axis="z", p=p,
            site_id="prep_first_ancX0",
        )
    )
    fault_sites.append(
        _add_prep_fault(
            first_state, flag_idx, axis="x", p=p,
            site_id="prep_first_flagZ0",
        )
    )

    n_gates = len(flag_qc.data)
    if n_gates % 4 != 0:
        raise ValueError("flag QASM gate count not divisible by 4")
    gates_per_stab = n_gates // 4

    for stab_i in range(t + 1):
        start_g = stab_i * gates_per_stab
        end_g = (stab_i + 1) * gates_per_stab
        _sym_run_gate_slice(
            flag_qc,
            first_state,
            start_g,
            end_g,
            base_p=p,
            gamma=gamma,
            phase=f"first_stab{stab_i}",
            fault_sites=fault_sites,
        )
        fault_sites.append(
            _add_1q_depol_fault(
                first_state, anc_idx,
                trigger_p=beta * p,
                site_id=f"meas_first_stab{stab_i}_ancX0",
            )
        )
        fault_sites.append(
            _add_1q_depol_fault(
                first_state, flag_idx,
                trigger_p=beta * p,
                site_id=f"meas_first_stab{stab_i}_flagZ0",
            )
        )
        meas_exprs[f"r_{stab_i}_ancX0"] = first_state.qubits[anc_idx].z
        meas_exprs[f"r_{stab_i}_flagZ0"] = first_state.qubits[flag_idx].x

    path_terms: List[BoolRef] = []
    for i in range(t):
        path_terms.append(Not(meas_exprs[f"r_{i}_ancX0"]))
        path_terms.append(Not(meas_exprs[f"r_{i}_flagZ0"]))
    path_terms.append(
        Or(meas_exprs[f"r_{t}_ancX0"], meas_exprs[f"r_{t}_flagZ0"])
    )
    path_condition = And(*path_terms) if len(path_terms) > 1 else path_terms[0]

    second_state = _copy_data_error_only_sym(first_state, flag_qc, raw_qc)
    idx_raw = _build_reg_index_map(raw_qc)
    raw_anc = [idx_raw[("ancX", j)] for j in range(4)]

    n_raw_gates = len(raw_qc.data)
    if n_raw_gates % 4 != 0:
        raise ValueError("raw QASM gate count not divisible by 4")
    raw_gates_per_stab = n_raw_gates // 4

    for j in range(4):
        start_g = j * raw_gates_per_stab
        end_g = (j + 1) * raw_gates_per_stab
        second_state.qubits[raw_anc[j]].x = BoolVal(False)
        second_state.qubits[raw_anc[j]].z = BoolVal(False)
        fault_sites.append(
            _add_prep_fault(
                second_state, raw_anc[j], axis="z", p=p,
                site_id=f"prep_second_stab{j}_ancX{j}",
            )
        )
        _sym_run_gate_slice(
            raw_qc,
            second_state,
            start_g,
            end_g,
            base_p=p,
            gamma=gamma,
            phase=f"second_stab{j}",
            fault_sites=fault_sites,
        )
        fault_sites.append(
            _add_1q_depol_fault(
                second_state, raw_anc[j],
                trigger_p=beta * p,
                site_id=f"meas_second_stab{j}_ancX{j}",
            )
        )
        meas_exprs[f"r_{t + 1}_ancX{j}"] = second_state.qubits[raw_anc[j]].z

    groups = detect_qubit_groups(flag_qc)
    data_idxs = sorted(groups["data"])
    data_x = [first_state.qubits[i].x for i in data_idxs]
    data_z = [first_state.qubits[i].z for i in data_idxs]

    return PathSymbolic(
        t=t,
        fault_sites=fault_sites,
        meas_exprs=meas_exprs,
        data_x=data_x,
        data_z=data_z,
        path_condition=path_condition,
    )


# ---------------------------------------------------------------------------
# Decoder + commute formula
# ---------------------------------------------------------------------------

def _load_pairs_from_file(path: Path) -> List[Tuple[str, str]]:
    pairs: List[Tuple[str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        xs, zs = line.split()
        pairs.append((xs, zs))
    return pairs


def _load_log_then_stab_pairs() -> List[Tuple[str, str]]:
    return _load_pairs_from_file(LOG_TXT) + _load_pairs_from_file(STAB_TXT)


def _parse_dec_formula(
    dec_name: str,
    sexpr: str,
    meas_names: Sequence[str],
) -> BoolRef:
    """Parse one decoder formula into a Z3 BoolRef over Bool(meas_name)."""
    decls = {m: Bool(m) for m in meas_names}
    decls[dec_name] = Bool(dec_name)
    smt = f"(assert (= {dec_name} {sexpr}))"
    asts = parse_smt2_string(smt, decls=decls)
    if not asts:
        raise ValueError(f"empty parse for {dec_name}")
    eq = asts[0]
    return eq.arg(1)


def build_failure_formula(
    path: PathSymbolic,
) -> Tuple[BoolRef, List[FaultSite]]:
    """Build the fail formula and return (formula, fault_sites_used)."""
    decoder_idx = 4 - path.t
    dec_file = decoder_path(ORIGIN_DIR, decoder_idx)
    meas_names, dec_formulas = parse_decoder_c_file(dec_file)

    subst_pairs: List[Tuple[BoolRef, BoolRef]] = []
    for m in meas_names:
        if m in path.meas_exprs:
            subst_pairs.append((Bool(m), path.meas_exprs[m]))
        else:
            mo = re.match(r"r_(\d+)_(ancX0|flagZ0)$", m)
            if mo is None:
                raise KeyError(f"decoder expects unknown meas var: {m}")
            r = int(mo.group(1))
            if r < path.t:
                subst_pairs.append((Bool(m), BoolVal(False)))
            else:
                raise KeyError(f"decoder expects meas var not built: {m}")

    n_data = len(path.data_x)
    fixed_x: List[BoolRef] = list(path.data_x)
    fixed_z: List[BoolRef] = list(path.data_z)
    for j in range(n_data):
        dx_name = f"dec{j}_x"
        dz_name = f"dec{j}_z"
        dx_expr = _parse_dec_formula(dx_name, dec_formulas[dx_name], meas_names)
        dz_expr = _parse_dec_formula(dz_name, dec_formulas[dz_name], meas_names)
        dx_expr = substitute(dx_expr, subst_pairs)
        dz_expr = substitute(dz_expr, subst_pairs)
        fixed_x[j] = Xor(fixed_x[j], dx_expr)
        fixed_z[j] = Xor(fixed_z[j], dz_expr)

    pairs = _load_log_then_stab_pairs()
    commute_terms: List[BoolRef] = []
    for xs, zs in pairs:
        if len(xs) != n_data or len(zs) != n_data:
            raise ValueError("generator length mismatch with data qubit count")
        acc = BoolVal(False)
        for j in range(n_data):
            if zs[j] == "1":
                acc = Xor(acc, fixed_x[j])
            if xs[j] == "1":
                acc = Xor(acc, fixed_z[j])
        commute_terms.append(Not(acc))
    if not commute_terms:
        all_commute = BoolVal(True)
    elif len(commute_terms) == 1:
        all_commute = commute_terms[0]
    else:
        all_commute = And(*commute_terms)

    fail = And(path.path_condition, Not(all_commute))
    return fail, path.fault_sites


# ---------------------------------------------------------------------------
# CNF + DIMACS export
# ---------------------------------------------------------------------------

def _cnf_from_formula(formula: BoolRef) -> List[BoolRef]:
    g = Goal()
    g.add(formula)
    t = Then("simplify", "tseitin-cnf")
    res = t(g)
    if len(res) != 1:
        raise RuntimeError(f"unexpected number of result goals: {len(res)}")
    ng = res[0]
    return list(ng)


def _atoms_in_clauses(clauses: Sequence[BoolRef]) -> List[str]:
    seen: Dict[str, None] = {}

    def walk(e: BoolRef) -> None:
        if is_const(e):
            if e.num_args() == 0 and not (is_true(e) or is_false(e)):
                seen.setdefault(str(e), None)
            return
        for ch in e.children():
            walk(ch)

    for cl in clauses:
        walk(cl)
    return list(seen.keys())


def _literal_from_clause_child(child: BoolRef) -> Tuple[bool, str]:
    """Return (is_positive, atom_name) for one clause literal."""
    if is_not(child):
        inner = child.arg(0)
        if not is_const(inner) or inner.num_args() != 0:
            raise ValueError(f"unexpected literal: {child}")
        return False, str(inner)
    if is_const(child) and child.num_args() == 0:
        return True, str(child)
    raise ValueError(f"unexpected clause child: {child}")


def _clause_literals(cl: BoolRef) -> List[Tuple[bool, str]]:
    if is_true(cl):
        return []
    if is_false(cl):
        return [(True, "__FALSE_CLAUSE__")]
    if cl.decl().name() == "or":
        return [_literal_from_clause_child(c) for c in cl.children()]
    return [_literal_from_clause_child(cl)]


def _fault_active_indicator(site: FaultSite) -> BoolRef:
    """Boolean that is True iff this fault site fires (any non-identity event)."""
    if site.kind in ("2q", "1q_depol"):
        if site.trigger is None:
            raise ValueError(f"site {site.site_id} missing trigger")
        return site.trigger
    if site.kind == "prep":
        if not site.free_bits:
            raise ValueError(f"site {site.site_id} missing flip bit")
        return site.free_bits[0]
    raise ValueError(f"unknown fault kind {site.kind}")


def _sinz_at_most_k(vars_: Sequence[BoolRef], k: int, prefix: str) -> BoolRef:
    """Sequential-counter encoding of at-most-k(vars_).

    Encodes via O(n*k) auxiliary registers. Returns a Z3 expression that is
    True iff at most k of the variables are True.
    """
    n = len(vars_)
    if k < 0:
        return BoolVal(False)
    if k >= n:
        return BoolVal(True)
    if k == 0:
        return And(*[Not(v) for v in vars_])
    s = [[Bool(f"{prefix}_s{i}_{j}") for j in range(k)] for i in range(n)]
    clauses: List[BoolRef] = []
    clauses.append(Or(Not(vars_[0]), s[0][0]))
    for j in range(1, k):
        clauses.append(Not(s[0][j]))
    for i in range(1, n):
        clauses.append(Or(Not(vars_[i]), s[i][0]))
        clauses.append(Or(Not(s[i - 1][0]), s[i][0]))
        for j in range(1, k):
            clauses.append(Or(Not(vars_[i]), Not(s[i - 1][j - 1]), s[i][j]))
            clauses.append(Or(Not(s[i - 1][j]), s[i][j]))
        clauses.append(Or(Not(vars_[i]), Not(s[i - 1][k - 1])))
    if not clauses:
        return BoolVal(True)
    return And(*clauses)


def emit_dimacs(
    formula: BoolRef,
    fault_sites: Sequence[FaultSite],
    out_path: Path,
    *,
    max_faults: Optional[int] = None,
) -> int:
    """Convert to CNF and write a MCC2021-style DIMACS file. Returns clause count.

    `max_faults`: if set, add an at-most-K cardinality constraint over the
    per-site `fault_active_indicator` Booleans. This truncates the WMC at the
    p^{K+1} order; for typical sim_513 parameters K=2 already gives errors
    below 1e-12.
    """
    constraint_terms: List[BoolRef] = []
    for site in fault_sites:
        c = site.constraint
        if not is_true(c):
            constraint_terms.append(c)
    if max_faults is not None:
        actives = [_fault_active_indicator(s) for s in fault_sites]
        card_expr = _sinz_at_most_k(actives, int(max_faults), prefix="card")
        if not is_true(card_expr):
            constraint_terms.append(card_expr)
    full = formula if not constraint_terms else And(formula, *constraint_terms)

    clauses = _cnf_from_formula(full)

    cleaned: List[BoolRef] = []
    for cl in clauses:
        if is_true(cl):
            continue
        cleaned.append(cl)

    atom_names = _atoms_in_clauses(cleaned)

    show_names: List[str] = []
    seen_show: set = set()
    for site in fault_sites:
        for v in site.all_show_vars():
            nm = str(v)
            if nm in seen_show:
                continue
            seen_show.add(nm)
            show_names.append(nm)

    weights_by_name: Dict[str, Tuple[float, float]] = {}
    for site in fault_sites:
        for nm, wts in site.weighted_lits.items():
            weights_by_name[nm] = wts

    ordered_atoms: List[str] = []
    seen_ord: set = set()
    for nm in show_names:
        if nm in atom_names and nm not in seen_ord:
            ordered_atoms.append(nm)
            seen_ord.add(nm)
    for nm in atom_names:
        if nm not in seen_ord:
            ordered_atoms.append(nm)
            seen_ord.add(nm)
    var_id: Dict[str, int] = {nm: i + 1 for i, nm in enumerate(ordered_atoms)}

    cnf_lines: List[str] = []
    for cl in cleaned:
        if is_false(cl):
            cnf_lines.append("0\n")
            continue
        lits = _clause_literals(cl)
        parts: List[str] = []
        for pos, nm in lits:
            if nm == "__FALSE_CLAUSE__":
                continue
            v = var_id[nm]
            parts.append(str(v if pos else -v))
        cnf_lines.append(" ".join(parts) + " 0\n")

    header_lines: List[str] = []
    header_lines.append("c [[5,1,3]] WMC instance (compact fault encoding)\n")
    show_vars = [str(var_id[nm]) for nm in show_names if nm in var_id]
    if show_vars:
        header_lines.append("c p show " + " ".join(show_vars) + " 0\n")
    for nm in show_names:
        if nm not in var_id:
            continue
        v = var_id[nm]
        wt, wf = weights_by_name.get(nm, (1.0, 1.0))
        header_lines.append(f"c p weight {v} {wt:.17g} 0\n")
        header_lines.append(f"c p weight {-v} {wf:.17g} 0\n")

    n_vars = len(ordered_atoms)
    n_clauses = len(cnf_lines)
    with out_path.open("w", encoding="utf-8") as f:
        f.write(f"p cnf {n_vars} {n_clauses}\n")
        for line in header_lines:
            f.write(line)
        for line in cnf_lines:
            f.write(line)
    return n_clauses


# ---------------------------------------------------------------------------
# gpmc invocation + parsing
# ---------------------------------------------------------------------------

_GPMC_EXACT_RE = re.compile(r"^c s exact\s+\S+\s+\S+\s+(\S+)\s*$")
_GPMC_LOG10_PREFIX = "c s log10-estimate"


def run_gpmc(
    cnf_path: Path,
    *,
    mode: int = 3,
    extra_args: Sequence[str] = (),
    timeout: Optional[float] = None,
) -> float:
    """Run gpmc and parse the exact weight from its stdout.

    Looks for the MCC-style line
        c s exact <double|arb> <prec-sci|int|float> <value>
    and returns the value as float. Falls back to nothing if not present.
    """
    cmd = [GPMC_BIN, f"-mode={mode}"]
    cmd.extend(extra_args)
    cmd.append(str(cnf_path))
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    text = proc.stdout + "\n" + proc.stderr
    weight: Optional[float] = None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith(_GPMC_LOG10_PREFIX):
            continue
        m = _GPMC_EXACT_RE.match(s)
        if m:
            try:
                weight = float(m.group(1))
            except ValueError:
                weight = None
    if weight is None:
        sys.stderr.write(text)
        raise RuntimeError(f"could not parse gpmc weight (rc={proc.returncode})")
    return weight


# ---------------------------------------------------------------------------
# High-level driver
# ---------------------------------------------------------------------------

@dataclass
class PathResult:
    t: int
    fail_prob: float
    n_vars: int
    n_clauses: int
    elapsed_sec: float
    cnf_path: Path


def compute_fail_prob_for_path(
    t: int,
    *,
    p: float,
    gamma: float,
    beta: float,
    cnf_path: Path,
    gpmc_extra_args: Sequence[str] = (),
    keep_cnf: bool = True,
    max_faults: Optional[int] = None,
) -> PathResult:
    start = time.time()
    sym = build_path_symbolic(t, p=p, gamma=gamma, beta=beta)
    fail_formula, sites = build_failure_formula(sym)
    n_clauses = emit_dimacs(fail_formula, sites, cnf_path, max_faults=max_faults)
    with cnf_path.open() as f:
        header = f.readline().strip().split()
        n_vars = int(header[2]) if len(header) >= 4 else -1
    prob = run_gpmc(cnf_path, extra_args=gpmc_extra_args)
    elapsed = time.time() - start
    if not keep_cnf:
        try:
            cnf_path.unlink()
        except FileNotFoundError:
            pass
    return PathResult(
        t=t,
        fail_prob=prob,
        n_vars=n_vars,
        n_clauses=n_clauses,
        elapsed_sec=elapsed,
        cnf_path=cnf_path,
    )


@dataclass
class WMCResult:
    p: float
    beta: float
    gamma: float
    paths: List[PathResult]
    ratio_fail: float
    total_elapsed_sec: float


def compute_ratio_fail(
    *,
    p: float,
    gamma: float,
    beta: float,
    work_dir: Optional[Path] = None,
    gpmc_extra_args: Sequence[str] = (),
    keep_cnf: bool = True,
    max_faults: Optional[int] = None,
) -> WMCResult:
    if work_dir is None:
        work_dir = Path(tempfile.mkdtemp(prefix="wmc_513_"))
    work_dir.mkdir(parents=True, exist_ok=True)
    start = time.time()
    paths: List[PathResult] = []
    for t in range(4):
        cnf_path = work_dir / f"path_t{t}_p{p:.8g}.cnf"
        pr = compute_fail_prob_for_path(
            t,
            p=p,
            gamma=gamma,
            beta=beta,
            cnf_path=cnf_path,
            gpmc_extra_args=gpmc_extra_args,
            keep_cnf=keep_cnf,
            max_faults=max_faults,
        )
        paths.append(pr)
    elapsed = time.time() - start
    return WMCResult(
        p=p,
        beta=beta,
        gamma=gamma,
        paths=paths,
        ratio_fail=sum(pp.fail_prob for pp in paths),
        total_elapsed_sec=elapsed,
    )


# ---------------------------------------------------------------------------
# CLI entry
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="WMC for [[5,1,3]] physical error rate")
    ap.add_argument("--p", type=float, required=True, help="2q gate fault rate")
    ap.add_argument("--gamma", type=float, default=1.0, help="idle factor (gamma)")
    ap.add_argument("--beta", type=float, default=1.0, help="meas factor (beta)")
    ap.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help="directory to write CNF files (defaults to tempdir)",
    )
    ap.add_argument(
        "--keep-cnf",
        action="store_true",
        help="keep CNF files after running gpmc",
    )
    ap.add_argument(
        "--no-keep-cnf",
        dest="keep_cnf",
        action="store_false",
        help="delete CNF files after running gpmc",
    )
    ap.set_defaults(keep_cnf=True)
    ap.add_argument(
        "--backend",
        choices=("cms", "dnf", "cadet"),
        default=None,
        help="Decoder backend subdirectory (default: DECODER_BACKEND env or cms)",
    )
    ap.add_argument(
        "--max-faults",
        type=int,
        default=None,
        help=(
            "Optional at-most-K cardinality constraint on the number of fired "
            "fault sites. Truncates WMC at the p^(K+1) order. Recommended K=2 "
            "for [[5,1,3]] (distance-3) physical error rates."
        ),
    )
    return ap.parse_args()


def main() -> None:
    args = _parse_args()
    if args.backend is not None:
        os.environ["DECODER_BACKEND"] = args.backend
    res = compute_ratio_fail(
        p=args.p,
        gamma=args.gamma,
        beta=args.beta,
        work_dir=args.work_dir,
        keep_cnf=args.keep_cnf,
        max_faults=args.max_faults,
    )
    print("[WMC] [[5,1,3]] result")
    for pr in res.paths:
        print(
            f"  path t={pr.t}: P(fail)={pr.fail_prob:.6g}  "
            f"vars={pr.n_vars}  clauses={pr.n_clauses}  "
            f"elapsed={pr.elapsed_sec:.2f}s"
        )
    print(
        f"[WMC] p={res.p:.8g} gamma={res.gamma:.8g} beta={res.beta:.8g} "
        f"ratio_fail={res.ratio_fail:.8g}  total={res.total_elapsed_sec:.2f}s"
    )


if __name__ == "__main__":
    main()
