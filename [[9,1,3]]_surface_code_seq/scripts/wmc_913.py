"""
PWMC pipeline for [[9,1,3]]_surface_code_seq (Chris d=3 flag protocol).

Computes logical error rate as:
    ratio_fail = sum_{path_idx in 0..3} P(fail | path_idx)

using learned decoders in decoder/<backend>/path_{path_idx}.txt and gpmc -mode=3.
Fault encoding and DIMACS export reuse [[5,1,3]] wmc_513.py.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from qiskit import QuantumCircuit
from z3 import And, Bool, BoolRef, BoolVal, Not, Or, Xor, parse_smt2_string, substitute

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

WMC_513_DIR = PROJECT_ROOT / "[[5,1,3]]_origin" / "scripts"
if str(WMC_513_DIR) not in sys.path:
    sys.path.insert(0, str(WMC_513_DIR))

import wmc_513 as wmc  # noqa: E402

from decoder_commute_verify import decoder_path, parse_decoder_c_file  # noqa: E402
from flag_analysis import CircuitXZ, apply_qasm_gate_into_state, detect_qubit_groups, new_clean_circuit_state  # noqa: E402

ORIGIN_DIR = Path(__file__).resolve().parents[1]
FLAG_QASM = ORIGIN_DIR / "[[9,1,3]]_flag_syndrome.qasm"
RAW_QASM = ORIGIN_DIR / "[[9,1,3]]_raw_syndrome.qasm"
STAB_TXT = ORIGIN_DIR / "[[9,1,3]].txt"
LOG_TXT = ORIGIN_DIR / "[[9,1,3]]_log_op.txt"

N_ANC = 4
N_FLAG = 2

FLAG_PREP = (
    ("ancX", "z"),
    ("ancZ", "x"),
    ("flagX", "x"),
    ("flagZ", "z"),
)
RAW_PREP = (("ancX", "z"), ("ancZ", "x"))


@dataclass
class PathSymbolic:
    path_idx: int
    fault_sites: List[wmc.FaultSite]
    meas_exprs: Dict[str, BoolRef]
    data_x: List[BoolRef]
    data_z: List[BoolRef]
    path_condition: BoolRef


def _build_reg_index_map(qc: QuantumCircuit) -> Dict[Tuple[str, int], int]:
    return wmc._build_reg_index_map(qc)


def _reset_qubits(state: CircuitXZ, qidxs: Sequence[int]) -> None:
    for qi in qidxs:
        state.qubits[qi].x = BoolVal(False)
        state.qubits[qi].z = BoolVal(False)


def _copy_data_errors(
    src: CircuitXZ, src_qc: QuantumCircuit, dst: CircuitXZ, dst_qc: QuantumCircuit
) -> None:
    src_g = detect_qubit_groups(src_qc)
    dst_g = detect_qubit_groups(dst_qc)
    src_data = sorted(src_g["data"])
    dst_data = sorted(dst_g["data"])
    if len(src_data) != len(dst_data):
        raise ValueError("data qubit count mismatch between circuits")
    for si, di in zip(src_data, dst_data):
        dst.qubits[di].x = src.qubits[si].x
        dst.qubits[di].z = src.qubits[si].z


def _new_state_carry_data(
    prev: CircuitXZ,
    prev_qc: QuantumCircuit,
    next_qc: QuantumCircuit,
    *,
    reset_groups: Tuple[str, ...],
) -> CircuitXZ:
    state = new_clean_circuit_state(next_qc.num_qubits)
    _copy_data_errors(prev, prev_qc, state, next_qc)
    groups = detect_qubit_groups(next_qc)
    reset_idxs: List[int] = []
    for g in reset_groups:
        reset_idxs.extend(groups.get(g, []))
    _reset_qubits(state, sorted(set(reset_idxs)))
    return state


def _sym_run_circuit(
    qc: QuantumCircuit,
    state: CircuitXZ,
    *,
    round_idx: int,
    base_p: float,
    gamma: float,
    prep_groups: Sequence[Tuple[str, str]],
) -> List[wmc.FaultSite]:
    fault_sites: List[wmc.FaultSite] = []
    groups = detect_qubit_groups(qc)

    for reg_name, axis in prep_groups:
        for j in groups.get(reg_name, []):
            fault_sites.append(
                wmc._add_prep_fault(
                    state,
                    j,
                    axis=axis,
                    p=base_p,
                    site_id=f"r{round_idx}_prep_{reg_name}{j}",
                )
            )

    participating: set = set()
    for instr, qargs, _ in qc.data:
        if instr.name not in ("barrier", "id", "reset", "measure"):
            participating.update(qc.find_bit(q).index for q in qargs)

    for gidx, (instr, qargs, _) in enumerate(qc.data):
        name = instr.name
        if name in ("barrier", "id", "reset", "measure"):
            continue
        qidxs = [qc.find_bit(q).index for q in qargs]
        apply_qasm_gate_into_state(state, name, qidxs)
        if name in ("cx", "cnot", "cz"):
            fault_sites.append(
                wmc._add_2q_fault(
                    state,
                    qidxs[0],
                    qidxs[1],
                    p=base_p,
                    site_id=f"r{round_idx}_g{gidx}_2q",
                )
            )
        active = set(qidxs)
        for qi in participating:
            if qi in active:
                continue
            fault_sites.append(
                wmc._add_1q_depol_fault(
                    state,
                    qi,
                    trigger_p=gamma * base_p,
                    site_id=f"r{round_idx}_g{gidx}_idle_q{qi}",
                )
            )
    return fault_sites


def _record_meas(
    state: CircuitXZ,
    qc: QuantumCircuit,
    round_idx: int,
    meas_exprs: Dict[str, BoolRef],
    *,
    include_ancX: bool,
    include_ancZ: bool,
    include_flagX: bool,
    include_flagZ: bool,
    beta: float,
    base_p: float,
    fault_sites: List[wmc.FaultSite],
) -> None:
    """Match proof_protocol: ancX/flagX use .z, ancZ/flagZ use .x."""
    idx = _build_reg_index_map(qc)
    tp = beta * base_p

    def _meas_one(qi: int, expr: BoolRef, key: str) -> None:
        fault_sites.append(
            wmc._add_1q_depol_fault(
                state,
                qi,
                trigger_p=tp,
                site_id=f"r{round_idx}_meas_{key}",
            )
        )
        meas_exprs[key] = expr

    for j in range(N_ANC):
        if include_ancX and ("ancX", j) in idx:
            qi = idx[("ancX", j)]
            _meas_one(qi, state.qubits[qi].z, f"r_{round_idx}_ancX{j}")
        if include_ancZ and ("ancZ", j) in idx:
            qi = idx[("ancZ", j)]
            _meas_one(qi, state.qubits[qi].x, f"r_{round_idx}_ancZ{j}")
    for j in range(N_FLAG):
        if include_flagX and ("flagX", j) in idx:
            qi = idx[("flagX", j)]
            _meas_one(qi, state.qubits[qi].z, f"r_{round_idx}_flagX{j}")
        if include_flagZ and ("flagZ", j) in idx:
            qi = idx[("flagZ", j)]
            _meas_one(qi, state.qubits[qi].x, f"r_{round_idx}_flagZ{j}")


def _meas(meas_exprs: Dict[str, BoolRef], key: str) -> BoolRef:
    if key not in meas_exprs:
        raise KeyError(f"missing measurement expr: {key}")
    return meas_exprs[key]


def _flags_all_zero(meas_exprs: Dict[str, BoolRef], round_idx: int) -> BoolRef:
    bits = (
        [_meas(meas_exprs, f"r_{round_idx}_flagX{j}") for j in range(N_FLAG)]
        + [_meas(meas_exprs, f"r_{round_idx}_flagZ{j}") for j in range(N_FLAG)]
    )
    return And(*[Not(b) for b in bits])


def _syndrome_all_zero(meas_exprs: Dict[str, BoolRef], round_idx: int) -> BoolRef:
    bits = (
        [_meas(meas_exprs, f"r_{round_idx}_ancX{j}") for j in range(N_ANC)]
        + [_meas(meas_exprs, f"r_{round_idx}_ancZ{j}") for j in range(N_ANC)]
    )
    return And(*[Not(b) for b in bits])


def _syndrome_equal(meas_exprs: Dict[str, BoolRef], a: int, b: int) -> BoolRef:
    terms: List[BoolRef] = []
    for j in range(N_ANC):
        terms.append(_meas(meas_exprs, f"r_{a}_ancX{j}") == _meas(meas_exprs, f"r_{b}_ancX{j}"))
        terms.append(_meas(meas_exprs, f"r_{a}_ancZ{j}") == _meas(meas_exprs, f"r_{b}_ancZ{j}"))
    return And(*terms)


def _run_flag_round(
    state: CircuitXZ | None,
    flag_qc: QuantumCircuit,
    round_idx: int,
    *,
    p: float,
    gamma: float,
    beta: float,
    fault_sites: List[wmc.FaultSite],
    meas_exprs: Dict[str, BoolRef],
    prev_qc: QuantumCircuit | None = None,
) -> CircuitXZ:
    if state is None:
        state = new_clean_circuit_state(flag_qc.num_qubits)
    else:
        state = _new_state_carry_data(
            state,
            prev_qc or flag_qc,
            flag_qc,
            reset_groups=("ancX", "ancZ", "flagX", "flagZ"),
        )
    fault_sites.extend(
        _sym_run_circuit(
            flag_qc, state, round_idx=round_idx, base_p=p, gamma=gamma, prep_groups=FLAG_PREP
        )
    )
    _record_meas(
        state,
        flag_qc,
        round_idx,
        meas_exprs,
        include_ancX=True,
        include_ancZ=True,
        include_flagX=True,
        include_flagZ=True,
        beta=beta,
        base_p=p,
        fault_sites=fault_sites,
    )
    return state


def _run_raw_round(
    state: CircuitXZ,
    flag_qc: QuantumCircuit,
    raw_qc: QuantumCircuit,
    round_idx: int,
    *,
    p: float,
    gamma: float,
    beta: float,
    fault_sites: List[wmc.FaultSite],
    meas_exprs: Dict[str, BoolRef],
) -> CircuitXZ:
    state = _new_state_carry_data(state, flag_qc, raw_qc, reset_groups=("ancX", "ancZ"))
    fault_sites.extend(
        _sym_run_circuit(
            raw_qc, state, round_idx=round_idx, base_p=p, gamma=gamma, prep_groups=RAW_PREP
        )
    )
    _record_meas(
        state,
        raw_qc,
        round_idx,
        meas_exprs,
        include_ancX=True,
        include_ancZ=True,
        include_flagX=False,
        include_flagZ=False,
        beta=beta,
        base_p=p,
        fault_sites=fault_sites,
    )
    return state


def _path_condition(path_idx: int, meas_exprs: Dict[str, BoolRef]) -> BoolRef:
    """Chris d=3 protocol branches (ter_1 .. ter_4)."""
    if path_idx == 0:
        return Not(_flags_all_zero(meas_exprs, 0))
    if path_idx == 1:
        return And(_flags_all_zero(meas_exprs, 0), Not(_flags_all_zero(meas_exprs, 1)))
    if path_idx == 2:
        return And(
            _flags_all_zero(meas_exprs, 0),
            _flags_all_zero(meas_exprs, 1),
            Not(_syndrome_equal(meas_exprs, 0, 1)),
        )
    if path_idx == 3:
        return And(
            _flags_all_zero(meas_exprs, 0),
            _flags_all_zero(meas_exprs, 1),
            _syndrome_equal(meas_exprs, 0, 1),
        )
    raise ValueError(f"path_idx must be 0..3, got {path_idx}")


def build_path_symbolic(
    path_idx: int,
    *,
    p: float,
    gamma: float,
    beta: float,
) -> PathSymbolic:
    if not (0 <= path_idx <= 3):
        raise ValueError("path_idx must be in 0..3")

    flag_qc = QuantumCircuit.from_qasm_file(str(FLAG_QASM))
    raw_qc = QuantumCircuit.from_qasm_file(str(RAW_QASM))
    fault_sites: List[wmc.FaultSite] = []
    meas_exprs: Dict[str, BoolRef] = {}

    state: CircuitXZ | None = None
    if path_idx == 0:
        state = _run_flag_round(None, flag_qc, 0, p=p, gamma=gamma, beta=beta, fault_sites=fault_sites, meas_exprs=meas_exprs)
        state = _run_raw_round(state, flag_qc, raw_qc, 1, p=p, gamma=gamma, beta=beta, fault_sites=fault_sites, meas_exprs=meas_exprs)
    else:
        state = _run_flag_round(None, flag_qc, 0, p=p, gamma=gamma, beta=beta, fault_sites=fault_sites, meas_exprs=meas_exprs)
        state = _run_flag_round(state, flag_qc, 1, p=p, gamma=gamma, beta=beta, fault_sites=fault_sites, meas_exprs=meas_exprs, prev_qc=flag_qc)
        state = _run_raw_round(state, flag_qc, raw_qc, 2, p=p, gamma=gamma, beta=beta, fault_sites=fault_sites, meas_exprs=meas_exprs)

    raw_groups = detect_qubit_groups(raw_qc)
    data_idxs = sorted(raw_groups["data"])
    data_x = [state.qubits[i].x for i in data_idxs]
    data_z = [state.qubits[i].z for i in data_idxs]

    return PathSymbolic(
        path_idx=path_idx,
        fault_sites=fault_sites,
        meas_exprs=meas_exprs,
        data_x=data_x,
        data_z=data_z,
        path_condition=_path_condition(path_idx, meas_exprs),
    )


def _load_pairs_from_file(path: Path) -> List[Tuple[str, str]]:
    pairs: List[Tuple[str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        xs, zs = line.split()
        pairs.append((xs, zs))
    return pairs


def _parse_dec_formula(dec_name: str, sexpr: str, meas_names: Sequence[str]) -> BoolRef:
    decls = {m: Bool(m) for m in meas_names}
    decls[dec_name] = Bool(dec_name)
    asts = parse_smt2_string(f"(assert (= {dec_name} {sexpr}))", decls=decls)
    if not asts:
        raise ValueError(f"empty parse for {dec_name}")
    return asts[0].arg(1)


def build_failure_formula(path: PathSymbolic) -> Tuple[BoolRef, List[wmc.FaultSite]]:
    dec_file = decoder_path(ORIGIN_DIR, path.path_idx)
    meas_names, dec_formulas = parse_decoder_c_file(dec_file)

    subst_pairs: List[Tuple[BoolRef, BoolRef]] = []
    for m in meas_names:
        if m in path.meas_exprs:
            subst_pairs.append((Bool(m), path.meas_exprs[m]))
        else:
            subst_pairs.append((Bool(m), BoolVal(False)))

    n_data = len(path.data_x)
    fixed_x = list(path.data_x)
    fixed_z = list(path.data_z)
    for j in range(n_data):
        dx_name = f"dec{j}_x"
        dz_name = f"dec{j}_z"
        dx_expr = substitute(
            _parse_dec_formula(dx_name, dec_formulas[dx_name], meas_names), subst_pairs
        )
        dz_expr = substitute(
            _parse_dec_formula(dz_name, dec_formulas[dz_name], meas_names), subst_pairs
        )
        fixed_x[j] = Xor(fixed_x[j], dx_expr)
        fixed_z[j] = Xor(fixed_z[j], dz_expr)

    pairs = _load_pairs_from_file(LOG_TXT) + _load_pairs_from_file(STAB_TXT)
    commute_terms: List[BoolRef] = []
    for xs, zs in pairs:
        if len(xs) != n_data or len(zs) != n_data:
            raise ValueError(f"generator length mismatch: expected {n_data}")
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


@dataclass
class PathResult:
    path_idx: int
    fail_prob: float
    n_vars: int
    n_clauses: int
    elapsed_sec: float
    cnf_path: Path


@dataclass
class WMCResult:
    p: float
    beta: float
    gamma: float
    paths: List[PathResult]
    ratio_fail: float
    total_elapsed_sec: float


def compute_fail_prob_for_path(
    path_idx: int,
    *,
    p: float,
    gamma: float,
    beta: float,
    cnf_path: Path,
    keep_cnf: bool = True,
    max_faults: Optional[int] = None,
) -> PathResult:
    start = time.time()
    sym = build_path_symbolic(path_idx, p=p, gamma=gamma, beta=beta)
    fail_formula, sites = build_failure_formula(sym)
    n_clauses = wmc.emit_dimacs(fail_formula, sites, cnf_path, max_faults=max_faults)
    with cnf_path.open() as f:
        header = f.readline().strip().split()
        n_vars = int(header[2]) if len(header) >= 4 else -1
    prob = wmc.run_gpmc(cnf_path)
    elapsed = time.time() - start
    if not keep_cnf:
        try:
            cnf_path.unlink()
        except FileNotFoundError:
            pass
    return PathResult(
        path_idx=path_idx,
        fail_prob=prob,
        n_vars=n_vars,
        n_clauses=n_clauses,
        elapsed_sec=elapsed,
        cnf_path=cnf_path,
    )


def compute_ratio_fail(
    *,
    p: float,
    gamma: float = 1.0,
    beta: float = 1.0,
    work_dir: Optional[Path] = None,
    keep_cnf: bool = True,
    max_faults: Optional[int] = 2,
) -> WMCResult:
    if work_dir is None:
        work_dir = Path(tempfile.mkdtemp(prefix="wmc_913_"))
    work_dir.mkdir(parents=True, exist_ok=True)
    start = time.time()
    paths: List[PathResult] = []
    for path_idx in range(4):
        cnf_path = work_dir / f"path_{path_idx}_p{p:.8g}.cnf"
        paths.append(
            compute_fail_prob_for_path(
                path_idx,
                p=p,
                gamma=gamma,
                beta=beta,
                cnf_path=cnf_path,
                keep_cnf=keep_cnf,
                max_faults=max_faults,
            )
        )
    return WMCResult(
        p=p,
        beta=beta,
        gamma=gamma,
        paths=paths,
        ratio_fail=sum(pp.fail_prob for pp in paths),
        total_elapsed_sec=time.time() - start,
    )


def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="WMC for [[9,1,3]]_surface_code_seq")
    ap.add_argument("--p", type=float, required=True)
    ap.add_argument("--gamma", type=float, default=1.0)
    ap.add_argument("--beta", type=float, default=1.0)
    ap.add_argument("--work-dir", type=Path, default=None)
    ap.add_argument("--keep-cnf", action="store_true")
    ap.add_argument("--max-faults", type=int, default=2)
    ap.add_argument("--path-idx", type=int, default=None, help="single path 0..3 only")
    return ap.parse_args()


def main() -> None:
    args = _parse_args()
    work_dir = args.work_dir or (ORIGIN_DIR / "wmc_results" / "_cnf_cache")
    work_dir.mkdir(parents=True, exist_ok=True)

    if args.path_idx is not None:
        cnf = work_dir / f"path_{args.path_idx}_p{args.p:.8g}.cnf"
        pr = compute_fail_prob_for_path(
            args.path_idx,
            p=args.p,
            gamma=args.gamma,
            beta=args.beta,
            cnf_path=cnf,
            keep_cnf=args.keep_cnf,
            max_faults=args.max_faults,
        )
        print(f"path_{pr.path_idx}: P(fail)={pr.fail_prob:.8g} vars={pr.n_vars} clauses={pr.n_clauses}")
        return

    res = compute_ratio_fail(
        p=args.p,
        gamma=args.gamma,
        beta=args.beta,
        work_dir=work_dir,
        keep_cnf=args.keep_cnf,
        max_faults=args.max_faults,
    )
    print("[WMC] [[9,1,3]]_surface_code_seq")
    for pr in res.paths:
        print(
            f"  path_{pr.path_idx}: P(fail)={pr.fail_prob:.6g}  "
            f"vars={pr.n_vars}  clauses={pr.n_clauses}  elapsed={pr.elapsed_sec:.2f}s"
        )
    print(f"  ratio_fail (logical error rate) = {res.ratio_fail:.8g}")


if __name__ == "__main__":
    main()
