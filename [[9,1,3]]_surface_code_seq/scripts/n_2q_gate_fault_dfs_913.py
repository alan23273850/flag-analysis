"""
Enumerate [[9,1,3]] surface_code_seq protocol trajectories with **at most n**
injected 2-qubit Pauli faults (after each CX/CZ: no fault or one of 15 non-identity
X/Z flips). Default n=1.

Protocol (Chris d=3 flag):
  - r_0: flag_syndrome
  - if flags != 0 -> path 0 -> raw at r_1
  - else r_1: flag_syndrome
  - if flags != 0 -> path 1 -> raw at r_2
  - elif syndrome r_1 != r_0 -> path 2 -> raw at r_2
  - else path 3 -> raw at r_2

Only 2q gate faults are enumerated (no prep / idle / meas noise), matching
`n_2q_gate_fault_dfs_513.py`.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple
import sys

from qiskit import QuantumCircuit
from z3 import BoolVal, is_true, simplify

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from circuit_op import load_qasm
from decoder_commute_verify import (
    DECODER_BACKENDS,
    decoder_path,
    verify_decoder_commute_913,
)
from flag_analysis import (
    CircuitXZ,
    apply_qasm_gate_into_state,
    detect_qubit_groups,
    load_symplectic_txt,
    new_clean_circuit_state,
    stabilizer_syndrome_from_data,
)

ORIGIN_DIR = Path(__file__).resolve().parents[1]
FLAG_QASM = ORIGIN_DIR / "[[9,1,3]]_flag_syndrome.qasm"
RAW_QASM = ORIGIN_DIR / "[[9,1,3]]_raw_syndrome.qasm"
STAB_TXT = ORIGIN_DIR / "[[9,1,3]].txt"
LOG_TXT = ORIGIN_DIR / "[[9,1,3]]_log_op.txt"

N_ANC = 4
N_FLAG = 2

_FAULT15: List[Tuple[bool, bool, bool, bool]] = [
    (a, b, c, d)
    for a in (False, True)
    for b in (False, True)
    for c in (False, True)
    for d in (False, True)
    if (a or b or c or d)
]

GateCont = Callable[[CircuitXZ, int, List["TwoQubitFaultPlacement"]], None]


@dataclass
class TwoQubitFaultPlacement:
    phase: str  # "flag" | "raw"
    round_idx: int
    gate_index: int
    gate_name: str
    ctrl: int
    targ: int
    flips: Tuple[bool, bool, bool, bool]

    def describe(self, qc: QuantumCircuit) -> str:
        return (
            f"{self.phase} round={self.round_idx} gate[{self.gate_index}]={self.gate_name} "
            f"({self._qubit_pretty(qc, self.ctrl)}->{self._qubit_pretty(qc, self.targ)}) "
            f"flips(cx,cz,tx,tz)={tuple(int(x) for x in self.flips)}"
        )

    @staticmethod
    def _qubit_pretty(qc: QuantumCircuit, idx: int) -> str:
        for qreg in qc.qregs:
            for j, bit in enumerate(qreg):
                if qc.find_bit(bit).index == idx:
                    return f"{qreg.name}[{j}]"
        return f"q[{idx}]"


@dataclass
class DfsStats:
    leaves_lut: int = 0
    leaves_syn_skip: int = 0
    decoder_pass: Dict[str, int] = field(default_factory=lambda: {b: 0 for b in DECODER_BACKENDS})
    decoder_fail: Dict[str, int] = field(default_factory=lambda: {b: 0 for b in DECODER_BACKENDS})


def _to_bit(expr) -> int:
    if isinstance(expr, bool):
        return int(expr)
    return int(is_true(simplify(expr)))


def _build_reg_index_map(qc: QuantumCircuit) -> Dict[Tuple[str, int], int]:
    m: Dict[Tuple[str, int], int] = {}
    g = 0
    for qreg in qc.qregs:
        for j in range(qreg.size):
            m[(qreg.name, j)] = g
            g += 1
    return m


def _flip_qubit_xz(state: CircuitXZ, qidx: int, flip_x: bool, flip_z: bool) -> None:
    x_now = bool(_to_bit(state.qubits[qidx].x))
    z_now = bool(_to_bit(state.qubits[qidx].z))
    if flip_x:
        state.qubits[qidx].x = BoolVal(not x_now)
    if flip_z:
        state.qubits[qidx].z = BoolVal(not z_now)


def _apply_2q_with_optional_fault(
    state: CircuitXZ,
    name: str,
    ctrl: int,
    targ: int,
    fault_choice: Optional[Tuple[bool, bool, bool, bool]],
    fault_count: int,
    max_faults: int,
) -> int:
    apply_qasm_gate_into_state(state, name, [ctrl, targ])
    if fault_choice is None:
        return fault_count
    if fault_count >= max_faults:
        return fault_count
    fc, fz, tc, tz = fault_choice
    _flip_qubit_xz(state, ctrl, fc, fz)
    _flip_qubit_xz(state, targ, tc, tz)
    return fault_count + 1


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


def _reset_qubits(state: CircuitXZ, qidxs: List[int]) -> None:
    for qi in qidxs:
        state.qubits[qi].x = BoolVal(False)
        state.qubits[qi].z = BoolVal(False)


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


def _data_xz_lists(state: CircuitXZ, qc: QuantumCircuit) -> Tuple[List[bool], List[bool]]:
    groups = detect_qubit_groups(qc)
    data_idxs = sorted(groups["data"])
    return (
        [bool(_to_bit(state.qubits[i].x)) for i in data_idxs],
        [bool(_to_bit(state.qubits[i].z)) for i in data_idxs],
    )


def _record_meas_deterministic(
    state: CircuitXZ,
    qc: QuantumCircuit,
    round_idx: int,
    meas: Dict[str, bool],
    *,
    include_ancX: bool,
    include_ancZ: bool,
    include_flagX: bool,
    include_flagZ: bool,
) -> None:
    idx = _build_reg_index_map(qc)

    def _meas_one(qi: int, read_z: bool, key: str) -> None:
        meas[key] = bool(_to_bit(state.qubits[qi].z if read_z else state.qubits[qi].x))

    for j in range(N_ANC):
        if include_ancX and ("ancX", j) in idx:
            _meas_one(idx[("ancX", j)], True, f"r_{round_idx}_ancX{j}")
        if include_ancZ and ("ancZ", j) in idx:
            _meas_one(idx[("ancZ", j)], False, f"r_{round_idx}_ancZ{j}")
    for j in range(N_FLAG):
        if include_flagX and ("flagX", j) in idx:
            _meas_one(idx[("flagX", j)], True, f"r_{round_idx}_flagX{j}")
        if include_flagZ and ("flagZ", j) in idx:
            _meas_one(idx[("flagZ", j)], False, f"r_{round_idx}_flagZ{j}")


def _flags_all_zero(meas: Dict[str, bool], round_idx: int) -> bool:
    for j in range(N_FLAG):
        if meas.get(f"r_{round_idx}_flagX{j}", False):
            return False
        if meas.get(f"r_{round_idx}_flagZ{j}", False):
            return False
    return True


def _syndrome_matches_data(
    meas: Dict[str, bool],
    raw_round_idx: int,
    data_x: List[bool],
    data_z: List[bool],
) -> bool:
    """Measured raw syndrome must match commutation syndrome from data (proof_protocol LUT constraint)."""
    gens = load_symplectic_txt(str(STAB_TXT))
    pred = stabilizer_syndrome_from_data(
        [BoolVal(x) for x in data_x],
        [BoolVal(z) for z in data_z],
        gens,
    )
    pred_b = [bool(is_true(simplify(p))) for p in pred]
    syn_meas = (
        [meas.get(f"r_{raw_round_idx}_ancX{j}", False) for j in range(N_ANC)]
        + [meas.get(f"r_{raw_round_idx}_ancZ{j}", False) for j in range(N_ANC)]
    )
    return syn_meas == pred_b


def _syndrome_equal(meas: Dict[str, bool], a: int, b: int) -> bool:
    for j in range(N_ANC):
        if meas.get(f"r_{a}_ancX{j}") != meas.get(f"r_{b}_ancX{j}"):
            return False
        if meas.get(f"r_{a}_ancZ{j}") != meas.get(f"r_{b}_ancZ{j}"):
            return False
    return True


def _dfs_gate_slice(
    qc: QuantumCircuit,
    state: CircuitXZ,
    lo: int,
    end: int,
    fault_count: int,
    max_faults: int,
    placements: List[TwoQubitFaultPlacement],
    phase: str,
    round_idx: int,
    cont: GateCont,
    stop_flag: Optional[List[bool]] = None,
) -> None:
    if stop_flag and stop_flag[0]:
        return
    if lo >= end:
        cont(state, fault_count, placements)
        return

    instr, qargs, _ = qc.data[lo]
    name = instr.name
    qidxs = [qc.find_bit(q).index for q in qargs]
    lname = name.lower()

    if lname in ("measure", "barrier", "reset", "id"):
        _dfs_gate_slice(
            qc, state, lo + 1, end, fault_count, max_faults, placements,
            phase, round_idx, cont, stop_flag,
        )
        return

    if lname in ("cx", "cnot", "cz"):
        c, t = qidxs[0], qidxs[1]
        choices: List[Optional[Tuple[bool, bool, bool, bool]]]
        if fault_count >= max_faults:
            choices = [None]
        else:
            choices = [None] + _FAULT15

        for fc in choices:
            if stop_flag and stop_flag[0]:
                return
            st = deepcopy(state)
            pls = list(placements)
            if fc is None:
                nfc = _apply_2q_with_optional_fault(st, name, c, t, None, fault_count, max_faults)
                npls = pls
            else:
                nfc = _apply_2q_with_optional_fault(st, name, c, t, fc, fault_count, max_faults)
                npls = pls + [
                    TwoQubitFaultPlacement(
                        phase=phase,
                        round_idx=round_idx,
                        gate_index=lo,
                        gate_name=name,
                        ctrl=c,
                        targ=t,
                        flips=fc,
                    )
                ]
            _dfs_gate_slice(
                qc, st, lo + 1, end, nfc, max_faults, npls,
                phase, round_idx, cont, stop_flag,
            )
        return

    apply_qasm_gate_into_state(state, name, qidxs)
    _dfs_gate_slice(
        qc, state, lo + 1, end, fault_count, max_faults, placements,
        phase, round_idx, cont, stop_flag,
    )


def dfs_protocol(
    *,
    max_2q_faults: int = 1,
    backends: Tuple[str, ...] = DECODER_BACKENDS,
    verbose: bool = False,
    stop_on_first_fail: bool = False,
) -> DfsStats:
    if max_2q_faults < 0:
        raise ValueError("max_2q_faults (n) must be >= 0")

    stats = DfsStats()
    stop_flag: Optional[List[bool]] = [False] if stop_on_first_fail else None

    flag_qc = load_qasm(str(FLAG_QASM))
    raw_qc = load_qasm(str(RAW_QASM))

    def verify_leaf(
        path_idx: int,
        raw_round_idx: int,
        state: CircuitXZ,
        meas: Dict[str, bool],
        placements: List[TwoQubitFaultPlacement],
        fault_count: int,
    ) -> None:
        stats.leaves_lut += 1
        data_x, data_z = _data_xz_lists(state, raw_qc)
        if not _syndrome_matches_data(meas, raw_round_idx, data_x, data_z):
            stats.leaves_syn_skip += 1
            if verbose:
                print(
                    f"SKIP syn_constraint path={path_idx} faults={fault_count} "
                    "(measured raw syndrome != data commutation syndrome)"
                )
            return
        for backend in backends:
            dec_file = decoder_path(ORIGIN_DIR, path_idx, backend=backend)
            ok, _ = verify_decoder_commute_913(
                meas_assign=meas,
                data_x=data_x,
                data_z=data_z,
                decoder_c_path=dec_file,
                log_path=LOG_TXT,
                stab_path=STAB_TXT,
            )
            if ok:
                stats.decoder_pass[backend] += 1
            else:
                stats.decoder_fail[backend] += 1
                if stop_flag is not None:
                    stop_flag[0] = True
                if verbose or stop_on_first_fail:
                    print(
                        f"FAIL backend={backend} path={path_idx} faults={fault_count} "
                        f"placements={[p.describe(flag_qc if p.phase == 'flag' else raw_qc) for p in placements]}"
                    )
            if verbose and ok:
                print(
                    f"PASS backend={backend} path={path_idx} faults={fault_count}"
                )

    def dfs_raw(
        raw_round_idx: int,
        path_idx: int,
        state: CircuitXZ,
        fault_count: int,
        placements: List[TwoQubitFaultPlacement],
        meas: Dict[str, bool],
    ) -> None:
        if stop_flag and stop_flag[0]:
            return

        def after_raw(st: CircuitXZ, fc: int, pls: List[TwoQubitFaultPlacement]) -> None:
            _record_meas_deterministic(
                st,
                raw_qc,
                raw_round_idx,
                meas,
                include_ancX=True,
                include_ancZ=True,
                include_flagX=False,
                include_flagZ=False,
            )
            verify_leaf(path_idx, raw_round_idx, st, meas, pls, fc)

        _dfs_gate_slice(
            raw_qc,
            state,
            0,
            len(raw_qc.data),
            fault_count,
            max_2q_faults,
            placements,
            "raw",
            raw_round_idx,
            after_raw,
            stop_flag,
        )

    def after_flag_round(
        round_idx: int,
        state: CircuitXZ,
        meas: Dict[str, bool],
        fault_count: int,
        placements: List[TwoQubitFaultPlacement],
    ) -> None:
        if stop_flag and stop_flag[0]:
            return

        if round_idx == 0:
            if not _flags_all_zero(meas, 0):
                raw_state = _new_state_carry_data(
                    state, flag_qc, raw_qc, reset_groups=("ancX", "ancZ")
                )
                dfs_raw(1, 0, raw_state, fault_count, placements, meas)
                return
            dfs_flag_round(1, state, fault_count, placements, meas)
            return

        if not _flags_all_zero(meas, 1):
            path_idx = 1
        elif not _syndrome_equal(meas, 0, 1):
            path_idx = 2
        else:
            path_idx = 3
        raw_state = _new_state_carry_data(
            state, flag_qc, raw_qc, reset_groups=("ancX", "ancZ")
        )
        dfs_raw(2, path_idx, raw_state, fault_count, placements, meas)

    def dfs_flag_round(
        round_idx: int,
        prev_state: Optional[CircuitXZ],
        fault_count: int,
        placements: List[TwoQubitFaultPlacement],
        meas: Dict[str, bool],
    ) -> None:
        if stop_flag and stop_flag[0]:
            return

        if prev_state is None:
            state0 = new_clean_circuit_state(flag_qc.num_qubits)
        else:
            state0 = _new_state_carry_data(
                prev_state,
                flag_qc,
                flag_qc,
                reset_groups=("ancX", "ancZ", "flagX", "flagZ"),
            )

        def after_flag(st: CircuitXZ, fc: int, pls: List[TwoQubitFaultPlacement]) -> None:
            _record_meas_deterministic(
                st,
                flag_qc,
                round_idx,
                meas,
                include_ancX=True,
                include_ancZ=True,
                include_flagX=True,
                include_flagZ=True,
            )
            after_flag_round(round_idx, st, meas, fc, pls)

        _dfs_gate_slice(
            flag_qc,
            state0,
            0,
            len(flag_qc.data),
            fault_count,
            max_2q_faults,
            placements,
            "flag",
            round_idx,
            after_flag,
            stop_flag,
        )

    dfs_flag_round(0, None, 0, [], {})
    return stats


def main() -> None:
    ap = argparse.ArgumentParser(
        description="DFS over [[9,1,3]] surface_code_seq with at most N 2q faults."
    )
    ap.add_argument(
        "n",
        nargs="?",
        type=int,
        default=1,
        metavar="N",
        help="Max 2q faults per trajectory (default: 1)",
    )
    ap.add_argument(
        "--backend",
        action="append",
        choices=DECODER_BACKENDS,
        help="Decoder backend to verify (repeatable; default: all three)",
    )
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument(
        "--stop-on-first-fail",
        action="store_true",
        help="Stop at first decoder commute FAIL (any backend).",
    )
    args = ap.parse_args()

    if args.n < 0:
        ap.error("N must be >= 0")

    backends = tuple(args.backend) if args.backend else DECODER_BACKENDS
    stats = dfs_protocol(
        max_2q_faults=args.n,
        backends=backends,
        verbose=args.verbose,
        stop_on_first_fail=args.stop_on_first_fail,
    )

    verified = stats.leaves_lut - stats.leaves_syn_skip
    print(f"max_2q_faults={args.n} | backends={','.join(backends)}")
    print(
        f"leaves LUT: {stats.leaves_lut}  |  syn_constraint skip: {stats.leaves_syn_skip}  |  "
        f"verified: {verified}"
    )
    for b in backends:
        p = stats.decoder_pass[b]
        f = stats.decoder_fail[b]
        print(f"  {b}: PASS {p}/{verified}" + ("" if f == 0 else f"  FAIL {f}"))


if __name__ == "__main__":
    main()
