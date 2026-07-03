"""Stochastic Monte Carlo simulator for [[9,1,3]] Chris d=3 flag protocol."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
import random
import sys

from qiskit import QuantumCircuit
from z3 import BoolVal, is_true, simplify

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from decoder_commute_verify import decoder_path, verify_decoder_commute_913
from flag_analysis import (
    CircuitXZ,
    apply_qasm_gate_into_state,
    detect_qubit_groups,
    new_clean_circuit_state,
)

ORIGIN_DIR = Path(__file__).resolve().parents[1]
FLAG_QASM = ORIGIN_DIR / "[[9,1,3]]_flag_syndrome.qasm"
RAW_QASM = ORIGIN_DIR / "[[9,1,3]]_raw_syndrome.qasm"
STAB_TXT = ORIGIN_DIR / "[[9,1,3]].txt"
LOG_TXT = ORIGIN_DIR / "[[9,1,3]]_log_op.txt"

N_ANC = 4
N_FLAG = 2
FLAG_PREP = (("ancX", "z"), ("ancZ", "x"), ("flagX", "x"), ("flagZ", "z"))
RAW_PREP = (("ancX", "z"), ("ancZ", "x"))

TWO_QUBIT_FAULT_P = 1e-4
IDLE_GAMMA = 1.0
MEAS_BETA = 1.0

_FLAG_QC_CACHE: Optional[QuantumCircuit] = None
_RAW_QC_CACHE: Optional[QuantumCircuit] = None


@dataclass
class ProtocolRecord:
    path_idx: int
    meas_assign: Dict[str, bool]


@dataclass
class FaultEvent:
    kind: str
    message: str


def preload_static_inputs() -> None:
    global _FLAG_QC_CACHE, _RAW_QC_CACHE
    if _FLAG_QC_CACHE is None:
        _FLAG_QC_CACHE = QuantumCircuit.from_qasm_file(str(FLAG_QASM))
    if _RAW_QC_CACHE is None:
        _RAW_QC_CACHE = QuantumCircuit.from_qasm_file(str(RAW_QASM))


def _to_bit(expr) -> int:
    if isinstance(expr, bool):
        return int(expr)
    return int(is_true(simplify(expr)))


def _build_reg_index_map(qc: QuantumCircuit) -> Dict[Tuple[str, int], int]:
    idx_map: Dict[Tuple[str, int], int] = {}
    global_idx = 0
    for qreg in qc.qregs:
        for local_idx in range(qreg.size):
            idx_map[(qreg.name, local_idx)] = global_idx
            global_idx += 1
    return idx_map


def _flip_qubit_xz(state: CircuitXZ, qidx: int, flip_x: bool, flip_z: bool) -> None:
    x_now = bool(_to_bit(state.qubits[qidx].x))
    z_now = bool(_to_bit(state.qubits[qidx].z))
    if flip_x:
        state.qubits[qidx].x = BoolVal(not x_now)
    if flip_z:
        state.qubits[qidx].z = BoolVal(not z_now)


def _apply_pauli_on_qubit(state: CircuitXZ, qidx: int, pauli: str) -> None:
    p = pauli.upper()
    if p == "X":
        _flip_qubit_xz(state, qidx, True, False)
    elif p == "Z":
        _flip_qubit_xz(state, qidx, False, True)
    elif p == "Y":
        _flip_qubit_xz(state, qidx, True, True)
    else:
        raise ValueError(f"Unsupported Pauli: {pauli}")


def _sample_single_qubit_depolarizing_pauli() -> str:
    return random.choice(("X", "Y", "Z"))


def _append_event(fault_events: List[FaultEvent], kind: str, message: str) -> None:
    fault_events.append(FaultEvent(kind=kind, message=message))


def _maybe_inject_prep_error(
    *,
    state: CircuitXZ,
    qidx: int,
    axis: str,
    base_p: float,
    phase: str,
    label: str,
    fault_events: List[FaultEvent],
) -> None:
    trigger_p = min(1.0, (2.0 * base_p) / 3.0)
    if random.random() >= trigger_p:
        return
    if axis == "x":
        _flip_qubit_xz(state, qidx, True, False)
        flips = "(1,0)"
    elif axis == "z":
        _flip_qubit_xz(state, qidx, False, True)
        flips = "(0,1)"
    else:
        raise ValueError("axis must be 'x' or 'z'")
    _append_event(fault_events, "prep", f"phase={phase}, qubit={label}, flips={flips}")


def _maybe_inject_single_qubit_depolarizing(
    *,
    state: CircuitXZ,
    qidx: int,
    trigger_p: float,
    phase: str,
    reason: str,
    label: str,
    fault_events: List[FaultEvent],
) -> None:
    p = min(1.0, trigger_p)
    if random.random() >= p:
        return
    pauli = _sample_single_qubit_depolarizing_pauli()
    _apply_pauli_on_qubit(state, qidx, pauli)
    if pauli == "X":
        flips = "(1,0)"
    elif pauli == "Z":
        flips = "(0,1)"
    else:
        flips = "(1,1)"
    _append_event(fault_events, reason, f"phase={phase}, qubit={label}, flips={flips}")


def _inject_two_qubit_stochastic_fault(
    state: CircuitXZ, ctrl: int, targ: int, p: float
) -> Optional[Tuple[bool, bool, bool, bool]]:
    if random.random() >= p:
        return None
    events = [
        (a, b, c, d)
        for a in (False, True)
        for b in (False, True)
        for c in (False, True)
        for d in (False, True)
        if (a or b or c or d)
    ]
    flip_cx, flip_cz, flip_tx, flip_tz = random.choice(events)
    _flip_qubit_xz(state, ctrl, flip_cx, flip_cz)
    _flip_qubit_xz(state, targ, flip_tx, flip_tz)
    return (flip_cx, flip_cz, flip_tx, flip_tz)


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


def _run_circuit_stochastic(
    qc: QuantumCircuit,
    state: CircuitXZ,
    *,
    round_idx: int,
    base_p: float,
    gamma: float,
    prep_groups: Sequence[Tuple[str, str]],
    fault_events: List[FaultEvent],
) -> None:
    groups = detect_qubit_groups(qc)
    for reg_name, axis in prep_groups:
        for j in groups.get(reg_name, []):
            _maybe_inject_prep_error(
                state=state,
                qidx=j,
                axis=axis,
                base_p=base_p,
                phase=f"r{round_idx}_prep_{reg_name}{j}",
                label=f"{reg_name}[{j}]",
                fault_events=fault_events,
            )

    participating: set[int] = set()
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
            _inject_two_qubit_stochastic_fault(state, qidxs[0], qidxs[1], p=base_p)
        active = set(qidxs)
        for qi in participating:
            if qi in active:
                continue
            _maybe_inject_single_qubit_depolarizing(
                state=state,
                qidx=qi,
                trigger_p=gamma * base_p,
                phase=f"r{round_idx}_g{gidx}_idle",
                reason="idle",
                label=f"q{qi}",
                fault_events=fault_events,
            )


def _record_meas_stochastic(
    state: CircuitXZ,
    qc: QuantumCircuit,
    round_idx: int,
    meas: Dict[str, bool],
    *,
    include_ancX: bool,
    include_ancZ: bool,
    include_flagX: bool,
    include_flagZ: bool,
    beta: float,
    base_p: float,
    fault_events: List[FaultEvent],
) -> None:
    idx = _build_reg_index_map(qc)
    tp = beta * base_p

    def _meas_one(qi: int, read_z: bool, key: str) -> None:
        _maybe_inject_single_qubit_depolarizing(
            state=state,
            qidx=qi,
            trigger_p=tp,
            phase=f"r{round_idx}_meas_{key}",
            reason="meas",
            label=key,
            fault_events=fault_events,
        )
        meas[key] = bool(_to_bit(state.qubits[qi].z if read_z else state.qubits[qi].x))

    for j in range(N_ANC):
        if include_ancX and ("ancX", j) in idx:
            qi = idx[("ancX", j)]
            _meas_one(qi, True, f"r_{round_idx}_ancX{j}")
        if include_ancZ and ("ancZ", j) in idx:
            qi = idx[("ancZ", j)]
            _meas_one(qi, False, f"r_{round_idx}_ancZ{j}")
    for j in range(N_FLAG):
        if include_flagX and ("flagX", j) in idx:
            qi = idx[("flagX", j)]
            _meas_one(qi, True, f"r_{round_idx}_flagX{j}")
        if include_flagZ and ("flagZ", j) in idx:
            qi = idx[("flagZ", j)]
            _meas_one(qi, False, f"r_{round_idx}_flagZ{j}")


def _flags_all_zero(meas: Dict[str, bool], round_idx: int) -> bool:
    for j in range(N_FLAG):
        if meas.get(f"r_{round_idx}_flagX{j}", False):
            return False
        if meas.get(f"r_{round_idx}_flagZ{j}", False):
            return False
    return True


def _syndrome_equal(meas: Dict[str, bool], a: int, b: int) -> bool:
    for j in range(N_ANC):
        if meas.get(f"r_{a}_ancX{j}") != meas.get(f"r_{b}_ancX{j}"):
            return False
        if meas.get(f"r_{a}_ancZ{j}") != meas.get(f"r_{b}_ancZ{j}"):
            return False
    return True


def _run_flag_round(
    state: CircuitXZ | None,
    flag_qc: QuantumCircuit,
    round_idx: int,
    *,
    base_p: float,
    gamma: float,
    beta: float,
    meas: Dict[str, bool],
    fault_events: List[FaultEvent],
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
    _run_circuit_stochastic(
        flag_qc,
        state,
        round_idx=round_idx,
        base_p=base_p,
        gamma=gamma,
        prep_groups=FLAG_PREP,
        fault_events=fault_events,
    )
    _record_meas_stochastic(
        state,
        flag_qc,
        round_idx,
        meas,
        include_ancX=True,
        include_ancZ=True,
        include_flagX=True,
        include_flagZ=True,
        beta=beta,
        base_p=base_p,
        fault_events=fault_events,
    )
    return state


def _run_raw_round(
    state: CircuitXZ,
    flag_qc: QuantumCircuit,
    raw_qc: QuantumCircuit,
    round_idx: int,
    *,
    base_p: float,
    gamma: float,
    beta: float,
    meas: Dict[str, bool],
    fault_events: List[FaultEvent],
) -> CircuitXZ:
    state = _new_state_carry_data(state, flag_qc, raw_qc, reset_groups=("ancX", "ancZ"))
    _run_circuit_stochastic(
        raw_qc,
        state,
        round_idx=round_idx,
        base_p=base_p,
        gamma=gamma,
        prep_groups=RAW_PREP,
        fault_events=fault_events,
    )
    _record_meas_stochastic(
        state,
        raw_qc,
        round_idx,
        meas,
        include_ancX=True,
        include_ancZ=True,
        include_flagX=False,
        include_flagZ=False,
        beta=beta,
        base_p=base_p,
        fault_events=fault_events,
    )
    return state


def _data_xz_lists(state: CircuitXZ, qc: QuantumCircuit) -> Tuple[List[bool], List[bool]]:
    groups = detect_qubit_groups(qc)
    data_idxs = sorted(groups["data"])
    return (
        [bool(_to_bit(state.qubits[i].x)) for i in data_idxs],
        [bool(_to_bit(state.qubits[i].z)) for i in data_idxs],
    )


def run_protocol_once() -> Tuple[ProtocolRecord, List[bool], List[bool], List[FaultEvent]]:
    if not (0.0 <= IDLE_GAMMA <= 1.0):
        raise ValueError(f"IDLE_GAMMA must be in [0,1], got {IDLE_GAMMA}")
    if not (1.0 <= MEAS_BETA <= 10.0):
        raise ValueError(f"MEAS_BETA must be in [1,10], got {MEAS_BETA}")

    preload_static_inputs()
    if _FLAG_QC_CACHE is None or _RAW_QC_CACHE is None:
        raise RuntimeError("QASM preload failed")
    flag_qc = _FLAG_QC_CACHE
    raw_qc = _RAW_QC_CACHE

    base_p = TWO_QUBIT_FAULT_P
    gamma = IDLE_GAMMA
    beta = MEAS_BETA
    fault_events: List[FaultEvent] = []
    meas: Dict[str, bool] = {}

    state = _run_flag_round(
        None, flag_qc, 0, base_p=base_p, gamma=gamma, beta=beta, meas=meas, fault_events=fault_events
    )

    if not _flags_all_zero(meas, 0):
        path_idx = 0
        state = _run_raw_round(
            state, flag_qc, raw_qc, 1, base_p=base_p, gamma=gamma, beta=beta, meas=meas, fault_events=fault_events
        )
    else:
        state = _run_flag_round(
            state,
            flag_qc,
            1,
            base_p=base_p,
            gamma=gamma,
            beta=beta,
            meas=meas,
            fault_events=fault_events,
            prev_qc=flag_qc,
        )
        if not _flags_all_zero(meas, 1):
            path_idx = 1
        elif not _syndrome_equal(meas, 0, 1):
            path_idx = 2
        else:
            path_idx = 3
        state = _run_raw_round(
            state, flag_qc, raw_qc, 2, base_p=base_p, gamma=gamma, beta=beta, meas=meas, fault_events=fault_events
        )

    data_x, data_z = _data_xz_lists(state, raw_qc)
    return ProtocolRecord(path_idx=path_idx, meas_assign=meas), data_x, data_z, fault_events


def main() -> None:
    rec, data_x, data_z, fault_events = run_protocol_once()
    print(f"Noise params: p={TWO_QUBIT_FAULT_P}, gamma={IDLE_GAMMA}, beta={MEAS_BETA}")
    print(f"path_idx={rec.path_idx}")
    dec_file = decoder_path(ORIGIN_DIR, rec.path_idx)
    ok, _ = verify_decoder_commute_913(
        meas_assign=rec.meas_assign,
        data_x=data_x,
        data_z=data_z,
        decoder_c_path=dec_file,
        log_path=LOG_TXT,
        stab_path=STAB_TXT,
    )
    print(f"path_{rec.path_idx}.txt commute check: {'PASS' if ok else 'FAIL'}")
    print(f"fault events: {len(fault_events)}")


if __name__ == "__main__":
    main()
