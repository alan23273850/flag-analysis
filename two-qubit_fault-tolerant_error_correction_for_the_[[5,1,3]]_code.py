from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import random

from z3 import BoolVal, is_true, simplify

from qiskit import QuantumCircuit

from decoder_commute_verify import decoder_file_index, verify_decoder_commute
from flag_analysis import (
    CircuitXZ,
    apply_qasm_gate_into_state,
    detect_qubit_groups,
    new_clean_circuit_state,
)


FLAG_QASM = Path("./[[5,1,3]]_origin/[[5,1,3]]_origin_flag_syndrome.qasm ")
RAW_QASM = Path("./[[5,1,3]]_origin/raw_[[5,1,3]]_origin_flag_syndrome.qasm")
ORIGIN_DIR = Path("./[[5,1,3]]_origin")
STAB_TXT = ORIGIN_DIR / "[[5,1,3]]_origin.txt"
LOG_TXT = ORIGIN_DIR / "[[5,1,3]]_log_op.txt"
TWO_QUBIT_FAULT_P = 1e-2
IDLE_GAMMA = 1e-2        # 0 <= gamma <= 1
MEAS_BETA = 1.0         # 1 <= beta <= 10


@dataclass
class LutRecord:
    first_stabilizer_index: int
    first_s: int
    first_f: int
    second_s_bits: List[int]
    bitstring6: str


@dataclass
class FaultEvent:
    kind: str
    message: str


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


def _qidx_to_name(qc: QuantumCircuit, qidx: int) -> str:
    for qreg in qc.qregs:
        for j, bit in enumerate(qreg):
            if qc.find_bit(bit).index == qidx:
                return f"{qreg.name}[{j}]"
    return f"q[{qidx}]"


def _apply_pauli_on_qubit(state: CircuitXZ, qidx: int, pauli: str) -> None:
    """Apply single-qubit Pauli to tracked (x,z) error bits."""
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
    # Uniform over X/Y/Z
    return random.choice(("X", "Y", "Z"))


def _append_event(fault_events: List[FaultEvent], kind: str, message: str) -> None:
    fault_events.append(FaultEvent(kind=kind, message=message))


def _maybe_inject_prep_error(
    *,
    state: CircuitXZ,
    qidx: int,
    axis: str,  # "x" or "z"
    base_p: float,
    phase: str,
    stab_or_seg: int,
    label: str,
    fault_events: List[FaultEvent],
) -> None:
    """
    Ancilla preparation error:
      trigger prob = 2p/3
      - axis 'x' => flip X
      - axis 'z' => flip Z
    """
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
    _append_event(
        fault_events,
        "prep",
        f"phase={phase}, stabilizer={stab_or_seg}, qubit={label}, flips={flips}",
    )


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
    """Inject 1q depolarizing error with trigger_p, then X/Y/Z uniformly."""
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
    _append_event(
        fault_events,
        reason,
        f"phase={phase}, qubit={label}, flips={flips}",
    )


def _run_gate_slice(
    qc: QuantumCircuit,
    state: CircuitXZ,
    start: int,
    end: int,
    *,
    base_p: float,
    gamma: float,
    phase: str,
    stab_or_seg: int,
    fault_events: List[FaultEvent],
) -> None:
    # Idle noise only applies to qubits that actually participate in this slice.
    participating_qubits: set[int] = set()
    for instr, qargs, _ in qc.data[start:end]:
        _ = instr
        participating_qubits.update(qc.find_bit(q).index for q in qargs)

    for gidx in range(start, end):
        instr, qargs, _ = qc.data[gidx]
        name = instr.name
        qidxs = [qc.find_bit(q).index for q in qargs]
        apply_qasm_gate_into_state(state, name, qidxs)
        if name in ("cx", "cnot", "cz"):
            flips = _inject_two_qubit_stochastic_fault(state, qidxs[0], qidxs[1], p=base_p)
            if flips is not None:
                fc, fz, tx, tz = flips
                ctrl_name = _qidx_to_name(qc, qidxs[0])
                targ_name = _qidx_to_name(qc, qidxs[1])
                _append_event(
                    fault_events,
                    "2q_gate",
                    f"phase={phase}, stabilizer={stab_or_seg}, gate_idx={gidx}, gate={name}, "
                    f"qubits={ctrl_name},{targ_name}, at=after_gate, flips=({int(fc)},{int(fz)})({int(tx)},{int(tz)})",
                )

        # Idle errors: one time step per processed gate, on qubits not touched by this gate.
        active = set(qidxs)
        for qi in participating_qubits:
            if qi in active:
                continue
            _maybe_inject_single_qubit_depolarizing(
                state=state,
                qidx=qi,
                trigger_p=gamma * base_p,
                phase=f"{phase}, gate_idx={gidx}, at=after_gate",
                reason="idle",
                label=_qidx_to_name(qc, qi),
                fault_events=fault_events,
            )


def _flip_qubit_xz(state: CircuitXZ, qidx: int, flip_x: bool, flip_z: bool) -> None:
    x_now = bool(_to_bit(state.qubits[qidx].x))
    z_now = bool(_to_bit(state.qubits[qidx].z))
    if flip_x:
        state.qubits[qidx].x = BoolVal(not x_now)
    if flip_z:
        state.qubits[qidx].z = BoolVal(not z_now)


def _inject_two_qubit_stochastic_fault(
    state: CircuitXZ, ctrl: int, targ: int, p: float = 1e-4
) -> Optional[Tuple[bool, bool, bool, bool]]:
    """
    Two-qubit stochastic Pauli-like fault after each CX/CZ gate.

    - With probability (1-p): no fault (I,I)
    - With probability p: pick uniformly from 15 non-identity events
      over control/target x/z flips.
    """
    if random.random() >= p:
        return None

    # Each tuple: (flip_cx, flip_cz, flip_tx, flip_tz), excluding all-False.
    events = [(a, b, c, d) for a in (False, True)
                           for b in (False, True)
                           for c in (False, True)
                           for d in (False, True)
                           if (a or b or c or d)]
    flip_cx, flip_cz, flip_tx, flip_tz = random.choice(events)
    _flip_qubit_xz(state, ctrl, flip_cx, flip_cz)
    _flip_qubit_xz(state, targ, flip_tx, flip_tz)
    return (flip_cx, flip_cz, flip_tx, flip_tz)


def _copy_data_error_only(src_state: CircuitXZ, src_qc: QuantumCircuit, dst_qc: QuantumCircuit) -> CircuitXZ:
    src_groups = detect_qubit_groups(src_qc)
    dst_groups = detect_qubit_groups(dst_qc)
    src_data = src_groups["data"]
    dst_data = dst_groups["data"]

    if len(src_data) != len(dst_data):
        raise ValueError("Data qubit count mismatch between first and second subround circuits.")

    dst_state = new_clean_circuit_state(dst_qc.num_qubits)
    for s_idx, d_idx in zip(src_data, dst_data):
        dst_state.qubits[d_idx].x = src_state.qubits[s_idx].x
        dst_state.qubits[d_idx].z = src_state.qubits[s_idx].z
    return dst_state


def _run_second_subround(
    first_state: CircuitXZ,
    first_qc: QuantumCircuit,
    raw_qc: QuantumCircuit,
    fault_events: List[FaultEvent],
    *,
    base_p: float,
    gamma: float,
    beta: float,
) -> List[int]:
    # Second subround ancillas are treated as fresh ancillas.
    # Data errors must carry over from first subround.
    second_state = _copy_data_error_only(first_state, first_qc, raw_qc)

    raw_idx = _build_reg_index_map(raw_qc)
    raw_anc = [raw_idx[("ancX", i)] for i in range(4)]

    n_gates = len(raw_qc.data)
    if n_gates % 4 != 0:
        raise ValueError(f"Raw circuit gate count {n_gates} cannot be split into 4 stabilizers.")
    gates_per_stab = n_gates // 4

    s_bits: List[int] = []
    for i in range(4):
        start = i * gates_per_stab
        end = (i + 1) * gates_per_stab
        # Explicitly initialize each second-subround ancX[i] as a fresh ancilla.
        second_state.qubits[raw_anc[i]].x = BoolVal(False)
        second_state.qubits[raw_anc[i]].z = BoolVal(False)
        # Prep error on |+>/|-> ancX => flip Z with prob 2p/3
        _maybe_inject_prep_error(
            state=second_state,
            qidx=raw_anc[i],
            axis="z",
            base_p=base_p,
        phase=f"second_subround at=before_gate_idx={start}",
            stab_or_seg=i,
            label=f"ancX[{i}]",
            fault_events=fault_events,
        )
        _run_gate_slice(
            raw_qc,
            second_state,
            start,
            end,
            base_p=base_p,
            gamma=gamma,
            phase="second_subround",
            stab_or_seg=i,
            fault_events=fault_events,
        )
        # Measurement error before reading ancX[i] in Z-basis
        _maybe_inject_single_qubit_depolarizing(
            state=second_state,
            qidx=raw_anc[i],
            trigger_p=beta * base_p,
            phase=f"second_subround seg={i} at=before_readout_after_gate_idx={end - 1}",
            reason="meas",
            label=f"ancX[{i}] (Z-basis readout)",
            fault_events=fault_events,
        )
        s_bits.append(_to_bit(second_state.qubits[raw_anc[i]].z))
    return s_bits


def _data_xz_lists(state: CircuitXZ, qc: QuantumCircuit) -> Tuple[List[bool], List[bool]]:
    groups = detect_qubit_groups(qc)
    data_idxs = sorted(groups["data"])
    return (
        [bool(_to_bit(state.qubits[i].x)) for i in data_idxs],
        [bool(_to_bit(state.qubits[i].z)) for i in data_idxs],
    )


def run_protocol_once() -> Tuple[Optional[LutRecord], Optional[List[bool]], Optional[List[bool]], List[FaultEvent]]:
    if not (0.0 <= IDLE_GAMMA <= 1.0):
        raise ValueError(f"IDLE_GAMMA must be in [0,1], got {IDLE_GAMMA}")
    if not (1.0 <= MEAS_BETA <= 10.0):
        raise ValueError(f"MEAS_BETA must be in [1,10], got {MEAS_BETA}")

    flag_qc = QuantumCircuit.from_qasm_file(str(FLAG_QASM))
    raw_qc = QuantumCircuit.from_qasm_file(str(RAW_QASM))

    first_state = new_clean_circuit_state(flag_qc.num_qubits)
    fault_events: List[FaultEvent] = []
    idx = _build_reg_index_map(flag_qc)
    anc_s_idx = idx[("ancX", 0)]
    flag_f_idx = idx[("flagZ", 0)]
    # First subround initial ancilla prep errors:
    # ancX[0] is |+>/|-> => prep error toggles Z with prob 2p/3
    _maybe_inject_prep_error(
        state=first_state,
        qidx=anc_s_idx,
        axis="z",
        base_p=TWO_QUBIT_FAULT_P,
        phase="first_subround at=before_gate_idx=0",
        stab_or_seg=0,
        label="ancX[0]",
        fault_events=fault_events,
    )
    # flagZ[0] is |0>/|1> => prep error toggles X with prob 2p/3
    _maybe_inject_prep_error(
        state=first_state,
        qidx=flag_f_idx,
        axis="x",
        base_p=TWO_QUBIT_FAULT_P,
        phase="first_subround at=before_gate_idx=0",
        stab_or_seg=0,
        label="flagZ[0]",
        fault_events=fault_events,
    )

    n_gates = len(flag_qc.data)
    if n_gates % 4 != 0:
        raise ValueError(f"Flag circuit gate count {n_gates} cannot be split into 4 stabilizers.")
    gates_per_stab = n_gates // 4

    for stab_i in range(4):
        start = stab_i * gates_per_stab
        end = (stab_i + 1) * gates_per_stab
        _run_gate_slice(
            flag_qc,
            first_state,
            start,
            end,
            base_p=TWO_QUBIT_FAULT_P,
            gamma=IDLE_GAMMA,
            phase="first_subround",
            stab_or_seg=stab_i,
            fault_events=fault_events,
        )

        # Measurement errors before first-subround readout
        _maybe_inject_single_qubit_depolarizing(
            state=first_state,
            qidx=anc_s_idx,
            trigger_p=MEAS_BETA * TWO_QUBIT_FAULT_P,
            phase=f"first_subround stab/seg={stab_i} at=before_readout",
            reason="meas",
            label="ancX[0] (Z-basis readout)",
            fault_events=fault_events,
        )
        _maybe_inject_single_qubit_depolarizing(
            state=first_state,
            qidx=flag_f_idx,
            trigger_p=MEAS_BETA * TWO_QUBIT_FAULT_P,
            phase=f"first_subround stab/seg={stab_i} at=before_readout",
            reason="meas",
            label="flagZ[0] (X-basis readout)",
            fault_events=fault_events,
        )
        s = _to_bit(first_state.qubits[anc_s_idx].z)
        f = _to_bit(first_state.qubits[flag_f_idx].x)

        # Enter second subround only when [s,f] != [0,0].
        if not (s == 0 and f == 0):
            second_s = _run_second_subround(
                first_state,
                flag_qc,
                raw_qc,
                fault_events,
                base_p=TWO_QUBIT_FAULT_P,
                gamma=IDLE_GAMMA,
                beta=MEAS_BETA,
            )
            bitstring6 = f"{s}{f}{''.join(str(b) for b in second_s)}"
            # Once second subround is triggered, protocol terminates at LUT.
            dx, dz = _data_xz_lists(first_state, flag_qc)
            return (
                LutRecord(
                    first_stabilizer_index=stab_i,
                    first_s=s,
                    first_f=f,
                    second_s_bits=second_s,
                    bitstring6=bitstring6,
                ),
                dx,
                dz,
                fault_events,
            )

    # If all first-subround checks are [0,0], there is no LUT output.
    return None, None, None, fault_events


def main() -> None:
    rec, data_x, data_z, fault_events = run_protocol_once()
    print(f"Noise params: p={TWO_QUBIT_FAULT_P}, gamma={IDLE_GAMMA}, beta={MEAS_BETA}")
    print("Injected error events:")
    if not fault_events:
        print("  (none)")
    else:
        for i, ev in enumerate(fault_events):
            print(f"  [{i}] kind={ev.kind}, {ev.message}")
    if rec is None:
        print("End of protocol with no LUT output.")
        return

    print("LUT record:")
    print(
        f"first_stab={rec.first_stabilizer_index} "
        f"(s,f)=({rec.first_s},{rec.first_f}) "
        f"second_s={rec.second_s_bits} "
        f"6bit={rec.bitstring6}"
    )

    dec_idx = decoder_file_index(rec.first_stabilizer_index)
    dec_path = Path(f"./decoder_C_{dec_idx}.txt")
    ok, _dec_vals = verify_decoder_commute(
        first_stabilizer_index=rec.first_stabilizer_index,
        bitstring6=rec.bitstring6,
        data_x=data_x or [],
        data_z=data_z or [],
        decoder_c_path=dec_path,
        log_path=LOG_TXT,
        stab_path=STAB_TXT,
    )
    print(f"decoder_C_{dec_idx}.txt commute check (fixed vs log+stab): {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
