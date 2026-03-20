from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple
import random

from z3 import BoolVal, is_true, simplify

from qiskit import QuantumCircuit

from flag_analysis import (
    CircuitXZ,
    apply_qasm_gate_into_state,
    detect_qubit_groups,
    new_clean_circuit_state,
)


FLAG_QASM = Path("./[[5,1,3]]_origin/[[5,1,3]]_origin_flag_syndrome.qasm ")
RAW_QASM = Path("./[[5,1,3]]_origin/raw_[[5,1,3]]_origin_flag_syndrome.qasm")
TWO_QUBIT_FAULT_P = 5e-2


@dataclass
class LutRecord:
    first_stabilizer_index: int
    first_s: int
    first_f: int
    second_s_bits: List[int]
    bitstring6: str


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


def _run_gate_slice(qc: QuantumCircuit, state: CircuitXZ, start: int, end: int) -> None:
    for instr, qargs, _ in qc.data[start:end]:
        name = instr.name
        qidxs = [qc.find_bit(q).index for q in qargs]
        apply_qasm_gate_into_state(state, name, qidxs)
        if name in ("cx", "cnot", "cz"):
            _inject_two_qubit_stochastic_fault(state, qidxs[0], qidxs[1], p=TWO_QUBIT_FAULT_P)


def _flip_qubit_xz(state: CircuitXZ, qidx: int, flip_x: bool, flip_z: bool) -> None:
    x_now = bool(_to_bit(state.qubits[qidx].x))
    z_now = bool(_to_bit(state.qubits[qidx].z))
    if flip_x:
        state.qubits[qidx].x = BoolVal(not x_now)
    if flip_z:
        state.qubits[qidx].z = BoolVal(not z_now)


def _inject_two_qubit_stochastic_fault(state: CircuitXZ, ctrl: int, targ: int, p: float = 1e-4) -> None:
    """
    Two-qubit stochastic Pauli-like fault after each CX/CZ gate.

    - With probability (1-p): no fault (I,I)
    - With probability p: pick uniformly from 15 non-identity events
      over control/target x/z flips.
    """
    if random.random() >= p:
        return

    # Each tuple: (flip_cx, flip_cz, flip_tx, flip_tz), excluding all-False.
    events = [(a, b, c, d) for a in (False, True)
                           for b in (False, True)
                           for c in (False, True)
                           for d in (False, True)
                           if (a or b or c or d)]
    flip_cx, flip_cz, flip_tx, flip_tz = random.choice(events)
    _flip_qubit_xz(state, ctrl, flip_cx, flip_cz)
    _flip_qubit_xz(state, targ, flip_tx, flip_tz)


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


def _run_second_subround(first_state: CircuitXZ, first_qc: QuantumCircuit, raw_qc: QuantumCircuit) -> List[int]:
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
        # Explicitly initialize each second-subround ancX[i] as a fresh ancilla.
        second_state.qubits[raw_anc[i]].x = BoolVal(False)
        second_state.qubits[raw_anc[i]].z = BoolVal(False)

        start = i * gates_per_stab
        end = (i + 1) * gates_per_stab
        _run_gate_slice(raw_qc, second_state, start, end)
        s_bits.append(_to_bit(second_state.qubits[raw_anc[i]].x))
    return s_bits


def run_protocol_once() -> List[LutRecord]:
    flag_qc = QuantumCircuit.from_qasm_file(str(FLAG_QASM))
    raw_qc = QuantumCircuit.from_qasm_file(str(RAW_QASM))

    first_state = new_clean_circuit_state(flag_qc.num_qubits)
    idx = _build_reg_index_map(flag_qc)
    anc_s_idx = idx[("ancX", 0)]
    flag_f_idx = idx[("flagZ", 0)]

    n_gates = len(flag_qc.data)
    if n_gates % 4 != 0:
        raise ValueError(f"Flag circuit gate count {n_gates} cannot be split into 4 stabilizers.")
    gates_per_stab = n_gates // 4

    lut_records: List[LutRecord] = []

    for stab_i in range(4):
        start = stab_i * gates_per_stab
        end = (stab_i + 1) * gates_per_stab
        _run_gate_slice(flag_qc, first_state, start, end)

        s = _to_bit(first_state.qubits[anc_s_idx].x)
        f = _to_bit(first_state.qubits[flag_f_idx].z)

        # Enter second subround only when [s,f] != [0,0].
        if not (s == 0 and f == 0):
            second_s = _run_second_subround(first_state, flag_qc, raw_qc)
            bitstring6 = f"{s}{f}{''.join(str(b) for b in second_s)}"
            lut_records.append(
                LutRecord(
                    first_stabilizer_index=stab_i,
                    first_s=s,
                    first_f=f,
                    second_s_bits=second_s,
                    bitstring6=bitstring6,
                )
            )
            # Once second subround is triggered, protocol terminates at LUT.
            return lut_records

    # If all first-subround checks are [0,0], there is no LUT output.
    return lut_records


def main() -> None:
    records = run_protocol_once()
    if not records:
        print("End of protocol with no LUT output.")
        return

    print("LUT records:")
    for rec in records:
        print(
            f"first_stab={rec.first_stabilizer_index} "
            f"(s,f)=({rec.first_s},{rec.first_f}) "
            f"second_s={rec.second_s_bits} "
            f"6bit={rec.bitstring6}"
        )


if __name__ == "__main__":
    main()
