"""
Enumerate [[5,1,3]] origin protocol trajectories with **at most n** injected
2-qubit Pauli faults (after each CX/CZ: either no fault, or one of 15 non-identity
X/Z flips on control/target). Default n=1 (same as former one_gate_fault_dfs_513).

Matches the flow of `proof_protocol.py` / `boolean_learning_output_C.txt` naming:
  - Rounds r_0..r_3: one stabilizer circuit per round (XZZXI, IXZZX, XIXZZ, ZXIXZ).
    Fault sites: r_k_faulty_gate{0..5}_{x0,z0,x1,z1} (gate index = position in that .qasm).
  - Round r_4: full raw syndrome QASM. Fault sites: r_4_faulty_gate{0..N-1}_*.

Readout (aligned with proof_protocol symbolic conventions):
  - Syndrome s (ancX): Z-basis  -> use q.z
  - Flag f (flagZ):   X-basis  -> use q.x
  - Raw ancX[i]:      Z-basis  -> use q.z

Config: `[[5,1,3]]_origin/[[5,1,3]]_origin_config.txt` (paths to stab-specific QASM + raw).
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import dataclass
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
    decoder_file_index,
    format_decoder_commute_fail_report,
    verify_decoder_commute,
)
from flag_analysis import (
    CircuitXZ,
    apply_qasm_gate_into_state,
    detect_qubit_groups,
    new_clean_circuit_state,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

ORIGIN_DIR = PROJECT_ROOT / "[[5,1,3]]_origin"
CONFIG_PATH = ORIGIN_DIR / "[[5,1,3]]_origin_config.txt"
STAB_TXT = ORIGIN_DIR / "[[5,1,3]]_origin.txt"
LOG_TXT = ORIGIN_DIR / "[[5,1,3]]_log_op.txt"

# Stabilizer order matches `protocols/origin_5_1_3_protocol.json` nodes.
STAB_INSTRS = ("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")


# ---------------------------------------------------------------------------
# Config loader (key = value lines)
# ---------------------------------------------------------------------------


def _load_config(path: Path) -> Dict[str, str]:
    cfg: Dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        cfg[k.strip()] = v.strip()
    return cfg


# ---------------------------------------------------------------------------
# Pauli / helpers
# ---------------------------------------------------------------------------


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


# 15 non-identity (ctrl_x, ctrl_z, targ_x, targ_z) flip patterns
_FAULT15: List[Tuple[bool, bool, bool, bool]] = [
    (a, b, c, d)
    for a in (False, True)
    for b in (False, True)
    for c in (False, True)
    for d in (False, True)
    if (a or b or c or d)
]


def _apply_2q_with_optional_fault(
    state: CircuitXZ,
    name: str,
    ctrl: int,
    targ: int,
    fault_choice: Optional[Tuple[bool, bool, bool, bool]],
    fault_count: int,
    max_faults: int,
) -> int:
    """Ideal 2q Clifford, then optional after-gate Pauli. Returns new fault_count."""
    apply_qasm_gate_into_state(state, name, [ctrl, targ])
    if fault_choice is None:
        return fault_count
    if fault_count >= max_faults:
        return fault_count
    fc, fz, tc, tz = fault_choice
    _flip_qubit_xz(state, ctrl, fc, fz)
    _flip_qubit_xz(state, targ, tc, tz)
    return fault_count + 1


def _copy_data_error_only(src_state: CircuitXZ, src_qc: QuantumCircuit, dst_qc: QuantumCircuit) -> CircuitXZ:
    src_groups = detect_qubit_groups(src_qc)
    dst_groups = detect_qubit_groups(dst_qc)
    src_data = src_groups["data"]
    dst_data = dst_groups["data"]
    if len(src_data) != len(dst_data):
        raise ValueError("Data qubit count mismatch between circuits.")
    dst = new_clean_circuit_state(dst_qc.num_qubits)
    for s_idx, d_idx in zip(src_data, dst_data):
        dst.qubits[d_idx].x = src_state.qubits[s_idx].x
        dst.qubits[d_idx].z = src_state.qubits[s_idx].z
    return dst


def _data_xz_lists(state: CircuitXZ, qc: QuantumCircuit) -> Tuple[List[bool], List[bool]]:
    groups = detect_qubit_groups(qc)
    data_idxs = sorted(groups["data"])
    return (
        [bool(_to_bit(state.qubits[i].x)) for i in data_idxs],
        [bool(_to_bit(state.qubits[i].z)) for i in data_idxs],
    )


# ---------------------------------------------------------------------------
# Fault placement (for reports + SMT block)
# ---------------------------------------------------------------------------


@dataclass
class TwoQubitFaultPlacement:
    """One injected 2q fault event (a trajectory may have up to n of these)."""

    phase: str  # "flag" | "raw"
    stab_or_seg: int  # stabilizer index 0..3, or raw segment 0..3
    gate_index: int  # index in qc.data for that circuit (global within that file)
    gate_name: str
    ctrl: int
    targ: int
    flips: Tuple[bool, bool, bool, bool]  # ctrl_x, ctrl_z, targ_x, targ_z

    def describe(self, qc: QuantumCircuit) -> str:
        return (
            f"{self.phase} stab/seg={self.stab_or_seg} gate[{self.gate_index}]={self.gate_name} "
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


def _fault_detail_lines_from_placements(
    placements: List[TwoQubitFaultPlacement],
    flag_qc: QuantumCircuit,
    raw_qc: QuantumCircuit,
) -> List[str]:
    """Structured lines for format_decoder_commute_fail_report (no leading indent)."""
    lines: List[str] = []
    for p in placements:
        qc = raw_qc if p.phase == "raw" else flag_qc
        if p.phase == "flag":
            lines.append(f"first_subround stabilizer_slice = {p.stab_or_seg}")
        else:
            lines.append(f"raw_subround segment_index = {p.stab_or_seg}")
        instr, qargs, _ = qc.data[p.gate_index]
        ci = qc.find_bit(qargs[0]).index
        ti = qc.find_bit(qargs[1]).index
        cp = TwoQubitFaultPlacement._qubit_pretty(qc, ci)
        tp = TwoQubitFaultPlacement._qubit_pretty(qc, ti)
        lines.append(f"qc.data[{p.gate_index}]  {instr.name}  {cp}, {tp}")
        fc, fz, tgx, tgz = p.flips
        lines.append(
            f"after_gate flips (ctrl then targ, x/z): ({int(fc)},{int(fz)}) ({int(tgx)},{int(tgz)})"
        )
    return lines


# ---------------------------------------------------------------------------
# DFS over one [start, end) slice of qc.data
# ---------------------------------------------------------------------------

GateCont = Callable[[CircuitXZ, int, List[TwoQubitFaultPlacement]], None]


def _dfs_gate_slice(
    qc: QuantumCircuit,
    state: CircuitXZ,
    start: int,
    end: int,
    lo: int,
    fault_count: int,
    max_faults: int,
    placements: List[TwoQubitFaultPlacement],
    round_tag: str,
    stab_or_seg: int,
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
            qc,
            state,
            start,
            end,
            lo + 1,
            fault_count,
            max_faults,
            placements,
            round_tag,
            stab_or_seg,
            cont,
            stop_flag,
        )
        if stop_flag and stop_flag[0]:
            return
        return

    if lname in ("cx", "cnot", "cz"):
        c, t = qidxs[0], qidxs[1]
        if fault_count >= max_faults:
            choices: List[Optional[Tuple[bool, bool, bool, bool]]] = [None]
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
                        phase=round_tag,
                        stab_or_seg=stab_or_seg,
                        gate_index=lo,
                        gate_name=name,
                        ctrl=c,
                        targ=t,
                        flips=fc,
                    )
                ]
            _dfs_gate_slice(
                qc,
                st,
                start,
                end,
                lo + 1,
                nfc,
                max_faults,
                npls,
                round_tag,
                stab_or_seg,
                cont,
                stop_flag,
            )
            if stop_flag and stop_flag[0]:
                return
        return

    # 1-qubit Clifford
    apply_qasm_gate_into_state(state, name, qidxs)
    _dfs_gate_slice(
        qc,
        state,
        start,
        end,
        lo + 1,
        fault_count,
        max_faults,
        placements,
        round_tag,
        stab_or_seg,
        cont,
        stop_flag,
    )
    if stop_flag and stop_flag[0]:
        return


# ---------------------------------------------------------------------------
# SMT fault asserts (paste into Z3 / output_C debugging)
# ---------------------------------------------------------------------------


def build_smt_fault_asserts_for_path(
    *,
    placements: List[TwoQubitFaultPlacement],
    gates_per_stab: int,
    raw_num_gates: int,
) -> str:
    """
    Emit (assert (= r_... true/false)) lines matching `boolean_learning_output_C.txt` fault var names.
    Multiple faults: each gate site's x0/z0/x1/z1 set from that event; others false.
    """
    vals: Dict[str, bool] = {}
    for k in range(4):
        for g in range(gates_per_stab):
            p = f"r_{k}_faulty_gate{g}"
            for nm in (f"{p}_x0", f"{p}_z0", f"{p}_x1", f"{p}_z1"):
                vals[nm] = False
    for g in range(raw_num_gates):
        p = f"r_4_faulty_gate{g}"
        for nm in (f"{p}_x0", f"{p}_z0", f"{p}_x1", f"{p}_z1"):
            vals[nm] = False

    for placement in placements:
        fc, fz, tc, tz = placement.flips
        g = placement.gate_index
        if placement.phase == "flag":
            rk = placement.stab_or_seg
            p = f"r_{rk}_faulty_gate{g}"
        else:
            p = f"r_4_faulty_gate{g}"
        vals[f"{p}_x0"] = fc
        vals[f"{p}_z0"] = fz
        vals[f"{p}_x1"] = tc
        vals[f"{p}_z1"] = tz

    lines = [f"(assert (= {nm} {'true' if v else 'false'}))" for nm, v in sorted(vals.items())]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Protocol DFS driver
# ---------------------------------------------------------------------------


@dataclass
class DfsStats:
    leaves_end: int = 0
    leaves_lut: int = 0
    decoder_pass: int = 0
    decoder_fail: int = 0


def _prepare_stab_state(prev_state: Optional[CircuitXZ], prev_qc: Optional[QuantumCircuit], qc: QuantumCircuit) -> CircuitXZ:
    if prev_state is None:
        return new_clean_circuit_state(qc.num_qubits)
    assert prev_qc is not None
    return _copy_data_error_only(prev_state, prev_qc, qc)


def dfs_protocol(
    cfg: Dict[str, str],
    *,
    max_2q_faults: int = 1,
    verbose: bool = False,
    print_fail_decoder_vars: bool = False,
    max_fail_prints: int = 5,
    stop_on_first_fail: bool = False,
) -> DfsStats:
    if max_2q_faults < 0:
        raise ValueError("max_2q_faults (n) must be >= 0")

    stats = DfsStats()
    fail_prints = 0
    # When True, set to True on first decoder commute FAIL; DFS checks this and unwinds.
    stop_flag: Optional[List[bool]] = [False] if stop_on_first_fail else None

    stab_paths = [Path(cfg[instr]) for instr in STAB_INSTRS]
    raw_path = Path(cfg["raw_syndrome"])
    stab_qcs = [load_qasm(str(p)) for p in stab_paths]
    raw_qc = load_qasm(str(raw_path))

    n_gates_flag = len(stab_qcs[0].data)
    if any(len(qc.data) != n_gates_flag for qc in stab_qcs):
        raise ValueError("All stabilizer QASM files must have the same gate list length.")
    gates_per_stab = n_gates_flag

    n_raw = len(raw_qc.data)
    if n_raw % 4 != 0:
        raise ValueError("Raw circuit gate count must be divisible by 4.")
    gates_per_raw_seg = n_raw // 4

    def run_raw_round(
        first_state: CircuitXZ,
        flag_qc_ref: QuantumCircuit,
        first_stab: int,
        placements_flag: List[TwoQubitFaultPlacement],
        fault_count_in: int,
    ) -> None:
        nonlocal fail_prints
        second_state = _copy_data_error_only(first_state, flag_qc_ref, raw_qc)
        raw_idx = _build_reg_index_map(raw_qc)
        raw_anc = [raw_idx[("ancX", i)] for i in range(4)]

        def finish_raw(st: CircuitXZ, fc: int, pl_raw: List[TwoQubitFaultPlacement]) -> None:
            nonlocal fail_prints
            sec_s = [_to_bit(st.qubits[raw_anc[j]].z) for j in range(4)]
            s0 = _to_bit(first_state.qubits[_build_reg_index_map(flag_qc_ref)[("ancX", 0)]].z)
            f0 = _to_bit(first_state.qubits[_build_reg_index_map(flag_qc_ref)[("flagZ", 0)]].x)
            bitstring6 = f"{s0}{f0}{''.join(str(b) for b in sec_s)}"
            dx, dz = _data_xz_lists(first_state, flag_qc_ref)
            dec_idx = decoder_file_index(first_stab)
            dec_path = ORIGIN_DIR / "decoder" / f"path_{dec_idx}.txt"
            ok, _ = verify_decoder_commute(
                first_stabilizer_index=first_stab,
                bitstring6=bitstring6,
                data_x=dx,
                data_z=dz,
                decoder_c_path=dec_path,
                log_path=LOG_TXT,
                stab_path=STAB_TXT,
            )
            stats.leaves_lut += 1
            if ok:
                stats.decoder_pass += 1
            else:
                is_first_decoder_fail = stats.decoder_fail == 0
                stats.decoder_fail += 1
                if stop_flag is not None:
                    stop_flag[0] = True
                should_print = (print_fail_decoder_vars and fail_prints < max_fail_prints) or (
                    stop_on_first_fail and is_first_decoder_fail
                )
                if should_print:
                    smt = build_smt_fault_asserts_for_path(
                        placements=pl_raw,
                        gates_per_stab=gates_per_stab,
                        raw_num_gates=n_raw,
                    )
                    fault_lines = _fault_detail_lines_from_placements(pl_raw, flag_qc_ref, raw_qc)
                    print(
                        format_decoder_commute_fail_report(
                            first_stabilizer_index=first_stab,
                            bitstring6=bitstring6,
                            data_x=dx,
                            data_z=dz,
                            decoder_c_path=dec_path,
                            log_path=LOG_TXT,
                            stab_path=STAB_TXT,
                            fault_detail_lines=fault_lines if fault_lines else None,
                            smt_fault_block=smt,
                            stopped_after_first_fail=bool(stop_on_first_fail and is_first_decoder_fail),
                        )
                    )
                    fail_prints += 1
            if verbose:
                print(
                    f"LUT stab={first_stab} 6bit={bitstring6} dec_C_{dec_idx}.txt "
                    f"commute={'PASS' if ok else 'FAIL'} fault_count={fc}"
                )

        def seg_loop(seg_i: int, st: CircuitXZ, fc: int, pls: List[TwoQubitFaultPlacement]) -> None:
            if stop_flag and stop_flag[0]:
                return
            if seg_i >= 4:
                finish_raw(st, fc, pls)
                return
            st.qubits[raw_anc[seg_i]].x = BoolVal(False)
            st.qubits[raw_anc[seg_i]].z = BoolVal(False)
            a = seg_i * gates_per_raw_seg
            b = (seg_i + 1) * gates_per_raw_seg

            def after_seg(st2: CircuitXZ, fc2: int, pls2: List[TwoQubitFaultPlacement]) -> None:
                seg_loop(seg_i + 1, st2, fc2, pls2)

            _dfs_gate_slice(
                raw_qc,
                st,
                a,
                b,
                a,
                fc,
                max_2q_faults,
                pls,
                "raw",
                seg_i,
                after_seg,
                stop_flag,
            )

        seg_loop(0, second_state, fault_count_in, list(placements_flag))

    def run_from_stab(
        stab_i: int,
        prev_state: Optional[CircuitXZ],
        prev_qc: Optional[QuantumCircuit],
        fault_count: int,
        placements_so_far: List[TwoQubitFaultPlacement],
    ) -> None:
        if stop_flag and stop_flag[0]:
            return
        if stab_i >= 4:
            stats.leaves_end += 1
            if verbose:
                print("protocol end (no LUT): all [s,f]==[0,0]")
            return

        qc = stab_qcs[stab_i]
        st0 = _prepare_stab_state(prev_state, prev_qc, qc)
        idx = _build_reg_index_map(qc)
        anc_s_idx = idx[("ancX", 0)]
        flag_f_idx = idx[("flagZ", 0)]

        def after_slice(st: CircuitXZ, fc: int, pls: List[TwoQubitFaultPlacement]) -> None:
            merged = placements_so_far + pls
            s = _to_bit(st.qubits[anc_s_idx].z)
            f = _to_bit(st.qubits[flag_f_idx].x)
            if s == 0 and f == 0:
                run_from_stab(stab_i + 1, st, qc, fc, merged)
            else:
                run_raw_round(st, qc, stab_i, merged, fc)

        _dfs_gate_slice(
            qc,
            st0,
            0,
            len(qc.data),
            0,
            fault_count,
            max_2q_faults,
            [],
            "flag",
            stab_i,
            after_slice,
            stop_flag,
        )

    run_from_stab(0, None, None, 0, [])
    return stats


def main() -> None:
    ap = argparse.ArgumentParser(
        description="DFS over [[5,1,3]] origin protocol with at most N 2q faults per trajectory."
    )
    ap.add_argument(
        "n",
        nargs="?",
        type=int,
        default=1,
        metavar="N",
        help="Max number of 2q faults allowed per run (default: 1)",
    )
    ap.add_argument("--config", type=Path, default=CONFIG_PATH, help="Path to origin config txt")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument(
        "--print-fail-decoder-vars",
        action="store_true",
        help="On decoder commute FAIL, print meas/dec/data and SMT fault assert block (limited count).",
    )
    ap.add_argument("--max-fail-print", type=int, default=5, help="Max FAIL reports when printing.")
    ap.add_argument(
        "--stop-on-first-fail",
        action="store_true",
        help="On first decoder commute FAIL, stop exploring further trajectories.",
    )
    args = ap.parse_args()

    if args.n < 0:
        ap.error("N must be >= 0")

    cfg = _load_config(args.config)
    stats = dfs_protocol(
        cfg,
        max_2q_faults=args.n,
        verbose=args.verbose,
        print_fail_decoder_vars=args.print_fail_decoder_vars,
        max_fail_prints=args.max_fail_print,
        stop_on_first_fail=args.stop_on_first_fail,
    )
    print(
        f"max_2q_faults={args.n} | stop_on_first_fail={args.stop_on_first_fail} | "
        f"leaves end (no LUT): {stats.leaves_end}  |  leaves LUT: {stats.leaves_lut}  |  "
        f"decoder commute PASS: {stats.decoder_pass}  FAIL: {stats.decoder_fail}"
    )


if __name__ == "__main__":
    main()
