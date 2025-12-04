from qiskit import QuantumCircuit
from qiskit.qasm2 import dumps as qasm2_dumps



import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
try:
    # Some 2.x installs expose a dedicated qasm2 loader
    from qiskit.qasm2 import loads as qasm2_loads  # type: ignore
    _HAS_QASM2 = True
except Exception:
    _HAS_QASM2 = False



def load_qasm(qasm_path: str) -> QuantumCircuit:
    """
    Loader shim for Qiskit 2.x: try the classic from_qasm_file first;
    if your build expects qasm2 loader, use it as a fallback.
    """
    try:
        return QuantumCircuit.from_qasm_file(qasm_path)
    except Exception:
        if _HAS_QASM2:
            with open(qasm_path, "r", encoding="utf-8") as f:
                txt = f.read()
            return qasm2_loads(txt)
        raise  # rethrow if no fallback available

def get_gate_only_indices(qc):
    """
    Return a list of indices in qc.data corresponding ONLY to real gates.
    Skips: measure, barrier, reset, id.
    """
    skip = {"measure", "barrier", "reset", "id"}
    gate_indices = []

    for i, inst in enumerate(qc.data):
        if inst.operation.name.lower() not in skip:
            gate_indices.append(i)

    return gate_indices


def _build_bit_lookup(qc: QuantumCircuit):
    """
    Return dict: { QubitObject -> (register_name, local_index) }
    Works across Qiskit 1.x/2.x because it iterates the circuit's qregs directly.
    """
    m = {}
    for reg in qc.qregs:
        # reg is a QuantumRegister; iterating gives Bit objects in order
        for j, bit in enumerate(reg):
            m[bit] = (reg.name, j)
    return m

def split_circuit_full_q_compact_flags(qc: QuantumCircuit):
    """
    Split `qc` by 'barrier'. In each slice:
      - keep full data register 'q' (same size),
      - create compact ancilla/flag regs containing only used bits for that slice,
      - ignore 'measure' ops; no classical registers created.
      - if the original slice ended at a barrier, append a barrier in the subcircuit.
    Returns: List[QuantumCircuit]
    """
    # 1) collect ops between barriers (skip measures), and remember if a slice ended at a barrier
    slices = []  # list of (ops, ended_by_barrier)
    cur = []
    for instr, qargs, cargs in qc.data:
        if instr.name == "barrier":
            if cur:
                slices.append((cur, True))  # this slice ended due to a barrier
                cur = []
        elif instr.name == "measure":
            continue
        else:
            cur.append((instr, qargs))
    if cur:
        slices.append((cur, False))  # last slice (no trailing barrier)

    # 2) map each Qubit object -> (register_name, local_index)
    bit_lookup = _build_bit_lookup(qc)

    # 3) original data register (kept full size)
    q_reg_orig = next((r for r in qc.qregs if r.name == "q"), None)
    if q_reg_orig is None:
        raise ValueError("Expected a data register named 'q'.")
    q_size = q_reg_orig.size

    subcircuits = []
    for k, (ops, ended_by_barrier) in enumerate(slices):
        # which anc/flag locals are used in this slice
        used_by_reg = {"ancX": set(), "ancZ": set(), "flagX": set(), "flagZ": set()}

        for instr, qargs in ops:
            for qb in qargs:
                rname, lidx = bit_lookup[qb]
                if rname in used_by_reg:
                    used_by_reg[rname].add(lidx)

        # build new regs: full 'q', compact anc/flag only if used
        q_new = QuantumRegister(q_size, "q")
        regs = [q_new]
        remap = {("q", j): q_new[j] for j in range(q_size)}

        def add_compact_reg(reg_name):
            idxs = sorted(used_by_reg[reg_name])
            if not idxs:
                return
            R = QuantumRegister(len(idxs), reg_name)
            regs.append(R)
            for new_i, old_i in enumerate(idxs):
                remap[(reg_name, old_i)] = R[new_i]

        for rn in ("ancX", "ancZ", "flagX", "flagZ"):
            add_compact_reg(rn)

        # build subcircuit and append remapped ops
        sub = QuantumCircuit(*regs, name=f"stab_{k}_fullq_compactflags")
        for instr, qargs in ops:
            new_qargs = []
            for qb in qargs:
                rname, lidx = bit_lookup[qb]
                key = (rname, lidx)
                if key not in remap:
                    raise KeyError(f"Unmapped bit {rname}[{lidx}] in slice {k}")
                new_qargs.append(remap[key])
            sub.append(instr, new_qargs, [])

        # append a barrier if the original block ended at a barrier
        if ended_by_barrier:
            sub.barrier(*sub.qubits)

        subcircuits.append(sub)

    return subcircuits


def save_qasm_full_slices(qc, prefix="stab"):
    """Split by barriers and save each stabilizer as a .qasm file (Qiskit 2.x compatible)."""
    subs = split_circuit_full_q_compact_flags(qc)
    for i, sub in enumerate(subs):
        filename = f"{prefix}_{i}.qasm"
        try:
            # Preferred (Qiskit 2.x)
            qasm_str = qasm2_dumps(sub)
        except Exception:
            # Fallback for Qiskit 1.x
            qasm_str = sub.qasm()
        with open(filename, "w", encoding="utf-8") as f:
            f.write(qasm_str)
        print(f"Saved: {filename}")

def remove_flag_gates(qasm_path: str, save_path: str = None):

    """
    Load a QASM file, remove all gates that act on flag qubits (flagX[...] or flagZ[...]),
    but preserve the original `barrier` gates.
    """
    qc = QuantumCircuit.from_qasm_file(qasm_path)
    new_qc = QuantumCircuit(*qc.qregs, *qc.cregs)

    # Map: global index -> register name
    regmap = {}
    idx = 0
    for qreg in qc.qregs:
        for _ in range(qreg.size):
            regmap[idx] = qreg.name.lower()
            idx += 1

    def is_flag_qubit(qbit):
        """Check if the given Qubit belongs to a flag register."""
        loc = qc.find_bit(qbit)
        reg_name = regmap.get(loc.index, "")
        return reg_name.startswith("flagx") or reg_name.startswith("flagz")

    # Filter out gates acting on flag qubits, but keep barriers
    for instr, qargs, cargs in qc.data:
        if instr.name == "barrier":
            # Always include barrier gates
            new_qc.append(instr, qargs, cargs)
        elif any(is_flag_qubit(q) for q in qargs):
            # Skip gates that act on flag qubits
            continue
        else:
            # Include all other gates
            new_qc.append(instr, qargs, cargs)

    # Optionally save the modified circuit
    if save_path:
        with open(save_path, "w") as f:
            f.write(qasm2_dumps(new_qc))  # Use qasm2_dumps to generate QASM string

    return new_qc