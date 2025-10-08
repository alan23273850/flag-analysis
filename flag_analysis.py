# flag_analysis.py
# Minimal Pauli-flow utilities focused on the FLAG qubit.
# Tested with Qiskit 2.x (with compatibility shims).
# Requires: pip install qiskit z3-solver

from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional


# --- Qiskit imports + loader shim ---
from qiskit import QuantumCircuit
try:
    # Some 2.x installs expose a dedicated qasm2 loader
    from qiskit.qasm2 import loads as qasm2_loads  # type: ignore
    _HAS_QASM2 = True
except Exception:
    _HAS_QASM2 = False

from z3 import BoolVal, Xor, Bool,simplify,substitute, And, Not,Or, PbLe, AtMost
from z3 import Solver, unsat, sat
# ---------------------------
# Pauli-flow data structures
# ---------------------------

@dataclass
class QubitXZ:
    x: object  # z3 BoolRef (or BoolVal)
    z: object

@dataclass
class CircuitXZ:
    qubits: List[QubitXZ]

# ---------------------------
# Boolean helpers
# ---------------------------

def bfalse(): return BoolVal(False)
def bxor(a, b): return Xor(a, b)

# ---------- Qiskit 2.2 regmap helper ----------
def _regmap_indices(qc: QuantumCircuit):
    """Return dict: global_qubit_index -> (qreg_name, local_index)."""
    regmap = {}
    idx = 0
    for qreg in qc.qregs:
        for j in range(qreg.size):
            regmap[idx] = (qreg.name, j)
            idx += 1
    return regmap

# ---------------------------
# Init
# ---------------------------

def new_clean_circuit_state(n_qubits: int) -> CircuitXZ:
    """Start with no errors anywhere: X=0, Z=0 per qubit."""
    return CircuitXZ([QubitXZ(bfalse(), bfalse()) for _ in range(n_qubits)])


def new_variable_circuit_state(qc: QuantumCircuit) -> CircuitXZ:
    """
    Initialize each qubit with named Bool variables based on qreg name + local index.
    Produces variables like: q0_x, q0_z, ancX1_x, ancX1_z, flagZ2_x, ...
    """
    regmap = _regmap_indices(qc)
    qubits: List[QubitXZ] = []
    for i in range(qc.num_qubits):
        regname, j = regmap[i]
        prefix = f"{regname}{j}"
        qubits.append(QubitXZ(x=Bool(f"{prefix}_x"), z=Bool(f"{prefix}_z")))
    return CircuitXZ(qubits)



# ---------------------------
# Group detection by register name
# ---------------------------

def detect_qubit_groups(qc: QuantumCircuit) -> Dict[str, List[int]]:
    """
    Group qubits by their register names for Qiskit 2.2 style.
    Each qreg has a name (like 'ancX', 'ancZ', 'flagX', 'flagZ', 'q'),
    and qubits are numbered globally across all registers.
    """
    groups = {'data': [], 'ancX': [], 'ancZ': [], 'flagX': [], 'flagZ': []}

    # Build mapping: global index → (register name, local index)
    idx = 0
    regmap = {}
    for qreg in qc.qregs:
        for j in range(qreg.size):
            regmap[idx] = qreg.name
            idx += 1

    # Classify each qubit index
    for i, reg in regmap.items():
        reg_l = reg.lower()
        if   reg_l.startswith("ancx"):  groups["ancX"].append(i)
        elif reg_l.startswith("ancz"):  groups["ancZ"].append(i)
        elif reg_l.startswith("flagx"): groups["flagX"].append(i)
        elif reg_l.startswith("flagz"): groups["flagZ"].append(i)
        else: groups["data"].append(i)

    return groups

# ---------------------------
# Clifford update rules
# ---------------------------

def apply_h(state: CircuitXZ, q: int) -> None:
    """Hadamard on qubit q: (x,z) <- (z,x)."""
    state.qubits[q].x, state.qubits[q].z = state.qubits[q].z, state.qubits[q].x

def apply_s(state: CircuitXZ, q: int) -> None:
    """Phase S on qubit q: (x,z) <- (x, x xor z)."""
    x, z = state.qubits[q].x, state.qubits[q].z
    state.qubits[q].z = bxor(x, z)

def apply_sdg(state: CircuitXZ, q: int) -> None:
    """S† on qubit q: (x,z) <- (x, z xor x).  (inverse of S)"""
    x, z = state.qubits[q].x, state.qubits[q].z
    state.qubits[q].z = bxor(z, x)

def apply_cnot(state: CircuitXZ, ctrl: int, targ: int) -> None:
    """
    CNOT(c->t):
      x_c' = x_c
      z_c' = z_c xor z_t
      x_t' = x_t xor x_c
      z_t' = z_t
    """
    xc, zc = state.qubits[ctrl].x, state.qubits[ctrl].z
    xt, zt = state.qubits[targ].x, state.qubits[targ].z
    state.qubits[ctrl].x = xc
    state.qubits[ctrl].z = bxor(zc, zt)
    state.qubits[targ].x = bxor(xt, xc)
    state.qubits[targ].z = zt

def apply_notnot(state: CircuitXZ, ctrl: int, targ: int) -> None:
    """
    NOTNOT(c->t):
      x_c' = x_c
      z_c' = z_c xor z_t
      x_t' = x_t xor x_c
      z_t' = z_t
    """
    xc, zc = state.qubits[ctrl].x, state.qubits[ctrl].z
    xt, zt = state.qubits[targ].x, state.qubits[targ].z
    state.qubits[ctrl].x = bxor(xc, zt)
    state.qubits[ctrl].z = zc
    state.qubits[targ].x = bxor(xt, zc)
    state.qubits[targ].z = zt

# ---------------------------
# Fault injection on FLAG
# ---------------------------

def inject_flag_error(state: CircuitXZ, flag_idx: int, kind: str) -> None:
    """
    Insert a Pauli error on the flag qubit at the *current time*.
    kind ∈ {'I','X','Z','Y'}.
    """
    k = kind.upper()
    if k == 'I':
        return
    if k in ('X', 'Y'):
        state.qubits[flag_idx].x = bxor(state.qubits[flag_idx].x, BoolVal(True))
    if k in ('Z', 'Y'):
        state.qubits[flag_idx].z = bxor(state.qubits[flag_idx].z, BoolVal(True))

def inject_flag_symbolic_one_axis(state, fidx: int, axis="z", prefix="ferr"):
    """
    Add a symbolic error variable on a flag qubit, restricted to X *or* Z axis.

    axis: "x" or "z"
    """
    var = Bool(f"{prefix}{fidx}_{axis}")

    if axis == "x":
        # Only X part of error
        state.qubits[fidx].x = Xor(state.qubits[fidx].x, var)
    elif axis == "z":
        # Only Z part of error
        state.qubits[fidx].z = Xor(state.qubits[fidx].z, var)
    else:
        raise ValueError("axis must be 'x' or 'z'")

    return var
def _inject_1q_fault_after(state, q, fault_kind=None, prefix="f"):
    """
    One-qubit Pauli on wire q after a 1q gate (or on ONE wire of a 2q gate):
      fault_kind: None -> symbolic; 'I'|'X'|'Z'|'Y' -> concrete.
    Returns dict {'fx','fz','act'}.
    """
    if fault_kind is None:
        fx = Bool(f"{prefix}_x"); fz = Bool(f"{prefix}_z")
    else:
        k = fault_kind.upper()
        fx = BoolVal(k in ("X","Y")); fz = BoolVal(k in ("Z","Y"))
    state.qubits[q].x = Xor(state.qubits[q].x, fx)
    state.qubits[q].z = Xor(state.qubits[q].z, fz)
    return {"fx": fx, "fz": fz, "act": Or(fx, fz)}

def _inject_2q_fault_after(state, q0: int, q1: int, fault_kind=None, prefix="f"):
    """
    Gate-agnostic 2-qubit Pauli injection on wires (q0, q1) *after* a 2q gate.
    It does not depend on the specific 2q gate; it simply toggles X/Z components.

    Args:
      state: CircuitXZ (your Pauli-flow state)
      q0, q1: qubit indices (in the circuit's global indexing)
      fault_kind:
        - None               -> symbolic on both wires
        - (k0, k1)           -> concrete per-wire, each in {'I','X','Z','Y'} (case-insensitive)
          e.g. ('Y','X') means inject Y on q0 and X on q1
      prefix: name prefix for z3 symbols (when symbolic)

    Returns:
      info: dict with z3 literals and helpers:
        {
          'fx0','fz0','fx1','fz1',   # per-wire Pauli indicator bits
          'act0','act1',             # per-wire activity (X or Z non-identity)
          'act'                      # any activity on either wire
        }
    """
    if fault_kind is None:
        fx0 = Bool(f"{prefix}_x0"); fz0 = Bool(f"{prefix}_z0")
        fx1 = Bool(f"{prefix}_x1"); fz1 = Bool(f"{prefix}_z1")
    else:
        k0, k1 = fault_kind
        k0 = k0.upper(); k1 = k1.upper()
        fx0 = BoolVal(k0 in ("X","Y")); fz0 = BoolVal(k0 in ("Z","Y"))
        fx1 = BoolVal(k1 in ("X","Y")); fz1 = BoolVal(k1 in ("Z","Y"))

    # Inject on both wires (gate-agnostic)
    state.qubits[q0].x = Xor(state.qubits[q0].x, fx0)
    state.qubits[q0].z = Xor(state.qubits[q0].z, fz0)
    state.qubits[q1].x = Xor(state.qubits[q1].x, fx1)
    state.qubits[q1].z = Xor(state.qubits[q1].z, fz1)

    act0 = Or(fx0, fz0); act1 = Or(fx1, fz1)
    info = {
        "fx0": fx0, "fz0": fz0, "fx1": fx1, "fz1": fz1,
        "act0": act0, "act1": act1,
        "act": Or(act0, act1),
        
    }
    return info

def add_fault_mode_constraints(solver, info, fault_mode="2q", fault_kind=None):
    """
    Add constraints to the solver according to the declared fault_mode.

    - "1q": at most one of {act_c, act_t} is true
    - "2q": both act_c and act_t are true
    - "either": no extra restriction
    """
    if fault_mode == "1q":
        solver.add(AtMost(info["act_c"], info["act_t"], 1))
    elif fault_mode == "2q":
        solver.add(Or(info["act_c"], info["act_t"]))
    # "either": do nothing
# ---------------------------
# QASM → Pauli-flow
# ---------------------------

def _qiskit_qubit_index(qc: QuantumCircuit, qobj) -> int:
    """
    Robustly get integer index from a Qiskit Qubit object in 1.x/2.x.
    """
    # find_bit returns a BitLocations with .index
    return qc.find_bit(qobj).index

def apply_qasm_gate_into_state(state: CircuitXZ, name: str, qidxs: List[int]) -> None:
    """Apply supported (Clifford) gates to Pauli-flow state."""
    if name == 'h':
        apply_h(state, qidxs[0])
    elif name == 's':
        apply_s(state, qidxs[0])
    elif name in ('sdg', 'sxdg'):  # sdg is the usual name; include alias just in case
        apply_sdg(state, qidxs[0])
    elif name in ('cx', 'cnot'):
        apply_cnot(state, qidxs[0], qidxs[1])

    elif name in ('notnot'):  # notnot is a common alias for ccx
        apply_notnot(state, qidxs[0], qidxs[1])

    elif name in ('id', 'barrier', 'reset', 'measure'):
        # ignored here; measurement is read via final Z/X bits directly
        pass
    else:
        raise NotImplementedError(f"Unsupported gate in Pauli-flow: {name}")

def _load_qasm(qasm_path: str) -> QuantumCircuit:
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

def build_state_from_qasm(qasm_path: str) -> Tuple[CircuitXZ, QuantumCircuit]:
    """Load QASM and walk the gates to produce final Pauli-flow state (no faults)."""
    qc = _load_qasm(qasm_path)
    st = new_clean_circuit_state(qc.num_qubits)
    for instr, qargs, _ in qc.data:
        name = instr.name
        qidxs = [_qiskit_qubit_index(qc, q) for q in qargs]
        apply_qasm_gate_into_state(st, name, qidxs)
    return st, qc

def build_variable_state_from_qasm(qasm_path: str) -> Tuple[CircuitXZ, QuantumCircuit, Dict[str, object]]:
    """
    Load QASM, build variable-initialized Pauli-flow state, propagate gates,
    and return (state, qc, varenv) where varenv maps variable names to z3 Bools.
    """
    qc = _load_qasm(qasm_path)

    # 1) variable-initialized state
    state = new_variable_circuit_state(qc)

    # 2) also return a varenv for easy substitutions/evaluation
    regmap = _regmap_indices(qc)
    varenv: Dict[str, object] = {}
    for i in range(qc.num_qubits):
        regname, j = regmap[i]
        prefix = f"{regname}{j}"
        varenv[f"{prefix}_x"] = state.qubits[i].x
        varenv[f"{prefix}_z"] = state.qubits[i].z

    # 3) walk the circuit (Clifford updates)
    for instr, qargs, _ in qc.data:
        name = instr.name
        qidxs = [_qiskit_qubit_index(qc, q) for q in qargs]
        apply_qasm_gate_into_state(state, name, qidxs)

    return state, qc, varenv


def build_state_with_fault_after_gate(qasm_path: str, gate_index: int, fault_mode="either", fault_kind=None):
    """
    Run circuit ideally; inject a fault right AFTER gate #gate_index only.
    fault_mode: '1q' | '2q' | 'either'  (for CNOTs)
    fault_kind:
      - None                -> symbolic
      - 'I'|'X'|'Z'|'Y'     -> for 1q gates, or for CNOT in mode='1q' (applied to one wire)
      - (kc,kt)             -> for CNOT in mode='2q'/'either' (concrete per-wire)
    Returns: (state, qc, site_info, groups)
    """
    qc = _load_qasm(qasm_path)
    state = new_clean_circuit_state(qc.num_qubits)
    groups = detect_qubit_groups(qc)
    site_info = None

    for i, (instr, qargs, _) in enumerate(qc.data):
        name = instr.name
        qidxs = [_qiskit_qubit_index(qc, q) for q in qargs]

        if name in ("h","s","sdg"):
            apply_qasm_gate_into_state(state, name, qidxs)
            if i == gate_index:
                info = _inject_1q_fault_after(
                    state, qidxs[0],
                    fault_kind=None if fault_kind is None else fault_kind,
                    prefix=f"f_site{i}"
                )
                site_info = {
                    "gate_index": i, "gate_name": name,
                    "qubits": (qidxs[0],),
                    "vars": info, "act": info["act"], "fault_mode": "1q"
                }

        elif name in ("cx","cnot", "notnot"):
            c, t = qidxs
            if (name == "cx"  or name == "cnot"):apply_cnot(state, c, t)
            else: apply_notnot(state, c, t)
            if i == gate_index:
                info = _inject_2q_fault_after(
                    state, c, t, 
                    fault_kind=None if fault_kind is None else fault_kind,
                    prefix=f"f_gate{i}"
                )
                site_info = {
                    "gate_index": i, "gate_name": "cx",
                    "qubits": (c, t),
                    "vars": info, "act": info["act"], "fault_mode": fault_mode
                }

        elif name in ("barrier","id","reset","measure"):
            if i == gate_index:
                site_info = {
                    "gate_index": i, "gate_name": name,
                    "qubits": tuple(qidxs),
                    "vars": {}, "act": BoolVal(False), "fault_mode": "none"
                }
        else:
            raise NotImplementedError(f"Unsupported gate: {name}")

    if site_info is None:
        raise IndexError(f"gate_index {gate_index} out of range (len={len(qc.data)})")

    return state, qc, site_info, groups


def build_state_with_faults_after_gates(qasm_path: str, gate_indices: list, fault_mode="either", fault_kind=None):
    """
    Run circuit ideally; inject faults right AFTER all gates in `gate_indices`.
    
    fault_mode: '1q' | '2q' | 'either'  (for 2-qubit gates)
    fault_kind:
      - None                -> symbolic
      - 'I'|'X'|'Z'|'Y'     -> for 1q gates, or for CNOT in mode='1q' (applied to one wire)
      - (kc,kt)             -> for 2q gates (concrete per-wire)
    
    Returns:
      state, qc, sites_info, groups
      where `sites_info` is a list of site_info dicts (one per injected fault)
    """
    qc = _load_qasm(qasm_path)
    state = new_clean_circuit_state(qc.num_qubits)
    groups = detect_qubit_groups(qc)
    sites_info = []

    gate_indices = set(gate_indices)  # so we can check membership quickly

    for i, (instr, qargs, _) in enumerate(qc.data):
        name = instr.name
        qidxs = [_qiskit_qubit_index(qc, q) for q in qargs]

        if name in ("h","s","sdg"):
            apply_qasm_gate_into_state(state, name, qidxs)
            if i in gate_indices:
                info = _inject_1q_fault_after(
                    state, qidxs[0],
                    fault_kind=None if fault_kind is None else fault_kind,
                    prefix=f"f_site{i}"
                )
                sites_info.append({
                    "gate_index": i, "gate_name": name,
                    "qubits": (qidxs[0],),
                    "vars": info, "act": info["act"], "fault_mode": "1q"
                })

        elif name in ("cx","cnot","notnot"):
            c, t = qidxs
            if name in ("cx","cnot"):
                apply_cnot(state, c, t)
            else:
                apply_notnot(state, c, t)
            if i in gate_indices:
                info = _inject_2q_fault_after(
                    state, c, t,
                    fault_kind=None if fault_kind is None else fault_kind,
                    prefix=f"faulty_gate{i}"
                )
                sites_info.append({
                    "gate_index": i, "gate_name": name,
                    "qubits": (c, t),
                    "vars": info, "act": info["act"], "fault_mode": fault_mode
                })

        elif name in ("barrier","id","reset","measure"):
            if i in gate_indices:
                sites_info.append({
                    "gate_index": i, "gate_name": name,
                    "qubits": tuple(qidxs),
                    "vars": {}, "act": BoolVal(False), "fault_mode": "none"
                })
        else:
            raise NotImplementedError(f"Unsupported gate: {name}")

    if not sites_info:
        raise IndexError(f"gate_indices {gate_indices} produced no injections (len={len(qc.data)})")

    return state, qc, sites_info, groups

def build_stab_equiv_errors(E_x, E_z, stab_txt_path, prefix="g"):
    """
    Construct stabilizer-equivalent errors.

    Args:
      E_x, E_z: lists of z3 Bool formulas for data qubits
      stab_txt_path: path to stabilizer .txt file
      prefix: name prefix for generator selector variables (default "g")

    Returns:
      (Epx, Epz, gsel)
        - Epx, Epz: new error expressions after applying all possible generator products
        - gsel: list of selector Bool variables, one per generator
    """
    # Load stabilizers from file
    gens = load_symplectic_txt(stab_txt_path)
    m = len(gens)     # number of generators
    n = len(E_x)      # number of data qubits

    # Create selector vars g0..g{m-1}
    gsel = [Bool(f"{prefix}{j}") for j in range(m)]

    # Collect which generators flip which qubit components
    addX = [[] for _ in range(n)]
    addZ = [[] for _ in range(n)]
    for j, (Sx, Sz) in enumerate(gens):
        gj = gsel[j]
        for i in range(n):
            if Sx[i]: addX[i].append(gj)   # Z on generator anticommutes with X error
            if Sz[i]: addZ[i].append(gj)   # X on generator anticommutes with Z error

    # Apply XOR modifications to each data qubit
    Epx, Epz = [], []
    for i in range(n):
        xi = E_x[i]
        for t in addX[i]:
            xi = Xor(xi, t)
        zi = E_z[i]
        for t in addZ[i]:
            zi = Xor(zi, t)
        Epx.append(xi)
        Epz.append(zi)

    return Epx, Epz, gsel
# ---------------------------
# Outputs we care about
# ---------------------------


def ancillas_Z(state: CircuitXZ, anc_idxs: List[int]):
    """Syndrome bits if ancillas are measured in Z basis (flips if X on ancilla)."""
    return [state.qubits[i].x for i in anc_idxs]

def ancillas_X(state: CircuitXZ, anc_idxs: List[int]):
    """Syndrome bits if ancillas are measured in X basis (flips if Z on ancilla)."""
    return [state.qubits[i].z for i in anc_idxs]

def flags_Z(state: CircuitXZ, flag_idxs: List[int]):
    """Flag measured in Z basis (flips if X on flag)."""
    return [state.qubits[i].x for i in flag_idxs]

def flags_X(state: CircuitXZ, flag_idxs: List[int]):
    """Flag measured in X basis (flips if Z on flag)."""
    return [state.qubits[i].z for i in flag_idxs]

def data_qubits(state: CircuitXZ, data_idxs: List[int]):
    """
    Return the residual Pauli error components (x,z) for each data qubit.
    Example:
      (False, True)  -> Z error
      (True, False)  -> X error
      (True, True)   -> Y error
      (False, False) -> no error
    """
    return [(state.qubits[i].x, state.qubits[i].z) for i in data_idxs]

def data_error_weight_literals(state, data_idxs):
    """
    Return a list of literals indicating whether each data qubit carries
    any non-trivial Pauli error (X or Z).
    Useful for counting error weight with PB constraints.
    """
    return [Or(state.qubits[i].x, state.qubits[i].z) for i in data_idxs]


def eval_under(boolexpr, assignment: dict, varenv: dict):
    """
    Evaluate a z3 Bool expression under a *partial assignment*.
    assignment: dict like {"q3_x": True, "ancX0_z": False}
    - Only these vars are substituted.
    - Any var not listed stays symbolic (instead of defaulting to False).
    """
    subs = []
    for name, val in assignment.items():
        if name in varenv:
            subs.append((varenv[name], BoolVal(val)))
    return simplify(substitute(boolexpr, subs))

def project_data_only(expr, varenv: dict):
    """
    Substitute anc/flag variables to False, keep all data variables symbolic.
    """
    subs = []
    for name, sym in varenv.items():
        if name.startswith(("ancX","ancZ","flagX","flagZ")):
            subs.append((sym, BoolVal(False)))
    return simplify(substitute(expr, subs))
# ---------------------------
# High-level helper:
#   What happens if the flag has an X/Z/Y error just before measurement?
# ---------------------------

def analyze_flag_errors_multi(
    qasm_path: str,
    anc_idxs: List[int],
    flag_idxs: List[int],
    flag_error_kinds: Dict[int, str],  # mapping: flag_idx -> 'I'|'X'|'Z'|'Y'
):
    """
    Build the Pauli-flow state from QASM, inject specified Pauli errors
    on one or more flag qubits (just before measurement), and return:
      - syn_flips:  List[BoolExpr]  (Z on each ancilla in anc_idxs)
      - flag_flips: List[BoolExpr]  (X-basis flip formula for each flag in flag_idxs)
    """
    state, _ = build_state_from_qasm(qasm_path)

    # Inject errors on the chosen flag qubits
    for fidx, kind in flag_error_kinds.items():
        inject_flag_error(state, fidx, kind)

    syn_flips  = ancillas_Z(state, anc_idxs)
    flag_flips = flags_X(state, flag_idxs)  # X-basis measurement assumed
    return syn_flips, flag_flips

# ---------------------------
# Stabilizers:
#  
# ---------------------------
def load_symplectic_txt(path: str):
    """
    Each line has 'XXXXXXX ZZZZZZZ' (0/1).
    Returns list of (Sx, Sz) where each is a list of 0/1.
    """
    gens = []
    with open(path, 'r') as f:
        for ln in f:
            if not ln.strip():
                continue
            xs, zs = ln.split()
            Sx = [int(c) for c in xs.strip()]
            Sz = [int(c) for c in zs.strip()]
            gens.append((Sx, Sz))
    return gens

def anticomm_formula(Sx, Sz, varenv):
    """
    Build z3 Bool formula:
      ⊕_i ( E_x[i]*S_z[i]  ⊕  E_z[i]*S_x[i] )
    where error vars are q{i}_x, q{i}_z (or data{i}_x/z).
    """
    acc = BoolVal(False)
    for i in range(len(Sx)):
        if Sz[i]:  # stabilizer has Z → anticommutes with X error
            acc = Xor(acc, varenv.get(f"q{i}_x", varenv.get(f"data{i}_x")))
        if Sx[i]:  # stabilizer has X → anticommutes with Z error
            #print("Adding Z term for qubit", i)
            acc = Xor(acc, varenv.get(f"q{i}_z", varenv.get(f"data{i}_z")))
        #print(f"Step {i}: acc = {acc}")
    return acc

# ---------------------------
# Ordered solver check: ancilla[i] ≡ stabilizer[i]
# ---------------------------

def _equiv(a, b) -> bool:
    """True iff a and b are logically equivalent (UNSAT of XOR)."""
    s = Solver()
    s.add(Xor(a, b))          # SAT means there exists an assignment where they differ
    return s.check() == unsat # UNSAT ⇒ no such assignment ⇒ equivalent

def _counterexample(a, b):
    """Return a counterexample model if a ≢ b, else None."""
    s = Solver()
    s.add(Xor(a, b))
    return s.model() if s.check() == sat else None


def check_ancillas_match_symplectic_ordered(qasm_path: str,
                                            stab_txt_path: str,
                                            order: str = "X-then-Z"):
    """
    Pairwise, ordered equivalence:
      ancilla[i]  ≡  anticommute_formula_from_txt_line[i]

    - Ancilla formulas are taken from the circuit and projected to data-only
      (anc/flag vars set False; data vars symbolic).
    - Stabilizer formulas come from the txt (file order preserved).
    - `order` tells how to concatenate ancillas from QASM registers:
         "X-then-Z" (default) means ancX first, then ancZ (both in QASM order).
         "Z-then-X" means ancZ first, then ancX.
      Choose the one that matches the line order in your .txt.
    """
    # Build circuit (symbolic) and detect groups
    state, qc, varenv = build_variable_state_from_qasm(qasm_path)
    groups = detect_qubit_groups(qc)

    # Ancilla flip formulas from circuit → project to data-only
    ancX = [project_data_only(e, varenv) for e in ancillas_X(state, groups["ancX"])]
    ancZ = [project_data_only(e, varenv) for e in ancillas_Z(state, groups["ancZ"])]

    ancillas = (ancX + ancZ) if order == "X-then-Z" else (ancZ + ancX)

    # Stabilizer anticommute formulas from txt (exact line order)
    gens = load_symplectic_txt(stab_txt_path)
    stabs = [anticomm_formula(Sx, Sz, varenv) for (Sx, Sz) in gens]

    if len(ancillas) != len(stabs):
        print(f"[COUNT MISMATCH] ancillas={len(ancillas)} vs stabs={len(stabs)}")
        return {"ok": False, "mismatches": list(range(min(len(ancillas), len(stabs))))}

    # Combine all equivalence checks into a single AND condition
    combined_condition = And(*[Xor(a, s) == False for a, s in zip(ancillas, stabs)])

    # Check if the combined condition is satisfied
    s = Solver()
    s.add(Not(combined_condition))  # Check if there exists a counterexample
    if s.check() == unsat:
        print("Overall ordered match: True")
        return {"ok": True, "mismatches": []}
    else:
        print("Overall ordered match: False")
        print("Counterexample:", s.model())
        mismatches = [i for i, (a, s) in enumerate(zip(ancillas, stabs)) if not _equiv(a, s)]
        return {"ok": False, "mismatches": mismatches}


def check_gate_k_with_fault(
    qasm_path: str,
    gate_index: int,
    fault_mode: str = "2q",   # '1q' or '2q' or 'either'
    fault_kind = None,            # None | 'X'|'Z'|'Y' | (kc,kt)
    w_min: int = 2
):
    """
    Inject ONE fault after gate k with the given fault_mode/kind.
    Return (ok, info) where ok=True if UNSAT (i.e., weight ≥ w_min for all assignments).
    """
    state, qc, site, groups = build_state_with_fault_after_gate(
        qasm_path, gate_index, fault_mode=fault_mode, fault_kind=fault_kind
    )
    b = data_error_weight_literals(state, groups["data"])

    s = Solver()
    # Enforce the chosen fault structure if it's a CNOT site
    add_fault_mode_constraints(s, site, fault_mode, fault_kind)

    # Look for a counterexample: data-weight ≤ w_min-1
    s.add(PbLe([(bi,1) for bi in b], w_min-1))

    if s.check() == unsat:
        return True, {"gate": qc.data[gate_index][0].name, "index": gate_index}
    mdl = s.model()
    # Report which wires (and which Pauli) got chosen under the model
    if site["gate_name"] == "cx":
        v = site["vars"]
        ctrl = (
            "Y" if (bool(mdl.eval(v["fxc"], True)) and bool(mdl.eval(v["fzc"], True)))
            else ("X" if bool(mdl.eval(v["fxc"], True))
            else ("Z" if bool(mdl.eval(v["fzc"], True)) else "I"))
        )
        targ = (
            "Y" if (bool(mdl.eval(v["fxt"], True)) and bool(mdl.eval(v["fzt"], True)))
            else ("X" if bool(mdl.eval(v["fxt"], True))
            else ("Z" if bool(mdl.eval(v["fzt"], True)) else "I"))
        )
        where = f"cx c={site['qubits'][0]}, t={site['qubits'][1]}  (ctrl={ctrl}, targ={targ})"
    else:
        v = site["vars"]
        k = "Y" if (bool(mdl.eval(v["fx"], True)) and bool(mdl.eval(v["fz"], True))) \
            else ("X" if bool(mdl.eval(v["fx"], True)) else ("Z" if bool(mdl.eval(v["fz"], True)) else "I"))
        where = f"{site['gate_name']} q={site['qubits'][0]}  ({k})"

    bvals = [bool(mdl.eval(e, True)) for e in b]
    return False, {
        "gate": qc.data[gate_index][0].name,
        "index": gate_index,
        "fault": where,
        "data_weight": sum(bvals),
        "data_bits_true": [groups["data"][i] for i,v in enumerate(bvals) if v],
    }
from z3 import Solver, ForAll, Exists, Or, Xor, PbLe

def forall_fault_exists_low_weight_per_gate(
    qasm_path: str,
    stab_txt_path: str,
    gate_indices=None,          # e.g. range(10) or [0,1,2]
    fault_mode: str = "2q",     # "2q" | "1q" | "either" (whatever your builder accepts)
    flag_axis: str = "z",       # inject flag error on this axis ("x" or "z")
    flag_prefix: str = "flagErr",
):
    """
    This is for checking 'bad loocation'
    For each gate in `gate_indices`:
      state, qc, site_info, groups = build_state_with_fault_after_gate(...)
      fault_vars = all 'f*' vars from site_info['vars'] (universally quantified)
      E' = stabilizer-equivalent data error
      b  = per-qubit error indicators
      Check:  ∀ fault_vars. ( Or(fault_vars) → ∃ gsel.  sum(b) ≤ 1 )
      And also: Xor(site_info['act'], flag_var)   (your extra constraint)

    Returns: dict {gate_index: {"result": sat/unsat, "num_fault_vars": int, "num_gens": int}}
    """
    results = {}
    unsat_gates = []
    # If user didn't pass indices, default to all gates
    if gate_indices is None:
        # Peek the circuit once to know how many gates
        _, qc, _, _ = build_state_with_fault_after_gate(qasm_path, gate_index=0, fault_mode=fault_mode)
        gate_indices = range(len(qc.data))

    for i in gate_indices:
        # 1) Build state with *symbolic* fault inserted after gate i
        state, qc, site_info, groups = build_state_with_fault_after_gate(
            qasm_path, gate_index=i, fault_mode=fault_mode
        )

        # 2) Collect the fault variables at this site (universally quantified)
        fault_vars = [v for k, v in site_info["vars"].items() if k.startswith("f")]
        if not fault_vars:
            # No fault DOFs at this site (e.g. a barrier/measure) → skip
            results[i] = {"result": "no-fault-vars", "num_fault_vars": 0, "num_gens": 0}
            continue

        # 3) Extract data error (E_x, E_z)
        data_idxs = groups["data"]
        E_x = [state.qubits[j].x for j in data_idxs]
        E_z = [state.qubits[j].z for j in data_idxs]

        

        # 4) Build stabilizer-equivalent errors E' using selector Booleans gsel
        Epx, Epz, gsel = build_stab_equiv_errors(E_x, E_z, stab_txt_path, prefix=f"g")

        # 5) Weight ≤ 1 predicate: sum over per-qubit indicators b_i = Or(E′x_i, E′z_i)
        b = [Or(xi, zi) for xi, zi in zip(Epx, Epz)]

        # 6) ∀ fault_vars: Or(fault_vars) → ∃ gsel: sum(b) ≤ 1
        s = Solver()
        body = Exists(gsel, PbLe([(bi, 1) for bi in b], 1))
        phi  = ForAll(fault_vars, Or(Or(fault_vars), False) == False)  # placeholder replaced below

        # Rebuild phi cleanly (the line above avoids z3py “no quantifier vars” edge cases if empty)
        phi = ForAll(fault_vars, 
                     Or(  # (¬any_fault) ∨ (∃ gsel: weight ≤ 1)
                        Or([v for v in fault_vars]) == False,
                        body
                     ))

        # 7) Add both the quantified property and your extra XOR constraint
        s.add(phi)
        

        res = s.check()
        results[i] = {
            "result": str(res),
            "num_fault_vars": len(fault_vars),
            "num_gens": len(gsel),
        }
        if str(res) == "unsat":
            unsat_gates.append(i)

    return results, unsat_gates
# ---------------------------
# Example CLI usage (optional)
# ---------------------------
if __name__ == "__main__":
    qasm_file = "my_flagged_round.qasm"

    # Suppose your layout has:
    # - ancillas at indices [8, 9, 10]  (three stabilizer measurements)
    # - flags at    indices [11, 12]    (two shared/parallel flags)
    anc_idxs  = [8, 9, 10]
    flag_idxs = [11, 12]

    # Inject a Z fault on flag 11 and a Y fault on flag 12 just before measurement
    errors = {11: "Z", 12: "Y"}

    syn, flg = analyze_flag_errors_multi(qasm_file, anc_idxs, flag_idxs, errors)

    print("Syndrome flips (per ancilla):")
    for i, expr in zip(anc_idxs, syn):
        print(f"  anc[{i}] Z -> {expr}")

    print("Flag X-basis flips (per flag):")
    for i, expr in zip(flag_idxs, flg):
        print(f"  flag[{i}] X-meas flip -> {expr}")