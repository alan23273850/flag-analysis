"""
Global coordinated learner for path_x (x from argv, 2<=x<=5): 10 learners (dataX0..dataZ4),
generalized syndrome = all round files from round_0, round_1, ... (in order), constraint condition.blif.

Usage: python decoder_learner_[[5,1,3]]_origin.py <x>   where 2 <= x <= 5

Final EQ formula (target to satisfy = find counterexample):
  1. PathCondition(E) == True   -- restrict error E (e.g. weight), path branch conditions
  2. C = L(Syndrome(E))         -- decoder output C (fresh) from syndrome of E; C_xi, C_zi = learner
  3. F_xi = E_xi ⊕ C_xi, F_zi = E_zi ⊕ C_zi  (residual = fix)
  4. For each line in [[5,1,3]]_origin.txt and [[5,1,3]]_log_op.txt (two 5-bit strings):
     Line parity = XOR_i [ (F_zi & s_z[i]) ⊕ (F_xi & s_x[i]) ]
     "Fixed" per line: Line parity == 0. Truly fixed: ALL lines hold.
  5. Target (logical error): at least ONE line has Line parity != 0
  + not_in_table. (No decoder-mismatch term.)
"""
from __future__ import annotations
import sys
import subprocess
import tempfile
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from blif_utils import (
    parse_blif,
    rename_blif,
    write_blif,
    combine_blifs,
    evaluate_blif,
    topological_sort_gates,
)

try:
    from z3 import (
        Bool,
        BoolVal,
        And,
        Or,
        Not,
        Xor,
        parse_smt2_string,
        substitute,
        simplify,
        Solver,
        sat,
        unsat,
        is_true,
        is_false,
    )
    _Z3_AVAILABLE = True
except ImportError:
    _Z3_AVAILABLE = False

# --- Config from argv ---
_PROJECT_ROOT = Path(__file__).resolve().parent


def _parse_path_arg() -> int:
    """Parse argv for path index x; 2 <= x <= 5."""
    if len(sys.argv) < 2:
        print("Usage: python decoder_learner_[[5,1,3]]_origin.py <x>   where 2 <= x <= 5", file=sys.stderr)
        sys.exit(1)
    try:
        x = int(sys.argv[1])
    except ValueError:
        print("Error: x must be an integer", file=sys.stderr)
        sys.exit(1)
    if not (2 <= x <= 5):
        print("Error: x must be in [2, 5]", file=sys.stderr)
        sys.exit(1)
    return x


def _discover_bi_paths(path_dir: Path) -> Tuple[List[Path], List[Path], List[str]]:
    """
    Scan path_dir/round_0, round_1, ... for *.blif; collect all as generalized syndrome.
    Returns (BI_BLIF_PATHS, BI_SMT2_PATHS, BI_LABELS).
    """
    blif_paths: List[Path] = []
    smt2_paths: List[Path] = []
    labels: List[str] = []
    r = 0
    while True:
        round_dir = path_dir / f"round_{r}"
        if not round_dir.is_dir():
            break
        # List *.blif in round, sorted by name
        blifs = sorted(round_dir.glob("*.blif"))
        for bp in blifs:
            name = bp.stem
            sp = round_dir / f"{name}.smt2"
            if sp.exists():
                blif_paths.append(bp)
                smt2_paths.append(sp)
                labels.append(f"{name}_r{r}")
        r += 1
    return blif_paths, smt2_paths, labels


def _init_config(x: int) -> None:
    """Set module-level config from path index x."""
    global PATH_DIR, BI_BLIF_PATHS, BI_SMT2_PATHS, BI_LABELS, NUM_BI, NUM_VARS
    global DATA_BLIF_PATHS, DATA_SMT2_PATHS, CONDITION_BLIF, CONDITION_SMT2

    PATH_DIR = _PROJECT_ROOT / "output_[[5,1,3]]_origin" / f"path_{x}"
    if not PATH_DIR.is_dir():
        print(f"Error: {PATH_DIR} does not exist", file=sys.stderr)
        sys.exit(1)

    BI_BLIF_PATHS, BI_SMT2_PATHS, BI_LABELS = _discover_bi_paths(PATH_DIR)
    if not BI_BLIF_PATHS:
        print("Error: no round BLIF files found", file=sys.stderr)
        sys.exit(1)
    NUM_BI = len(BI_BLIF_PATHS)

    # var_names.txt line count
    var_file = PATH_DIR / "var_names.txt"
    NUM_VARS = 0
    if var_file.exists():
        NUM_VARS = len(var_file.read_text().strip().splitlines())

    # 10 data formulas (fixed)
    DATA_BLIF_PATHS = [
        PATH_DIR / "dataX0.blif", PATH_DIR / "dataX1.blif", PATH_DIR / "dataX2.blif",
        PATH_DIR / "dataX3.blif", PATH_DIR / "dataX4.blif",
        PATH_DIR / "dataZ0.blif", PATH_DIR / "dataZ1.blif", PATH_DIR / "dataZ2.blif",
        PATH_DIR / "dataZ3.blif", PATH_DIR / "dataZ4.blif",
    ]
    DATA_SMT2_PATHS = [p.with_suffix(".smt2") for p in DATA_BLIF_PATHS]
    CONDITION_BLIF = PATH_DIR / "condition.blif"
    CONDITION_SMT2 = PATH_DIR / "condition.smt2"


# Initialize on import (after argv is available)
_PATH_X = _parse_path_arg()
_init_config(_PATH_X)

CODE_ORIGIN_DIR = _PROJECT_ROOT / "[[5,1,3]]_origin"
ORIGIN_TXT = CODE_ORIGIN_DIR / "[[5,1,3]]_origin.txt"
LOG_OP_TXT = CODE_ORIGIN_DIR / "[[5,1,3]]_log_op.txt"
_ABC_DIR = _PROJECT_ROOT / "abc"
ABC_BIN = _ABC_DIR / "abc" if (_ABC_DIR / "abc").exists() else Path("abc")

NUM_DATA = 10
NUM_DATA_QUBITS = 5


def _xor_chain(exprs: List[object]) -> object:
    """XOR of all expressions (parity). Empty -> False."""
    if not exprs:
        return BoolVal(False)
    acc = exprs[0]
    for e in exprs[1:]:
        acc = Xor(acc, e)
    return acc


def _load_stab_log_lines() -> List[Tuple[List[int], List[int]]]:
    """Load all lines from [[5,1,3]]_origin.txt and [[5,1,3]]_log_op.txt."""
    result: List[Tuple[List[int], List[int]]] = []
    for path in (ORIGIN_TXT, LOG_OP_TXT):
        if not path.exists():
            continue
        for line in path.read_text(encoding="utf-8").strip().splitlines():
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2 or len(parts[0]) != 5 or len(parts[1]) != 5:
                continue
            Sx = [1 if c == "1" else 0 for c in parts[0]]
            Sz = [1 if c == "1" else 0 for c in parts[1]]
            result.append((Sx, Sz))
    return result


def _build_logical_error_constraint(
    data_forms: List[object],
    h_formulas: List[object],
) -> Optional[object]:
    """Build Z3 formula: (E⊕C) is logical error per stabilizer/logical-op lines."""
    if not _Z3_AVAILABLE or len(data_forms) != NUM_DATA or len(h_formulas) != NUM_DATA:
        return None
    lines = _load_stab_log_lines()
    if not lines:
        return None
    F_x = [Xor(data_forms[i], h_formulas[i]) for i in range(NUM_DATA_QUBITS)]
    F_z = [Xor(data_forms[NUM_DATA_QUBITS + i], h_formulas[NUM_DATA_QUBITS + i]) for i in range(NUM_DATA_QUBITS)]
    line_fail_formulas: List[object] = []
    for Sx, Sz in lines:
        terms: List[object] = []
        for i in range(NUM_DATA_QUBITS):
            t_x = And(F_x[i], BoolVal(True)) if Sz[i] else BoolVal(False)
            t_z = And(F_z[i], BoolVal(True)) if Sx[i] else BoolVal(False)
            terms.append(t_x)
            terms.append(t_z)
        line_parity = _xor_chain(terms)
        line_fail_formulas.append(line_parity)
    return Or(*line_fail_formulas)


def _run_abc_sat(blif_path: Path, timeout_sec: int = 60) -> Optional[Dict[str, int]]:
    """Run ABC on BLIF; return assignment dict for primary inputs if SAT, else None."""
    import shutil
    with tempfile.NamedTemporaryFile(mode="w", suffix=".blif", delete=False) as f:
        pass
    out_blif = Path(f.name)
    out_cex = out_blif.with_suffix(".cex")
    try:
        shutil.copy(blif_path, out_blif)
        work = out_blif.parent
        cex_name = out_cex.name
        cmd = [str(ABC_BIN), "-c", f"read_blif {out_blif.name}; sat; write_cex {cex_name}"]
        r = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout_sec, cwd=str(work)
        )
        if "UNSATISFIABLE" in (r.stdout or "") or "UNSATISFIABLE" in (r.stderr or ""):
            return None
        if "SATISFIABLE" not in (r.stdout or "") and "SATISFIABLE" not in (r.stderr or ""):
            return None
        cex_path = work / cex_name
        if not cex_path.exists():
            return None
        lines = cex_path.read_text().splitlines()
        if not lines:
            return None
        bits = []
        for c in lines[0]:
            if c in "01":
                bits.append(int(c))
            elif c == "#":
                break
        inp, _, _ = parse_blif(blif_path)
        if len(bits) < len(inp):
            return None
        return {inp[i]: bits[i] for i in range(len(inp))}
    finally:
        out_blif.unlink(missing_ok=True)
        if out_cex.exists():
            out_cex.unlink(missing_ok=True)


def _build_mq_blif(bi: Tuple[int, ...], work_dir: Path) -> Path:
    """Build BLIF: condition & (b0==BI[0]) & ... & (bN==BI[N]), output 'valid'."""
    assert len(bi) == NUM_BI
    c_inp, c_out, c_gates = parse_blif(CONDITION_BLIF)
    c_inp, c_out, c_gates = rename_blif(c_inp, c_out, c_gates, "c0_", keep_inputs=True)
    shared = c_inp

    specs: List[Tuple[List[str], List[str], List]] = [(c_inp, c_out, c_gates)]
    prefixes = ["c0_"]
    for i, bp in enumerate(BI_BLIF_PATHS):
        inp, out, gates = parse_blif(bp)
        inp, out, gates = rename_blif(inp, out, gates, f"b{i}_", keep_inputs=True)
        specs.append((inp, out, gates))
        prefixes.append(f"b{i}_")
    all_inputs, all_outputs, all_gates = combine_blifs(specs, prefixes, shared)

    wire_c0 = "c0_out"
    wire_b = [f"b{i}_out" for i in range(NUM_BI)]
    for i in range(NUM_BI):
        const_name = f"const{i}"
        eq_name = f"eq{i}"
        all_gates.append(([], const_name, [["1"]] if bi[i] else [["0"]]))
        if bi[i]:
            all_gates.append(([wire_b[i], const_name], eq_name, [["11", "1"]]))
        else:
            all_gates.append(([wire_b[i], const_name], eq_name, [["00", "1"]]))
    and_inputs = [wire_c0] + [f"eq{i}" for i in range(NUM_BI)]
    all_gates.append((and_inputs, "valid", [["1" * len(and_inputs), "1"]]))

    out_path = work_dir / "mq.blif"
    write_blif(out_path, all_inputs, ["valid"], all_gates)
    return out_path


def _evaluate_condition0(assignment: Dict[str, int]) -> int:
    inp, _, gates = parse_blif(CONDITION_BLIF)
    val = evaluate_blif(inp, gates, assignment)
    return val.get("out", 0)


def _evaluate_bi(assignment: Dict[str, int]) -> Tuple[int, ...]:
    result = []
    for bp in BI_BLIF_PATHS:
        inp, _, gates = parse_blif(bp)
        val = evaluate_blif(inp, gates, assignment)
        result.append(val.get("out", 0))
    return tuple(result)


def _evaluate_data(assignment: Dict[str, int]) -> Tuple[int, ...]:
    result = []
    for dp in DATA_BLIF_PATHS:
        inp, _, gates = parse_blif(dp)
        val = evaluate_blif(inp, gates, assignment)
        result.append(val.get("out", 0))
    return tuple(result)


class PathXGlobalCoordinator:
    """Global table and query handlers for path_x. BI = N-tuple of 0/1. Data = 10-tuple."""

    def __init__(self, work_dir: Optional[Path] = None):
        self.table: Dict[Tuple[int, ...], Tuple[int, ...]] = {}
        self.work_dir = work_dir or Path(tempfile.mkdtemp())

    def membership_query(self, learner_idx: int, bi: Tuple[int, ...]) -> Optional[int]:
        if len(bi) != NUM_BI:
            raise ValueError(f"BI must be length {NUM_BI}, got {len(bi)}")
        bi = tuple(int(b) for b in bi)
        if bi in self.table:
            return self.table[bi][learner_idx]
        phi = _build_mq_z3(bi)
        if phi is None:
            raise RuntimeError("MQ: Z3 unified formulas not available")
        print("Running Z3 SAT for MQ")
        assignment = _run_z3_sat(phi)
        if assignment is None:
            phi_any = _build_any_valid_syndrome_z3()
            if phi_any is None:
                return None
            print("MQ: requested BI UNSAT, trying any valid syndrome")
            assignment = _run_z3_sat(phi_any)
            if assignment is None:
                return None
            bi_forms = [_load_unified_formula(p) for p in BI_SMT2_PATHS]
            data_forms = [_load_unified_formula(p) for p in DATA_SMT2_PATHS]
            bi_vals = [_eval_unified_formula(f, assignment) for f in bi_forms]
            data_vals_list = [_eval_unified_formula(f, assignment) for f in data_forms]
            if None in bi_vals or None in data_vals_list:
                return None
            bi_found = tuple(bi_vals)
            data_vals = tuple(data_vals_list)
            self.table[bi_found] = data_vals
            return data_vals[learner_idx]
        data_forms = [_load_unified_formula(p) for p in DATA_SMT2_PATHS]
        data_vals_list = [_eval_unified_formula(f, assignment) for f in data_forms]
        if None in data_vals_list:
            raise RuntimeError("MQ: Z3 unified formula eval failed for data")
        data_vals = tuple(data_vals_list)
        self.table[bi] = data_vals
        return data_vals[learner_idx]

    def _build_eq_z3(
        self,
        conjectures: List[Callable[[Tuple[int, ...]], int]],
        use_table_block: bool = True,
    ) -> Optional[object]:
        cond = _load_unified_formula(CONDITION_SMT2)
        bi_forms = [_load_unified_formula(p) for p in BI_SMT2_PATHS]
        data_forms = [_load_unified_formula(p) for p in DATA_SMT2_PATHS]
        if cond is None or any(f is None for f in bi_forms) or any(f is None for f in data_forms):
            return None
        h_formulas = []
        for i in range(NUM_DATA):
            ones = [bi for bi, row in self.table.items() if row[i] == 1]
            if not ones:
                h_formulas.append(BoolVal(False))
            else:
                terms = [And(*(bi_forms[j] if r[j] else Not(bi_forms[j]) for j in range(NUM_BI))) for r in ones]
                h_formulas.append(Or(*terms))
        if use_table_block and self.table:
            not_in_table = And([
                Not(And(*(bi_forms[j] if r[j] else Not(bi_forms[j]) for j in range(NUM_BI))))
                for r in self.table
            ])
        else:
            not_in_table = BoolVal(True)
        logical_error = _build_logical_error_constraint(data_forms, h_formulas)
        if logical_error is not None:
            return And(cond, not_in_table, logical_error)
        return And(cond, not_in_table)

    def equivalence_query(
        self,
        conjectures: List[Callable[[Tuple[int, ...]], int]],
        timeout_sec: int = 120,
    ) -> Optional[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
        if len(conjectures) != NUM_DATA:
            raise ValueError(f"Need {NUM_DATA} conjectures, got {len(conjectures)}")
        phi = self._build_eq_z3(conjectures, use_table_block=True)
        if phi is None:
            return None
        print("Running Z3 SAT for EQ")
        assignment = _run_z3_sat(phi, timeout_sec=timeout_sec)
        if assignment is None:
            return None
        bi_forms = [_load_unified_formula(p) for p in BI_SMT2_PATHS]
        data_forms = [_load_unified_formula(p) for p in DATA_SMT2_PATHS]
        bi_cex_list = [_eval_unified_formula(f, assignment) for f in bi_forms]
        data_vals_list = [_eval_unified_formula(f, assignment) for f in data_forms]
        if None in bi_cex_list or None in data_vals_list:
            raise RuntimeError("EQ: Z3 unified formula eval failed for counterexample")
        bi_cex = tuple(bi_cex_list)
        data_vals = tuple(data_vals_list)
        data_eval = data_vals
        if bi_cex in self.table:
            return None
        self.table[bi_cex] = data_vals
        for i in range(NUM_DATA):
            if self.table[bi_cex][i] != data_eval[i]:
                raise AssertionError(
                    f"Learner {i}: learned != formula output for BI={bi_cex}"
                )
        print_verification_detail(assignment, bi_cex, data_vals)
        return (bi_cex, data_vals)

    def _get_bi_from_eq_model(self, assignment: Dict[str, int], eq_blif_path: Path) -> Optional[Tuple[int, ...]]:
        inp, _, gates = parse_blif(eq_blif_path)
        val = evaluate_blif(inp, gates, assignment)
        bi = []
        for i in range(NUM_BI):
            w = f"b{i}_out"
            bi.append(val.get(w, 0))
        return tuple(bi)

    def _build_eq_blif(
        self,
        conjectures: List[Callable[[Tuple[int, ...]], int]],
        use_table_block: bool = True,
    ) -> Optional[Path]:
        c_inp, c_out, c_gates = parse_blif(CONDITION_BLIF)
        shared = c_inp
        specs: List[Tuple[List[str], List[str], List]] = [(c_inp, c_out, c_gates)]
        prefixes = ["c0_"]
        for i, bp in enumerate(BI_BLIF_PATHS):
            inp, out, gates = parse_blif(bp)
            specs.append((inp, out, gates))
            prefixes.append(f"b{i}_")
        for i, dp in enumerate(DATA_BLIF_PATHS):
            inp, out, gates = parse_blif(dp)
            specs.append((inp, out, gates))
            prefixes.append(f"d{i}_")
        all_inputs, all_outputs, all_gates = combine_blifs(specs, prefixes, shared)
        wire_b = [f"b{i}_out" for i in range(NUM_BI)]
        wire_d = [f"d{i}_out" for i in range(NUM_DATA)]

        bi_range = 2 ** NUM_BI
        for i in range(NUM_DATA):
            ones = []
            for bits in range(bi_range):
                bi_tuple = tuple((bits >> j) & 1 for j in range(NUM_BI))
                if conjectures[i](bi_tuple) == 1:
                    ones.append("".join(str((bits >> j) & 1) for j in range(NUM_BI)))
            hi_out = f"h{i}"
            if ones:
                all_gates.append((wire_b, hi_out, [[pat, "1"] for pat in ones]))
            else:
                all_gates.append(([], hi_out, [["0"]]))

        for i in range(NUM_DATA):
            all_gates.append(([f"h{i}", wire_d[i]], f"xor{i}", [["10", "1"], ["01", "1"]]))
        mismatch_rows = [
            ["".join("1" if k == i else "0" for k in range(NUM_DATA)), "1"]
            for i in range(NUM_DATA)
        ]
        all_gates.append(([f"xor{i}" for i in range(NUM_DATA)], "mismatch", mismatch_rows))

        not_in_table_inputs: List[str] = []
        if not use_table_block or not self.table:
            all_gates.append(([], "one", [["1"]]))
            not_in_table_inputs.append("one")
        else:
            for idx, row in enumerate(self.table):
                eq_name = f"eq_row_{idx}"
                pat = "".join(str(t) for t in row)
                all_gates.append((wire_b, eq_name, [[pat, "1"]]))
                neq_name = f"neq_{idx}"
                all_gates.append(([eq_name], neq_name, [["0", "1"]]))
                not_in_table_inputs.append(neq_name)
        not_in_table_val = not_in_table_inputs[0] if len(not_in_table_inputs) == 1 else "not_in_table"
        if len(not_in_table_inputs) > 1:
            all_gates.append((
                not_in_table_inputs,
                not_in_table_val,
                [["1" * len(not_in_table_inputs), "1"]],
            ))
        miter_inputs = ["c0_out", "mismatch", not_in_table_val]
        all_gates.append((miter_inputs, "miter", [["1" * len(miter_inputs), "1"]]))
        all_gates = topological_sort_gates(all_inputs, all_gates)
        out_path = self.work_dir / "eq.blif"
        write_blif(out_path, all_inputs, ["miter"], all_gates)
        return out_path

    def verify_learned_formulas(self, default_unseen: int = 0) -> Tuple[bool, Optional[Dict[str, int]]]:
        phi = self._build_eq_z3(
            [lambda bi, i=j: self.table[bi][i] if bi in self.table else default_unseen for j in range(NUM_DATA)],
            use_table_block=False,
        )
        if phi is None:
            return False, None
        assignment = _run_z3_sat(phi, timeout_sec=120)
        if assignment is None:
            return True, None
        return False, assignment


DATA_LABELS = [
    "dataX0", "dataX1", "dataX2", "dataX3", "dataX4",
    "dataZ0", "dataZ1", "dataZ2", "dataZ3", "dataZ4",
]


def _load_var_original_names(path_dir: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    var_names_file = path_dir / "var_names.txt"
    if not var_names_file.exists():
        return out
    try:
        lines = var_names_file.read_text(encoding="utf-8").strip().splitlines()
        for i, line in enumerate(lines):
            out[f"v{i + 1}"] = line.strip() or f"v{i + 1}"
    except OSError:
        pass
    return out


_z3_formula_cache: Dict[Path, Optional[object]] = {}
_unified_symbols: List[object] = []
_unified_by_name: Dict[str, object] = {}
_unified_formula_cache: Dict[Path, Optional[object]] = {}


def _ensure_unified_symbols() -> bool:
    global _unified_symbols, _unified_by_name
    if _unified_symbols:
        return True
    if not _Z3_AVAILABLE:
        return False
    path = PATH_DIR / "var_names.txt"
    if not path.exists():
        return False
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    _unified_symbols = []
    _unified_by_name = {}
    for i, line in enumerate(lines):
        name = line.strip() or f"v{i + 1}"
        sym = Bool(name)
        _unified_symbols.append(sym)
        _unified_by_name[name] = sym
    return True


def _load_unified_formula(smt2_path: Path):
    if not _Z3_AVAILABLE or not _ensure_unified_symbols():
        return None
    if smt2_path in _unified_formula_cache:
        return _unified_formula_cache[smt2_path]
    raw = _load_z3_formula(smt2_path)
    if raw is None:
        _unified_formula_cache[smt2_path] = None
        return None
    consts = _get_bool_constants_z3(raw)
    subst = [(c, _unified_by_name[name]) for c in consts
              if (name := c.decl().name()) in _unified_by_name]
    if not subst:
        _unified_formula_cache[smt2_path] = raw
        return raw
    unified = substitute(raw, subst)
    _unified_formula_cache[smt2_path] = unified
    return unified


def _run_z3_sat(phi, timeout_sec: int = 120) -> Optional[Dict[str, int]]:
    if not _Z3_AVAILABLE or phi is None or not _ensure_unified_symbols():
        return None
    s = Solver()
    s.set("timeout", timeout_sec * 1000)
    s.add(phi)
    if s.check() != sat:
        return None
    model = s.model()
    return {
        f"v{i + 1}": 1 if is_true(model.eval(sym)) else 0
        for i, sym in enumerate(_unified_symbols)
    }


def _eval_unified_formula(phi, assignment: Dict[str, int]) -> Optional[int]:
    if not _Z3_AVAILABLE or phi is None or not _ensure_unified_symbols():
        return None
    try:
        subst = [
            (_unified_symbols[i], BoolVal(bool(assignment.get(f"v{i + 1}", 0))))
            for i in range(len(_unified_symbols))
        ]
        sub = substitute(phi, subst)
        simplified = simplify(sub)
        if is_true(simplified):
            return 1
        if is_false(simplified):
            return 0
        return None
    except Exception:
        return None


def _build_mq_z3(bi: Tuple[int, ...]) -> Optional[object]:
    cond = _load_unified_formula(CONDITION_SMT2)
    if cond is None:
        return None
    bi_forms = [_load_unified_formula(p) for p in BI_SMT2_PATHS]
    if any(f is None for f in bi_forms):
        return None
    clauses = [cond]
    for i in range(NUM_BI):
        clauses.append(bi_forms[i] if bi[i] else Not(bi_forms[i]))
    return And(*clauses)


def _build_any_valid_syndrome_z3() -> Optional[object]:
    return _load_unified_formula(CONDITION_SMT2)


def _get_bool_constants_z3(expr) -> List:
    out: List = []
    seen = set()
    def visit(e):
        if id(e) in seen:
            return
        seen.add(id(e))
        try:
            if e.num_args() == 0 and str(e.sort()) == "Bool":
                out.append(e)
            for c in e.children():
                visit(c)
        except Exception:
            pass
    visit(expr)
    return out


def _load_z3_formula(smt2_path: Path):
    if smt2_path in _z3_formula_cache:
        return _z3_formula_cache[smt2_path]
    if not _Z3_AVAILABLE or not smt2_path.exists():
        _z3_formula_cache[smt2_path] = None
        return None
    try:
        content = smt2_path.read_text(encoding="utf-8")
        parsed = parse_smt2_string(content)
        assertions = parsed[1] if isinstance(parsed, tuple) and len(parsed) > 1 else parsed
        formula = assertions[0] if assertions else None
        _z3_formula_cache[smt2_path] = formula
        return formula
    except Exception:
        _z3_formula_cache[smt2_path] = None
        return None


def _print_sat_assignment(assignment: Dict[str, int]) -> None:
    var_display = _load_var_original_names(PATH_DIR)
    def vkey(x: str) -> Tuple[int, int]:
        return (0, int(x[1:])) if x.startswith("v") and x[1:].isdigit() else (1, 0)
    fault_names = sorted(assignment.keys(), key=vkey)
    print("  SAT assignment (boolean values):")
    for name in fault_names:
        print(f"    {var_display.get(name, name)} = {assignment[name]}  ({'true' if assignment[name] else 'false'})")
    bits = "".join(str(assignment.get(n, "?")) for n in fault_names)
    print(f"  Compact (v1..v{len(fault_names)}): {bits}")


def print_verification_detail(
    assignment: Dict[str, int],
    bi_eval: Tuple[int, ...],
    data_eval: Tuple[int, ...],
) -> None:
    var_display = _load_var_original_names(PATH_DIR)
    fault_names = sorted(assignment.keys(), key=lambda x: (len(x), x))
    print("    Fault variables:")
    for name in fault_names:
        print(f"      {var_display.get(name, name)} = {'true' if assignment[name] else 'false'}")
    print("    BI formulas:", NUM_BI)
    for i in range(NUM_BI):
        lab = BI_LABELS[i] if i < len(BI_LABELS) else f"b{i}"
        print(f"      {lab} => {bi_eval[i]}  ({'true' if bi_eval[i] else 'false'})")
    print("    Data formulas:", NUM_DATA)
    for i in range(NUM_DATA):
        lab = DATA_LABELS[i] if i < len(DATA_LABELS) else f"d{i}"
        print(f"      {lab} => {data_eval[i]}  ({'true' if data_eval[i] else 'false'})")


def learner_formula_dnf(bi: Tuple[int, ...], labels: List[str]) -> str:
    terms = []
    for j, b in enumerate(bi):
        lab = labels[j] if j < len(labels) else f"b{j}"
        terms.append(lab if b else f"¬{lab}")
    return "(" + " ∧ ".join(terms) + ")"


def print_learner_formulas(coord: PathXGlobalCoordinator) -> None:
    for i in range(NUM_DATA):
        ones = [bi for bi, row in coord.table.items() if row[i] == 1]
        formula = "0" if not ones else " ∨ ".join(learner_formula_dnf(bi, BI_LABELS) for bi in ones)
        print(f"  {DATA_LABELS[i]}: {formula}")


def run_learning_loop(coord: PathXGlobalCoordinator, default_unseen: int = 0) -> int:
    round_count = 0
    while True:
        def make_conjecture(i: int):
            def H_i(bi: Tuple[int, ...]) -> int:
                return coord.table[bi][i] if bi in coord.table else default_unseen
            return H_i
        conjectures = [make_conjecture(i) for i in range(NUM_DATA)]
        result = coord.equivalence_query(conjectures)
        if result is None:
            break
        bi_cex, data_vals = result
        print(f"  EQ round {round_count}: cex BI={bi_cex} data={data_vals}")
        round_count += 1
    return round_count


if __name__ == "__main__":
    import io
    import sys
    output_path = _PROJECT_ROOT / "output_[[5,1,3]]_origin" / f"output_{_PATH_X}.txt"
    buf = io.StringIO()

    class TeeWriter:
        def __init__(self, stream, buf):
            self.stream = stream
            self.buf = buf
        def write(self, s):
            self.stream.write(s)
            self.buf.write(s)
        def flush(self):
            self.stream.flush()
            self.buf.flush()

    orig_stdout = sys.stdout
    sys.stdout = TeeWriter(orig_stdout, buf)
    try:
        print(f"Path{_PATH_X} global coordinator – learning loop until equivalent (BI={NUM_BI}, vars={NUM_VARS})")
        coord = PathXGlobalCoordinator()
        n_rounds = run_learning_loop(coord, default_unseen=0)
        print(f"Equivalent after {n_rounds} EQ rounds. Table size: {len(coord.table)}")
        print("\nLearner boolean formulas (DNF over BI):")
        print_learner_formulas(coord)
        print("\nVerification: PathCondition(E) & (C=L(Syn(E))) & logical_error (stab/log) → UNSAT?")
        verified, sat_assignment = coord.verify_learned_formulas(default_unseen=0)
        if verified:
            print("  UNSAT: learned decoder consistent (path + logical-error condition).")
        else:
            print("  SAT: exists (v) satisfying path + logical-error parity.")
            if sat_assignment is not None:
                _print_sat_assignment(sat_assignment)
    finally:
        sys.stdout = orig_stdout
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(buf.getvalue(), encoding="utf-8")
    print(f"(Output also written to {output_path})", file=orig_stdout)
