"""
boolformula.py

Python data structure mirroring bull's boolformula_t (AND/OR/XOR/literal)
and conversion from Z3 BoolRef so formulas can be printed in the same
format as boolformula_print in C: { a | b } for OR, { a & b } for AND,
{ a ^ b } for XOR, and literals as signed integers.
"""

from __future__ import annotations
from typing import Union, List, Dict, Optional
from dataclasses import dataclass

# Z3: only needed when converting
try:
    from z3 import (
        BoolRef, BoolVal, Not, Xor as Z3Xor, is_true, is_false, is_and, is_or, is_not, is_const,
    )
    _Z3_AVAILABLE = True
except ImportError:
    _Z3_AVAILABLE = False


# ---------------------------
# BoolFormula tree (mirrors boolformula_t)
# ---------------------------

@dataclass(frozen=True)
class Literal:
    """Literal: signed variable id (positive = var, negative = not var)."""
    lit: int

@dataclass
class And:
    """Conjunction (AND of children)."""
    children: List["BoolFormula"]

@dataclass
class Or:
    """Disjunction (OR of children)."""
    children: List["BoolFormula"]

@dataclass
class Xor:
    """Exclusive disjunction (XOR of children)."""
    children: List["BoolFormula"]

BoolFormula = Union[Literal, And, Or, Xor]


def format_boolformula(f: BoolFormula) -> str:
    """
    Print in the same style as C boolformula_print:
    - Literal: signed integer
    - And: { ... & ... & ... }
    - Or:  { ... | ... | ... }
    - Xor: { ... ^ ... ^ ... }
    - Empty And -> { T }, empty Or -> { F }, empty Xor -> { 0 }
    """
    if isinstance(f, Literal):
        return str(f.lit)
    if isinstance(f, And):
        if not f.children:
            return "{ T }"
        parts = [format_boolformula(c) for c in f.children]
        return "{ " + " & ".join(parts) + " }"
    if isinstance(f, Or):
        if not f.children:
            return "{ F }"
        parts = [format_boolformula(c) for c in f.children]
        return "{ " + " | ".join(parts) + " }"
    if isinstance(f, Xor):
        if not f.children:
            return "{ 0 }"
        parts = [format_boolformula(c) for c in f.children]
        return "{ " + " ^ ".join(parts) + " }"
    raise TypeError(f"Not a BoolFormula: {type(f)}")


def to_extended_dimacs(bf: BoolFormula, num_vars: int) -> str:
    """
    Convert a BoolFormula to extended DIMACS (supporting 'x' lines for XOR).
    Assumes bf is a Conjunction of (Literal | Or | Xor) or just a single (Literal | Or | Xor).
    Raises ValueError if formula contains unsupported nesting (e.g. Or containing Xor).
    """
    clauses = []
    
    def get_lit(node: BoolFormula) -> int:
        if isinstance(node, Literal):
            return node.lit
        elif isinstance(node, Or) and not node.children:
            return 0  # { F } acts as identity in XOR / is ignored in OR
        raise ValueError(f"Cannot use {type(node)} as literal")

    def add_clause(c: BoolFormula):
        if isinstance(c, And):
            for child in c.children:
                add_clause(child)
        elif isinstance(c, Literal):
            clauses.append([c.lit])
        elif isinstance(c, Or):
            lits = []
            for child in c.children:
                if isinstance(child, And) and not child.children:
                    return # { T } inside OR means the whole clause is True, ignore it
                l = get_lit(child)
                if l != 0: lits.append(l)
            if not lits and c.children:
                pass # All children were { F } -> empty clause
            clauses.append(lits)
        elif isinstance(c, Xor):
            lits = []
            parity = 0  # 1 means we need to negate the XOR
            for child in c.children:
                if isinstance(child, And) and not child.children:
                    parity ^= 1
                    continue
                l = get_lit(child)
                if l != 0: lits.append(l)
            if parity == 1:
                if not lits:
                    # XOR is just { T } -> True clause -> ignore
                    return
                # Negate the first literal to flip the XOR parity
                lits[0] = -lits[0]
            if not lits:
                clauses.append(('x', [])) # empty XOR -> False
            else:
                clauses.append(('x', lits))
        else:
            raise ValueError(f"Unsupported clause type: {type(c)}")

    add_clause(bf)

    lines = [f"p cnf {num_vars} {len(clauses)}"]
    for cl in clauses:
        if isinstance(cl, tuple) and cl[0] == 'x':
            if not cl[1]:
                lines.append("0") # empty XOR is false
            else:
                lines.append("x " + " ".join(str(l) for l in cl[1]) + " 0")
        else:
            if not cl:
                lines.append("0")
            else:
                lines.append(" ".join(str(l) for l in cl) + " 0")

    return "\n".join(lines) + "\n"


def to_pure_cnf(bf: BoolFormula) -> List[List[int]]:
    """
    Recursively converts a BoolFormula into pure CNF (list of clauses)
    without introducing auxiliary variables. Note: this can cause exponential 
    blowup in size for highly nested formulas.
    """
    if isinstance(bf, Literal):
        return [[bf.lit]]
    elif isinstance(bf, And):
        clauses = []
        for c in bf.children:
            clauses.extend(to_pure_cnf(c))
        return clauses
    elif isinstance(bf, Or):
        if not bf.children:
            return [[]]  # False
        if len(bf.children) == 1:
            return to_pure_cnf(bf.children[0])
        import itertools
        child_cnfs = [to_pure_cnf(c) for c in bf.children]
        clauses = []
        # Protect against massive exponential blowup
        total_combos = 1
        for cnfs in child_cnfs:
            total_combos *= len(cnfs)
        if total_combos > 1000000:
            raise ValueError(f"Exponential blowup in pure CNF conversion (requires {total_combos} combinations). Tseitin transformation is required for this formula size, but it was disabled.")
            
        for combo in itertools.product(*child_cnfs):
            new_clause = []
            for cl in combo:
                new_clause.extend(cl)
            new_clause = list(set(new_clause))
            # If clause contains both A and -A, it's a tautology, we can skip it
            if not any(-lit in new_clause for lit in new_clause):
                clauses.append(new_clause)
        return clauses
    elif isinstance(bf, Xor):
        if not bf.children:
            return [[]]  # Empty XOR is False
        if len(bf.children) == 1:
            return to_pure_cnf(bf.children[0])
        lits = []
        parity = 0
        for c in bf.children:
            if isinstance(c, Literal):
                lits.append(c.lit)
            elif isinstance(c, And) and not c.children:
                parity ^= 1
            elif isinstance(c, Or) and not c.children:
                pass
            else:
                raise ValueError(f"Xor contains complex child which cannot be trivially flattened: {c}")
        
        if parity == 1:
            if not lits:
                return []  # True
            lits[0] = -lits[0]
            
        if len(lits) > 15:
            raise ValueError(f"XOR clause with {len(lits)} literals would produce 2^{len(lits)-1} clauses in pure CNF. Tseitin transformation is required for this size, but it was disabled.")
            
        import itertools
        clauses = []
        n = len(lits)
        for signs in itertools.product([1, -1], repeat=n):
            if sum(1 for s in signs if s == -1) % 2 == 0:
                clause = [s * l for s, l in zip(signs, lits)]
                clauses.append(clause)
        return clauses


def to_pure_dimacs(bf: BoolFormula, num_vars: int) -> str:
    """
    Convert a BoolFormula to pure DIMACS CNF (no 'x' lines, no auxiliary variables).
    """
    clauses = to_pure_cnf(bf)
    lines = [f"p cnf {num_vars} {len(clauses)}"]
    for cl in clauses:
        if not cl:
            lines.append("0")
        else:
            lines.append(" ".join(str(l) for l in cl) + " 0")
    return "\n".join(lines) + "\n"


def _is_bool_var(e: "BoolRef") -> bool:
    """True iff e is an uninterpreted boolean variable (not True/False, not compound)."""
    if not is_const(e):
        return False
    if is_true(e) or is_false(e):
        return False
    return True


def _flatten_xor(e: "BoolRef") -> List["BoolRef"]:
    """Flatten nested Xor(a, Xor(b, c)) into [a, b, c]."""
    if not _Z3_AVAILABLE or not isinstance(e, BoolRef):
        return [e]
    if e.decl().name() == "xor":
        out: List[BoolRef] = []
        for c in e.children():
            out.extend(_flatten_xor(c))
        return out
    return [e]


def z3_to_boolformula(
    expr: "BoolRef",
    var_to_id: Optional[Dict["BoolRef", int]] = None,
) -> BoolFormula:
    """
    Convert a Z3 BoolRef to our BoolFormula tree (boolformula_t style).

    - var_to_id: mapping from Z3 BoolRef variable to integer id (positive).
      If None, we collect all variables in expr and assign ids 1, 2, 3, ...
      in deterministic order (sorted by decl name).
    """
    if not _Z3_AVAILABLE:
        raise RuntimeError("z3 is required for z3_to_boolformula")

    if var_to_id is None:
        from dimacs_bridge import collect_bool_symbols
        syms = collect_bool_symbols(expr)
        var_to_id = {v: (i + 1) for i, v in enumerate(sorted(syms, key=lambda s: s.decl().name()))}

    def lit_for(e: BoolRef) -> int:
        if is_true(e) or is_false(e):
            raise ValueError("BoolVal True/False as literal not supported in boolformula_t style")
        if _is_bool_var(e):
            return var_to_id[e]
        if is_not(e) and _is_bool_var(e.children()[0]):
            return -var_to_id[e.children()[0]]
        raise ValueError(f"Not a literal: {e}")

    def go(e: BoolRef) -> BoolFormula:
        if is_true(e):
            return And([])   # { T }
        if is_false(e):
            return Or([])    # { F }
        if _is_bool_var(e):
            return Literal(lit_for(e))
        if is_not(e):
            sub = e.children()[0]
            if is_not(sub):
                return go(sub.children()[0])  # Not(Not(x)) = x
            if is_true(sub):
                return Or([]) # Not(True) = False
            if is_false(sub):
                return And([]) # Not(False) = True
            if _is_bool_var(sub):
                return Literal(lit_for(e))
            # Push Not into the tree (boolformula_t has no NOT node):
            # Not(And(a,b)) = Or(Not(a), Not(b)); Not(Or(a,b)) = And(Not(a), Not(b));
            # Not(Xor(a,b,c)) = Xor(Not(a), b, c).
            if is_and(sub):
                return Or([go(Not(c)) for c in sub.children()])
            if is_or(sub):
                return And([go(Not(c)) for c in sub.children()])
            if sub.decl().name() == "xor":
                flat = _flatten_xor(sub)
                if not flat:
                    return Or([])
                first_neg = go(Not(flat[0]))
                rest = [go(c) for c in flat[1:]]
                return Xor([first_neg] + rest)
            # Not(Iff(a,b)) = Xor(a,b)
            if sub.decl().name() == "=" and sub.num_args() == 2:
                return go(Z3Xor(sub.children()[0], sub.children()[1]))
            raise ValueError(f"Unhandled Z3 expression in Not: {e}")
        if is_and(e):
            return And([go(c) for c in e.children()])
        if is_or(e):
            return Or([go(c) for c in e.children()])
        if e.decl().name() == "xor":
            flat = _flatten_xor(e)
            return Xor([go(c) for c in flat])
        # Z3 Iff (a == b) for booleans: Iff(a,b) = Not(Xor(a,b)) = Xor(Not(a), b)
        if e.decl().name() == "=" and e.num_args() == 2:
            a, b = e.children()[0], e.children()[1]
            return go(Not(Z3Xor(a, b)))
        raise ValueError(f"Unhandled Z3 expression: {e}")

    return go(expr)


def z3_to_boolformula_str(
    expr: "BoolRef",
    var_to_id: Optional[Dict["BoolRef", int]] = None,
) -> str:
    """Convert Z3 formula to BoolFormula and return its string (boolformula_print style)."""
    bf = z3_to_boolformula(expr, var_to_id)
    return format_boolformula(bf)


def formula_to_boolformula_str(
    expr,
    var_to_id: Optional[Dict["BoolRef", int]] = None,
) -> str:
    """
    Convert a Z3 formula to boolformula_t-style string.
    If var_to_id is provided, use it for literal numbering; otherwise ids are auto-assigned.
    If conversion fails (e.g. quantifiers, pseudo-boolean), return str(expr).
    """
    if not _Z3_AVAILABLE:
        return str(expr)
    try:
        from z3 import BoolRef
        if isinstance(expr, BoolRef):
            return z3_to_boolformula_str(expr, var_to_id)
    except (ValueError, TypeError, AttributeError):
        pass
    return str(expr)
