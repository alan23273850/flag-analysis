"""
Utilities for parsing, combining, and evaluating BLIF files.
Used by the path_5 global coordinated learner.
"""
from __future__ import annotations
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

BlifGate = Tuple[List[str], str, List[List[str]]]

# ---------------------------------------------------------------------------
# BLIF parsing
# ---------------------------------------------------------------------------

def parse_blif(path: Path) -> Tuple[List[str], List[str], List[BlifGate]]:
    """
    Parse a BLIF file. Returns (inputs, outputs, gates).
    Each gate is (input_names, output_name, rows) where rows are truth table rows
    (e.g. ["10", "1"] for XOR).
    """
    text = path.read_text()
    inputs: List[str] = []
    outputs: List[str] = []
    gates: List[BlifGate] = []

    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line == ".end":
            continue
        if line.startswith(".model"):
            continue
        parts = line.split()
        if parts[0] == ".inputs":
            inputs = parts[1:]
        elif parts[0] == ".outputs":
            outputs = parts[1:]
        elif parts[0] == ".names":
            # .names in1 in2 ... out  then rows like "10 1" or "1"
            gate_inputs = parts[1:-1]
            gate_output = parts[-1]
            gates.append((gate_inputs, gate_output, []))
        else:
            # continuation of previous .names: truth table row
            if gates and isinstance(gates[-1][2], list):
                row = line.split()
                if row:
                    gates[-1][2].append(row)
    return inputs, outputs, gates


def rename_blif(
    inputs: List[str],
    outputs: List[str],
    gates: List[Tuple[List[str], str, List[List[str]]]],
    prefix: str,
    keep_inputs: bool = True,
) -> Tuple[List[str], List[str], List[BlifGate]]:
    """Rename all internal names with prefix; optionally keep input names as-is."""
    name_map: Dict[str, str] = {}
    for inp in inputs:
        name_map[inp] = inp if keep_inputs else f"{prefix}{inp}"
    for _, out, _ in gates:
        name_map[out] = f"{prefix}{out}"
    for inps, out, rows in gates:
        for inp in inps:
            if inp not in name_map:
                name_map[inp] = f"{prefix}{inp}"

    def map_name(n: str) -> str:
        return name_map.get(n, n)

    new_inputs = [map_name(i) for i in inputs]
    new_outputs = [map_name(o) for o in outputs]
    new_gates = [
        ([map_name(i) for i in gi], map_name(go), rows)
        for (gi, go, rows) in gates
    ]
    return new_inputs, new_outputs, new_gates


def write_blif(
    out_path: Path,
    inputs: List[str],
    outputs: List[str],
    gates: List[Tuple[List[str], str, List[List[str]]]],
    model_name: str = "formula",
) -> None:
    """Write a single BLIF file."""
    lines = [f".model {model_name}", f".inputs {' '.join(inputs)}", f".outputs {' '.join(outputs)}"]
    for gi, go, rows in gates:
        line = f".names {' '.join(gi)} {go}"
        lines.append(line)
        for row in rows:
            lines.append(" ".join(row))
    lines.append(".end")
    out_path.write_text("\n".join(lines) + "\n")


def evaluate_blif(
    inputs: List[str],
    gates: List[Tuple[List[str], str, List[List[str]]]],
    assignment: Dict[str, int],
) -> Dict[str, int]:
    """
    Given input assignment (name -> 0 or 1), propagate through gates and return values for all nodes.
    """
    val: Dict[str, int] = dict(assignment)
    # Topological order: gates may reference earlier outputs
    for gi, go, rows in gates:
        if go in val:
            continue
        in_vals = [val.get(x, 0) for x in gi]
        out_val = _eval_names(gi, rows, in_vals)
        val[go] = out_val
    return val


def _eval_names(
    gate_inputs: List[str],
    rows: List[List[str]],
    in_vals: List[int],
) -> int:
    """Evaluate one .names gate. rows are e.g. ['10', '1'] or ['1']."""
    if not gate_inputs:
        # 0-input: constant from first row
        if rows and rows[0]:
            return 1 if rows[0][0] == "1" else 0
        return 0
    n = len(gate_inputs)
    for row in rows:
        if len(row) < 2:
            continue
        pattern = row[0]
        result = row[-1]
        if len(pattern) != n:
            continue
        match = True
        for i, c in enumerate(pattern):
            if c == "-":
                continue
            if int(c) != in_vals[i]:
                match = False
                break
        if match:
            return 1 if result == "1" else 0
    return 0


def topological_sort_gates(
    primary_inputs: List[str],
    gates: List[BlifGate],
) -> List[BlifGate]:
    """Order gates so every gate's inputs are defined (primary or earlier output)."""
    defined = set(primary_inputs)
    result: List[BlifGate] = []
    remaining = list(gates)
    while remaining:
        made_progress = False
        for i, (gi, go, rows) in enumerate(remaining):
            if go in defined:
                continue
            if all(x in defined for x in gi):
                result.append((gi, go, rows))
                defined.add(go)
                remaining.pop(i)
                made_progress = True
                break
        if not made_progress:
            break
    result.extend(remaining)
    return result


# ---------------------------------------------------------------------------
# Combine multiple BLIFs with shared inputs
# ---------------------------------------------------------------------------

def combine_blifs(
    blif_specs: List[Tuple[List[str], List[str], List[Tuple[List[str], str, List[List[str]]]]]],
    prefixes: List[str],
    shared_inputs: List[str],
) -> Tuple[List[str], List[str], List[BlifGate]]:
    """
    Combine multiple BLIFs that share the same inputs. Each spec is (inputs, outputs, gates).
    prefixes[i] is used for the i-th spec. shared_inputs is the common input list (e.g. v1..v88).
    """
    all_inputs = list(shared_inputs)
    all_outputs: List[str] = []
    all_gates: List[Tuple[List[str], str, List[List[str]]]] = []

    for (inp, out, gates), prefix in zip(blif_specs, prefixes):
        ninp, nout, ngates = rename_blif(inp, out, gates, prefix, keep_inputs=True)
        all_outputs.extend(nout)
        all_gates.extend(ngates)

    return all_inputs, all_outputs, all_gates
