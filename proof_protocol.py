from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
from pathlib  import Path

from circuit_op import *
from protocol import *

from qiskit import QuantumCircuit


from z3 import BoolVal, Xor, Bool,simplify,substitute, And, Not,Or, PbLe, AtMost,ForAll, Implies, Exists, PbGe, AtLeast

from z3 import Solver, unsat, sat

from flag_analysis import *
from circuit_op import *
from protocol import *

from copy import deepcopy


from dataclasses import dataclass
from typing import List, Dict, Any

@dataclass
class PathRecord:
    nodes: List[str]                   # node IDs from root → leaf
    conditions: List[Dict]             # branch conditions taken
    instr_steps: List[Dict[str, Any]]  # [{"instr": name}, ...]
    state: Dict[str, List[Any]]        # final symbolic state {"data":..., "ancX":..., ...}


from copy import deepcopy

from typing import List, Dict, Any



from typing import List, Dict

def proof_protocol(protocol,
                  start_node: str,
                  init_state,
                  config: Dict, 
                  t: int,
                 stab_txt_path: str
                  ):

    all_paths = []

    def dfs(round_idx: int,
            node_id: str,
            cur_state,
            cur_path: List[Dict]):

        node = protocol[node_id]

        # -------------------------------
        # Execute instruction (if exists)
        # -------------------------------
        instr = node.instructions[0] if node.instructions else None
        state_after = cur_state
        site_info = []   # <--- always create a local holder

        if instr is not None  and node.branches:
            if instr not in config:
                raise KeyError(f"Instruction '{instr}' not found in config")
            
            qasm_path = config[instr]
            qc = load_qasm(qasm_path)
            gate_list = get_gate_only_indices(qc)
            groups = detect_qubit_groups(qc)

            print("round:", round_idx, "node:", node_id, "instr:", instr)
            state_after, site_info = symbolic_execution_of_state(
                qasm_path,
                cur_state,
                round_idx,
                fault_gate=gate_list,
                track_steps=False
            )

        # -------------------------------
        # Leaf node (no branches)
        # -------------------------------
        if not node.branches  :
           
            step = {
                "round": round_idx,
                "node": node_id,
                "next": None,
                "instruction": instr,
                "condition": None,
                "state": { a :  cur_path[-1]["state"][a] for a in ["data"] if a in cur_path[-1]["state"] } ,
                "site_info": site_info   # <-- store entire list
            }
            full_path = cur_path + [step]

            all_paths.append(full_path)
            if instr == 'Break': 
                return
            elif instr.startswith("LUT_") :
                all_condition = [  step["condition"] for step in full_path if step["condition"] is not None ]
                gen_syn =  parse_lut_instr(instr)
                proof_path(full_path, t, gen_syn ,all_condition  ,stab_txt_path)
                return

        # -------------------------------
        # Branching
        # -------------------------------
        for br in node.branches:
            cond_dict = br.condition.to_dict() if br.condition is not None else None

            full_state = [s["state"] for s in cur_path if s["state"] is not None] + [state_after]
            print("condition dict:", cond_dict)
            z3_condition = condition_to_z3(cond_dict, full_state, groups)

            
            
            step = {
                "round": round_idx,
                "node": node_id,
                "next": br.target,
                "instruction": instr,
                "condition": z3_condition,
                "state": state_to_raw_expr_dict(state_after, groups),
                "site_info": site_info   # <-- store entire list
            }

            dfs(
                round_idx + 1,
                br.target,
                state_after,
                cur_path + [step]
            )

    dfs(0, start_node, init_state, [])
    return all_paths
from z3 import BoolVal, And, Or, Not

# -----------------------------
# Parse "s_1", "f_3", etc.
# -----------------------------

def parse_var(name: str):
    """
    Convert 's_1' -> ('s', 1)
            'f_3' -> ('f', 3)
    """
    if "_" not in name:
        raise ValueError(f"Bad variable format (expected like 's_1'): {name}")
    group, idx_str = name.split("_", 1)
    return group, int(idx_str)


def read_state_variable(q_type: str, index: int, state: dict, groups : Dict):
    """
    Map (group, index) into your protocol state structure.

    Expect state like:
        state["syn"]  : List[BoolRef]   # syndrome ancilla bits
        state["flag"] : List[BoolRef]   # flag qubit bits

    's_i' -> state["syn"][i]
    'f_i' -> state["flag"][i]
    """

    state =  state_to_raw_expr_dict(state, groups)
    
    
    if q_type == "s":
        syn = [ q.z for q in state["ancX"]] + [q.z for q in state["ancZ"]]
        if index < 0 or index >= len(syn):
            raise IndexError(f"s_{index} out of range (len syn = {len(syn)})")
        return syn[index]

    if q_type == "f":
       
        flags = [ q.z for q in state["flagX"]] + [q.z for q in state["flagZ"]]
        if index < 0 or index >= len(flags):
            raise IndexError(f"f_{index} out of range (len flag = {len(flags)})")
        return flags[index]

    raise ValueError(f"Unknown variable group in condition: {q_type!r}")


def read_operand(x, full_state: dict , groups:Dict):
    """
    Turn an operand from the condition dict into a z3 Bool expression.

    Supports:
      - 0, 1              → False / True
      - 's_i', 'f_j'      → read from state via (group, index)
    """
    # numeric constants
    if x == 0:
        return False
    if x == 1:
        return True

    if isinstance(x, str):
        q_type, idx = parse_var(x)

        print("read_operand:", q_type, idx)
        if q_type == 's' :
            return ([ q.z for q in state_to_raw_expr_dict(full_state[idx],groups)["ancX"] if state_to_raw_expr_dict(full_state[idx],groups)["ancX"] != [] ] 
                    + [q.x for q in state_to_raw_expr_dict(full_state[idx], groups)["ancZ"]  if state_to_raw_expr_dict(full_state[idx],groups)["ancZ"] != []])
         
        elif q_type == 'f' :   
            return [ q.z for q in state_to_raw_expr_dict(full_state[idx], groups)["flagX"]] + [q.x for q in state_to_raw_expr_dict(full_state[idx], groups)["flagZ"]]

    raise ValueError(f"Unsupported operand in condition: {x!r}")

def parse_lut_instr(name: str):
    """
    Parse something like 'LUT_s_0_f_0_s_1' into:
      [('s', 0), ('f', 0), ('s', 1)]
    """
    if not name.startswith("LUT_"):
        raise ValueError(f"Not a LUT instruction: {name}")

    tokens = name[4:].split("_")  # drop 'LUT_' and split
    if len(tokens) % 2 != 0:
        raise ValueError(f"Bad LUT format: {name}")

    pairs = []
    for kind, idx_str in zip(tokens[0::2], tokens[1::2]):
        if kind not in ("s", "f"):
            raise ValueError(f"Unknown LUT kind '{kind}' in {name}")
        pairs.append((kind, int(idx_str)))
    return pairs
# -----------------------------------
# Main translator: condition → z3 expr
# -----------------------------------

def condition_to_z3(cond: dict | None, full_state: dict, groups:Dict) -> Bool:
    """
    Convert a protocol condition dict into a z3 BoolRef
    using the given `state`, where:

        state["syn"]  : List[BoolRef]  # syndromes (ancilla-based)
        state["flag"] : List[BoolRef]  # flags

    cond JSON forms:
      - {"type":"and", "operands":[cond1, cond2, ...]}
      - {"type":"or",  "operands":[cond1, cond2, ...]}
      - {"type":"not", "operands":[cond1]}
      - {"type":"equal",     "left":..., "right":...}
      - {"type":"not_equal", "left":..., "right":...}
    """
    if cond is None:
        return BoolVal(True)  # no condition → always true

    t = cond["type"]

    if t == "and":
        return And(*(condition_to_z3(c, full_state,  groups) for c in cond["operands"]))

    if t == "or":
        return Or(*(condition_to_z3(c,full_state , groups) for c in cond["operands"]))

    if t == "not":
        sub = cond["operands"][0]  # your JSON uses 'operands' even for NOT
        return Not(condition_to_z3(sub, full_state , groups))

    if t in ("equal"):
        L = read_operand(cond["left"],  full_state , groups)
        R = read_operand(cond["right"], full_state , groups)
        if  isinstance(R, bool):
            return And(*[L[i] == BoolVal(False) for i in range(len(L))])
        else :
            return And(*[L[i] == R[i] for i in range(len(L))])  

    raise ValueError(f"Unknown condition type: {t}")


def proof_path(path : list[dict], t : int , gen_syn : list ,all_condtion : list, stab_txt_path: str) :
    """
    Given a path (list of steps with conditions and states), build a z3 formula
    that encodes the conditions along the path.

    Each step in the path is a dict with keys:
      - "round": int
      - "node": str
      - "next": str | None
      - "instruction": str | None
      - "condition": z3 BoolRef | None
      - "state": dict

    The resulting formula is the AND of all step conditions.
    """
    
    vars = [v for step in path for info in step["site_info"]  for v in info["vars"].values()]
    faults =  [info["act"] for step in path for info in step["site_info"]]
    gen_syn_z3 = []
    for type, idx in gen_syn:
        if type == 's' :
            syn =  [ anc.z for anc in path[idx]["state"]["ancX"]] + [ anc.x for anc in path[idx]["state"]["ancZ"]]
            gen_syn_z3 += syn 
        elif type == 'f' :
            flag =  [ flag.z for flag in path[idx]["state"]["flagX"]] + [ flag.x for flag in path[idx]["state"]["flagZ"] ]
            gen_syn_z3 += flag 

    at_most_t_faults = [AtMost( *faults , t)]

    return uniqness_proof(vars, at_most_t_faults,all_condtion,  gen_syn_z3, path[-1]["state"]["data"],stab_txt_path)
    
    
    

    # Add fault constraints if needed
    

    return simplify(path_formula)