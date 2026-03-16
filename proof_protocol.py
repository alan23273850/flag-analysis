from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
from pathlib  import Path

from circuit_op import *
from protocol import *

from qiskit import QuantumCircuit


from z3 import Bool, BoolVal, Xor, Bool,simplify,substitute, And, Not,Or, PbLe, AtMost,ForAll, Implies, Exists, PbGe, AtLeast,PbEq

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
from typing import Dict, List

from typing import Dict, List

from typing import Dict, List

def to_human_readable(expr):
    """將 Z3 布林方程式轉換為傳統數學邏輯符號，並智慧省略多餘括號"""
    if is_const(expr):
        return str(expr)

    kind = expr.decl().kind()
    children = expr.children()

    # 內部輔助函數：由父節點判斷是否要幫子節點加括號
    def format_child(child, parent_kind):
        if is_const(child):
            return to_human_readable(child)

        child_kind = child.decl().kind()

        # NOT 節點自己有處理邏輯，直接回傳
        if child_kind == Z3_OP_NOT:
            return to_human_readable(child)

        # 如果子節點的運算子跟父節點一樣 (例如 AND 裡面包 AND)，就不加括號
        if child_kind == parent_kind and child_kind in (Z3_OP_AND, Z3_OP_OR, Z3_OP_XOR):
            return to_human_readable(child)

        # 其他情況 (例如 OR 裡面包 AND)，才需要把子節點括起來
        return f"({to_human_readable(child)})"

    # --- 主要運算子處理 ---
    if kind == Z3_OP_NOT:
        child = children[0]
        # NOT 的子節點如果是單一變數或另一個 NOT，不加括號
        if is_const(child) or child.decl().kind() == Z3_OP_NOT:
            return f"¬{to_human_readable(child)}"
        else:
            return f"¬({to_human_readable(child)})"

    elif kind == Z3_OP_AND:
        return " ∧ ".join(format_child(c, kind) for c in children)

    elif kind == Z3_OP_OR:
        return " ∨ ".join(format_child(c, kind) for c in children)

    elif kind == Z3_OP_XOR:
        return " ⊕ ".join(format_child(c, kind) for c in children)

    elif kind == Z3_OP_EQ:
        return f"{format_child(children[0], kind)} ↔ {format_child(children[1], kind)}"

    else:
        return str(expr)

def proof_protocol(protocol,
                  start_node: str,
                  init_state,
                  config: Dict,
                  t: int):

    all_paths = []

    def dfs(round_idx: int,
            node_id: str,
            cur_state,          # CircuitXZ
            cur_groups,         # dict or None
            cur_path: List[Dict]):

        #if len(all_paths) == 2 :  return  # stop exploring more branches

        node = protocol[node_id]

        # -------------------------------
        # Execute instruction (if exists)
        # -------------------------------
        instr = node.instructions[0] if node.instructions else None
        print("Current node:", node_id, "Instruction:", instr)
        state_after = cur_state
        site_info = []
        groups =  cur_groups

        if instr is not None and node.branches:
            # This is a circuit instruction (e.g. flag_syndrome, raw_syndrome)
            if instr not in config:
                raise KeyError(f"Instruction '{instr}' not found in config")

            qasm_path = config[instr]
            qc = load_qasm(qasm_path)
            gate_list = get_gate_only_indices(qc)
            groups = detect_qubit_groups(qc)   # new groups for this circuit

            print("round:", round_idx, "node:", node_id, "instr:", instr)


            state_after, site_info = symbolic_execution_of_state(
                qasm_path,
                cur_state,
                round_idx,
                fault_gate=gate_list,
                track_steps=False
            )
            '''
            for i, q in enumerate(state_after.qubits):
              #print(f"  q[{i}]: X = {q.x}, Z = {q.z}")
            '''

        # -------------------------------
        # Build dict-view of state_after
        # -------------------------------
        if groups is not None:

            state_dict = state_to_raw_expr_dict(state_after, groups)

        elif instr is None or instr == 'Break' or instr.startswith("LUT_"):
            print("in condition round:", round_idx, "node:", node_id, "instr:", instr)
            # Break instruction with no anc/flag structure; keep only data as a list
            #print("state after:", state_after)
            state_dict = state_to_raw_expr_dict(state_after, groups)



        else:
            # no anc/flag structure yet; keep only data as a list
            state_dict = state_to_raw_expr_dict(state_after, groups)

        if f'{instr}' + '_flag_group' in config:
            import json
            with open(config[f'{instr}' + '_flag_group'], "r") as f:
                flag_group= json.load(f)

            state_dict['flagX'] = [[state_dict.copy()['flagX'][i] for i in g ] for g in flag_group['flagX']]
            state_dict['flagZ'] = [[state_dict.copy()['flagZ'][i] for i in g ] for g in flag_group['flagZ']]

        # -------------------------------
        # Leaf node (no branches)
        # -------------------------------
        if not node.branches:
            step = {
                "round": round_idx,
                "node": node_id,
                "next": None,
                "instruction": instr,
                "condition": None,
                "state": state_dict,   # always dict
                "site_info": site_info
            }
            full_path = cur_path + [step]

            all_paths.append(full_path)

            # Leaf behavior
            if instr == 'Break':
                return
            elif instr and instr.startswith("LUT_"):  # FIX: Added null check
                print("LUT", instr)
                all_condition = [
                    s["condition"] for s in full_path
                    if s["condition"] is not None
                ]
                gen_syn = parse_lut_instr(instr)
                print("gen_syn:", gen_syn)

                ####check the last syndrome is the correct syndrome due to the data qubits
                # final data error
                E_x = [dq.x for dq in full_path[-1]["state"]["data"]]
                E_z = [dq.z for dq in full_path[-1]["state"]["data"]]


                gens, syn_measured = last_ancilla_formulas(full_path, config)
                pred_syn = stabilizer_syndrome_from_data(E_x, E_z, gens)



                # enforce: measured syndrome == commutation syndrome
                syn_constraint = And(*[
                    s_m == s_p for s_m, s_p in zip(syn_measured, pred_syn)
                ])



                all_condition = all_condition + [syn_constraint]


                print("len all paths:", len(all_paths))

                proof_path(full_path, t, gen_syn, all_condition, config['stab_txt_path'], config['log_txt_path'])

                return
            else:
                # leaf with some other instruction, but nothing to prove
                return

        # -------------------------------
        # Branching
        # -------------------------------
        for br in node.branches:
            cond_dict = br.condition.to_dict() if br.condition is not None else None

            # all states in full_state are dicts now
            full_state = [s["state"] for s in cur_path] + [state_dict]
            z3_condition = condition_to_z3(cond_dict, full_state, groups)

            step = {
                "round": round_idx,
                "node": node_id,
                "next": br.target,
                "instruction": instr,
                "condition": z3_condition,
                "state": state_dict,       # dict
                "site_info": site_info
            }
            next_groups = data_only_groups_from_state_dict(state_dict)

            dfs(
                round_idx + 1,
                br.target,
                state_after,   # still CircuitXZ
                next_groups,        # carry the same groups forward
                cur_path + [step]
            )


    # initial call: no groups yet, clean data state



    dfs(0, start_node, init_state, None, [])
    return all_paths



def proof_protocol_boolean(protocol,
                  start_node: str,
                  init_state,
                  config: Dict,
                  t: int):

    all_paths = []
    all_path_data = []  # Store collected data for each path

    def dfs(round_idx: int,
            node_id: str,
            cur_state,          # CircuitXZ
            cur_groups,         # dict or None
            cur_path: List[Dict]):

        #if len(all_paths) == 2 :  return  # stop exploring more branches

        node = protocol[node_id]

        # -------------------------------
        # Execute instruction (if exists)
        # -------------------------------
        instr = node.instructions[0] if node.instructions else None
        # print("Current node:", node_id, "Instruction:", instr)
        state_after = cur_state
        site_info = []
        groups =  cur_groups

        if instr is not None and node.branches:
            # This is a circuit instruction (e.g. flag_syndrome, raw_syndrome)
            if instr not in config:
                raise KeyError(f"Instruction '{instr}' not found in config")

            qasm_path = config[instr]
            qc = load_qasm(qasm_path)
            gate_list = get_gate_only_indices(qc)
            groups = detect_qubit_groups(qc)   # new groups for this circuit

            # print("round:", round_idx, "node:", node_id, "instr:", instr)


            state_after, site_info = symbolic_execution_of_state(
                qasm_path,
                cur_state,
                round_idx,
                fault_gate=gate_list,
                track_steps=False
            )
            '''
            for i, q in enumerate(state_after.qubits):
              #print(f"  q[{i}]: X = {q.x}, Z = {q.z}")
            '''

        # -------------------------------
        # Build dict-view of state_after
        # -------------------------------
        if groups is not None:

            state_dict = state_to_raw_expr_dict(state_after, groups)

        elif instr is None or instr == 'Break' or instr.startswith("LUT_"):
            # print("in condition round:", round_idx, "node:", node_id, "instr:", instr)
            # Break instruction with no anc/flag structure; keep only data as a list
            #print("state after:", state_after)
            state_dict = state_to_raw_expr_dict(state_after, groups)

        else:
            # no anc/flag structure yet; keep only data as a list
            state_dict = state_to_raw_expr_dict(state_after, groups)

        if f'{instr}' + '_flag_group' in config:
            import json
            with open(config[f'{instr}' + '_flag_group'], "r") as f:
                flag_group= json.load(f)

            state_dict['flagX'] = [[state_dict.copy()['flagX'][i] for i in g ] for g in flag_group['flagX']]
            state_dict['flagZ'] = [[state_dict.copy()['flagZ'][i] for i in g ] for g in flag_group['flagZ']]

        # -------------------------------
        # Leaf node (no branches)
        # -------------------------------
        if not node.branches:
            step = {
                "round": round_idx,
                "node": node_id,
                "next": None,
                "instruction": instr,
                "condition": None,
                "state": state_dict,   # always dict
                "site_info": site_info
            }
            full_path = cur_path + [step]

            all_paths.append(full_path)

            # Collect the three pieces of information for this path
            # 1. Last round data qubit formulas
            last_data = full_path[-1]["state"]["data"]

            # 2. All ancilla and flag formulas in each round
            anc_flag_per_round = []
            for step_info in full_path:
                # Extract ancilla formulas (syndrome measurements)
                ancX_list = step_info["state"].get("ancX", [])
                ancZ_list = step_info["state"].get("ancZ", [])
                flagX_list = step_info["state"].get("flagX", [])
                flagZ_list = step_info["state"].get("flagZ", [])

                # For ancX: measure in Z basis (use .z)
                ancX_formulas = [q.z for q in ancX_list]

                # For ancZ: measure in X basis (use .x)
                ancZ_formulas = [q.x for q in ancZ_list]

                # For flagX: measure in Z basis (use .z), handle nested lists
                flagX_formulas = []
                for item in flagX_list:
                    if isinstance(item, list):
                        flagX_formulas.append([q.z for q in item])
                    else:
                        flagX_formulas.append(item.z)

                # For flagZ: measure in X basis (use .x), handle nested lists
                flagZ_formulas = []
                for item in flagZ_list:
                    if isinstance(item, list):
                        flagZ_formulas.append([q.x for q in item])
                    else:
                        flagZ_formulas.append(item.x)

                round_anc_flag = {
                    "round": step_info["round"],
                    "ancX_formulas": ancX_formulas,  # Z-basis measurements
                    "ancZ_formulas": ancZ_formulas,  # X-basis measurements
                    "flagX_formulas": flagX_formulas,  # Z-basis measurements (may be nested)
                    "flagZ_formulas": flagZ_formulas   # X-basis measurements (may be nested)
                }
                anc_flag_per_round.append(round_anc_flag)

            # 3. All conditions along the path
            path_conditions = [
                s["condition"] for s in full_path
                if s["condition"] is not None
            ]

            # Store collected data
            faults = [info["act"] for step in full_path for info in step["site_info"]]
            path_data = {
                "last_data": last_data,
                "anc_flag_per_round": anc_flag_per_round,
                "conditions": path_conditions,
                "faults": faults
            }
            all_path_data.append(path_data)

            # Leaf behavior
            if instr == 'Break':
                return
            elif instr and instr.startswith("LUT_"):  # FIX: Added null check
                # print("LUT", instr)
                gen_syn = parse_lut_instr(instr)
                # print("gen_syn:", gen_syn)

                # print("len all paths:", len(all_paths))

                # Print collected data for this path
                # print(f"Path {len(all_paths)} data collected:")
                # print(f"  - Last data qubits: {len(last_data)} qubits")
                # print(f"  - Rounds with anc/flag: {len(anc_flag_per_round)}")
                # print(f"  - Path conditions: {len(path_conditions)}")

                return
            else:
                # leaf with some other instruction, but nothing to prove
                return

        # -------------------------------
        # Branching
        # -------------------------------
        for br in node.branches:
            cond_dict = br.condition.to_dict() if br.condition is not None else None

            # all states in full_state are dicts now
            full_state = [s["state"] for s in cur_path] + [state_dict]
            z3_condition = condition_to_z3(cond_dict, full_state, groups)

            step = {
                "round": round_idx,
                "node": node_id,
                "next": br.target,
                "instruction": instr,
                "condition": z3_condition,
                "state": state_dict,       # dict
                "site_info": site_info
            }
            next_groups = data_only_groups_from_state_dict(state_dict)

            dfs(
                round_idx + 1,
                br.target,
                state_after,   # still CircuitXZ
                next_groups,        # carry the same groups forward
                cur_path + [step]
            )


    # initial call: no groups yet, clean data state



    dfs(0, start_node, init_state, None, [])

    # Print all collected path data
    print("\n" + "="*80)
    print(f"COLLECTED DATA FROM {len(all_path_data)} PATHS")
    print("="*80)
    for i, path_data in enumerate(all_path_data):
        if i == 0: continue # 1 <= i <= 4
        print(f"\n--- PATH {i+1} ---")
        print(f"\n0. Fault variables (count = {len(path_data['faults'])}):")
        for f_idx, f in enumerate(path_data["faults"]):
            print(f"   f[{f_idx}]: {f}")

        # Print last round data qubits
        print(f"\n1. Last Round Data Qubits ({len(path_data['last_data'])} qubits):")
        for idx, dq in enumerate(path_data['last_data']):
            print(f"   Data[{idx}]: X={dq.x}, Z={dq.z}")

        # Print ancilla and flag formulas for each round
        print(f"\n2. Ancilla & Flag Formulas per Round ({len(path_data['anc_flag_per_round'])} rounds):")
        for round_data in path_data['anc_flag_per_round']:
            print(f"   Round {round_data['round']}:")

            # Print ancX formulas (Z-basis syndrome measurements)
            ancX_formulas = round_data['ancX_formulas']
            print(f"     ancX (Z-basis): {len(ancX_formulas)} measurements")
            for idx, formula in enumerate(ancX_formulas):
                print(f"       ancX[{idx}]: {formula}")

            # Print ancZ formulas (X-basis syndrome measurements)
            ancZ_formulas = round_data['ancZ_formulas']
            print(f"     ancZ (X-basis): {len(ancZ_formulas)} measurements")
            for idx, formula in enumerate(ancZ_formulas):
                print(f"       ancZ[{idx}]: {formula}")

            # Print flagX formulas (Z-basis flag measurements)
            flagX_formulas = round_data['flagX_formulas']
            print(f"     flagX (Z-basis): {len(flagX_formulas)} groups/measurements")
            for idx, formula in enumerate(flagX_formulas):
                if isinstance(formula, list):
                    print(f"       flagX[{idx}] (group of {len(formula)}):")
                    for sub_idx, sub_formula in enumerate(formula):
                        print(f"         [{sub_idx}]: {sub_formula}")
                else:
                    print(f"       flagX[{idx}]: {formula}")

            # Print flagZ formulas (X-basis flag measurements)
            flagZ_formulas = round_data['flagZ_formulas']
            print(f"     flagZ (X-basis): {len(flagZ_formulas)} groups/measurements")
            for idx, formula in enumerate(flagZ_formulas):
                if isinstance(formula, list):
                    print(f"       flagZ[{idx}] (group of {len(formula)}):")
                    for sub_idx, sub_formula in enumerate(formula):
                        print(f"         [{sub_idx}]: {sub_formula}")
                else:
                    print(f"       flagZ[{idx}]: {formula}")

        # Print path conditions
        print(f"\n3. Path Conditions ({len(path_data['conditions'])} conditions):")
        for cond_idx, cond in enumerate(path_data['conditions']):
            print(f"   Condition {cond_idx}: {cond}")

        # We want to learn the decoder in this section.
        print(f"\n4. Decoder Boolean formulas (from truth table, DNF):")
        solver = Solver()
        #############################

        # Build fresh Boolean variables for every fault in
        # path_data['faults'], and constrain them equal to
        # the original formulas.
        fault_bools = []
        for fault_idx, fault in enumerate(path_data["faults"]):
            fresh = Bool(f"fault_{fault_idx}")
            solver.add(fresh == fault)
            fault_bools.append(fresh)

        # Require that **exactly one** fault occurs (i.e. among all
        # fault_bools, there is one and only one that is True).
        # If there are no faults, we leave the solver unconstrained here.
        if len(fault_bools) == 1:
            solver.add(fault_bools[0])
        elif len(fault_bools) > 1:
            # Sum of fault_bools treated as 0/1 must equal 1.
            solver.add(PbEq([(b, 1) for b in fault_bools], 1))

        # For each condition in path_data["conditions"], create a fresh
        # Boolean variable condition_i and equate it to that condition.
        condition_bools = []
        for cond_idx, cond in enumerate(path_data["conditions"]):
            c = Bool(f"condition_{cond_idx}")
            solver.add(c == cond)
            condition_bools.append(c)

        # conditionAll is true iff all condition_i are true.
        conditionAll = Bool("conditionAll")
        if condition_bools:
            if len(condition_bools) == 1:
                solver.add(conditionAll == condition_bools[0])
            else:
                solver.add(conditionAll == And(*condition_bools))
            solver.add(conditionAll)

        # For each final data qubit in path_data["last_data"], create
        # fresh variables data{i}_x and data{i}_z for its X/Z formulas.
        # Then create corresponding decoder variables dec{i}_x, dec{i}_z,
        # and fixed variables fixed{i}_x, fixed{i}_z such that
        # fixed{i}_x == Xor(data{i}_x, dec{i}_x) and similarly for Z.
        n_data = len(path_data["last_data"])
        decoder_vars = []
        for data_idx, dq in enumerate(path_data["last_data"]):
            dx = Bool(f"data{data_idx}_x")
            dz = Bool(f"data{data_idx}_z")
            solver.add(dx == dq.x)
            solver.add(dz == dq.z)

            decx = Bool(f"dec{data_idx}_x")
            decz = Bool(f"dec{data_idx}_z")
            decoder_vars.extend([decx, decz])

            fixedx = Bool(f"fixed{data_idx}_x")
            fixedz = Bool(f"fixed{data_idx}_z")

            solver.add(fixedx == Xor(dx, decx))
            solver.add(fixedz == Xor(dz, decz))

        # Read logical operators and stabilizers from config paths.
        # Each line has two bitstrings xs, zs. For each pair, enforce
        # that its length matches the number of data qubits, and build
        # a commute_i variable that is the XOR over j of
        #   (fixed_j_x & zs[j]) XOR (fixed_j_z & xs[j]).
        bit_pairs: list[tuple[str, str]] = []
        for txt_key in ["log_txt_path", "stab_txt_path"]:
            txt_path = str(config[txt_key])
            with open(txt_path, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    xs, zs = line.split()
                    if len(xs) != n_data or len(zs) != n_data:
                        raise ValueError(
                            f"Length mismatch in {txt_key}: "
                            f"len(xs)={len(xs)}, len(zs)={len(zs)}, expected {n_data}"
                        )
                    bit_pairs.append((xs, zs))

        commute_bools = []
        for row_idx, (xs, zs) in enumerate(bit_pairs):
            terms = []
            for j in range(n_data):
                fx = Bool(f"fixed{j}_x")
                fz = Bool(f"fixed{j}_z")

                # Build (fixed_j_x & zs[j]) XOR (fixed_j_z & xs[j])
                per_terms = []
                if zs[j] == "1":
                    per_terms.append(fx)
                if xs[j] == "1":
                    per_terms.append(fz)

                if len(per_terms) == 1:
                    terms.append(per_terms[0])
                elif len(per_terms) == 2:
                    terms.append(Xor(per_terms[0], per_terms[1]))
                # if both bits are '0', this contributes nothing

            if not terms:
                row_parity = BoolVal(False)
            elif len(terms) == 1:
                row_parity = terms[0]
            else:
                acc = terms[0]
                for t in terms[1:]:
                    acc = Xor(acc, t)
                row_parity = acc

            commute = Bool(f"commute_{row_idx}")
            solver.add(commute == Not(row_parity))
            commute_bools.append(commute)

        # Fresh variable that is the AND of all commute_i.
        if commute_bools:
            all_commute = Bool("all_commute")
            if len(commute_bools) == 1:
                solver.add(all_commute == commute_bools[0])
            else:
                solver.add(all_commute == And(*commute_bools))
        # We will assert that all commute constraints hold later.

        # Build fresh Boolean variables for every ancilla/flag formula
        # in path_data['anc_flag_per_round'], and constrain them equal
        # to the original formulas. Collect these measurement variables
        # (r_...) as decoder inputs.
        meas_vars = []
        for round_data in sorted(path_data["anc_flag_per_round"], key=lambda x: x['round']):
            for key in ["ancX_formulas", "ancZ_formulas", "flagX_formulas", "flagZ_formulas"]:
                formulas = round_data.get(key, [])
                key = key[:-len('_formulas')]
                idx = 0
                for f in formulas:
                    # Handle nested groups (lists of formulas) as in flagX/flagZ
                    if isinstance(f, list):
                        idx2 = 0
                        for sub_f in f:
                            v = Bool(f"r_{round_data['round']}_{key}{idx}_{idx2}")
                            solver.add(v == sub_f)
                            meas_vars.append(v)
                            idx2 += 1
                    else:
                        v = Bool(f"r_{round_data['round']}_{key}{idx}")
                        solver.add(v == f)
                        meas_vars.append(v)
                    idx += 1

        # Decoder Learning in C: pass SMT2 + variable name lists to C routine
        #   decoder_learning_in_C. Set _use_c_decoder_learning to True when
        #   libdecoder_learning.so is built and available.
        _use_c_decoder_learning = True
        if _use_c_decoder_learning:
            learned_formulas_dict = {}

            import json
            import ctypes
            try:
                here = os.path.dirname(os.path.abspath(__file__))
                so_path = os.path.join(here, "bull", "trunk", "Src", "c_examples", "libdecoder_learning.so")
                lib = ctypes.CDLL(so_path)

                # 1. 準備純字串 List
                smt2_str = solver.to_smt2()
                meas_var_names = [str(v) for v in meas_vars]
                decoder_var_names = [str(v) for v in decoder_vars]
                all_commute_name = "all_commute"

                # 2. 將字串轉為 bytes 列表
                meas_bytes = [name.encode("utf-8") for name in meas_var_names]
                dec_bytes = [name.encode("utf-8") for name in decoder_var_names]

                # 3. 建立 ctypes 的字串陣列 (char**)
                MeasArrayType = ctypes.c_char_p * len(meas_bytes)
                meas_c_array = MeasArrayType(*meas_bytes)

                DecArrayType = ctypes.c_char_p * len(dec_bytes)
                dec_c_array = DecArrayType(*dec_bytes)

                # 4. 設定 C 函數的參數型態 (注意：多了長度參數)
                lib.decoder_learning_in_C.argtypes = [
                    ctypes.c_char_p,                 # smt2_str
                    ctypes.POINTER(ctypes.c_char_p), # meas_names (char**)
                    ctypes.c_int,                    # num_meas (陣列長度)
                    ctypes.POINTER(ctypes.c_char_p), # decoder_names (char**)
                    ctypes.c_int,                    # num_decoders (陣列長度)
                    ctypes.c_char_p                  # all_commute_name
                ]
                lib.decoder_learning_in_C.restype = ctypes.c_void_p

                # 5. 呼叫 C 函數，直接傳入陣列與長度
                ptr = lib.decoder_learning_in_C(
                    smt2_str.encode("utf-8"),
                    meas_c_array,
                    len(meas_bytes),
                    dec_c_array,
                    len(dec_bytes),
                    all_commute_name.encode("utf-8")
                )

                try:
                    # 4. 將 void 指標轉型為 char 指標，並讀取 C 語言產生的字串
                    result_str = ctypes.cast(ptr, ctypes.c_char_p).value.decode("utf-8")

                    # 5. 直接在 Python 記憶體中解析 JSON
                    c_formulas_json = json.loads(result_str)

                    # 建立翻譯字典：{ 'r_0_ancX0': r_0_ancX0的底層宣告, ... }
                    meas_decls = {v.decl().name(): v.decl() for v in meas_vars}
                    # 將 SMT2 字串解析回 Z3 AST (接續你原本的邏輯)
                    for dec_name, smt2_expr in c_formulas_json.items():
                        parsed_ast_vector = parse_smt2_string(f"(assert {smt2_expr})", decls=meas_decls)
                        learned_formulas_dict[dec_name] = parsed_ast_vector[0]

                finally:
                    # 6. 【極度重要】不論 Python 解析 JSON 是否出錯，都必須呼叫 C 釋放記憶體
                    if ptr:
                        lib.free_c_string(ptr)

            except Exception as e:
                _use_c_decoder_learning = False
                print("decoder_learning_in_C failed", e)
                quit()

        # Synthesize an explicit Boolean formula for each decoder variable
        # from the truth table (DNF: one term per row where output is True).
        def truth_table_to_formula(dec_var, meas_vars, decoder_tables):
            # DNF: OR over all rows where dec_var is True.
            # For each such row_key, term = AND over j: (meas_vars[j] if row_key[j] else Not(meas_vars[j])).
            terms = []
            for row_key, row_entry in decoder_tables.items():
                if row_entry[dec_var] is not True:
                    continue
                literals = []
                for j, v in enumerate(meas_vars):
                    if row_key[j]:
                        literals.append(v)
                    else:
                        literals.append(Not(v))
                terms.append(And(*literals) if literals else BoolVal(True))
            if not terms:
                return BoolVal(False)
            if len(terms) == 1:
                return terms[0]
            return Or(*terms)

        # Decoder Learning in Python:
        #   Build a truth table mapping the measurement variables
        #   (meas_vars) to each decoder variable decoder{i}_x / decoder{i}_z.
        #   Loop:
        #     1. Call solver.check().
        #     2. If SAT, read the model, project it to (meas_vars, decoder_vars)
        #        and add a new row to the truth table.
        #        Then add a blocking clause on meas_vars so that the next
        #        model must use a different measurement pattern.
        #     3. If UNSAT, stop; at this point the accumulated truth table
        #        implicitly defines Boolean formulas for each decoder.
        #   We organize it as:
        #       decoder_tables: Dict[row_key: Tuple[bool,...],
        #                            Dict[decoder_var, bool]]
        if not _use_c_decoder_learning:
            decoder_tables: dict[tuple[bool, ...], dict] = {}

            while True:
                # Construct decoder formulas from the table.
                solver.push() ##########################################################
                solver.add(Not(all_commute))
                # print(decoder_tables)
                # print('Table Size:', len(decoder_tables))
                # print("---- start decoder formulas ----")
                for dec_var in decoder_vars:
                    formula = truth_table_to_formula(dec_var, meas_vars, decoder_tables)
                    formula_flat = simplify(formula)
                    set_option(max_width=10000)
                    # print(f"  {dec_var}: {formula_flat}")
                    solver.add(dec_var == formula_flat)
                # print("----- end decoder formulas -----\n")
                ############################################

                # print(solver.to_smt2())
                res = solver.check()
                if res == unsat:
                    solver.pop() ##########################################################
                    break
                if res != sat:
                    print('The solver cannot solve the following smt2.')
                    print(solver.to_smt2())
                    quit()

                model = solver.model()
                solver.pop() ##########################################################
                # Current measurement pattern (row key)
                row_key = tuple(bool(model.eval(v)) for v in meas_vars)

                # Fill table row: map each decoder output under this measurement
                # pattern. One row_key maps to a dict {decoder_var -> bool}.
                if row_key in decoder_tables:
                    print("The row_key should not exist in decoder_tables.")
                    quit()

                # Now we want to find a feasible decoder under this row_key.
                solver.push() ##########################################################
                solver.add(all_commute)
                for i, meaV in enumerate(meas_vars):
                    solver.add(meaV == BoolVal(row_key[i]))
                res = solver.check()
                if res != sat:
                    print('So strange!')
                    quit()
                model = solver.model()
                solver.pop() ##########################################################

                # Fill the decoder truth table for this row_key using the decoder
                # assignments in the model. If any value is not decided by the model,
                # bool(model.eval(...)) will raise, which is what we want for debugging.
                decoder_tables[row_key] = {
                    decV: bool(model.eval(decV))
                    for decV in decoder_vars
                }

        # At this point, solver is UNSAT under additional blocking clauses
        # on meas_vars, and decoder_tables holds the learned truth tables.
        # We verify the learned decoders without blocking clauses again.
        solver.add(Not(all_commute))
        # Write decoder formulas for this path (machine-readable for Python)
        decoder_basename = f"decoder_C_{i}.txt" if _use_c_decoder_learning else f"decoder_Python_{i}.txt"
        decoder_out_path = Path(__file__).parent / decoder_basename
        with open(decoder_out_path, "w", encoding="utf-8") as dec_file:
            # Format for Python: line1 = meas_var_names (comma-sep); then per line: dec_name TAB sexpr.
            # Rebuild formula_flat: decls = {name: Bool(name) for name in meas_names}; add dec names to decls;
            # parsed = parse_smt2_string("(assert (= " + dec_name + " " + sexpr + "))", decls=decls); formula_flat = parsed[0].children()[1]
            dec_file.write("# decoder path " + str(i) + ": meas_var_names then dec_name\\tsexpr per line (parse (assert (= name sexpr)) with decls)\n")
            dec_file.write("meas_var_names:\t" + ",".join(str(v) for v in meas_vars) + "\n")
            # Each following line: decoder_var_name TAB formula_flat.sexpr()
            # To reconstruct: parse_smt2_string("(assert (= " + dec_name + " " + sexpr + "))", decls=decls)
            # then take the right-hand side of the equality as formula_flat.
            for dec_var in decoder_vars:
                if not _use_c_decoder_learning:
                    formula = truth_table_to_formula(dec_var, meas_vars, decoder_tables)
                else:
                    formula = learned_formulas_dict[str(dec_var)]
                formula_flat = simplify(formula)
                dec_file.write(f"{dec_var}\t{formula_flat.sexpr()}\n")
                solver.add(dec_var == formula_flat)
        print(f"  Decoder formulas written to {decoder_out_path}")
        if not _use_c_decoder_learning:
            print('Table Size:', len(decoder_tables))

        print(f"\n5. Verified SMT formula:")
        print(solver.to_smt2())
        res = solver.check()
        if res == unsat:
            print('Verified!')
        else:
            print('Unverified!')
            quit()

    print("\n" + "="*80)

    return all_paths, all_path_data
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


from z3 import BoolVal, BoolRef


def read_operand(x, full_state , group):
    """
    Interpret operands:
      - int / bool
      - 's_k' (syndrome bit from round k)
      - 'f_k' (flag bit from round k)
    """

    # ints / bools
    if isinstance(x, int) or isinstance(x, bool):
        return x

    # already Z3
    if isinstance(x, BoolRef):
        return x

    # strings like 's_0', 'f_1', ...
    if isinstance(x, str):
        parts = x.split("_")
        if len(parts) != 2:
            raise ValueError(f"Bad variable name in condition: {x!r}")
        kind, idx_str = parts
        round_idx = int(idx_str)  # This is the round index, not array index

        if not full_state:
            raise ValueError("full_state is empty in read_operand")

        # Check if we have enough rounds
        if round_idx >= len(full_state):
            raise IndexError(f"{kind}_{round_idx} refers to round {round_idx}, but we only have {len(full_state)} rounds")

        state_dict = full_state[round_idx]   # Get state from the specific round
        ancX_list = state_dict.get("ancX", [])
        ancZ_list = state_dict.get("ancZ", [])
        flagX_list = state_dict.get("flagX", [])
        flagZ_list = state_dict.get("flagZ", [])

        if kind == "s"  :

            # ancX.z followed by ancZ.x - take all syndrome bits from this round
            syn_bits = [q.z for q in ancX_list] + [q.x for q in ancZ_list]
            if not syn_bits:
                raise IndexError(f"s_{round_idx} refers to round {round_idx}, but no syndrome bits found in that round")

            # Return the OR of all syndrome bits (this represents "any syndrome fired")

            return syn_bits

        elif kind == "f"  :



            flag_bits = [[q.z for q in (g if isinstance(g, list) else [g])]  for g in flagX_list ] + \
                        [[q.x for q in (g if isinstance(g, list) else [g])]  for g in flagZ_list ]


            # flagX.z followed by flagZ.x - take all flag bits from this round
            if not flag_bits:
                raise IndexError(f"f_{round_idx} refers to round {round_idx}, but no flag bits found in that round")

            # Return the OR of all flag bits (this represents "any flag fired")

            return flag_bits

        else:
            raise ValueError(f"Unknown condition variable kind: {kind!r}")

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

def parse_pauli_instruction(instr: str):
    """
    Example:
      'XIXZZ_s IYXXY_f' → [('XIXZZ','s'), ('IYXXY','f')]
    """
    out = []
    for tok in instr.split():
        if "_" not in tok:
            continue
        p, tag = tok.split("_")
        out.append((p, tag))
    return out
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

    if t == "equal":
        L = read_operand(cond["left"],  full_state, groups)
        R = read_operand(cond["right"], full_state, groups)

        # Normalize Python bools to Z3


        # --- list vs bool (e.g. s_0 == 0) ---
        if isinstance(L, list) and isinstance(R, bool):
            L = [x for sub in (L if isinstance(L, list) else [L])
         for x in (sub if isinstance(sub, list) else [sub])]
            if R:
                # every bit in L must be true
                return And(*[li == BoolVal(True) for li in L])
            else:
                # every bit in L must be false
                return And(*[li == BoolVal(False) for li in L])

        if isinstance(R, list) and isinstance(L, bool):
            R = [x for sub in (R if isinstance(R, list) else [R])
         for x in (sub if isinstance(sub, list) else [sub])]
            if L:
                # every bit in L must be true
                return And(*[ri == BoolVal(True) for ri in L])
            else:
                # every bit in L must be false
                return And(*[ri == BoolVal(False) for ri in L])

        if isinstance(L, list) and isinstance(R, int):
            if R < 0:
                raise ValueError(f"Negative integer in equality: {R}")
            if R > len(L):
                raise ValueError(f"Integer in equality exceeds list length: {R} > {len(L)}")

            else:
                from z3 import PbEq
                return PbEq([(Or(li), 1) for li in L], R)

        if isinstance(R, list) and isinstance(L, int):
            if L < 0:
                raise ValueError(f"Negative integer in equality: {L}")
            if L > len(R):
                raise ValueError(f"Integer in equality exceeds list length: {L} > {len(R)}")

            else:
                from z3 import PbEq
                return PbEq([(Or(ri), 1) for ri in R], L)





        # --- list vs list: pairwise equality ---
        if isinstance(L, list) and isinstance(R, list):
            if len(L) != len(R):
                raise ValueError(
                    f"List lengths differ in equality: {len(L)} vs {len(R)}"
                )
            return And(*[li == ri for li, ri in zip(L, R)])

        # --- scalar vs scalar (Z3 or int/bool) ---
        # By here, L and R should both be scalar Z3 expressions or ints.
        # If ints 0/1 sneak in, you can map them to BoolVal as well:
        if isinstance(L, int):
            L = BoolVal(bool(L))
        if isinstance(R, int):
            R = BoolVal(bool(R))

        # Now they are both scalar Z3 expressions → simple equality.
        return L == R


def proof_path(path : list[dict], t : int , gen_syn : list ,all_condtion : list, stab_txt_path: str,  log_txt_path: str) :
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
        # --- flatten helper (no function) ---
        ancX = path[idx]["state"]["ancX"]
        ancZ = path[idx]["state"]["ancZ"]
        flagX = path[idx]["state"]["flagX"]
        flagZ = path[idx]["state"]["flagZ"]

        ancX = [q for g in ancX for q in (g if isinstance(g, list) else [g])]
        ancZ = [q for g in ancZ for q in (g if isinstance(g, list) else [g])]
        flagX = [q for g in flagX for q in (g if isinstance(g, list) else [g])]
        flagZ = [q for g in flagZ for q in (g if isinstance(g, list) else [g])]

        if type == "s":
            syn = [a.z for a in ancX] + [a.x for a in ancZ]
            gen_syn_z3 += syn

        elif type == "f":
            flg = [q.z for q in flagX] + [q.x for q in flagZ]
            gen_syn_z3 += flg


    #    print("gen_syn_z3:", gen_syn_z3)
    at_most_t_faults = [PbEq( [(f,1) for f in faults ], t)]






    #return uniqness_proof(vars, at_most_t_faults,all_condtion,  gen_syn_z3, path[-1]["state"]["data"],stab_txt_path, log_txt_path)
   # return  uniqueness_build_goal(vars, at_most_t_faults,all_condtion,  gen_syn_z3, path[-1]["state"]["data"],stab_txt_path, log_txt_path)
    #return uniqueness_solve_with_cryptominisat(vars, at_most_t_faults,all_condtion,  gen_syn_z3, path[-1]["state"]["data"],stab_txt_path, log_txt_path)

    status, model_lits, out = uniqueness_solve_with_cryptominisat(vars, at_most_t_faults,all_condtion,  gen_syn_z3, path[-1]["state"]["data"],stab_txt_path, log_txt_path)



    return 0

    # Add fault constraints if needed


    return simplify(path_formula)