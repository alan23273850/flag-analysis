import json
import os
from typing import Optional, List, Dict

class Condition:
    def __init__(self, cond_type: str, left=None, right=None, operand=None, operands=None):
        self.type = cond_type
        self.left = left          # For binary conditions
        self.right = right        # For binary conditions
        self.operand = operand    # For unary conditions (like 'not')
        self.operands = operands  # For n-ary conditions (like 'and', 'or')

    def to_dict(self):
        d = {"type": self.type}
        if self.type in ("equal", "not_equal"):
            d["left"] = self.left
            d["right"] = self.right
        elif self.type == "not":
            d["operand"] = self.operand.to_dict()
        elif self.type in ("and", "or"):
            d["operands"] = [op.to_dict() for op in self.operands]
        else:
            raise ValueError(f"Unknown condition type: {self.type}")
        return d


class Branch:
    def __init__(self, target: str, condition: Optional[Condition] = None):
        self.target = target
        self.condition = condition

    def to_dict(self):
        return {
            "condition": self.condition.to_dict() if self.condition else None,
            "target": self.target
        }


class Node:
    def __init__(self, node_id: str, instructions: List[str], branches: List[Branch]):
        self.node_id = node_id
        self.instructions = instructions
        self.branches = branches

    def to_dict(self):
        return {
            "instructions": self.instructions,
            "branches": [b.to_dict() for b in self.branches]
        }


class Protocol:
    def __init__(self, start_node: str):
        self.start_node = start_node
        self.nodes: Dict[str, Node] = {}

    def add_node(self, node: Node):
        self.nodes[node.node_id] = node

    def to_dict(self):
        return {
            "start_node": self.start_node,
            "nodes": {nid: node.to_dict() for nid, node in self.nodes.items()}
        }

    def save_to_file(self, filepath: str):
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, "w") as f:
            json.dump(self.to_dict(), f, indent=4)
        print(f"Protocol saved to {filepath}")
        
        
def build_protocol_d_5_lai() -> Protocol:
    protocol = Protocol(start_node="root")

    # Root node with unconditional branch to flag_measure
    root_node = Node(
        node_id="root",
        instructions=["flagged_syndrome"],
        branches=[Branch(target="all_zero"),
                  Branch(target = "not_all_zero")]
    )
    protocol.add_node(root_node)

    #condition s_1 == 0 
    condition_s_1_all_zero = Condition(
        cond_type="equal",
        left="s_1",
        right=0
    )
    
    #condition flag are all zero
    condition_f_1_all_zero = Condition(cond_type="equal", left="f_1", right=0)
    
    #node all zero in first round
    round_1_all_zero = Node("all_zero",[], [])
    
    protocol.add_node(round_1_all_zero)
    
    #node not all zero in first round
    round_1_not_all_zero = Node("not_all_zero" ,["flagged_syndrome"],branches=[Branch(target="flag_all_zero"),
                  Branch(target = "flag_not_all_zero")])

    protocol.add_node(round_1_not_all_zero)
    return protocol