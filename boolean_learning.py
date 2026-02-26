from proof_protocol import *
from protocol import *
from flag_analysis import *
from circuit_op import *
start_node, protocol = load_protocol("protocols/origin_5_1_3_protocol.json")

config = read_config("[[5,1,3]]_origin/[[5,1,3]]_origin_config.txt")
gen = load_symplectic_txt(str(config['stab_txt_path']))

init_state = new_clean_circuit_state(len(gen[0][0]))
t = 1

# Output directory for per-path CNF files (path_1/, path_2/, ... with Xi.cnf, Zi.cnf, round_*/ancX.cnf, ancZ.cnf, condition*.cnf)
path = proof_protocol_boolean(protocol, start_node, init_state, config, t, cnf_output_dir="path_blif_output")
