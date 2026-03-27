from proof_protocol import *
from protocol import *
from flag_analysis import *
from circuit_op import *

origin_dir = "[[5,1,3]]_origin"
protocol_path = "protocols/origin_5_1_3_protocol.json"
config_path = f"{origin_dir}/[[5,1,3]]_origin_config.txt"

start_node, protocol = load_protocol(protocol_path)

config = read_config(config_path)
gen = load_symplectic_txt(str(config['stab_txt_path']))

init_state = new_clean_circuit_state(len(gen[0][0]))
t = 1

path = proof_protocol_boolean(protocol, start_node, init_state, config, t, origin_dir=origin_dir)
