from __future__ import annotations
import io
import sys
from pathlib import Path

from proof_protocol import *
from protocol import *
from flag_analysis import *
from circuit_op import *

_PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_PATH = _PROJECT_ROOT / "output_[[5,1,3]]_origin" / "output.txt"


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


def main():
    start_node, protocol = load_protocol("protocols/origin_5_1_3_protocol.json")
    config = read_config("[[5,1,3]]_origin/[[5,1,3]]_origin_config.txt")
    gen = load_symplectic_txt(str(config['stab_txt_path']))
    init_state = new_clean_circuit_state(len(gen[0][0]))
    t = 1
    path = proof_protocol_boolean(protocol, start_node, init_state, config, t, cnf_output_dir="output_[[5,1,3]]_origin")
    return path


if __name__ == "__main__":
    buf = io.StringIO()
    orig_stdout = sys.stdout
    sys.stdout = TeeWriter(orig_stdout, buf)
    try:
        main()
    finally:
        sys.stdout = orig_stdout
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(buf.getvalue(), encoding="utf-8")
    print(f"(Output also written to {OUTPUT_PATH})", file=orig_stdout)
