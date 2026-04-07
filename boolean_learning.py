from proof_protocol import *
from protocol import *
from flag_analysis import *
from circuit_op import *
import argparse
import re
import sys
from pathlib import Path

# [[5,1,3]] origin
# origin_dir = "[[5,1,3]]_origin"
# protocol_path = "protocols/origin_5_1_3_protocol.json"
# config_path = f"{origin_dir}/[[5,1,3]]_origin_config.txt"
# at_most_t_faults = 1

# [[19,1,5]]_[2,2,2,1,1,1]_T
origin_dir = "[[19,1,5]]_[2,2,2,1,1,1]_T"
protocol_path = "protocols/d_5_lai_protocol.json"
config_path = f"{origin_dir}/[[19,1,5]]_[2,2,2,1,1,1]_T_lai_d_5_protocol_config.txt"
at_most_t_faults = 2
enable_decoder = True

ap = argparse.ArgumentParser(description="Run boolean learning and save console output to file.")
ap.add_argument(
    "--output-dir",
    type=Path,
    default=Path(origin_dir),
    help="Directory for per-path output logs",
)
args = ap.parse_args()

args.output_dir.mkdir(parents=True, exist_ok=True)


class PathSplitWriter:
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self._line_buf = ""
        self._current_file = None
        self._path_header_re = re.compile(r"^--- PATH (\d+) ---\s*$")

    def write(self, s: str) -> int:
        if not s:
            return 0
        self._line_buf += s
        while True:
            nl_idx = self._line_buf.find("\n")
            if nl_idx < 0:
                break
            line = self._line_buf[: nl_idx + 1]
            self._line_buf = self._line_buf[nl_idx + 1 :]
            self._consume_line(line)
        return len(s)

    def _consume_line(self, line: str) -> None:
        m = self._path_header_re.match(line.strip())
        if m:
            path_no = m.group(1)
            if self._current_file is not None:
                self._current_file.flush()
                self._current_file.close()
            out_path = self.output_dir / f"boolean_learning_output_C_path_{path_no}.txt"
            self._current_file = open(out_path, "w", encoding="utf-8")
        if self._current_file is not None:
            self._current_file.write(line)
            self._current_file.flush()

    def flush(self) -> None:
        if self._current_file is not None:
            self._current_file.flush()

    def close(self) -> None:
        if self._line_buf:
            self._consume_line(self._line_buf)
            self._line_buf = ""
        if self._current_file is not None:
            self._current_file.flush()
            self._current_file.close()
            self._current_file = None


_stdout_orig = sys.stdout
_stderr_orig = sys.stderr
_path_writer = PathSplitWriter(args.output_dir)
sys.stdout = _path_writer
sys.stderr = _path_writer

try:
    start_node, protocol = load_protocol(protocol_path)

    config = read_config(config_path)
    gen = load_symplectic_txt(str(config['stab_txt_path']))

    init_state = new_clean_circuit_state(len(gen[0][0]))
    path = proof_protocol_boolean(
        protocol,
        start_node,
        init_state,
        config,
        at_most_t_faults,
        origin_dir=origin_dir,
        enable_decoder=enable_decoder,
    )
finally:
    sys.stdout.flush()
    sys.stderr.flush()
    sys.stdout = _stdout_orig
    sys.stderr = _stderr_orig
    _path_writer.close()
