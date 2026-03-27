from __future__ import annotations

import argparse
import csv
import concurrent.futures as cf
import importlib.util
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from decoder_commute_verify import (
    decoder_file_index,
    preload_decoder_assets,
    verify_decoder_commute,
)


SIM_PATH = Path(__file__).resolve().parent / "sim_513.py"
ORIGIN_DIR = PROJECT_ROOT / "[[5,1,3]]_origin"
STAB_TXT = ORIGIN_DIR / "[[5,1,3]]_origin.txt"
LOG_TXT = ORIGIN_DIR / "[[5,1,3]]_log_op.txt"

OUT_DIR = ORIGIN_DIR / "mc_results"

TARGET_FAILS = 1600
BETA = 1.0
GAMMA = 1.0
P_POINTS = [0.00014 + i * 1e-6 for i in range(1, 5)]  # 1e-5 ... 1e-4
DEFAULT_PROCESSES = 128
DEFAULT_CHUNK_RUNS = 200


@dataclass
class PointResult:
    p: float
    beta: float
    gamma: float
    target_fails: int
    decoder_pass: int
    decoder_fail: int
    ratio_fail: float
    no_lut_runs: int
    total_runs: int
    elapsed_sec: float


_SIM_MOD_WORKER: Any | None = None


def _p_token(p: float) -> str:
    # Keep filename readable while supporting scientific notation.
    return f"{p:.8g}".replace("+", "")


def _output_paths() -> tuple[Path, Path]:
    p_min = min(P_POINTS)
    p_max = max(P_POINTS)
    suffix = f"_p{_p_token(p_min)}-p{_p_token(p_max)}"
    return (
        OUT_DIR / f"mc_table_513{suffix}.csv",
        OUT_DIR / f"mc_table_513{suffix}.txt",
    )


def _load_sim_module(path: Path) -> Any:
    spec = importlib.util.spec_from_file_location("sim513_mod", str(path.resolve()))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load simulator module from {path}")
    mod = importlib.util.module_from_spec(spec)
    # Ensure dataclass/type resolution can find this module by name.
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _worker_init(sim_path_str: str) -> None:
    global _SIM_MOD_WORKER
    _SIM_MOD_WORKER = _load_sim_module(Path(sim_path_str))
    _SIM_MOD_WORKER.preload_static_inputs()
    preload_decoder_assets(
        decoder_paths=[ORIGIN_DIR / "decoder" / f"path_{i}.txt" for i in range(1, 5)],
        log_path=LOG_TXT,
        stab_path=STAB_TXT,
    )


def _run_chunk(sim_mod: Any, p: float, runs: int) -> tuple[int, int, int, int]:
    sim_mod.TWO_QUBIT_FAULT_P = p
    sim_mod.IDLE_GAMMA = GAMMA
    sim_mod.MEAS_BETA = BETA

    dec_pass = 0
    dec_fail = 0
    no_lut = 0
    total = 0

    for _ in range(runs):
        total += 1
        rec, data_x, data_z, _fault_events = sim_mod.run_protocol_once()
        if rec is None:
            no_lut += 1
            continue

        dec_idx = decoder_file_index(rec.first_stabilizer_index)
        dec_path = ORIGIN_DIR / "decoder" / f"path_{dec_idx}.txt"
        ok, _ = verify_decoder_commute(
            first_stabilizer_index=rec.first_stabilizer_index,
            bitstring6=rec.bitstring6,
            data_x=data_x or [],
            data_z=data_z or [],
            decoder_c_path=dec_path,
            log_path=LOG_TXT,
            stab_path=STAB_TXT,
        )
        if ok:
            dec_pass += 1
        else:
            dec_fail += 1
    return dec_pass, dec_fail, no_lut, total


def _worker_run_chunk(p: float, runs: int) -> tuple[int, int, int, int]:
    if _SIM_MOD_WORKER is None:
        raise RuntimeError("Worker simulator module not initialized")
    return _run_chunk(_SIM_MOD_WORKER, p, runs)


def _run_one_point_serial(sim_mod: Any, p: float) -> PointResult:
    start = time.time()
    dec_pass = 0
    dec_fail = 0
    no_lut = 0
    total = 0
    last_report_fails = -1

    while dec_fail < TARGET_FAILS:
        c_pass, c_fail, c_no_lut, c_total = _run_chunk(sim_mod, p, runs=DEFAULT_CHUNK_RUNS)
        dec_pass += c_pass
        dec_fail += c_fail
        no_lut += c_no_lut
        total += c_total

        if dec_fail // 100 > last_report_fails // 100:
            last_report_fails = dec_fail
            print(
                f"[MC] p={p:.8g} progress fail={dec_fail}/{TARGET_FAILS} "
                f"pass={dec_pass} no_lut={no_lut} total={total}"
            )

    elapsed = time.time() - start
    # ratio_fail = decoder_fail / total_runs (includes no-LUT runs; those are not decoder_fail)
    ratio = dec_fail / total if total > 0 else 0.0
    return PointResult(
        p=p,
        beta=BETA,
        gamma=GAMMA,
        target_fails=TARGET_FAILS,
        decoder_pass=dec_pass,
        decoder_fail=dec_fail,
        ratio_fail=ratio,
        no_lut_runs=no_lut,
        total_runs=total,
        elapsed_sec=elapsed,
    )


def _run_one_point_parallel(p: float, processes: int, chunk_runs: int) -> PointResult:
    start = time.time()
    dec_pass = 0
    dec_fail = 0
    no_lut = 0
    total = 0
    last_report_fails = -1

    with cf.ProcessPoolExecutor(
        max_workers=processes,
        initializer=_worker_init,
        initargs=(str(SIM_PATH.resolve()),),
    ) as ex:
        while dec_fail < TARGET_FAILS:
            futures = [ex.submit(_worker_run_chunk, p, chunk_runs) for _ in range(processes)]
            for fut in cf.as_completed(futures):
                c_pass, c_fail, c_no_lut, c_total = fut.result()
                dec_pass += c_pass
                dec_fail += c_fail
                no_lut += c_no_lut
                total += c_total

                if dec_fail // 100 > last_report_fails // 100:
                    last_report_fails = dec_fail
                    print(
                        f"[MC] p={p:.8g} progress fail={dec_fail}/{TARGET_FAILS} "
                        f"pass={dec_pass} no_lut={no_lut} total={total}"
                    )

                if total % 20000 < c_total:
                    print(
                        f"[MC] p={p:.8g} heartbeat total={total} "
                        f"fail={dec_fail} pass={dec_pass} no_lut={no_lut}"
                    )

    elapsed = time.time() - start
    ratio = dec_fail / total if total > 0 else 0.0
    return PointResult(
        p=p,
        beta=BETA,
        gamma=GAMMA,
        target_fails=TARGET_FAILS,
        decoder_pass=dec_pass,
        decoder_fail=dec_fail,
        ratio_fail=ratio,
        no_lut_runs=no_lut,
        total_runs=total,
        elapsed_sec=elapsed,
    )


def _write_outputs(results: list[PointResult]) -> tuple[Path, Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path, txt_path = _output_paths()

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "p",
                "beta",
                "gamma",
                "target_fails",
                "decoder_pass",
                "decoder_fail",
                "ratio_fail",  # decoder_fail / total_runs
                "no_lut_runs",
                "total_runs",
                "elapsed_sec",
            ]
        )
        for r in results:
            writer.writerow(
                [
                    f"{r.p:.8g}",
                    f"{r.beta:.8g}",
                    f"{r.gamma:.8g}",
                    r.target_fails,
                    r.decoder_pass,
                    r.decoder_fail,
                    f"{r.ratio_fail:.8g}",
                    r.no_lut_runs,
                    r.total_runs,
                    f"{r.elapsed_sec:.3f}",
                ]
            )

    lines: list[str] = []
    lines.append("Monte Carlo table for [[5,1,3]] two-qubit simulator")
    lines.append(f"Target fails per point: {TARGET_FAILS}")
    lines.append(f"beta={BETA}, gamma={GAMMA}")
    lines.append("ratio_fail = decoder_fail / total_runs  (all protocol runs, incl. no LUT)")
    lines.append("")
    lines.append(
        "p,beta,gamma,target_fails,decoder_pass,decoder_fail,ratio_fail,no_lut_runs,total_runs,elapsed_sec"
    )
    for r in results:
        lines.append(
            f"{r.p:.8g},{r.beta:.8g},{r.gamma:.8g},{r.target_fails},{r.decoder_pass},"
            f"{r.decoder_fail},{r.ratio_fail:.8g},{r.no_lut_runs},{r.total_runs},{r.elapsed_sec:.3f}"
        )
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return csv_path, txt_path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Monte Carlo table generator for [[5,1,3]]")
    parser.add_argument(
        "--processes",
        "-j",
        type=int,
        default=DEFAULT_PROCESSES,
        help=f"Worker processes per p point (default: {DEFAULT_PROCESSES})",
    )
    parser.add_argument(
        "--chunk-runs",
        type=int,
        default=DEFAULT_CHUNK_RUNS,
        help=f"Trials per worker task (default: {DEFAULT_CHUNK_RUNS})",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    if args.processes < 1:
        raise ValueError("--processes must be >= 1")
    if args.chunk_runs < 1:
        raise ValueError("--chunk-runs must be >= 1")

    use_parallel = args.processes > 1
    sim_mod = None if use_parallel else _load_sim_module(SIM_PATH)
    if sim_mod is not None:
        sim_mod.preload_static_inputs()
    preload_decoder_assets(
        decoder_paths=[ORIGIN_DIR / "decoder" / f"path_{i}.txt" for i in range(1, 5)],
        log_path=LOG_TXT,
        stab_path=STAB_TXT,
    )
    results: list[PointResult] = []

    for p in P_POINTS:
        mode = f"{args.processes} proc x chunk {args.chunk_runs}" if use_parallel else "serial"
        print(f"[MC] start p={p:.8g} mode={mode}")
        if use_parallel:
            r = _run_one_point_parallel(p, processes=args.processes, chunk_runs=args.chunk_runs)
        else:
            r = _run_one_point_serial(sim_mod, p)
        results.append(r)
        print(
            f"[MC] done p={p:.8g} ratio={r.ratio_fail:.8g} "
            f"(pass={r.decoder_pass}, fail={r.decoder_fail}, no_lut={r.no_lut_runs}, total={r.total_runs})"
        )

    csv_path, txt_path = _write_outputs(results)
    print(f"[MC] wrote: {csv_path}")
    print(f"[MC] wrote: {txt_path}")


if __name__ == "__main__":
    main()

