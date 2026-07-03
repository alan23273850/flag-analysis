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

from decoder_commute_verify import decoder_path, preload_decoder_assets, verify_decoder_commute_913

SIM_PATH = Path(__file__).resolve().parent / "sim_913.py"
ORIGIN_DIR = Path(__file__).resolve().parents[1]
STAB_TXT = ORIGIN_DIR / "[[9,1,3]].txt"
LOG_TXT = ORIGIN_DIR / "[[9,1,3]]_log_op.txt"
OUT_DIR = ORIGIN_DIR / "mc_results"

TARGET_FAILS = 800
BETA = 1.0
GAMMA = 1.0
P_POINTS = [1e-4, 1e-5, 1e-6]
DEFAULT_PROCESSES = 64
DEFAULT_CHUNK_RUNS = 100

_SIM_MOD_WORKER: Any | None = None


@dataclass
class PointResult:
    p: float
    beta: float
    gamma: float
    target_fails: int
    decoder_pass: int
    decoder_fail: int
    ratio_fail: float
    total_runs: int
    elapsed_sec: float


def _p_token(p: float) -> str:
    return f"{p:.8g}".replace("+", "")


def _backend_name() -> str:
    import os
    return (os.environ.get("DECODER_BACKEND", "cms") or "cms").lower()


def _output_paths() -> tuple[Path, Path]:
    p_min = min(P_POINTS)
    p_max = max(P_POINTS)
    suffix = f"_p{_p_token(p_min)}-p{_p_token(p_max)}_{_backend_name()}"
    return (
        OUT_DIR / f"mc_table_913{suffix}.csv",
        OUT_DIR / f"mc_table_913{suffix}.txt",
    )


def _load_sim_module(path: Path) -> Any:
    spec = importlib.util.spec_from_file_location("sim913_mod", str(path.resolve()))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load simulator module from {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _worker_init(sim_path_str: str) -> None:
    global _SIM_MOD_WORKER
    _SIM_MOD_WORKER = _load_sim_module(Path(sim_path_str))
    _SIM_MOD_WORKER.preload_static_inputs()
    preload_decoder_assets(
        decoder_paths=[decoder_path(ORIGIN_DIR, i) for i in range(4)],
        log_path=LOG_TXT,
        stab_path=STAB_TXT,
    )


def _run_chunk(sim_mod: Any, p: float, runs: int) -> tuple[int, int, int]:
    sim_mod.TWO_QUBIT_FAULT_P = p
    sim_mod.IDLE_GAMMA = GAMMA
    sim_mod.MEAS_BETA = BETA

    dec_pass = 0
    dec_fail = 0
    total = 0

    for _ in range(runs):
        total += 1
        rec, data_x, data_z, _ = sim_mod.run_protocol_once()
        dec_file = decoder_path(ORIGIN_DIR, rec.path_idx)
        ok, _ = verify_decoder_commute_913(
            meas_assign=rec.meas_assign,
            data_x=data_x or [],
            data_z=data_z or [],
            decoder_c_path=dec_file,
            log_path=LOG_TXT,
            stab_path=STAB_TXT,
        )
        if ok:
            dec_pass += 1
        else:
            dec_fail += 1
    return dec_pass, dec_fail, total


def _worker_run_chunk(p: float, runs: int) -> tuple[int, int, int]:
    if _SIM_MOD_WORKER is None:
        raise RuntimeError("Worker simulator module not initialized")
    return _run_chunk(_SIM_MOD_WORKER, p, runs)


def _run_one_point_parallel(p: float, processes: int, chunk_runs: int) -> PointResult:
    start = time.time()
    dec_pass = 0
    dec_fail = 0
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
                c_pass, c_fail, c_total = fut.result()
                dec_pass += c_pass
                dec_fail += c_fail
                total += c_total

                if dec_fail // 100 > last_report_fails // 100:
                    last_report_fails = dec_fail
                    print(
                        f"[MC913] p={p:.8g} fail={dec_fail}/{TARGET_FAILS} "
                        f"pass={dec_pass} total={total} ratio={dec_fail/total:.6g}"
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
                "ratio_fail",
                "ratio_over_p",
                "total_runs",
                "elapsed_sec",
            ]
        )
        for r in results:
            ratio_over_p = r.ratio_fail / r.p if r.p > 0 else 0.0
            writer.writerow(
                [
                    f"{r.p:.8g}",
                    f"{r.beta:.8g}",
                    f"{r.gamma:.8g}",
                    r.target_fails,
                    r.decoder_pass,
                    r.decoder_fail,
                    f"{r.ratio_fail:.8g}",
                    f"{ratio_over_p:.8g}",
                    r.total_runs,
                    f"{r.elapsed_sec:.3f}",
                ]
            )

    lines = [
        "Monte Carlo table for [[9,1,3]] Chris d=3 protocol",
        f"Target fails per point: {TARGET_FAILS}",
        f"beta={BETA}, gamma={GAMMA}",
        "ratio_fail = decoder_fail / total_runs",
        "ratio_over_p = ratio_fail / p  (compare to PWMC ~249)",
        "",
        "p,beta,gamma,target_fails,decoder_pass,decoder_fail,ratio_fail,ratio_over_p,total_runs,elapsed_sec",
    ]
    for r in results:
        ratio_over_p = r.ratio_fail / r.p if r.p > 0 else 0.0
        lines.append(
            f"{r.p:.8g},{r.beta:.8g},{r.gamma:.8g},{r.target_fails},{r.decoder_pass},"
            f"{r.decoder_fail},{r.ratio_fail:.8g},{ratio_over_p:.8g},{r.total_runs},{r.elapsed_sec:.3f}"
        )
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return csv_path, txt_path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Monte Carlo table for [[9,1,3]]")
    parser.add_argument("-j", "--processes", type=int, default=DEFAULT_PROCESSES)
    parser.add_argument("--chunk-runs", type=int, default=DEFAULT_CHUNK_RUNS)
    parser.add_argument("--target-fails", type=int, default=TARGET_FAILS)
    parser.add_argument(
        "--backend",
        choices=("cms", "dnf", "cadet"),
        default=None,
        help="Decoder backend subdirectory (default: DECODER_BACKEND env or cms)",
    )
    return parser.parse_args()


def main() -> None:
    global TARGET_FAILS
    args = _parse_args()
    if args.backend is not None:
        import os
        os.environ["DECODER_BACKEND"] = args.backend
    TARGET_FAILS = args.target_fails
    if args.processes < 1:
        raise ValueError("--processes must be >= 1")

    results: list[PointResult] = []
    for p in P_POINTS:
        print(f"[MC913] start p={p:.8g} ({args.processes} workers)")
        r = _run_one_point_parallel(p, processes=args.processes, chunk_runs=args.chunk_runs)
        rop = r.ratio_fail / r.p if r.p > 0 else 0.0
        results.append(r)
        print(
            f"[MC913] done p={r.p:.8g} ratio_fail={r.ratio_fail:.8g} "
            f"ratio/p={rop:.6g} total={r.total_runs} elapsed={r.elapsed_sec:.1f}s"
        )

    csv_path, txt_path = _write_outputs(results)
    print(f"[MC913] wrote: {csv_path}")
    print(f"[MC913] wrote: {txt_path}")


if __name__ == "__main__":
    main()
