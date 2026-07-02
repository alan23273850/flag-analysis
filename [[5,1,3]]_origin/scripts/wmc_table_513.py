"""
Driver for the projected weighted model counter pipeline (wmc_513).

Mirrors monte_carlo_table_513.py's CLI shape and CSV output layout, but
computes `ratio_fail` exactly (modulo the at-most-K cardinality truncation)
via `gpmc -mode=3` instead of Monte Carlo sampling.

Output schema (CSV columns):
  p,beta,gamma,max_faults,
  ratio_fail,
  path_t0_fail,path_t1_fail,path_t2_fail,path_t3_fail,
  total_vars,total_clauses,elapsed_sec

`max_faults` defaults to 2 (sufficient for distance-3 [[5,1,3]] at typical
physical error rates). Set to None / 0 / a larger value via CLI as needed.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import wmc_513 as wmc  # noqa: E402


ORIGIN_DIR = wmc.ORIGIN_DIR
OUT_DIR = ORIGIN_DIR / "wmc_results"

DEFAULT_P_POINTS = (1e-5, 2e-5, 5e-5, 1e-4)
DEFAULT_MAX_FAULTS = 2
DEFAULT_BETA = 1.0
DEFAULT_GAMMA = 1.0


@dataclass
class PointRow:
    p: float
    beta: float
    gamma: float
    max_faults: Optional[int]
    ratio_fail: float
    path_fail: List[float]
    total_vars: int
    total_clauses: int
    elapsed_sec: float


def _p_token(p: float) -> str:
    return f"{p:.8g}".replace("+", "")


def _output_path(p_points: List[float], suffix_extra: str) -> Path:
    p_min = min(p_points)
    p_max = max(p_points)
    suffix = f"_p{_p_token(p_min)}-p{_p_token(p_max)}{suffix_extra}"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUT_DIR / f"wmc_table_513{suffix}.csv"


def _run_one_point(
    p: float,
    *,
    beta: float,
    gamma: float,
    max_faults: Optional[int],
    work_dir: Path,
    keep_cnf: bool,
) -> PointRow:
    res = wmc.compute_ratio_fail(
        p=p,
        gamma=gamma,
        beta=beta,
        work_dir=work_dir,
        keep_cnf=keep_cnf,
        max_faults=max_faults,
    )
    return PointRow(
        p=p,
        beta=beta,
        gamma=gamma,
        max_faults=max_faults,
        ratio_fail=res.ratio_fail,
        path_fail=[pp.fail_prob for pp in res.paths],
        total_vars=sum(pp.n_vars for pp in res.paths),
        total_clauses=sum(pp.n_clauses for pp in res.paths),
        elapsed_sec=res.total_elapsed_sec,
    )


def _write_outputs(rows: List[PointRow], csv_path: Path) -> None:
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "p",
                "beta",
                "gamma",
                "max_faults",
                "ratio_fail",
                "path_t0_fail",
                "path_t1_fail",
                "path_t2_fail",
                "path_t3_fail",
                "total_vars",
                "total_clauses",
                "elapsed_sec",
            ]
        )
        for r in rows:
            mf = "" if r.max_faults is None else str(r.max_faults)
            paths = r.path_fail + [float("nan")] * (4 - len(r.path_fail))
            w.writerow(
                [
                    f"{r.p:.8g}",
                    f"{r.beta:.8g}",
                    f"{r.gamma:.8g}",
                    mf,
                    f"{r.ratio_fail:.8g}",
                    f"{paths[0]:.8g}",
                    f"{paths[1]:.8g}",
                    f"{paths[2]:.8g}",
                    f"{paths[3]:.8g}",
                    r.total_vars,
                    r.total_clauses,
                    f"{r.elapsed_sec:.3f}",
                ]
            )


def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="WMC table generator for [[5,1,3]] (gpmc-based)"
    )
    ap.add_argument(
        "--p-points",
        type=float,
        nargs="+",
        default=list(DEFAULT_P_POINTS),
        help="List of physical error rates to evaluate.",
    )
    ap.add_argument("--beta", type=float, default=DEFAULT_BETA)
    ap.add_argument("--gamma", type=float, default=DEFAULT_GAMMA)
    ap.add_argument(
        "--max-faults",
        type=int,
        default=DEFAULT_MAX_FAULTS,
        help="At-most-K cardinality truncation on number of fired fault sites.",
    )
    ap.add_argument(
        "--no-max-faults",
        action="store_true",
        help="Disable the cardinality truncation (exact, but much slower).",
    )
    ap.add_argument(
        "--processes",
        "-j",
        type=int,
        default=1,
        help="Parallel processes (one WMC point per worker).",
    )
    ap.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help="Directory to write CNF files (defaults to a temp dir per run).",
    )
    ap.add_argument(
        "--keep-cnf",
        action="store_true",
        help="Keep CNF files (default: drop after parsing).",
    )
    return ap.parse_args()


def main() -> None:
    args = _parse_args()
    max_faults = None if args.no_max_faults else args.max_faults
    p_points: List[float] = list(args.p_points)

    base_work_dir = args.work_dir
    if base_work_dir is None:
        base_work_dir = OUT_DIR / "_cnf_cache"
    base_work_dir.mkdir(parents=True, exist_ok=True)

    suffix_extra = f"_K{max_faults}" if max_faults is not None else "_Knone"
    csv_path = _output_path(p_points, suffix_extra)

    start = time.time()
    rows: List[PointRow] = []
    if args.processes <= 1:
        for p in p_points:
            sub_dir = base_work_dir / f"p_{_p_token(p)}"
            sub_dir.mkdir(parents=True, exist_ok=True)
            print(
                f"[WMC] start p={p:.8g} max_faults={max_faults} "
                f"work_dir={sub_dir}"
            )
            row = _run_one_point(
                p,
                beta=args.beta,
                gamma=args.gamma,
                max_faults=max_faults,
                work_dir=sub_dir,
                keep_cnf=args.keep_cnf,
            )
            rows.append(row)
            paths_str = ", ".join(f"t{i}={v:.4g}" for i, v in enumerate(row.path_fail))
            print(
                f"[WMC] done  p={p:.8g} ratio_fail={row.ratio_fail:.6g} "
                f"({paths_str}) elapsed={row.elapsed_sec:.2f}s"
            )
    else:
        def _job(p: float) -> PointRow:
            sub_dir = base_work_dir / f"p_{_p_token(p)}"
            sub_dir.mkdir(parents=True, exist_ok=True)
            return _run_one_point(
                p,
                beta=args.beta,
                gamma=args.gamma,
                max_faults=max_faults,
                work_dir=sub_dir,
                keep_cnf=args.keep_cnf,
            )

        with cf.ProcessPoolExecutor(max_workers=args.processes) as ex:
            futures = {ex.submit(_job, p): p for p in p_points}
            for fut in cf.as_completed(futures):
                row = fut.result()
                rows.append(row)
                paths_str = ", ".join(f"t{i}={v:.4g}" for i, v in enumerate(row.path_fail))
                print(
                    f"[WMC] done p={row.p:.8g} ratio_fail={row.ratio_fail:.6g} "
                    f"({paths_str}) elapsed={row.elapsed_sec:.2f}s"
                )

    rows.sort(key=lambda r: r.p)
    _write_outputs(rows, csv_path)
    print(f"[WMC] wrote: {csv_path}")
    print(f"[WMC] total wall: {time.time() - start:.2f}s")


if __name__ == "__main__":
    main()
