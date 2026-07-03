#!/usr/bin/env python3
"""Sweep decreasing p (K=2) until ratio_fail <= p."""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import wmc_913 as wmc  # noqa: E402

OUT_DIR = wmc.ORIGIN_DIR / "wmc_results"


def _p_token(p: float) -> str:
    return f"{p:.8g}".replace("+", "")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start-p", type=float, default=1e-7)
    ap.add_argument("--factor", type=float, default=0.1, help="multiply p each step")
    ap.add_argument("--min-p", type=float, default=1e-20)
    ap.add_argument("--max-faults", type=int, default=2)
    ap.add_argument("--keep-cnf", action="store_true")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "wmc_sweep_small_p_K2.csv"
    txt_path = OUT_DIR / "wmc_sweep_small_p_K2.txt"

    rows: list[dict] = []
    if csv_path.is_file():
        with csv_path.open(newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))

    p = args.start_p
    crossed = False
    while p >= args.min_p:
        already = {float(r["p"]) for r in rows}
        if p in already or any(abs(float(r["p"]) - p) / max(p, 1e-30) < 1e-3 for r in rows):
            print(f"[skip] p={p:.8g} already in table")
            p *= args.factor
            continue

        work_dir = OUT_DIR / "_cnf_cache" / f"sweep_p_{_p_token(p)}"
        work_dir.mkdir(parents=True, exist_ok=True)
        print(f"[sweep] start p={p:.8g} max_faults={args.max_faults}", flush=True)
        t0 = time.time()
        res = wmc.compute_ratio_fail(
            p=p,
            max_faults=args.max_faults,
            work_dir=work_dir,
            keep_cnf=args.keep_cnf,
        )
        rf = res.ratio_fail
        ratio_over_p = rf / p if p > 0 else float("inf")
        row = {
            "p": f"{p:.8g}",
            "ratio_fail": f"{rf:.8g}",
            "ratio_over_p": f"{ratio_over_p:.8g}",
            "path_0": f"{res.paths[0].fail_prob:.8g}",
            "path_1": f"{res.paths[1].fail_prob:.8g}",
            "path_2": f"{res.paths[2].fail_prob:.8g}",
            "path_3": f"{res.paths[3].fail_prob:.8g}",
            "elapsed_sec": f"{time.time() - t0:.1f}",
        }
        rows.append(row)
        print(
            f"[sweep] done  p={p:.8g} ratio_fail={rf:.8g} ratio/p={ratio_over_p:.6g} "
            f"elapsed={row['elapsed_sec']}s",
            flush=True,
        )

        with csv_path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(
                f,
                fieldnames=[
                    "p",
                    "ratio_fail",
                    "ratio_over_p",
                    "path_0",
                    "path_1",
                    "path_2",
                    "path_3",
                    "elapsed_sec",
                ],
            )
            w.writeheader()
            w.writerows(rows)

        txt_path.write_text(
            "p,ratio_fail,ratio_over_p,path_0,path_1,path_2,path_3,elapsed_sec\n"
            + "\n".join(
                f"{r['p']},{r['ratio_fail']},{r['ratio_over_p']},"
                f"{r['path_0']},{r['path_1']},{r['path_2']},{r['path_3']},{r['elapsed_sec']}"
                for r in rows
            )
            + "\n",
            encoding="utf-8",
        )

        if rf <= p:
            print(f"[sweep] STOP: ratio_fail <= p at p={p:.8g}", flush=True)
            crossed = True
            break

        p *= args.factor

    if not crossed:
        print(f"[sweep] reached min_p={args.min_p:.8g} without ratio_fail <= p", flush=True)
    print(f"[sweep] wrote {csv_path}", flush=True)


if __name__ == "__main__":
    main()
