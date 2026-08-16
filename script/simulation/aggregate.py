#!/usr/bin/env python3
"""
Aggregate per-simulation result JSONs into a single CSV.

Usage:
  python3 aggregate.py --mode one_wgd
  python3 aggregate.py --mode one_wgd --name myrun
"""
import argparse
import csv
import glob
import json
import os
import re
import sys

from run_single import MODES

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def main():
    parser = argparse.ArgumentParser(description="Aggregate simulation results into a CSV")
    parser.add_argument("--mode", required=True, choices=list(MODES.keys()))
    parser.add_argument("--name", default="", help="Label used when running run_single.py")
    args = parser.parse_args()

    cfg    = MODES[args.mode]
    suffix = f"_{args.name}" if args.name else ""

    sim_outputs = os.path.join(SCRIPT_DIR, cfg["sim_outputs_dir"] + suffix)
    pattern     = os.path.join(sim_outputs, "sim_*/result.json")
    files       = sorted(
        glob.glob(pattern),
        key=lambda p: int(re.search(r"sim_(\d+)", p).group(1)),
    )

    if not files:
        print(f"No result files found matching: {pattern}", file=sys.stderr)
        sys.exit(1)

    results = []
    errors  = []
    for f in files:
        with open(f) as fp:
            r = json.load(fp)
        if "error" in r:
            errors.append(r["sim_id"])
        results.append(r)

    if errors:
        print(f"Warning: {len(errors)} simulation(s) had errors: {errors}", file=sys.stderr)

    # CSV column order from the first successful result
    first_ok = next((r for r in results if "error" not in r), results[0])
    fields   = list(first_ok.keys())

    out_csv = os.path.join(SCRIPT_DIR, cfg["results_csv"].replace(".csv", f"{suffix}.csv"))
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            writer.writerow({k: "" if v is None else v for k, v in r.items()})

    print(f"Aggregated {len(results)} simulations ({len(errors)} errors) → {out_csv}")


if __name__ == "__main__":
    main()
