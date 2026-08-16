#!/usr/bin/env python3
"""
Single-simulation runner — designed to be called as a Slurm array task.

Usage:
  python3 run_single.py --mode one_wgd --sim_id 1
  python3 run_single.py --mode one_wgd --sim_id $SLURM_ARRAY_TASK_ID --name myrun

Each invocation runs one simulation + inference and writes its result to:
  sim_outputs_<mode>[_<name>]/sim_<sim_id>/result.json

Run aggregate.py afterwards to collect all results into a single CSV.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import time

# ── Mode configurations ────────────────────────────────────────────────────────

MODES = {
    "standard": {
        "description":     "Rate parameter recovery across multiple trees",
        "sim_trees": [
            "test_data/tiley2016/Brassicaceae/tree.newick",
            "test_data/tiley2016/Eudicots/tree.newick",
            "test_data/tiley2016/LandPlants/tree.newick",
            "test_data/tiley2016/Monocots/tree.newick",
        ],
        "infer_tree":      None,
        "true_qs":         [],
        "root_lambda":     -1,
        "sim_root_lambda": 1.5,
        "rate_init":       "generic",
        "wgd_threshold":   10,
        "target_events":   6.0,
        "rates": {
            "gain":        1.4,
            "loss":        1.2,
            "innovation":  0.7,
            "elimination": 0.5,
        },
        "results_csv":     "simulation_study_results.csv",
        "sim_outputs_dir": "sim_outputs_standard",
    },
    "one_wgd": {
        "description":     "WGD detection — single WGD event",
        "sim_trees": [
            "test_data/tiley2016/Monocots/tree_wgd.newick",
        ],
        "infer_tree":      "test_data/tiley2016/Monocots/tree.newick",
        "true_qs":         [0.3],
        "root_lambda":     -1,
        "sim_root_lambda": 1.5,
        "rate_init":       "generic",
        "wgd_threshold":   10,
        "target_events":   5.0,
        "rates": {
            "gain":        1.7,
            "loss":        1.0,
            "innovation":  0.2,
            "elimination": 0.15,
        },
        "results_csv":     "wgd_detection_results_one.csv",
        "sim_outputs_dir": "sim_outputs_one_wgd",
    },
    "two_wgds": {
        "description":     "WGD detection — two WGD events",
        "sim_trees": [
            "test_data/tiley2016/Monocots/tree_wgds.newick",
        ],
        "infer_tree":      "test_data/tiley2016/Monocots/tree.newick",
        "true_qs":         [0.5, 0.5],
        "root_lambda":     1.0,
        "rate_init":       "generic",
        "wgd_threshold":   10,
        "target_events":   5.0,
        "rates": {
            "gain":        1.3,
            "loss":        1.0,
            "innovation":  0.1,
            "elimination": 0.05,
        },
        "results_csv":     "wgd_detection_results_two.csv",
        "sim_outputs_dir": "sim_outputs_two_wgds",
    },
    "four_wgds": {
        "description":     "WGD detection — four WGD events",
        "sim_trees": [
            "test_data/tiley2016/Monocots/tree_4wgds.newick",
        ],
        "infer_tree":      "test_data/tiley2016/Monocots/tree.newick",
        "true_qs":         [0.5, 0.5, 0.5, 0.5],
        "root_lambda":     1.0,
        "rate_init":       "generic",
        "wgd_threshold":   10,
        "target_events":   5.0,
        "rates": {
            "gain":        1.3,
            "loss":        1.0,
            "innovation":  0.1,
            "elimination": 0.05,
        },
        "results_csv":     "wgd_detection_results_four.csv",
        "sim_outputs_dir": "sim_outputs_four_wgds",
    },
}

# ── Shared constants ───────────────────────────────────────────────────────────

NUM_SITES = 1000
MAX_STATE = 50

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
REPO_ROOT  = os.path.dirname(os.path.dirname(SCRIPT_DIR))

# ── Tree helpers ───────────────────────────────────────────────────────────────

def compute_tree_length(path):
    with open(path) as f:
        content = f.read()
    lengths = re.findall(r':([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)', content)
    return sum(float(x) for x in lengths)

# ── Param file writers ─────────────────────────────────────────────────────────

def write_sim_param_file(path, cfg, tree, branch_mul, seed):
    rates = cfg["rates"]
    with open(path, "w") as f:
        f.write(f"_treePath        = {tree}\n")
        f.write(f"_outputFasta     = sim_output/simulated.fasta\n")
        f.write(f"_gain            = {rates['gain']}\n")
        f.write(f"_loss            = {rates['loss']}\n")
        f.write(f"_innovation      = {rates['innovation']}\n")
        f.write(f"_elimination     = {rates['elimination']}\n")
        f.write(f"_gainFunc        = CONST\n")
        f.write(f"_lossFunc        = CONST\n")
        f.write(f"_innovationFunc  = CONST\n")
        f.write(f"_eliminationFunc = CONST\n")
        f.write(f"_maxState        = {MAX_STATE}\n")
        f.write(f"_branchMul       = {branch_mul:.6f}\n")
        sim_lambda = cfg["sim_root_lambda"] if cfg["root_lambda"] <= 0 else cfg["root_lambda"]
        f.write(f"_rootLambda      = {sim_lambda}\n")
        f.write(f"_numSites        = {NUM_SITES}\n")
        f.write(f"_seed            = {seed}\n")
        f.write(f"_wgdMode         = disable\n")
        if cfg["true_qs"]:
            f.write(f"_fixedWgdQ       = {','.join(str(q) for q in cfg['true_qs'])}\n")


def write_genevol_param_file(path, cfg, tree, branch_mul):
    rates = cfg["rates"]
    n     = len(cfg["true_qs"])
    with open(path, "w") as f:
        f.write(f"_treePath        = {tree}\n")
        f.write(f"_dataPath        = sim_output/simulated.fasta\n")
        if cfg["rate_init"] == "true":
            f.write(f"_gain            = {rates['gain']}\n")
            f.write(f"_loss            = {rates['loss']}\n")
            f.write(f"_innovation      = {rates['innovation']}\n")
            f.write(f"_elimination     = {rates['elimination']}\n")
        else:
            f.write(f"_gain            = 0.5\n")
            f.write(f"_loss            = 0.5\n")
            f.write(f"_innovation      = 0.01\n")
            f.write(f"_elimination     = 0.01\n")
        f.write(f"_gainFunc        = CONST\n")
        f.write(f"_lossFunc        = CONST\n")
        f.write(f"_innovationFunc  = CONST\n")
        f.write(f"_eliminationFunc = CONST\n")
        f.write(f"_maxState        = {MAX_STATE}\n")
        f.write(f"_branchMul       = {branch_mul:.6f}\n")
        if cfg["root_lambda"] > 0:
            f.write(f"_rootLambda      = {cfg['root_lambda']}\n")
        if n > 0:
            f.write(f"_wgdMode         = detect\n")
            f.write(f"_wgdThreshold    = {cfg['wgd_threshold']}\n")
            f.write(f"_modelCriterion  = AIC\n")
        else:
            f.write(f"_wgdMode         = disable\n")

# ── Runners ────────────────────────────────────────────────────────────────────

def run_simulator(param_file_abs, sim_output_abs):
    subprocess.run([
        os.path.join(REPO_ROOT, "Simulator/simulator"),
        f"param={param_file_abs}",
    ], check=True, cwd=sim_output_abs)


def run_genevol(param_file_abs, sim_output_abs):
    result = subprocess.run([
        os.path.join(REPO_ROOT, "GenEvol/genEvol"),
        f"param={param_file_abs}",
    ], capture_output=True, text=True, cwd=sim_output_abs)
    return result.stdout + result.stderr

# ── Output parsers ─────────────────────────────────────────────────────────────

def parse_recovered(log):
    recovered = {}
    for line in log.splitlines():
        m = re.match(r"\s+GeneCount\.(\w+?)(\d+)_\d+\s*=\s*([\d.eE+\-]+)", line)
        if m:
            name, val = m.group(1), float(m.group(3))
            if name in ("gain", "loss", "innovation", "elimination"):
                recovered[name] = val
        m2 = re.match(r"\s+Poisson\.lambda\s*=\s*([\d.eE+\-]+)", line)
        if m2:
            recovered["lambda"] = float(m2.group(1))
    return recovered


def parse_wgd_results(log):
    block_match = re.search(r"=== WGD Detection Results ===(.*?)(?:===|\Z)", log, re.DOTALL)
    if not block_match or "No WGD events detected." in block_match.group(1):
        return []
    rows = []
    for m in re.finditer(
        r"^\s+\d+\s+\d+\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)", block_match.group(1), re.MULTILINE
    ):
        rows.append({"q": float(m.group(1)), "delta_aic": float(m.group(2))})
    return rows

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Run a single simulation + inference")
    parser.add_argument("--mode",   required=True, choices=list(MODES.keys()))
    parser.add_argument("--sim_id", required=True, type=int,
                        help="Simulation index (use $SLURM_ARRAY_TASK_ID in batch jobs)")
    parser.add_argument("--name",   default="",
                        help="Optional label appended to output directory name")
    args = parser.parse_args()

    cfg       = MODES[args.mode]
    rates     = cfg["rates"]
    n_wgds    = len(cfg["true_qs"])
    opt_lambda = cfg["root_lambda"] <= 0

    sim_tree   = cfg["sim_trees"][(args.sim_id - 1) % len(cfg["sim_trees"])]
    infer_tree = cfg.get("infer_tree") or sim_tree
    tree_label = os.path.basename(os.path.dirname(sim_tree))

    raw        = compute_tree_length(os.path.join(REPO_ROOT, sim_tree))
    target_len = cfg["target_events"] / max(rates.values())
    bm         = target_len / raw

    suffix     = f"_{args.name}" if args.name else ""
    sim_output = os.path.join(SCRIPT_DIR, cfg["sim_outputs_dir"] + suffix,
                              f"sim_{args.sim_id:04d}")
    os.makedirs(sim_output, exist_ok=True)

    print(f"[sim {args.sim_id}] mode={args.mode} tree={tree_label} branchMul={bm:.6f}")

    t_start = time.time()
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            sim_param     = os.path.join(tmpdir, "sim_param.txt")
            genevol_param = os.path.join(tmpdir, "genevol_param.txt")

            write_sim_param_file(sim_param, cfg, sim_tree, bm, args.sim_id)
            write_genevol_param_file(genevol_param, cfg, infer_tree, bm)

            run_simulator(os.path.abspath(sim_param), sim_output)
            log = run_genevol(os.path.abspath(genevol_param), sim_output)

        with open(os.path.join(sim_output, "genevol_log.txt"), "w") as f:
            f.write(log)

        rec      = parse_recovered(log)
        wgd_rows = parse_wgd_results(log) if n_wgds > 0 else []
        elapsed  = time.time() - t_start

        result = {}
        result["sim_id"]  = args.sim_id
        result["seed"]    = args.sim_id
        result["mode"]    = args.mode
        result["tree"]    = tree_label

        for i, q in enumerate(cfg["true_qs"]):
            result[f"true_q_{i+1}"] = q

        result["true_gain"]        = rates["gain"]
        result["true_loss"]        = rates["loss"]
        result["true_innovation"]  = rates["innovation"]
        result["true_elimination"] = rates["elimination"]
        if opt_lambda:
            result["true_lambda"] = cfg.get("sim_root_lambda")

        if n_wgds > 0:
            result["n_detected"] = len(wgd_rows)
            for i in range(n_wgds):
                result[f"detected_q_{i+1}"] = wgd_rows[i]["q"]   if i < len(wgd_rows) else None
                result[f"delta_aic_{i+1}"]  = wgd_rows[i]["delta_aic"] if i < len(wgd_rows) else None
            result["inferred_gain"]        = rec.get("gain")
            result["inferred_loss"]        = rec.get("loss")
            result["inferred_innovation"]  = rec.get("innovation")
            result["inferred_elimination"] = rec.get("elimination")
        else:
            result["recovered_gain"]        = rec.get("gain")
            result["recovered_loss"]        = rec.get("loss")
            result["recovered_innovation"]  = rec.get("innovation")
            result["recovered_elimination"] = rec.get("elimination")
            if opt_lambda:
                result["recovered_lambda"] = rec.get("lambda")
            for p in ("gain", "loss", "innovation", "elimination"):
                t, r = rates[p], rec.get(p)
                if r is not None:
                    result[f"diff_{p}"]    = abs(r - t)
                    result[f"err_pct_{p}"] = abs(r - t) / t * 100 if t != 0 else None
                else:
                    result[f"diff_{p}"]    = None
                    result[f"err_pct_{p}"] = None

        result["runtime_sec"] = round(elapsed, 2)

    except Exception as e:
        elapsed = time.time() - t_start
        print(f"[sim {args.sim_id}] ERROR: {e}", file=sys.stderr)
        result = {"sim_id": args.sim_id, "seed": args.sim_id, "mode": args.mode,
                  "tree": tree_label, "error": str(e), "runtime_sec": round(elapsed, 2)}

    with open(os.path.join(sim_output, "result.json"), "w") as f:
        json.dump(result, f, indent=2)

    print(f"[sim {args.sim_id}] done in {elapsed:.1f}s → {sim_output}/result.json")


if __name__ == "__main__":
    main()
