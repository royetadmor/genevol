#!/usr/bin/env python3
"""
Unified simulation study.

Usage:
  python3 run_study.py --mode standard   # rate parameter recovery across multiple trees
  python3 run_study.py --mode one_wgd    # WGD detection — single WGD event
  python3 run_study.py --mode two_wgds   # WGD detection — two WGD events

Adding a new mode: add an entry to MODES below — no other code changes needed.

Tree notes:
  sim_trees   — list of trees used for simulation, cycled across sims (round-robin)
  infer_tree  — tree used for genevol inference; None means same tree as simulation

Root lambda notes:
  root_lambda > 0  — fixed lambda used for both simulation and inference (no optimization)
  root_lambda <= 0 — simulation uses sim_root_lambda; inference optimizes from empirical mean
  sim_root_lambda  — required when root_lambda <= 0; the true lambda used to generate data

WGD notes:
  true_qs — list of q values per WGD node (DFS left-to-right order in sim_trees[0])
            [] disables WGD; non-empty enables detect mode and adds q columns to CSV
"""
import argparse
import csv
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
        "infer_tree":      None,        # same as sim tree
        "true_qs":         [],
        "root_lambda":     -1,          # -1 = optimize from empirical mean during inference
        "sim_root_lambda": 1.5,         # true lambda used to generate simulated data
        "rate_init":       "generic",   # genevol starts rates from generic values
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
            "test_data/tiley2016/Brassicaceae/tree_wgd.newick",
        ],
        "infer_tree":      "test_data/tiley2016/Brassicaceae/tree.newick",
        "true_qs":         [0.15],
        "root_lambda":     1.0,         # fixed lambda for both simulation and inference
        "rate_init":       "true",      # genevol starts rates from true values
        "wgd_threshold":   10,
        "target_events":   5.0,
        "rates": {
            "gain":        1.3,
            "loss":        1.0,
            "innovation":  0.1,
            "elimination": 0.05,
        },
        "results_csv":     "wgd_detection_results_one.csv",
        "sim_outputs_dir": "sim_outputs_one_wgd",
    },
    "two_wgds": {
        "description":     "WGD detection — two WGD events",
        "sim_trees": [
            "test_data/tiley2016/Brassicaceae/tree_wgds.newick",
        ],
        "infer_tree":      "test_data/tiley2016/Brassicaceae/tree.newick",
        "true_qs":         [0.4, 0.2],
        "root_lambda":     1.0,         # fixed lambda for both simulation and inference
        "rate_init":       "true",      # genevol starts rates from true values
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
}

# ── Shared constants ───────────────────────────────────────────────────────────

NUM_SIMS  = 5
NUM_SITES = 1000
MAX_STATE = 50

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
REPO_ROOT  = os.path.dirname(os.path.dirname(SCRIPT_DIR))

# ── CSV fields ─────────────────────────────────────────────────────────────────

def build_csv_fields(cfg):
    n      = len(cfg["true_qs"])
    multi  = len(cfg["sim_trees"]) > 1
    opt_lambda  = cfg["root_lambda"] <= 0

    fields = ["sim_id", "seed"]
    if multi:
        fields.append("tree")

    if n > 0:
        fields += [f"true_q_{i+1}" for i in range(n)]

    fields += ["true_gain", "true_loss", "true_innovation", "true_elimination"]
    if opt_lambda:
        fields.append("true_lambda")

    if n > 0:
        fields += ["n_detected"]
        fields += [f"detected_q_{i+1}" for i in range(n)]
        fields += [f"delta_aic_{i+1}"  for i in range(n)]
        fields += ["inferred_gain", "inferred_loss", "inferred_innovation", "inferred_elimination"]
    else:
        fields += ["recovered_gain", "recovered_loss", "recovered_innovation", "recovered_elimination"]
        if opt_lambda:
            fields.append("recovered_lambda")
        fields += ["diff_gain", "diff_loss", "diff_innovation", "diff_elimination"]
        fields += ["err_pct_gain", "err_pct_loss", "err_pct_innovation", "err_pct_elimination"]

    fields.append("runtime_sec")
    return fields

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
        # else: no _rootLambda → genevol infers lambda from empirical mean and optimizes it
        if n > 0:
            f.write(f"_wgdMode         = detect\n")
            f.write(f"_wgdThreshold    = {cfg['wgd_threshold']}\n")
            f.write(f"_modelCriterion  = AIC\n")
        else:
            f.write(f"_wgdMode         = disable\n")

# ── Docker runners ─────────────────────────────────────────────────────────────

def run_simulator(param_file_abs, sim_output_abs):
    param_filename = os.path.basename(param_file_abs)
    subprocess.run([
        "docker", "run", "--rm",
        "-v", f"{param_file_abs}:/app/genevol/{param_filename}:ro",
        "-v", f"{sim_output_abs}:/app/genevol/sim_output",
        "genevol",
        "./Simulator/simulator", f"param=/app/genevol/{param_filename}",
    ], check=True)


def run_genevol(param_file_abs, sim_output_abs):
    param_filename = os.path.basename(param_file_abs)
    result = subprocess.run([
        "docker", "run", "--rm",
        "-v", f"{param_file_abs}:/app/genevol/{param_filename}:ro",
        "-v", f"{sim_output_abs}:/app/genevol/sim_output:ro",
        "genevol",
        "./GenEvol/genEvol", f"param=/app/genevol/{param_filename}",
    ], capture_output=True, text=True)
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
    """Parse the === WGD Detection Results === summary block.

    Returns list of {q, delta_aic} in detection order, or [] if none detected.
    """
    if "No WGD events detected." in log:
        return []
    rows = []
    for m in re.finditer(
        r"^\s+\d+\s+\d+\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)", log, re.MULTILINE
    ):
        rows.append({"q": float(m.group(1)), "delta_aic": float(m.group(2))})
    return rows

# ── Row builders ───────────────────────────────────────────────────────────────

def fmt(v):
    return f"{v:.6f}" if v is not None else ""


def build_row(cfg, csv_fields, sim_id, tree_label, seed, rec, wgd_rows, elapsed):
    rates     = cfg["rates"]
    n         = len(cfg["true_qs"])
    n_det     = len(wgd_rows)
    opt_lambda     = cfg["root_lambda"] <= 0
    true_lambda    = cfg.get("sim_root_lambda") if opt_lambda else cfg["root_lambda"]

    row = {"sim_id": sim_id, "seed": seed}
    if "tree" in csv_fields:
        row["tree"] = tree_label

    for i, q in enumerate(cfg["true_qs"]):
        row[f"true_q_{i+1}"] = fmt(q)

    row.update({
        "true_gain":        fmt(rates["gain"]),
        "true_loss":        fmt(rates["loss"]),
        "true_innovation":  fmt(rates["innovation"]),
        "true_elimination": fmt(rates["elimination"]),
    })
    if opt_lambda:
        row["true_lambda"] = fmt(true_lambda)

    if n > 0:
        row["n_detected"] = n_det
        for i in range(n):
            if i < n_det:
                row[f"detected_q_{i+1}"] = f"{wgd_rows[i]['q']:.6f}"
                row[f"delta_aic_{i+1}"]  = f"{wgd_rows[i]['delta_aic']:.4f}"
            else:
                row[f"detected_q_{i+1}"] = ""
                row[f"delta_aic_{i+1}"]  = ""
        row.update({
            "inferred_gain":        fmt(rec.get("gain")),
            "inferred_loss":        fmt(rec.get("loss")),
            "inferred_innovation":  fmt(rec.get("innovation")),
            "inferred_elimination": fmt(rec.get("elimination")),
        })
    else:
        row.update({
            "recovered_gain":        fmt(rec.get("gain")),
            "recovered_loss":        fmt(rec.get("loss")),
            "recovered_innovation":  fmt(rec.get("innovation")),
            "recovered_elimination": fmt(rec.get("elimination")),
        })
        if opt_lambda:
            row["recovered_lambda"] = fmt(rec.get("lambda"))
        for p in ("gain", "loss", "innovation", "elimination"):
            t, r = rates[p], rec.get(p)
            if r is not None:
                d = abs(r - t)
                row[f"diff_{p}"]    = fmt(d)
                row[f"err_pct_{p}"] = fmt(d / t * 100) if t != 0 else ""
            else:
                row[f"diff_{p}"] = row[f"err_pct_{p}"] = ""

    row["runtime_sec"] = f"{elapsed:.1f}"
    return row


def build_error_row(cfg, csv_fields, sim_id, tree_label, seed, elapsed):
    rates = cfg["rates"]
    n     = len(cfg["true_qs"])
    opt_lambda  = cfg["root_lambda"] <= 0
    true_lambda = cfg.get("sim_root_lambda") if opt_lambda else cfg["root_lambda"]

    row = {"sim_id": sim_id, "seed": seed}
    if "tree" in csv_fields:
        row["tree"] = tree_label
    for i, q in enumerate(cfg["true_qs"]):
        row[f"true_q_{i+1}"] = fmt(q)
    row.update({
        "true_gain":        fmt(rates["gain"]),
        "true_loss":        fmt(rates["loss"]),
        "true_innovation":  fmt(rates["innovation"]),
        "true_elimination": fmt(rates["elimination"]),
    })
    if opt_lambda:
        row["true_lambda"] = fmt(true_lambda)
    if n > 0:
        row["n_detected"] = ""
        for i in range(n):
            row[f"detected_q_{i+1}"] = ""
            row[f"delta_aic_{i+1}"]  = ""
        row.update({
            "inferred_gain": "", "inferred_loss": "",
            "inferred_innovation": "", "inferred_elimination": "",
        })
    else:
        row.update({
            "recovered_gain": "", "recovered_loss": "",
            "recovered_innovation": "", "recovered_elimination": "",
        })
        if opt_lambda:
            row["recovered_lambda"] = ""
        for p in ("gain", "loss", "innovation", "elimination"):
            row[f"diff_{p}"] = row[f"err_pct_{p}"] = ""
    row["runtime_sec"] = f"{elapsed:.1f}"
    return row

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Unified simulation study")
    parser.add_argument("--mode", required=True, choices=list(MODES.keys()),
                        help="Study mode: " + ", ".join(MODES.keys()))
    parser.add_argument("--name", default="",
                        help="Optional label appended to output filenames, e.g. 'opt_lambda'")
    args = parser.parse_args()

    cfg        = MODES[args.mode]
    rates      = cfg["rates"]
    n_wgds     = len(cfg["true_qs"])
    target_len = cfg["target_events"] / max(rates.values())

    suffix      = f"_{args.name}" if args.name else ""
    base_csv    = cfg["results_csv"].replace(".csv", f"{suffix}.csv")
    results_csv = os.path.join(SCRIPT_DIR, base_csv)
    sim_outputs = os.path.join(SCRIPT_DIR, cfg["sim_outputs_dir"] + suffix)
    csv_fields  = build_csv_fields(cfg)

    # Precompute branch multipliers for each unique sim tree.
    branch_muls = {}
    for tree_path in cfg["sim_trees"]:
        raw = compute_tree_length(os.path.join(REPO_ROOT, tree_path))
        branch_muls[tree_path] = (os.path.basename(os.path.dirname(tree_path)), raw, target_len / raw)

    print(f"Mode:  {args.mode} — {cfg['description']}")
    print(f"Rates: {rates}")
    print(f"Trees:")
    for tp, (label, raw, bm) in branch_muls.items():
        print(f"  {label}: raw_length={raw:.1f}  branchMul={bm:.6f}  scaled={raw * bm:.2f}")
    if n_wgds:
        print(f"True q values:  {cfg['true_qs']}  ({n_wgds} WGD(s))")
        print(f"Inference tree: {cfg['infer_tree']}")
    print(f"Running {NUM_SIMS} simulations. Results -> {results_csv}\n")

    os.makedirs(sim_outputs, exist_ok=True)

    with open(results_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_fields)
        writer.writeheader()

        for sim_id in range(1, NUM_SIMS + 1):
            sim_tree              = cfg["sim_trees"][(sim_id - 1) % len(cfg["sim_trees"])]
            infer_tree            = cfg["infer_tree"] or sim_tree
            tree_label, _, bm     = branch_muls[sim_tree]
            seed                  = sim_id

            print(f"[{sim_id:02d}/{NUM_SIMS}] tree={tree_label}  seed={seed}  branchMul={bm:.6f}")

            sim_output = os.path.join(sim_outputs, f"sim_{sim_id:02d}")
            os.makedirs(sim_output, exist_ok=True)

            t_start = time.time()
            try:
                with tempfile.TemporaryDirectory() as tmpdir:
                    sim_param     = os.path.join(tmpdir, "sim_param.txt")
                    genevol_param = os.path.join(tmpdir, "genevol_param.txt")

                    write_sim_param_file(sim_param, cfg, sim_tree, bm, seed)
                    write_genevol_param_file(genevol_param, cfg, infer_tree, bm)

                    run_simulator(os.path.abspath(sim_param), sim_output)
                    log = run_genevol(os.path.abspath(genevol_param), sim_output)

                with open(os.path.join(sim_output, "genevol_log.txt"), "w") as lf:
                    lf.write(log)

                rec      = parse_recovered(log)
                wgd_rows = parse_wgd_results(log) if n_wgds > 0 else []
                elapsed  = time.time() - t_start

                row = build_row(cfg, csv_fields, sim_id, tree_label, seed, rec, wgd_rows, elapsed)
                writer.writerow(row)
                csvfile.flush()

                if n_wgds > 0:
                    n_det = len(wgd_rows)
                    q_str = "  ".join(
                        f"q{i+1}={wgd_rows[i]['q']:.4f}(ΔAIC={wgd_rows[i]['delta_aic']:.2f})"
                        for i in range(n_det)
                    )
                    print(f"         {n_det}/{n_wgds} detected  {q_str}")
                print(f"         gain={rec.get('gain','?'):.4f}  loss={rec.get('loss','?'):.4f}  "
                      f"innov={rec.get('innovation','?'):.4f}  elim={rec.get('elimination','?'):.4f}")

            except Exception as e:
                elapsed = time.time() - t_start
                print(f"         ERROR: {e}", file=sys.stderr)
                row = build_error_row(cfg, csv_fields, sim_id, tree_label, seed, elapsed)
                writer.writerow(row)
                csvfile.flush()

    print(f"\nDone. Results written to {results_csv}")


if __name__ == "__main__":
    main()
