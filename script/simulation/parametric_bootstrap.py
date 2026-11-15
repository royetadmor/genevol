#!/usr/bin/env python3
"""
Usage:
  # Sequential (local / Docker):
  python3 parametric_bootstrap.py --param-file param_file.txt

  # HPC (Run in parallel on a cluster with Slurm):
  python3 parametric_bootstrap.py --param-file param_file.txt --hpc [--n-replicates 100]
"""

import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time

SCRIPT_DIR   = os.path.dirname(os.path.realpath(__file__))
REPO_ROOT    = os.path.dirname(os.path.dirname(SCRIPT_DIR))
GENEVOL_BIN  = os.path.join(REPO_ROOT, 'GenEvol', 'genEvol')
SIM_BIN      = os.path.join(REPO_ROOT, 'Simulator', 'simulator')
SLURM_SCRIPT = os.path.join(SCRIPT_DIR, 'bootstrap_rep.slurm')


# ── Param file I/O ─────────────────────────────────────────────────────────────

def parse_param_file(path):
    """Return list of (key, value) pairs; blank/comment lines stored as (None, raw)."""
    entries = []
    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                entries.append((None, stripped))
                continue
            m = re.match(r'^(\S+)\s*=\s*(.+)$', stripped)
            if m:
                entries.append((m.group(1), m.group(2).strip()))
            else:
                entries.append((None, stripped))
    return entries


def entries_to_dict(entries):
    return {k: v for k, v in entries if k is not None}


def set_param(entries, key, value):
    """Update existing key in-place, or append if absent."""
    for i, (k, v) in enumerate(entries):
        if k == key:
            entries[i] = (key, str(value))
            return
    entries.append((key, str(value)))


def remove_params(entries, keys):
    keys = set(keys)
    return [(k, v) for k, v in entries if k not in keys]


def write_param_file(path, entries):
    with open(path, 'w') as f:
        for k, v in entries:
            if k is None:
                f.write(v + '\n')
            else:
                f.write(f"{k} = {v}\n")


def resolve_path(raw, base_dir):
    """Resolve a param-file path relative to base_dir if it isn't absolute."""
    if not raw or os.path.isabs(raw):
        return raw
    return os.path.realpath(os.path.join(base_dir, raw))


# ── Run helpers ────────────────────────────────────────────────────────────────

def run_binary(binary, param_path, extra_args=None):
    """Run a genevol-family binary; return combined stdout+stderr. Raises on failure."""
    cmd    = [binary, f"param={param_path}"] + (extra_args or [])
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] {os.path.basename(binary)} failed (exit {result.returncode})", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"{os.path.basename(binary)} exited with code {result.returncode}")
    return result.stdout + result.stderr


# ── Output parsers ─────────────────────────────────────────────────────────────

def parse_ml_params(stdout):
    """Extract the 'Parameter values:' block from genEvol output.
    Returns {bpp_param_name: float}."""
    params   = {}
    in_block = False
    for line in stdout.split('\n'):
        if 'Parameter values:' in line:
            in_block = True
            continue
        if in_block:
            m = re.match(r'\s+([\w.]+)\s*=\s*([\d.eE+\-]+)', line)
            if m:
                params[m.group(1)] = float(m.group(2))
            elif line.strip() == '':
                break
    return params


def build_sim_param_substitutions(ml_params, pdict):
    """Build all param-file substitutions needed to seed the simulator with ML values.
    Handles multi-subparam functions (e.g. LINEAR: lossK for K=0,1 → '_loss=v0,v1').
    Also updates rate_distribution when Gamma ASRV was used."""
    subst = {}
    for ptype, key in [('loss',        '_loss'),
                       ('gain',        '_gain'),
                       ('innovation',  '_innovation'),
                       ('elimination', '_elimination')]:
        sub = {}
        for name, val in ml_params.items():
            m = re.match(rf'GeneCount\.{ptype}(\d+)_\d+', name)
            if m:
                sub[int(m.group(1))] = val
        if sub:
            subst[key] = ','.join(str(sub[i]) for i in sorted(sub))
    if 'Poisson.lambda' in ml_params:
        subst['_rootLambda'] = str(ml_params['Poisson.lambda'])
    if 'Gamma.alpha' in ml_params:
        rd = pdict.get('rate_distribution', '')
        m  = re.search(r'n\s*=\s*(\d+)', rd)
        if m:
            subst['rate_distribution'] = f'Gamma(n={m.group(1)}, alpha={ml_params["Gamma.alpha"]})'
    return subst


def parse_bootstrap_delta_aic(stdout):
    """Extract ΔAIC from a 1-WGD-capped genEvol detect run.
    Returns 0.0 if no WGD was detected."""
    if 'No WGD events detected' in stdout:
        return 0.0
    # Detection table row: "  1   <childId>   <q>   <deltaAIC>"
    m = re.search(r'^\s*1\s+\d+\s+[\d.eE+\-]+\s+([\d.eE+\-]+)', stdout, re.MULTILINE)
    return float(m.group(1)) if m else 0.0


# ── Dataset helpers ────────────────────────────────────────────────────────────

def count_gene_families(data_path):
    """Count gene family columns in the TSV data file (header minus species column)."""
    with open(data_path) as f:
        header = f.readline()
    return max(1, len(header.strip().split('\t')) - 1)


# ── ASCII histogram ────────────────────────────────────────────────────────────

def ascii_histogram(values, title, n_bins=12, bar_width=32):
    if not values:
        return
    lo, hi = min(values), max(values)
    if lo == hi:
        hi = lo + 1.0
    bw     = (hi - lo) / n_bins
    counts = [0] * n_bins
    for v in values:
        counts[min(int((v - lo) / bw), n_bins - 1)] += 1
    peak = max(counts) or 1
    print(f"\n{title}")
    for i in range(n_bins):
        left  = lo + i * bw
        right = left + bw
        bar   = '█' * int(counts[i] / peak * bar_width)
        print(f"  {left:7.2f} – {right:7.2f}  |{bar:<{bar_width}}  {counts[i]}")


# ── Progress / Results ─────────────────────────────────────────────────────────

def report_progress(i, n_reps, last_delta_aic, delta_aics, t_start):
    elapsed = time.time() - t_start
    eta     = elapsed / i * (n_reps - i)
    cur_pct = sorted(delta_aics)[int(len(delta_aics) * 0.95) - 1]
    print(f'  [{i:>{len(str(n_reps))}}/{n_reps}]  '
          f'ΔAIC={last_delta_aic:.2f}  '
          f'running 95th pct={cur_pct:.2f}  '
          f'ETA={eta:.0f}s')



def print_results(delta_aics, out_dir):
    n_reps  = len(delta_aics)
    pct     = 0.95
    sorted_ = sorted(delta_aics)
    threshold = sorted_[max(0, int(n_reps * pct) - 1)]
    mean_val  = sum(delta_aics) / n_reps
    std_val   = (sum((x - mean_val) ** 2 for x in delta_aics) / n_reps) ** 0.5
    n_zero    = sum(1 for v in delta_aics if v == 0.0)

    ascii_histogram(delta_aics, f'Bootstrap ΔAIC distribution ({n_reps} replicates):')

    dist_path = os.path.join(out_dir, 'bootstrap_delta_aic.txt')
    with open(dist_path, 'w') as f:
        f.write('delta_aic\n')
        for v in delta_aics:
            f.write(f'{v}\n')
    print(f'\nFull distribution written to: {dist_path}')
    print(f'(Replicates with no WGD detected: {n_zero}/{n_reps})')

    w   = 58
    sep = '═' * w
    print(f'\n╔{sep}╗')
    print(f'║  Bootstrap ΔAIC threshold (95th pct, {n_reps} reps): {threshold:<8.4f}            ║')
    print(f'║  Distribution: mean={mean_val:.2f}  std={std_val:.2f}  max={max(delta_aics):.2f}{" " * (w - 38 - len(f"{mean_val:.2f}{std_val:.2f}{max(delta_aics):.2f}"))}║')
    print(f'║{" " * w}║')
    print(f'║  Add these lines to your param file and re-run genevol:   ║')
    print(f'║    _wgdMode      = detect                                 ║')
    print(f'║    _wgdThreshold = {threshold:<8.4f}                              ║')
    print(f'╚{sep}╝')


# ── Bootstrap phases ───────────────────────────────────────────────────────────

def _run_setup(args):
    """Steps 1+2: fit base model, write sim_template.txt and inf_template.txt to out_dir.
    Returns (sim_base, orig_entries, tree_path) for use by the calling loop."""
    param_file = os.path.realpath(args.param_file)
    param_dir  = os.path.dirname(param_file)
    out_dir    = os.path.realpath(args.output_dir)

    print('=' * 60)
    print('Step 1: Fitting base model (no WGD) on real data...')
    print('=' * 60)

    orig_entries = parse_param_file(param_file)
    pdict        = entries_to_dict(orig_entries)

    tree_path = resolve_path(pdict.get('_treePath', ''), param_dir)
    data_path = resolve_path(pdict.get('_dataPath', ''), param_dir)

    base_entries = list(orig_entries)
    set_param(base_entries, '_wgdMode', 'disable')
    set_param(base_entries, '_treePath', tree_path)
    if data_path:
        set_param(base_entries, '_dataPath', data_path)
    os.makedirs(out_dir, exist_ok=True)
    base_param_path = os.path.join(out_dir, 'base_model_params.txt')
    write_param_file(base_param_path, base_entries)

    print(f"  Running: {GENEVOL_BIN}")
    base_stdout = run_binary(GENEVOL_BIN, base_param_path)
    with open(os.path.join(out_dir, 'base_model_log.txt'), 'w') as f:
        f.write(base_stdout)

    ml_params = parse_ml_params(base_stdout)
    if not ml_params:
        print('[ERROR] Could not parse ML parameters from genEvol output.', file=sys.stderr)
        print('        Check base_model_log.txt for details.', file=sys.stderr)
        sys.exit(1)

    print('  ML parameters extracted:')
    for name, val in sorted(ml_params.items()):
        print(f'    {name} = {val}')

    sim_subst = build_sim_param_substitutions(ml_params, pdict)

    print('\nStep 2: Preparing simulation templates...')

    num_sites = count_gene_families(data_path) if data_path else 1000
    print(f'  _numSites = {num_sites}')

    # sim_template: ML rates applied; _seed/_outputFasta/_outputTree added per-rep
    sim_base = list(orig_entries)
    for key, val in sim_subst.items():
        set_param(sim_base, key, val)
    set_param(sim_base, '_treePath', tree_path)
    set_param(sim_base, '_dataPath', data_path)
    set_param(sim_base, '_wgdMode',  'disable')
    set_param(sim_base, '_numSites', str(num_sites))
    sim_base = remove_params(sim_base, ['_fixedWgdQ', '_wgdThreshold', '_maxWgdEvents',
                                        '_seed', '_outputFasta', '_outputTree'])
    write_param_file(os.path.join(out_dir, 'sim_template.txt'), sim_base)

    # inf_template: detection settings applied; _dataPath added per-rep
    inf_base = list(orig_entries)
    set_param(inf_base, '_treePath',     tree_path)
    set_param(inf_base, '_wgdMode',      'detect')
    set_param(inf_base, '_maxWgdEvents', '1')
    set_param(inf_base, '_wgdThreshold', '0')
    inf_base = remove_params(inf_base, ['_fixedWgdQ', '_dataPath'])
    write_param_file(os.path.join(out_dir, 'inf_template.txt'), inf_base)

    with open(os.path.join(out_dir, 'bootstrap_config.json'), 'w') as f:
        json.dump({'seed': args.seed, 'n_reps': args.n_replicates}, f)

    print(f'  Templates written to {out_dir}/')
    return sim_base, orig_entries, tree_path


def _run_rep(args, rep_id=None):
    """Run one replicate. Reads templates from out_dir, writes delta_aic.txt, returns the value."""
    out_dir = os.path.realpath(args.output_dir)
    if rep_id is None:
        rep_id = args.rep_id

    with open(os.path.join(out_dir, 'bootstrap_config.json')) as f:
        config = json.load(f)

    rep_dir  = os.path.join(out_dir, f'rep_{rep_id:04d}')
    os.makedirs(rep_dir, exist_ok=True)

    sim_template = parse_param_file(os.path.join(out_dir, 'sim_template.txt'))
    inf_template = parse_param_file(os.path.join(out_dir, 'inf_template.txt'))

    sim_fasta = os.path.join(rep_dir, 'simulated.fasta')
    sim_nwk   = os.path.join(rep_dir, 'simulated.nwk')

    rep_sim = list(sim_template)
    set_param(rep_sim, '_seed',        str(config['seed'] + rep_id))
    set_param(rep_sim, '_outputFasta', sim_fasta)
    set_param(rep_sim, '_outputTree',  sim_nwk)
    sim_param_path = os.path.join(rep_dir, 'sim_params.txt')
    write_param_file(sim_param_path, rep_sim)
    run_binary(SIM_BIN, sim_param_path)

    rep_inf = list(inf_template)
    set_param(rep_inf, '_dataPath', sim_fasta)
    inf_param_path = os.path.join(rep_dir, 'inf_params.txt')
    write_param_file(inf_param_path, rep_inf)

    inf_stdout = run_binary(GENEVOL_BIN, inf_param_path)
    with open(os.path.join(rep_dir, 'genevol_log.txt'), 'w') as f:
        f.write(inf_stdout)

    delta_aic = parse_bootstrap_delta_aic(inf_stdout)
    with open(os.path.join(rep_dir, 'delta_aic.txt'), 'w') as f:
        f.write(f'{delta_aic}\n')

    return delta_aic


def _run_collect(out_dir):
    """Glob rep_*/delta_aic.txt, read all values, and call print_results."""
    files = sorted(glob.glob(os.path.join(out_dir, 'rep_*/delta_aic.txt')))
    if not files:
        print('[ERROR] No replicate results found in ' + out_dir, file=sys.stderr)
        sys.exit(1)
    delta_aics = []
    for path in files:
        with open(path) as f:
            delta_aics.append(float(f.read().strip()))
    print_results(delta_aics, out_dir)


def _run_full(args):
    """Sequential bootstrap: setup then loop over all replicates."""
    _run_setup(args)

    out_dir = os.path.realpath(args.output_dir)
    n_reps  = args.n_replicates

    print(f'\nStep 3: Running {n_reps} bootstrap replicates...')
    print('=' * 60)

    delta_aics = []
    t_loop = time.time()

    for i in range(1, n_reps + 1):
        delta_aic = _run_rep(args, rep_id=i)
        delta_aics.append(delta_aic)
        if i % 10 == 0 or i == n_reps:
            report_progress(i, n_reps, delta_aic, delta_aics, t_loop)

    print_results(delta_aics, out_dir)

    if not args.keep_tmp:
        shutil.rmtree(out_dir)
        print('\n(Temporary files removed. Use --keep-tmp to retain them.)')


def _run_hpc(args):
    """HPC bootstrap: setup, submit SLURM array, poll until done, collect."""
    _run_setup(args)

    out_dir    = os.path.realpath(args.output_dir)
    param_file = os.path.realpath(args.param_file)
    n_reps     = args.n_replicates

    logs_dir = os.path.join(out_dir, 'logs')
    os.makedirs(logs_dir, exist_ok=True)

    print(f'\nSubmitting SLURM array job ({n_reps} tasks)...')
    result = subprocess.run(
        [
            'sbatch', '--parsable',
            f'--array=1-{n_reps}',
            f'--output={os.path.join(logs_dir, "bootstrap_%A_%a.out")}',
            f'--error={os.path.join(logs_dir, "bootstrap_%A_%a.err")}',
            SLURM_SCRIPT, out_dir, param_file, os.path.join(SCRIPT_DIR, 'parametric_bootstrap.py'),
        ],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print('[ERROR] sbatch failed:', result.stderr.strip(), file=sys.stderr)
        sys.exit(1)

    job_id = result.stdout.strip().split(';')[0]  # --parsable may append cluster name
    print(f'  Submitted job array {job_id}')
    print('  Waiting for replicates to finish (checking every 30 s)...')

    while True:
        sq = subprocess.run(['squeue', '-j', job_id, '-h'], capture_output=True, text=True)
        running = [l for l in sq.stdout.strip().splitlines() if l.strip()]
        if not running:
            break
        print(f'  {len(running)} task(s) still running...', end='\r', flush=True)
        time.sleep(30)

    print('\n  All tasks complete. Collecting results...')
    _run_collect(out_dir)

    if not args.keep_tmp:
        shutil.rmtree(out_dir)
        print('\n(Temporary files removed. Use --keep-tmp to retain them.)')


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Parametric bootstrap to calibrate WGD detection ΔAIC threshold.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--param-file',       required=True,
                        help='genEvol param file for the real dataset')
    parser.add_argument('--n-replicates', type=int, default=100, metavar='N',
                        help='Number of bootstrap replicates')
    parser.add_argument('--output-dir',   default='bootstrap_tmp', metavar='DIR',
                        help='Working directory for bootstrap files')
    parser.add_argument('--seed',         type=int, default=42,
                        help='Base seed; replicate i uses seed+i')
    parser.add_argument('--keep-tmp',     action='store_true',
                        help='Retain output-dir after completion')
    parser.add_argument('--hpc',          action='store_true',
                        help='Run replicates as a SLURM array job')
    # Internal flags: called back by bootstrap_rep.slurm on each compute node, not by the user
    parser.add_argument('--run-rep', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--rep-id',  type=int,            help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.run_rep:
        if args.rep_id is None:
            print('[ERROR] --run-rep requires --rep-id', file=sys.stderr)
            sys.exit(1)
        _run_rep(args)
        return 0

    if args.hpc:
        _run_hpc(args)
    else:
        _run_full(args)
    return 0


if __name__ == '__main__':
    sys.exit(main())
