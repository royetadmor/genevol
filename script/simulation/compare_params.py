#!/usr/bin/env python3
"""
Compare true simulation parameters to parameters recovered by genevol.
Usage: python3 compare_params.py simulation_param_file.txt < genevol_log.txt
"""
import sys
import re

PARAM_NAMES = ['gain', 'loss', 'innovation', 'elimination']

def parse_true_params(path):
    """Returns dict: param_name -> list of float values (comma-separated in file)."""
    true = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or '=' not in line:
                continue
            key, val = line.split('=', 1)
            key = key.strip().lstrip('_')
            if key in PARAM_NAMES:
                try:
                    true[key] = [float(v.strip()) for v in val.split(',')]
                except ValueError:
                    pass
    return true

def parse_recovered_params(log):
    """Returns dict: param_name -> {index -> float}, where index is the param index."""
    recovered = {}
    for line in log.splitlines():
        m = re.match(r'\s+GeneCount\.(\w+?)(\d+)_\d+\s*=\s*([\d.eE+\-]+)', line)
        if m:
            name, idx, val = m.group(1), int(m.group(2)), float(m.group(3))
            if name in PARAM_NAMES:
                recovered.setdefault(name, {})[idx] = val
    return recovered

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: compare_params.py <simulation_param_file>", file=sys.stderr)
        sys.exit(1)

    true = parse_true_params(sys.argv[1])
    recovered = parse_recovered_params(sys.stdin.read())

    print(f"{'Parameter':<18} {'True':>10} {'Recovered':>12} {'Abs diff':>12} {'Error%':>10}")
    print('-' * 66)
    for p in PARAM_NAMES:
        true_vals = true.get(p, [])
        rec_map   = recovered.get(p, {})
        n = max(len(true_vals), len(rec_map))
        for i in range(n):
            label = f"{p}[{i}]"
            t = true_vals[i] if i < len(true_vals) else None
            r = rec_map.get(i)
            if t is not None and r is not None:
                diff = abs(r - t)
                err  = diff / t * 100 if t != 0 else float('inf')
                print(f"{label:<18} {t:>10.4f} {r:>12.4f} {diff:>12.4f} {err:>9.1f}%")
            else:
                print(f"{label:<18} {t!s:>10} {r!s:>12}")
