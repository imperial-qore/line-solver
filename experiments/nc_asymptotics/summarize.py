"""Aggregate grid_results.jsonl into the comparison tables."""
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
METHODS = ['kt', 'le', 'ble', 'nrl', 'nrp', 'nre']


def load():
    rows = []
    with open(os.path.join(HERE, 'grid_results.jsonl')) as fh:
        for line in fh:
            rows.append(json.loads(line))
    return rows


def agg(vals):
    v = np.array([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
    n_bad = len([x for x in vals if x is None or not np.isfinite(x)])
    if v.size == 0:
        return None, None, n_bad
    return float(np.median(v)), float(np.max(v)), n_bad


def table(rows, keyfn, keylabel, field, title, fmt='{:9.2e}', stat='median'):
    buckets = defaultdict(lambda: defaultdict(list))
    for r in rows:
        buckets[keyfn(r)][r['method']].append(r[field])
    keys = sorted(buckets)
    out = [title, '-' * len(title)]
    out.append(f"{keylabel:>10s} | " + " ".join(f"{m:>9s}" for m in METHODS))
    for k in keys:
        cells = []
        for m in METHODS:
            med, mx, bad = agg(buckets[k][m])
            val = med if stat == 'median' else mx
            cells.append('     n/a ' if val is None else fmt.format(val))
        out.append(f"{str(k):>10s} | " + " ".join(cells))
    return "\n".join(out)


def main():
    rows = load()
    blocks = {
        'A': 'Block A - pure queueing network (Z=0), random demands L~U(0.1,1)',
        'B': 'Block B - with think time (Z=5 per class), random demands',
        'C': 'Block C - near-balanced demands (L=1+-1%), Z=0',
    }
    lines = []
    lines.append(f"models: {len({(r['block'], r['M'], r['R'], r['Ntot'], r['rep']) for r in rows})}"
                 f"   method evaluations: {len(rows)}")
    fails = defaultdict(int)
    for r in rows:
        if r['error']:
            fails[(r['method'], r['error'])] += 1
        elif not np.isfinite(r['lG']):
            fails[(r['method'], 'non-finite lG')] += 1
    lines.append("failures: " + (", ".join(f"{m}/{e}={n}" for (m, e), n in sorted(fails.items()))
                                 or "none"))
    for b, desc in blocks.items():
        sub = [r for r in rows if r['block'] == b]
        if not sub:
            continue
        lines += ['', '=' * 78, desc, '=' * 78]
        for field, stat, what in [('abs_dlG', 'median', 'median |log G_hat - log G|'),
                                  ('abs_dlG', 'max', 'worst-case |log G_hat - log G|'),
                                  ('tput_err', 'median', 'median max_r rel. err. of X_r = G(N-e_r)/G(N)'),
                                  ('tput_err', 'max', 'worst-case rel. err. of X_r')]:
            lines.append('')
            lines.append(table(sub, lambda r: r['M'], 'stations M', field,
                               f'{what}, by station count', stat=stat))
            lines.append('')
            lines.append(table(sub, lambda r: r['R'], 'classes R', field,
                               f'{what}, by class count', stat=stat))
            lines.append('')
            lines.append(table(sub, lambda r: r['Ntot'], 'total N', field,
                               f'{what}, by total population', stat=stat))
        lines.append('')
        lines.append(table(sub, lambda r: r['Ntot'], 'total N', 'secs',
                           'median seconds per model (all G evaluations, 8-way loaded host)',
                           fmt='{:9.3f}'))
    # overall ranking on block A
    lines += ['', '=' * 78, 'Overall (Block A)', '=' * 78]
    for field, what in [('abs_dlG', '|log G| error'), ('tput_err', 'throughput rel. error')]:
        lines.append('')
        lines.append(f'{what}: median / p90 / max over all Block A models')
        for m in METHODS:
            v = np.array([r[field] for r in rows
                          if r['block'] == 'A' and r['method'] == m
                          and r[field] is not None and np.isfinite(r[field])])
            if v.size:
                lines.append(f"  {m:>5s}  n={v.size:4d}  {np.median(v):.3e}  "
                             f"{np.percentile(v, 90):.3e}  {np.max(v):.3e}")
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(HERE, 'summary.txt'), 'w') as fh:
        fh.write(text + "\n")


if __name__ == '__main__':
    main()
