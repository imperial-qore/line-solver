"""Emit the comparison report as a single self-contained HTML page."""
import html
import json
import math
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'nc_asymptotics_report.html')

METHODS = ['kt', 'le', 'ble', 'nrl', 'nrp', 'nre']
LI = ['kt', 'le', 'ble']
NR = ['nrl', 'nrp', 'nre']
COLOR = {'kt': 'c1', 'le': 'c2', 'ble': 'c3', 'nrl': 'c4', 'nrp': 'c5', 'nre': 'c6'}
FLOOR = 1e-8
CEIL = 10.0

rows = [json.loads(l) for l in open(os.path.join(HERE, 'grid_results.jsonl'))]
A = [r for r in rows if r['block'] == 'A']


def med(vals):
    v = [x for x in vals if x is not None and np.isfinite(x)]
    return float(np.median(v)) if v else None


def stat(block, method, field, key=None, val=None, stat='median'):
    sel = [r[field] for r in rows if r['block'] == block and r['method'] == method
           and (key is None or r[key] == val)]
    sel = [x for x in sel if x is not None and np.isfinite(x)]
    if not sel:
        return None
    if stat == 'median':
        return float(np.median(sel))
    if stat == 'p90':
        return float(np.percentile(sel, 90))
    return float(np.max(sel))


# ---------------------------------------------------------------- log scale
def domain(series):
    """Decade-aligned log domain covering the drawn values, floored at FLOOR."""
    vals = [max(v, FLOOR) for _, ys in series for v in ys if v is not None]
    lo = math.floor(math.log10(min(vals)))
    hi = math.ceil(math.log10(max(vals)))
    if hi - lo < 2:                       # never squeeze a panel below two decades
        hi = lo + 2
    return lo, hi


def panel(title, xs, xlabel, series, w=452, h=248, pad=(46, 74, 34, 12)):
    """One log-y line panel. series = [(method, [values aligned with xs])]"""
    pl, pr, pb, pt = pad
    x0, x1, y0, y1 = pl, w - pr, pt + 18, h - pb
    lo, hi = domain(series)

    def ylog(v):
        v = max(min(v if v and v > 0 else FLOOR, 10.0 ** hi), 10.0 ** lo)
        f = (math.log10(v) - lo) / (hi - lo)
        return y1 - f * (y1 - y0)

    out = [f'<svg class="panel" viewBox="0 0 {w} {h}" role="img" '
           f'aria-label="{html.escape(title)}">']
    out.append(f'<text class="ptitle" x="0" y="11">{html.escape(title)}</text>')
    step = 1 if hi - lo <= 5 else 2
    for e in range(lo, hi + 1, step):
        y = ylog(10.0 ** e)
        out.append(f'<line class="grid" x1="{x0}" y1="{y:.1f}" x2="{x1}" y2="{y:.1f}"/>')
        out.append(f'<text class="tick ty" x="{x0 - 7}" y="{y + 3.5:.1f}">1e{e}</text>')
    n = len(xs)

    def px(i):
        return x0 if n == 1 else x0 + i * (x1 - x0) / (n - 1)

    for i, xv in enumerate(xs):
        out.append(f'<text class="tick tx" x="{px(i):.1f}" y="{y1 + 16}">{xv}</text>')
    out.append(f'<line class="axis" x1="{x0}" y1="{y1}" x2="{x1}" y2="{y1}"/>')
    out.append(f'<text class="axlabel" x="{(x0 + x1) / 2:.0f}" y="{h - 6}">'
               f'{html.escape(xlabel)}</text>')

    labels = []
    for m, vals in series:
        pts, clamped = [], []
        for i, v in enumerate(vals):
            if v is None:
                continue
            pts.append((px(i), ylog(v)))
            clamped.append(v < 10.0 ** lo)
        if not pts:
            continue
        d = ' '.join(f'{"M" if k == 0 else "L"}{x:.1f} {y:.1f}'
                     for k, (x, y) in enumerate(pts))
        out.append(f'<path class="line s-{COLOR[m]}" d="{d}"/>')
        for (x, y), cl in zip(pts, clamped):
            cls = 'dot hollow' if cl else 'dot'
            out.append(f'<circle class="{cls} s-{COLOR[m]}" cx="{x:.1f}" cy="{y:.1f}" '
                       f'r="3.4"/>')
        labels.append([m, pts[-1][0], pts[-1][1], pts[-1][1]])

    # push overlapping end labels apart, keeping them inside the panel, and draw a
    # leader wherever a label had to move off its own line
    labels.sort(key=lambda t: t[3])
    for k in range(1, len(labels)):
        if labels[k][3] - labels[k - 1][3] < 12:
            labels[k][3] = labels[k - 1][3] + 12
    if labels and labels[-1][3] > h - pb - 2:
        shift = labels[-1][3] - (h - pb - 2)
        for lab in labels:
            lab[3] -= shift
    for m, lx, ly, ty in labels:
        if abs(ty - ly) > 1.5:
            out.append(f'<path class="leader s-{COLOR[m]}" d="M{lx + 4:.1f} {ly:.1f} '
                       f'L{lx + 9:.1f} {ty:.1f}"/>')
        out.append(f'<text class="endlab s-{COLOR[m]}" x="{lx + 11:.1f}" '
                   f'y="{ty + 3.5:.1f}">{m}</text>')
    out.append('</svg>')
    return '\n'.join(out)


def rankchart(field, title, w=520, h=222):
    """Ranked median / p90 / max range plot over Block A. Identity is the row label."""
    data = []
    for m in METHODS:
        data.append((m, stat('A', m, field), stat('A', m, field, stat='p90'),
                     stat('A', m, field, stat='max')))
    data.sort(key=lambda t: t[1])
    pl, pr, pt, pb = 54, 62, 24, 30
    x0, x1 = pl, w - pr
    rowh = (h - pt - pb) / len(data)
    out = [f'<svg class="panel" viewBox="0 0 {w} {h}" role="img" '
           f'aria-label="{html.escape(title)}">']
    out.append(f'<text class="ptitle" x="0" y="11">{html.escape(title)}</text>')

    def px(v):
        v = max(min(v, CEIL), FLOOR)
        f = (math.log10(v) - math.log10(FLOOR)) / (math.log10(CEIL) - math.log10(FLOOR))
        return x0 + f * (x1 - x0)

    for e in range(-8, 2, 2):
        x = px(10.0 ** e)
        out.append(f'<line class="grid" x1="{x:.1f}" y1="{pt - 6}" x2="{x:.1f}" '
                   f'y2="{h - pb}"/>')
        out.append(f'<text class="tick tx" x="{x:.1f}" y="{h - pb + 15}">1e{e}</text>')
    for i, (m, mn, p90, mx) in enumerate(data):
        y = pt + rowh * (i + 0.5)
        out.append(f'<line class="rangebar" x1="{px(mn):.1f}" y1="{y:.1f}" '
                   f'x2="{px(mx):.1f}" y2="{y:.1f}"/>')
        out.append(f'<line class="p90tick" x1="{px(p90):.1f}" y1="{y - 5:.1f}" '
                   f'x2="{px(p90):.1f}" y2="{y + 5:.1f}"/>')
        out.append(f'<circle class="dot s-{COLOR[m]}" cx="{px(mn):.1f}" cy="{y:.1f}" r="4.6"/>')
        out.append(f'<text class="rowlab" x="{x0 - 10}" y="{y + 4:.1f}">{m}</text>')
        out.append(f'<text class="rowval" x="{x1 + 8}" y="{y + 4:.1f}">{mn:.1e}</text>')
    out.append('</svg>')
    return '\n'.join(out)


# --------------------------------------------------------------- panel data
Ms, Rs, Ns = [2, 3, 5, 10, 20], [1, 2, 3], [2, 5, 10, 20, 50, 100]
Ns_t = [2, 5, 10, 20, 50]


def sers(fam, field, key, keys):
    return [(m, [stat('A', m, field, key, k) for k in keys]) for m in fam]


figs = {}
for fam, tag in ((LI, 'li'), (NR, 'nr')):
    label = 'Load-independent expansions' if tag == 'li' else 'Norlund-Rice family'
    figs[f'M_{tag}'] = panel(label, Ms, 'queueing stations  M',
                             sers(fam, 'abs_dlG', 'M', Ms))
    figs[f'R_{tag}'] = panel(label, Rs, 'job classes  R', sers(fam, 'abs_dlG', 'R', Rs))
    figs[f'N_{tag}'] = panel(label, Ns, 'total population  N',
                             sers(fam, 'abs_dlG', 'Ntot', Ns))
    figs[f'T_{tag}'] = panel(label, Ns_t, 'total population  N',
                             sers(fam, 'tput_err', 'Ntot', Ns_t))
    figs[f'leg_{tag}'] = fam
bench = json.load(open(os.path.join(HERE, 'bench_nrl_nre.json')))


def bsel(sweep, method, var):
    r = sorted([b for b in bench if b['sweep'] == sweep and b['method'] == method],
               key=lambda b: b[var])
    return [b[var] for b in r], [b['cpu'] for b in r], r


def benchpanel(sweep, var, xlabel, title):
    xs, ya, ra = bsel(sweep, 'nrl', var)
    _, yb, rb = bsel(sweep, 'nre', var)
    return panel(title, xs, xlabel, [('nrl', ya), ('nre', yb)])


def callpanel():
    """kernel calls against class count -- the whole of the cost difference"""
    xs, _, ra = bsel('R', 'nrl', 'R')
    _, _, rb = bsel('R', 'nre', 'R')
    ca = [float(r['gld'] or r['gldsingle']) for r in ra]
    cb = [float(r['gldsingle']) for r in rb]
    return panel('kernel evaluations per constant', xs, 'job classes  R',
                 [('nrl', ca), ('nre', cb)])


figs['b_N'] = benchpanel('N', 'Ntot', 'total population  N  (M=5, R=2)',
                         'CPU seconds per log G')
figs['b_M'] = benchpanel('M', 'M', 'queueing stations  M  (R=2, N=50)',
                         'CPU seconds per log G')
figs['b_R'] = benchpanel('R', 'R', 'job classes  R  (M=5, N=48)',
                         'CPU seconds per log G')
figs['b_C'] = callpanel()

brow = []
for sw, var, lab in (('N', 'Ntot', 'N'), ('M', 'M', 'M'), ('R', 'R', 'R')):
    xs, _, ra = bsel(sw, 'nrl', var)
    _, _, rb = bsel(sw, 'nre', var)
    for p_, q_ in zip(ra, rb):
        na = p_['gld'] or p_['gldsingle']
        nb = q_['gldsingle']
        brow.append(
            f'<tr><th scope="row">{lab} = {p_[var]}</th>'
            f'<td>{p_["cpu"]:.3f}</td><td>{q_["cpu"]:.3f}</td>'
            f'<td>{q_["cpu"] / p_["cpu"]:.2f}&times;</td>'
            f'<td>{na}</td><td>{nb}</td>'
            f'<td>{p_["kernel_s"] / na:.4f}</td><td>{q_["kernel_s"] / nb:.4f}</td></tr>')
brow = ''.join(brow)
r1 = [b for b in bench if b['sweep'] == 'R' and b['R'] == 1]
nrl_r1 = next(b for b in r1 if b['method'] == 'nrl')['cpu']
nre_r1 = next(b for b in r1 if b['method'] == 'nre')['cpu']

figs['rank_lg'] = rankchart('abs_dlG', 'median, p90 and worst |log G - log G|')
figs['rank_x'] = rankchart('tput_err', 'median, p90 and worst relative error in X')


# ------------------------------------------------------------------- tables
def table(field, keyname, keys, caption, stat_='median'):
    head = ''.join(f'<th>{m}</th>' for m in METHODS)
    body = []
    for k in keys:
        cells = []
        best = min((v for v in (stat('A', m, field, keyname, k) for m in METHODS)
                    if v is not None), default=None)
        for m in METHODS:
            v = stat('A', m, field, keyname, k, stat_)
            if v is None:
                cells.append('<td class="na">-</td>')
            else:
                cls = ' class="best"' if best is not None and v <= best * 1.0000001 else ''
                cells.append(f'<td{cls}>{v:.2e}</td>')
        body.append(f'<tr><th scope="row">{k}</th>{"".join(cells)}</tr>')
    return (f'<figure class="tblwrap"><table><caption>{html.escape(caption)}</caption>'
            f'<thead><tr><th scope="col">{keyname}</th>{head}</tr></thead>'
            f'<tbody>{"".join(body)}</tbody></table></figure>')


overall = []
for m in METHODS:
    overall.append((m,
                    stat('A', m, 'abs_dlG'), stat('A', m, 'abs_dlG', stat='p90'),
                    stat('A', m, 'abs_dlG', stat='max'),
                    stat('A', m, 'tput_err'), stat('A', m, 'tput_err', stat='p90'),
                    stat('A', m, 'tput_err', stat='max'),
                    stat('A', m, 'secs', 'Ntot', 50)))
overall.sort(key=lambda t: t[1])
orow = ''.join(
    f'<tr><th scope="row"><span class="swatch s-{COLOR[m]}"></span>{m}</th>'
    f'<td>{a:.2e}</td><td>{b:.2e}</td><td>{c:.2e}</td>'
    f'<td>{d:.2e}</td><td>{e:.2e}</td><td>{f:.2e}</td><td>{g:.2f}</td></tr>'
    for m, a, b, c, d, e, f, g in overall)

# structural identities, recomputed here so the page never quotes a stale number
by = defaultdict(dict)
for r in rows:
    by[(r['block'], r['M'], r['R'], r['Ntot'], r['rep'])][r['method']] = r
c_ble = 1 - math.log(2 * math.pi) / 2
lediff = {M: float(np.mean([by[k]['le']['lG'] - by[k]['ble']['lG']
                            for k in by if k[0] == 'A' and k[1] == M])) for M in Ms}
nrdiff1 = float(np.mean([by[k]['nrl']['lG'] - by[k]['nrp']['lG']
                         for k in by if k[0] == 'A' and k[2] == 1]))
logistic_lap = math.log(0.25 * math.sqrt(4 * math.pi))
xdiff = max(abs(by[k]['le']['tput_err'] - by[k]['ble']['tput_err']) for k in by
            if by[k]['le']['tput_err'] is not None
            and np.isfinite(by[k]['le']['tput_err']))
blow = [r for r in rows if r['method'] in NR and abs(r['dlG']) > 1]
nrtot = len([r for r in rows if r['method'] in NR])
nmodels = len({(r['block'], r['M'], r['R'], r['Ntot'], r['rep']) for r in rows})

lerow = ' '.join(f'<span class="kv"><b>M={M}</b>{lediff[M]:+.6f}</span>' for M in Ms)
blowrow = ''.join(
    f'<tr><td>{r["block"]}</td><td>{r["M"]}</td><td>{r["R"]}</td><td>{r["Ntot"]}</td>'
    f'<td><span class="swatch s-{COLOR[r["method"]]}"></span>{r["method"]}</td>'
    f'<td>{r["dlG"]:+.3f}</td></tr>'
    for r in sorted(blow, key=lambda r: -abs(r['dlG'])))

CSS = """
:root{
  color-scheme: light;
  --plane:#f2f3f5; --surface:#fbfbfc; --raise:#ffffff;
  --ink:#0d0e10; --ink2:#4b5058; --muted:#7d838c;
  --rule:#dcdfe4; --grid:#e6e8ec; --axis:#c2c6cd;
  --c1:#2a78d6; --c2:#eb6834; --c3:#1baf7a;
  --c4:#4a3aa7; --c5:#e34948; --c6:#008300;
  --good:#0ca30c; --crit:#d03b3b;
  --range:#c9ced6; --best:rgba(42,120,214,.10);
}
:root:not([data-theme="light"]){}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    color-scheme: dark;
    --plane:#101113; --surface:#17181a; --raise:#1e2023;
    --ink:#f0f1f3; --ink2:#b6bcc4; --muted:#868c95;
    --rule:#2b2e33; --grid:#26282c; --axis:#3a3e44;
    --c1:#3987e5; --c2:#d95926; --c3:#199e70;
    --c4:#9085e9; --c5:#e66767; --c6:#008300;
    --good:#0ca30c; --crit:#e66767;
    --range:#3c4149; --best:rgba(57,135,229,.16);
  }
}
:root[data-theme="dark"]{
  color-scheme: dark;
  --plane:#101113; --surface:#17181a; --raise:#1e2023;
  --ink:#f0f1f3; --ink2:#b6bcc4; --muted:#868c95;
  --rule:#2b2e33; --grid:#26282c; --axis:#3a3e44;
  --c1:#3987e5; --c2:#d95926; --c3:#199e70;
  --c4:#9085e9; --c5:#e66767; --c6:#008300;
  --good:#0ca30c; --crit:#e66767;
  --range:#3c4149; --best:rgba(57,135,229,.16);
}

*{box-sizing:border-box}
body{
  margin:0; background:var(--plane); color:var(--ink);
  font-family:"IBM Plex Sans", ui-sans-serif, system-ui, sans-serif;
  font-size:15px; line-height:1.62; -webkit-font-smoothing:antialiased;
}
.wrap{max-width:1140px; margin:0 auto; padding:44px 26px 84px;
      display:flex; flex-direction:column; gap:44px}
.prose{max-width:69ch}
h1,h2,h3{font-family:"IBM Plex Serif", Georgia, serif; text-wrap:balance;
         margin:0; font-weight:600; letter-spacing:-.012em}
h1{font-size:clamp(28px,4.1vw,42px); line-height:1.12}
h2{font-size:23px; line-height:1.25}
h3{font-size:16.5px; line-height:1.3; font-weight:600}
p{margin:0 0 .85em}
p:last-child{margin-bottom:0}
a{color:var(--c1)}
code,.mono{font-family:"IBM Plex Mono", ui-monospace, Menlo, monospace}
code{font-size:.885em; background:var(--raise); border:1px solid var(--rule);
     border-radius:3px; padding:.08em .34em}

.eyebrow{font-family:"IBM Plex Mono", monospace; font-size:11.5px;
  letter-spacing:.13em; text-transform:uppercase; color:var(--muted); margin:0 0 12px}
.lede{font-size:17.5px; line-height:1.55; color:var(--ink2); margin-top:14px}

header .meta{display:flex; flex-wrap:wrap; gap:0 26px; margin-top:22px;
  padding-top:16px; border-top:1px solid var(--rule);
  font-family:"IBM Plex Mono", monospace; font-size:12.5px; color:var(--muted)}
header .meta b{color:var(--ink2); font-weight:500}

section{display:flex; flex-direction:column; gap:18px}
.sechead{display:flex; align-items:baseline; gap:14px; border-bottom:1px solid var(--rule);
  padding-bottom:9px}
.sechead .n{font-family:"IBM Plex Mono", monospace; font-size:12px; color:var(--muted)}

.figgrid{display:grid; grid-template-columns:repeat(auto-fit,minmax(330px,1fr)); gap:20px}
.card{background:var(--surface); border:1px solid var(--rule); border-radius:5px;
      padding:16px 18px 12px}
.panel{width:100%; height:auto; display:block; overflow:visible}
.ptitle{font-family:"IBM Plex Sans",sans-serif; font-size:11.5px; font-weight:600;
        fill:var(--ink2); letter-spacing:.02em}
.grid{stroke:var(--grid); stroke-width:1}
.axis{stroke:var(--axis); stroke-width:1}
.tick{font-family:"IBM Plex Mono",monospace; font-size:9.5px; fill:var(--muted)}
.ty{text-anchor:end}
.tx{text-anchor:middle}
.axlabel{font-family:"IBM Plex Sans",sans-serif; font-size:10.5px; fill:var(--muted);
         text-anchor:middle}
.line{fill:none; stroke-width:2; stroke-linejoin:round; stroke-linecap:round}
.dot{stroke:var(--surface); stroke-width:1.6}
.dot.hollow{fill:var(--surface); stroke-width:2}
.leader{fill:none; stroke-width:1.2; opacity:.75}
.endlab{font-family:"IBM Plex Mono",monospace; font-size:11px; font-weight:500}
.rowlab{font-family:"IBM Plex Mono",monospace; font-size:12px; fill:var(--ink);
        text-anchor:end}
.rowval{font-family:"IBM Plex Mono",monospace; font-size:10.5px; fill:var(--muted)}
.rangebar{stroke:var(--range); stroke-width:5; stroke-linecap:round}
.p90tick{stroke:var(--ink2); stroke-width:1.6}
.s-c1{stroke:var(--c1)} .s-c2{stroke:var(--c2)} .s-c3{stroke:var(--c3)}
.s-c4{stroke:var(--c4)} .s-c5{stroke:var(--c5)} .s-c6{stroke:var(--c6)}
circle.s-c1{fill:var(--c1)} circle.s-c2{fill:var(--c2)} circle.s-c3{fill:var(--c3)}
circle.s-c4{fill:var(--c4)} circle.s-c5{fill:var(--c5)} circle.s-c6{fill:var(--c6)}
circle.hollow.s-c1,circle.hollow.s-c2,circle.hollow.s-c3,
circle.hollow.s-c4,circle.hollow.s-c5,circle.hollow.s-c6{fill:var(--surface)}
text.s-c1{fill:var(--c1); stroke:none} text.s-c2{fill:var(--c2); stroke:none}
text.s-c3{fill:var(--c3); stroke:none} text.s-c4{fill:var(--c4); stroke:none}
text.s-c5{fill:var(--c5); stroke:none} text.s-c6{fill:var(--c6); stroke:none}

.legend{display:flex; flex-wrap:wrap; gap:6px 18px; margin-top:4px;
  font-family:"IBM Plex Mono",monospace; font-size:11.5px; color:var(--ink2)}
.legend span{display:inline-flex; align-items:center; gap:7px}
.swatch{width:10px; height:10px; border-radius:2px; display:inline-block;
  flex:0 0 auto; margin-right:7px; vertical-align:-1px}
.swatch.s-c1{background:var(--c1)} .swatch.s-c2{background:var(--c2)}
.swatch.s-c3{background:var(--c3)} .swatch.s-c4{background:var(--c4)}
.swatch.s-c5{background:var(--c5)} .swatch.s-c6{background:var(--c6)}
.legend .swatch{margin-right:0}

.tblwrap{margin:0; overflow-x:auto; background:var(--surface);
  border:1px solid var(--rule); border-radius:5px}
table{border-collapse:collapse; width:100%; font-variant-numeric:tabular-nums;
  font-family:"IBM Plex Mono",monospace; font-size:12.5px}
caption{caption-side:top; text-align:left; padding:13px 16px 9px; color:var(--ink2);
  font-family:"IBM Plex Sans",sans-serif; font-size:12.5px; font-weight:500}
th,td{padding:6px 13px; text-align:right; white-space:nowrap}
thead th{color:var(--muted); font-weight:500; border-bottom:1px solid var(--rule);
  padding-bottom:8px}
tbody th{text-align:left; color:var(--ink); font-weight:500}
tbody tr+tr td, tbody tr+tr th{border-top:1px solid var(--grid)}
td.best{background:var(--best); color:var(--ink); font-weight:500}
td.na{color:var(--muted)}
.tblnote{font-size:12.5px; color:var(--muted); padding:9px 16px 12px;
  border-top:1px solid var(--grid); font-family:"IBM Plex Sans",sans-serif}

.findings{display:flex; flex-direction:column; gap:0;
  border:1px solid var(--rule); border-radius:5px; background:var(--surface);
  overflow:hidden}
.finding{display:grid; grid-template-columns:minmax(0,1fr) minmax(0,1.05fr);
  gap:18px 26px; padding:17px 20px; align-items:start}
.finding+.finding{border-top:1px solid var(--rule)}
.finding h3{margin-bottom:5px}
.finding p{font-size:14px; color:var(--ink2); margin:0}
.ident{font-family:"IBM Plex Mono",monospace; font-size:12.5px; color:var(--ink);
  background:var(--raise); border:1px solid var(--rule); border-radius:4px;
  padding:11px 13px; line-height:1.75; overflow-x:auto}
.ident b{font-weight:600}
.kv{display:inline-block; margin-right:14px; white-space:nowrap; color:var(--ink2)}
.kv b{color:var(--muted); font-weight:500; margin-right:5px}

.rec{display:grid; grid-template-columns:repeat(auto-fit,minmax(250px,1fr)); gap:16px}
.rec > div{background:var(--surface); border:1px solid var(--rule); border-radius:5px;
  padding:15px 17px; display:flex; flex-direction:column; gap:7px}
.rec .verdict{font-family:"IBM Plex Mono",monospace; font-size:14px; font-weight:600;
  display:flex; align-items:center; gap:8px}
.rec p{font-size:13.5px; color:var(--ink2); margin:0}
.pill{font-family:"IBM Plex Mono",monospace; font-size:10.5px; letter-spacing:.05em;
  text-transform:uppercase; padding:2px 7px; border-radius:99px; align-self:flex-start;
  border:1px solid currentColor}
.pill.use{color:var(--good)} .pill.avoid{color:var(--crit)} .pill.ok{color:var(--muted)}

footer{border-top:1px solid var(--rule); padding-top:18px; color:var(--muted);
  font-size:13px}
footer code{font-size:12px}
@media (max-width:720px){ .finding{grid-template-columns:1fr} }
@media (prefers-reduced-motion:no-preference){}
"""


def fig(inner, note=None, legend=None):
    leg = ''
    if legend:
        leg = ('<div class="legend">' + ''.join(
            f'<span><i class="swatch s-{COLOR[m]}"></i>{m}</span>' for m in legend)
            + '</div>')
    nt = f'<div class="tblnote" style="border:0;padding:6px 0 0">{note}</div>' if note else ''
    return f'<div class="card">{inner}{leg}{nt}</div>'


HTML = f"""<title>Six Asymptotic Constants</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Serif:wght@500;600&display=swap">
<style>{CSS}</style>

<div class="wrap">
<header>
  <p class="eyebrow">LINE / SolverNC &middot; normalizing constant</p>
  <h1>Six asymptotic normalizing constants on load-independent closed networks</h1>
  <p class="lede prose"><code>kt</code>, <code>le</code>, <code>ble</code>,
  <code>nrl</code>, <code>nrp</code> and <code>nre</code> measured against exact
  convolution over a grid of queueing stations, job classes and populations.
  <b>nre wins outright</b> &mdash; three to four decades below every other method
  on <code>log G</code> and two on the throughput ratio &mdash; and it is the only
  one of the six that never blows up. It is also the slowest, by about 20x.</p>
  <div class="meta">
    <span><b>models</b> {nmodels}</span>
    <span><b>evaluations</b> {len(rows)}</span>
    <span><b>reference</b> pfqn_ca (exact)</span>
    <span><b>codebase</b> python/line_solver</span>
    <span><b>failures</b> none</span>
  </div>
</header>

<section>
  <div class="sechead"><span class="n">01</span><h2>How the six rank</h2></div>
  <div class="prose"><p>Every model is a closed product-form network with
  single-server load-independent queues, random demands
  <code>L ~ U(0.1, 1)</code> and no think time. Error is measured two ways: on
  <code>log G(N)</code> itself, and on the throughput ratio
  <code>X&#8203;<sub>r</sub> = G(N &minus; e<sub>r</sub>) / G(N)</code>, which is
  what the solver actually reports and where an additive bias in
  <code>log G</code> cancels. Dot is the median, the vertical tick p90, the bar
  runs to the worst case over 425 models.</p></div>
  <div class="figgrid">
    {fig(figs['rank_lg'])}
    {fig(figs['rank_x'], note='Ratio error is measured on the 350 models with N &le; 50, where all R+1 constants were evaluated exactly.')}
  </div>
  <figure class="tblwrap">
    <table>
      <caption>Block A summary, sorted by median <code>|log G</code> error<code>|</code>.
      Seconds are for one model at N = 50 on an eight-way loaded host, indicative only.</caption>
      <thead><tr><th scope="col">method</th>
        <th colspan="3">|log G error|</th><th colspan="3">relative error in X</th>
        <th>s @ N=50</th></tr>
        <tr><th scope="col"></th><th>median</th><th>p90</th><th>worst</th>
        <th>median</th><th>p90</th><th>worst</th><th></th></tr></thead>
      <tbody>{orow}</tbody>
    </table>
  </figure>
</section>

<section>
  <div class="sechead"><span class="n">02</span><h2>What each axis does</h2></div>
  <div class="prose"><p>Median <code>|log G</code> error<code>|</code> on a log
  scale, each panel cut to the decades its own data occupies and floored at
  <code>1e-8</code>; a hollow marker sits on that floor and means the value is
  below it. The two families are drawn apart because they fail in different
  places, not because they are measured differently &mdash; the axes, the models
  and the reference are identical, so a value in the left panel is directly
  comparable with one in the right.</p></div>
  <div class="figgrid">{fig(figs['M_li'], legend=LI)}{fig(figs['M_nr'], legend=NR)}</div>
  <div class="figgrid">{fig(figs['R_li'], legend=LI)}{fig(figs['R_nr'], legend=NR)}</div>
  <div class="figgrid">{fig(figs['N_li'], legend=LI)}{fig(figs['N_nr'], legend=NR)}</div>
  <div class="prose"><p>Only <code>le</code> depends on the station count, and it
  depends on it linearly &mdash; that is the whole of its disadvantage against
  <code>ble</code> (&sect;03). Only the Norlund-Rice methods depend on the class
  count, and <code>nrl</code>/<code>nrp</code> lose about half a decade per class
  added. Population is the axis on which <code>kt</code> earns its keep: it is an
  asymptotic in <code>N</code> and is the worst of the six at
  <code>N = 2</code>.</p></div>
  <div class="figgrid">
    {fig(figs['T_li'], legend=LI, note='Throughput ratio, not the constant. le and ble coincide exactly.')}
    {fig(figs['T_nr'], legend=NR, note='nre is the only one of the three that improves monotonically with population.')}
  </div>
</section>

<section>
  <div class="sechead"><span class="n">03</span><h2>Four exact identities behind the numbers</h2></div>
  <div class="findings">

    <div class="finding">
      <div><h3>le and ble differ by a constant, and it is the published one</h3>
      <p>The gap between the two is <code>&minus;(M&minus;1)(1 &minus; log 2&pi; / 2)</code>
      to six decimals at every station count, and nothing else. So
      <code>le</code>'s error growing to 1.6 at <code>M = 20</code> is not the
      expansion degrading &mdash; it is the missing correction term, and
      <code>ble</code> is flat in <code>M</code> once it is restored.</p></div>
      <div class="ident">measured mean of <b>log G(le) &minus; log G(ble)</b><br>{lerow}
      <br><span class="kv"><b>&minus;(M&minus;1)&middot;{c_ble:.6f}</b>reproduces every row</span></div>
    </div>

    <div class="finding">
      <div><h3>That constant cancels in every mean value</h3>
      <p>Because the offset is additive in <code>log G</code> and independent of
      <code>N</code>, it divides out of <code>G(N &minus; e<sub>r</sub>)/G(N)</code>.
      Across all 350 ratio models the two methods do not merely agree closely,
      they agree to floating point. Choosing <code>le</code> over <code>ble</code>
      changes only <code>getProbNormConstAggr</code>.</p></div>
      <div class="ident">max over 350 models of<br>
      <b>|X&#8203;<sub>err</sub>(le) &minus; X&#8203;<sub>err</sub>(ble)|</b>
      &nbsp;=&nbsp; {xdiff:.1e}</div>
    </div>

    <div class="finding">
      <div><h3>nrl is nrp minus the Laplace error of the logistic density</h3>
      <p>Both methods push the Norlund-Rice contour through a transform of the
      unit torus &mdash; a logit for <code>nrl</code>, a probit for
      <code>nrp</code> &mdash; then Laplace-approximate. At <code>R = 1</code> the
      torus collapses (<code>t &minus; t&#772; &equiv; 0</code>) and the integral
      reduces to the transform's own density, which integrates to one. Laplace is
      <em>exact</em> on a Gaussian, so <code>nrp</code> is exact there; on the
      logistic it is not, and the shortfall is a pure constant.</p></div>
      <div class="ident">measured mean at R = 1 of<br>
      <b>log G(nrl) &minus; log G(nrp)</b> &nbsp;=&nbsp; {nrdiff1:+.6f}<br>
      <b>log(&frac14;&radic;(4&pi;))</b> &nbsp;=&nbsp; {logistic_lap:+.6f}
      <br><span class="kv"><b>consequence</b>nrp &lt; nrl at every R</span></div>
    </div>

    <div class="finding">
      <div><h3>The Laplace step occasionally fails outright</h3>
      <p><code>nrl</code> and <code>nrp</code> share a numerical Hessian, and on
      {len(blow)//2} of the {nmodels} models it returns a curvature that puts
      <code>log G</code> up to eight nats <em>above</em> the truth &mdash; both
      methods, same models, always at <code>R &ge; 2</code>. That is not a large
      error, it is a wrong answer, and nothing in the returned value flags it.
      <code>nre</code>, which tilts the contour to the saddle instead, never
      exceeds 0.19 anywhere in the grid.</p></div>
      <div class="ident"><table style="font-size:11.5px">
        <thead><tr><th>blk</th><th>M</th><th>R</th><th>N</th><th>method</th>
        <th>&Delta; log G</th></tr></thead>
        <tbody>{blowrow}</tbody></table>
        <div style="color:var(--muted);padding-top:7px">{len(blow)} of {nrtot}
        Norlund-Rice evaluations exceed one nat</div></div>
    </div>

  </div>
</section>

<section>
  <div class="sechead"><span class="n">04</span><h2>Full grids</h2></div>
  {table('abs_dlG', 'M', Ms, 'Median |log G error| by station count, Block A. Best in each row shaded.')}
  {table('abs_dlG', 'R', Rs, 'Median |log G error| by class count. At R = 1 both nrp (4.1e-08) and nre (2.2e-16) are exact up to their own arithmetic.')}
  {table('abs_dlG', 'Ntot', Ns, 'Median |log G error| by total population.')}
  {table('tput_err', 'Ntot', Ns_t, 'Median relative error in the throughput ratio X_r, by total population. N = 100 is omitted: the R+1 reference constants were not evaluated there.')}
  {table('secs', 'Ntot', Ns, 'Median WALL seconds per model from the accuracy grid: eight workers on a loaded host, and N <= 50 sums R+1 evaluations where N = 100 is one. Kept for the ordering only -- the clean per-evaluation CPU cost is in section 05.')}
</section>

<section>
  <div class="sechead"><span class="n">05</span><h2>What nre costs, and why</h2></div>
  <div class="prose"><p>Measured separately from the accuracy grid, because the
  host it ran on carries a load average near 80 and wall clock there says nothing.
  These are <b>CPU seconds for one <code>log G</code> evaluation</b>, best of
  three, single-threaded BLAS. Both methods spend 97&ndash;99% of that inside the
  load-dependent kernel <code>pfqn_gld</code> / <code>pfqn_gldsingle</code>, so
  the comparison reduces to two questions: what one kernel call costs, and how
  many each method makes.</p></div>
  <div class="figgrid">{fig(figs['b_N'], legend=['nrl', 'nre'])}{fig(figs['b_M'], legend=['nrl', 'nre'])}</div>
  <div class="figgrid">{fig(figs['b_R'], legend=['nrl', 'nre'])}{fig(figs['b_C'], legend=['nrl', 'nre'], note='Counts, not seconds. nrl is 4R&sup2;-R+1 exactly; nre is adaptive.')}</div>
  <div class="prose"><p><b>Per call they cost the same</b> &mdash; within 20%,
  scaling as <code>M<sup>1.05</sup>N<sup>1.9</sup></code> for both, flat in
  <code>R</code>. So the entire runtime difference is the call count, and the call
  counts have completely different characters. <code>nrl</code> makes exactly
  <code>4R&sup2; &minus; R + 1</code> of them &mdash; a numerical Hessian over
  <code>R</code> dimensions, three evaluations per diagonal entry and four per
  off-diagonal, plus one at the mode &mdash; and that count does not move with
  <code>M</code> or <code>N</code>. <code>nre</code>'s is set by a Newton root
  find and is adaptive: 1, 36, 114, 244 at <code>R = 1&hellip;4</code> here, and
  27 to 75 across the <code>R = 2</code> runs alone.</p>
  <p>Two things follow. <b>At one class <code>nre</code> is the faster method</b>,
  by {nrl_r1 / nre_r1:.1f}&times; ({nre_r1:.3f}s against {nrl_r1:.3f}s): the torus collapses to
  dimension zero, so it returns <code>pfqn_gldsingle</code> after a single call
  while <code>nrl</code> still builds a Hessian. It is also exact there and
  <code>nrl</code> is 0.12 off, so that case has no trade-off in it at all. From
  two classes up the ratio is roughly 2&times; and climbs with the class count
  &mdash; 1.8&times;, 2.6&times;, 3.2&times; at <code>R = 2, 3, 4</code>. The
  spread you see along the <code>M</code> and <code>N</code> sweeps (1.7&times; to
  3.6&times;) is not a size effect: it is <code>nre</code>'s iteration count
  moving, and it disappears once the cost is divided by calls.</p></div>
  <figure class="tblwrap">
    <table>
      <caption>CPU seconds per constant, kernel calls made, and seconds per call.
      Sweeps vary one axis at a time about M = 5, R = 2, N = 50.</caption>
      <thead><tr><th scope="col"></th>
        <th colspan="3">CPU s per log G</th>
        <th colspan="2">kernel calls</th><th colspan="2">s per call</th></tr>
        <tr><th scope="col"></th><th>nrl</th><th>nre</th><th>ratio</th>
        <th>nrl</th><th>nre</th><th>nrl</th><th>nre</th></tr></thead>
      <tbody>{brow}</tbody>
    </table>
    <div class="tblnote">Read the ratio column down each block, not across blocks:
    the three sweeps overlap only at M = 5, R = 2.</div>
  </figure>
</section>

<section>
  <div class="sechead"><span class="n">06</span><h2>Which to use</h2></div>
  <div class="rec">
    <div><span class="pill use">default</span>
      <div class="verdict"><span class="swatch s-c6"></span>nre</div>
      <p>Best on both metrics at every station count, class count and population
      in the grid, exact at <code>R = 1</code>, and the only method with no
      catastrophic case. Pay for it in time: ~20x <code>ble</code> at
      <code>N = 50</code> and it grows with <code>N</code>, so it is a poor fit
      inside a fixed-point loop.</p></div>
    <div><span class="pill use">when speed matters</span>
      <div class="verdict"><span class="swatch s-c3"></span>ble</div>
      <p>Ratio error 1e-3 median, 1.9e-2 worst, flat in <code>M</code>, and about
      as cheap as the six get. The natural choice when the constant feeds mean
      values rather than being reported itself.</p></div>
    <div><span class="pill ok">large N only</span>
      <div class="verdict"><span class="swatch s-c1"></span>kt</div>
      <p>The only method whose ratio error falls monotonically with population,
      from 7.8e-2 at <code>N = 2</code> to 7.5e-4 at <code>N = 50</code>. Below
      <code>N &asymp; 10</code> it is the worst of the six on ratios; above it,
      competitive.</p></div>
    <div><span class="pill ok">reporting log G only</span>
      <div class="verdict"><span class="swatch s-c2"></span>le</div>
      <p>Identical to <code>ble</code> for every mean value. Prefer it only when
      the published Cas17 eq. 34 form is what is wanted; on <code>log G</code>
      itself it is worse by <code>(M&minus;1)&middot;0.081</code>.</p></div>
    <div><span class="pill avoid">not on these models</span>
      <div class="verdict"><span class="swatch s-c5"></span>nrp</div>
      <p>Exact at <code>R = 1</code>, then 3.3e-1 on <code>log G</code> and 9e-2
      on ratios at <code>R = 3</code>, with the shared Hessian failure. Its reason
      to exist is genuine load dependence, which this grid does not exercise.</p></div>
    <div><span class="pill avoid">not on these models</span>
      <div class="verdict"><span class="swatch s-c4"></span>nrl</div>
      <p>Strictly <code>nrp</code> shifted by {abs(nrdiff1):.3f}, so it is behind
      <code>nrp</code> everywhere, including at <code>R = 1</code> where
      <code>nrp</code> is exact and it is not. Its one advantage is a fixed
      <code>4R&sup2;&minus;R+1</code> kernel budget (&sect;05).</p></div>
  </div>
</section>

<footer class="prose">
  <p><b>Method.</b> <code>kt</code>, <code>le</code> and <code>ble</code> are
  called through <code>pfqn_nc</code>. <code>nrl</code>, <code>nrp</code> and
  <code>nre</code> are load-dependent evaluators, so they are driven the way
  <code>SolverNC</code> drives them on a load-independent model &mdash;
  <code>pfqn_ncld</code> with <code>mu &equiv; 1</code> &mdash; which was checked
  to reproduce <code>SolverNC(model, 'nrl'|'nrp'|'nre')</code> exactly. The
  reference <code>pfqn_ca</code> was checked against brute-force enumeration of
  the product-form state space to machine precision.</p>
  <p>Two further blocks were run and are in <code>summary.txt</code>: a think-time
  block (<code>Z = 5</code> per class) and a near-balanced-demand block
  (<code>L = 1 &plusmn; 1%</code>). Neither changes the ranking. The balanced
  block sharpens it: <code>ble</code>'s ratio error falls to 2.5e-7 at
  <code>M = 20</code> while <code>kt</code> stays at 8.7e-4, and
  <code>nrl</code>/<code>nrp</code> reach 3.0e-2 at <code>R = 3</code>. The
  think-time block produces the worst Norlund-Rice failure in the study, +8.0
  nats at <code>M = 2, R = 3, N = 50</code>.</p>
  <p>Grid, raw results and summary:
  <code>experiments/nc_asymptotics/</code>.</p>
</footer>
</div>
"""

with open(OUT, 'w') as fh:
    fh.write(HTML)
print('wrote', OUT, len(HTML), 'bytes')
