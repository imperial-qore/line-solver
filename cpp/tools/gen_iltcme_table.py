#!/usr/bin/env python3
"""Regenerate cpp/src/api/mam/iltcme_table.cpp from the vendored iltcme.json.

Emits only the entries a C++ caller can actually select. There are two callers
and they select differently:

  matlab_ilt takes "the entry of smallest cv2 whose n+1 <= maxFnEvals", and its
  only caller derives maxFnEvals = min(1000, max(11, iter_max)), so sweeping
  11..1000 enumerates its whole reachable set. Entry 0 is always reachable
  because matlab_ilt seeds the search with it before testing any bound.

  cme_table_entry (api/mam/cme.h, added 2026-08-01) takes the entry of smallest
  cv2 among those with a GIVEN n, for any n at all -- dist_fit_me walks the
  orders from the smallest upwards and stops at the first whose reach covers the
  target SCV, so the low-n tail matlab_ilt can never select is exactly what it
  uses most. Restricting the table to the matlab_ilt set left it unable to fit
  below 22 phases what MATLAB fits in 4.

Run from the repository root:  python3 cpp/tools/gen_iltcme_table.py
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, 'matlab/lib/thirdparty/iltcme/iltcme.json')
DST = os.path.join(ROOT, 'cpp/src/api/mam/iltcme_table.cpp')


def reachable_ilt(table):
    sel = set()
    for mfe in range(11, 1001):
        bi = 0
        for i, e in enumerate(table[1:], 1):
            if e['cv2'] < table[bi]['cv2'] and e['n'] + 1 <= mfe:
                bi = i
        sel.add(bi)
    return sel


def reachable_cme(table):
    best = {}
    for i, e in enumerate(table):
        n = e['n']
        if n not in best or e['cv2'] < table[best[n]]['cv2']:
            best[n] = i
    return set(best.values())


def reachable(table):
    return sorted(reachable_ilt(table) | reachable_cme(table))


def main():
    d = json.load(open(SRC))
    idx = reachable(d)
    f = repr
    out = ['''/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * GENERATED FILE -- do not edit by hand.
 * Regenerate with cpp/tools/gen_iltcme_table.py from
 * matlab/lib/thirdparty/iltcme/iltcme.json.
 *
 * VENDORED THIRD-PARTY DATA. These are the concentrated-matrix-exponential
 * (CME) coefficients of the ILT-CME project, the same table MATLAB's
 * matlab_ilt reads from iltcme.json. Their licence terms are NOT STATED
 * upstream; see THIRD-PARTY-NOTICES.md, which records the decision to vendor
 * with that known.
 *
 * ONLY THE SELECTABLE ENTRIES ARE VENDORED, as the union of what the two
 * consumers can reach: matlab_ilt takes the smallest cv2 with n+1 <= maxFnEvals
 * over the legal maxFnEvals range 11..1000, and cme_table_entry takes the
 * smallest cv2 at each distinct n. The rest can never be chosen by any input.
 */
#include "line/api/mam/iltcme_table.h"

namespace line {
namespace mam {
namespace iltcme {

namespace {
''']
    for k, i in enumerate(idx):
        e = d[i]
        out.append('const double A%d[] = {%s};\n' % (k, ','.join(f(float(v)) for v in e['a'])))
        out.append('const double B%d[] = {%s};\n' % (k, ','.join(f(float(v)) for v in e['b'])))
    out.append('\n}  // namespace\n\nconst CmeEntry kTable[] = {\n')
    for k, i in enumerate(idx):
        e = d[i]
        out.append('    {%d, %s, %s, %s, %s, A%d, B%d, %d},\n' % (
            e['n'], f(float(e['c'])), f(float(e['omega'])), f(float(e['mu1'])),
            f(float(e['cv2'])), k, k, len(e['a'])))
    out.append('};\n\nconst std::size_t kTableSize = sizeof(kTable) / sizeof(kTable[0]);\n\n')
    out.append('}  // namespace iltcme\n}  // namespace mam\n}  // namespace line\n')
    open(DST, 'w').write(''.join(out))
    print('wrote %s: %d of %d entries' % (DST, len(idx), len(d)))


if __name__ == '__main__':
    main()
