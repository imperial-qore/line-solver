"""
Decision-diagram reachable set and Kronecker rate descriptor of a stochastic
Petri net, so that mdd_mcd can analyse it.

Levels are of two kinds:

* place levels: one per (Place, class) pair, holding a token count. A multiclass
  net therefore has P*R of them, ordered place-major, so level p*R+k is class k
  in place p.
* phase levels: one per mode whose firing time has more than one phase, holding
  the phase the running server occupies.

The rate structure factorises exactly under single-server firing semantics: a
mode fires at a constant rate whenever every input level holds its enabling
multiplicity and no inhibitor level has reached its threshold, so

    W_l^e[i, i + fire(l) - enab(l)] = 1     for enab(l) <= i < inhib(l)

at every place level. A phase-type mode contributes two event families on its
phase level, the internal phase changes D0 (marking unchanged) and the firings
D1 (marking moved), each gated by the same per-level enabling indicators. Both
are products of per-level terms, which is what Eq. 1 of the paper requires.

Phase-type firing and the memory policy
---------------------------------------
LINE discards a running server's phase when its mode becomes disabled, i.e.
preemptive repeat. Resetting a mode's phase is then triggered by a JOINT
condition on the place levels, which is not a product of per-level terms and has
no Kronecker form. What this descriptor encodes is preemptive resume: a disabled
mode's phase freezes and continues when the mode is re-enabled. The two policies
coincide exactly when a mode is never disabled while running, so reachability
records, for free, whether any phase-type mode was ever found disabled.
``phmemory='exact'`` (the default) errors when one was; ``'resume'`` proceeds
deliberately with the resume semantics.

Other restrictions, each an error and never a silent approximation: no immediate
transitions (they make vanishing states, which must be eliminated before a
Kronecker rate descriptor exists) and no marking-dependent firing rates. A
multi-server mode is accepted only when its enabling touches ONE level, because
the enabling degree min_l floor(m(l)/enab(l)) is otherwise not a product of
per-level terms.

See also: mdd_mcd, mdd_descriptor, mdd_reachset.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp

from ..io.logging import line_error, line_printf
from ..sn.network_struct import NodeType
from ..mdd import MDD


def _param(nparam, name, default=None):
    """Read a nodeparam field, which may be a dataclass or a dict."""
    if nparam is None:
        return default
    if isinstance(nparam, dict):
        return nparam.get(name, default)
    return getattr(nparam, name, default)


def _timing_name(t) -> str:
    if t is None:
        return ''
    if hasattr(t, 'name'):
        return str(t.name).upper()
    return str(t).upper().split('.')[-1]


def _arcvec(mat, places, R, nnodes, fillval):
    """(nnodes x nclasses) arc matrix -> place-major length-(P*R) level vector."""
    P = len(places)
    v = np.full(P * R, fillval, dtype=float)
    if mat is None:
        return v
    m = np.asarray(mat, dtype=float)
    if m.ndim == 1:
        m = m.reshape(nnodes, R)
    for pp in range(P):
        for k in range(R):
            x = m[places[pp], k]
            if np.isfinite(fillval) and fillval == 0:
                x = max(0.0, x)
            v[pp * R + k] = x
    return v


def spn_mdd(model, options: Optional[Dict] = None) -> Tuple[object, Dict, Dict]:
    """Build the reachable set and Kronecker descriptor of a stochastic Petri net.

    Parameters
    ----------
    model : a Network holding Places and Transitions
    options : dict with optional keys
        bound    - per-place-level token bound (scalar or length P*R); inferred
                   from a place invariant when omitted
        phmemory - 'exact' (default) or 'resume', see the module docstring
        descriptor - build the Kronecker rate descriptor (default True). Pass
                   False for the MDD-rec route, which reads only the reachable
                   set: the restrictions that exist purely because a Kronecker
                   form must factorise per level (marking-dependent firing
                   rates, multi-server modes drawing from several places) are
                   then lifted, in exchange for the firing times having to be
                   exponential
        verbose  - print the net summary (default False)

    Returns
    -------
    (mdds, desc, info) with mdds the MDD.to_struct export, desc the descriptor
    for mdd_mcd (carrying levelkind and, when one exists, the place invariant)
    and info the descriptive metadata.
    """
    if options is None:
        options = {}
    verbose = bool(options.get('verbose', False))
    bound_opt = options.get('bound')
    phmemory = str(options.get('phmemory', 'exact'))
    descriptor = bool(options.get('descriptor', True))

    sn = model.getStruct()
    R = int(sn.nclasses)
    nodetype = np.ravel(np.asarray(sn.nodetype, dtype=int))
    places = [int(i) for i in np.nonzero(nodetype == int(NodeType.PLACE))[0]]
    transitions = [int(i) for i in np.nonzero(nodetype == int(NodeType.TRANSITION))[0]]
    if not places or not transitions:
        line_error('spn_mdd', 'the model holds no Place or no Transition node')
    P = len(places)
    L = P * R                                  # place levels, place-major
    nnodes = int(sn.nnodes)

    # ---- collect the (transition, mode) pairs
    md: List[Dict] = []
    for ind in transitions:
        nparam = sn.nodeparam[ind] if sn.nodeparam is not None else None
        nmodes = int(_param(nparam, 'nmodes', 0) or 0)
        timing = _param(nparam, 'timingstrategies', None) or []
        firingdep = _param(nparam, 'firingdep', None) or []
        firingproc = _param(nparam, 'firingproc', None) or []
        firingpie = _param(nparam, 'firingpie', None) or []
        enabling = _param(nparam, 'enabling', None) or []
        inhibiting = _param(nparam, 'inhibiting', None) or []
        firing = _param(nparam, 'firing', None) or []
        nmodeservers = np.asarray(_param(nparam, 'nmodeservers', np.ones(nmodes)), dtype=float)
        for m in range(nmodes):
            if m < len(timing) and _timing_name(timing[m]) == 'IMMEDIATE':
                line_error('spn_mdd',
                           'mode %d of node %d is IMMEDIATE; vanishing states must be eliminated '
                           'before the net has a Kronecker rate descriptor' % (m + 1, ind + 1))
            dep = firingdep[m] if m < len(firingdep) else None
            if dep is not None and descriptor:
                line_error('spn_mdd',
                           'mode %d of node %d has a marking-dependent firing rate; g(marking) is '
                           'not a product of per-level terms' % (m + 1, ind + 1))
            proc = firingproc[m] if m < len(firingproc) else None
            if proc is None or len(proc) < 2 or proc[0] is None:
                line_error('spn_mdd',
                           'mode %d of node %d has no Markovian firing process; a general '
                           'distribution has no finite phase level' % (m + 1, ind + 1))
            D0 = np.atleast_2d(np.asarray(proc[0], dtype=float))
            D1 = np.atleast_2d(np.asarray(proc[1], dtype=float))
            nph = D0.shape[0]
            pv = firingpie[m] if m < len(firingpie) else None
            if pv is None or np.size(pv) == 0:
                pv = np.zeros(nph)
                pv[0] = 1.0
            pv = np.asarray(pv, dtype=float).ravel()
            md.append({
                'trans': ind,
                'mode': m,
                'enab': _arcvec(enabling[m] if m < len(enabling) else None, places, R, nnodes, 0.0),
                'inhib': _arcvec(inhibiting[m] if m < len(inhibiting) else None, places, R,
                                 nnodes, np.inf),
                'fire': _arcvec(firing[m] if m < len(firing) else None, places, R, nnodes, 0.0),
                'D0': D0,
                'D1': D1,
                'pie': pv / pv.sum(),
                'nph': nph,
                'srv': float(nmodeservers[m]) if m < nmodeservers.size else 1.0,
                'dep': dep,
            })
            if not descriptor and nph > 1:
                line_error('spn_mdd',
                           'mode %d of node %d has a phase-type firing time; the '
                           'reachable-set-only mode carries no phase level, and a product-form '
                           'marking process must be memoryless in the marking alone'
                           % (m + 1, ind + 1))
    E = len(md)
    for e in range(E) if descriptor else []:
        nz = int(np.count_nonzero(md[e]['enab']))
        if md[e]['srv'] != 1 and nz > 1:
            line_error('spn_mdd',
                       'mode %d of node %d has %g servers and draws from %d levels; the enabling '
                       'degree min_l floor(m(l)/enab(l)) is then not a product of per-level terms '
                       'and admits no Kronecker form'
                       % (md[e]['mode'] + 1, md[e]['trans'] + 1, md[e]['srv'], nz))
        if md[e]['srv'] != 1 and nz == 0:
            line_error('spn_mdd',
                       'mode %d of node %d has %g servers but consumes from no place, so its '
                       'enabling degree is unbounded and its firing rate undefined'
                       % (md[e]['mode'] + 1, md[e]['trans'] + 1, md[e]['srv']))

    # ---- phase levels for the multi-phase modes
    phaseof = np.zeros(E, dtype=int)           # 0 = no phase level (1-based when set)
    for e in range(E) if descriptor else []:
        if md[e]['nph'] > 1:
            phaseof[e] = L + int(np.count_nonzero(phaseof)) + 1
    Q = int(np.count_nonzero(phaseof))
    K = L + Q

    netm = np.zeros((E, L))
    for e in range(E):
        netm[e, :] = md[e]['fire'] - md[e]['enab']

    # ---- initial marking and per-level bounds
    init0 = _init_marking(model, sn, places, R)
    winv, vinv = _place_invariant(netm, init0)
    if bound_opt is not None and np.size(bound_opt) > 0:
        bound = np.asarray(bound_opt, dtype=float).ravel()
        if bound.size == 1:
            bound = np.full(L, float(bound[0]))
    elif winv is not None:
        bound = np.zeros(L)
        for l in range(L):
            bound[l] = np.floor(vinv / winv[l]) if winv[l] > 0 else vinv
    else:
        line_error('spn_mdd',
                   'the net has no place invariant with positive weights, so the marking is not '
                   'bounded a priori; pass options["bound"]')

    domain = np.zeros(K, dtype=int)
    domain[:L] = (bound + 1).astype(int)
    for e in range(E):
        if phaseof[e] > 0:
            domain[phaseof[e] - 1] = md[e]['nph']

    init = np.zeros(K, dtype=int)
    init[:L] = init0.astype(int)
    for e in range(E):
        if phaseof[e] > 0:
            init[phaseof[e] - 1] = int(np.nonzero(md[e]['pie'] > 0)[0][0])

    # ---- reachable set; the closure records which modes were ever disabled, so
    # the phase-memory question is answered without a second pass over |S|
    ever_disabled = np.zeros(E, dtype=bool)
    mdd = MDD(domain)
    init_t = tuple(int(v) for v in init)
    mdd.insert(init_t)
    frontier = [init_t]
    head = 0
    while head < len(frontier):
        s = frontier[head]
        head += 1
        succ, dis = _next(s, md, netm, phaseof, domain, L)
        ever_disabled |= dis
        for t in succ:
            if not mdd.member(t):
                mdd.insert(t)
                frontier.append(t)
        if head > 1024 and 2 * head > len(frontier):
            frontier = frontier[head:]
            head = 0
    mdd.compact()
    mdds = mdd.to_struct()

    badph = [e for e in range(E) if phaseof[e] > 0 and ever_disabled[e]]
    if badph and descriptor and phmemory.lower() != 'resume':
        line_error('spn_mdd',
                   'mode %d of node %d has a phase-type firing time AND is disabled in some '
                   'reachable marking. LINE discards the phase on disabling (preemptive repeat) '
                   'but that reset is a joint condition on the place levels and has no Kronecker '
                   'form, so this descriptor would encode preemptive resume instead and disagree '
                   'with SolverCTMC. Pass phmemory="resume" to accept the resume semantics.'
                   % (md[badph[0]]['mode'] + 1, md[badph[0]]['trans'] + 1))

    # ---- Kronecker event matrices
    events = []
    for e in range(E) if descriptor else []:
        gate = [l for l in range(L)
                if md[e]['enab'][l] > 0 or np.isfinite(md[e]['inhib'][l])]
        move = [l for l in range(L) if netm[e, l] != 0]
        touched = sorted(set(gate) | set(move))
        degl = -1
        if md[e]['srv'] != 1:
            nzl = np.nonzero(md[e]['enab'] > 0)[0]
            if nzl.size:
                degl = int(nzl[0])                       # level carrying the degree
        if phaseof[e] == 0:
            # single-phase mode: one firing event, rate on the first touched level
            if not touched:
                continue
            lev = []
            W = []
            for l in touched:
                lev.append(l)
                W.append(_placemat(l, md[e], netm[e, l], int(domain[l]), l == degl))
            W[0] = float(md[e]['D1'][0, 0]) * W[0]       # scalar rate
            events.append({'a': md[e]['trans'], 'b': md[e]['mode'], 'lev': lev, 'W': W})
        else:
            q = int(phaseof[e]) - 1
            # (1) internal phase changes: marking unchanged, gated by enabling
            D0off = md[e]['D0'] - np.diag(np.diag(md[e]['D0']))
            if np.count_nonzero(D0off) > 0:
                lev = []
                W = []
                for l in gate:
                    lev.append(l)
                    W.append(_placemat(l, md[e], 0.0, int(domain[l]), False))
                lev.append(q)
                W.append(sp.csr_matrix(D0off))
                events.append({'a': md[e]['trans'], 'b': md[e]['mode'], 'lev': lev, 'W': W})
            # (2) firings: marking moved, phase redrawn through D1
            lev = []
            W = []
            for l in touched:
                lev.append(l)
                W.append(_placemat(l, md[e], netm[e, l], int(domain[l]), l == degl))
            lev.append(q)
            W.append(sp.csr_matrix(md[e]['D1']))
            events.append({'a': md[e]['trans'], 'b': md[e]['mode'], 'lev': lev, 'W': W})

    levelkind = np.concatenate([np.ones(L, dtype=int), 2 * np.ones(Q, dtype=int)])
    desc = {
        'K': K,
        'domain': domain,
        'events': events,
        'levelkind': levelkind,
    }
    if winv is not None:
        desc['invariant'] = {
            'weights': np.concatenate([winv, np.zeros(Q)]),
            'value': float(vinv),
        }

    # ---- descriptive information
    placenames = [str(sn.nodenames[places[pp]]) for pp in range(P)]
    classnames = [str(sn.classnames[k]) for k in range(R)]
    levelname = [''] * K
    for pp in range(P):
        for k in range(R):
            levelname[pp * R + k] = '%s.%s' % (placenames[pp], classnames[k])
    for e in range(E):
        if phaseof[e] > 0:
            levelname[phaseof[e] - 1] = 'phase(%s.m%d)' % (
                str(sn.nodenames[md[e]['trans']]), md[e]['mode'] + 1)
    info = {
        'places': places,
        'placenames': placenames,
        'classnames': classnames,
        'levelkind': levelkind,
        'levelname': levelname,
        'modes': md,
        'nnodes': nnodes,
        'nclasses': R,
        'descriptor': descriptor,
        'init': init,
        'mdd': mdd,
        'phaseof': phaseof,
        'everDisabled': ever_disabled,
        'nplacelevels': L,
    }

    if verbose:
        line_printf('\nSPN -> MDD: %d places x %d classes = %d place levels, %d phase levels\n'
                    % (P, R, L, Q))
        line_printf('  modes = %d, |S| = %d, bounds = %s\n'
                    % (E, mdd.cardinality(), np.array2string(bound.astype(int))))
    return mdds, desc, info


# ---------------------------------------------------------------------------
def _placemat(l, mde, net, d, applydegree):
    """Local matrix of one mode at place level l.

    Move the count by net, but only from local states that satisfy this level's
    enabling and inhibition. applydegree scales each row by the enabling degree
    min(floor(m/enab), srv), i.e. the number of concurrently firing servers. It
    is carried by the single enabling level of a multi-server mode; with one
    server the degree is 1 and the indicator below is already the whole story.
    """
    i = np.arange(d)
    ok = (i >= mde['enab'][l]) & (i < mde['inhib'][l])
    j = i + int(net)
    ok = ok & (j >= 0) & (j <= d - 1)
    val = np.ones(d)
    if applydegree and mde['enab'][l] > 0:
        val = np.minimum(np.floor(i / mde['enab'][l]), mde['srv'])
    return sp.coo_matrix((val[ok], (i[ok], j[ok])), shape=(d, d)).tocsr()


# ---------------------------------------------------------------------------
def _init_marking(model, sn, places, R):
    """Tokens per (place, class) at time zero.

    Taken from the model state when set, and from the reference station
    otherwise.
    """
    P = len(places)
    init = np.zeros(P * R)
    anyset = False
    nodes = model.get_nodes()   # python spells getNodeByIndex as a list accessor
    for pp in range(P):
        node = nodes[places[pp]]
        st = getattr(node, '_state', None)
        if st is None and hasattr(node, 'getState'):
            st = node.getState()
        if st is not None and np.size(st) > 0:
            anyset = True
            stv = np.ravel(np.asarray(st, dtype=float))
            for k in range(min(R, stv.size)):
                init[pp * R + k] = stv[k]
    if not anyset or np.all(init == 0):
        njobs = np.ravel(np.asarray(sn.njobs, dtype=float))
        refstat = np.ravel(np.asarray(sn.refstat, dtype=int))
        node_to_station = np.ravel(np.asarray(sn.nodeToStation, dtype=int))
        for k in range(R):
            ref = int(refstat[k]) if k < refstat.size else 0
            for pp in range(P):
                if node_to_station[places[pp]] == ref:
                    init[pp * R + k] = njobs[k]
    return init


# ---------------------------------------------------------------------------
def _place_invariant(netm, init):
    """A strictly positive place invariant w, i.e. w >= 0 with netm @ w = 0.

    Scanning a null-space BASIS for a positive vector is not enough, and the
    fork-join net is the counterexample: its two minimal-support invariants are
    (1,1,0,1) and (1,0,1,1), neither positive, while their sum (2,1,1,2) is.
    Which basis a codebase's null space routine returns then decides whether the
    net is accepted, which is how MATLAB and the JAR came to disagree on it. So
    the non-negative generators are computed directly, by Farkas' algorithm on
    [netm.T | I] -- the same construction spn_sinvariants uses on the incidence
    matrix -- and summed.
    """
    L = netm.shape[1]
    if np.all(np.abs(netm.sum(axis=1)) < 1e-12):
        return np.ones(L), float(np.sum(init))
    if netm.size == 0:
        return None, None
    E = netm.shape[0]
    M = netm.T.copy()                      # L x E, row l is level l's column
    B = np.eye(L)                          # the combination each row stands for
    for e in range(E):
        keep = np.nonzero(np.abs(M[:, e]) < 1e-12)[0]
        pos = np.nonzero(M[:, e] > 1e-12)[0]
        neg = np.nonzero(M[:, e] < -1e-12)[0]
        rowsM = [M[keep, :]]
        rowsB = [B[keep, :]]
        for a in pos:
            for b in neg:
                c = -M[b, e] * M[a, :] + M[a, e] * M[b, :]
                d = -M[b, e] * B[a, :] + M[a, e] * B[b, :]
                g = _gcd_vec(np.concatenate([c, d]))
                if g > 0:
                    c = c / g
                    d = d / g
                rowsM.append(c[None, :])
                rowsB.append(d[None, :])
        M = np.vstack(rowsM)
        B = np.vstack(rowsB)
        M, B = _minimal_support(M, B)
    if B.shape[0] == 0:
        return None, None
    s = B.sum(axis=0)
    if np.all(s > 1e-9):
        w = s / s[s > 1e-9].min()
        return w, float(w @ init)
    return None, None


def _gcd_vec(v):
    """gcd of the integral entries, 0 when any entry is not an integer."""
    from math import gcd
    g = 0
    for x in v:
        if abs(x - round(x)) > 1e-9:
            return 0
        g = gcd(g, abs(int(round(x))))
    return g


def _minimal_support(M, B):
    """Drop every row whose support strictly contains another's, which is what
    leaves the minimal supports and stops the pair expansion from blowing up."""
    n = B.shape[0]
    supp = np.abs(B) > 1e-12
    drop = np.zeros(n, dtype=bool)
    for i in range(n):
        if drop[i]:
            continue
        for j in range(n):
            if i == j or drop[j]:
                continue
            if np.all(supp[j] <= supp[i]) and np.any(supp[i] > supp[j]):
                drop[i] = True
                break
    keep = ~drop
    return M[keep, :], B[keep, :]


# ---------------------------------------------------------------------------
def _next(s, md, netm, phaseof, domain, L):
    """Successors of s, plus a flag per mode saying it was disabled here."""
    E = len(md)
    m = np.asarray(s[:L], dtype=float)
    T = []
    disabled = np.zeros(E, dtype=bool)
    for e in range(E):
        if np.any(m < md[e]['enab']) or np.any(m >= md[e]['inhib']):
            disabled[e] = True
            continue
        if phaseof[e] == 0:
            t = list(s)
            newm = m + netm[e, :]
            if np.all(newm >= 0) and np.all(newm <= domain[:L] - 1):
                for l in range(L):
                    t[l] = int(newm[l])
                T.append(tuple(t))
        else:
            q = int(phaseof[e]) - 1
            ph = int(s[q])
            D0row = md[e]['D0'][ph, :].copy()
            D0row[ph] = 0.0
            for j in np.nonzero(D0row != 0)[0]:
                t = list(s)
                t[q] = int(j)
                T.append(tuple(t))
            newm = m + netm[e, :]
            if np.all(newm >= 0) and np.all(newm <= domain[:L] - 1):
                for j in np.nonzero(md[e]['D1'][ph, :] != 0)[0]:
                    t = list(s)
                    for l in range(L):
                        t[l] = int(newm[l])
                    t[q] = int(j)
                    T.append(tuple(t))
    return T, disabled
