"""Exact delayed-hit rate of a cache, as a transition reward over the generator.

Shared by the CTMC analyzer and by the ENV state-vector blend, which solves the
same generators one stage at a time. Port of MATLAB
matlab/src/solvers/CTMC/cache_delayed_hit_rate.m, twin of the JAR's static
SolverCTMC.delayedHitRate.
"""

import numpy as np

from ...sn import NodeType


def cache_delayed_hit_rate(sn, isf, infgen, state_space, pi, col_span=None):
    """Delayed-hit rate per ORIGINATING class of the cache at stateful index isf.

    A fetch of item i completes on exactly the transitions that clear block A
    bit i, and each such transition releases the block-B counts of item i as
    delayed hits, so the rate is a TRANSITION reward over the generator, not a
    state reward (the alternative arrival-rate identity lambda_i*phi_i is only
    PASTA-exact).

    ``col_span`` is the (start, end) column range of this stateful node inside
    the concatenated state vector. Pass it when the caller already knows it;
    otherwise it is derived from ``sn.space``, which a released struct no longer
    carries.

    Returns a length-nclasses vector, all zeros when the node carries no
    retrieval sub-system.
    """
    K = int(sn.nclasses)
    out = np.zeros(K)
    if infgen is None or state_space is None or pi is None:
        return out
    ind = int(sn.statefulToNode[isf])
    if sn.nodetype[ind] != NodeType.CACHE:
        return out
    npar = sn.nodeparam[ind] if sn.nodeparam is not None and ind in sn.nodeparam else None
    if npar is None:
        return out
    if int(getattr(npar, 'retrieval_system_capacity', 0) or 0) <= 0:
        return out

    from ...state.ctmc_ssg import cache_retrieval_class_map
    _, rc_items, rc_orig = cache_retrieval_class_map(sn, ind)
    if not rc_items:
        return out

    space = np.atleast_2d(np.asarray(state_space, dtype=float))
    pi = np.ravel(np.asarray(pi, dtype=float))
    if space.shape[0] != pi.shape[0]:
        return out
    Q = np.asarray(infgen.todense(), dtype=float) if hasattr(infgen, 'todense') \
        else np.asarray(infgen, dtype=float)
    if Q.shape[0] != space.shape[0]:
        return out

    # Column span of this stateful node inside the concatenated state vector.
    if col_span is not None:
        col_off, col_end = int(col_span[0]), int(col_span[1])
        width = col_end - col_off
    else:
        if sn.space is None or isf not in sn.space or sn.space[isf] is None:
            return out
        col_off = 0
        for jsf in range(int(sn.nstateful)):
            if jsf == isf:
                break
            w = 0
            if sn.space is not None and jsf in sn.space and sn.space[jsf] is not None:
                w = int(np.atleast_2d(np.asarray(sn.space[jsf])).shape[1])
            col_off += w
        width = int(np.atleast_2d(np.asarray(sn.space[isf])).shape[1])

    nitems = int(getattr(npar, 'nitems', 0))
    tcc = int(getattr(npar, 'total_cache_capacity', 0))
    lvs = width - (tcc + nitems + len(rc_items))
    if lvs < 0:
        return out
    a0 = col_off + lvs + tcc
    b0 = a0 + nitems
    if b0 + len(rc_items) > space.shape[1]:
        return out

    offdiag = Q - np.diag(np.diag(Q))
    for j, item in enumerate(rc_items):
        i = item - 1
        rows = np.where((space[:, a0 + i] != 0) & (space[:, b0 + j] > 0))[0]
        for rr in rows:
            nz = np.where(offdiag[rr] != 0)[0]
            completes = nz[space[nz, a0 + i] == 0]
            if completes.size == 0:
                continue
            out[rc_orig[j]] += (pi[rr] * space[rr, b0 + j]
                                * float(np.sum(offdiag[rr, completes])))
    return out
