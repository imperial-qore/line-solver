"""Exact product-form analysis under Krzesinski state-dependent routing.

Krzesinski, A. E., "Multiclass Queueing Networks with State-Dependent Routing",
Performance Evaluation 7(2):125-143, 1987. The joint distribution is eq. (16);
the coefficients xi are those of Section 3.2, obtained from the
state-independent part of the routing matrix.
"""

import time

import numpy as np

from ...api.pfqn.sdr import pfqn_sdr, pfqn_sdrmva, pfqn_sdrvisits
from ...api.sn.network_struct import SchedStrategy

__all__ = ['solver_nc_sdr_analyzer']


def solver_nc_sdr_analyzer(sn, options=None):
    """Return (QN, UN, RN, TN, CN, XN, lG, runtime, iter, method)."""
    tstart = time.time()
    M = int(sn.nstations)
    R = int(sn.nclasses)
    sdr = sn.sdr

    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    if np.any(np.isinf(njobs)):
        raise ValueError('state-dependent routing is defined for closed networks '
                         'only; the model has an open class')
    if int(getattr(sn, 'nchains', R)) != R:
        raise ValueError('state-dependent routing does not support class switching: '
                         'the product form of Krzesinski (1987) is stated over closed '
                         'chains whose customers keep their class')
    if int(sn.nstateful) != M:
        raise ValueError('state-dependent routing requires every stateful node to be '
                         'a station: the product form is over queue lengths, and a '
                         'stateless node holds none')

    rates = np.asarray(sn.rates, dtype=float)
    S = np.zeros((M, R))
    ok = np.isfinite(rates) & (rates > 0)
    S[ok] = 1.0 / rates[ok]

    Ntot = int(np.sum(njobs))
    alpha = np.ones((M, max(1, Ntot)))
    for i in range(M):
        sched = int(sn.sched[i]) if sn.sched is not None else -1
        if sched == int(SchedStrategy.INF):
            alpha[i, :] = np.arange(1, alpha.shape[1] + 1)
        else:
            c = float(sn.nservers[i]) if sn.nservers is not None else 1.0
            if np.isfinite(c) and c > 1:
                alpha[i, :] = np.minimum(np.arange(1, alpha.shape[1] + 1), c)
    lld = getattr(sn, 'lldscaling', None)
    if lld is not None and np.size(lld) > 0:
        lld = np.atleast_2d(np.asarray(lld, dtype=float))
        for i in range(min(M, lld.shape[0])):
            for k in range(min(alpha.shape[1], lld.shape[1])):
                if lld[i, k] > 0:
                    alpha[i, k] *= lld[i, k]

    # A BCMP center served FCFS must hold one rate for every chain; eq. (16)
    # admits chain-dependent rates only at the symmetric disciplines
    for i in range(M):
        if sn.sched is not None and int(sn.sched[i]) == int(SchedStrategy.FCFS):
            act = [r for r in range(R) if njobs[r] > 0 and S[i, r] > 0]
            if len(set(np.round(S[i, act], 12))) > 1:
                raise ValueError('station %d is FCFS with chain-dependent service times, '
                                 'which has no BCMP product form; use PS, LCFSPR or INF, '
                                 'or equalize the service times' % i)

    rt = np.asarray(sn.rt, dtype=float)
    P = np.zeros((M, M, R))
    for r in range(R):
        for i in range(M):
            isf = int(sn.stationToStateful[i])
            for j in range(M):
                jsf = int(sn.stationToStateful[j])
                P[i, j, r] = rt[isf * R + r, jsf * R + r]
    xi = pfqn_sdrvisits(sdr, P)

    # 'sdr' evaluates the product form (16) exactly by state enumeration, which
    # is general in the branch topology; 'sdr.mva' runs the paper's Section 4
    # MVA and convolution, which costs O(J T M (V_1...V_J)^2) instead of the
    # state-space size but requires single-centre branches. Both are exact.
    method = 'sdr'
    if options is not None and str(getattr(options, 'method', '')).lower() == 'sdr.mva':
        method = 'sdr.mva'
    if method == 'sdr.mva':
        QN, TN, UN, RN, lG = pfqn_sdrmva(S, xi, njobs.astype(int), sdr, alpha)
    else:
        QN, TN, UN, RN, _, lG, _, _ = pfqn_sdr(S, xi, njobs.astype(int), sdr, alpha)

    XN = np.zeros(R)
    CN = np.zeros(R)
    for r in range(R):
        ref = int(sn.refstat[r]) if sn.refstat is not None else 0
        XN[r] = TN[ref, r]
        if XN[r] > 0:
            CN[r] = njobs[r] / XN[r]

    return QN, UN, RN, TN, CN, XN, lG, time.time() - tstart, 1, method
