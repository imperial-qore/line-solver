"""
Exact MAP/MAP/1 solution for a single-class, single-server, open
Source -> FCFS Queue -> Sink model.

The decomposition methods (dec.source/dec.mmap) approximate this queue: they fit
the arrival to a simpler process or treat the service as a renewal phase-type,
discarding autocorrelation. When the arrival or the service is a genuinely
correlated (non-renewal) MAP, this routine returns the exact mean queue length
via the matrix-geometric MAP/MAP/1 solver (q_ct_map_map_1).

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np

from line_solver.lib.thirdparty.qmam import q_ct_map_map_1
from line_solver.api.mam.map_analysis import map_lambda
from . import MAMResult
from .transient_qbd import _find_source_queue, _proc_ph_entry, _read_map, _is_renewal_map


def _is_markovian_map(D0, D1):
    """True when (D0, D1) is a genuine MAP rather than a RAP or an ME process.

    A MAP has non-negative off-diagonal rates in D0, non-negative rates in D1,
    and (D0 + D1) is an infinitesimal generator (zero row sums). A RAP or an ME
    violates the sign conditions while still defining a valid point process.
    """
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    ns = D0.shape[0]
    tol = 1e-9 * max(1.0, float(np.max(np.abs(np.concatenate((D0.flatten(), D1.flatten()))))))
    off_diag = D0[~np.eye(ns, dtype=bool)]
    return (bool(np.all(off_diag >= -tol))
            and bool(np.all(D1 >= -tol))
            and bool(np.all(np.abs(np.sum(D0 + D1, axis=1)) <= tol)))


def solver_mam_mapmap1_exact(sn):
    """Return an exact MAMResult when sn is a single-class single-server open
    MAP/MAP/1 queue with a correlated (non-renewal) arrival or service;
    otherwise return None so the caller falls back to the decomposition methods.
    """
    if sn.nclasses != 1:
        return None
    njobs = np.asarray(sn.njobs).flatten()
    if not np.isinf(njobs[0]):
        return None

    source_idx, queue_idx = _find_source_queue(sn)
    if source_idx is None or queue_idx is None:
        return None
    # The arrival seen by the queue equals the source MAP only when the queue is
    # fed directly by the source: require source and queue to be the only stations.
    for i in range(sn.nstations):
        if i not in (source_idx, queue_idx):
            return None
    nservers = np.asarray(sn.nservers).flatten()
    if int(nservers[queue_idx]) != 1:
        return None

    Da0, Da1 = _read_map(_proc_ph_entry(sn, source_idx))
    Ds0, Ds1 = _read_map(_proc_ph_entry(sn, queue_idx))
    if Da1 is None or Ds1 is None:
        return None

    # see _kb/06-solver-catalog.md (MAM: "MAP/MAP/1 exact fast-path") --
    # genuine MAPs only; RAP/ME fall back to the RAP/RAP/1 QBD
    if not _is_markovian_map(Da0, Da1) or not _is_markovian_map(Ds0, Ds1):
        return None

    # Only fire when the decomposition methods are inexact: a genuinely
    # correlated (non-renewal) MAP arrival or service.
    if _is_renewal_map(Da0, Da1) and _is_renewal_map(Ds0, Ds1):
        return None

    lam = map_lambda(Da0, Da1)
    mu = map_lambda(Ds0, Ds1)
    if not (lam < mu):
        return None   # unstable or degenerate: leave to the fallback path

    res = q_ct_map_map_1(Da0, Da1, Ds0, Ds1)
    ql = np.asarray(res.queue_length).flatten()
    EN = float(np.sum(np.arange(len(ql)) * ql))
    rho = lam / mu

    M = sn.nstations
    K = 1
    QN = np.zeros((M, K)); UN = np.zeros((M, K)); RN = np.zeros((M, K)); TN = np.zeros((M, K))
    TN[source_idx, 0] = lam
    TN[queue_idx, 0] = lam
    QN[queue_idx, 0] = EN
    UN[queue_idx, 0] = rho
    RN[queue_idx, 0] = EN / lam
    CN = np.array([[EN / lam]])
    XN = np.array([[lam]])
    return MAMResult(QN=QN, UN=UN, RN=RN, TN=TN, CN=CN, XN=XN,
                     totiter=1, method="exact.mapmap1", runtime=0.0)
