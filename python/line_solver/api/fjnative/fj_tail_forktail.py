"""
ForkTail black-box tail-latency approximation for fork-join requests.

Port of matlab/src/api/fj/fj_tail_forktail.m and
matlab/src/api/fj/fj_mg1_respt_moments.m; JAR twin
jline.api.fj.FJ_tail_forktail.

Reference:
    M. Nguyen, S. Alesawi, N. Li, H. Che, H. Jiang, "ForkTail: A Black-Box
    Fork-Join Tail Latency Prediction Model for User-Facing Datacenter
    Workloads", ACM HPDC 2018, pp. 206-217.
"""

import numpy as np
from scipy.optimize import brentq
from scipy.special import digamma, polygamma

FINE_TOL = 1e-8
COARSE_TOL = 1e-3


def fj_mg1_respt_moments(lambda_, ES, ES2, ES3):
    """
    Mean and variance of the M/G/1 FCFS response time, the white-box inputs of
    ForkTail, from the first three moments of the service time:

        E[T] = E[S]*(1 + rho/(1-rho) * (1+SCV_S)/2)
        V[T] = E[W]^2 + lambda*E[S^3]/(3*(1-rho)) + E[S^2] - E[S]^2

    with rho = lambda*E[S] and E[W] = lambda*E[S^2]/(2*(1-rho)), the
    Pollaczek-Khinchine mean waiting time.

    Args:
        lambda_: arrival rate at the branch
        ES: first moment of the service time
        ES2: second moment of the service time
        ES3: third moment of the service time

    Returns:
        (ET, VT): response time mean and variance

    Raises:
        ValueError: if the branch is unstable
    """
    rho = lambda_ * ES
    if rho >= 1.0:
        raise ValueError(
            f"The branch is unstable (rho = {rho:g} >= 1); "
            "the response time moments do not exist.")
    if not np.isfinite(ES3):
        raise ValueError(
            "The service law has no finite third moment, so the ForkTail response "
            "time variance is undefined (a Pareto branch with shape <= 3, for "
            "instance). Use a service law with three finite moments.")
    scv_s = (ES2 - ES**2) / ES**2
    ET = ES * (1.0 + rho/(1.0-rho) * (1.0 + scv_s)/2.0)
    EW = lambda_ * ES2 / (2.0*(1.0-rho))
    VT = EW**2 + lambda_*ES3/(3.0*(1.0-rho)) + ES2 - ES**2
    return ET, VT


def ge_fit(ET, VT):
    """
    Match a generalized exponential law F_T(x) = (1-exp(-x/beta))^alpha on a
    mean and a variance.

    The squared coefficient of variation depends on the shape alone and
    decreases monotonically in it, so the shape is recovered by a scalar
    root-find on a logarithmic scale and the scale then follows in closed
    form. SCV = 1 is the exponential case alpha = 1, kept exact.

    Args:
        ET: mean of the branch response time
        VT: variance of the branch response time

    Returns:
        (alpha, beta): fitted shape and scale
    """
    if ET <= 0 or VT <= 0:
        raise ValueError("The task response time mean and variance must be positive.")
    scv = VT / ET**2
    if abs(scv - 1.0) < FINE_TOL:
        alpha = 1.0
    else:
        def residual(la):
            a = np.exp(la)
            m = digamma(a + 1.0) - digamma(1.0)
            return (polygamma(1, 1.0) - polygamma(1, a + 1.0)) / m**2 - scv

        lo, hi = -30.0, 30.0
        while residual(lo) < 0 and lo > -700.0:
            lo -= 30.0    # smaller shape -> larger SCV
        while residual(hi) > 0 and hi < 700.0:
            hi += 30.0    # larger shape -> smaller SCV
        alpha = float(np.exp(brentq(residual, lo, hi, xtol=1e-14, rtol=1e-14)))
    beta = ET / (digamma(alpha + 1.0) - digamma(1.0))
    return alpha, beta


def fj_tail_forktail(ET, VT, K=1, p=99, P=None):
    """
    Predict the p-th percentile of a fork-join request response time.

    Each branch is fitted with a generalized exponential law matched on the
    mean and variance of its task response time, and the request response time
    is the maximum over the branches, taken as the product of the branch CDFs
    (exact only for independent branches):

        F_X(x) = prod_i (1 - exp(-x/beta_i))^alpha_i
        homogeneous: x_p = -beta*log(1 - p^(1/(K*alpha)))
        random fanout: F_X(x) = sum_i P_i * (1 - exp(-x/beta))^(K_i*alpha)

    The approximation rests on the central limit theorem for G/G/m queues in
    heavy traffic, so it is a HIGH-LOAD result and under-predicts at low load,
    the more so the more variable the service. Prefer the FJ_codes route
    (fj_is_homogeneous) on the homogeneous MAP/PH/1 class, which is more
    accurate there; ForkTail covers the heterogeneous branches and mixed
    service laws that route rejects.

    Args:
        ET: mean task response time, a scalar (homogeneous branches) or one
            entry per branch
        VT: variance of the task response time, same shape as ET
        K: number of branches, or a sequence of distinct fanouts when the
            fanout is random; ignored when ET is a sequence
        p: percentile, a fraction in (0,1) or a percentage in (0,100)
        P: probabilities of the fanouts in K, required when K is a sequence

    Returns:
        (xp, alpha, beta): the predicted percentile and the fitted parameters
    """
    p = float(p)
    if p > 1.0:
        p = p / 100.0
    if p <= 0.0 or p >= 1.0:
        raise ValueError("The percentile must lie strictly between 0 and 1 (or 0 and 100).")

    ET = np.atleast_1d(np.asarray(ET, dtype=float)).ravel()
    VT = np.atleast_1d(np.asarray(VT, dtype=float)).ravel()
    if ET.size != VT.size:
        raise ValueError("ET and VT must have the same number of entries.")

    fits = [ge_fit(ET[i], VT[i]) for i in range(ET.size)]
    alpha = np.array([f[0] for f in fits])
    beta = np.array([f[1] for f in fits])

    Kv = np.atleast_1d(np.asarray(K)).ravel()

    if ET.size == 1 and Kv.size > 1:
        # random fanout: mix the homogeneous request laws over the fanout
        # distribution and invert numerically
        if P is None:
            raise ValueError("A vector of fanouts K needs a probability vector P of the same length.")
        Pv = np.atleast_1d(np.asarray(P, dtype=float)).ravel()
        if Pv.size != Kv.size:
            raise ValueError("A vector of fanouts K needs a probability vector P of the same length.")
        if np.any(Pv < 0) or abs(Pv.sum() - 1.0) > COARSE_TOL:
            raise ValueError("The fanout probabilities P must be non-negative and sum to one.")
        a, b = alpha[0], beta[0]

        def mixres(x):
            return float(np.sum(Pv * (1.0 - np.exp(-x/b))**(Kv*a))) - p

        xlo = -b * np.log(1.0 - p**(1.0/(Kv.min()*a)))
        xhi = -b * np.log(1.0 - p**(1.0/(Kv.max()*a)))
        if xlo == xhi:
            return float(xlo), a, b
        xp = brentq(mixres, min(xlo, xhi), max(xlo, xhi), xtol=1e-14, rtol=1e-14)
        return float(xp), a, b

    if ET.size == 1:
        # homogeneous: the product of K identical CDFs raises the shape to
        # K*alpha and inverts in closed form
        a, b = alpha[0], beta[0]
        xp = -b * np.log(1.0 - p**(1.0/(int(Kv[0])*a)))
        return float(xp), a, b

    # inhomogeneous: solve prod_i (1-exp(-x/beta_i))^alpha_i = p. F_X is
    # bounded above by any single branch CDF, so the request percentile is at
    # least the largest branch percentile
    logp = np.log(p)

    def residual(x):
        return float(np.sum(alpha * np.log1p(-np.exp(-x/beta)))) - logp

    xlo = float(np.max(-beta * np.log(1.0 - p**(1.0/alpha))))
    xhi = xlo
    while residual(xhi) < 0:
        xhi *= 2.0
        if not np.isfinite(xhi):
            raise RuntimeError("Could not bracket the ForkTail percentile.")
    while residual(xlo) > 0 and xlo > np.finfo(float).tiny:
        xlo /= 2.0
    xp = brentq(residual, xlo, xhi, xtol=1e-14, rtol=1e-14)
    return float(xp), alpha, beta


def forktail_percentiles(solver, percentiles, jobclass=None):
    """
    Fork-join request tail latency by ForkTail, driven from a solved model.

    Mirrors the MATLAB entry point @NetworkSolver/getPerctRespT.m with
    method='forktail'. The branch arrival rate is read from the solved
    throughputs and the service moments from the node distributions.

    The topology gate is strict, because the approximation is defined on a
    per-branch M/G/1 task response time: exactly one fork with a matching join,
    and every branch a single station feeding that join.

    Args:
        solver: a solved (or solvable) network solver exposing getAvg() and .model
        percentiles: percentiles, fractions in (0,1) or percentages in (0,100)
        jobclass: optional class index (0-based) or name to restrict the output

    Returns:
        (PercRT, PercTable): a list of per-class dicts and a pandas DataFrame
    """
    import pandas as pd
    from ...lang.base import NodeType

    model = solver.model
    sn = model.get_struct()
    nodetype = [int(nt) for nt in np.asarray(sn.nodetype).ravel()]
    forks = [i for i, nt in enumerate(nodetype) if nt == int(NodeType.FORK)]
    if not forks:
        raise ValueError("The forktail method requires a model with a Fork node.")
    if len(forks) > 1:
        raise ValueError(
            f"The forktail method supports a single fork-join pair; this model has {len(forks)} forks.")
    f = forks[0]

    fj = np.atleast_2d(np.asarray(sn.fj))
    join_candidates = np.where(fj[f, :] > 0)[0]
    if join_candidates.size == 0:
        raise ValueError("The fork node has no matching join; the request response time is undefined.")
    join_idx = int(join_candidates[0])

    conn = np.atleast_2d(np.asarray(sn.connmatrix))
    branches = [int(j) for j in np.where(conn[f, :] > 0)[0]]
    names = list(sn.nodenames) if getattr(sn, 'nodenames', None) is not None else [str(i) for i in branches]
    for b in branches:
        if not sn.isstation[b]:
            raise ValueError(f"Branch node {names[b]} is not a station; "
                             "the forktail method needs one queueing station per branch.")
        if conn[b, join_idx] == 0:
            raise ValueError(f"Branch station {names[b]} does not feed the join directly; "
                             "the forktail method needs one station per branch.")

    QN, UN, RN, TN = solver.getAvg()[:4]
    TN = np.atleast_2d(np.asarray(TN, dtype=float))
    UN = np.atleast_2d(np.asarray(UN, dtype=float))

    classnames = list(sn.classnames)
    if jobclass is None:
        classes = list(range(int(sn.nclasses)))
    elif isinstance(jobclass, str):
        classes = [classnames.index(jobclass)]
    else:
        classes = [int(jobclass)]

    pcts = np.atleast_1d(np.asarray(percentiles, dtype=float)).ravel()

    PercRT = []
    rows = []
    for r in classes:
        ET = np.zeros(len(branches))
        VT = np.zeros(len(branches))
        rho = np.zeros(len(branches))
        traverses = True
        for bi, b in enumerate(branches):
            ist = int(sn.nodeToStation[b])
            lam = TN[ist, r]
            if lam <= FINE_TOL:
                traverses = False
                break
            svc = model.get_nodes()[b].get_service(model.get_classes()[r])
            ES = svc.getMean()
            VS = svc.getVar()
            ES2 = VS + ES**2
            ES3 = svc.getSkew()*VS**1.5 + 3*ES*ES2 - 2*ES**3
            ET[bi], VT[bi] = fj_mg1_respt_moments(lam, ES, ES2, ES3)
            rho[bi] = UN[ist, r]
        if not traverses:
            continue
        if rho.max() < 0.5:
            import warnings
            warnings.warn(
                "ForkTail is a heavy-traffic approximation; the busiest branch is at utilization "
                f"{rho.max():.2f}, so the tail is likely under-predicted.")
        values = np.array([fj_tail_forktail(ET, VT, None, p)[0] for p in pcts])
        PercRT.append({'class': classnames[r],
                       'percentiles': pcts / 100.0 if np.any(pcts > 1) else pcts,
                       'values': values,
                       'method': 'forktail'})
        for pi, p in enumerate(pcts):
            rows.append({'JobClass': classnames[r],
                         'Percentile': p/100.0 if p > 1 else p,
                         'ResponseTime': values[pi]})
    return PercRT, pd.DataFrame(rows)
