"""
Near-immediate feedback elimination for the robust queueing network analyzer.

Native Python twin of matlab/src/api/npfqn/npfqn_feedback_elim.m, implementing
Section 4 of W. Whitt and W. You (2022), A robust queueing network analyzer
based on indices of dispersion, Naval Research Logistics 69(1), 36-56.
"""

from typing import Any, Dict, Optional, Sequence

import numpy as np


def npfqn_feedback_elim(P, rho: Sequence[float], cs2: Optional[Sequence[float]] = None,
                        lambda_: Optional[Sequence[float]] = None,
                        immediateOnly: bool = False) -> Dict[str, Any]:
    """
    Eliminate near-immediate feedback from an open queueing network.

    WHY FEEDBACK BREAKS DECOMPOSITION. A parametric decomposition treats the
    arrival stream at each station as if it were renewal. Feedback destroys that
    badly: a customer that leaves a busy station and comes straight back arrives
    exactly when the station is busy, so the flow is strongly correlated with the
    queue it feeds. The fix is not to model the correlation but to REMOVE the
    feedback, by folding the repeated visits into the service time.

    THE TRANSFORMATION. With feedback probability ``p`` at a station, a customer
    is served a geometric number of times, so the effective service is
    ``S_p = sum_{i=1}^{N} S_i`` with ``N`` geometric of mean ``1/(1-p)``. Hence

    * effective mean service ``E[S]/(1-p)``,
    * effective service SCV ``p + (1-p)cs^2`` (eq. 37 and the line after it),
    * fresh arrival rate ``lambda(1-p)``,
    * per-visit waiting time ``= (1-p)`` times the wait in the modified system.

    The modified system has the same heavy-traffic limits for queue length,
    workload, waiting time and external departures, so the elimination is
    asymptotically exact rather than merely plausible.

    NEAR-IMMEDIATE, NOT JUST IMMEDIATE. Feedback rarely returns a customer in one
    hop. What matters is whether it returns WITHOUT PASSING A BUSIER STATION: a
    detour through a station of lower traffic intensity is fast on the time scale
    of the busy station, so it behaves like immediate feedback. The probability
    computed here is therefore the probability of returning to station ``i``
    through stations of strictly smaller ``rho`` only.

    Args:
        P: routing matrix, substochastic, ``P[i][j]`` from station i to j
        rho: traffic intensity of each station, which fixes what counts as
            "near-immediate"
        cs2: service SCV of each station; the modified SCVs are returned when given
        lambda_: arrival rate of each station; the modified rates are returned
            when given
        immediateOnly: keep only the self-loops ``P[i][i]``, i.e. immediate
            feedback in the strict sense of Section 4.1

    Returns:
        Dict with ``feedbackProb`` (p-hat per station), ``modifiedScv``,
        ``modifiedRates``, ``modifiedRouting`` (the routing with the eliminated
        feedback removed and the remaining rows renormalized), and
        ``visitInflation`` (``1/(1-p)``, the mean visits per customer).

    References:
        W. Whitt, W. You (2022). A robust queueing network analyzer based on
        indices of dispersion. Naval Research Logistics 69(1), 36-56, Section 4.
    """
    P = np.asarray(P, dtype=float)
    m = P.shape[0]
    if P.shape != (m, m):
        raise ValueError('The routing matrix must be square.')
    if np.any(P < -1e-12) or np.any(P.sum(axis=1) > 1 + 1e-9):
        raise ValueError('The routing matrix must be substochastic.')
    rho = np.asarray(rho, dtype=float)
    if rho.size != m:
        raise ValueError('One traffic intensity per station is required.')

    phat = np.zeros(m)
    for i in range(m):
        if immediateOnly:
            phat[i] = P[i, i]
            continue
        # Stations a customer may pass through on a near-immediate return: those
        # NOT MORE loaded than i. A detour through a busier station is not fast
        # on the time scale of station i, so it is not near-immediate; one
        # through a station of equal load is, which is why the test is <= and
        # not <. This is the cloud of eqs. (3.8)-(3.9) with H = {i}, and the
        # same one solver_rqna applies -- the two must not drift.
        idx = [j for j in range(m) if j != i and rho[j] <= rho[i] + 1e-9]
        ret = P[i, i]
        if idx:
            Q = P[np.ix_(idx, idx)]
            r = P[np.ix_(idx, [i])].ravel()
            # (I-Q)^-1 r is the probability of eventually reaching i from each
            # allowed station without leaving the allowed set.
            reach = np.linalg.solve(np.eye(len(idx)) - Q, r)
            ret += float(P[i, idx] @ reach)
        phat[i] = min(max(ret, 0.0), 1.0 - 1e-12)

    result: Dict[str, Any] = {
        'feedbackProb': phat,
        'visitInflation': 1.0 / (1.0 - phat),
    }
    if cs2 is not None:
        cs2 = np.asarray(cs2, dtype=float)
        if cs2.size != m:
            raise ValueError('One service SCV per station is required.')
        # eq. (37): the geometric sum of service times.
        result['modifiedScv'] = phat + (1.0 - phat) * cs2
    if lambda_ is not None:
        lam = np.asarray(lambda_, dtype=float)
        if lam.size != m:
            raise ValueError('One arrival rate per station is required.')
        result['modifiedRates'] = lam * (1.0 - phat)

    # The reduced network. For IMMEDIATE feedback the reduction is exact and
    # unambiguous: drop the self-loop and renormalize the rest of the row, since
    # a customer that does not feed back goes where it would have gone anyway.
    # For NEAR-IMMEDIATE feedback the return path runs through other stations, so
    # there is no such row-local reduction; the elimination then applies to the
    # SERVICE description at the station, which is what modifiedScv and
    # modifiedRates carry, and this field is returned as the immediate-feedback
    # reduction only, for reference.
    Pmod = P.copy()
    for i in range(m):
        loop = P[i, i]
        if loop <= 0:
            continue
        Pmod[i, i] = 0.0
        rest = Pmod[i].sum()
        if rest > 0:
            Pmod[i] *= (P[i].sum() - loop) / rest
    result['modifiedRouting'] = Pmod
    result['reductionExact'] = bool(immediateOnly or np.allclose(phat, np.diag(P)))
    return result
