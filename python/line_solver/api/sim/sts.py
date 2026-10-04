"""
Standardized time series areas of the batched quantile process.

The prefix quantiles are exact order statistics, obtained from a Fenwick tree
over the within-batch ranks that is advanced across all batches simultaneously,
so the cost is O(b m log m) with the m loop carrying only vectorized numpy
statements. This mirrors the MATLAB twin; the JAR uses a scalar Fenwick tree per
batch. All three compute the same exact order statistics, so results agree
exactly.

References:
    Original MATLAB: matlab/src/api/sim/sim_sts_quantile_areas.m
    C. Alexopoulos, D. Goldsman, A. Lolos, K. D. Dingec, J. R. Wilson,
    "Steady-State Quantile Estimation Using Standardized Time Series", 2020/2023.
    A. Lolos et al., Proc. Winter Simulation Conference, 2023, theorems 1 to 3.
"""

from math import ceil, sqrt
from typing import Dict, Sequence

import numpy as np

__all__ = ['sts_quantile_areas', 'DEFAULT_WEIGHT']

#: Constant STS weight function that normalizes the Brownian bridge area. The
#: requirement is that int_0^1 w(t)B(t)dt be standard normal for a standard
#: Brownian bridge B; for a constant w = c that variance is c^2/12.
DEFAULT_WEIGHT = sqrt(12.0)


def sts_quantile_areas(y: Sequence[float], b: int, m: int, p: float,
                       weight: float = DEFAULT_WEIGHT) -> Dict[str, object]:
    """
    Signed STS areas and variance-parameter estimators of the quantile process.

    Splits the ``b*m`` observations in ``y`` into b nonoverlapping batches of
    size m. With ``yhat_p(j,m)`` the empirical p-quantile of batch j and
    ``yhat_p(j,k)`` that of its first k observations, the STS process of batch j
    is

        T_{j,m}(k/m) = (k/sqrt(m)) (yhat_p(j,m) - yhat_p(j,k)),

    its signed area is ``A_p(w;j,m) = m^-1 sum_k w(k/m) T_{j,m}(k/m)``, and the
    three estimators of ``sigma_p^2 = lim n Var(ytilde_p(n))`` are

        A_p(w;b,m) = b^-1 sum_j A_p(w;j,m)^2                        (STS area)
        N_p(b,m)   = (b-1)^-1 m sum_j (yhat_p(j,m)-ytilde_p(n))^2   (NBQ)
        V_p(w;b,m) = [b A_p(w;b,m) + (b-1) N_p(b,m)] / (2b-1)       (combined)

    with ``ytilde_p(n)`` the full-sample empirical p-quantile. The first two have
    limiting chi-square laws on b and b-1 degrees of freedom and are
    asymptotically independent, so the combined estimator carries 2b-1 degrees of
    freedom and is about sqrt(2) less variable than either component.

    On i.i.d. Exp(1) data, where ``sigma_p^2 = p(1-p)/f(y_p)^2`` is exact, both
    A_p and N_p are unbiased to within 7%.

    Args:
        y: Exactly ``b*m`` observations, in time order
        b: Batch count, at least 1
        m: Batch size, at least 1
        p: Quantile probability in (0,1)
        weight: Constant STS weight function, nonzero

    Returns:
        Dict with keys ``areas``, ``bqe``, ``quantile``, ``Ap``, ``Np``, ``Vp``,
        ``b``, ``m``, ``n`` and ``analyzer``. ``Np`` and ``Vp`` are NaN when
        ``b < 2``, there being no between-batch degrees of freedom.

    Raises:
        ValueError: If the arguments are out of range or ``y`` is the wrong size.
    """
    if int(b) != b or b < 1:
        raise ValueError("The batch count b must be a positive integer, got %r" % (b,))
    if int(m) != m or m < 1:
        raise ValueError("The batch size m must be a positive integer, got %r" % (m,))
    if not 0.0 < p < 1.0:
        raise ValueError("p must lie in (0,1), got %r" % (p,))
    if weight == 0.0 or not np.isfinite(weight):
        raise ValueError("weight must be nonzero and finite, got %r" % (weight,))

    b = int(b)
    m = int(m)
    n = b * m
    v = np.asarray(y, dtype=float).ravel()
    if v.size != n:
        raise ValueError("y must hold exactly b*m = %d observations, got %d"
                         % (n, v.size))
    if not np.all(np.isfinite(v)):
        raise ValueError("The sample path must be finite")

    # column j is batch j in time order
    ym = v.reshape(m, b, order='F')
    order = np.argsort(ym, axis=0, kind='stable')
    sorted_vals = np.take_along_axis(ym, order, axis=0)

    # rank[i, j] is the 1-based position of ym[i, j] among its batch's values
    rank = np.empty((m, b), dtype=np.int64)
    rows = np.arange(1, m + 1, dtype=np.int64)[:, None]
    np.put_along_axis(rank, order, np.broadcast_to(rows, (m, b)), axis=0)

    bqe = sorted_vals[ceil(m * p) - 1, :].copy()

    # Fenwick tree per batch over the within-batch ranks, advanced in lockstep
    fen = np.zeros((m + 1, b), dtype=np.int64)
    cols = np.arange(b)
    top = 1 << (m.bit_length() - 1)
    acc = np.zeros(b)

    for k in range(1, m + 1):
        pos = rank[k - 1, :].copy()
        while True:
            active = pos <= m
            if not active.any():
                break
            idx = pos[active]
            fen[idx, cols[active]] += 1
            pos[active] = idx + (idx & (-idx))

        target = ceil(p * k)
        pos = np.zeros(b, dtype=np.int64)
        remaining = np.full(b, target, dtype=np.int64)
        step = top
        while step >= 1:
            cand = pos + step
            ok = cand <= m
            if ok.any():
                okcols = cols[ok]
                values = fen[cand[ok], okcols]
                move = values < remaining[ok]
                if move.any():
                    sel = okcols[move]
                    remaining[sel] -= values[move]
                    pos[sel] = cand[sel]
            step >>= 1

        acc += k * (bqe - sorted_vals[pos, cols])

    areas = weight * acc / (m * sqrt(m))
    quantile = float(np.partition(v, ceil(n * p) - 1)[ceil(n * p) - 1])

    ap = float(np.mean(areas ** 2))
    if b >= 2:
        np_est = float(m * np.sum((bqe - quantile) ** 2) / (b - 1))
        vp = (b * ap + (b - 1) * np_est) / (2 * b - 1)
    else:
        # a single batch carries no between-batch degrees of freedom; areas and
        # bqe stay valid and firquest pools them across replications instead
        np_est = float('nan')
        vp = float('nan')

    return {'areas': areas, 'bqe': bqe, 'quantile': quantile, 'Ap': ap,
            'Np': np_est, 'Vp': vp, 'b': b, 'm': m, 'n': n,
            'analyzer': 'sts_quantile_areas'}
