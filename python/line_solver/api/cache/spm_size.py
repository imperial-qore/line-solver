"""Ray (WKB) asymptotic expansion of the cost-capped cache normalizing constant.

Native-Python port of matlab/src/api/cache/cache_spm_size.m (and
jar/.../jline/api/cache/Cache_spm_size.java).

Approximates what ``cache_erec(gamma, m, sigma, k)`` computes exactly, in the
SAME normalization, so the two are interchangeable.  This is the item-size
extension of ``retrieval_rayint``, which carries the size-free expansion; call
that one when there are no storage costs.

Writing ``E = prod_j m_j! * H``, the size-free recursion::

    E(m,n) = E(m,n-1) + sum_j gamma_{n,j} m_j E(m-1_j,n-1),  E(0,0)=1

relaxes to ``H ~ exp(phi/eps)`` with ``n = y/eps``, ``m_j = x_j/eps``, whose
eikonal ``e^{phi_y} = 1 + sum_j gamma_j(y) e^{-phi_j}`` carries the ray
constants ``xi_j = e^{-phi_j}``.  With per-item storage costs ``sigma_i`` and
per-list cost caps ``k_j`` the recursion gains the cost coordinate::

    E(m,k) = E_i(m,k) + sum_j m_j gamma_ij E_i(m-1_j, k-sigma_i 1_j),

so the shift ``1_j`` becomes ``e_j(y) = (1_j, s(y) 1_j)`` in the enlarged space
``X = (x,kappa)`` and the eikonal picks up the size tilt::

    e^{phi_y} = 1 + sum_j gamma_j(y) e^{-phi_{x_j} - s(y) phi_{kappa_j}},

with the second family of ray constants ``zeta_j = e^{-phi_{kappa_j}}``.  The
rays integrate to the discrete saddle point of the product generating function::

    sum_{m,k} H(m,k) prod_j z_j^{m_j} w_j^{k_j}
        = prod_i ( 1 + sum_j gamma_ij z_j w_j^{sigma_i} ),

namely, with ``D_i = 1 + sum_j gamma_ij xi_j zeta_j^{sigma_i}`` and
``Psi = sum_i log D_i``::

    m_j = sum_i gamma_ij xi_j zeta_j^{sigma_i} / D_i
    k_j = sum_i sigma_i gamma_ij xi_j zeta_j^{sigma_i} / D_i
    log H(m,k) ~ Psi - sum_j m_j log xi_j - sum_j k_j log zeta_j
                 - (d/2) log(2 pi) - (1/2) log det grad^2 Psi

where ``d`` is the number of saddle coordinates and, with
``pi_ij = gamma_ij xi_j zeta_j^{sigma_i} / D_i`` and
``Q^i_{jl} = delta_{jl} pi_ij - pi_ij pi_il``::

    grad^2 Psi = sum_i [1; sigma_i] [1; sigma_i]' (x) Q^i

Setting ``zeta_j = 1`` recovers the size-free expansion exactly.

CAPS ARE CUMULATIVE.  ``cache_erec`` sums over the states of cost AT MOST
``k_j``, so this function does the same by default (``'atmost'``).  The shadow
price ``eta_j = log zeta_j <= 0`` obeys complementary slackness: a list whose
unconstrained mean cost already meets its cap is SLACK, keeps ``zeta_j = 1`` and
drops out of the saddle, which then degenerates continuously to the size-free
expansion; a list whose cap BINDS sits at ``eta_j < 0``, and the states below
the boundary decay geometrically with ratio ``zeta_j``, contributing the
amplitude factor ``1/(1-zeta_j)``.  Pass ``'exact'`` to obtain instead the
constant resolving the cost exactly at ``k_j``, which is the raw Laplace
formula above with no such factor.

SIZE DIVERSITY IS REQUIRED.  The Hessian integrand
``[1;sigma_i][1;sigma_i]' (x) Q^i`` has rank ``h``, not ``2h``, so
``grad^2 Psi`` is nonsingular only if the sizes actually vary.  This is not an
artefact: with a single item size the cost of list j is ``sigma*m_j``
identically and the cap carries no information.  That case is detected and
answered exactly rather than passed to a singular saddle.  If the sizes share a
common divisor the cost lives on a sublattice; the sizes and caps are divided
through by their gcd, which is an exact reduction and removes the corresponding
lattice factor.

OCCUPANCY.  ``out.pij`` is the saddle occupancy
``pi_il = gamma_il xi_l zeta_l^{sigma_i} / D_i`` and ``out.K`` its per-list cost.
These are EXACT-COST quantities: the saddle conditions are
``sum_i pi_ij = m_j`` and ``sum_i sigma_i pi_ij = k_j``, so ``out.K`` equals the
cap exactly on every binding list.  Under cumulative caps the true mean cost is
strictly below the cap; use ``cache_cost`` and ``cache_prob_erec`` for that.

ACCURACY.  The expansion is ``O(1/n)`` at fixed occupancy.  With a
well-separated cap the observed error in ``log E`` is around ``1e-2`` at
``n = 200`` and halves at each doubling of ``n``.  It degrades as a binding
``zeta_j`` approaches 1, i.e. in the transition between the binding and slack
regimes, where the geometric resummation ``1/(1-zeta_j)`` is no longer sharp;
``out.zeta`` and ``out.binding`` report where the saddle sits and a warning is
raised inside that region.
"""
import warnings
from math import gcd
from types import SimpleNamespace

import numpy as np
from scipy.special import gammaln


def cache_spm_size(gamma, m, sigma, k, costmode='atmost'):
    """Ray expansion of the cost-capped cache normalizing constant.

    Parameters
    ----------
    gamma : ndarray (n, h)
        Item popularity probabilities (access factors).
    m : ndarray (h,)
        Cache list capacities, non-negative integers.
    sigma : ndarray (n,)
        Item storage costs (sizes), positive integers.
    k : ndarray (h,)
        Per-list storage cost caps, integers.
    costmode : {'atmost', 'exact'}
        ``'atmost'`` (default) matches ``cache_erec``; ``'exact'`` resolves the
        cost exactly at ``k``.

    Returns
    -------
    E : float
        Normalizing constant, same normalization as ``cache_erec`` (may
        overflow; use ``logE``).
    logE : float
        Natural logarithm of ``E``.
    out : SimpleNamespace
        ``xi``, ``zeta``, ``binding``, ``pij``, ``K``, ``phi``,
        ``logdetSigma``, ``span``, ``method``, ``relerrEst``, ``iter``.
    """
    gamma = np.asarray(gamma, dtype=float)
    if gamma.ndim != 2 or gamma.size == 0:
        raise ValueError("cache_spm_size: the access factors must be a non-empty n x h matrix.")
    n0, h0 = gamma.shape
    m = np.asarray(m, dtype=float).ravel()
    if m.size != h0:
        raise ValueError("cache_spm_size: the capacity vector must have one entry per cache "
                         "list (%d given, %d expected)." % (m.size, h0))
    if np.any(m < 0) or np.any(np.abs(m - np.rint(m)) > 0):
        raise ValueError("cache_spm_size: list capacities must be non-negative integers.")
    if sigma is None or k is None:
        raise ValueError("cache_spm_size: the item sizes and the cost caps are both required. "
                         "Use retrieval_rayint for the size-free expansion.")
    sigma = np.asarray(sigma, dtype=float).ravel()
    k = np.asarray(k, dtype=float).ravel()
    if sigma.size == 0 or k.size == 0:
        raise ValueError("cache_spm_size: the item sizes and the cost caps are both required. "
                         "Use retrieval_rayint for the size-free expansion.")
    if sigma.size != n0:
        raise ValueError("cache_spm_size: the item size vector must have one entry per item.")
    if k.size != h0:
        raise ValueError("cache_spm_size: the cost cap vector must have one entry per cache list.")
    if np.any(sigma <= 0) or np.any(np.abs(sigma - np.rint(sigma)) > 0):
        raise ValueError("cache_spm_size: item sizes must be positive integers.")
    if np.any(np.abs(k - np.rint(k)) > 0):
        raise ValueError("cache_spm_size: storage cost caps must be integers.")
    costmode = str(costmode).lower()
    if costmode not in ('atmost', 'exact'):
        raise ValueError("cache_spm_size: the cost mode must be 'atmost' or 'exact' "
                         "('%s' given)." % costmode)
    capped = True   # cleared below when a single item size makes the cap uninformative

    out = SimpleNamespace(xi=np.zeros(h0), zeta=np.ones(h0), binding=np.zeros(h0, dtype=bool),
                          pij=np.hstack([np.ones((n0, 1)), np.zeros((n0, h0))]),
                          K=np.zeros(h0), phi=np.nan, logdetSigma=np.nan, span=1,
                          method='', relerrEst=np.nan, iter=0)

    # --- boundaries, matching cache_erec ---
    if m.sum() > n0 or np.any(k < 0):
        out.method = 'boundary'
        return 0.0, -np.inf, out
    if m.sum() == 0:
        out.method = 'boundary'
        if costmode == 'exact' and np.any(k > 0):
            return 0.0, -np.inf, out
        return 1.0, 0.0, out

    # --- items that can never be cached and lists of zero capacity drop out ---
    alive = gamma.sum(axis=1) > 0
    live = m > 0
    G = gamma[np.ix_(alive, live)]
    mk = m[live]
    hk = int(mk.size)
    n = int(G.shape[0])
    sg = sigma[alive]
    kk = k[live]
    if mk.sum() > n:
        out.method = 'boundary'
        return 0.0, -np.inf, out
    if mk.sum() == n:
        raise ValueError("cache_spm_size: the expansion requires sum(m) < n; at sum(m) = n the "
                         "saddle point escapes to infinity. Use cache_erec for a full cache.")

    # --- exact reductions on the cost lattice ---
    span = 0
    for s in sg:
        span = gcd(span, int(s))
    out.span = span
    if costmode == 'exact' and np.any(np.mod(kk, span) != 0):
        out.method = 'lattice'      # unreachable off the sublattice
        return 0.0, -np.inf, out
    sg = sg / span
    kk = np.floor(kk / span)
    # per-list feasibility: the m_j cheapest (dearest) reachable items bound the cost
    for j in range(hk):
        idx = np.flatnonzero(G[:, j] > 0)
        if idx.size < mk[j]:
            out.method = 'boundary'
            return 0.0, -np.inf, out
        srt = np.sort(sg[idx])
        mj = int(mk[j])
        if srt[:mj].sum() > kk[j]:
            out.method = 'boundary'
            return 0.0, -np.inf, out
        if costmode == 'exact' and srt[srt.size - mj:].sum() < kk[j]:
            out.method = 'boundary'
            return 0.0, -np.inf, out
    # a single item size makes the cost of list j equal to sigma*m_j identically,
    # so the cap carries no information and the 2h saddle is singular (rank h)
    if np.all(sg == sg[0]):
        cost = sg[0] * mk
        feasible = np.all(cost == kk) if costmode == 'exact' else np.all(cost <= kk)
        if not feasible:
            out.method = 'uniform-size'
            return 0.0, -np.inf, out
        capped = False              # fall through to the size-free expansion
        out.method = 'uniform-size'

    # --- saddle point ---
    if capped:
        th, et, bind, it, P, D = _saddle_cost(G, mk, sg, kk, costmode == 'atmost')
        ix = np.flatnonzero(bind)
        dof = hk + int(ix.size)
        phi = float(np.sum(np.log(D)) - mk @ th - kk[ix] @ et[ix])
        logdet = _logdet(_hessian(P, sg, ix))
        logH = phi - 0.5 * dof * np.log(2 * np.pi) - 0.5 * logdet
        if costmode == 'atmost' and ix.size > 0:
            logH -= float(np.sum(np.log1p(-np.exp(et[ix]))))   # resummation below the cap
        if not out.method:
            out.method = 'spm-size'
    else:
        th, it, P, D = _saddle(G, mk)
        et = np.zeros(hk)
        bind = np.zeros(hk, dtype=bool)
        phi = float(np.sum(np.log(D)) - mk @ th)
        logdet = _logdet(_hessian(P, np.zeros(n), np.empty(0, dtype=int)))
        logH = phi - 0.5 * hk * np.log(2 * np.pi) - 0.5 * logdet
        if not out.method:
            out.method = 'spm'

    logE = float(logH + gammaln(m + 1.0).sum())     # back to the cache_erec normalization
    with np.errstate(over='ignore'):
        E = float(np.exp(logE))

    # --- ray quantities, reported on the original item and list indexing ---
    out.xi[live] = np.exp(th)
    out.zeta[live] = np.exp(et / span)
    out.binding[live] = bind
    pij = np.zeros((n0, h0 + 1))
    pij[np.ix_(alive, np.hstack([[False], live]))] = P
    pij[:, 0] = 1.0 - pij[:, 1:].sum(axis=1)
    out.pij = pij
    out.K = sigma @ pij[:, 1:]
    out.phi = phi
    out.logdetSigma = float(logdet)
    out.iter = int(it)
    out.relerrEst = float(0.14 * (1.0 / mk.min() + 1.0 / (n - mk.sum())))
    if min(mk.min(), n - mk.sum()) < 2:
        warnings.warn("cache_spm_size: the smallest occupancy is %d, so the expansion is only "
                      "qualitative here (estimated relative error %.0f%%); cache_erec is exact."
                      % (int(min(mk.min(), n - mk.sum())), 100 * out.relerrEst), RuntimeWarning)
    if costmode == 'atmost' and bind.any() and np.max(np.exp(et[bind])) > 0.8:
        warnings.warn("cache_spm_size: a binding cost cap has zeta = %.3f, i.e. it sits in the "
                      "transition between the binding and slack regimes where the geometric "
                      "resummation below the cap is not sharp; the reported relative error does "
                      "not cover it." % float(np.max(np.exp(et[bind]))), RuntimeWarning)
    return E, logE, out


def _occupancy(G, sg, th, et):
    """pi_ij = gamma_ij xi_j zeta_j^{sigma_i} / D_i, D_i = 1 + sum_j (that numerator)."""
    A = G * np.exp(th[None, :] + np.outer(sg, et))
    D = 1.0 + A.sum(axis=1)
    return A / D[:, None], D


def _theta0(G, tgt):
    n = G.shape[0]
    slack = max(1.0 - tgt.sum() / n, 1e-9)
    return np.log(np.maximum(tgt, 1e-12) / np.maximum(G.sum(axis=0) * slack, 1e-12))


def _obj(G, sg, th, et, tgtm, tgtk, bind):
    _, D = _occupancy(G, sg, th, et)
    return float(np.sum(np.log(D)) - tgtm @ th - tgtk[bind] @ et[bind])


def _hessian(P, sg, ix):
    """grad^2 Psi in (theta, eta), restricted to the free eta coordinates ix.

    Each block is sum_i w_i Q^i with Q^i = diag(pi_i) - pi_i pi_i' and
    w = 1, sigma, sigma^2.
    """
    def qw(w):
        return np.diag((w[:, None] * P).sum(axis=0)) - P.T @ (w[:, None] * P)

    htt = qw(np.ones(P.shape[0]))
    if ix.size == 0:
        return htt
    hte = qw(np.asarray(sg, dtype=float))
    hee = qw(np.asarray(sg, dtype=float) ** 2)
    return np.block([[htt, hte[:, ix]], [hte[ix, :], hee[np.ix_(ix, ix)]]])


def _newton_step(H, g):
    d = -np.linalg.solve(H, g)
    if not np.all(np.isfinite(d)):
        raise ValueError("cache_spm_size: the saddle-point Newton step is not finite. With item "
                         "sizes this is the rank-h degeneracy of the size-tilted Hessian: the "
                         "sizes must genuinely vary for the cost coordinate to carry information.")
    return d


def _damp(d):
    step = 1.0
    while np.max(np.abs(step * d)) > 2.0:   # keep the tilts within a factor e^2 per iteration
        step /= 2.0
    return step


def _logdet(H):
    try:
        L = np.linalg.cholesky(0.5 * (H + H.T))
    except np.linalg.LinAlgError:
        raise ValueError("cache_spm_size: the saddle-point Hessian is not positive definite; "
                         "the ray map is singular here. With item sizes this happens when the "
                         "sizes do not vary over the items the cache can hold, in which case the "
                         "cost cap carries no information.")
    return 2.0 * float(np.sum(np.log(np.diag(L))))


def _saddle(G, tgt):
    """Size-free saddle: Newton on theta = log xi for sum_i gamma_ij xi_j / D_i = m_j."""
    n, h = G.shape
    th = _theta0(G, tgt)
    zeros_n = np.zeros(n)
    zeros_h = np.zeros(h)
    noix = np.empty(0, dtype=int)
    it = 0
    for it in range(1, 201):
        P, D = _occupancy(G, zeros_n, th, zeros_h)
        g = P.sum(axis=0) - tgt
        if np.max(np.abs(g)) <= 1e-12 * max(1.0, float(np.max(np.abs(tgt)))):
            break
        d = _newton_step(_hessian(P, zeros_n, noix), g)
        th = th + _damp(d) * d
    P, D = _occupancy(G, zeros_n, th, zeros_h)
    return th, it, P, D


def _saddle_cost(G, tgtm, sg, tgtk, cumulative):
    """Cost-constrained saddle.

    Minimises the convex dual f(theta,eta) = sum_i log D_i - m.theta - k.eta over
    eta <= 0 when the caps are cumulative, so that complementary slackness selects
    the binding lists; over all of R^{2h} when the cost is resolved exactly.
    """
    n, h = G.shape
    th = _theta0(G, tgtm)
    et = np.zeros(h)
    bind = np.ones(h, dtype=bool)
    tol = 1e-12 * max(1.0, float(np.max(np.abs(tgtm))), float(np.max(np.abs(tgtk))))
    it = 0
    for it in range(1, 201):
        P, D = _occupancy(G, sg, th, et)
        gth = P.sum(axis=0) - tgtm
        get = sg @ P - tgtk
        if cumulative:
            # at eta_j = 0 the cap binds when the mean cost exceeds it
            bind = (et < 0) | (get > 0)
        ix = np.flatnonzero(bind)
        g = np.concatenate([gth, get[ix]])
        if np.max(np.abs(g)) <= tol:
            break
        d = _newton_step(_hessian(P, sg, ix), g)
        step = _damp(d)
        fcur = _obj(G, sg, th, et, tgtm, tgtk, bind)
        # Backtrack until the dual decreases. The slack is essential, not cosmetic:
        # Newton reaches the floating-point floor of f in a handful of steps, and a
        # strict test then rejects every step and halves to zero without converging.
        ftol = 1e-12 * (1.0 + abs(fcur))
        for _ in range(40):
            thn = th + step * d[:h]
            etn = et.copy()
            etn[ix] = et[ix] + step * d[h:]
            if cumulative:
                etn = np.minimum(etn, 0.0)
            if _obj(G, sg, thn, etn, tgtm, tgtk, bind) <= fcur + ftol:
                break
            step /= 2.0
        moved = max(float(np.max(np.abs(thn - th))), float(np.max(np.abs(etn - et))))
        th, et = thn, etn
        if moved <= 1e-13:          # the iterate can no longer move: at the floor
            break
    P, D = _occupancy(G, sg, th, et)
    if cumulative:
        bind = et < 0
    return th, et, bind, it, P, D
