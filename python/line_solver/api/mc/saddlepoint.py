"""
Saddlepoint approximation of the counting process of a Markovian arrival process.

Pr{N(t)=k}, the probability that the MAP (D0,D1) records exactly k events in
(0,t], obtained by steepest-descent inversion of the counting generating
function instead of by forming the k-th superdiagonal block of expm(t*X).

The counting generating function is the matrix exponential

    sum_k P(k,t) z^k = expm(t*(D0 + z*D1)),

so the cumulant generating function of N(t) is

    eta(theta) = spectral abscissa of A(theta) = D0 + exp(theta)*D1,

the Perron root of an irreducible Metzler matrix: real, simple, strictly convex
in theta, with eta(0)=0 and eta'(0)=lambda. Inverting by steepest descent gives
Daniels (1954),

    Pr{N(t)=k} ~ g(theta*) * exp(t*eta(theta*) - k*theta*)
                           / sqrt(2*pi*t*eta''(theta*)),

with the saddle theta* solving eta'(theta*) = k/t and g the amplitude of the
Perron projection, g(theta) = (pi0*v)*(u*1), u and v the left and right Perron
vectors normalised by u*v = 1.

THE EXPANSION PARAMETER IS K2 = t*eta''(theta*), THE VARIANCE OF THE COUNT, not
its mean and not t. Measured error laws, constants flat to two digits over
Erlang orders 1..8 and horizons 10..160:

    err('daniels') = 0.083 / K2        err('daniels2') = 0.017 / K2**2

For a renewal Erlang(r) the count variance rate is lambda/r, so K2 = lambda*t/r
and an Erlang-4 at t=50 is as accurate as a Poisson at t=12.5: low variability
shrinks the parameter, it does not break the method. Below K2 = 5 the expansion
is out of its regime and the call warns.

This is an asymptotic method, not a quadrature: use it for rare-event and
large-deviation coefficients, where k/t is away from lambda or where the
probability underflows. For the bulk of the transient distribution, i.e. every
block k=0..N-1 at once at moderate t, uniformization (ctmc_uniformization,
ctmc_foxglynn) is both exact and faster.

ATTRIBUTION. The first-order form is Daniels (1954). The amplitude g and the
whole 'daniels2' bracket are NOT a rederivation: they are Jensen, "Saddlepoint
Expansions for Sums of Markov Dependent Variables on a Continuous State Space",
Probab. Th. Rel. Fields 89, 1991, Eq. (4.4) with the coefficients on p.191. His
gamma_0(s) = (sum_i c_i)(sum_i r_i P(Y_0=i)) is exactly g under his own
normalisation sum_i r_i c_i = 1, and expanding his
alpha_0 + (1/n){-alpha_3/2 + alpha_4/8 - 5*alpha_5/24} reproduces
g*(1 + lam4/8 - 5*lam3^2/24) - g''/(2*K2) + g'*K3/(2*K2^2) term for term; his
Theorem 4.1 gives the O(n^-2) error measured here as 0.017/K2^2. Jensen works
with discrete-n sums over a Markov chain, so the continuous-time MAP counting
process is that result transcribed, n -> t and the kernel eigenvalue -> the
Perron root of D0+exp(theta)*D1.

Key algorithms:
    ctmc_saddlepoint: Pr{N(t)=k} by the first- or second-order saddlepoint
"""

import math
import numpy as np
from scipy.linalg import expm
from typing import Optional, Sequence, Tuple

from line_solver.api.io import line_error, line_warning

# Below this value of K2 = t*eta''(theta*) the expansion is out of its regime.
# Do NOT threshold on lambda*t: for Erlang(r) the count variance rate is
# lambda/r, so K2 = lambda*t/r, and lambda*t over-warns on Poisson-like
# processes while under-warning on low-variability ones.
K2_MIN = 5.0

_METHODS = {
    'daniels2': (2, True), 'sp2': (2, True),
    'daniels': (1, True), 'sp1': (1, True),
    'plain': (1, False), 'bare': (1, False),
}


class PerronState:
    """Perron root of A(theta) with its first two derivatives and the amplitude."""

    __slots__ = ('eta', 'deta', 'd2eta', 'ampl')

    def __init__(self, eta: float, deta: float, d2eta: float, ampl: float):
        self.eta = eta
        self.deta = deta
        self.d2eta = d2eta
        self.ampl = ampl


def _perronstate(D0: np.ndarray, D1: np.ndarray, pi0: np.ndarray, th: float) -> PerronState:
    """
    Perron root of A(th) = D0 + exp(th)*D1 with deta, d2eta and the amplitude.

    Only EIGENVALUES are taken from the eigensolver; the Perron vectors come
    from bordered solves, the idiom ctmc_solve already uses. That keeps all four
    codebases on one algorithm: neither the JAR (commons-math hands back Schur
    blocks, not eigenvectors, as soon as a complex pair appears) nor C++
    (eig.h exposes values only) can supply a left eigenvector.
    """
    n = D0.shape[0]
    W = math.exp(th) * D1                     # A'(th) = A''(th) = exp(th)*D1
    A = D0 + W
    eta = float(np.max(np.linalg.eigvals(A).real))
    Ashift = A - eta * np.eye(n)
    rhs = np.zeros(n)
    rhs[n - 1] = 1.0
    # (A-eta*I)v = 0 with the last row replaced by sum(v)=1. A row may be
    # dropped because A-eta*I is a singular irreducible M-matrix, every proper
    # principal submatrix of which is nonsingular.
    M = Ashift.copy()
    M[n - 1, :] = 1.0
    v = np.linalg.solve(M, rhs)
    # u(A-eta*I) = 0 by the same construction on the transpose
    Mt = Ashift.T.copy()
    Mt[n - 1, :] = 1.0
    u = np.linalg.solve(Mt, rhs)
    u = u / (u @ v)                           # u*v = 1 fixes the residual scale
    deta = float(u @ W @ v)
    # First-order eigenvector perturbation (A-eta*I)v' = (eta'*I-W)v taken with
    # u*v'=0; the bordered system is nonsingular because the Perron root of an
    # irreducible Metzler matrix is simple.
    B = np.zeros((n + 1, n + 1))
    B[:n, :n] = Ashift
    B[:n, n] = v
    B[n, :n] = u
    r = np.zeros(n + 1)
    r[:n] = (deta * np.eye(n) - W) @ v
    vp = np.linalg.solve(B, r)[:n]
    d2eta = deta + 2.0 * float(u @ W @ vp)
    ampl = float((pi0 @ v) * (u @ np.ones(n)))
    return PerronState(eta, deta, d2eta, ampl)


def _solvesaddle(D0, D1, pi0, r: float, th0: float, thmin: float, thmax: float) -> Tuple[float, int]:
    """
    Saddle of the counting cumulant generating function at rate r, the root of
    eta'(th) = r. eta' is continuous and strictly increasing from 0 to +inf, so
    the root exists and is unique for every r>0; it is bracketed by geometric
    expansion from th0 and refined by Newton on log(eta'), safeguarded by
    bisection.
    """
    TOL = 1e-13
    MAXIT = 200

    def deriv(x):
        return _perronstate(D0, D1, pi0, x).deta

    th = min(max(th0, thmin), thmax)
    d1 = deriv(th)
    lo = hi = th
    dlo = dhi = d1
    step = 1.0
    while dlo > r:
        hi, dhi = lo, dlo
        lo = lo - step
        if lo <= thmin:
            lo = thmin
            dlo = deriv(lo)
            if dlo > r:
                line_error('ctmc_saddlepoint',
                           'The rate k/t=%g is below the representable range of eta\'.' % r)
            break
        dlo = deriv(lo)
        step *= 2.0
    step = 1.0
    while dhi < r:
        lo, dlo = hi, dhi
        hi = hi + step
        if hi >= thmax:
            hi = thmax
            dhi = deriv(hi)
            if dhi < r:
                line_error('ctmc_saddlepoint',
                           'The rate k/t=%g is above the representable range of eta\'.' % r)
            break
        dhi = deriv(hi)
        step *= 2.0
    th = min(max(th, lo), hi)
    logr = math.log(r)
    iters = 0
    for it in range(1, MAXIT + 1):
        iters = it
        si = _perronstate(D0, D1, pi0, th)
        f = math.log(si.deta) - logr
        if abs(f) <= TOL:
            break
        if f > 0.0:
            hi = th
        else:
            lo = th
        thn = th - f * si.deta / si.d2eta
        if not math.isfinite(thn) or thn <= lo or thn >= hi:
            thn = 0.5 * (lo + hi)
        if abs(thn - th) <= TOL * max(1.0, abs(th)):
            th = thn
            break
        th = thn
    return th, iters


def ctmc_saddlepoint(D0, D1=None, t=None, k=None, method: Optional[str] = None,
                     pi0: Optional[Sequence[float]] = None):
    """
    Saddlepoint approximation of Pr{N(t)=k} for the MAP counting process.

    Parameters
    ----------
    D0 : array (K,K), or the MAP pair (D0,D1) / [D0,D1], in which case the
        remaining arguments shift left by one.
    D1 : array (K,K), nonnegative. D0+D1 must be an irreducible generator.
    t : float or array. Time horizon, broadcast against k.
    k : int or array. Event count, a nonnegative integer, broadcast against t.
    method : 'daniels2' (default) second-order saddlepoint, error O(1/K2**2);
        'daniels' first order with the Perron amplitude, error O(1/K2);
        'plain' the bare first-order form with the amplitude set to 1.
    pi0 : array (K,). Initial phase distribution; the stationary distribution
        of D0+D1 if None.

    Returns
    -------
    p : ndarray. Approximation of Pr{N(t)=k}.
    logp : ndarray. Its natural logarithm, evaluated without forming p, so it
        stays accurate below the smallest positive double.
    theta : ndarray. The saddle theta*, -inf where k=0.
    info : dict with per-point arrays eta, deta, d2eta, d3eta, d4eta, ampl,
        corr, k2 (the expansion parameter), iter, exact, and the scalar lambda.
    """
    # MAP pair form (D0,D1): every later argument sits one position to the left.
    # method defaults to None rather than to the method name precisely so that
    # this single rebinding is unambiguous.
    if isinstance(D0, (tuple, list)) and len(D0) == 2 and np.ndim(D0[0]) == 2:
        D0, D1, t, k, method, pi0 = D0[0], D0[1], D1, t, k, method
    if method is None:
        method = 'daniels2'

    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    nph = D0.shape[0]
    if D0.shape[1] != nph or D1.shape != (nph, nph):
        line_error('ctmc_saddlepoint', 'D0 and D1 must be square matrices of the same order.')
    if np.any(D1 < 0):
        line_error('ctmc_saddlepoint', 'D1 must be nonnegative.')
    Q = D0 + D1
    if np.max(np.abs(Q.sum(axis=1))) > 1e-8 * max(1.0, float(np.max(np.abs(Q)))):
        line_error('ctmc_saddlepoint',
                   'D0+D1 must be an infinitesimal generator (zero row sums).')
    maxrate = float(np.max(np.abs(D1)))
    if maxrate <= 0.0:
        line_error('ctmc_saddlepoint',
                   'D1 has no counted transitions, the counting process is identically zero.')

    key = str(method).strip().lower()
    if key not in _METHODS:
        line_error('ctmc_saddlepoint',
                   "Unknown method '%s', expected daniels2, daniels or plain." % key)
    order, useampl = _METHODS[key]

    if pi0 is None:
        # stationary distribution of D0+D1
        M = np.vstack([Q.T, np.ones(nph)])
        b = np.zeros(nph + 1)
        b[-1] = 1.0
        pi0 = np.linalg.lstsq(M, b, rcond=None)[0]
    pi0 = np.asarray(pi0, dtype=float).reshape(-1)
    if pi0.size != nph:
        line_error('ctmc_saddlepoint', 'pi0 must have one entry per phase.')
    if abs(pi0.sum() - 1.0) > 1e-8:
        line_error('ctmc_saddlepoint', 'pi0 must sum to one.')

    t_arr = np.atleast_1d(np.asarray(t, dtype=float))
    k_arr = np.atleast_1d(np.asarray(k, dtype=float))
    if t_arr.size == 1 and k_arr.size > 1:
        t_arr = np.repeat(t_arr, k_arr.size).reshape(k_arr.shape)
    elif k_arr.size == 1 and t_arr.size > 1:
        k_arr = np.repeat(k_arr, t_arr.size).reshape(t_arr.shape)
    elif t_arr.shape != k_arr.shape:
        line_error('ctmc_saddlepoint', 't and k must be scalars or arrays of the same size.')
    if np.any(t_arr < 0):
        line_error('ctmc_saddlepoint', 'The horizon t must be nonnegative.')
    if np.any(k_arr < 0) or np.any(k_arr != np.round(k_arr)):
        line_error('ctmc_saddlepoint', 'The count k must be a nonnegative integer.')

    shape = t_arr.shape
    p = np.zeros(shape)
    logp = np.full(shape, -np.inf)
    theta = np.full(shape, -np.inf)
    info = {nm: np.full(shape, np.nan) for nm in
            ('eta', 'deta', 'd2eta', 'd3eta', 'd4eta', 'ampl', 'corr', 'k2')}
    info['iter'] = np.zeros(shape, dtype=int)
    info['exact'] = np.zeros(shape, dtype=bool)

    onesvec = np.ones(nph)
    # exp(theta) multiplies D1, so the saddle is confined to the range over
    # which A(theta) is representable; never active for a feasible k/t
    thmax = math.log(np.finfo(float).max / 1e6) - math.log(maxrate)
    thmin = math.log(np.finfo(float).tiny * 1e6) - math.log(maxrate)
    info['lambda'] = _perronstate(D0, D1, pi0, 0.0).deta

    tf = t_arr.reshape(-1)
    kf = k_arr.reshape(-1)
    pf = p.reshape(-1)
    logpf = logp.reshape(-1)
    thetaf = theta.reshape(-1)
    flat = {nm: info[nm].reshape(-1) for nm in info if isinstance(info[nm], np.ndarray)}

    # Sorting by the rate k/t lets each Newton solve warm-start from the
    # previous saddle, the saddle being a monotone function of that rate alone
    rate = np.where(tf > 0, kf / np.where(tf > 0, tf, 1.0), 0.0)
    worst_k2 = np.inf
    worst_at = (0.0, 0)
    thprev = 0.0

    for i in np.argsort(rate, kind='stable'):
        ti = float(tf[i])
        ki = int(kf[i])
        if ti == 0.0:
            # No time has elapsed, so the count is zero with probability one
            flat['exact'][i] = True
            if ki == 0:
                pf[i] = 1.0
                logpf[i] = 0.0
            continue
        if ki == 0:
            # The saddle runs off to -inf; the exact value is one matrix
            # exponential of the taboo generator and costs no more than a step
            # of the approximation itself
            flat['exact'][i] = True
            pf[i] = float(pi0 @ expm(ti * D0) @ onesvec)
            logpf[i] = math.log(pf[i]) if pf[i] > 0 else -np.inf
            continue

        th, iters = _solvesaddle(D0, D1, pi0, ki / ti, thprev, thmin, thmax)
        thprev = th
        thetaf[i] = th
        flat['iter'][i] = iters
        s = _perronstate(D0, D1, pi0, th)
        flat['eta'][i] = s.eta
        flat['deta'][i] = s.deta
        flat['d2eta'][i] = s.d2eta

        K2 = ti * s.d2eta
        flat['k2'][i] = K2
        if K2 < worst_k2:
            worst_k2 = K2
            worst_at = (ti, ki)
        if not (K2 > 0.0):
            line_error('ctmc_saddlepoint',
                       'The cumulant generating function is not strictly convex at the saddle '
                       "(t=%g, k=%d): eta''=%g. D0+D1 is probably reducible." % (ti, ki, s.d2eta))
        base = ti * s.eta - ki * th - 0.5 * math.log(2.0 * math.pi * K2)
        ampl = s.ampl if useampl else 1.0
        flat['ampl'][i] = s.ampl

        if order == 1:
            corr = ampl
        else:
            # The higher cumulants and the derivatives of the amplitude come
            # from central differences of the analytic eta'' and g, both of
            # which carry full precision at each evaluation point
            h = 1e-3 * max(1.0, abs(th))
            sp = _perronstate(D0, D1, pi0, th + h)
            sm = _perronstate(D0, D1, pi0, th - h)
            d3 = (sp.d2eta - sm.d2eta) / (2.0 * h)
            d4 = (sp.d2eta - 2.0 * s.d2eta + sm.d2eta) / (h * h)
            flat['d3eta'][i] = d3
            flat['d4eta'][i] = d4
            K3 = ti * d3
            K4 = ti * d4
            lam3sq = K3 * K3 / (K2 ** 3)
            lam4 = K4 / (K2 * K2)
            if useampl:
                gp = (sp.ampl - sm.ampl) / (2.0 * h)
                gpp = (sp.ampl - 2.0 * s.ampl + sm.ampl) / (h * h)
            else:
                gp = 0.0
                gpp = 0.0
            # Steepest descent to O(1/K2), Jensen (1991) Eq. (4.4): the Daniels
            # bracket on the amplitude, plus the two terms the amplitude
            # contributes through its own curvature along the contour
            corr = ampl * (1.0 + lam4 / 8.0 - 5.0 * lam3sq / 24.0) \
                - gpp / (2.0 * K2) + gp * K3 / (2.0 * K2 * K2)
            if corr <= 0.0:
                line_warning('ctmc_saddlepoint',
                             'The second-order correction is nonpositive at t=%g, k=%d; the '
                             'expansion has broken down, returning the first-order value.'
                             % (ti, ki))
                corr = ampl
        flat['corr'][i] = corr
        logpf[i] = base + math.log(corr)
        pf[i] = math.exp(logpf[i])

    if worst_k2 < K2_MIN:
        # Once per call, not once per point: a vectorised call spans hundreds of
        # counts and the caller needs the worst one, not a page of repetitions
        line_warning('ctmc_saddlepoint',
                     "K2 = t*eta''(theta*) = %.2f at t=%g, k=%d is below %g, so the saddlepoint "
                     'expansion is outside its asymptotic regime there and the result is '
                     'unreliable (expect a relative error near %.0e). K2 is the variance of the '
                     'count, not its mean: a low-variability process needs a longer horizon than '
                     'its rate suggests. Take the exact value from the block chain instead.'
                     % (worst_k2, worst_at[0], worst_at[1], K2_MIN, 0.017 / worst_k2 ** 2))

    return p, logp, theta, info
