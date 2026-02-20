"""
KPC-Toolbox Fitting Functions.

Native Python implementations of the KPC-Toolbox fitting algorithms for
Markovian Arrival Processes (MAPs) and Phase-Type (PH) distributions.

Based on the Kronecker Product Composition (KPC) method.

References:
    [1] G.Casale, E.Z.Zhang, E.Smirni. Trace Data Characterization and Fitting
        for Markov Modeling, Elsevier Performance Evaluation, 67(2):61-79,
        Feb 2010.
    [2] G.Casale, E.Z.Zhang, E.Smirni. KPC-Toolbox: Simple Yet Effective Trace
        Fitting Using Markovian Arrival Processes. in Proc. of QEST 2008,
        83-92, St.Malo, France, IEEE Press, September 2008.
"""

import numpy as np
import time
import warnings
from typing import Union, Tuple, Optional, List, Dict, Any
from dataclasses import dataclass, field
from scipy.optimize import minimize
from math import factorial

ArrayLike = Union[np.ndarray, list]

# Default tolerance for KPC fitting
KPCFIT_TOL = 1e-10


# ---------------------------------------------------------------------------
# Helper: scale a MAP to a desired mean (MATLAB map_scale semantics)
# ---------------------------------------------------------------------------
def _map_scale_to_mean(D0, D1, new_mean):
    """Scale a MAP (D0, D1) so that its mean inter-arrival time equals new_mean.

    This reproduces the MATLAB semantics of map_scale(MAP, NEWMEAN):
        ratio = map_mean(MAP) / NEWMEAN
        D0 *= ratio; D1 *= ratio; normalize
    """
    from ..mam import map_mean, map_normalize
    current_mean = map_mean(D0, D1)
    if current_mean <= 0 or not np.isfinite(current_mean) or new_mean <= 0:
        return D0.copy(), D1.copy()
    ratio = current_mean / new_mean
    D0_s = D0 * ratio
    D1_s = D1 * ratio
    return map_normalize(D0_s, D1_s)


@dataclass
class KpcfitTraceData:
    """Data structure for KPC fitting."""
    S: np.ndarray  # Original trace
    E: np.ndarray  # Moments E[X^k]
    AC: np.ndarray  # Autocorrelation coefficients for fitting
    ACFull: np.ndarray  # Full autocorrelation coefficients
    ACLags: np.ndarray  # Lags for AC
    BC: np.ndarray  # Bicovariance coefficients
    BCGridLags: np.ndarray  # Grid lags for BC
    BCLags: np.ndarray  # Lag combinations for BC


@dataclass
class KpcfitPhOptions:
    """Options for PH distribution fitting.

    Matches MATLAB kpcfit_ph_options fields:
        Verbose, Runs, MinNumStates, MaxNumStates, MinExactMom
    """
    verbose: bool = True
    runs: int = 5
    min_num_states: int = 2
    max_num_states: int = 32
    min_exact_mom: int = 3


@dataclass
class KpcfitResult:
    """Result of KPC fitting."""
    MAP: Optional[Tuple[np.ndarray, np.ndarray]]
    fac: float  # Autocorrelation fitting objective
    fbc: float  # Bicovariance fitting objective
    sub_maps: List[Tuple[np.ndarray, np.ndarray]]


def kpcfit_tol() -> float:
    """Return default tolerance for KPC fitting."""
    return KPCFIT_TOL


def kpcfit_version() -> str:
    """Return KPC-Toolbox version string."""
    return '0.3.2'


def logspacei(start: int, end: int, n: int) -> np.ndarray:
    """
    Generate logarithmically spaced integers.

    Args:
        start: Starting value
        end: Ending value
        n: Number of points

    Returns:
        Array of unique integers, logarithmically spaced
    """
    if start <= 0:
        start = 1
    if end <= start:
        return np.array([start], dtype=int)
    log_vals = np.logspace(np.log10(start), np.log10(end), n)
    return np.unique(np.round(log_vals).astype(int))


def kpcfit_init(trace: ArrayLike, ac_lags: ArrayLike = None,
                bc_grid_lags: ArrayLike = None, smooth: int = 0,
                max_moment: int = 3) -> KpcfitTraceData:
    """
    Prepare trace data for KPC fitting.

    Args:
        trace: Array of inter-arrival times
        ac_lags: Lags for autocorrelation (default: logarithmically spaced)
        bc_grid_lags: Lags for bicovariance grid (default: logarithmically spaced)
        smooth: Smoothing window size (0 = no smoothing)
        max_moment: Maximum moment order to compute

    Returns:
        KpcfitTraceData structure with preprocessed data
    """
    from ..trace import trace_acf, trace_joint, trace_bicov

    S = np.asarray(trace, dtype=np.float64).ravel()
    n = len(S)

    n_min_support_ac = 10

    # Default AC lags
    if ac_lags is None:
        ac_lags = logspacei(1, n // n_min_support_ac, 500)
    ac_lags = np.asarray(ac_lags, dtype=int)

    # Default BC grid lags
    if bc_grid_lags is None:
        bc_grid_lags = logspacei(1, max(ac_lags) if len(ac_lags) > 0 else 10, 5)
    bc_grid_lags = np.asarray(bc_grid_lags, dtype=int)

    # Compute moments
    E = np.zeros(max_moment)
    for j in range(max_moment):
        E[j] = np.mean(S ** (j + 1))

    # Compute autocorrelations
    print("init: computing moments from the trace")
    AC = trace_acf(S, ac_lags)
    ACFull = trace_acf(S, np.arange(1, n // n_min_support_ac + 1))

    # Apply smoothing if requested
    if smooth > 0:
        print("init: computing smoothed autocorrelations from the trace")
        try:
            from scipy.ndimage import uniform_filter1d
            AC = uniform_filter1d(AC, size=smooth, mode='nearest')
            ACFull = uniform_filter1d(ACFull, size=smooth, mode='nearest')
        except ImportError:
            pass  # Skip smoothing if scipy not available
    else:
        print("init: computing autocorrelations from the trace")

    # Truncate AC where values become negligible
    small_idx = np.where(np.abs(AC) < 1e-6)[0]
    if len(small_idx) > 0:
        posmax = ac_lags[small_idx[0]]
        keep_mask = ac_lags <= posmax
        ac_lags = ac_lags[keep_mask]
        AC = AC[keep_mask]
        bc_grid_lags = bc_grid_lags[bc_grid_lags <= posmax]

    print(f"Using {len(ac_lags)} AC Lags")

    # Compute bicovariances
    print("init: computing bicovariances from the trace")
    BC, BCLags = trace_bicov(S, bc_grid_lags)

    return KpcfitTraceData(
        S=S, E=E, AC=AC, ACFull=ACFull, ACLags=ac_lags,
        BC=BC, BCGridLags=bc_grid_lags, BCLags=BCLags
    )


def kpcfit_sub_eval_acfit(SCV: np.ndarray, G2: np.ndarray,
                           acf_lags: np.ndarray) -> Tuple[float, np.ndarray]:
    """
    Evaluate autocorrelation fit for given SCV and G2 parameters.

    Computes the autocorrelation function of a MAP composed of J MAP(2)s
    with the given SCVs and G2 decay rates, evaluated at the specified lags.

    Args:
        SCV: Squared coefficient of variation for each composing MAP
        G2: Autocorrelation decay rate for each composing MAP
        acf_lags: Lags at which to evaluate

    Returns:
        Tuple of (final_SCV, acf_coefficients)
    """
    J = len(G2)
    SCV = np.asarray(SCV).ravel()
    G2 = np.asarray(G2).ravel()
    acf_lags = np.asarray(acf_lags).ravel()

    SCVj = SCV[0]
    acf_coeff = 0.5 * (1 - 1 / SCV[0]) * (G2[0] ** acf_lags)

    for j in range(1, J):
        SCVj_1 = SCVj
        SCVj = (1 + SCVj) * (1 + SCV[j]) / 2 - 1
        r0j = 0.5 * (1 - 1 / SCV[j])
        X = SCV[j] * r0j * (G2[j] ** acf_lags)
        acf_coeff = (X + SCVj_1 * acf_coeff * (1 + X)) / SCVj

    return SCVj, acf_coeff


# ============================================================================
# Autocorrelation fitting (kpcfit_sub_acfit)
# ============================================================================

def kpcfit_sub_acfit(E: np.ndarray, SA: np.ndarray, SAlags: np.ndarray,
                     J: int, MaxIterAC: int = 300, MaxRunsAC: int = 50,
                     MaxResAC: int = 10,
                     AnimateAC: bool = False
                     ) -> Tuple[List[np.ndarray], List[np.ndarray], np.ndarray]:
    """
    Fit autocorrelation structure via multi-start optimization.

    Finds SCV and G2 parameters for J MAP(2) components that best reproduce
    the empirical autocorrelation SA at the specified lags.

    Args:
        E: Moments of the trace (E[0]=E[X], E[1]=E[X^2], ...)
        SA: Empirical autocorrelation values at SAlags
        SAlags: Lags at which SA is evaluated
        J: Number of MAP(2) components (log2 of total states)
        MaxIterAC: Maximum iterations per optimization run
        MaxRunsAC: Number of random restarts
        MaxResAC: Maximum number of results to return
        AnimateAC: (ignored in Python, MATLAB animation flag)

    Returns:
        Tuple of (resSCV, resG2, fobjAC) where each resSCV[i] / resG2[i]
        is an ndarray of length J.
    """
    MaxTimeAC = 10  # seconds
    MaxResAC = min(MaxResAC, MaxRunsAC)

    SCV = (E[1] - E[0] ** 2) / E[0] ** 2
    NSA = np.linalg.norm(SA, 2)  # normalization constant
    SAlags = np.asarray(SAlags).ravel()

    TOL = KPCFIT_TOL

    # -- helper: x -> (SCVj, G2j) ------
    def xtopar(x):
        SCVj = x[:J]
        G2j = x[J:]
        return SCVj, G2j

    # -- objective function ------
    def objfun(x):
        SCVj, G2j = xtopar(x)
        SCVJ, acfCoeff = kpcfit_sub_eval_acfit(SCVj, G2j, SAlags)
        f = np.linalg.norm(SA - acfCoeff, 1) / NSA + (SCVJ - SCV) ** 2 / SCV ** 2
        return f

    # -- inequality constraints  c(x) <= 0 ------
    def ineq_constraints(x):
        SCVj, G2j = xtopar(x)
        c = []
        # SCVj(1) >= 0.5
        c.append((0.5 - TOL) - SCVj[0])
        # SCVj(j) >= 1 for j>=2
        for j in range(1, J):
            c.append((1 + TOL) - SCVj[j])
        # G2j(j) < 1 and G2j(j) > 0
        for j in range(J):
            c.append(G2j[j] - (1 - TOL))
            c.append(TOL - G2j[j])
        return np.array(c)

    # -- stagnation callback ------
    class _StagnationCallback:
        def __init__(self, max_iter, max_time, tstart, f_best):
            self.stagnval = 0.0
            self.stagniter = 0
            self.max_iter = max_iter
            self.max_time = max_time
            self.tstart = tstart
            self.f_best = f_best
            self.niter = 0

        def __call__(self, xk):
            self.niter += 1
            fval = objfun(xk)
            # time check
            elapsed = time.time() - self.tstart
            if elapsed > self.max_time and self.niter > self.max_iter:
                raise StopIteration
            # stagnation check
            if self.stagnval == 0:
                self.stagnval = fval
            delta = abs(fval - self.stagnval) / max(abs(self.stagnval), 1e-30)
            if delta < 0.01:
                self.stagniter += 1
                if self.stagniter >= 100:
                    raise StopIteration
            else:
                self.stagniter = 0
            self.stagnval = fval
            # exceed-best check
            if self.niter % self.max_iter == 0 and self.niter > 1:
                if fval > self.f_best:
                    raise StopIteration

    # -- multi-start optimization ------
    fset = []
    xparamset = []
    f_best = np.inf

    lb = np.full(2 * J, TOL)
    ub = np.full(2 * J, np.inf)
    bounds = list(zip(lb, ub))

    constraints = {'type': 'ineq', 'fun': lambda x: -ineq_constraints(x)}

    for i in range(MaxRunsAC):
        x0 = np.concatenate([1 + np.random.rand(J), np.random.rand(J)])
        tstart = time.time()
        cb = _StagnationCallback(MaxIterAC, MaxTimeAC, tstart, f_best)

        print(f"acfit: run {i+1} of {MaxRunsAC} ", end='')
        t0 = time.time()
        try:
            res = minimize(objfun, x0, method='SLSQP',
                           bounds=bounds,
                           constraints=constraints,
                           options={'maxiter': MaxIterAC, 'ftol': TOL,
                                    'disp': False},
                           callback=cb)
            x = res.x
            f = res.fun
        except StopIteration:
            x = cb.__dict__.get('_last_x', x0)
            f = objfun(x0)
            try:
                f = objfun(x)
            except Exception:
                pass

        dt = time.time() - t0
        print(f"- objfun: {f:.6f} ({dt:.3f} sec) ", end='')
        xparamset.append(x)
        fset.append(f)
        if f < f_best:
            print("**best**", end='')
            f_best = f
        print()

    fset = np.array(fset)
    ind = np.argsort(fset)

    resSCV = []
    resG2 = []
    fobjAC = []
    for i in range(min(MaxResAC, len(ind))):
        SCVj, G2j = xtopar(xparamset[ind[i]])
        G2j = np.clip(G2j, TOL, 1 - TOL)
        resSCV.append(SCVj)
        resG2.append(G2j)
        fobjAC.append(fset[ind[i]])

    return resSCV, resG2, np.array(fobjAC)


# ============================================================================
# Bicovariance fitting (kpcfit_sub_bcfit)
# ============================================================================

def kpcfit_sub_bcfit(E: np.ndarray, SCVj: np.ndarray, G2j: np.ndarray,
                     BC: np.ndarray, BCLags: np.ndarray,
                     MaxIterBC: int = 10, MaxRunsBC: int = 30
                     ) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Fit bicovariance by optimizing E1 and E3 for each MAP component.

    Given fixed SCV and G2 from autocorrelation fitting, finds E1 and E3
    for each MAP(2) that best reproduce the empirical bicovariances.

    Args:
        E: Moments of the trace
        SCVj: SCV values for each MAP (fixed from AC fitting)
        G2j: G2 values for each MAP (fixed from AC fitting)
        BC: Empirical bicovariance values
        BCLags: Lag combinations for BC
        MaxIterBC: Maximum iterations per optimization run
        MaxRunsBC: Number of random restarts

    Returns:
        Tuple of (E1j, E3j, f_best) optimal first and third moments
    """
    from ..mam import map_joint, map_scale, map_normalize

    MaxTimeBC = 100  # seconds
    NumMAPs = len(SCVj)
    TOL = 1e-9
    EPSTOL = 10 * TOL

    SCVj = np.asarray(SCVj).ravel()
    G2j = np.asarray(G2j).ravel()
    BC = np.asarray(BC).ravel()
    BCLags = np.asarray(BCLags)

    # Normalizing constant
    NBC = np.linalg.norm(BC, 2)
    if NBC < 1e-30:
        NBC = 1.0

    # Truncate bicovariance entries with decreasing lags
    keep = []
    for idx in range(BCLags.shape[0]):
        if np.any(np.diff(BCLags[idx, :]) < 0):
            continue
        keep.append(idx)
    keep = np.array(keep)
    if len(keep) < len(BC):
        BCLags = BCLags[keep]
        BC = BC[keep]

    # Initial point
    E1j_init = E[0] ** (1.0 / NumMAPs) * np.ones(NumMAPs)
    E2j = np.zeros(NumMAPs)
    E3j_init = np.zeros(NumMAPs)
    t = E[0] * np.random.rand()
    for idx in range(NumMAPs):
        E2j[idx] = (1 + SCVj[idx]) * E1j_init[idx] ** 2
        E3j_init[idx] = (3.0 / 2.0 + t) * E2j[idx] ** 2 / E1j_init[idx]

    x0base = np.concatenate([E1j_init, E3j_init])

    fold = [0.0]  # mutable container for last-good objective

    def xtopar(x):
        E1j = x[:NumMAPs].copy()
        E3j = x[NumMAPs:].copy()
        # Force E1j(1) = E(1) / prod(E1j(2:end))
        prod_rest = np.prod(E1j[1:])
        if prod_rest > 0:
            E1j[0] = E[0] / prod_rest
        return E1j, E3j

    def nnlcon(x):
        E1j, E3j = xtopar(x)
        c = []
        for j in range(1, NumMAPs):
            E2j_loc = (1 + SCVj[j]) * E1j[j] ** 2
            # E2j > 2*E1j^2
            c.append((2 + EPSTOL) * E1j[j] ** 2 - E2j_loc)
            # E3j > (3/2)*E2j^2/E1j
            c.append((3.0 / 2.0 + EPSTOL) * E2j_loc ** 2 / E1j[j] - E3j[j])
        if SCVj[0] > 1:
            E2j_0 = (1 + SCVj[0]) * E1j[0] ** 2
            c.append((2 + EPSTOL) * E1j[0] ** 2 - E2j_0)
            c.append((3.0 / 2.0 + EPSTOL) * E2j_0 ** 2 / E1j[0] - E3j[0])
        # E3 product constraint
        temp = np.prod(E3j)
        temp = temp / (factorial(3)) ** (NumMAPs - 1)
        if len(E) > 2:
            c.append(temp / E[2] - 2)
            c.append(0.5 - temp / E[2])
        # MAP(2) feasibility: G2 bounds from E3
        for j in range(1, NumMAPs):
            E2j_loc = (1 + SCVj[j]) * E1j[j] ** 2
            denom = E2j_loc - 2 * E1j[j] ** 2
            if abs(denom) < 1e-30:
                continue
            val = (1.0 / 3.0) * E1j[j] / denom * E3j[j] - 0.5 * E2j_loc ** 2 / denom
            c.append(val - 1.0 / 1e16)
            c.append(1e-16 - val)
        return np.array(c)

    def objfun(x):
        E1j, E3j = xtopar(x)
        compMAP, subMaps, err = kpcfit_sub_compose(E1j, SCVj, E3j, G2j)
        if err != 0 or compMAP is None:
            f = max(2 * fold[0], 1e6)
            return f
        # Scale MAP to match E(1)
        compMAP = _map_scale_to_mean(compMAP[0], compMAP[1], E[0])
        BCj = np.ones(BCLags.shape[0])
        for indexL in range(BCLags.shape[0]):
            try:
                BCj[indexL] = map_joint(compMAP[0], compMAP[1],
                                        BCLags[indexL, :], np.array([1, 1, 1]))
            except Exception:
                BCj[indexL] = 0.0
        f = np.linalg.norm(BC - BCj, 2) / NBC
        if np.isnan(f):
            f = 2 * fold[0]
        else:
            fold[0] = f
        return f

    # Multi-start optimization
    xparamset = []
    fset = []
    f_best = np.inf

    lb = np.full(2 * NumMAPs, EPSTOL)
    bounds = [(EPSTOL, None)] * (2 * NumMAPs)
    constraints = {'type': 'ineq', 'fun': lambda x: -nnlcon(x)}

    x0 = x0base.copy()
    for ind in range(MaxRunsBC):
        tstart = time.time()

        try:
            res = minimize(objfun, x0, method='SLSQP',
                           bounds=bounds,
                           constraints=constraints,
                           options={'maxiter': MaxIterBC, 'ftol': TOL,
                                    'disp': False})
            x = res.x
            f = res.fun
        except Exception:
            x = x0.copy()
            f = np.inf

        xparamset.append(x)
        fset.append(f)
        if f < f_best:
            f_best = f
            xparam_best = x.copy()

        # Perturb x0 for next run
        x0 = x0base.copy()
        x0[:NumMAPs] *= (0.25 + 1.75 * np.random.rand(NumMAPs))
        t = np.random.rand() * E[0]
        for idx in range(NumMAPs):
            E2j_loc = (1 + SCVj[idx]) * x0[idx] ** 2
            x0[NumMAPs + idx] = (3.0 / 2.0 + t) * E2j_loc ** 2 / x0[idx]

    E1j, E3j = xtopar(xparam_best)
    prod_rest = np.prod(E1j[1:])
    if prod_rest > 0:
        E1j[0] = E[0] / prod_rest
    return E1j, E3j, f_best


# ============================================================================
# BIC order selection (kpcfit_sub_bic)
# ============================================================================

def kpcfit_sub_bic(SA: np.ndarray, orders: ArrayLike) -> int:
    """
    Select MAP order using BIC (Schwarz Bayesian Criterion).

    Performs AR model fitting at different orders and selects the order
    that minimizes the SBC criterion.

    Args:
        SA: Autocorrelation sequence (ACFull from kpcfit_init)
        orders: Candidate orders (e.g., [2, 4, 8, 16, 32, 64, 128])

    Returns:
        Recommended number of MAPs (log2 of states)
    """
    orders = np.asarray(orders, dtype=int)
    nlags = len(SA)

    # nlagsend = point of first decay below threshold
    nlagsend_idx = np.where(SA < 1e-6)[0]
    ordermax = int(max(orders))

    if len(nlagsend_idx) == 0:
        nlagsend = nlags
    else:
        nlagsend = int(nlagsend_idx[0]) - 1 + ordermax

    NLAGSMAX = 10000

    if nlagsend > NLAGSMAX:
        SAlags_base = logspacei(1, nlagsend - ordermax, NLAGSMAX)
        # MATLAB adds linear fill for initial range
        SAlags = SAlags_base
    else:
        if nlagsend > ordermax:
            SAlags = np.arange(1, nlagsend - ordermax + 1)
        else:
            ordermax = max(1, nlagsend - 2)
            SAlags = np.arange(1, max(1, nlagsend - ordermax) + 1)
            orders = orders[orders <= ordermax]

    if len(SAlags) == 0 or len(orders) == 0:
        return 3

    n_samples = len(SAlags)
    SAlags = np.round(SAlags).astype(int)
    SAlags_y = SAlags.copy()

    # Ensure indices are within bounds (1-based MATLAB -> 0-based Python)
    max_idx = len(SA) - 1
    SAlags_y = SAlags_y[SAlags_y <= max_idx]
    if len(SAlags_y) == 0:
        return 3

    Y = SA[SAlags_y - 1]
    n_samples = len(Y)

    SBC = []
    for order_val in orders:
        order_val = int(order_val)

        # Check we have enough data for this order
        max_lag_x = SAlags_y + order_val
        if np.any(max_lag_x > max_idx + 1):
            valid = (SAlags_y + order_val) <= max_idx + 1
            if not np.any(valid):
                SBC.append(np.inf)
                continue
            valid_lags = SAlags_y[valid]
        else:
            valid_lags = SAlags_y

        # Build X matrix (AR design matrix)
        X = np.zeros((len(valid_lags), order_val))
        for i in range(order_val):
            lags_x = valid_lags + (i + 1)
            X[:, i] = SA[lags_x - 1]

        Y_fit = SA[valid_lags - 1]

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                b, residuals, rank, s = np.linalg.lstsq(X, Y_fit, rcond=None)
            r = Y_fit - X @ b
            sse = np.sum(r ** 2)
            n = len(Y_fit)
            sbc = n * np.log(max(sse, 1e-300)) - n * np.log(n) + np.log(n) * order_val
            SBC.append(sbc)
        except Exception:
            SBC.append(np.inf)

    if len(SBC) == 0 or all(np.isinf(SBC)):
        return 3

    best_idx = int(np.argmin(SBC))
    Nstar = int(np.log2(orders[best_idx]))

    # Optional warning for poor fit
    order = 2 ** Nstar
    valid_lags = SAlags_y[(SAlags_y + order) <= max_idx + 1]
    if len(valid_lags) > 0:
        X = np.zeros((len(valid_lags), order))
        for i in range(order):
            X[:, i] = SA[valid_lags + i]
        Y_check = SA[valid_lags - 1]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                b = np.linalg.lstsq(X, Y_check, rcond=None)[0]
            r = Y_check - X @ b
            rssrho = np.sum(r ** 2) / max(np.sum(Y_check ** 2), 1e-30)
            if rssrho > 0.05:
                print("Warning! The residual sum of square indicates that the "
                      "trace may not be well characterized by a MAP with less "
                      "than 64 states!")
        except Exception:
            pass

    return Nstar


# ============================================================================
# MAP composition (kpcfit_sub_compose)
# ============================================================================

def kpcfit_sub_compose(E1j: np.ndarray, SCVj: np.ndarray,
                        E3j: np.ndarray, G2j: np.ndarray,
                        verbose: int = 0) -> Tuple[Optional[Tuple[np.ndarray, np.ndarray]],
                                                    List[Tuple[np.ndarray, np.ndarray]], int]:
    """
    Compose a large MAP from smaller MAP(2)s using Kronecker products.

    The first MAP is an MMPP(2); remaining MAPs are feasible MAP(2) blocks.

    Args:
        E1j: First moments for each MAP
        SCVj: SCVs for each MAP
        E3j: Third moments for each MAP
        G2j: Autocorrelation decay rates for each MAP
        verbose: Verbosity level

    Returns:
        Tuple of (composed_MAP, list_of_sub_MAPs, error_code)
    """
    from ..mam import (map_kpc, map_isfeasible, map_normalize, map_erlang,
                       map_exponential, map_feasblock, map_embedded, map2_fit)
    from ..kpctoolbox import mmpp2_fit3

    J = len(G2j)
    E1j = np.asarray(E1j).ravel()
    SCVj = np.asarray(SCVj).ravel()
    E3j = np.asarray(E3j).ravel()
    G2j = np.asarray(G2j).ravel()

    sub_maps = []

    # First MAP is an MMPP2
    try:
        E2_1 = (1 + SCVj[0]) * E1j[0] ** 2
        kpc_map = mmpp2_fit3(E1j[0], E2_1, E3j[0], G2j[0])
    except Exception:
        kpc_map = None

    kpc_feasible = False
    if kpc_map is not None:
        try:
            # Check for large imaginary parts
            if (np.any(np.abs(np.imag(kpc_map[0])) > 1e-4) or
                    np.any(np.abs(np.imag(kpc_map[1])) > 1e-4)):
                kpc_feasible = False
            else:
                kpc_map = (np.real(kpc_map[0]), np.real(kpc_map[1]))
                kpc_feasible = map_isfeasible(kpc_map[0], kpc_map[1])
        except Exception:
            kpc_feasible = False

    if not kpc_feasible:
        if SCVj[0] < 0.5:
            kpc_map = map_erlang(E1j[0], 2)
            if verbose > 0:
                print("MAP 1 is erlang-2")
        else:
            if verbose > 0:
                print(f"MAP 1 has presumably infeasible E3: {E3j[0]:.6f}")
            E2_1 = (1 + SCVj[0]) * E1j[0] ** 2
            kpc_map, fit_err = map2_fit(E1j[0], E2_1, -1, G2j[0])
            if kpc_map is None or fit_err != 0:
                if verbose > 0:
                    print(f"MAP 1 has infeasible G2: {G2j[0]:.6f} (err={fit_err})")
                kpc_map, fit_err = map2_fit(E1j[0], E2_1, -1, 0)
                if fit_err != 0 or kpc_map is None:
                    return None, [], 1

    sub_maps.append(kpc_map)

    for j in range(1, J):
        E2_j = (1 + SCVj[j]) * E1j[j] ** 2
        map_j = map_feasblock(E1j[j], E2_j, E3j[j], G2j[j])

        try:
            feasible = map_isfeasible(map_j[0], map_j[1])
        except Exception:
            if verbose > 0:
                print(f"Could not check feasibility of MAP-{j+1}")
            feasible = False

        if not feasible:
            if SCVj[j] < 1:
                if verbose > 0:
                    print(f"MAP {j+1} has low variability")
                map_j = map_exponential(E1j[j])
            else:
                if verbose > 0:
                    print(f"MAP {j+1} has presumably infeasible E3")
                map_j, fit_err = map2_fit(E1j[j], E2_j, -1, G2j[j])
                if fit_err != 0 or map_j is None:
                    map_j, fit_err = map2_fit(E1j[j], E2_j, -1, 0)
                    if verbose > 0:
                        print(f"MAP {j+1} has infeasible G2: {G2j[j]:.6f}")
                    if fit_err != 0 or map_j is None:
                        return None, sub_maps, 5

                # MATLAB: eigenvalue normalization for high-SCV blocks
                if map_j is not None:
                    try:
                        lam = np.linalg.eigvals(-np.linalg.inv(map_j[0]))
                        D0 = np.diag(-1.0 / lam)
                        P = map_embedded(map_j[0], map_j[1])
                        p = min(np.real(np.linalg.eigvals(P)))
                        p = max(0, min(1, p))
                        D1 = -D0 @ np.array([[p, 1 - p], [p, 1 - p]])
                        map_j = map_normalize(np.real(D0), np.real(D1))
                    except Exception:
                        pass

        if map_j is None:
            if verbose > 0:
                print(f"Replacing MAP {j+1} with an exponential")
            map_j = map_exponential(E1j[j])

        sub_maps.append(map_j)
        kpc_map = map_kpc([kpc_map, map_j])

    # Check feasibility of all sub-MAPs
    for j, smap in enumerate(sub_maps):
        if not map_isfeasible(smap[0], smap[1]):
            if verbose > 0:
                print(f"MAP {j+1} is infeasible")
            return None, sub_maps, 10

    kpc_map = map_normalize(kpc_map[0], kpc_map[1])
    return kpc_map, sub_maps, 0


# ============================================================================
# Manual KPC fitting (kpcfit_manual)
# ============================================================================

def kpcfit_manual(NumMAPs: int, E: np.ndarray, AC: np.ndarray,
                  ACLags: np.ndarray, BC: np.ndarray, BCLags: np.ndarray,
                  MaxIterAC: int = 300, MaxRunsAC: int = 50,
                  MaxResAC: int = 10, MaxRunsBC: int = 30,
                  MaxIterBC: int = 10, OnlyAC: bool = False,
                  AnimateAC: bool = False, MaxRetMAPs: int = 1,
                  ParallelBC: int = 0
                  ) -> Tuple[Optional[Tuple[np.ndarray, np.ndarray]],
                             float, float,
                             List[Tuple[np.ndarray, np.ndarray]],
                             List[Tuple[np.ndarray, np.ndarray]],
                             List[float], List[float],
                             List[List[Tuple[np.ndarray, np.ndarray]]]]:
    """
    Manual KPC fitting with user-specified parameters.

    Performs two-phase fitting:
    1) Autocorrelation fitting to determine SCV and G2 for each MAP component
    2) Bicovariance fitting to determine E1 and E3 for each MAP component

    Args:
        NumMAPs: Number of MAP(2)s to compose (log2 of total states)
        E: Moments of the trace
        AC: Autocorrelation values at ACLags
        ACLags: Autocorrelation lags
        BC: Bicovariance values
        BCLags: Bicovariance lag combinations
        MaxIterAC: Max iterations per AC fitting run
        MaxRunsAC: Number of AC fitting restarts
        MaxResAC: Number of AC results to keep for BC fitting
        MaxRunsBC: Number of BC fitting restarts
        MaxIterBC: Max iterations per BC fitting run
        OnlyAC: If True, skip bicovariance fitting
        AnimateAC: (unused in Python)
        MaxRetMAPs: Maximum number of MAPs returned
        ParallelBC: (unused in Python, parallelization flag)

    Returns:
        Tuple of (bestMAP, fac, fbc, subMAPs,
                  otherMAPs, otherFACs, otherFBCs, otherSubMAPs)
    """
    from ..mam import (map_acf, map_scv, map_skew, map_joint,
                       map_moment)

    E = np.asarray(E).ravel()
    AC = np.asarray(AC).ravel()
    ACLags = np.asarray(ACLags).ravel()
    BC = np.asarray(BC).ravel()
    BCLags = np.asarray(BCLags)

    MaxResAC = min(MaxResAC, MaxRunsAC)

    # --- Phase 1: fit autocorrelations ---
    print()
    resSCV, resG2, fobjAC = kpcfit_sub_acfit(
        E, AC, ACLags, NumMAPs, MaxIterAC, MaxRunsAC, MaxResAC, AnimateAC)

    # --- Phase 2: fit bicovariances ---
    print()
    resE1 = [None] * len(resSCV)
    resE3 = [None] * len(resSCV)
    fobjBC = np.zeros(len(resSCV))

    if not OnlyAC:
        for i in range(len(resSCV)):
            print(f"fitting: processing acfit subresult {i+1} of top {len(resSCV)}")
            E1j, E3j, foBC = kpcfit_sub_bcfit(
                E, resSCV[i], resG2[i], BC, BCLags, MaxIterBC, MaxRunsBC)
            resE1[i] = E1j
            resE3[i] = E3j
            fobjBC[i] = foBC
    else:
        for i in range(len(resSCV)):
            resE1[i] = np.ones(NumMAPs)
            resE3[i] = (3.0 / 2.0 + 0.01) * (1 + resSCV[i]) ** 2
            fobjBC[i] = -1

    # --- Sort by bicovariance objective ---
    ind = np.argsort(fobjBC)

    # Truncate BCLags with decreasing lags (same as MATLAB)
    keep = []
    for index in range(BCLags.shape[0]):
        if np.any(np.diff(BCLags[index, :]) < 0):
            continue
        keep.append(index)
    if len(keep) < len(BC):
        keep = np.array(keep)
        BCLags = BCLags[keep]
        BC = BC[keep]

    # --- Helper: evaluate objective functions for a composed MAP ---
    def evaluate_obj_function(map_tuple):
        D0, D1 = map_tuple
        tSCV = (E[1] - E[0] ** 2) / E[0] ** 2
        acf_vals = map_acf(D0, D1, ACLags)
        objAC = (np.linalg.norm(AC - acf_vals, 1) / np.linalg.norm(AC, 2) +
                 (map_scv(D0, D1) - tSCV) ** 2 / tSCV ** 2)

        mapBC = np.ones(BCLags.shape[0])
        for j_idx in range(BCLags.shape[0]):
            try:
                mapBC[j_idx] = map_joint(D0, D1, BCLags[j_idx, :],
                                          np.array([1, 1, 1]))
            except Exception:
                mapBC[j_idx] = 0.0
        objBC = np.linalg.norm(BC - mapBC, 2) / max(np.linalg.norm(BC, 2), 1e-30)
        return objAC, objBC

    # --- Compose intermediate results into final MAPs ---
    composedMAPs = 0
    MAPs = []
    subs = []
    FACs = []
    FBCs = []

    for k in range(len(ind)):
        index = ind[k]
        COMPOSE_VERBOSE = 1
        newMAP, newSubMAPs, errorCode = kpcfit_sub_compose(
            resE1[index], resSCV[index], resE3[index], resG2[index],
            COMPOSE_VERBOSE)
        if errorCode != 0 or newMAP is None:
            print(f"[{k+1}] Discarded (errorCode={errorCode}).")
            continue

        # Scale MAP to match mean exactly
        newMAP = _map_scale_to_mean(newMAP[0], newMAP[1], E[0])

        # Compute objective functions
        newfAC, newfBC = evaluate_obj_function(newMAP)

        composedMAPs += 1
        MAPs.append(newMAP)
        FACs.append(newfAC)
        FBCs.append(newfBC)
        subs.append(newSubMAPs)

        if COMPOSE_VERBOSE > 0:
            print(f"[{k+1}] OK - fac= {newfAC:.6f}; fbc= {newfBC:.6f}")

        if composedMAPs == MaxRetMAPs:
            break

    if composedMAPs == 0:
        print("+++ KPC FAILED +++")
        return None, np.inf, np.inf, [], [], [], [], []

    bestMAP = MAPs[0]
    subMAPs = subs[0]
    fbc = FBCs[0]
    fac = FACs[0]

    D0, D1 = bestMAP
    print(f"Returned {composedMAPs} MAPs:")
    print(f"1) fAC={fac:.6f}, fBC={fbc:.6f}, SCV={map_scv(D0, D1):.6f}, "
          f"ACF(1)={float(map_acf(D0, D1, 1)):.6f}, SKEW={map_skew(D0, D1):.6f}")

    otherMAPs = []
    otherFACs = []
    otherFBCs = []
    otherSubMAPs = []

    for i in range(1, composedMAPs):
        otherMAPs.append(MAPs[i])
        otherFACs.append(FACs[i])
        otherFBCs.append(FBCs[i])
        otherSubMAPs.append(subs[i])
        D0o, D1o = MAPs[i]
        print(f"{i+1}) fAC={FACs[i]:.6f}, fBC={FBCs[i]:.6f}, "
              f"SCV={map_scv(D0o, D1o):.6f}, "
              f"ACF(1)={float(map_acf(D0o, D1o, 1)):.6f}, "
              f"SKEW={map_skew(D0o, D1o):.6f}")

    print()
    return bestMAP, fac, fbc, subMAPs, otherMAPs, otherFACs, otherFBCs, otherSubMAPs


# ============================================================================
# Automatic KPC fitting (kpcfit_auto)
# ============================================================================

def kpcfit_auto(trace_data: KpcfitTraceData,
                OnlyAC: bool = False,
                NumStates: int = None,
                NumMAPs: int = None,
                MaxIterAC: int = 300,
                MaxIterBC: int = 10,
                MaxRunsAC: int = 50,
                AnimateAC: bool = False,
                MaxRunsBC: int = 30,
                MaxResAC: int = None,
                MaxRetMAPs: int = 1,
                ParallelBC: int = 0
                ) -> Tuple[Optional[Tuple[np.ndarray, np.ndarray]],
                           float, float,
                           List[Tuple[np.ndarray, np.ndarray]],
                           List[Tuple[np.ndarray, np.ndarray]],
                           List[float], List[float],
                           List[List[Tuple[np.ndarray, np.ndarray]]]]:
    """
    Automatic MAP fitting using the KPC method.

    Orchestrates the full KPC fitting pipeline:
    1) BIC-based order selection (if NumMAPs not specified)
    2) Autocorrelation fitting
    3) Bicovariance fitting (unless OnlyAC=True)
    4) MAP composition and validation

    Args:
        trace_data: Preprocessed trace data from kpcfit_init
        OnlyAC: If True, fit only moments and autocorrelations
        NumStates: Number of states (must be power of 2, overrides NumMAPs)
        NumMAPs: Number of MAP(2)s = log2(NumStates). Auto-selected if None.
        MaxIterAC: Max iterations per AC fitting run
        MaxIterBC: Max iterations per BC fitting run
        MaxRunsAC: Number of AC fitting restarts
        AnimateAC: (unused in Python)
        MaxRunsBC: Number of BC fitting restarts
        MaxResAC: AC results to pass to BC fitting (default: min(MaxRunsAC, 10))
        MaxRetMAPs: Maximum number of MAPs returned
        ParallelBC: (unused in Python)

    Returns:
        Tuple of (bestMAP, fac, fbc, kpcMAPs,
                  otherMAPs, otherFACs, otherFBCs, otherSubMAPs)
    """
    from ..mam import map_moment, map_acf

    if MaxResAC is None:
        MaxResAC = min(MaxRunsAC, 10)

    if NumStates is not None:
        if NumStates > 0 and (NumStates & (NumStates - 1)) != 0:
            raise ValueError(
                f"Requested {NumStates} states, but kpc-toolbox can return "
                f"only a number of states that is a power of 2")
        if NumMAPs is None:
            NumMAPs = int(np.ceil(np.log2(NumStates)))

    print("** KPC fitting algorithm initialized **")

    # Order selection
    if NumMAPs is None:
        print("init: performing order selection")
        try:
            NumMAPs = kpcfit_sub_bic(trace_data.ACFull,
                                     np.array([2, 4, 8, 16, 32, 64, 128]))
        except Exception:
            print("BIC order selection failed. Using 3 MAPs")
            NumMAPs = 3

    print(f"init: kpc-toolbox will search for a MAP with "
          f"2^{NumMAPs}={2**NumMAPs} states")

    # Fitting
    print("fitting: running KPC-based fitting script")
    result = kpcfit_manual(
        NumMAPs, trace_data.E, trace_data.AC, trace_data.ACLags,
        trace_data.BC, trace_data.BCLags,
        MaxIterAC=MaxIterAC, MaxRunsAC=MaxRunsAC,
        MaxResAC=MaxResAC, MaxRunsBC=MaxRunsBC,
        MaxIterBC=MaxIterBC, OnlyAC=OnlyAC,
        AnimateAC=AnimateAC, MaxRetMAPs=MaxRetMAPs,
        ParallelBC=ParallelBC)

    bestMAP, fac, fbc, kpcMAPs = result[0], result[1], result[2], result[3]
    otherMAPs, otherFACs, otherFBCs, otherSubMAPs = result[4], result[5], result[6], result[7]

    # Display final moments and autocorrelations
    print("** KPC fitting algorithm completed ** ")
    print()
    if bestMAP is not None:
        D0, D1 = bestMAP
        print("                  Moments Comparison       ")
        print("          Original Trace          Fitted MAP")
        for k in range(len(trace_data.E)):
            mk = map_moment(D0, D1, k + 1)
            print(f"   {trace_data.E[k]:.10e}   {mk:.10e}")
        print()
        print("                             Autocorrelation Comparison       ")
        print("          Lag                     Original Trace          Fitted MAP")
        n_show = min(10, len(trace_data.ACLags))
        for k in range(n_show):
            acf_k = float(map_acf(D0, D1, int(trace_data.ACLags[k])))
            print(f"   {trace_data.ACLags[k]:10d}   {trace_data.AC[k]:.10e}   {acf_k:.10e}")

    return bestMAP, fac, fbc, kpcMAPs, otherMAPs, otherFACs, otherFBCs, otherSubMAPs


# ============================================================================
# Characteristic polynomial and Prony's method
# ============================================================================

def kpcfit_hyper_charpoly(E: np.ndarray, n: int) -> np.ndarray:
    """
    Compute characteristic polynomial coefficients for hyperexponential fitting.

    Solves a system of equations relating moments to polynomial coefficients
    used in Prony's method for hyperexponential fitting.

    Args:
        E: Moments E[X^k]
        n: Number of phases

    Returns:
        Polynomial coefficients
    """
    E_full = np.concatenate([[1], np.asarray(E).ravel()])
    f = np.array([factorial(i) for i in range(2 * n)])

    A = np.zeros((n + 1, n + 1))
    for i in range(n):
        for j in range(n + 1):
            idx = n + i - j
            if 0 <= idx < len(E_full) and 0 <= idx < len(f):
                A[i, j] = E_full[idx] / f[idx]
    A[n, 0] = 1

    b = np.zeros(n + 1)
    b[n] = 1

    try:
        m = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        m = np.linalg.lstsq(A, b, rcond=None)[0]

    return m


def kpcfit_ph_prony(E: np.ndarray, n: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit a hyperexponential PH distribution using Prony's method.

    Args:
        E: Moments E[X^k]
        n: Number of phases

    Returns:
        Tuple (D0, D1) representing the PH distribution as a MAP
    """
    E = np.asarray(E).ravel()
    f = np.array([factorial(i) for i in range(2 * n)])

    m = kpcfit_hyper_charpoly(E, n)
    theta = np.roots(m[::-1])  # Roots of characteristic polynomial

    # Filter real positive roots
    theta = np.real(theta[np.abs(np.imag(theta)) < 1e-10])
    theta = theta[theta > 0]

    if len(theta) < n:
        # Pad with positive values
        theta = np.concatenate([theta, np.abs(theta[-1]) * np.ones(n - len(theta))])
    theta = theta[:n]

    # Compute entry probabilities
    C = np.zeros((n, n))
    for i in range(n):
        C[i, :] = f[i + 1] * (theta ** (i + 1))

    try:
        M = np.linalg.solve(C, E[:n])
    except np.linalg.LinAlgError:
        M = np.linalg.lstsq(C, E[:n], rcond=None)[0]

    D0 = np.diag(-1 / theta)
    D1 = -D0 @ (np.ones((n, 1)) @ M.reshape(1, -1))

    return D0, D1


# ============================================================================
# PH fitting options
# ============================================================================

def kpcfit_ph_options(E: np.ndarray, **kwargs) -> KpcfitPhOptions:
    """
    Create options for PH fitting.

    Args:
        E: Moments vector
        **kwargs: Option overrides matching MATLAB names:
            - Verbose / verbose: Print progress (default: True)
            - Runs / runs: Number of optimization runs (default: 5)
            - MinNumStates / min_num_states: Minimum states (default: 2)
            - MaxNumStates / max_num_states: Maximum states (default: 32)
            - MinExactMom / min_exact_mom: Minimum moments to fit exactly (default: 3)

    Returns:
        KpcfitPhOptions dataclass
    """
    import math

    options = KpcfitPhOptions()
    E = np.asarray(E)

    # Map MATLAB-style option names to Python field names
    name_map = {
        'verbose': 'verbose',
        'runs': 'runs',
        'minnumstates': 'min_num_states',
        'maxnumstates': 'max_num_states',
        'minexactmom': 'min_exact_mom',
        'min_num_states': 'min_num_states',
        'max_num_states': 'max_num_states',
        'min_exact_mom': 'min_exact_mom',
    }

    for key, value in kwargs.items():
        key_lower = key.lower().replace('_', '')
        # Direct field match
        if hasattr(options, key):
            setattr(options, key, value)
        elif key_lower in name_map:
            setattr(options, name_map[key_lower], value)

    # Ensure states are powers of 2
    if options.min_num_states > 0:
        if options.min_num_states != 2 ** round(math.log2(options.min_num_states)):
            options.min_num_states = 2 ** math.ceil(math.log2(options.min_num_states))
    if options.max_num_states > 0:
        if options.max_num_states != 2 ** round(math.log2(options.max_num_states)):
            options.max_num_states = 2 ** math.ceil(math.log2(options.max_num_states))

    # Swap if min > max
    if options.min_num_states > options.max_num_states:
        options.min_num_states, options.max_num_states = (
            options.max_num_states, options.min_num_states)

    # Check moment requirements
    if 2 * options.max_num_states - 1 > len(E):
        raise ValueError(f"MaxNumStates of {options.max_num_states} requires "
                        f"at least {2 * options.max_num_states - 1} moments.")

    return options


# ============================================================================
# PH exact fitting
# ============================================================================

def kpcfit_ph_exact(E: np.ndarray, options: KpcfitPhOptions
                    ) -> List[Tuple[np.ndarray, np.ndarray]]:
    """
    Exact PH fitting methods (Prony's method for hyperexp, APH for low variability).

    Args:
        E: Moments E[X^k] (typically normalized so E[X]=1)
        options: Fitting options

    Returns:
        List of fitted PH distributions (D0, D1)
    """
    from ..mam import map_isfeasible, map_erlang, map_moment, map2_fit
    from ..kpctoolbox import aph_fit

    E = np.asarray(E).ravel()
    PH_EXACT = []

    SCV = (E[1] - E[0] ** 2) / E[0] ** 2

    if SCV > 1:
        # Higher variability - use hyperexponential fitting
        if options.verbose:
            print(f"kpcfit_ph: HIGHER variability than an exponential "
                  f"(var/mean^2 = {SCV:.6f})")
            print("kpcfit_ph: starting exact hyper-exponential fitting method "
                  "(Prony's method)")

        for n in range(2, options.max_num_states + 1):
            if len(E) < 2 * n - 1:
                if options.verbose:
                    print(f"kpcfit_ph: not enough moments given in input to "
                          f"fit hyper-exp({n})")
                break

            PH = kpcfit_ph_prony(E, n)
            if map_isfeasible(PH[0], PH[1]):
                PH_EXACT.append(PH)
                if options.verbose:
                    print(f"\t\t\thyper-exp({n}): feasible, matched exactly "
                          f"{2*n-1} moments. result saved.")
            else:
                if options.verbose:
                    print(f"\t\t\thyper-exp({n}): infeasible to fit exactly.")
                break

    elif SCV < 1:
        # Lower variability
        if options.verbose:
            print(f"kpcfit_ph: LOWER variability than an exponential "
                  f"(var/mean^2 = {SCV:.6f})")

        n = 1
        while 1 / n > SCV:
            n += 1

        if options.verbose:
            print(f"kpcfit_ph: exact fitting of E[X^2] requires at least "
                  f"{n} states")

        if options.min_exact_mom >= 2 and options.max_num_states < n:
            if options.verbose:
                print(f"kpcfit_ph: impossible to fit exactly E[X^2] with "
                      f"MaxNumStates = {options.max_num_states}, "
                      f"increasing to MaxNumStates = {n}.")
            return PH_EXACT

        if n == 2:
            if options.verbose:
                print("kpcfit_ph: attempting PH(2) fitting method")
            PH, err = map2_fit(E[0], E[1], E[2] if len(E) > 2 else -1, 0)
            if PH is None:
                PH, err = map2_fit(E[0], E[1], -1, 0)
                if PH is not None and map_isfeasible(PH[0], PH[1]):
                    PH_EXACT.append(PH)
                    if options.verbose:
                        print(f"\t\t\tph(2): feasible, matched exactly 2 "
                              f"moments. result saved.")
            else:
                PH_EXACT.append(PH)
                if options.verbose:
                    print(f"\t\t\tph(2): feasible, matched exactly 3 "
                          f"moments. result saved.")
        elif abs(SCV - 1 / n) < KPCFIT_TOL:
            # Erlang case
            ERL = map_erlang(E[0], n)
            ERL_moments = np.array([map_moment(ERL[0], ERL[1], k + 1)
                                    for k in range(len(E))])
            if np.linalg.norm(E - ERL_moments) < KPCFIT_TOL:
                if options.verbose:
                    print(f"kpcfit_ph: erlang moment set. fitted erlang-{n}. "
                          f"result saved.")
                PH_EXACT.append(ERL)
        else:
            # APH fitting
            maxorder = options.max_num_states
            if options.verbose:
                print(f"kpcfit_ph: fitting APH distribution (best effort, "
                      f"max order = {maxorder}).")
            try:
                PH = aph_fit(E[0], E[1], E[2] if len(E) > 2 else E[1] * 1.5,
                             maxorder)
                if map_isfeasible(PH[0], PH[1]):
                    aph_matched = 0
                    n_ph = PH[0].shape[0]
                    for k in range(2 * n_ph - 1):
                        if k < len(E):
                            mk = map_moment(PH[0], PH[1], k + 1)
                            if abs(E[k] - mk) < KPCFIT_TOL * mk:
                                aph_matched += 1
                    if options.verbose:
                        print(f"\t\t\t      aph({n_ph}): feasible, matched "
                              f"exactly {aph_matched} moments. result saved.")
                    PH_EXACT.append(PH)
                else:
                    if options.verbose:
                        print("kpcfit_ph: cannot fit APH distribution.")
            except Exception:
                if options.verbose:
                    print("kpcfit_ph: cannot fit APH distribution.")
    else:
        # SCV == 1 (exponential)
        if options.verbose:
            print(f"kpcfit_ph: SAME variability than an exponential "
                  f"(var/mean^2 = {SCV:.6f})")
        EXP = (np.array([[-1 / E[0]]]), np.array([[1 / E[0]]]))
        EXP_moments = np.array([map_moment(EXP[0], EXP[1], k + 1)
                                for k in range(len(E))])
        if np.linalg.norm(E - EXP_moments) < KPCFIT_TOL:
            if options.verbose:
                print("kpcfit_ph: exponential moment set. fitted exponential. "
                      "result saved.")
        PH_EXACT.append(EXP)

    return PH_EXACT


# ============================================================================
# PH search (optimization-based approximate PH fitting)
# ============================================================================

def kpcfit_ph_search(E: np.ndarray, J: int, options: KpcfitPhOptions,
                     x0: np.ndarray = None, max_aph_order: int = None
                     ) -> Tuple[Tuple[np.ndarray, np.ndarray], float, np.ndarray]:
    """
    Optimization-based search for PH distributions via KPC composition.

    Searches the moment space of J PH(2) distributions to find a composition
    whose moments match the target moments E. Uses fmincon-equivalent
    constrained optimization.

    Args:
        E: Target moments (typically normalized to E[X]=1)
        J: Number of PH(2)s to compose
        options: Fitting options
        x0: Initial point (None for automatic)
        max_aph_order: Maximum APH order (default: ceil(1/SCV))

    Returns:
        Tuple of (PH_distribution, score, x_solution)
    """
    from ..mam import (map_isfeasible, map_moment, map_feasblock,
                       map_kpc, map_normalize)
    from ..kpctoolbox import aph_fit

    E = np.asarray(E).ravel()
    SCV = (E[1] - E[0] ** 2) / E[0] ** 2

    if max_aph_order is None:
        max_aph_order = max(2, int(np.ceil(1.0 / max(SCV, 1e-10))))

    order = [2] * J  # Current version assumes all PH(2)s
    K = max(2 * max(order) - 1, len(E))
    F = np.array([factorial(k) for k in range(1, K + 1)])

    TOL = KPCFIT_TOL

    # Weight function
    w = np.power(np.log(E[-1]) ** (1.0 / len(E)), -np.arange(1, len(E) + 1))

    # Initialize PH(2) components
    PH = [None] * J
    logEtable = np.zeros((J, K))

    for j in range(J):
        PH[j] = map_feasblock(np.random.rand(), 1000 * np.random.rand(), -1, 0)
        for k in range(K):
            mk = map_moment(PH[j][0], PH[j][1], k + 1)
            logEtable[j, k] = np.log(max(mk, 1e-300))

    if x0 is None or len(x0) == 0:
        x0 = logEtable.ravel()
    else:
        x0 = np.asarray(x0).ravel()
        if len(x0) == J * K:
            pass
        else:
            x0 = logEtable.ravel()

    # Stagnation tracking
    stagnval = [0.0]
    stagniter = [0]
    lgkx = [x0.copy()]

    def _callback(xk):
        fval = objfun(xk)
        if not np.isnan(fval) and np.sum(np.isnan(xk)) == 0:
            lgkx[0] = xk.copy()
        else:
            raise StopIteration
        if stagnval[0] == 0:
            stagnval[0] = fval
        delta = abs(fval - stagnval[0]) / max(abs(stagnval[0]), 1e-30)
        if delta < 0.01:
            stagniter[0] += 1
            if stagniter[0] >= 100:
                raise StopIteration
        else:
            stagniter[0] = 0
        stagnval[0] = fval

    def _compute_logEtable_full(logEtable_in):
        """Compute higher-order moments from first 3 via charpoly."""
        logEtable_out = logEtable_in.copy()
        if SCV < 1:
            # APH for first component
            j = 0
            try:
                ph = aph_fit(np.exp(logEtable_in[j, 0]),
                             np.exp(logEtable_in[j, 1]),
                             np.exp(logEtable_in[j, 2]),
                             max_aph_order)
                ph = _map_scale_to_mean(ph[0], ph[1], np.exp(logEtable_in[j, 0]))
                for k in range(K):
                    mk = map_moment(ph[0], ph[1], k + 1)
                    logEtable_out[j, k] = np.log(max(mk, 1e-300))
            except Exception:
                pass
            for j in range(1, J):
                try:
                    ph_j = map_feasblock(np.exp(logEtable_in[j, 0]),
                                         np.exp(logEtable_in[j, 1]),
                                         np.exp(logEtable_in[j, 2]), 0)
                    ph_j = _map_scale_to_mean(ph_j[0], ph_j[1],
                                              np.exp(logEtable_in[j, 0]))
                    for k in range(K):
                        mk = map_moment(ph_j[0], ph_j[1], k + 1)
                        logEtable_out[j, k] = np.log(max(mk, 1e-300))
                except Exception:
                    pass
        else:
            # Hyperexponential: use charpoly for higher moments
            for j in range(1, J):
                try:
                    moments_j = np.exp(logEtable_in[j, :]) / F
                    m_poly = kpcfit_hyper_charpoly(moments_j, 2)
                    if np.sum(np.isnan(m_poly)) > 0:
                        # Exponential fallback
                        exp_map = (np.array([[-1]]), np.array([[1]]))
                        for k in range(3, K):
                            idxs = np.arange(k - 1, k - (2 * order[j] - 2) - 1, -1)
                            valid = (idxs >= 0) & (idxs < K)
                            if not np.all(valid):
                                continue
                            Ej = np.array([map_moment(exp_map[0], exp_map[1], int(idx_val) + 1)
                                           for idx_val in idxs])
                            val = -F[k] * m_poly[1:] @ (Ej / F[idxs])
                            if val > 0:
                                logEtable_out[j, k] = np.log(val)
                    else:
                        for k in range(3, K):
                            idxs = np.arange(k - 1, k - (2 * order[j] - 2) - 1, -1)
                            valid = (idxs >= 0) & (idxs < K)
                            if not np.all(valid):
                                continue
                            Ej = np.exp(logEtable_in[j, idxs[valid]])
                            if np.any(np.isnan(Ej)):
                                exp_map = (np.array([[-1]]), np.array([[1]]))
                                Ej = np.array([map_moment(exp_map[0], exp_map[1], int(idx_val) + 1)
                                               for idx_val in idxs[valid]])
                            val = -F[k] * m_poly[1:len(idxs[valid]) + 1] @ (Ej / F[idxs[valid]])
                            if val > 0:
                                logEtable_out[j, k] = np.log(val)
                except Exception:
                    pass
        return logEtable_out

    def nnlcon(x):
        logET = x.reshape(J, K)
        logET = _compute_logEtable_full(logET)
        logEcur = np.sum(logET, axis=0) - (J - 1) * np.log(F)
        c = []
        j = 0
        if SCV < 1:
            # APH moment constraints
            c.append(-(logET[j, 1] - np.log(max_aph_order) + np.log(max_aph_order - 1)))
            c.append(-(logET[j, 2] - 2 * logET[j, 1] + logET[j, 0] -
                       np.log(max_aph_order + 1) + np.log(max_aph_order)))
        else:
            # SCV > 1
            c.append(-(logET[j, 1] - 2 * logET[j, 0] - np.log(2)) - 10 * TOL)
            c.append(-(logET[j, 2] - (np.log(3.0 / 2.0) + 2 * logET[j, 1] - logET[j, 0])))

        for j in range(1, J):
            c.append(-(logET[j, 1] - 2 * logET[j, 0] - np.log(2)) - 10 * TOL)
            c.append(-(logET[j, 2] - (np.log(3.0 / 2.0) + 2 * logET[j, 1] - logET[j, 0])))

        c = np.array(c)
        if np.any(np.isnan(c)):
            c = 1e10 * np.ones_like(c)

        # Equality constraints: first MinExactMom moments match
        ceq = w[:options.min_exact_mom] * (
            np.log(E[:options.min_exact_mom]) - logEcur[:options.min_exact_mom])
        if np.any(np.isnan(ceq)):
            ceq = 1e10 * np.ones_like(ceq)

        return c, ceq

    def objfun(x):
        if np.any(np.isnan(x)):
            return 1e10
        logET = x.reshape(J, K)
        if SCV < 1:
            logET = _compute_logEtable_full(logET)
        else:
            logET = _compute_logEtable_full(logET)
        logEcur = np.sum(logET, axis=0) - (J - 1) * np.log(F)
        f = w @ np.abs(np.log(E) - logEcur[:len(E)])
        if np.isnan(f):
            f = 1e10
        return f

    # Build constraints for scipy
    def ineq_fun(x):
        c, ceq = nnlcon(x)
        return -c  # scipy uses >= 0

    def eq_fun(x):
        c, ceq = nnlcon(x)
        return ceq

    bounds = [(-200, 200)] * len(x0)  # Matching MATLAB bounds

    constraints_list = [
        {'type': 'ineq', 'fun': ineq_fun},
        {'type': 'eq', 'fun': eq_fun},
    ]

    try:
        res = minimize(objfun, x0, method='SLSQP',
                       bounds=bounds,
                       constraints=constraints_list,
                       options={'maxiter': 200, 'ftol': TOL, 'disp': False},
                       callback=_callback)
        x = res.x
        score = res.fun
    except StopIteration:
        x = lgkx[0]
        score = objfun(x)

    # Reconstruct PH from solution
    logEtable = x.reshape(J, K)

    j = 0
    if SCV < 1:
        try:
            PH[0] = aph_fit(np.exp(logEtable[0, 0]),
                            np.exp(logEtable[0, 1]),
                            np.exp(logEtable[0, 2]),
                            max_aph_order)
            PH[0] = _map_scale_to_mean(PH[0][0], PH[0][1],
                                        np.exp(logEtable[0, 0]))
        except Exception:
            PH[0] = map_feasblock(np.exp(logEtable[0, 0]),
                                  np.exp(logEtable[0, 1]),
                                  np.exp(logEtable[0, 2]), 0)
    else:
        PH[0] = map_feasblock(np.exp(logEtable[0, 0]),
                              np.exp(logEtable[0, 1]),
                              np.exp(logEtable[0, 2]), 0)

    for j in range(1, J):
        PH[j] = map_feasblock(np.exp(logEtable[j, 0]),
                              np.exp(logEtable[j, 1]),
                              np.exp(logEtable[j, 2]), 0)

    KPC_PH = PH[0]
    for j in range(1, J):
        KPC_PH = map_kpc([KPC_PH, PH[j]])

    return KPC_PH, score, x


# ============================================================================
# Parameter-space fitting (psfit_hyperexp_weighted, psfit_aph_weighted)
# ============================================================================

def psfit_hyperexp_weighted(E: np.ndarray, n: int
                            ) -> Tuple[Tuple[np.ndarray, np.ndarray],
                                       np.ndarray, np.ndarray]:
    """
    Fit hyperexponential distribution by parameter-space optimization.

    Optimizes entry probabilities and rates of an n-state hyperexponential
    to minimize weighted moment distance. First 3 moments are matched
    as equality constraints.

    Args:
        E: Moments E[X^k]
        n: Number of phases

    Returns:
        Tuple of (MAP, Eapx, x) where MAP is (D0, D1),
        Eapx are the fitted moments, x is the solution vector
    """
    from ..mam import map_normalize

    E = np.asarray(E).ravel()

    # Scale moments to E[X]=1
    Eunscaled = E.copy()
    Escale = np.array([E[0] ** k for k in range(1, len(E) + 1)])
    E = E / Escale

    w = np.power(np.log(E[-1]) ** (1.0 / len(E)), -np.arange(1, len(E) + 1))

    # Initial point: alpha (entry probs) and theta (service times)
    x0 = np.random.rand(2 * n)
    x0[:n] = x0[:n] / np.sum(x0[:n])

    TOL = 1e-8
    EPSTOL = 100 * TOL
    MAXITER = 100

    def topar(x):
        alpha = x[:n]
        T = np.diag(-1.0 / x[n:])
        return alpha, T

    def _compute_moments(alpha, T):
        n_local = len(alpha)
        invT = np.linalg.inv(-T)
        e = np.ones(n_local)
        Eapx = np.zeros(len(E))
        invT_k = invT.copy()
        for j in range(len(E)):
            Eapx[j] = factorial(j + 1) * alpha @ invT_k @ e
            invT_k = invT_k @ invT
        return Eapx

    def nnlcon(x):
        alpha, T = topar(x)
        Eapx = _compute_moments(alpha, T)
        c = -alpha  # alpha >= 0
        ceq = np.zeros(3 + n)
        ceq[:3] = w[:3] * np.abs(np.log(E[:3]) - np.log(np.maximum(Eapx[:3], 1e-300)))
        ceq[3:] = np.sum(alpha) - 1  # sum(alpha) = 1
        return c, ceq

    def objfun(x):
        alpha, T = topar(x)
        Eapx = _compute_moments(alpha, T)
        Eapx = np.maximum(Eapx, 1e-300)
        return w @ np.abs(np.log(E) - np.log(Eapx))

    def ineq_fun(x):
        c, ceq = nnlcon(x)
        return -c

    def eq_fun(x):
        c, ceq = nnlcon(x)
        return ceq

    bounds = [(EPSTOL, None)] * (2 * n)

    try:
        res = minimize(objfun, x0, method='SLSQP',
                       bounds=bounds,
                       constraints=[
                           {'type': 'ineq', 'fun': ineq_fun},
                           {'type': 'eq', 'fun': eq_fun},
                       ],
                       options={'maxiter': MAXITER, 'ftol': TOL, 'disp': False})
        x = res.x
    except Exception:
        x = x0

    alpha, T = topar(x)
    D0 = T
    D1 = -T @ (np.ones((n, 1)) @ alpha.reshape(1, -1))
    MAP = _map_scale_to_mean(D0, D1, Eunscaled[0])
    MAP = map_normalize(MAP[0], MAP[1])
    Eapx = _compute_moments(alpha, T) * Escale
    return MAP, Eapx, x


def psfit_aph_weighted(E: np.ndarray, n: int
                       ) -> Tuple[Tuple[np.ndarray, np.ndarray],
                                  np.ndarray, np.ndarray]:
    """
    Fit APH (acyclic phase-type) distribution by parameter-space optimization.

    Optimizes entry probabilities and rates of an n-state APH distribution
    (bidiagonal generator) to minimize weighted moment distance.

    Args:
        E: Moments E[X^k]
        n: Number of phases

    Returns:
        Tuple of (MAP, Eapx, x) where MAP is (D0, D1),
        Eapx are the fitted moments, x is the solution vector
    """
    from ..mam import map_normalize

    E = np.asarray(E).ravel()

    # Scale moments to E[X]=1
    Eunscaled = E.copy()
    Escale = np.array([E[0] ** k for k in range(1, len(E) + 1)])
    E = E / Escale

    w = np.power(np.log(E[-1]) ** (1.0 / len(E)), -np.arange(1, len(E) + 1))

    # Initial point
    x0 = np.random.rand(2 * n)
    x0[:n] = x0[:n] / np.sum(x0[:n])

    TOL = 1e-8
    EPSTOL = 100 * TOL
    MAXITER = 100

    def topar(x):
        alpha = x[:n]
        # APH: bidiagonal T (diagonal + superdiagonal)
        T = np.diag(-1.0 / x[n:])
        if n > 1:
            T += np.diag(1.0 / x[n:n + n - 1], 1)
        return alpha, T

    def _compute_moments(alpha, T):
        n_local = len(alpha)
        try:
            invT = np.linalg.inv(-T)
        except np.linalg.LinAlgError:
            invT = np.linalg.pinv(-T)
        e = np.ones(n_local)
        Eapx = np.zeros(len(E))
        invT_k = invT.copy()
        for j in range(len(E)):
            Eapx[j] = factorial(j + 1) * alpha @ invT_k @ e
            invT_k = invT_k @ invT
        return Eapx

    def nnlcon(x):
        alpha, T = topar(x)
        Eapx = _compute_moments(alpha, T)
        c = -alpha  # alpha >= 0
        ceq = np.zeros(3 + n)
        ceq[:3] = w[:3] * np.abs(np.log(E[:3]) - np.log(np.maximum(Eapx[:3], 1e-300)))
        ceq[3:] = np.sum(alpha) - 1
        return c, ceq

    def objfun(x):
        alpha, T = topar(x)
        Eapx = _compute_moments(alpha, T)
        Eapx = np.maximum(Eapx, 1e-300)
        return w @ np.abs(np.log(E) - np.log(Eapx))

    def ineq_fun(x):
        c, ceq = nnlcon(x)
        return -c

    def eq_fun(x):
        c, ceq = nnlcon(x)
        return ceq

    bounds = [(EPSTOL, None)] * (2 * n)

    try:
        res = minimize(objfun, x0, method='SLSQP',
                       bounds=bounds,
                       constraints=[
                           {'type': 'ineq', 'fun': ineq_fun},
                           {'type': 'eq', 'fun': eq_fun},
                       ],
                       options={'maxiter': MAXITER, 'ftol': TOL, 'disp': False})
        x = res.x
    except Exception:
        x = x0

    alpha, T = topar(x)
    D0 = T
    D1 = -T @ (np.ones((n, 1)) @ alpha.reshape(1, -1))
    MAP = _map_scale_to_mean(D0, D1, Eunscaled[0])
    MAP = map_normalize(MAP[0], MAP[1])
    Eapx = _compute_moments(alpha, T) * Escale
    return MAP, Eapx, x


# ============================================================================
# PH auto fitting (complete version with approximate methods)
# ============================================================================

def kpcfit_ph_auto(E: np.ndarray, options: KpcfitPhOptions = None
                   ) -> List[Tuple[Tuple[np.ndarray, np.ndarray], float, Any, str]]:
    """
    Automatic PH distribution fitting using exact and approximate methods.

    Performs three phases:
    1) Exact fitting (Prony's method for hyperexp, APH for low variability)
    2) Approximate fitting in moment space (KPC-based, kpcfit_ph_search)
    3) Approximate fitting in parameter space (hyperexp or APH direct)

    Args:
        E: Moments E[X^k]
        options: Fitting options (default: auto-generated)

    Returns:
        List of (PH_distribution, distance, params, method) tuples
        sorted by fitting quality. PH_distribution is a (D0, D1) tuple.
    """
    from ..mam import map_isfeasible, map_moment
    import math

    E = np.asarray(E).ravel()

    if options is None:
        options = kpcfit_ph_options(E)

    if options.verbose:
        print(f"kpcfit_ph: version {kpcfit_version()}")
        print("kpcfit_ph: type 'help kpcfit_ph_options' for options")

    # Scale moments to E[X]=1
    E_unscaled = E.copy()
    E_scale = np.array([E[0] ** (k + 1) for k in range(len(E))])
    E_normalized = E / E_scale

    # --- distance function ---
    def dist_fun(E_ref, E_apx):
        E_ref = np.asarray(E_ref)
        E_apx = np.asarray(E_apx)
        # Avoid log of zero/negative
        E_ref_safe = np.maximum(np.abs(E_ref), 1e-300)
        E_apx_safe = np.maximum(np.abs(E_apx), 1e-300)
        w = np.power(np.log10(E_ref_safe[-1]) ** (1.0 / len(E_ref)),
                     -np.arange(1, len(E_ref) + 1))
        return np.dot(w, np.abs(np.log10(E_ref_safe) - np.log10(E_apx_safe)))

    # ==================== Phase 1: Exact fitting ====================
    if options.verbose:
        print("\nkpcfit_ph: starting exact fitting methods")
    PH_EXACT = kpcfit_ph_exact(E_normalized, options)

    # ==================== Phase 2: Approximate fitting - moment space ====================
    if options.verbose:
        print("\nkpcfit_ph: starting approximate fitting method -- "
              "moment space (KPC-based)")
        print("\t\t\tRUN\t\tORD\tDIST\t\tRESULT", end='')

    PH_APX_MS = []  # Each entry: [MAP, dist, x]
    t0 = time.time()

    min_J = max(2, int(math.ceil(math.log2(options.min_num_states))))
    max_J = int(math.ceil(math.log2(options.max_num_states)))

    for J in range(min_J, max_J + 1):
        mindist = np.inf
        best = 0
        bestidx = 0
        randomruns = int(math.ceil(2 * options.runs / 3))

        for run in range(1, options.runs + 1):
            if options.verbose:
                suffix = 'r' if run > randomruns else ''
                print(f"\n\t\t\t{run}{suffix}\t", end='')
                print(f"\t{2**J}", end='')

            x0_val = None
            if run > randomruns and best > 0 and bestidx > 0:
                x0_val = PH_APX_MS[bestidx - 1][2]

            try:
                APX, score_val, x_val = kpcfit_ph_search(
                    E_normalized, J, options, x0_val)

                if map_isfeasible(APX[0], APX[1]):
                    # Compute dist
                    apx_moments = np.array([map_moment(APX[0], APX[1], k + 1)
                                            for k in range(len(E_normalized))])
                    dist = dist_fun(E_normalized, apx_moments)

                    if np.isnan(dist):
                        if options.verbose:
                            print(f"\t-------- \tx no solution found.", end='')
                    else:
                        if options.verbose:
                            print(f"\t{dist:.6f}", end='')
                        if dist < mindist:
                            if best == 0:
                                PH_APX_MS.append([APX, dist, x_val])
                                if options.verbose:
                                    n_states = APX[0].shape[0]
                                    print(f"\to PH({n_states}) saved. "
                                          f"\t elapsed: {time.time()-t0:.3f}s",
                                          end='')
                            else:
                                PH_APX_MS[-1] = [APX, dist, x_val]
                                if options.verbose:
                                    n_states = APX[0].shape[0]
                                    print(f"\t+ PH({n_states}) updated. "
                                          f"\t elapsed: {time.time()-t0:.3f}s",
                                          end='')
                            best = run
                            bestidx = len(PH_APX_MS)
                            mindist = dist
                        else:
                            if options.verbose:
                                print(f"\t\t\t\t\t\t elapsed: "
                                      f"{time.time()-t0:.3f}s", end='')
                else:
                    if options.verbose:
                        print(f"\t\t\t\tx infeasible result.", end='')
            except Exception:
                if options.verbose:
                    print(f"\t\t\t\tx optimization solver failed.", end='')

    # ==================== Phase 3: Approximate fitting - parameter space ====================
    SCV = (E_normalized[1] - E_normalized[0] ** 2) / E_normalized[0] ** 2
    PH_APX_PS = []

    if SCV >= 1:
        if options.verbose:
            print("\n\nkpcfit_ph: starting approximate fitting method -- "
                  "hyper-exponential parameter space (non-KPC-based)")
            print("\t\t\tRUN\t\tORD\tDIST\t\tRESULT", end='')

        for J in range(min_J, max_J + 1):
            mindist = np.inf
            best = 0
            bestidx = 0
            randomruns = int(math.ceil(2 * options.runs / 3))

            for run in range(1, options.runs + 1):
                if options.verbose:
                    print(f"\n\t\t\t{run}\t", end='')
                    print(f"\t{2**J}", end='')
                try:
                    APX, Eapx, x_val = psfit_hyperexp_weighted(
                        E_normalized, 2 ** J)
                    if map_isfeasible(APX[0], APX[1]):
                        apx_moments = np.array([
                            map_moment(APX[0], APX[1], k + 1)
                            for k in range(len(E_normalized))])
                        dist = dist_fun(E_normalized, apx_moments)
                        if np.isnan(dist):
                            if options.verbose:
                                print(f"\t-------- \tx no solution found.",
                                      end='')
                        else:
                            if options.verbose:
                                print(f"\t{dist:.6f}", end='')
                            if dist < mindist:
                                if best == 0:
                                    PH_APX_PS.append([APX, dist, x_val])
                                else:
                                    PH_APX_PS[-1] = [APX, dist, x_val]
                                best = run
                                bestidx = len(PH_APX_PS)
                                mindist = dist
                                if options.verbose:
                                    print(f"\to saved.", end='')
                    else:
                        if options.verbose:
                            print(f"\t\t\t\tinfeasible result.", end='')
                except Exception:
                    if options.verbose:
                        print(f"\t\t\t\tx optimization solver failed.", end='')
    else:
        if options.verbose:
            print("\n\nkpcfit_ph: starting approximate fitting method -- "
                  "aph parameter space (non-KPC-based)")
            print("\t\t\tRUN\t\tORD\tDIST\t\tRESULT", end='')

        for J in range(min_J, max_J + 1):
            mindist = np.inf
            best = 0
            bestidx = 0
            randomruns = int(math.ceil(2 * options.runs / 3))

            for run in range(1, options.runs + 1):
                if options.verbose:
                    print(f"\n\t\t\t{run}\t", end='')
                    print(f"\t{2**J}", end='')
                try:
                    APX, Eapx, x_val = psfit_aph_weighted(
                        E_normalized, 2 ** J)
                    if map_isfeasible(APX[0], APX[1]):
                        apx_moments = np.array([
                            map_moment(APX[0], APX[1], k + 1)
                            for k in range(len(E_normalized))])
                        dist = dist_fun(E_normalized, apx_moments)
                        if np.isnan(dist):
                            if options.verbose:
                                print(f"\t-------- \tx no solution found.",
                                      end='')
                        else:
                            if options.verbose:
                                print(f"\t{dist:.6f}", end='')
                            if dist < mindist:
                                if best == 0:
                                    PH_APX_PS.append([APX, dist, x_val])
                                else:
                                    PH_APX_PS[-1] = [APX, dist, x_val]
                                best = run
                                bestidx = len(PH_APX_PS)
                                mindist = dist
                                if options.verbose:
                                    print(f"\to saved.", end='')
                    else:
                        if options.verbose:
                            print(f"\t\t\t\tinfeasible result.", end='')
                except Exception:
                    if options.verbose:
                        print(f"\t\t\t\tx optimization solver failed.", end='')

    # ==================== Collect all results ====================
    nExactResults = len(PH_EXACT)
    nApproxMS = len(PH_APX_MS)
    nApproxPS = len(PH_APX_PS)

    results = []

    for j, ph in enumerate(PH_EXACT):
        ph_scaled = _map_scale_to_mean(ph[0], ph[1], E_unscaled[0])
        try:
            ph_moments = np.array([map_moment(ph[0], ph[1], k + 1)
                                   for k in range(len(E_normalized))])
            dist = dist_fun(E_normalized, ph_moments)
        except Exception:
            dist = np.inf
        results.append((ph_scaled, dist, None, 'exact'))

    for j, entry in enumerate(PH_APX_MS):
        ph = entry[0]
        ph_scaled = _map_scale_to_mean(ph[0], ph[1], E_unscaled[0])
        try:
            ph_moments = np.array([map_moment(ph[0], ph[1], k + 1)
                                   for k in range(len(E_normalized))])
            dist = dist_fun(E_normalized, ph_moments)
        except Exception:
            dist = np.inf
        results.append((ph_scaled, dist, entry[2], 'approx_moment_space'))

    for j, entry in enumerate(PH_APX_PS):
        ph = entry[0]
        ph_scaled = _map_scale_to_mean(ph[0], ph[1], E_unscaled[0])
        try:
            ph_moments = np.array([map_moment(ph[0], ph[1], k + 1)
                                   for k in range(len(E_normalized))])
            dist = dist_fun(E_normalized, ph_moments)
        except Exception:
            dist = np.inf
        results.append((ph_scaled, dist, entry[2], 'approx_param_space'))

    if options.verbose:
        print(f"\n\nReturned {nExactResults} PH distribution(s) by exact "
              f"moment matching.")
        print(f"Returned {nApproxMS} PH distribution(s) by approximate "
              f"moment matching (moment space).")
        print(f"Returned {nApproxPS} PH distribution(s) by approximate "
              f"moment matching (parameter space).")
        print()

    return results


# ============================================================================
# Manual PH fitting (kpcfit_ph_manual)
# ============================================================================

def kpcfit_ph_manual(E: np.ndarray,
                     NumStates: int = 4,
                     NumPHs: int = 2,
                     NumMoms: int = 6,
                     MaxIterMM: int = 100,
                     WantHyper: bool = False
                     ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Manual PH fitting via KPC composition with user-specified parameters.

    Fits a PH distribution by composing NumPHs PH(2) distributions using
    Kronecker products, optimizing to match the given moments.

    Args:
        E: Moments E[X^k]
        NumStates: Number of states (must be even, default: 4)
        NumPHs: Number of PH(2)s to compose (default: 2)
        NumMoms: Number of moments to consider (default: 6)
        MaxIterMM: Maximum iterations (default: 100)
        WantHyper: If True, prefer hyperexponential structure

    Returns:
        Tuple (D0, D1) representing the fitted PH distribution
    """
    from ..mam import (map_normalize, map_moment, map_feasblock, map_kpc,
                       map_rand, map_isfeasible)
    from ..kpctoolbox import mmpp2_fit3

    E = np.asarray(E).ravel()

    if NumStates % 2 != 0:
        raise ValueError("kpcfit_ph: number of states must be multiple of 2")

    TOL = 1e-10
    EPSTOL = 100 * TOL
    n = len(E)
    nexact = 2

    # Initialize PH components
    PH = [None] * NumPHs
    x0 = np.zeros(3 * NumPHs)

    for j in range(NumPHs):
        if j == 0:
            # Random MAP
            D0_r = np.random.rand(2, 2)
            D1_r = np.random.rand(2, 2)
            from ..mam import map_normalize as mn
            PH[j] = mn(D0_r, D1_r)
        else:
            PH[j] = map_feasblock(np.random.rand(),
                                  1000 * np.random.rand(), -1, 0)
        for i in range(3):
            x0[i * NumPHs + j] = map_moment(PH[j][0], PH[j][1], i + 1)

    # Stagnation tracking
    stagnval = [0.0]
    stagniter = [0]
    lgkx = [x0.copy()]

    def xtopar(x):
        x_abs = np.abs(x) + EPSTOL
        E1j = x_abs[:NumPHs]
        SCVj = x_abs[NumPHs:2 * NumPHs]
        E2j = (1 + SCVj) * E1j ** 2
        E3j = x_abs[2 * NumPHs:3 * NumPHs]

        # Fit first PH
        nonlocal PH
        if WantHyper:
            PH[0] = map_feasblock(E1j[0], E2j[0], E3j[0], np.random.rand())
        else:
            try:
                PH[0] = mmpp2_fit3(E1j[0], E2j[0], E3j[0], np.random.rand())
                PH[0] = (np.real(PH[0][0]), np.real(PH[0][1]))
            except Exception:
                PH[0] = map_feasblock(E1j[0], E2j[0], E3j[0], np.random.rand())
        PH[0] = _map_scale_to_mean(PH[0][0], PH[0][1], E1j[0])

        result = PH[0]
        for j_idx in range(1, NumPHs - 1):
            PHk = map_feasblock(E1j[j_idx], E2j[j_idx], E3j[j_idx],
                                np.random.rand())
            result = map_kpc([result, PHk])

        # Adjust last PH to match overall moments
        m1_cur = map_moment(result[0], result[1], 1)
        m2_cur = map_moment(result[0], result[1], 2)
        if m1_cur > 0:
            E1j[NumPHs - 1] = E[0] / m1_cur
        if m2_cur > 0:
            E2j[NumPHs - 1] = E[1] * 2 / m2_cur
        if n > 2:
            m3_cur = map_moment(result[0], result[1], 3)
            if m3_cur > 0:
                E3j[NumPHs - 1] = E[2] * 6 / m3_cur

        last_ph = map_feasblock(E1j[NumPHs - 1], E2j[NumPHs - 1],
                                E3j[NumPHs - 1], 0)
        result = map_kpc([result, last_ph])

        # Make diagonal D0
        D0 = np.diag(np.diag(result[0]))
        D1 = result[1]
        # Re-normalize
        row_sums = D0 + D1
        for i_row in range(D0.shape[0]):
            s = np.sum(row_sums[i_row, :])
            if s != 0:
                D1[i_row, :] *= (-D0[i_row, i_row]) / s if s > 0 else 1

        PH_final = (D0, D1)

        Eest = np.zeros(n)
        for k in range(n):
            Eest[k] = map_moment(PH_final[0], PH_final[1], k + 1)
        return Eest

    def objfun(x):
        Eest = xtopar(x)
        f = 0.0
        for k in range(nexact, n):
            if Eest[k] > 0:
                f += abs(np.log10(E[k]) - np.log10(Eest[k]))
        return f

    def eq_con(x):
        Eest = xtopar(x)
        ceq = []
        for k in range(2, nexact):
            ceq.append(np.linalg.norm(E[k] - Eest[k], 1))
        return np.array(ceq) if ceq else np.array([0.0])

    def _callback(xk):
        fval = objfun(xk)
        if not np.isnan(fval) and np.sum(np.isnan(xk)) == 0:
            lgkx[0] = xk.copy()
        else:
            raise StopIteration
        if stagnval[0] == 0:
            stagnval[0] = fval
        delta = abs(fval - stagnval[0]) / max(abs(stagnval[0]), 1e-30)
        if delta < 0.01:
            stagniter[0] += 1
            if stagniter[0] >= 100:
                raise StopIteration
        else:
            stagniter[0] = 0
        stagnval[0] = fval

    bounds = [(EPSTOL, None)] * len(x0)
    constraints_list = []
    if nexact > 2:
        constraints_list.append({'type': 'eq', 'fun': eq_con})

    try:
        res = minimize(objfun, x0, method='SLSQP',
                       bounds=bounds,
                       constraints=constraints_list,
                       options={'maxiter': MaxIterMM, 'ftol': TOL,
                                'disp': False},
                       callback=_callback)
    except StopIteration:
        pass

    # Reconstruct from best known solution
    xtopar(lgkx[0])

    # Get the current PH and compose
    # Re-run xtopar to set PH properly, then compose
    x_best = lgkx[0]
    x_abs = np.abs(x_best) + EPSTOL
    E1j = x_abs[:NumPHs]
    SCVj = x_abs[NumPHs:2 * NumPHs]
    E2j = (1 + SCVj) * E1j ** 2
    E3j = x_abs[2 * NumPHs:3 * NumPHs]

    if WantHyper:
        PH[0] = map_feasblock(E1j[0], E2j[0], E3j[0], np.random.rand())
    else:
        try:
            PH[0] = mmpp2_fit3(E1j[0], E2j[0], E3j[0], np.random.rand())
            PH[0] = (np.real(PH[0][0]), np.real(PH[0][1]))
        except Exception:
            PH[0] = map_feasblock(E1j[0], E2j[0], E3j[0], np.random.rand())
    PH[0] = _map_scale_to_mean(PH[0][0], PH[0][1], E1j[0])

    result = PH[0]
    for j_idx in range(1, NumPHs - 1):
        PHk = map_feasblock(E1j[j_idx], E2j[j_idx], E3j[j_idx],
                            np.random.rand())
        result = map_kpc([result, PHk])

    m1_cur = map_moment(result[0], result[1], 1)
    m2_cur = map_moment(result[0], result[1], 2)
    if m1_cur > 0:
        E1j[NumPHs - 1] = E[0] / m1_cur
    if m2_cur > 0:
        E2j[NumPHs - 1] = E[1] * 2 / m2_cur
    if n > 2:
        m3_cur = map_moment(result[0], result[1], 3)
        if m3_cur > 0:
            E3j[NumPHs - 1] = E[2] * 6 / m3_cur

    last_ph = map_feasblock(E1j[NumPHs - 1], E2j[NumPHs - 1],
                            E3j[NumPHs - 1], 0)
    result = map_kpc([result, last_ph])

    D0 = np.diag(np.diag(result[0]))
    D1 = result[1]

    from ..mam import map_normalize as mn
    PH_result = mn(D0, D1)
    PH_result = _map_scale_to_mean(PH_result[0], PH_result[1], E[0])

    if map_isfeasible(PH_result[0], PH_result[1]):
        pass
    else:
        warnings.warn("PH distribution may not be valid, check for negative "
                      "rates in D1 or in off-diagonal of D0")

    PH_result = _map_scale_to_mean(PH_result[0], PH_result[1], E[0])
    return PH_result


# ============================================================================
# PH summary (kpcfit_ph_summary)
# ============================================================================

def kpcfit_ph_summary(PH_list: List[Tuple[Tuple[np.ndarray, np.ndarray],
                                           float, Any, str]],
                      E: np.ndarray) -> None:
    """
    Print summary of PH fitting results.

    Displays moment comparison table for exact and approximate fits.

    Args:
        PH_list: Results from kpcfit_ph_auto
        E: Target moments
    """
    from ..mam import map_moment

    E = np.asarray(E).ravel()

    exact = [(i, ph) for i, ph in enumerate(PH_list) if ph[3] == 'exact']
    apxms = [(i, ph) for i, ph in enumerate(PH_list)
             if ph[3] == 'approx_moment_space']
    apxps = [(i, ph) for i, ph in enumerate(PH_list)
             if ph[3] == 'approx_param_space']

    # Header
    header = "Log10\tTRACE        "
    for j, (idx, _) in enumerate(exact):
        header += f"\tEXMM (idx={j+1})   "
    for j, (idx, _) in enumerate(apxms):
        header += f"\tAPXMS (idx={len(exact)+j+1})  "
    for j, (idx, _) in enumerate(apxps):
        header += f"\tAPXPS (idx={len(exact)+len(apxms)+j+1})  "
    print(header)

    for k in range(len(E)):
        if k == 0:
            line = "E[X]"
        else:
            line = f"E[X^{k+1}]"
        line += f"\t{np.log10(max(E[k], 1e-300)):13.6f}"
        for j, (idx, ph) in enumerate(exact):
            mk = map_moment(ph[0][0], ph[0][1], k + 1)
            line += f"\t{np.log10(max(mk, 1e-300)):13.6f}"
        for j, (idx, ph) in enumerate(apxms):
            mk = map_moment(ph[0][0], ph[0][1], k + 1)
            line += f"\t{np.log10(max(mk + 1e-300, 1e-300)):13.6f}"
        for j, (idx, ph) in enumerate(apxps):
            mk = map_moment(ph[0][0], ph[0][1], k + 1)
            line += f"\t{np.log10(max(mk + 1e-300, 1e-300)):13.6f}"
        print(line)


__all__ = [
    'KPCFIT_TOL',
    'KpcfitTraceData',
    'KpcfitPhOptions',
    'KpcfitResult',
    'kpcfit_tol',
    'kpcfit_version',
    'logspacei',
    'kpcfit_init',
    'kpcfit_sub_eval_acfit',
    'kpcfit_sub_acfit',
    'kpcfit_sub_bcfit',
    'kpcfit_sub_bic',
    'kpcfit_sub_compose',
    'kpcfit_hyper_charpoly',
    'kpcfit_ph_prony',
    'kpcfit_ph_options',
    'kpcfit_ph_exact',
    'kpcfit_ph_search',
    'kpcfit_ph_auto',
    'kpcfit_ph_manual',
    'kpcfit_ph_summary',
    'kpcfit_auto',
    'kpcfit_manual',
    'psfit_hyperexp_weighted',
    'psfit_aph_weighted',
    '_map_scale_to_mean',
]
