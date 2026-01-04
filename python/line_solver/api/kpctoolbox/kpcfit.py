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
from typing import Union, Tuple, Optional, List, Dict, Any
from dataclasses import dataclass, field
from scipy.optimize import minimize, NonlinearConstraint
from math import factorial

ArrayLike = Union[np.ndarray, list]

# Default tolerance for KPC fitting
KPCFIT_TOL = 1e-10


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
    """Options for PH distribution fitting."""
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
    AC = trace_acf(S, ac_lags)
    ACFull = trace_acf(S, np.arange(1, n // n_min_support_ac + 1))

    # Apply smoothing if requested
    if smooth > 0:
        try:
            from scipy.ndimage import uniform_filter1d
            AC = uniform_filter1d(AC, size=smooth, mode='nearest')
            ACFull = uniform_filter1d(ACFull, size=smooth, mode='nearest')
        except ImportError:
            pass  # Skip smoothing if scipy not available

    # Truncate AC where values become negligible
    small_idx = np.where(np.abs(AC) < 1e-6)[0]
    if len(small_idx) > 0:
        posmax = ac_lags[small_idx[0]]
        keep_mask = ac_lags <= posmax
        ac_lags = ac_lags[keep_mask]
        AC = AC[keep_mask]
        bc_grid_lags = bc_grid_lags[bc_grid_lags <= posmax]

    # Compute bicovariances
    BC, BCLags = trace_bicov(S, bc_grid_lags)

    return KpcfitTraceData(
        S=S, E=E, AC=AC, ACFull=ACFull, ACLags=ac_lags,
        BC=BC, BCGridLags=bc_grid_lags, BCLags=BCLags
    )


def kpcfit_sub_eval_acfit(SCV: np.ndarray, G2: np.ndarray,
                           acf_lags: np.ndarray) -> Tuple[float, np.ndarray]:
    """
    Evaluate autocorrelation fit for given SCV and G2 parameters.

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


def kpcfit_sub_bic(SA: np.ndarray, orders: ArrayLike) -> int:
    """
    Select MAP order using BIC criterion.

    Args:
        SA: Autocorrelation sequence
        orders: Candidate orders (powers of 2)

    Returns:
        Recommended number of MAPs (log2 of states)
    """
    orders = np.asarray(orders)
    nlags = len(SA)

    # Find where AC decays below threshold
    nlagsend_idx = np.where(SA < 1e-6)[0]
    ordermax = max(orders)

    if len(nlagsend_idx) == 0:
        nlagsend = nlags
    else:
        nlagsend = nlagsend_idx[0] - 1 + ordermax

    NLAGSMAX = 10000

    if nlagsend > NLAGSMAX:
        SA_lags = logspacei(1, nlagsend - ordermax, NLAGSMAX)
    else:
        if nlagsend > ordermax:
            SA_lags = np.arange(1, nlagsend - ordermax + 1)
        else:
            ordermax = max(1, nlagsend - 2)
            SA_lags = np.arange(1, max(1, nlagsend - ordermax) + 1)
            orders = orders[orders <= ordermax]

    if len(SA_lags) == 0:
        return int(np.log2(orders[0])) if len(orders) > 0 else 2

    n_samples = len(SA_lags)
    SA_lags_y = SA_lags.copy()

    # Ensure indices are within bounds
    max_idx = len(SA) - 1
    SA_lags_y = SA_lags_y[SA_lags_y <= max_idx]
    if len(SA_lags_y) == 0:
        return int(np.log2(orders[0])) if len(orders) > 0 else 2

    Y = SA[SA_lags_y - 1] if SA_lags_y[0] > 0 else SA[SA_lags_y]

    SBC = []
    for order in orders:
        order = int(order)

        # Determine valid range for this order
        max_valid_lag = max_idx - order
        if max_valid_lag < 1:
            SBC.append(np.inf)
            continue

        valid_lags = SA_lags_y[SA_lags_y <= max_valid_lag]
        if len(valid_lags) == 0:
            SBC.append(np.inf)
            continue

        # Build X matrix with consistent dimensions
        X = np.zeros((len(valid_lags), order))
        for i in range(order):
            lags_x = valid_lags + (i + 1)
            X[:, i] = SA[lags_x - 1] if valid_lags[0] > 0 else SA[lags_x]

        Y_fit = SA[valid_lags - 1] if valid_lags[0] > 0 else SA[valid_lags]

        try:
            b, residuals, rank, s = np.linalg.lstsq(X, Y_fit, rcond=None)
            r = Y_fit - X @ b
            sse = np.sum(r ** 2)
            n = len(Y_fit)
            sbc = n * np.log(max(sse, 1e-300)) - n * np.log(n) + np.log(n) * order
            SBC.append(sbc)
        except:
            SBC.append(np.inf)

    if len(SBC) == 0 or all(np.isinf(SBC)):
        return 3  # Default

    best_idx = np.argmin(SBC)
    return int(np.log2(orders[best_idx]))


def kpcfit_sub_compose(E1j: np.ndarray, SCVj: np.ndarray,
                        E3j: np.ndarray, G2j: np.ndarray,
                        verbose: int = 0) -> Tuple[Optional[Tuple[np.ndarray, np.ndarray]],
                                                    List[Tuple[np.ndarray, np.ndarray]], int]:
    """
    Compose a large MAP from smaller MAP(2)s using Kronecker products.

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
    except:
        kpc_map = None

    if kpc_map is None or not map_isfeasible(kpc_map):
        if SCVj[0] < 0.5:
            kpc_map = map_erlang(E1j[0], 2)
            if verbose > 0:
                print("MAP 1 is erlang-2")
        else:
            if verbose > 0:
                print(f"MAP 1 has presumably infeasible E3: {E3j[0]}")
            E2_1 = (1 + SCVj[0]) * E1j[0] ** 2
            kpc_map, fit_err = map2_fit(E1j[0], E2_1, -1, G2j[0])
            if kpc_map is None or fit_err != 0:
                if verbose > 0:
                    print(f"MAP 1 has infeasible G2: {G2j[0]}")
                kpc_map, fit_err = map2_fit(E1j[0], E2_1, -1, 0)
                if fit_err != 0 or kpc_map is None:
                    return None, [], 1

    sub_maps.append(kpc_map)

    for j in range(1, J):
        E2_j = (1 + SCVj[j]) * E1j[j] ** 2
        map_j = map_feasblock(E1j[j], E2_j, E3j[j], G2j[j])

        try:
            feasible = map_isfeasible(map_j)
        except:
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
                        print(f"MAP {j+1} has infeasible G2: {G2j[j]}")
                    if fit_err != 0 or map_j is None:
                        return None, sub_maps, 5

        if map_j is None:
            if verbose > 0:
                print(f"Replacing MAP {j+1} with exponential")
            map_j = map_exponential(E1j[j])

        sub_maps.append(map_j)
        kpc_map = map_kpc([kpc_map, map_j])

    # Check feasibility of all sub-MAPs
    for j, smap in enumerate(sub_maps):
        if not map_isfeasible(smap):
            if verbose > 0:
                print(f"MAP {j+1} is infeasible")
            return None, sub_maps, 10

    kpc_map = map_normalize(kpc_map)
    return kpc_map, sub_maps, 0


def kpcfit_hyper_charpoly(E: np.ndarray, n: int) -> np.ndarray:
    """
    Compute characteristic polynomial coefficients for hyperexponential fitting.

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


def kpcfit_ph_options(E: np.ndarray, **kwargs) -> KpcfitPhOptions:
    """
    Create options for PH fitting.

    Args:
        E: Moments vector
        **kwargs: Option overrides
            - verbose: Print progress (default: True)
            - runs: Number of optimization runs (default: 5)
            - min_num_states: Minimum states (default: 2)
            - max_num_states: Maximum states (default: 32)
            - min_exact_mom: Minimum moments to fit exactly (default: 3)

    Returns:
        KpcfitPhOptions dataclass
    """
    options = KpcfitPhOptions()

    for key, value in kwargs.items():
        key_lower = key.lower().replace('_', '')
        if hasattr(options, key):
            setattr(options, key, value)
        elif key_lower == 'minnumstates':
            options.min_num_states = value
        elif key_lower == 'maxnumstates':
            options.max_num_states = value
        elif key_lower == 'minexactmom':
            options.min_exact_mom = value

    # Ensure states are powers of 2
    import math
    if options.min_num_states > 0:
        options.min_num_states = 2 ** math.ceil(math.log2(options.min_num_states))
    if options.max_num_states > 0:
        options.max_num_states = 2 ** math.ceil(math.log2(options.max_num_states))

    # Check moment requirements
    E = np.asarray(E)
    if 2 * options.max_num_states - 1 > len(E):
        raise ValueError(f"MaxNumStates of {options.max_num_states} requires "
                        f"at least {2 * options.max_num_states - 1} moments.")

    return options


def kpcfit_ph_exact(E: np.ndarray, options: KpcfitPhOptions
                    ) -> List[Tuple[np.ndarray, np.ndarray]]:
    """
    Exact PH fitting methods (Prony's method for hyperexp, APH for low variability).

    Args:
        E: Moments E[X^k]
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
            print(f"kpcfit_ph: HIGHER variability than exponential (SCV = {SCV:.4f})")
            print("kpcfit_ph: starting exact hyper-exponential fitting (Prony's method)")

        for n in range(2, options.max_num_states + 1):
            if len(E) < 2 * n - 1:
                if options.verbose:
                    print(f"kpcfit_ph: not enough moments for hyper-exp({n})")
                break

            PH = kpcfit_ph_prony(E, n)
            if map_isfeasible(PH):
                PH_EXACT.append(PH)
                if options.verbose:
                    print(f"\t\thyper-exp({n}): feasible, matched {2*n-1} moments. saved.")
            else:
                if options.verbose:
                    print(f"\t\thyper-exp({n}): infeasible")
                break

    elif SCV < 1:
        # Lower variability
        if options.verbose:
            print(f"kpcfit_ph: LOWER variability than exponential (SCV = {SCV:.4f})")

        n = 1
        while 1 / n > SCV:
            n += 1

        if options.verbose:
            print(f"kpcfit_ph: exact E[X^2] fitting requires at least {n} states")

        if n == 2:
            if options.verbose:
                print("kpcfit_ph: attempting PH(2) fitting")
            PH, err = map2_fit(E[0], E[1], E[2] if len(E) > 2 else -1, 0)
            if PH is None:
                PH, err = map2_fit(E[0], E[1], -1, 0)
            if PH is not None and map_isfeasible(PH):
                PH_EXACT.append(PH)
                if options.verbose:
                    print("\t\tph(2): feasible. saved.")
        elif abs(SCV - 1/n) < KPCFIT_TOL:
            # Erlang case
            ERL = map_erlang(E[0], n)
            if options.verbose:
                print(f"kpcfit_ph: erlang moment set. fitted erlang-{n}. saved.")
            PH_EXACT.append(ERL)
        else:
            # APH fitting
            maxorder = options.max_num_states
            if options.verbose:
                print(f"kpcfit_ph: fitting APH (max order = {maxorder})")
            try:
                PH = aph_fit(E[0], E[1], E[2] if len(E) > 2 else E[1] * 1.5, maxorder)
                if map_isfeasible(PH):
                    PH_EXACT.append(PH)
                    if options.verbose:
                        print(f"\t\taph({len(PH[0])}): feasible. saved.")
            except:
                if options.verbose:
                    print("kpcfit_ph: cannot fit APH distribution")
    else:
        # SCV == 1 (exponential)
        if options.verbose:
            print(f"kpcfit_ph: SAME variability as exponential (SCV = {SCV:.4f})")
        EXP = (np.array([[-1/E[0]]]), np.array([[1/E[0]]]))
        PH_EXACT.append(EXP)
        if options.verbose:
            print("kpcfit_ph: exponential. saved.")

    return PH_EXACT


def kpcfit_ph_auto(E: np.ndarray, options: KpcfitPhOptions = None
                   ) -> List[Tuple[Tuple[np.ndarray, np.ndarray], float, Any, str]]:
    """
    Automatic PH distribution fitting using exact and approximate methods.

    Args:
        E: Moments E[X^k]
        options: Fitting options (default: auto-generated)

    Returns:
        List of (PH_distribution, distance, params, method) tuples
    """
    from ..mam import map_isfeasible, map_scale, map_moment

    E = np.asarray(E).ravel()

    if options is None:
        options = kpcfit_ph_options(E)

    if options.verbose:
        print(f"kpcfit_ph: version 1.0 (Python)")

    # Scale moments to E[X]=1
    E_unscaled = E.copy()
    E_scale = np.array([E[0] ** (k + 1) for k in range(len(E))])
    E_normalized = E / E_scale

    # Exact fitting
    if options.verbose:
        print("\nkpcfit_ph: starting exact fitting methods")
    PH_EXACT = kpcfit_ph_exact(E_normalized, options)

    # Collect results
    results = []

    def dist_fun(E_ref, E_apx):
        E_ref = np.asarray(E_ref)
        E_apx = np.asarray(E_apx)
        w = np.power(np.log10(E_ref[-1]) ** (1 / len(E_ref)),
                     -np.arange(1, len(E_ref) + 1))
        return np.dot(w, np.abs(np.log10(E_ref) - np.log10(E_apx)))

    for ph in PH_EXACT:
        ph_scaled = map_scale(ph, E_unscaled[0])
        try:
            ph_moments = map_moment(ph, list(range(1, len(E) + 1)))
            dist = dist_fun(E_normalized, ph_moments / E_scale)
        except:
            dist = np.inf
        results.append((ph_scaled, dist, None, 'exact'))

    return results


__all__ = [
    'KPCFIT_TOL',
    'KpcfitTraceData',
    'KpcfitPhOptions',
    'KpcfitResult',
    'kpcfit_tol',
    'logspacei',
    'kpcfit_init',
    'kpcfit_sub_eval_acfit',
    'kpcfit_sub_bic',
    'kpcfit_sub_compose',
    'kpcfit_hyper_charpoly',
    'kpcfit_ph_prony',
    'kpcfit_ph_options',
    'kpcfit_ph_exact',
    'kpcfit_ph_auto',
]
