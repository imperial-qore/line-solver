"""
KPCToolbox library functions - comprehensive utilities:
- Trace analysis (single and multi-trace statistics)
- MMPP (Markov Modulated Poisson Process) fitting
- KPC fitting (Kronecker Product Composition)
- MVPH (Multivariate Phase-Type) analysis
- CTMC/DTMC utilities
- Advanced probability computations

Native Python implementation using numpy/scipy. No JVM dependencies.
"""

import numpy as np


# ========== TRACE ANALYSIS (SINGLE) ==========

def lib_kpc_trace_mean(trace):
    """Compute mean of trace.

    Delegates to trace_mean from api/kpctoolbox/trace.
    """
    try:
        from line_solver.api.kpctoolbox.trace import trace_mean
        return trace_mean(trace)
    except Exception as e:
        raise RuntimeError(f"trace_mean failed: {str(e)}")


def lib_kpc_trace_var(trace):
    """Compute variance of trace (unbiased, ddof=1).

    Delegates to trace_var from api/kpctoolbox/trace.
    """
    try:
        from line_solver.api.kpctoolbox.trace import trace_var
        return trace_var(trace)
    except Exception as e:
        raise RuntimeError(f"trace_var failed: {str(e)}")


def lib_kpc_trace_scv(trace):
    """Compute squared coefficient of variation.

    Delegates to trace_scv from api/kpctoolbox/trace.
    """
    try:
        from line_solver.api.kpctoolbox.trace import trace_scv
        return trace_scv(trace)
    except Exception as e:
        raise RuntimeError(f"trace_scv failed: {str(e)}")


def lib_kpc_trace_acf(trace, max_lag):
    """Compute autocorrelation function at lags 1..max_lag.

    Delegates to trace_acf from api/kpctoolbox/trace.

    Args:
        trace: Array of observations.
        max_lag: Maximum lag to compute.

    Returns:
        numpy array of autocorrelation values at lags 1..max_lag.
    """
    try:
        from line_solver.api.kpctoolbox.trace import trace_acf
        max_lag = int(max_lag)
        if max_lag <= 0:
            return np.array([])
        lags = np.arange(1, max_lag + 1)
        return trace_acf(trace, lags)
    except Exception as e:
        raise RuntimeError(f"trace_acf failed: {str(e)}")


def lib_kpc_trace_idc(trace, time_window=None):
    """Compute index of dispersion for counts.

    Delegates to trace_idc from api/kpctoolbox/trace.
    The time_window parameter is accepted for backward compatibility
    but is not used by the proper MATLAB-matching implementation,
    which computes IDC via the index of dispersion for intervals (IDI).

    Args:
        trace: Array of inter-arrival times.
        time_window: Ignored (kept for backward compatibility).

    Returns:
        float: The IDC value.
    """
    try:
        from line_solver.api.kpctoolbox.trace import trace_idc
        return trace_idc(trace)
    except Exception as e:
        raise RuntimeError(f"trace_idc failed: {str(e)}")


def lib_kpc_trace_iat2counts(trace, num_bins):
    """Convert inter-arrival times to counting process.

    Delegates to trace_iat2counts from api/kpctoolbox/trace.
    The num_bins parameter is converted to a time scale by dividing
    the total trace time by the number of bins.

    Args:
        trace: Array of inter-arrival times.
        num_bins: Number of bins. Converted to time scale = total_time / num_bins.

    Returns:
        numpy array of counts.
    """
    try:
        from line_solver.api.kpctoolbox.trace import trace_iat2counts
        trace_arr = np.asarray(trace, dtype=np.float64).ravel()
        num_bins = int(num_bins)
        if num_bins <= 0:
            raise ValueError("num_bins must be positive")
        total_time = np.sum(trace_arr)
        scale = total_time / num_bins
        return trace_iat2counts(trace_arr, scale)
    except Exception as e:
        raise RuntimeError(f"trace_iat2counts failed: {str(e)}")


def lib_kpc_trace_summary(trace, max_lag=10):
    """Compute comprehensive trace summary.

    Delegates to trace_summary from api/kpctoolbox/trace for the full
    MATLAB-matching summary, but returns a simplified dict for backward
    compatibility.

    Args:
        trace: Array of observations.
        max_lag: Maximum lag for autocorrelation (default: 10).

    Returns:
        dict with keys: mean, var, scv, acf (array of acf values at lags 1..max_lag).
    """
    try:
        from line_solver.api.kpctoolbox.trace import trace_mean, trace_var, trace_scv, trace_acf
        lags = np.arange(1, max_lag + 1)
        return {
            'mean': trace_mean(trace),
            'var': trace_var(trace),
            'scv': trace_scv(trace),
            'acf': trace_acf(trace, lags),
        }
    except Exception as e:
        raise RuntimeError(f"trace_summary failed: {str(e)}")


# ========== MMPP FITTING ==========

def lib_kpc_mmpp2_fit(moments, acf_lag1=None):
    """Fit 2-state MMPP from first three moments and optional lag-1 autocorrelation.

    Args:
        moments: Array-like with at least 3 moments [E1, E2, E3].
        acf_lag1: Lag-1 autocorrelation (optional, default 0).

    Returns:
        dict with keys 'D0' and 'D1' (numpy arrays).
    """
    try:
        from line_solver.api.kpctoolbox.mmpp import mmpp2_fit
        moments = np.asarray(moments, dtype=np.float64).flatten()
        E1, E2, E3 = float(moments[0]), float(moments[1]), float(moments[2])
        acf = float(acf_lag1) if acf_lag1 is not None else 0.0
        D0, D1 = mmpp2_fit(E1, E2, E3, acf)
        return {'D0': np.asarray(D0), 'D1': np.asarray(D1)}
    except Exception as e:
        raise RuntimeError(f"MMPP2 fitting failed: {str(e)}")


def lib_kpc_mmpp2_fit_moments(moments):
    """Fit MMPP2 from moments only (no autocorrelation).

    Args:
        moments: Array-like with at least 3 moments [E1, E2, E3].

    Returns:
        dict with keys 'D0' and 'D1' (numpy arrays).
    """
    try:
        from line_solver.api.kpctoolbox.mmpp import mmpp2_fit
        moments = np.asarray(moments, dtype=np.float64).flatten()
        E1, E2, E3 = float(moments[0]), float(moments[1]), float(moments[2])
        D0, D1 = mmpp2_fit(E1, E2, E3, 0.0)
        return {'D0': np.asarray(D0), 'D1': np.asarray(D1)}
    except Exception as e:
        raise RuntimeError(f"MMPP2 fit_moments failed: {str(e)}")


def lib_kpc_mmpp2_fit_with_acf(moments, acf_lag1, acf_lag2=None):
    """Fit MMPP2 with autocorrelation information.

    If acf_lag2 is provided, uses mmpp2_fit3 with the G2 = acf_lag2/acf_lag1
    decay ratio. Otherwise uses mmpp2_fit with acf_lag1.

    Args:
        moments: Array-like with at least 3 moments [E1, E2, E3].
        acf_lag1: Lag-1 autocorrelation.
        acf_lag2: Lag-2 autocorrelation (optional).

    Returns:
        dict with keys 'D0' and 'D1' (numpy arrays).
    """
    try:
        moments = np.asarray(moments, dtype=np.float64).flatten()
        E1, E2, E3 = float(moments[0]), float(moments[1]), float(moments[2])

        if acf_lag2 is not None:
            from line_solver.api.kpctoolbox.mmpp import mmpp2_fit3
            # G2 = ratio of consecutive autocorrelations
            if abs(acf_lag1) < 1e-14:
                G2 = 0.0
            else:
                G2 = float(acf_lag2) / float(acf_lag1)
            D0, D1 = mmpp2_fit3(E1, E2, E3, G2)
        else:
            from line_solver.api.kpctoolbox.mmpp import mmpp2_fit
            D0, D1 = mmpp2_fit(E1, E2, E3, float(acf_lag1))

        return {'D0': np.asarray(D0), 'D1': np.asarray(D1)}
    except Exception as e:
        raise RuntimeError(f"MMPP2 fitting with ACF failed: {str(e)}")


# ========== MVPH (MULTIVARIATE PHASE-TYPE) ==========

def lib_kpc_mvph_mean_x(alpha, S, T, D):
    """Compute mean of first variable (X) in a bivariate PH distribution.

    Delegates to mvph_mean_x from api/kpctoolbox/mvph.

    Args:
        alpha: Initial probability vector (1D array).
        S: First phase generator matrix (for X).
        T: Second phase generator matrix (for Y).
        D: Transition matrix between phases.

    Returns:
        float: Mean of X.
    """
    try:
        from line_solver.api.kpctoolbox.mvph import mvph_mean_x
        return mvph_mean_x(alpha, S, T, D)
    except Exception as e:
        raise RuntimeError(f"mvph_mean_x failed: {str(e)}")


def lib_kpc_mvph_mean_y(alpha, S, T, D):
    """Compute mean of second variable (Y) in a bivariate PH distribution.

    Delegates to mvph_mean_y from api/kpctoolbox/mvph.

    Args:
        alpha: Initial probability vector (1D array).
        S: First phase generator matrix (for X).
        T: Second phase generator matrix (for Y).
        D: Transition matrix between phases.

    Returns:
        float: Mean of Y.
    """
    try:
        from line_solver.api.kpctoolbox.mvph import mvph_mean_y
        return mvph_mean_y(alpha, S, T, D)
    except Exception as e:
        raise RuntimeError(f"mvph_mean_y failed: {str(e)}")


def lib_kpc_mvph_cov(alpha, S, T, D):
    """Compute covariance between X and Y in a bivariate PH distribution.

    Delegates to mvph_cov from api/kpctoolbox/mvph.

    Args:
        alpha: Initial probability vector.
        S: First phase generator matrix (for X).
        T: Second phase generator matrix (for Y).
        D: Transition matrix between phases.

    Returns:
        float: Covariance of X and Y.
    """
    try:
        from line_solver.api.kpctoolbox.mvph import mvph_cov
        return mvph_cov(alpha, S, T, D)
    except Exception as e:
        raise RuntimeError(f"mvph_cov failed: {str(e)}")


def lib_kpc_mvph_corr(alpha, S, T, D):
    """Compute correlation between X and Y in a bivariate PH distribution.

    Delegates to mvph_corr from api/kpctoolbox/mvph.

    Args:
        alpha: Initial probability vector.
        S: First phase generator matrix (for X).
        T: Second phase generator matrix (for Y).
        D: Transition matrix between phases.

    Returns:
        float: Correlation of X and Y.
    """
    try:
        from line_solver.api.kpctoolbox.mvph import mvph_corr
        return mvph_corr(alpha, S, T, D)
    except Exception as e:
        raise RuntimeError(f"mvph_corr failed: {str(e)}")


def lib_kpc_mvph_joint(alpha, S, T, D, n1, n2):
    """Compute joint moment E[X^n1 * Y^n2] of a bivariate PH distribution.

    Delegates to mvph_joint from api/kpctoolbox/mvph.

    Args:
        alpha: Initial probability vector.
        S: First phase generator matrix (for X).
        T: Second phase generator matrix (for Y).
        D: Transition matrix between phases.
        n1: Power for first variable X (non-negative integer).
        n2: Power for second variable Y (non-negative integer).

    Returns:
        float: Joint moment E[X^n1 * Y^n2].
    """
    try:
        from line_solver.api.kpctoolbox.mvph import mvph_joint
        return mvph_joint(alpha, S, T, D, n1, n2)
    except Exception as e:
        raise RuntimeError(f"mvph_joint failed: {str(e)}")


# ========== MULTI-TRACE ANALYSIS ==========

def lib_kpc_mtrace_mean(traces):
    """Compute per-trace means for multiple traces.

    Args:
        traces: List of trace arrays.

    Returns:
        List of floats, one mean per trace.
    """
    try:
        return [float(np.mean(np.asarray(t, dtype=np.float64))) for t in traces]
    except Exception as e:
        raise RuntimeError(f"mtrace_mean failed: {str(e)}")


def lib_kpc_mtrace_var(traces):
    """Compute per-trace variances (unbiased) for multiple traces.

    Args:
        traces: List of trace arrays.

    Returns:
        List of floats, one variance per trace.
    """
    try:
        return [float(np.var(np.asarray(t, dtype=np.float64), ddof=1)) for t in traces]
    except Exception as e:
        raise RuntimeError(f"mtrace_var failed: {str(e)}")


def lib_kpc_mtrace_summary(traces, max_lag=10):
    """Compute summary statistics for each trace.

    Args:
        traces: List of trace arrays.
        max_lag: Maximum lag for autocorrelation (default: 10).

    Returns:
        dict with keys: means, vars, num_traces, summaries (list of per-trace dicts).
    """
    try:
        summaries = [lib_kpc_trace_summary(t, max_lag) for t in traces]
        return {
            'means': lib_kpc_mtrace_mean(traces),
            'vars': lib_kpc_mtrace_var(traces),
            'num_traces': len(traces),
            'summaries': summaries,
        }
    except Exception as e:
        raise RuntimeError(f"mtrace_summary failed: {str(e)}")


# ========== BASIC UTILITIES ==========

def lib_kpc_eye(n):
    """Create identity matrix.

    Args:
        n: Size of the identity matrix.

    Returns:
        numpy array of shape (n, n).
    """
    try:
        return np.eye(int(n))
    except Exception as e:
        raise RuntimeError(f"eye failed: {str(e)}")


def lib_kpc_ones(m, n):
    """Create matrix of all ones.

    Args:
        m: Number of rows.
        n: Number of columns.

    Returns:
        numpy array of shape (m, n).
    """
    try:
        return np.ones((int(m), int(n)))
    except Exception as e:
        raise RuntimeError(f"ones failed: {str(e)}")


def lib_kpc_zeros(m, n):
    """Create matrix of zeros.

    Args:
        m: Number of rows.
        n: Number of columns.

    Returns:
        numpy array of shape (m, n).
    """
    try:
        return np.zeros((int(m), int(n)))
    except Exception as e:
        raise RuntimeError(f"zeros failed: {str(e)}")


def lib_kpc_maxpos(data):
    """Find position of maximum value.

    Args:
        data: Array-like data.

    Returns:
        int: Index of maximum value.
    """
    try:
        data = np.asarray(data, dtype=np.float64).ravel()
        return int(np.argmax(data))
    except Exception as e:
        raise RuntimeError(f"maxpos failed: {str(e)}")


def lib_kpc_minpos(data):
    """Find position of minimum value.

    Args:
        data: Array-like data.

    Returns:
        int: Index of minimum value.
    """
    try:
        data = np.asarray(data, dtype=np.float64).ravel()
        return int(np.argmin(data))
    except Exception as e:
        raise RuntimeError(f"minpos failed: {str(e)}")


# ========== KPC FITTING ==========

def lib_kpc_fit_auto(trace, num_maps=None):
    """Automatically fit a PH distribution from trace data using KPC-Toolbox.

    Delegates to kpcfit_ph_auto from api/kpctoolbox/kpcfit.py.

    Args:
        trace: Array of inter-arrival times or moments.
        num_maps: Number of MAP components (optional).

    Returns:
        List of fitting results: (MAP, distance, params, method) tuples.
    """
    try:
        from line_solver.api.kpctoolbox.kpcfit import kpcfit_ph_auto, kpcfit_ph_options
        E = np.asarray(trace, dtype=np.float64).ravel()

        options = None
        if num_maps is not None:
            num_states = 2 ** int(num_maps)
            options = kpcfit_ph_options(E, max_num_states=num_states)

        return kpcfit_ph_auto(E, options)
    except Exception as e:
        raise RuntimeError(f"KPC auto fitting failed: {str(e)}")


def lib_kpc_fit_manual(trace, num_maps, structure='general'):
    """Manually fit KPC with specified number of components.

    Delegates to kpcfit_ph_auto with constrained options.

    Args:
        trace: Array of moments to fit.
        num_maps: Number of KPC components (log2 of states).
        structure: Structure type ('general', 'mmpp', 'ph', etc.).

    Returns:
        List of fitting results.
    """
    try:
        from line_solver.api.kpctoolbox.kpcfit import kpcfit_ph_auto, kpcfit_ph_options
        E = np.asarray(trace, dtype=np.float64).ravel()
        num_states = 2 ** int(num_maps)
        options = kpcfit_ph_options(
            E,
            min_num_states=num_states,
            max_num_states=num_states,
        )
        return kpcfit_ph_auto(E, options)
    except Exception as e:
        raise RuntimeError(f"KPC manual fitting failed: {str(e)}")


# ========== CTMC/DTMC MARKOV CHAIN METHODS ==========

def lib_kpc_ctmc_steady_state(Q):
    """Compute steady-state distribution for continuous-time Markov chain.

    Delegates to ctmc_solve from api/kpctoolbox/mc.py.

    Args:
        Q (array-like): Generator matrix (infinitesimal matrix).

    Returns:
        ndarray: Steady-state probability distribution.
    """
    try:
        from line_solver.api.kpctoolbox.mc import ctmc_solve
        Q = np.asarray(Q, dtype=np.float64)
        return ctmc_solve(Q)
    except Exception as e:
        raise RuntimeError(f"CTMC steady-state computation failed: {str(e)}")


def lib_kpc_dtmc_steady_state(P):
    """Compute steady-state distribution for discrete-time Markov chain.

    Delegates to dtmc_solve from api/kpctoolbox/mc.py.

    Args:
        P (array-like): Transition probability matrix.

    Returns:
        ndarray: Steady-state probability distribution.
    """
    try:
        from line_solver.api.kpctoolbox.mc import dtmc_solve
        P = np.asarray(P, dtype=np.float64)
        return dtmc_solve(P)
    except Exception as e:
        raise RuntimeError(f"DTMC steady-state computation failed: {str(e)}")


def lib_kpc_ctmc_transient(Q, pi0, t):
    """Compute transient distribution for CTMC at time t using uniformization.

    Delegates to ctmc_uniformization from api/kpctoolbox/mc.py.

    Args:
        Q (array-like): Generator matrix.
        pi0 (array-like): Initial probability distribution.
        t (float): Time at which to evaluate distribution.

    Returns:
        ndarray: Transient probability distribution at time t.
    """
    try:
        from line_solver.api.kpctoolbox.mc import ctmc_uniformization
        Q = np.asarray(Q, dtype=np.float64)
        pi0 = np.asarray(pi0, dtype=np.float64).ravel()
        result, _ = ctmc_uniformization(pi0, Q, float(t))
        return result
    except Exception as e:
        raise RuntimeError(f"CTMC transient computation failed: {str(e)}")


# ========== APH (ACYCLIC PHASE-TYPE) FUNCTIONS ==========

def lib_kpc_aph_from_moments(moments):
    """Fit acyclic phase-type distribution from moments.

    Delegates to aph_fit from api/kpctoolbox/aph.py.

    Args:
        moments (array-like): First few moments [E1, E2, E3].

    Returns:
        dict: APH parameters with keys 'alpha', 'A', 'D0', 'D1'.
    """
    try:
        from line_solver.api.kpctoolbox.aph import aph_fit
        moments = np.asarray(moments, dtype=np.float64).flatten()
        e1 = float(moments[0])
        e2 = float(moments[1])
        e3 = float(moments[2]) if len(moments) > 2 else e2 * 1.5

        MAP, is_exact = aph_fit(e1, e2, e3)
        D0 = MAP['D0']
        D1 = MAP['D1']

        # Extract alpha from D1: D1 = exit_rates * alpha^T (rank-1 for renewal)
        # For a renewal MAP: alpha = pie (stationary distribution at departures)
        # alpha = row of D1 normalized: D1[i,:] / sum(D1[i,:])
        n = D0.shape[0]
        row_sums = np.sum(D1, axis=1)
        nonzero_rows = np.where(row_sums > 1e-14)[0]
        if len(nonzero_rows) > 0:
            alpha = D1[nonzero_rows[0], :] / row_sums[nonzero_rows[0]]
        else:
            alpha = np.ones(n) / n

        return {
            'alpha': alpha,
            'A': D0,
            'D0': D0,
            'D1': D1,
        }
    except Exception as e:
        raise RuntimeError(f"APH fitting from moments failed: {str(e)}")


def lib_kpc_aph_bounds(moments):
    """Compute APH feasibility bounds for given moments.

    Checks whether the moment set [E1, E2, E3] can be represented
    by an acyclic phase-type distribution, and returns normalized
    moment bounds.

    Args:
        moments (array-like): First few moments [E1, E2, E3].

    Returns:
        dict: Feasibility information with keys:
            - 'n2': normalized second moment E2/(E1^2)
            - 'n3': normalized third moment E3/(E1*E2)
            - 'n2_feasible': bool, whether n2 is in APH range
            - 'n3_lb': lower bound on n3
            - 'n3_ub': upper bound on n3
            - 'min_order': minimum APH order needed
    """
    try:
        moments = np.asarray(moments, dtype=np.float64).flatten()
        e1 = float(moments[0])
        e2 = float(moments[1])
        e3 = float(moments[2]) if len(moments) > 2 else e2 * 1.5

        n2 = e2 / (e1 * e1)
        n3 = e3 / (e1 * e2)

        # Find minimum feasible order
        n2_feas = False
        n3_lb = float('inf')
        n3_ub = float('-inf')
        min_order = -1
        nmax = 100

        for n in range(2, nmax + 1):
            # Compute bounds for this order
            lb_n2 = (n + 1.0) / n
            ub_n2 = n / (n - 1.0) if n > 1 else float('inf')

            if n2 >= lb_n2:
                n2_feas = True

                # Lower bound on n3
                if n2 <= (n + 4.0) / (n + 1):
                    pn = ((n + 1) * (n2 - 2) / (3 * n2 * (n - 1))) * \
                         (-2 * np.sqrt(n + 1.0) / np.sqrt(max(1e-14, 4.0 * (n + 1) - 3 * n * n2)) - 1)
                    an = (n2 - 2) / (pn * (1 - n2) + np.sqrt(max(1e-14, pn * pn + pn * n * (n2 - 2) / (n - 1))))
                    ln = ((3 + an) * (n - 1) + 2 * an) / ((n - 1) * (1 + an * pn)) - \
                         (2 * an * (n + 1)) / (2 * (n - 1) + an * pn * (n * an + 2 * n - 2))
                    n3_lb = ln
                else:
                    n3_lb = n2 * (n + 1) / n

                # Upper bound on n3
                if n2 <= ub_n2:
                    un = (1.0 / (n * n * n2)) * (2 * (n - 2) * (n * n2 - n - 1) *
                         np.sqrt(max(1e-14, 1 + n * (n2 - 2) / (n - 1))) +
                         (n + 2) * (3 * n * n2 - 2 * n - 2))
                    n3_ub = un
                else:
                    n3_ub = float('inf')

                if n3 >= n3_lb and n3 <= n3_ub:
                    min_order = n
                    break

        return {
            'n2': n2,
            'n3': n3,
            'n2_feasible': n2_feas,
            'n3_lb': n3_lb,
            'n3_ub': n3_ub,
            'min_order': min_order,
        }
    except Exception as e:
        raise RuntimeError(f"APH bounds computation failed: {str(e)}")
