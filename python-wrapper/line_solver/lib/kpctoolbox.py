"""
KPCToolbox library functions - comprehensive utilities:
- Trace analysis (single and multi-trace statistics)
- MMPP (Markov Modulated Poisson Process) fitting
- KPC fitting (Kronecker Product Composition)
- MVPH (Multivariate Phase-Type) analysis
- CTMC/DTMC utilities
- Advanced probability computations
"""

import jpype
import numpy as np
from line_solver import jlineMatrixFromArray, jlineMatrixToArray


# ========== TRACE ANALYSIS (SINGLE) ==========

def lib_kpc_trace_mean(trace):
    """Compute mean of trace."""
    trace = np.asarray(trace)
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.trace_mean(
            jlineMatrixFromArray(trace)
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"trace_mean failed: {str(e)}")


def lib_kpc_trace_var(trace):
    """Compute variance of trace."""
    trace = np.asarray(trace)
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.trace_var(
            jlineMatrixFromArray(trace)
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"trace_var failed: {str(e)}")


def lib_kpc_trace_scv(trace):
    """Compute squared coefficient of variation."""
    trace = np.asarray(trace)
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.trace_scv(
            jlineMatrixFromArray(trace)
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"trace_scv failed: {str(e)}")


def lib_kpc_trace_acf(trace, max_lag):
    """Compute autocorrelation function."""
    trace = np.asarray(trace)
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.trace_acf(
            jlineMatrixFromArray(trace), int(max_lag)
        )
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"trace_acf failed: {str(e)}")


def lib_kpc_trace_idc(trace, time_window):
    """Compute index of dispersion for counts."""
    trace = np.asarray(trace)
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.trace_idc(
            jlineMatrixFromArray(trace), float(time_window)
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"trace_idc failed: {str(e)}")


def lib_kpc_trace_iat2counts(trace, num_bins):
    """Convert inter-arrival times to count histogram."""
    trace = np.asarray(trace)
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.trace_iat2counts(
            jlineMatrixFromArray(trace), int(num_bins)
        )
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"trace_iat2counts failed: {str(e)}")


def lib_kpc_trace_summary(trace):
    """Compute comprehensive trace summary."""
    trace = np.asarray(trace)
    return {
        'mean': lib_kpc_trace_mean(trace),
        'var': lib_kpc_trace_var(trace),
        'scv': lib_kpc_trace_scv(trace),
        'acf_lag1': lib_kpc_trace_acf(trace, 1)[0] if lib_kpc_trace_acf(trace, 1).size > 0 else 0
    }


# ========== MMPP FITTING ==========

def lib_kpc_mmpp2_fit(moments, acf_lag1=None):
    """Fit 2-state MMPP."""
    moments = np.asarray(moments).flatten()
    try:
        MMPP = jpype.JPackage('jline').lib.kpctoolbox.MMPP2
        if acf_lag1 is not None:
            result = MMPP.fit(jlineMatrixFromArray(moments), float(acf_lag1))
        else:
            result = MMPP.fit(jlineMatrixFromArray(moments))
        return {'parameters': jlineMatrixToArray(result) if hasattr(result, 'toArray') else result}
    except Exception as e:
        raise RuntimeError(f"MMPP2 fitting failed: {str(e)}")


def lib_kpc_mmpp2_fit_moments(moments):
    """Fit MMPP2 from moments only."""
    moments = np.asarray(moments).flatten()
    try:
        MMPP = jpype.JPackage('jline').lib.kpctoolbox.MMPP2
        result = MMPP.fit1(jlineMatrixFromArray(moments))
        return {'parameters': jlineMatrixToArray(result) if hasattr(result, 'toArray') else result}
    except Exception as e:
        raise RuntimeError(f"MMPP2 fit1 failed: {str(e)}")


def lib_kpc_mmpp2_fit_with_acf(moments, acf_lag1, acf_lag2=None):
    """Fit MMPP2 with autocorrelation."""
    moments = np.asarray(moments).flatten()
    try:
        MMPP = jpype.JPackage('jline').lib.kpctoolbox.MMPP2
        if acf_lag2 is None:
            result = MMPP.fit2(jlineMatrixFromArray(moments), float(acf_lag1))
        else:
            result = MMPP.fit3(jlineMatrixFromArray(moments), float(acf_lag1), float(acf_lag2))
        return {'parameters': jlineMatrixToArray(result) if hasattr(result, 'toArray') else result}
    except Exception as e:
        raise RuntimeError(f"MMPP2 fitting failed: {str(e)}")


# ========== MVPH (MULTIVARIATE PHASE-TYPE) ==========

def lib_kpc_mvph_mean_x(alpha, A, B):
    """Compute mean of first output (X)."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    B = np.asarray(B)
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.mvph_mean_x(
            jlineMatrixFromArray(alpha), jlineMatrixFromArray(A), jlineMatrixFromArray(B)
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"mvph_mean_x failed: {str(e)}")


def lib_kpc_mvph_mean_y(alpha, A, C):
    """Compute mean of second output (Y)."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    C = np.asarray(C)
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.mvph_mean_y(
            jlineMatrixFromArray(alpha), jlineMatrixFromArray(A), jlineMatrixFromArray(C)
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"mvph_mean_y failed: {str(e)}")


def lib_kpc_mvph_cov(alpha, A, B, C):
    """Compute covariance between outputs."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    B = np.asarray(B)
    C = np.asarray(C)
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.mvph_cov(
            jlineMatrixFromArray(alpha), jlineMatrixFromArray(A),
            jlineMatrixFromArray(B), jlineMatrixFromArray(C)
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"mvph_cov failed: {str(e)}")


def lib_kpc_mvph_corr(alpha, A, B, C):
    """Compute correlation between outputs."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    B = np.asarray(B)
    C = np.asarray(C)
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.mvph_corr(
            jlineMatrixFromArray(alpha), jlineMatrixFromArray(A),
            jlineMatrixFromArray(B), jlineMatrixFromArray(C)
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"mvph_corr failed: {str(e)}")


def lib_kpc_mvph_joint(alpha, A, B, C, x_vals, y_vals):
    """Compute joint distribution of outputs."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    B = np.asarray(B)
    C = np.asarray(C)
    x_vals = np.asarray(x_vals)
    y_vals = np.asarray(y_vals)
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.mvph_joint(
            jlineMatrixFromArray(alpha), jlineMatrixFromArray(A),
            jlineMatrixFromArray(B), jlineMatrixFromArray(C),
            jlineMatrixFromArray(x_vals), jlineMatrixFromArray(y_vals)
        )
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"mvph_joint failed: {str(e)}")


# ========== MULTI-TRACE ANALYSIS ==========

def lib_kpc_mtrace_mean(traces):
    """Compute mean across multiple traces."""
    try:
        traces_jline = [jlineMatrixFromArray(np.asarray(t)) for t in traces]
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.mtrace_mean(traces_jline)
        return float(result)
    except Exception as e:
        raise RuntimeError(f"mtrace_mean failed: {str(e)}")


def lib_kpc_mtrace_var(traces):
    """Compute variance across multiple traces."""
    try:
        traces_jline = [jlineMatrixFromArray(np.asarray(t)) for t in traces]
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.mtrace_var(traces_jline)
        return float(result)
    except Exception as e:
        raise RuntimeError(f"mtrace_var failed: {str(e)}")


def lib_kpc_mtrace_summary(traces):
    """Compute summary statistics across multiple traces."""
    return {
        'mean': lib_kpc_mtrace_mean(traces),
        'var': lib_kpc_mtrace_var(traces),
        'num_traces': len(traces)
    }


# ========== BASIC UTILITIES ==========

def lib_kpc_eye(n):
    """Create identity matrix."""
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.eye(int(n))
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"eye failed: {str(e)}")


def lib_kpc_ones(m, n):
    """Create matrix of all ones."""
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.ones(int(m), int(n))
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"ones failed: {str(e)}")


def lib_kpc_zeros(m, n):
    """Create matrix of zeros."""
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.zeros(int(m), int(n))
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"zeros failed: {str(e)}")


def lib_kpc_maxpos(data):
    """Find position of maximum value."""
    data = np.asarray(data)
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.maxpos(
            jlineMatrixFromArray(data)
        )
        return int(result)
    except Exception as e:
        raise RuntimeError(f"maxpos failed: {str(e)}")


def lib_kpc_minpos(data):
    """Find position of minimum value."""
    data = np.asarray(data)
    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.minpos(
            jlineMatrixFromArray(data)
        )
        return int(result)
    except Exception as e:
        raise RuntimeError(f"minpos failed: {str(e)}")


# ========== KPC FITTING ==========

def lib_kpc_fit_auto(moments, acf_lag1=None, acf_lag2=None):
    """Automatically fit KPC distribution from moments and autocorrelation.

    Args:
        moments (array-like): First few moments
        acf_lag1 (float): Lag-1 autocorrelation (optional)
        acf_lag2 (float): Lag-2 autocorrelation (optional)

    Returns:
        dict: KPC fitting results
    """
    moments = np.asarray(moments).flatten()
    try:
        KPCFit = jpype.JPackage('jline').lib.kpctoolbox.KPCFit
        if acf_lag2 is not None:
            result = KPCFit.fitAuto(
                jlineMatrixFromArray(moments),
                float(acf_lag1), float(acf_lag2)
            )
        elif acf_lag1 is not None:
            result = KPCFit.fitAuto(
                jlineMatrixFromArray(moments),
                float(acf_lag1)
            )
        else:
            result = KPCFit.fitAuto(jlineMatrixFromArray(moments))
        return jlineMatrixToArray(result) if hasattr(result, 'toArray') else result
    except Exception as e:
        raise RuntimeError(f"KPC auto fitting failed: {str(e)}")


def lib_kpc_fit_manual(moments, n_components, structure='general'):
    """Manually fit KPC with specified number of components.

    Args:
        moments (array-like): Moments to fit
        n_components (int): Number of KPC components
        structure (str): Structure type ('general', 'mmpp', 'ph', etc.)

    Returns:
        dict: KPC parameters
    """
    moments = np.asarray(moments).flatten()
    try:
        KPCFit = jpype.JPackage('jline').lib.kpctoolbox.KPCFit
        result = KPCFit.fitManual(
            jlineMatrixFromArray(moments),
            int(n_components),
            str(structure)
        )
        return jlineMatrixToArray(result) if hasattr(result, 'toArray') else result
    except Exception as e:
        raise RuntimeError(f"KPC manual fitting failed: {str(e)}")


# ========== CTMC/DTMC MARKOV CHAIN METHODS ==========

def lib_kpc_ctmc_steady_state(Q):
    """Compute steady-state distribution for continuous-time Markov chain.

    Args:
        Q (array-like): Generator matrix (infinitesimal matrix)

    Returns:
        ndarray: Steady-state probability distribution
    """
    Q = np.asarray(Q)
    try:
        MC = jpype.JPackage('jline').lib.kpctoolbox.MarkovChain
        result = MC.ctmcSteadyState(jlineMatrixFromArray(Q))
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"CTMC steady-state computation failed: {str(e)}")


def lib_kpc_dtmc_steady_state(P):
    """Compute steady-state distribution for discrete-time Markov chain.

    Args:
        P (array-like): Transition probability matrix

    Returns:
        ndarray: Steady-state probability distribution
    """
    P = np.asarray(P)
    try:
        MC = jpype.JPackage('jline').lib.kpctoolbox.MarkovChain
        result = MC.dtmcSteadyState(jlineMatrixFromArray(P))
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"DTMC steady-state computation failed: {str(e)}")


def lib_kpc_ctmc_transient(Q, t):
    """Compute transient distribution for CTMC at time t.

    Args:
        Q (array-like): Generator matrix
        t (float): Time at which to evaluate distribution

    Returns:
        ndarray: Transient probability distribution
    """
    Q = np.asarray(Q)
    try:
        MC = jpype.JPackage('jline').lib.kpctoolbox.MarkovChain
        result = MC.ctmcTransient(jlineMatrixFromArray(Q), float(t))
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"CTMC transient computation failed: {str(e)}")


# ========== APH (ACYCLIC PHASE-TYPE) FUNCTIONS ==========

def lib_kpc_aph_from_moments(moments):
    """Fit acyclic phase-type distribution from moments.

    Args:
        moments (array-like): First few moments

    Returns:
        dict: APH parameters {'alpha': alpha, 'A': A}
    """
    moments = np.asarray(moments).flatten()
    try:
        APH = jpype.JPackage('jline').lib.kpctoolbox.APH
        result = APH.fromMoments(jlineMatrixFromArray(moments))
        alpha = jlineMatrixToArray(result.first)
        A = jlineMatrixToArray(result.second)
        return {'alpha': alpha, 'A': A}
    except Exception as e:
        raise RuntimeError(f"APH fitting from moments failed: {str(e)}")


def lib_kpc_aph_bounds(m1, m2, m3, order=2):
    """Compute moment bounds for acyclic phase-type distributions.

    Args:
        m1 (float): First moment (mean)
        m2 (float): Second moment
        m3 (float): Third moment
        order (int): Phase order (2 or 3)

    Returns:
        dict: Moment bounds for higher orders
    """
    try:
        APH = jpype.JPackage('jline').lib.kpctoolbox.APH
        if order == 2:
            result = APH.bounds2ndOrder(float(m1), float(m2), float(m3))
        elif order == 3:
            result = APH.bounds3rdOrder(float(m1), float(m2), float(m3))
        else:
            raise ValueError(f"Order {order} not supported")
        return jlineMatrixToArray(result) if hasattr(result, 'toArray') else float(result)
    except Exception as e:
        raise RuntimeError(f"APH bounds computation failed: {str(e)}")
