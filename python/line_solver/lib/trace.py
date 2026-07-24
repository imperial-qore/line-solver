"""
Low-level trace analysis utilities.

Provides direct access to trace processing functions for:
- Autocorrelation and cross-covariance computation
- Index of dispersion calculations
- Moment and lag computation
- Trace statistics and summaries
- Multiple trace processing

These are complementary to api.trace functions, providing more granular control
over trace analysis.
"""

import numpy as np
from line_solver import native_to_array


def lib_trace_acf(trace, max_lag):
    """
    Compute autocorrelation function of a trace.

    Args:
        trace (array-like): Inter-arrival times or other sequence
        max_lag (int): Maximum lag to compute

    Returns:
        ndarray: Autocorrelation values for lags 1 to max_lag
    """
    trace = np.asarray(trace)
    trace_arr = native_to_array(trace)

    try:
            trace_arr, max_lag
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"trace_acf failed: {str(e)}")


def lib_trace_bicov(trace1, trace2, max_lag):
    """
    Compute bivariate covariance between two traces.

    Args:
        trace1 (array-like): First trace
        trace2 (array-like): Second trace
        max_lag (int): Maximum lag

    Returns:
        ndarray: Cross-covariance values
    """
    trace1 = np.asarray(trace1)
    trace2 = np.asarray(trace2)
    trace1_arr = native_to_array(trace1)
    trace2_arr = native_to_array(trace2)

    try:
            trace1_arr, trace2_arr, max_lag
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"trace_bicov failed: {str(e)}")


def lib_trace_gamma(trace, max_lag):
    """
    Compute lag gamma (correlation coefficient) of a trace.

    Args:
        trace (array-like): Trace data
        max_lag (int): Maximum lag

    Returns:
        ndarray: Gamma values for each lag
    """
    trace = np.asarray(trace)
    trace_arr = native_to_array(trace)

    try:
            trace_arr, max_lag
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"trace_gamma failed: {str(e)}")


def lib_trace_iat2bins(trace, num_bins):
    """
    Bin inter-arrival times from a trace.

    Args:
        trace (array-like): Trace of arrival times
        num_bins (int): Number of bins

    Returns:
        ndarray: Histogram of inter-arrival times
    """
    trace = np.asarray(trace)
    trace_arr = native_to_array(trace)

    try:
            trace_arr, num_bins
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"trace_iat2bins failed: {str(e)}")


def lib_trace_iat2counts(trace, num_bins):
    """
    Convert inter-arrival times to count histogram.

    Args:
        trace (array-like): Trace data
        num_bins (int): Number of bins

    Returns:
        ndarray: Counts in each time bin
    """
    trace = np.asarray(trace)
    trace_arr = native_to_array(trace)

    try:
            trace_arr, num_bins
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"trace_iat2counts failed: {str(e)}")


def lib_trace_idc(trace, time_window):
    """
    Compute index of dispersion for counts (IDC).

    Args:
        trace (array-like): Trace data
        time_window (float): Time window for analysis

    Returns:
        float: Index of dispersion for counts
    """
    trace = np.asarray(trace)
    trace_arr = native_to_array(trace)

    try:
            trace_arr, float(time_window)
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"trace_idc failed: {str(e)}")


def lib_trace_idi(trace, num_intervals):
    """
    Compute index of dispersion for intervals (IDI).

    Args:
        trace (array-like): Trace data
        num_intervals (int): Number of intervals

    Returns:
        float: Index of dispersion for intervals
    """
    trace = np.asarray(trace)
    trace_arr = native_to_array(trace)

    try:
            trace_arr, num_intervals
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"trace_idi failed: {str(e)}")


def lib_trace_mean(trace):
    """
    Compute mean of trace values.

    Args:
        trace (array-like): Trace data

    Returns:
        float: Mean value
    """
    trace = np.asarray(trace)
    trace_arr = native_to_array(trace)

    try:
            trace_arr
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"trace_mean failed: {str(e)}")


def lib_trace_scv(trace):
    """
    Compute squared coefficient of variation (SCV) of trace.

    Args:
        trace (array-like): Trace data

    Returns:
        float: Squared coefficient of variation
    """
    trace = np.asarray(trace)
    trace_arr = native_to_array(trace)

    try:
            trace_arr
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"trace_scv failed: {str(e)}")


def lib_trace_var(trace):
    """
    Compute variance of trace values.

    Args:
        trace (array-like): Trace data

    Returns:
        float: Variance
    """
    trace = np.asarray(trace)
    trace_arr = native_to_array(trace)

    try:
            trace_arr
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"trace_var failed: {str(e)}")


def lib_trace_skew(trace):
    """
    Compute skewness of trace values.

    Args:
        trace (array-like): Trace data

    Returns:
        float: Skewness coefficient
    """
    trace = np.asarray(trace)
    trace_arr = native_to_array(trace)

    try:
            trace_arr
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"trace_skew failed: {str(e)}")


def lib_trace_summary(trace):
    """
    Compute comprehensive summary statistics of a trace.

    Args:
        trace (array-like): Trace data

    Returns:
        dict: Dictionary with keys: mean, var, skew, scv, idc, idi
    """
    trace = np.asarray(trace)
    mean = lib_trace_mean(trace)
    var = lib_trace_var(trace)
    skew = lib_trace_skew(trace)
    scv = lib_trace_scv(trace)
    idc = lib_trace_idc(trace, 1.0)
    idi = lib_trace_idi(trace, 100)

    return {
        'mean': mean,
        'var': var,
        'skew': skew,
        'scv': scv,
        'idc': idc,
        'idi': idi
    }


def lib_mtrace_mean(traces):
    """
    Compute mean across multiple traces.

    Args:
        traces (list of array-like): Multiple traces

    Returns:
        float: Mean across all traces
    """
    traces = [np.asarray(t) for t in traces]
    traces_arr = [native_to_array(t) for t in traces]

    try:
            traces_arr
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"mtrace_mean failed: {str(e)}")


def lib_mtrace_summary(traces):
    """
    Compute summary statistics across multiple traces.

    Args:
        traces (list of array-like): Multiple traces

    Returns:
        dict: Summary statistics with confidence intervals
    """
    mean = lib_mtrace_mean(traces)
    # Additional statistics would be computed here
    return {
        'mean': mean,
        'traces_count': len(traces)
    }
