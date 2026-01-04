
"""
Trace data analysis functions.

This module provides functions for analyzing empirical trace data,
including statistical measures and transformations of time series
and event traces.

These functions support trace-driven modeling and empirical analysis
of system measurements and workload characterization.
"""

import jpype
import numpy as np
from line_solver import jlineMatrixToArray, jlineMatrixFromArray


def trace_mean(trace):
    """
    Calculate the mean of a trace.

    Args:
        trace: Array of trace values

    Returns:
        float: Mean value of the trace
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)

    result = jpype.JPackage('jline').api.trace.Trace_meanKt.trace_mean(java_trace)
    return float(result)


def trace_var(trace):
    """
    Calculate the variance of a trace.

    Args:
        trace: Array of trace values

    Returns:
        float: Variance of the trace
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)

    result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_var(java_trace)
    return float(result)


def mtrace_mean(trace, ntypes, type_array):
    """
    Calculate mean of multi-type trace data.

    Computes the mean for each type in a multi-type trace dataset.

    Args:
        trace: Array of trace values
        ntypes: Number of different types in the trace
        type_array: Array indicating the type of each trace element

    Returns:
        numpy.ndarray: Mean values for each type
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)
    java_type_array = jpype.JArray(jpype.JInt)(type_array)

    result = jpype.JPackage('jline').api.trace.Mtrace_meanKt.mtrace_mean(
        java_trace, jpype.JInt(ntypes), java_type_array
    )

    return jlineMatrixToArray(result)


def trace_scv(trace):
    """
    Calculate the squared coefficient of variation of a trace.

    Args:
        trace: Array of trace values

    Returns:
        float: Squared coefficient of variation (variance/mean²)
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)
    result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_scv(java_trace)
    return float(result)


def trace_acf(trace, lags=None):
    """
    Calculate the autocorrelation function of a trace.

    Args:
        trace: Array of trace values
        lags: List of lag values to compute autocorrelation for (default: [1])

    Returns:
        numpy.ndarray: Autocorrelation values at specified lags
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)

    if lags is None:
        lags = [1]
    java_lags = jpype.JArray(jpype.JInt)(lags)

    result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_acf(java_trace, java_lags)
    return np.array([float(x) for x in result])


def trace_gamma(trace, limit=1000):
    """
    Fit gamma distribution parameters to trace data.

    Args:
        trace: Array of trace values
        limit: Maximum number of iterations for fitting (default: 1000)

    Returns:
        tuple: (shape, scale, rate) parameters of fitted gamma distribution
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)
    result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_gamma(java_trace, jpype.JInt(limit))
    return (float(result[0]), float(result[1]), float(result[2]))


def trace_iat2counts(trace, scale):
    """
    Convert inter-arrival times to count data.

    Args:
        trace: Array of inter-arrival times
        scale: Time scale for binning

    Returns:
        numpy.ndarray: Count values in each time bin
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)
    result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_iat2counts(java_trace, jpype.JDouble(scale))
    return np.array([int(x) for x in result])


def trace_idi(trace, kset, option=None, n=1):
    """
    Calculate index of dispersion for intervals.

    Measures the variability of intervals in the trace data.

    Args:
        trace: Array of trace values
        kset: Array of interval lengths to analyze
        option: Analysis option (optional)
        n: Parameter for analysis method (default: 1)

    Returns:
        tuple: (idi_values, support_values) - IDI values and their support
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)
    java_kset = jpype.JArray(jpype.JInt)(kset)

    if option is None:
        result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_idi(java_trace, java_kset)
    else:
        result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_idi(
            java_trace, java_kset, jpype.JString(option), jpype.JInt(n)
        )

    idi_values = np.array([float(x) for x in result.getFirst()])
    support_values = np.array([int(x) for x in result.getSecond()])
    return (idi_values, support_values)


def trace_idc(trace):
    """
    Calculate index of dispersion for counts.

    Args:
        trace: Array of trace values

    Returns:
        float: Index of dispersion for counts
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)
    result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_idc(java_trace)
    return float(result)


def trace_pmf(data):
    """
    Calculate probability mass function of discrete data.

    Args:
        data: Array of discrete data values

    Returns:
        tuple: (pmf_values, unique_values) - PMF and corresponding unique values
    """
    java_data = jpype.JArray(jpype.JInt)([int(x) for x in data])
    result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_pmf(java_data)

    pmf_values = np.array([float(x) for x in result.getFirst()])
    unique_values = np.array([int(x) for x in result.getSecond()])
    return (pmf_values, unique_values)


def trace_shuffle(trace):
    """
    Randomly shuffle trace data.

    Args:
        trace: Array of trace values

    Returns:
        numpy.ndarray: Shuffled trace data
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)
    result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_shuffle(java_trace)
    return np.array([float(x) for x in result])


def trace_bicov(trace, grid):
    """
    Calculate bicovariance of trace data.

    Args:
        trace: Array of trace values
        grid: Grid of lag values for bicovariance calculation

    Returns:
        tuple: (bicov_values, lag_combinations) - Bicovariance values and lag pairs
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)
    java_grid = jpype.JArray(jpype.JInt)(grid)

    result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_bicov(java_trace, java_grid)

    bicov_values = np.array([float(x) for x in result.getFirst()])
    lag_combinations = []
    for lag_array in result.getSecond():
        lag_combinations.append([int(x) for x in lag_array])

    return (bicov_values, lag_combinations)


def trace_iat2bins(trace, scale):
    """
    Convert inter-arrival times to histogram bins.

    Args:
        trace: Array of inter-arrival times
        scale: Time scale for binning

    Returns:
        tuple: (counts, bin_membership) - Bin counts and membership array
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)
    result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_iat2bins(java_trace, jpype.JDouble(scale))

    counts = np.array([int(x) for x in result.getFirst()])
    bin_membership = np.array([int(x) for x in result.getSecond()])
    return (counts, bin_membership)


def trace_joint(trace, lag, order):
    """
    Calculate joint statistics of trace data.

    Args:
        trace: Array of trace values
        lag: Array of lag values
        order: Array of moment orders

    Returns:
        float: Joint statistic value
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)
    java_lag = jpype.JArray(jpype.JInt)(lag)
    java_order = jpype.JArray(jpype.JInt)(order)

    result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_joint(java_trace, java_lag, java_order)
    return float(result)


def trace_summary(trace):
    """
    Calculate comprehensive summary statistics for trace data.

    Args:
        trace: Array of trace values

    Returns:
        dict: Dictionary containing various statistics:
            - mean, scv, mad, skew, kurt: Basic statistics
            - q25, q50, q75, p95: Quantiles and percentiles
            - min, max, iqr: Range statistics
            - acf1-acf4: Autocorrelation at lags 1-4
            - idc_scv_ratio: IDC to SCV ratio
    """
    java_trace = jpype.JArray(jpype.JDouble)(trace)
    result = jpype.JPackage('jline').api.trace.Trace_varKt.trace_summary(java_trace)

    result_array = [float(x) for x in result]

    return {
        'mean': result_array[0],
        'scv': result_array[1],
        'mad': result_array[2],
        'skew': result_array[3],
        'kurt': result_array[4],
        'q25': result_array[5],
        'q50': result_array[6],
        'q75': result_array[7],
        'p95': result_array[8],
        'min': result_array[9],
        'max': result_array[10],
        'iqr': result_array[11],
        'acf1': result_array[12],
        'acf2': result_array[13],
        'acf3': result_array[14],
        'acf4': result_array[15],
        'idc_scv_ratio': result_array[16]
    }


def trace_cdf(trace, x_values=None):
    """
    Calculate cumulative distribution function of trace data.

    Args:
        trace: Array of trace values
        x_values: X values to evaluate CDF at (optional, defaults to unique sorted trace values)

    Returns:
        dict: Dictionary with 'x' and 'y' arrays for CDF plot
    """
    import numpy as np

    trace = np.array(trace).flatten()

    if x_values is None:
        x_values = np.unique(np.sort(trace))
    else:
        x_values = np.array(x_values)

    n = len(trace)
    cdf_values = np.zeros_like(x_values, dtype=float)

    for i, x in enumerate(x_values):
        cdf_values[i] = np.sum(trace <= x) / n

    return {
        'x': x_values,
        'y': cdf_values
    }


def trace_pdf(trace, x_values=None, bandwidth=None):
    """
    Calculate probability density function using kernel density estimation.

    Args:
        trace: Array of trace values
        x_values: X values to evaluate PDF at (optional)
        bandwidth: KDE bandwidth (optional, uses default if None)

    Returns:
        dict: Dictionary with 'x' and 'y' arrays for PDF plot
    """
    import numpy as np
    from scipy.stats import gaussian_kde

    trace = np.array(trace).flatten()

    trace = trace[np.isfinite(trace)]

    if len(trace) == 0:
        raise ValueError("No finite values in trace data")

    if x_values is None:
        x_min, x_max = np.min(trace), np.max(trace)
        x_range = x_max - x_min
        x_values = np.linspace(x_min - 0.1*x_range, x_max + 0.1*x_range, 100)
    else:
        x_values = np.array(x_values)

    kde = gaussian_kde(trace)
    if bandwidth is not None:
        kde.set_bandwidth(bandwidth)

    pdf_values = kde(x_values)

    return {
        'x': x_values,
        'y': pdf_values
    }


def trace_hist(trace, bins=None):
    """
    Calculate histogram of trace data.

    Args:
        trace: Array of trace values
        bins: Number of bins or bin specification (optional, defaults to 'auto')

    Returns:
        dict: Dictionary with 'counts', 'bin_edges', and 'bin_centers'
    """
    import numpy as np

    trace = np.array(trace).flatten()
    trace = trace[np.isfinite(trace)]

    if bins is None:
        bins = 'auto'

    counts, bin_edges = np.histogram(trace, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    return {
        'counts': counts,
        'bin_edges': bin_edges,
        'bin_centers': bin_centers
    }


def trace_moment(trace, k):
    """
    Calculate the k-th moment of trace data.

    Args:
        trace: Array of trace values
        k: Moment order

    Returns:
        float: k-th moment value
    """
    import numpy as np

    trace = np.array(trace).flatten()
    trace = trace[np.isfinite(trace)]

    if len(trace) == 0:
        return 0.0

    return float(np.mean(trace ** k))


def trace_skew(trace):
    """
    Calculate skewness of trace data.

    Args:
        trace: Array of trace values

    Returns:
        float: Skewness value
    """
    import numpy as np
    from scipy.stats import skew

    trace = np.array(trace).flatten()
    trace = trace[np.isfinite(trace)]

    if len(trace) < 3:
        return 0.0

    return float(skew(trace))


def trace_kurt(trace):
    """
    Calculate kurtosis of trace data.

    Args:
        trace: Array of trace values

    Returns:
        float: Excess kurtosis value (Fisher's definition)
    """
    import numpy as np
    from scipy.stats import kurtosis

    trace = np.array(trace).flatten()
    trace = trace[np.isfinite(trace)]

    if len(trace) < 4:
        return 0.0

    return float(kurtosis(trace, fisher=True))


def trace_fit_gamma(trace):
    """
    Fit gamma distribution to trace data using method of moments.

    Args:
        trace: Array of positive trace values

    Returns:
        dict: Gamma distribution parameters:
            - 'shape': Shape parameter (alpha)
            - 'scale': Scale parameter (beta)
            - 'rate': Rate parameter (1/beta)
            - 'mean': Sample mean
            - 'var': Sample variance

    Raises:
        ValueError: If no positive finite values or zero variance
    """
    import numpy as np

    trace = np.array(trace).flatten()
    trace = trace[np.isfinite(trace) & (trace > 0)]

    if len(trace) == 0:
        raise ValueError("No positive finite values in trace data")

    sample_mean = np.mean(trace)
    sample_var = np.var(trace, ddof=1)

    if sample_var <= 0:
        raise ValueError("Sample variance must be positive for gamma fitting")

    scale = sample_var / sample_mean
    shape = sample_mean / scale
    rate = 1.0 / scale

    return {
        'shape': float(shape),
        'scale': float(scale),
        'rate': float(rate),
        'mean': float(sample_mean),
        'var': float(sample_var)
    }


def mtrace_corr(trace, ntypes, type_array, lags=None):
    """
    Calculate cross-correlation matrix for multi-type trace data.

    Args:
        trace: Array of trace values
        ntypes: Number of different types
        type_array: Array indicating the type of each trace element
        lags: Array of lag values (default: [0,1,2,3,4])

    Returns:
        dict: Dictionary with:
            - 'correlations': 3D array (ntypes x ntypes x nlags)
            - 'lags': Array of lag values used

    Raises:
        ValueError: If trace and type_array have different lengths
    """
    import numpy as np

    trace = np.array(trace).flatten()
    type_array = np.array(type_array).flatten()

    if len(trace) != len(type_array):
        raise ValueError("trace and type_array must have same length")

    if lags is None:
        lags = [0, 1, 2, 3, 4]
    lags = np.array(lags)

    nlags = len(lags)
    correlations = np.zeros((ntypes, ntypes, nlags))

    for i in range(ntypes):
        for j in range(ntypes):
            for k, lag in enumerate(lags):
                if lag == 0:
                    mask_i = (type_array == i)
                    mask_j = (type_array == j)
                    if np.sum(mask_i) > 0 and np.sum(mask_j) > 0:
                        if i == j:
                            correlations[i, j, k] = 1.0
                        else:
                            correlations[i, j, k] = np.corrcoef(
                                trace[mask_i], trace[mask_j]
                            )[0, 1] if len(trace[mask_i]) > 1 and len(trace[mask_j]) > 1 else 0.0
                else:
                    if len(trace) > lag:
                        trace1 = trace[:-lag]
                        trace2 = trace[lag:]
                        type1 = type_array[:-lag]
                        type2 = type_array[lag:]

                        mask1 = (type1 == i)
                        mask2 = (type2 == j)

                        if np.sum(mask1) > 0 and np.sum(mask2) > 0:
                            correlations[i, j, k] = np.corrcoef(
                                trace1[mask1], trace2[mask2]
                            )[0, 1] if len(trace1[mask1]) > 1 and len(trace2[mask2]) > 1 else 0.0

    return {
        'correlations': correlations,
        'lags': lags
    }


def mtrace_cov(trace, ntypes, type_array, lags=None):
    """
    Calculate cross-covariance matrix for multi-type trace data.

    Args:
        trace: Array of trace values
        ntypes: Number of different types
        type_array: Array indicating the type of each trace element
        lags: Array of lag values (default: [0,1,2,3,4])

    Returns:
        dict: Dictionary with:
            - 'covariances': 3D array (ntypes x ntypes x nlags)
            - 'lags': Array of lag values used

    Raises:
        ValueError: If trace and type_array have different lengths
    """
    import numpy as np

    trace = np.array(trace).flatten()
    type_array = np.array(type_array).flatten()

    if len(trace) != len(type_array):
        raise ValueError("trace and type_array must have same length")

    if lags is None:
        lags = [0, 1, 2, 3, 4]
    lags = np.array(lags)

    nlags = len(lags)
    covariances = np.zeros((ntypes, ntypes, nlags))

    for i in range(ntypes):
        for j in range(ntypes):
            for k, lag in enumerate(lags):
                if lag == 0:
                    mask_i = (type_array == i)
                    mask_j = (type_array == j)
                    if np.sum(mask_i) > 0 and np.sum(mask_j) > 0:
                        if i == j:
                            covariances[i, j, k] = np.var(trace[mask_i], ddof=1)
                        else:
                            covariances[i, j, k] = np.cov(
                                trace[mask_i], trace[mask_j], ddof=1
                            )[0, 1] if len(trace[mask_i]) > 1 and len(trace[mask_j]) > 1 else 0.0
                else:
                    if len(trace) > lag:
                        trace1 = trace[:-lag]
                        trace2 = trace[lag:]
                        type1 = type_array[:-lag]
                        type2 = type_array[lag:]

                        mask1 = (type1 == i)
                        mask2 = (type2 == j)

                        if np.sum(mask1) > 0 and np.sum(mask2) > 0:
                            covariances[i, j, k] = np.cov(
                                trace1[mask1], trace2[mask2], ddof=1
                            )[0, 1] if len(trace1[mask1]) > 1 and len(trace2[mask2]) > 1 else 0.0

    return {
        'covariances': covariances,
        'lags': lags
    }

# Additional MTRACE functions
def mtrace_backward_moment(traces, order):
    """
    Computes backward moments for multiple traces.

    Args:
        traces: List of trace data
        order: Moment order

    Returns:
        numpy.ndarray: Backward moments
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_backward_momentKt.mtrace_backward_moment(
        jlineMatrixFromArray(traces), int(order)
    )
    return jlineMatrixToArray(result)


def mtrace_bootstrap(traces, n_bootstrap=100, statistic='mean'):
    """
    Performs bootstrap analysis on multiple traces.

    Args:
        traces: List of trace data
        n_bootstrap: Number of bootstrap samples
        statistic: Statistic to compute ('mean', 'var', etc.)

    Returns:
        dict: Bootstrap results
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_bootstrapKt.mtrace_bootstrap(
        jlineMatrixFromArray(traces), int(n_bootstrap), jpype.JString(statistic)
    )
    return {
        'estimates': jlineMatrixToArray(result.getEstimates()),
        'confidence_interval': jlineMatrixToArray(result.getConfidenceInterval()),
        'std_error': float(result.getStdError())
    }


def mtrace_count(traces):
    """
    Counts events in multiple traces.

    Args:
        traces: List of trace data

    Returns:
        numpy.ndarray: Event counts per trace
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_countKt.mtrace_count(
        jlineMatrixFromArray(traces)
    )
    return jlineMatrixToArray(result)


def mtrace_cross_moment(trace1, trace2, order1, order2):
    """
    Computes cross moments between two traces.

    Args:
        trace1: First trace data
        trace2: Second trace data
        order1: Moment order for trace1
        order2: Moment order for trace2

    Returns:
        float: Cross moment value
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_cross_momentKt.mtrace_cross_moment(
        jlineMatrixFromArray(trace1), jlineMatrixFromArray(trace2),
        int(order1), int(order2)
    )
    return float(result)


def mtrace_forward_moment(traces, order):
    """
    Computes forward moments for multiple traces.

    Args:
        traces: List of trace data
        order: Moment order

    Returns:
        numpy.ndarray: Forward moments
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_forward_momentKt.mtrace_forward_moment(
        jlineMatrixFromArray(traces), int(order)
    )
    return jlineMatrixToArray(result)


def mtrace_iat2counts(traces, interval_length):
    """
    Converts inter-arrival times to counts for multiple traces.

    Args:
        traces: List of inter-arrival time traces
        interval_length: Length of counting interval

    Returns:
        dict: Count data for each trace
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_iat2countsKt.mtrace_iat2counts(
        jlineMatrixFromArray(traces), float(interval_length)
    )
    return {
        'counts': jlineMatrixToArray(result.getCounts()),
        'intervals': jlineMatrixToArray(result.getIntervals())
    }


def mtrace_joint(traces, bins=None):
    """
    Computes joint distribution for multiple traces.

    Args:
        traces: List of trace data
        bins: Number of bins for histogram (optional)

    Returns:
        dict: Joint distribution information
    """
    if bins is not None:
        result = jpype.JPackage('jline').api.trace.Mtrace_jointKt.mtrace_joint(
            jlineMatrixFromArray(traces), int(bins)
        )
    else:
        result = jpype.JPackage('jline').api.trace.Mtrace_jointKt.mtrace_joint(
            jlineMatrixFromArray(traces)
        )

    return {
        'joint_distribution': jlineMatrixToArray(result.getJointDistribution()),
        'marginals': jlineMatrixToArray(result.getMarginals()),
        'bins': jlineMatrixToArray(result.getBins()) if hasattr(result, 'getBins') else None
    }


def mtrace_merge(traces):
    """
    Merges multiple traces into a single trace.

    Args:
        traces: List of trace data to merge

    Returns:
        numpy.ndarray: Merged trace
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_mergeKt.mtrace_merge(
        jlineMatrixFromArray(traces)
    )
    return jlineMatrixToArray(result)


def mtrace_moment(traces, order):
    """
    Computes moments of specified order for multiple traces.

    Args:
        traces: List of trace data
        order: Moment order

    Returns:
        numpy.ndarray: Moments for each trace
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_momentKt.mtrace_moment(
        jlineMatrixFromArray(traces), int(order)
    )
    return jlineMatrixToArray(result)


def mtrace_moment_simple(traces):
    """
    Computes simple moments (mean, variance) for multiple traces.

    Args:
        traces: List of trace data

    Returns:
        dict: Simple moments for each trace
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_moment_simpleKt.mtrace_moment_simple(
        jlineMatrixFromArray(traces)
    )
    return {
        'means': jlineMatrixToArray(result.getMeans()),
        'variances': jlineMatrixToArray(result.getVariances())
    }


def mtrace_pc(traces):
    """
    Computes principal components for multiple traces.

    Args:
        traces: List of trace data

    Returns:
        dict: Principal component analysis results
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_pcKt.mtrace_pc(
        jlineMatrixFromArray(traces)
    )
    return {
        'components': jlineMatrixToArray(result.getComponents()),
        'explained_variance': jlineMatrixToArray(result.getExplainedVariance()),
        'loadings': jlineMatrixToArray(result.getLoadings())
    }


def mtrace_sigma(traces):
    """
    Computes sigma parameter for multiple traces.

    Args:
        traces: List of trace data

    Returns:
        numpy.ndarray: Sigma values for each trace
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_sigmaKt.mtrace_sigma(
        jlineMatrixFromArray(traces)
    )
    return jlineMatrixToArray(result)


def mtrace_sigma2(traces):
    """
    Computes squared sigma parameter for multiple traces.

    Args:
        traces: List of trace data

    Returns:
        numpy.ndarray: Squared sigma values for each trace
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_sigma2Kt.mtrace_sigma2(
        jlineMatrixFromArray(traces)
    )
    return jlineMatrixToArray(result)


def mtrace_split(trace, n_parts):
    """
    Splits a trace into multiple parts.

    Args:
        trace: Trace data to split
        n_parts: Number of parts to split into

    Returns:
        list: List of trace parts
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_splitKt.mtrace_split(
        jlineMatrixFromArray(trace), int(n_parts)
    )
    return [jlineMatrixToArray(part) for part in result]


def mtrace_summary(traces):
    """
    Computes summary statistics for multiple traces.

    Args:
        traces: List of trace data

    Returns:
        dict: Summary statistics
    """
    result = jpype.JPackage('jline').api.trace.Mtrace_summaryKt.mtrace_summary(
        jlineMatrixFromArray(traces)
    )
    return {
        'count': jlineMatrixToArray(result.getCount()),
        'mean': jlineMatrixToArray(result.getMean()),
        'std': jlineMatrixToArray(result.getStd()),
        'min': jlineMatrixToArray(result.getMin()),
        'percentiles': jlineMatrixToArray(result.getPercentiles()),
        'max': jlineMatrixToArray(result.getMax())
    }

