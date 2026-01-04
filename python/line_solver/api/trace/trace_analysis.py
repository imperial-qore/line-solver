"""
Trace Analysis Functions.

Native Python implementations for statistical analysis of empirical trace data.
Provides functions for computing means, variances, correlations, and other
statistical measures from measurement traces.
"""

import numpy as np
from typing import Union, Tuple, Optional

ArrayLike = Union[np.ndarray, list]


def trace_mean(trace: ArrayLike) -> float:
    """
    Compute the arithmetic mean of trace data.

    Args:
        trace: Array of trace values.

    Returns:
        Mean value of the trace.

    Example:
        >>> trace_mean([1.0, 2.0, 3.0, 4.0, 5.0])
        3.0
    """
    trace = np.asarray(trace, dtype=np.float64).ravel()
    return float(np.mean(trace))


def trace_var(trace: ArrayLike) -> float:
    """
    Compute the variance of trace data.

    Uses population variance (ddof=0) for consistency with Kotlin.

    Args:
        trace: Array of trace values.

    Returns:
        Variance of the trace.

    Example:
        >>> trace_var([1.0, 2.0, 3.0, 4.0, 5.0])
        2.0
    """
    trace = np.asarray(trace, dtype=np.float64).ravel()
    e1 = np.mean(trace)
    e2 = np.mean(trace ** 2)
    return float(e2 - e1 * e1)


def trace_scv(trace: ArrayLike) -> float:
    """
    Compute the squared coefficient of variation (SCV).

    SCV = Var(X) / E[X]^2

    Args:
        trace: Array of trace values.

    Returns:
        Squared coefficient of variation.

    Example:
        >>> trace_scv([1.0, 2.0, 3.0])  # Var=0.667, Mean=2, SCV=0.167
    """
    mean = trace_mean(trace)
    var = trace_var(trace)
    return float(var / (mean * mean)) if mean != 0 else 0.0


def _autocov(X: np.ndarray) -> np.ndarray:
    """Compute autocovariance using direct method."""
    M = len(X)
    acv = np.zeros(M - 1)

    for p in range(M - 1):
        acv[p] = np.mean(X[:M - p] * X[p:])

    return acv


def trace_acf(trace: ArrayLike, lags: ArrayLike = None) -> np.ndarray:
    """
    Compute the autocorrelation function at specified lags.

    Args:
        trace: Array of trace values.
        lags: Array of lag values (default: [1]).

    Returns:
        Array of autocorrelation values at each lag.

    Example:
        >>> trace = np.random.randn(100)
        >>> acf = trace_acf(trace, [1, 2, 3])
    """
    trace = np.asarray(trace, dtype=np.float64).ravel()
    if lags is None:
        lags = np.array([1])
    lags = np.asarray(lags, dtype=np.int32).ravel()

    max_lag = np.max(lags)

    if max_lag > len(trace) - 2:
        # Filter out lags that are too large
        valid_lags = lags[(lags <= len(trace) - 2) & (lags > 0)]
        if len(valid_lags) == 0:
            return np.array([])
        lags = valid_lags

    mean = trace_mean(trace)
    centered = trace - mean

    autocov = _autocov(centered)
    rho = np.zeros(len(lags))

    for i, lag in enumerate(lags):
        if lag < len(autocov):
            rho[i] = autocov[lag] / autocov[0] if autocov[0] != 0 else 0.0

    return rho


def trace_gamma(trace: ArrayLike, limit: int = 1000) -> np.ndarray:
    """
    Estimate the autocorrelation decay rate of a trace.

    Args:
        trace: Array of trace values.
        limit: Maximum lag considered.

    Returns:
        Array containing [GAMMA, RHO0, RESIDUALS].

    Example:
        >>> gamma, rho0, residuals = trace_gamma(trace_data)
    """
    trace = np.asarray(trace, dtype=np.float64).ravel()

    M1 = trace_mean(trace)
    M2 = np.mean(trace ** 2)

    max_lag = min(limit, len(trace) - 1)
    lag = np.arange(1, max_lag + 1)
    rho = trace_acf(trace, lag)

    VAR = M2 - M1 * M1
    SCV = VAR / (M1 * M1) if M1 != 0 else 1.0
    RHO0 = 0.5 * (1.0 - 1.0 / SCV) if SCV != 0 else 0.0

    # Grid search for best gamma
    best_gamma = 0.99
    min_residuals = np.inf

    for gamma_int in range(990, 1000):
        g = gamma_int / 1000.0
        expected = RHO0 * (g ** lag[:len(rho)])
        residuals = np.sum((rho - expected) ** 2)
        if residuals < min_residuals:
            min_residuals = residuals
            best_gamma = g

    return np.array([best_gamma, RHO0, min_residuals])


def trace_iat2counts(trace: ArrayLike, scale: float) -> np.ndarray:
    """
    Compute the counting process from inter-arrival times.

    Args:
        trace: Array of inter-arrival times.
        scale: Time scale for counting.

    Returns:
        Array of counts after `scale` units of time from each arrival.

    Example:
        >>> iat = [0.5, 0.3, 0.8, 0.2, 0.4]
        >>> counts = trace_iat2counts(iat, 1.0)
    """
    S = np.asarray(trace, dtype=np.float64).ravel()
    n = len(S)

    # Cumulative sum with 0 at start
    CS = np.zeros(n + 1)
    CS[1:] = np.cumsum(S)

    C = []
    for i in range(n - 1):
        cur = i
        while cur + 1 < n and CS[cur + 1] - CS[i] <= scale:
            cur += 1
        C.append(cur - i)
        if cur == n - 1:
            break

    return np.array(C, dtype=np.int32)


def trace_idi(trace: ArrayLike, kset: ArrayLike, option: str = None, n: int = 1
              ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the Index of Dispersion for Intervals.

    Args:
        trace: Array of trace values.
        kset: Set of k values to compute IDI for.
        option: Aggregation option (None, 'aggregate', 'aggregate-mix').
        n: Aggregation parameter.

    Returns:
        Tuple of (IDI values, support values).

    Example:
        >>> idi, support = trace_idi(trace_data, [10, 20, 50])
    """
    S = np.asarray(trace, dtype=np.float64).ravel()
    kset = np.asarray(kset, dtype=np.int32).ravel()

    IDIk = []
    support = []

    for k in kset:
        if option is None:
            support_val = len(S) - k - 1
            if support_val <= 0:
                continue
            support.append(support_val)

            Sk = np.zeros(len(S) - k)
            for t in range(len(S) - k - 1):
                Sk[t] = np.sum(S[t:t + k])

            variance = trace_var(Sk)
            mean = trace_mean(Sk)
            IDIk.append(k * variance / (mean * mean) if mean != 0 else 0.0)

        elif option == 'aggregate':
            keff = k // n
            if keff <= 0:
                continue
            support_val = len(S) // keff
            support.append(support_val)

            Sk = np.zeros(len(S) - keff)
            for t in range(len(S) - keff - 1):
                Sk[t] = np.sum(S[t:t + keff])

            variance = trace_var(Sk)
            mean = trace_mean(Sk)
            IDIk.append(k * variance / (mean * mean) if mean != 0 else 0.0)

    return np.array(IDIk), np.array(support, dtype=np.int32)


def trace_idc(trace: ArrayLike) -> float:
    """
    Compute the Index of Dispersion for Counts.

    Asymptotically equal to IDI.

    Args:
        trace: Array of trace values.

    Returns:
        IDC value.

    Example:
        >>> idc = trace_idc(inter_arrival_times)
    """
    S = np.asarray(trace, dtype=np.float64).ravel()
    k = min(1000, len(S) // 30)
    if k < 1:
        k = 1

    idi, _ = trace_idi(S, [k])
    return float(idi[0]) if len(idi) > 0 else 0.0


def trace_pmf(X: ArrayLike) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the probability mass function of discrete data.

    Args:
        X: Array of discrete values.

    Returns:
        Tuple of (PMF values, unique values).

    Example:
        >>> pmf, values = trace_pmf([1, 2, 2, 3, 3, 3])
    """
    X = np.asarray(X, dtype=np.int32).ravel()
    unique_values, counts = np.unique(X, return_counts=True)
    pmf = counts / len(X)
    return pmf, unique_values


def trace_shuffle(trace: ArrayLike) -> np.ndarray:
    """
    Shuffle trace data randomly.

    Args:
        trace: Array of trace values.

    Returns:
        Shuffled trace array.

    Example:
        >>> shuffled = trace_shuffle([1, 2, 3, 4, 5])
    """
    trace = np.asarray(trace, dtype=np.float64).ravel()
    result = trace.copy()
    np.random.shuffle(result)
    return result


def trace_joint(trace: ArrayLike, lag: ArrayLike, order: ArrayLike) -> float:
    """
    Compute joint moments E[X^{k_1}_{i} * X^{k_2}_{i+j} * ...].

    Args:
        trace: Array of trace values.
        lag: Cumulative lag values.
        order: Moment orders.

    Returns:
        Joint moment value.

    Example:
        >>> jm = trace_joint(trace, [0, 1], [1, 1])  # E[X_i * X_{i+1}]
    """
    S = np.asarray(trace, dtype=np.float64).ravel()
    lag = np.asarray(lag, dtype=np.int32).ravel()
    order = np.asarray(order, dtype=np.int32).ravel()

    sorted_lag = np.sort(lag)
    K = len(sorted_lag)
    base_lag = sorted_lag[0]
    adjusted_lag = sorted_lag - base_lag

    max_lag = np.max(adjusted_lag)
    valid_length = len(S) - max_lag

    if valid_length <= 0:
        return 0.0

    result = 0.0
    for i in range(valid_length):
        product = 1.0
        for j, ord_val in enumerate(order):
            idx = i + adjusted_lag[min(j, len(adjusted_lag) - 1)]
            if idx < len(S):
                product *= S[idx] ** ord_val
        result += product

    return float(result / valid_length)


def trace_iat2bins(trace: ArrayLike, scale: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute counts in bins with specified timescale.

    Args:
        trace: Array of inter-arrival times.
        scale: Bin timescale.

    Returns:
        Tuple of (counts per bin, bin membership for each element).

    Example:
        >>> counts, bins = trace_iat2bins(iat_data, 1.0)
    """
    S = np.asarray(trace, dtype=np.float64).ravel()
    n = len(S)

    CS = np.zeros(n + 1)
    CS[1:] = np.cumsum(S)

    num_bins = int(np.ceil((CS[n] - CS[0]) / scale))
    C = np.zeros(num_bins, dtype=np.int32)
    bC = []

    cur = 0
    last = 0

    for i in range(num_bins):
        if cur >= n - 1:
            break

        while cur < n - 1 and CS[cur + 1] <= (i + 1) * scale:
            cur += 1

        C[i] = cur - last
        for _ in range(cur - last):
            bC.append(i)
        last = cur

    return C, np.array(bC, dtype=np.int32)


def _percentile(sorted_data: np.ndarray, p: float) -> float:
    """Compute percentile of sorted data."""
    index = (p / 100.0) * (len(sorted_data) - 1)
    lower = int(np.floor(index))
    upper = int(np.ceil(index))

    if lower == upper:
        return float(sorted_data[lower])

    weight = index - lower
    return float(sorted_data[lower] * (1 - weight) + sorted_data[upper] * weight)


def trace_summary(trace: ArrayLike) -> np.ndarray:
    """
    Compute comprehensive summary statistics for a trace.

    Args:
        trace: Array of trace values.

    Returns:
        Array containing [MEAN, SCV, MAD, SKEW, KURT, Q25, Q50, Q75, P95,
        MIN, MAX, IQR, ACF1, ACF2, ACF3, ACF4, IDC_SCV_RATIO].

    Example:
        >>> summary = trace_summary(trace_data)
        >>> print(f"Mean: {summary[0]}, SCV: {summary[1]}")
    """
    m = np.asarray(trace, dtype=np.float64).ravel()

    mean = trace_mean(m)
    scv = trace_scv(m)
    sorted_m = np.sort(m)

    # Percentiles
    q25 = _percentile(sorted_m, 25.0)
    q50 = _percentile(sorted_m, 50.0)  # median
    q75 = _percentile(sorted_m, 75.0)
    p95 = _percentile(sorted_m, 95.0)
    min_val = float(sorted_m[0])
    max_val = float(sorted_m[-1])
    iqr = q75 - q25

    # MAD (median absolute deviation)
    mad = float(np.median(np.abs(m - q50)))

    # Skewness and kurtosis
    variance = trace_var(m)
    std = np.sqrt(variance) if variance > 0 else 1.0
    z = (m - mean) / std
    skew = float(np.mean(z ** 3))
    kurt = float(np.mean(z ** 4) - 3.0)  # excess kurtosis

    # ACF for lags 1-4
    acf = trace_acf(m, [1, 2, 3, 4])
    if len(acf) < 4:
        acf = np.concatenate([acf, np.zeros(4 - len(acf))])

    # IDC/SCV ratio
    idc = trace_idc(m)
    idc_scv_ratio = idc / scv if scv != 0 else 0.0

    return np.array([
        mean, scv, mad, skew, kurt, q25, q50, q75, p95,
        min_val, max_val, iqr, acf[0], acf[1], acf[2], acf[3], idc_scv_ratio
    ])


def mtrace_mean(trace: ArrayLike, ntypes: int, types: ArrayLike) -> np.ndarray:
    """
    Compute the mean of a trace, divided by types.

    Args:
        trace: Array of trace values.
        ntypes: Number of different types.
        types: Array indicating the type of each element.

    Returns:
        Array containing the mean values for each type.

    Example:
        >>> trace = [1.0, 2.0, 3.0, 4.0]
        >>> types = [0, 1, 0, 1]
        >>> means = mtrace_mean(trace, 2, types)
    """
    trace = np.asarray(trace, dtype=np.float64).ravel()
    types = np.asarray(types, dtype=np.int32).ravel()

    mean = np.zeros(ntypes)
    for c in range(ntypes):
        mask = types == c
        if np.any(mask):
            mean[c] = np.mean(trace[mask])
        else:
            mean[c] = np.nan

    return mean


def mtrace_var(trace: ArrayLike, ntypes: int, types: ArrayLike) -> np.ndarray:
    """
    Compute the variance of a trace, divided by types.

    Args:
        trace: Array of trace values.
        ntypes: Number of different types.
        types: Array indicating the type of each element.

    Returns:
        Array containing the variance values for each type.

    Example:
        >>> var_by_type = mtrace_var(trace, 2, types)
    """
    trace = np.asarray(trace, dtype=np.float64).ravel()
    types = np.asarray(types, dtype=np.int32).ravel()

    var = np.zeros(ntypes)
    for c in range(ntypes):
        mask = types == c
        if np.sum(mask) > 1:
            class_trace = trace[mask]
            e1 = np.mean(class_trace)
            e2 = np.mean(class_trace ** 2)
            var[c] = e2 - e1 * e1
        else:
            var[c] = np.nan

    return var


def mtrace_count(trace: ArrayLike, ntypes: int, types: ArrayLike) -> np.ndarray:
    """
    Count elements per type in a trace.

    Args:
        trace: Array of trace values.
        ntypes: Number of different types.
        types: Array indicating the type of each element.

    Returns:
        Array containing counts for each type.

    Example:
        >>> counts = mtrace_count(trace, 2, types)
    """
    types = np.asarray(types, dtype=np.int32).ravel()

    counts = np.zeros(ntypes, dtype=np.int32)
    for c in range(ntypes):
        counts[c] = np.sum(types == c)

    return counts


def mtrace_sigma(T: ArrayLike, L: ArrayLike) -> np.ndarray:
    """
    Compute one-step class transition probabilities from a marked trace.

    Computes P(C_k = j | C_{k-1} = i) empirically from the trace.

    Args:
        T: Array of inter-arrival times (unused, for API consistency)
        L: Array of class labels

    Returns:
        C x C matrix where element (i,j) is the probability of observing
        class j after class i.
    """
    L = np.asarray(L).ravel()
    marks = np.unique(L)
    C = len(marks)

    sigma = np.zeros((C, C))

    for i in range(C):
        for j in range(C):
            count = np.sum((L[:-1] == marks[i]) & (L[1:] == marks[j]))
            sigma[i, j] = count / (len(L) - 1)

    return sigma


def mtrace_sigma2(T: ArrayLike, L: ArrayLike) -> np.ndarray:
    """
    Compute two-step class transition probabilities from a marked trace.

    Computes P(C_k = h | C_{k-1} = j, C_{k-2} = i) empirically.

    Args:
        T: Array of inter-arrival times (unused, for API consistency)
        L: Array of class labels

    Returns:
        C x C x C 3D array of transition probabilities.
    """
    L = np.asarray(L).ravel()
    marks = np.unique(L)
    C = len(marks)

    sigma = np.zeros((C, C, C))

    for i in range(C):
        for j in range(C):
            for h in range(C):
                count = np.sum((L[:-2] == marks[i]) & (L[1:-1] == marks[j]) & (L[2:] == marks[h]))
                sigma[i, j, h] = count / (len(L) - 2)

    return sigma


def mtrace_cross_moment(T: ArrayLike, L: ArrayLike, k: int) -> np.ndarray:
    """
    Compute the k-th order moment of inter-arrival times between class pairs.

    Args:
        T: Array of inter-arrival times
        L: Array of class labels
        k: Order of the moment

    Returns:
        C x C matrix where element (i,j) is E[T^k | C_{t-1} = i, C_t = j]
    """
    T = np.asarray(T, dtype=np.float64).ravel()
    L = np.asarray(L).ravel()
    marks = np.unique(L)
    C = len(marks)

    MC = np.zeros((C, C))
    count = np.zeros((C, C))

    for t in range(1, len(T)):
        for i in range(C):
            for j in range(C):
                if L[t - 1] == marks[i] and L[t] == marks[j]:
                    MC[i, j] += T[t] ** k
                    count[i, j] += 1

    # Avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        MC = np.where(count > 0, MC / count, np.nan)

    return MC


def mtrace_forward_moment(T: ArrayLike, A: ArrayLike, orders: ArrayLike,
                          norm: bool = True) -> np.ndarray:
    """
    Compute forward moments of a marked trace.

    Forward moment for class c is E[T_{k+1}^order | C_k = c].

    Args:
        T: Array of inter-arrival times
        A: Array of class labels
        orders: Array of moment orders to compute
        norm: If True, normalize by class probability

    Returns:
        Matrix of shape (C, len(orders)) containing forward moments.
    """
    T = np.asarray(T, dtype=np.float64).ravel()
    A = np.asarray(A).ravel()
    orders = np.asarray(orders).ravel()
    marks = np.unique(A)
    C = len(marks)

    M = np.zeros((C, len(orders)))

    for j, k in enumerate(orders):
        for c in range(C):
            mask = A[:-1] == marks[c]
            if np.any(mask):
                M[c, j] = np.mean((T[1:] ** k) * mask)
                if norm:
                    M[c, j] = M[c, j] * (len(T) - 1) / np.sum(mask)
            else:
                M[c, j] = np.nan

    return M


def mtrace_backward_moment(T: ArrayLike, A: ArrayLike, orders: ArrayLike,
                           norm: bool = True) -> np.ndarray:
    """
    Compute backward moments of a marked trace.

    Backward moment for class c is E[T_k^order | C_k = c].

    Args:
        T: Array of inter-arrival times
        A: Array of class labels
        orders: Array of moment orders to compute
        norm: If True, normalize by class probability

    Returns:
        Matrix of shape (C, len(orders)) containing backward moments.
    """
    T = np.asarray(T, dtype=np.float64).ravel()
    A = np.asarray(A).ravel()
    orders = np.asarray(orders).ravel()
    marks = np.unique(A)
    C = len(marks)

    M = np.zeros((C, len(orders)))

    for j, k in enumerate(orders):
        for c in range(C):
            mask = A == marks[c]
            if np.any(mask):
                M[c, j] = np.mean((T ** k) * mask)
                if norm:
                    M[c, j] = M[c, j] * len(T) / np.sum(mask)
            else:
                M[c, j] = np.nan

    return M


def mtrace_cov(T: ArrayLike, A: ArrayLike) -> np.ndarray:
    """
    Compute lag-1 covariance between classes in a marked trace.

    Args:
        T: Array of inter-arrival times
        A: Array of class labels

    Returns:
        C x C matrix of 2x2 covariance matrices.
    """
    T = np.asarray(T, dtype=np.float64).ravel()
    A = np.asarray(A).ravel()
    C = int(np.max(A))
    N = len(A)

    # Return list of covariance matrices
    COV = [[None for _ in range(C)] for _ in range(C)]

    for c1 in range(C):
        for c2 in range(C):
            X0c1v = np.zeros(N - 1)
            X1c2v = np.zeros(N - 1)

            for i in range(N - 1):
                if A[i] == c1 + 1:  # 1-indexed classes
                    X0c1v[i] = T[i]
                if A[i + 1] == c2 + 1:
                    X1c2v[i] = T[i + 1]

            COV[c1][c2] = np.cov(X0c1v, X1c2v)

    return COV


def mtrace_pc(T: ArrayLike, L: ArrayLike) -> np.ndarray:
    """
    Compute the probability of arrival for each class.

    Args:
        T: Array of inter-arrival times (unused, for API consistency)
        L: Array of class labels

    Returns:
        Array of class probabilities.
    """
    L = np.asarray(L).ravel()
    labels = np.unique(L)
    m = len(labels)

    pc = np.zeros(m)
    for i in range(m):
        pc[i] = np.sum(L == labels[i]) / len(L)

    return pc


def mtrace_summary(T: ArrayLike, L: ArrayLike) -> dict:
    """
    Compute comprehensive summary statistics for a marked trace.

    Args:
        T: Array of inter-arrival times
        L: Array of class labels

    Returns:
        Dictionary containing:
        - M: First 5 aggregate moments
        - ACF: Autocorrelation for lags 1-100
        - F1, F2: Forward moments of order 1 and 2
        - B1, B2: Backward moments of order 1 and 2
        - C1, C2: Cross moments of order 1 and 2
        - Pc: Class probabilities
        - Pab: One-step transition probabilities
    """
    T = np.asarray(T, dtype=np.float64).ravel()
    L = np.asarray(L).ravel()

    summary = {
        'M': np.array([np.mean(T ** k) for k in range(1, 6)]),
        'ACF': trace_acf(T, np.arange(1, 101)),
        'F1': mtrace_forward_moment(T, L, [1]),
        'F2': mtrace_forward_moment(T, L, [2]),
        'B1': mtrace_backward_moment(T, L, [1]),
        'B2': mtrace_backward_moment(T, L, [2]),
        'C1': mtrace_cross_moment(T, L, 1),
        'C2': mtrace_cross_moment(T, L, 2),
        'Pc': mtrace_pc(T, L),
        'Pab': mtrace_sigma(T, L),
    }

    return summary


def mtrace_split(T: ArrayLike, L: ArrayLike) -> list:
    """
    Split a multi-class trace into per-class traces.

    For each class, computes inter-arrival times between consecutive
    events of that class.

    Args:
        T: Array of inter-arrival times
        L: Array of class labels

    Returns:
        List of arrays, one per class, containing inter-arrival times.
    """
    T = np.asarray(T, dtype=np.float64).ravel()
    L = np.asarray(L).ravel()
    labels = np.unique(L)
    C = len(labels)

    TCUM = np.cumsum(T)

    TL = []
    for c in range(C):
        mask = L == labels[c]
        cum_times = np.concatenate([[0], TCUM[mask]])
        TL.append(np.diff(cum_times))

    return TL


def mtrace_merge(t1: ArrayLike, t2: ArrayLike) -> Tuple[np.ndarray, np.ndarray]:
    """
    Merge two traces into a single marked (multi-class) trace.

    Args:
        t1: Inter-arrival times of the first trace
        t2: Inter-arrival times of the second trace

    Returns:
        Tuple of (T, L) where:
        - T: Merged inter-arrival times
        - L: Class labels (1 for t1, 2 for t2)
    """
    t1 = np.asarray(t1, dtype=np.float64).ravel()
    t2 = np.asarray(t2, dtype=np.float64).ravel()

    # Compute cumulative times
    cum1 = np.cumsum(t1)
    cum2 = np.cumsum(t2)

    # Combine and sort
    all_times = np.concatenate([[0], cum1, cum2])
    sorted_indices = np.argsort(all_times)
    sorted_times = all_times[sorted_indices]

    # Compute inter-arrival times
    T = np.diff(sorted_times)

    # Assign labels
    IDX = sorted_indices[1:]  # Skip the leading 0
    L = np.zeros(len(T), dtype=np.int32)

    # Label 1 for t1 events (indices 1 to len(t1))
    # Label 2 for t2 events (indices len(t1)+1 onwards)
    L[(IDX >= 1) & (IDX <= len(t1))] = 1
    L[(IDX > len(t1))] = 2

    return T, L


def mtrace_joint(T: ArrayLike, A: ArrayLike, i: ArrayLike) -> np.ndarray:
    """
    Compute class-dependent joint moments.

    Computes E[(X^(a)_j)^i[0] * (X^(a)_{j+1})^i[1]] for all classes a.

    Args:
        T: Array of inter-arrival times
        A: Array of class labels
        i: Array of exponents [i0, i1]

    Returns:
        Array of joint moments, one per class.
    """
    T = np.asarray(T, dtype=np.float64).ravel()
    A = np.asarray(A).ravel()
    i = np.asarray(i).ravel()
    C = int(np.max(A))
    N = len(A)

    JM = np.zeros(C)

    for a in range(1, C + 1):
        # Events of class a, excluding first and last
        Ta = T[(A[1:-1] == a)]
        Na = len(Ta)

        if Na == 0:
            JM[a - 1] = np.nan
            continue

        tmp = 0.0
        for j in range(N - 2):
            if A[j + 1] == a:
                tmp += (T[j] ** i[0]) * (T[j + 1] ** i[1])

        JM[a - 1] = tmp / Na

    return JM


def mtrace_moment(T: ArrayLike, A: ArrayLike, orders: ArrayLike,
                  after: bool = False, norm: bool = False) -> np.ndarray:
    """
    Compute class-dependent moments of a multi-class trace.

    Args:
        T: Array of inter-arrival times
        A: Array of class labels
        orders: Array of moment orders to compute
        after: If True, compute moments of Bucholz variables (forward)
               If False, compute moments of Horvath variables (backward)
        norm: If True, normalize by class probability

    Returns:
        Matrix of shape (C, len(orders)) containing moments per class.
    """
    T = np.asarray(T, dtype=np.float64).ravel()
    A = np.asarray(A).ravel()
    orders = np.asarray(orders).ravel()
    marks = np.unique(A)
    C = len(marks)

    M = np.zeros((C, len(orders)))

    for j, k in enumerate(orders):
        for c in range(C):
            if after:
                mask = A[:-1] == marks[c]
                if np.any(mask):
                    M[c, j] = np.mean((T[1:] ** k) * mask)
                    if norm:
                        M[c, j] = M[c, j] * (len(T) - 1) / np.sum(mask)
            else:
                mask = A == marks[c]
                if np.any(mask):
                    M[c, j] = np.mean((T ** k) * mask)
                    if norm:
                        M[c, j] = M[c, j] * len(T) / np.sum(mask)

    return M


def mtrace_moment_simple(T: ArrayLike, A: ArrayLike, k: int) -> np.ndarray:
    """
    Simple interface to compute k-th order moments per class.

    Args:
        T: Array of inter-arrival times
        A: Array of class labels
        k: Order of the moment

    Returns:
        Array of k-th order moments, one per class.
    """
    return mtrace_moment(T, A, [k], after=False, norm=True)[:, 0]


def mtrace_bootstrap(T: ArrayLike, A: ArrayLike, n_samples: int = 100,
                     seed: int = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate bootstrap samples from a marked trace.

    Args:
        T: Array of inter-arrival times
        A: Array of class labels
        n_samples: Number of bootstrap samples to generate
        seed: Random seed

    Returns:
        Tuple of (T_boot, A_boot) containing bootstrap samples.
    """
    if seed is not None:
        np.random.seed(seed)

    T = np.asarray(T, dtype=np.float64).ravel()
    A = np.asarray(A).ravel()
    n = len(T)

    indices = np.random.choice(n, size=n_samples, replace=True)

    return T[indices], A[indices]


def mtrace_iat2counts(T: ArrayLike, L: ArrayLike, scale: float) -> np.ndarray:
    """
    Compute counting process from marked inter-arrival times.

    For each class, counts arrivals in windows of specified scale.

    Args:
        T: Array of inter-arrival times
        L: Array of class labels
        scale: Time scale for counting

    Returns:
        Array of counts per class per window.
    """
    T = np.asarray(T, dtype=np.float64).ravel()
    L = np.asarray(L).ravel()
    labels = np.unique(L)
    C = len(labels)

    # Split by class and compute counts
    TL = mtrace_split(T, L)

    counts = []
    for c in range(C):
        if len(TL[c]) > 0:
            counts.append(trace_iat2counts(TL[c], scale))
        else:
            counts.append(np.array([]))

    return counts


def trace_bicov(trace: ArrayLike, grid: ArrayLike) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute bicovariance of a trace.

    Bicovariance measures third-order correlation structure at
    multiple lag combinations.

    Args:
        trace: Array of trace values.
        grid: Array of lags to form the grid (e.g., [1, 2, 3, 4, 5]).

    Returns:
        Tuple of (bicov, bicov_lags) where:
        - bicov: Array of bicovariance values
        - bicov_lags: 2D array of lag combinations (each row is [1, i, j])

    Example:
        >>> trace = np.random.randn(1000)
        >>> bicov, lags = trace_bicov(trace, [1, 2, 3, 4, 5])
    """
    trace = np.asarray(trace, dtype=np.float64).ravel()
    grid = np.asarray(grid).ravel()

    bicov_lags = []
    for i in grid:
        for j in grid:
            bicov_lags.append([1, int(i), int(j)])

    bicov_lags = np.array(bicov_lags)
    bicov = np.zeros(len(bicov_lags))

    for k, lags in enumerate(bicov_lags):
        bicov[k] = trace_joint(trace, lags, [1, 1, 1])

    return bicov, bicov_lags


__all__ = [
    # Single trace analysis
    'trace_mean',
    'trace_var',
    'trace_scv',
    'trace_acf',
    'trace_gamma',
    'trace_iat2counts',
    'trace_idi',
    'trace_idc',
    'trace_pmf',
    'trace_shuffle',
    'trace_joint',
    'trace_iat2bins',
    'trace_summary',
    'trace_bicov',
    # Multi-class trace analysis
    'mtrace_mean',
    'mtrace_var',
    'mtrace_count',
    'mtrace_sigma',
    'mtrace_sigma2',
    'mtrace_cross_moment',
    'mtrace_forward_moment',
    'mtrace_backward_moment',
    'mtrace_cov',
    'mtrace_pc',
    'mtrace_summary',
    'mtrace_split',
    'mtrace_merge',
    'mtrace_joint',
    'mtrace_moment',
    'mtrace_moment_simple',
    'mtrace_bootstrap',
    'mtrace_iat2counts',
]
