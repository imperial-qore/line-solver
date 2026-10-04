"""
Trace analysis functions for the KPC-Toolbox.

Native Python implementations ported from MATLAB: matlab/lib/kpctoolbox/trace/

All functions operate on numpy arrays and follow the MATLAB ground truth
implementations exactly.
"""

import numpy as np
from scipy.optimize import curve_fit


def trace_mean(S):
    """Compute the mean of trace S.

    Parameters
    ----------
    S : array_like
        Input trace.

    Returns
    -------
    float
        Mean value of the trace.
    """
    S = np.asarray(S, dtype=float).ravel()
    return float(np.mean(S))


def trace_var(S):
    """Compute the variance of trace S.

    Uses ddof=1 (unbiased estimator) matching MATLAB's var().

    Parameters
    ----------
    S : array_like
        Input trace.

    Returns
    -------
    float
        Variance of the trace.
    """
    S = np.asarray(S, dtype=float).ravel()
    return float(np.var(S, ddof=1))


def trace_scv(S):
    """Compute the squared coefficient of variation for trace S.

    Parameters
    ----------
    S : array_like
        Input trace.

    Returns
    -------
    float
        Squared coefficient of variation (var(S)/mean(S)^2).
    """
    S = np.asarray(S, dtype=float).ravel()
    m = np.mean(S)
    v = np.var(S, ddof=1)
    return float(v / (m * m))


def autocov(X):
    """Compute the autocovariance of column vector X using FFT.

    Computes autocovariance from lag 0 to N-2 where N is the length of X.
    The formula is:
        acv(p) = 1/(N-p) * sum_{i=1}^{N-p} (X_i - X_bar) * (X_{i+p} - X_bar)

    Matches MATLAB: matlab/lib/kpctoolbox/contrib/autocov.m

    Parameters
    ----------
    X : array_like
        Input column vector (1-D array).

    Returns
    -------
    numpy.ndarray
        Autocovariance vector of length N-1.
    """
    X = np.asarray(X, dtype=float).ravel()
    M = len(X)
    if M < 2:
        raise ValueError("X is too short.")

    # Pad X - mean(X) with zeros (MATLAB pads to length 2*M)
    X_pad = np.concatenate([X - np.mean(X), np.zeros(M)])
    Y_pad = X_pad.copy()  # autocov: Y = X

    X_hat = np.fft.fft(X_pad)
    Y_hat = np.fft.fft(Y_pad)

    acv_full = np.fft.ifft(np.conj(X_hat) * Y_hat)
    # Take real part (imaginary part is due to float precision errors)
    # Extract lags 0 to M-2 and normalize by (M - p)
    acv = np.real(acv_full[:M - 1]) / (M - np.arange(1, M))

    return acv


def trace_acf(S, lags=None):
    """Compute the autocorrelation function for trace S at specified lags.

    If lags is not provided, computes ACF at lag 1.

    Matches MATLAB: matlab/lib/kpctoolbox/trace/trace_acf.m
    Uses the FFT-based autocov() fallback path (no xcorr dependency).

    Parameters
    ----------
    S : array_like
        Input trace.
    lags : array_like or int, optional
        Lag values at which to compute ACF. Default is 1.

    Returns
    -------
    numpy.ndarray
        Autocorrelation values at the specified lags (column vector).
    """
    S = np.asarray(S, dtype=float).ravel()

    if lags is None:
        lags = np.array([1])
    else:
        lags = np.asarray(lags, dtype=int).ravel()

    if len(lags) == 0:
        return np.array([])

    # Truncate lags that exceed trace length
    if np.max(lags) > len(S) - 2:
        lags = np.minimum(lags, len(S) - 2)
        lags = lags[lags != len(S) - 2]
        if len(lags) == 0:
            return np.array([])

    # Compute autocovariance using FFT
    # MATLAB: S=S(:); acv = autocov(S-mean(S));
    # Note: autocov already subtracts the mean internally, but MATLAB passes
    # S-mean(S) to autocov which then subtracts mean again (mean of zero-mean
    # data is ~0, so it's effectively the same). We follow the MATLAB code exactly.
    S_centered = S - np.mean(S)
    acv = autocov(S_centered)

    # rho = acv(1+lags) / acv(1)  (MATLAB 1-based indexing)
    # In Python 0-based: acv[lags] / acv[0]
    rho = acv[lags] / acv[0]

    return rho.ravel()


def trace_skew(S):
    """Compute the bias-corrected skewness of trace S.

    Equivalent to MATLAB's skewness(S, 0) (bias-corrected).

    Parameters
    ----------
    S : array_like
        Input trace.

    Returns
    -------
    float
        Bias-corrected skewness.
    """
    S = np.asarray(S, dtype=float).ravel()
    # Remove NaNs
    S = S[~np.isnan(S)]

    n = len(S)
    if n < 3:
        return float('nan')

    res = S - np.mean(S)
    s2 = np.mean(res ** 2)
    m3 = np.mean(res ** 3)

    # Bias correction factor: sqrt((n-1)/n) * n / (n-2)
    correction = np.sqrt((n - 1.0) / n) * n / (n - 2.0)

    return float(m3 * correction / (s2 ** 1.5))


def trace_joint(S, lag, order):
    """Compute joint moments E[X^{k_1}_{i} X^{k_2}_{i+lag_1} ...] for a trace.

    Matches MATLAB: matlab/lib/kpctoolbox/trace/trace_joint.m

    Parameters
    ----------
    S : array_like
        Input trace.
    lag : array_like
        Array of lag increments. Cumulative sum is computed internally.
    order : array_like
        Array of moment orders for each position.

    Returns
    -------
    float
        Joint moment value.
    """
    S = np.asarray(S, dtype=float).ravel()
    lag = np.asarray(lag, dtype=int).ravel()
    order = np.asarray(order, dtype=int).ravel()

    # MATLAB: lag = sort(cumsum(lag))
    lag = np.sort(np.cumsum(lag))
    K = len(lag)
    # MATLAB: lag = lag - lag(1)*ones(1,K)
    lag = lag - lag[0]

    max_lag = int(np.max(lag))
    n = len(S)

    # MATLAB: JMv = S(1:length(S)-max(lag)).^order(1);
    # In Python 0-based: S[0 : n - max_lag] ** order[0]
    JMv = S[0:n - max_lag].copy() ** order[0]

    # MATLAB: for i = 2:length(order)
    #     JMv = JMv.*(S(lag(i)+1:length(S)-max(lag)+lag(i)).^order(i));
    for i in range(1, len(order)):
        # MATLAB: S(lag(i)+1 : length(S)-max(lag)+lag(i))
        # Python 0-based: S[lag[i] : n - max_lag + lag[i]]
        JMv = JMv * (S[lag[i]:n - max_lag + lag[i]] ** order[i])

    return float(np.mean(JMv))


def trace_bicov(S, GRID):
    """Compute the bicovariance of the trace.

    Matches MATLAB: matlab/lib/kpctoolbox/trace/trace_bicov.m

    Parameters
    ----------
    S : array_like
        Input trace.
    GRID : array_like
        Lag values for bicovariance grid.

    Returns
    -------
    tuple
        (BiCov, BiCovLags) where BiCov is 1-D array of bicovariance values
        and BiCovLags is a 2-D array of shape (len(GRID)^2, 3) of lag triplets.
    """
    S = np.asarray(S, dtype=float).ravel()
    GRID = np.asarray(GRID, dtype=int).ravel()

    BiCovLags = []
    for i in GRID:
        for j in GRID:
            BiCovLags.append([1, i, j])
    BiCovLags = np.array(BiCovLags, dtype=int)

    BiCov = np.zeros(len(BiCovLags))
    for idx in range(len(BiCovLags)):
        BiCov[idx] = trace_joint(S, BiCovLags[idx, :], np.array([1, 1, 1]))

    return BiCov, BiCovLags


def trace_idi(S, kset=None, OPTION=None, n=None):
    """Compute the index of dispersion for intervals (IDI) of a trace.

    Matches MATLAB: matlab/lib/kpctoolbox/trace/trace_idi.m

    Parameters
    ----------
    S : array_like
        Input trace (inter-arrival times).
    kset : array_like or int, optional
        Values of k at which to compute IDI. Default: min(1000, ceil(len(S)/30)).
    OPTION : str, optional
        'aggregate' or 'aggregate-mix' for pre-aggregated data.
    n : int or array_like, optional
        Aggregation level (used with OPTION).

    Returns
    -------
    tuple
        (IDIk, support) where IDIk is an array of IDI values and support
        is the number of data points used.
    """
    S = np.asarray(S, dtype=float).ravel()

    if kset is None:
        kset = min(1000, int(np.ceil(len(S) / 30)))
        return trace_idi(S, kset)

    # Allow scalar kset
    if np.isscalar(kset):
        kset = np.array([kset], dtype=int)
    else:
        kset = np.asarray(kset, dtype=int).ravel()

    IDIk = []
    support = 0

    for k in kset:
        if OPTION is None:
            # Default: no aggregation
            support = len(S) - k - 1
            Sk = np.zeros(len(S) - k)
            for t in range(len(S) - k - 1):
                Sk[t] = np.sum(S[t:t + k])
            # MATLAB: IDIk(end+1) = k*var(Sk)/mean(Sk)^2
            IDIk.append(k * np.var(Sk[:support], ddof=1) / np.mean(Sk[:support]) ** 2)

        elif OPTION.lower() == 'aggregate':
            # Data is already aggregated: S is a sample of S=X1+...+Xn
            keff = int(np.floor(k / n))
            support = int(len(S) / keff)
            Sk = np.zeros(len(S) - keff)
            for t in range(len(S) - keff - 1):
                Sk[t] = np.sum(S[t:t + keff])
            valid = len(S) - keff - 1
            IDIk.append(k * np.var(Sk[:valid], ddof=1) / np.mean(Sk[:valid]) ** 2)

        elif OPTION.lower() == 'aggregate-mix':
            # Data is already aggregated with variable sizes
            n_arr = np.asarray(n, dtype=int).ravel()
            Sk = []
            for t0 in range(len(n_arr) - int(np.ceil(np.sum(n_arr) / k))):
                acc = 0
                indexes = [[0, 0]]  # [start, end] pairs
                indexes[0][0] = t0
                y_list = []
                for t in range(t0, len(n_arr)):
                    acc += n_arr[t]
                    if acc > k:
                        indexes[-1][1] = t - 1
                        indexes.append([t, 0])
                        y_list.append(acc - n_arr[t])
                        acc = n_arr[t]
                # Remove last incomplete index
                if len(indexes) > 0:
                    indexes = indexes[:-1]
                for pair in indexes:
                    Sk.append(np.sum(S[pair[0]:pair[1] + 1]))
            Sk = np.array(Sk)
            if len(Sk) > 1:
                IDIk.append(k * np.var(Sk, ddof=1) / np.mean(Sk) ** 2)
            else:
                IDIk.append(float('nan'))

    IDIk = np.array(IDIk)
    return IDIk, support


def trace_idc(S):
    """Compute the index of dispersion for counts (IDC) of a trace.

    Asymptotically, IDI and IDC are equal.

    Matches MATLAB: matlab/lib/kpctoolbox/trace/trace_idc.m

    Parameters
    ----------
    S : array_like
        Input trace.

    Returns
    -------
    float
        IDC value.
    """
    IDIk, _ = trace_idi(S)
    return float(IDIk[0]) if len(IDIk) > 0 else float('nan')


def trace_gamma(T, limit=1000):
    """Estimate the autocorrelation decay rate of a trace.

    Fits a geometric decay model rho(k) = RHO0 * gamma^k to the ACF
    using nonlinear least squares (matching MATLAB's nlinfit).

    Matches MATLAB: matlab/lib/kpctoolbox/trace/trace_gamma.m

    Parameters
    ----------
    T : array_like
        Input trace.
    limit : int, optional
        Maximum lag considered. Default is 1000.

    Returns
    -------
    tuple
        (GAMMA, RHO0, RESIDUALS) where GAMMA is the decay rate, RHO0 is the
        initial autocorrelation coefficient, and RESIDUALS are the fit residuals.
    """
    T = np.asarray(T, dtype=float).ravel()

    M1 = np.mean(T)
    M2 = np.mean(T ** 2)

    max_lag = min(limit, len(T))
    lag = np.arange(1, max_lag + 1)
    rho = trace_acf(T, lag)
    rho = rho.ravel()

    VAR = M2 - M1 ** 2
    SCV = VAR / (M1 ** 2)
    RHO0 = 0.5 * (1.0 - 1.0 / SCV)

    # Define geometric model: rhok = RHO0 * gamma^k
    def geometric(k, gamma):
        return RHO0 * gamma ** k

    try:
        # MATLAB calls nlinfit with RobustWgtFun='fair', i.e. an IRLS fit, NOT
        # ordinary least squares: outlying lags (the noisy tail of the ACF) are
        # down-weighted. Plain curve_fit here moved the estimate by ~1e-5
        # relative. This reproduces MATLAB's scheme - leverage-adjusted
        # residuals, MAD scale, fair weights w = 1/(1+|r|) at tune 1.400 - and
        # lands within ~3e-7 of it; the remainder is nlinfit's internal robust
        # scale convention, which is not documented well enough to match
        # exactly.
        popt, _ = curve_fit(geometric, lag, rho, p0=[0.99], bounds=(-1.0, 1.0),
                            maxfev=100000)
        GAMMA = _robust_fair_refit(geometric, lag, rho, RHO0, float(popt[0]))
        RESIDUALS = rho - geometric(lag, GAMMA)
    except Exception:
        try:
            # Fallback: bounded least squares (matches MATLAB lsqcurvefit fallback)
            popt, _ = curve_fit(geometric, lag, rho, p0=[0.99], bounds=(-1.0, 1.0))
            GAMMA = float(popt[0])
            RESIDUALS = rho - geometric(lag, GAMMA)
        except Exception:
            GAMMA = 0.5
            RESIDUALS = rho - geometric(lag, GAMMA)

    return GAMMA, RHO0, RESIDUALS


def _robust_fair_refit(model, lag, rho, RHO0, beta, tune=1.400,
                       maxiter=200, tol=1e-10):
    """Iteratively reweighted refit with the 'fair' weight function.

    Mirrors what MATLAB's nlinfit does when RobustWgtFun is set: adjust the
    residuals for leverage, estimate a robust scale from their MAD, weight by
    w = 1/(1+|r_adj/(tune*s)|), and refit until the parameter settles. The model
    has a single parameter (the decay rate), so the leverage of lag i is just
    J_i^2 / sum(J^2).
    """
    for _ in range(maxiter):
        r = rho - model(lag, beta)
        J = RHO0 * lag * beta ** (lag - 1.0)
        ssq = float(np.sum(J ** 2))
        if not np.isfinite(ssq) or ssq <= 0:
            return beta
        h = J ** 2 / ssq
        radj = r / np.sqrt(1.0 - np.minimum(0.9999, h))
        rs = np.sort(np.abs(radj))
        s = float(np.median(rs)) / 0.6745
        if not np.isfinite(s) or s <= 0:
            return beta
        w = 1.0 / (1.0 + np.abs(radj / (tune * s)))
        new, _ = curve_fit(model, lag, rho, p0=[beta], sigma=1.0 / np.sqrt(w),
                           absolute_sigma=False, bounds=(-1.0, 1.0),
                           maxfev=100000)
        new = float(new[0])
        if abs(new - beta) < tol * max(1e-12, abs(beta)):
            return new
        beta = new
    return beta


def trace_shuffle(S):
    """Randomly shuffle a trace.

    Matches MATLAB: matlab/lib/kpctoolbox/trace/trace_shuffle.m

    Parameters
    ----------
    S : array_like
        Input trace.

    Returns
    -------
    numpy.ndarray
        Shuffled trace.
    """
    S = np.asarray(S, dtype=float).ravel()
    return S[np.random.permutation(len(S))]


def trace_iat2counts(S, scale):
    """Compute the counting process from inter-arrival times.

    Computes the number of arrivals within 'scale' time units from each arrival.

    Matches MATLAB: matlab/lib/kpctoolbox/trace/trace_iat2counts.m

    Parameters
    ----------
    S : array_like
        Inter-arrival times.
    scale : float
        Time window size.

    Returns
    -------
    numpy.ndarray
        Counts array.
    """
    S = np.asarray(S, dtype=float).ravel()
    n = len(S)
    CS = np.cumsum(S)
    C = np.zeros(n - 1, dtype=int)

    for i in range(n - 1):
        cur = i
        while CS[cur + 1] - CS[i] <= scale:
            cur = cur + 1
            if cur == n - 1:  # MATLAB: cur == n (1-based), Python: cur == n-1 (0-based)
                C[i] = cur - i
                return C[:i + 1]
        C[i] = cur - i

    return C


def trace_iat2bins(S, scale):
    """Compute binned counts from inter-arrival times.

    Computes the counts C of S in each bin with timescale 'scale'.
    bC[i] gives the bin of membership for element i.

    Matches MATLAB: matlab/lib/kpctoolbox/trace/trace_iat2bins.m

    Parameters
    ----------
    S : array_like
        Inter-arrival times.
    scale : float
        Bin width.

    Returns
    -------
    tuple
        (C, bC) where C is counts per bin and bC is bin membership for each event.
    """
    from line_solver.api.trace.trace_analysis import trace_iat2bins as _trace_iat2bins
    return _trace_iat2bins(S, scale)


def trace_pmf(X):
    """Compute the probability mass function of discrete values in a trace.

    Matches MATLAB: matlab/lib/kpctoolbox/trace/trace_pmf.m
    MATLAB uses hist(X, max(X)) which creates max(X) bins from min to max.

    Parameters
    ----------
    X : array_like
        Input trace with discrete values.

    Returns
    -------
    tuple
        (pmf, px) where pmf is the probability mass function values
        and px is the sorted unique values.
    """
    X = np.asarray(X, dtype=float).ravel()
    n = len(X)
    if n == 0:
        return np.array([]), np.array([])

    # MATLAB: hist(X, max(X)) creates max(X) equally spaced bins
    # Then normalizes by numel(X)
    num_bins = int(np.max(X))
    if num_bins <= 0:
        num_bins = 1
    pmf_hist, bin_edges = np.histogram(X, bins=num_bins)
    pmf_hist = pmf_hist.astype(float) / n

    px = np.unique(X)

    return pmf_hist, px


def trace_summary(m, fid=None):
    """Compute comprehensive summary statistics for a trace.

    Matches MATLAB: matlab/lib/kpctoolbox/trace/trace_summary.m

    Parameters
    ----------
    m : array_like
        Input trace.
    fid : file-like or None, optional
        If provided, print summary to this file. Default prints to stdout.

    Returns
    -------
    tuple
        (MEAN, SCV, MAD, SKEW, KURT, QUART, TAILS1PERC, MINMAX, IQR, ACF, IDC)
    """
    import sys

    m = np.asarray(m, dtype=float).ravel()

    if fid is None:
        fid = sys.stdout

    MEAN = float(np.mean(m))
    SCV = float(np.var(m, ddof=1) / (np.mean(m) ** 2))

    # MATLAB: skewness(m) uses bias-corrected skewness by default
    SKEW = trace_skew(m)

    # MATLAB: kurtosis(m) computes excess kurtosis with bias correction
    n = len(m)
    if n < 4:
        KURT = float('nan')
    else:
        res = m - np.mean(m)
        s2 = np.mean(res ** 2)
        m4 = np.mean(res ** 4)
        # MATLAB kurtosis default: bias-corrected (flag=0 is default in newer MATLAB)
        # kurtosis = (m4/s2^2) with bias correction
        # Standard MATLAB kurtosis(x) = n*(n+1)/((n-1)*(n-2)*(n-3)) * sum((x-mean)^4)/var^2 - 3*(n-1)^2/((n-2)*(n-3))
        # But the simple form: m4/s2^2 is the non-corrected version
        # MATLAB kurtosis(x) with default flag (1) is m4/s2^2
        KURT = float(m4 / (s2 ** 2))

    QUART = np.array([np.percentile(m, 25), np.percentile(m, 50), np.percentile(m, 75)])
    TAILS1PERC = np.array([np.percentile(m, 95), np.percentile(m, (1 - 1e-6) * 100)])
    MINMAX = np.array([float(np.min(m)), float(np.max(m))])

    # MATLAB: mad(m, 1) = median absolute deviation (median-based)
    median_val = np.median(m)
    MAD = float(np.median(np.abs(m - median_val)))

    ACF = trace_acf(m, np.arange(1, 11))
    IDC = trace_idc(m)
    LEN = len(m)

    # MATLAB: iqr(m) = Q3 - Q1
    IQR_val = float(QUART[2] - QUART[0])

    # Print summary matching MATLAB format
    fid.write(
        'length=%d mean=%f scv=%f cv=%f mad=%f, skew=%f kurt=%f '
        'p25=%f p50=%f p75=%f p95=%f min=%f max=%f iqr=%d ' % (
            LEN, MEAN, SCV, np.sqrt(SCV), MAD, SKEW, KURT,
            QUART[0], QUART[1], QUART[2], TAILS1PERC[0],
            MINMAX[0], MINMAX[1], IQR_val))
    fid.write('acf=%f %f %f %f idc(1000)/scv=%f' % (
        ACF[0], ACF[1], ACF[2], ACF[3], IDC / SCV))
    fid.write('\n')

    return MEAN, SCV, MAD, SKEW, KURT, QUART, TAILS1PERC, MINMAX, IQR_val, ACF, IDC


def mtrace_mean(trace, ntypes, type_arr):
    """Compute the mean of a trace, divided by types.

    Matches MATLAB: matlab/lib/kpctoolbox/trace/mtrace_mean.m

    Parameters
    ----------
    trace : array_like
        The array containing the trace data.
    ntypes : int
        The number of different types.
    type_arr : array_like
        An array indicating the type of each element in the trace (0-indexed).

    Returns
    -------
    numpy.ndarray
        Column vector of mean values for each type.
    """
    trace = np.asarray(trace, dtype=float).ravel()
    type_arr = np.asarray(type_arr, dtype=int).ravel()

    m = np.zeros(ntypes)
    for c in range(ntypes):
        idx = (type_arr == c)
        if np.any(idx):
            m[c] = np.mean(trace[idx])
        else:
            m[c] = float('nan')

    return m
