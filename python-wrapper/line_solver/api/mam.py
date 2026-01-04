
"""
Matrix-Analytic Methods (MAM) for MAP/PH distributions.

This module provides functions for analyzing Markovian Arrival Processes (MAPs),
Phase-Type (PH) distributions, and related matrix-analytic methods. It includes
fitting algorithms, moment calculations, and various transformations.

Key function categories:
- MAP analysis: map_pie, map_mean, map_var, map_scv, map_skew
- MAP fitting: map2_fit, mmpp2_fit, aph_fit, aph2_fit
- PH distributions: Phase-type analysis and fitting
- Transformations: map_scale, map_normalize, map_timereverse
- QBD methods: qbd_R, qbd_mapmap1, qbd_raprap1
- Compression: compress_adaptive, compress_spectral

These functions support advanced stochastic modeling with correlated
arrivals and non-exponential service times.
"""

import jpype
import numpy as np
from line_solver import jlineMatrixToArray, jlineMatrixFromArray


def map_pie(D0, D1=None):
    """
    Compute equilibrium distribution of the embedded discrete-time process.
    
    Calculates the stationary distribution of the discrete-time Markov chain
    embedded at departure instants for a Markovian Arrival Process (MAP).
    The embedded process has transition matrix P = (-D0)^(-1) * D1.
    
    Args:
        D0: Generator matrix for non-arrival transitions, or MAP as container {D0,D1}.
        D1: Generator matrix for arrival transitions (optional if D0 is container).
        
    Returns:
        numpy.ndarray: Equilibrium distribution of the discrete-time Markov chain
                      embedded at departure instants.
    """
    if hasattr(D0, 'get') and hasattr(D0, 'length'):
        if D0.length() == 2:
            return jlineMatrixToArray(
                jpype.JPackage('jline').api.mam.Map_pieKt.map_pie(D0)
            )
        else:
            raise ValueError("D0 container must have exactly 2 elements")

    elif D1 is None:
        if isinstance(D0, (list, np.ndarray)):
            D0_array = np.array(D0)
            if D0_array.ndim == 3 and D0_array.shape[0] == 2:
                D0_mat = jlineMatrixFromArray(D0_array[0])
                D1_mat = jlineMatrixFromArray(D0_array[1])
                return jlineMatrixToArray(
                    jpype.JPackage('jline').api.mam.Map_pieKt.map_pie(D0_mat, D1_mat)
                )
            else:
                raise ValueError("When D1 is None, D0 must be a 3D array with shape (2, n, n)")
        else:
            raise ValueError("Invalid input type for D0 when D1 is None")

    else:
        return jlineMatrixToArray(
            jpype.JPackage('jline').api.mam.Map_pieKt.map_pie(
                jlineMatrixFromArray(D0),
                jlineMatrixFromArray(D1)
            )
        )


def map_mean(D0, D1=None):
    """
    Calculate mean inter-arrival time of a MAP.

    Computes the expected inter-arrival time for a Markovian Arrival Process.

    Args:
        D0: Generator matrix for non-arrival transitions, or MAP container
        D1: Generator matrix for arrival transitions (optional if D0 is container)

    Returns:
        float: Mean inter-arrival time
    """
    if hasattr(D0, 'get') and hasattr(D0, 'length'):
        if D0.length() == 2:
            return jpype.JPackage('jline').api.mam.Map_meanKt.map_mean(D0)
        else:
            raise ValueError("D0 container must have exactly 2 elements")

    elif D1 is None:
        if isinstance(D0, (list, np.ndarray)):
            D0_array = np.array(D0)
            if D0_array.ndim == 3 and D0_array.shape[0] == 2:
                D0_mat = jlineMatrixFromArray(D0_array[0])
                D1_mat = jlineMatrixFromArray(D0_array[1])
                return jpype.JPackage('jline').api.mam.Map_meanKt.map_mean(D0_mat, D1_mat)
            else:
                raise ValueError("When D1 is None, D0 must be a 3D array with shape (2, n, n)")
        else:
            raise ValueError("Invalid input type for D0 when D1 is None")

    else:
        return jpype.JPackage('jline').api.mam.Map_meanKt.map_mean(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1)
        )


def map_var(D0, D1=None):
    """
    Calculate variance of inter-arrival times of a MAP.

    Computes the variance of inter-arrival times for a Markovian Arrival Process.

    Args:
        D0: Generator matrix for non-arrival transitions, or MAP container
        D1: Generator matrix for arrival transitions (optional if D0 is container)

    Returns:
        float: Variance of inter-arrival times
    """
    if D1 is None:
        return jpype.JPackage('jline').api.mam.Map_varKt.map_var(D0)
    else:
        return jpype.JPackage('jline').api.mam.Map_varKt.map_var(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1)
        )


def map_scv(D0, D1=None):
    """
    Calculate squared coefficient of variation of a MAP.

    Computes the SCV (variance/mean²) of inter-arrival times.

    Args:
        D0: Generator matrix for non-arrival transitions, or MAP container
        D1: Generator matrix for arrival transitions (optional if D0 is container)

    Returns:
        float: Squared coefficient of variation
    """
    if D1 is None:
        return jpype.JPackage('jline').api.mam.Map_scvKt.map_scv(D0)
    else:
        return jpype.JPackage('jline').api.mam.Map_scvKt.map_scv(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1)
        )


def map_skew(D0, D1=None):
    """
    Calculate skewness of inter-arrival times of a MAP.

    Computes the third central moment normalized by variance^(3/2).

    Args:
        D0: Generator matrix for non-arrival transitions, or MAP container
        D1: Generator matrix for arrival transitions (optional if D0 is container)

    Returns:
        float: Skewness of inter-arrival times
    """
    if D1 is None:
        return jpype.JPackage('jline').api.mam.Map_skewKt.map_skew(D0)
    else:
        return jpype.JPackage('jline').api.mam.Map_skewKt.map_skew(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1)
        )


def map_moment(D0, D1, order):
    """
    Compute (power) moments of interarrival times of the specified order.
    
    Calculates the k-th moment E[X^k] where X is the interarrival time
    random variable of the MAP.

    Args:
        D0: Generator matrix for non-arrival transitions
        D1: Generator matrix for arrival transitions  
        order: Moment order to calculate (1=>E[X], 2=>E[X^2], ...)

    Returns:
        float: The k-th power moment of interarrival times
        
    Examples:
        - map_moment(D0, D1, 1) returns the first moment E[X]
        - map_moment(D0, D1, 2) returns the second moment E[X^2]
    """
    return jpype.JPackage('jline').api.mam.Map_momentKt.map_moment(
        jlineMatrixFromArray(D0),
        jlineMatrixFromArray(D1),
        jpype.JInt(order)
    )


def map_lambda(D0, D1=None):
    """
    Calculate arrival rate of a MAP.

    Computes the long-run arrival rate (λ) of the Markovian Arrival Process.

    Args:
        D0: Generator matrix for non-arrival transitions, or MAP container
        D1: Generator matrix for arrival transitions (optional if D0 is container)

    Returns:
        float: Arrival rate λ
    """
    if D1 is None:
        return jpype.JPackage('jline').api.mam.Map_lambdaKt.map_lambda(D0)
    else:
        return jpype.JPackage('jline').api.mam.Map_lambdaKt.map_lambda(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1)
        )


def map_acf(D0, D1, lags):
    """
    Compute autocorrelation coefficients of interarrival times.
    
    Calculates the lag-k autocorrelation coefficients of the interarrival
    times for a Markovian Arrival Process. The autocorrelation coefficient
    at lag k measures the correlation between interarrival times that are
    k arrivals apart.

    Args:
        D0: Generator matrix for non-arrival transitions
        D1: Generator matrix for arrival transitions
        lags: Array of positive integers specifying the lags of the
              autocorrelation coefficients to compute

    Returns:
        numpy.ndarray: Autocorrelation coefficients returned in the same
                      order as the lags vector
                      
    Examples:
        - map_acf(D0, D1, [1]) returns the lag-1 autocorrelation coefficient
        - map_acf(D0, D1, [1, 2, 3, 4, 5]) returns the first five autocorrelation coefficients
        - map_acf(D0, D1, np.logspace(0, 4, 5)) returns five logarithmically spaced
          autocorrelation coefficients in [1e0, 1e4]
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Map_acfKt.map_acf(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1),
            jlineMatrixFromArray(lags)
        )
    )


def map_acfc(D0, D1, lags, u):
    """
    Compute autocorrelation of a MAP's counting process at specified lags.
    
    Calculates the autocorrelation function of the counting process N(t)
    representing the number of arrivals in time interval [0,t], evaluated
    at discrete time points separated by timeslot length u.

    Args:
        D0: Generator matrix for non-arrival transitions
        D1: Generator matrix for arrival transitions
        lags: Array of lag values (positive integers)
        u: Length of timeslot (timescale parameter)

    Returns:
        numpy.ndarray: Autocorrelation values of the counting process
                      at the specified lags
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Map_acfcKt.map_acfc(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1),
            jlineMatrixFromArray(lags),
            jpype.JDouble(u)
        )
    )


def map_idc(D0, D1=None):
    """
    Compute the asymptotic index of dispersion.
    
    Calculates I = SCV(1 + 2*sum_{k=1}^∞ ρ_k) where SCV is the squared
    coefficient of variation and ρ_k is the lag-k autocorrelation coefficient
    of inter-arrival times. I is also the limiting value of the index of
    dispersion for counts and of the index of dispersion for intervals.

    Args:
        D0: Generator matrix for non-arrival transitions, or MAP container {D0,D1}
        D1: Generator matrix for arrival transitions (optional if D0 is container)

    Returns:
        float: Asymptotic index of dispersion
        
    Examples:
        - For a renewal process: map_idc(map_renewal(MAP)) equals the SCV of the MAP
    """
    if D1 is None:
        return jpype.JPackage('jline').api.mam.Map_idcKt.map_idc(D0)
    else:
        return jpype.JPackage('jline').api.mam.Map_idcKt.map_idc(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1)
        )


def map_gamma(D0, D1=None):
    """
    Estimate the auto-correlation decay rate of a MAP.
    
    For MAPs of order higher than 2, performs an approximation of the 
    autocorrelation function (ACF) curve using non-linear least squares fitting.
    For second-order MAPs, returns the geometric ACF decay rate. For Poisson
    processes (order 1), returns 0.

    Args:
        D0: Generator matrix for non-arrival transitions, or MAP container {D0,D1}
        D1: Generator matrix for arrival transitions (optional if D0 is container)

    Returns:
        float: Autocorrelation decay rate (GAMMA parameter)
    """
    if D1 is None:
        return jpype.JPackage('jline').api.mam.Map_gammaKt.map_gamma(D0)
    else:
        return jpype.JPackage('jline').api.mam.Map_gammaKt.map_gamma(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1)
        )


def map_gamma2(MAP):
    """
    Calculate the second largest eigenvalue of the embedded transition matrix.
    
    Computes the second-order gamma parameter (γ₂) which is the second largest
    eigenvalue (in absolute value) of the embedded discrete-time transition matrix
    P = (-D0)^(-1) * D1. This parameter characterizes the autocorrelation structure
    of the MAP beyond the dominant eigenvalue.

    Args:
        MAP: MAP structure or container {D0, D1}

    Returns:
        complex: Second largest eigenvalue of the embedded transition matrix
                (complex number, but typically real for feasible MAPs)
                
    Note:
        The eigenvalues are sorted by descending absolute value, where the first
        eigenvalue is 1 (for irreducible MAPs) and γ₂ is the second eigenvalue.
        This parameter is important for analyzing higher-order correlation properties.
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Map_gamma2Kt.map_gamma2(MAP)
    )


def map_cdf(D0, D1, points):
    """
    Calculate cumulative distribution function of a MAP.

    Evaluates the CDF of inter-arrival times at specified points.

    Args:
        D0: Generator matrix for non-arrival transitions
        D1: Generator matrix for arrival transitions
        points: Array of points to evaluate CDF at

    Returns:
        numpy.ndarray: CDF values at specified points
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Map_cdfKt.map_cdf(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1),
            jlineMatrixFromArray(points)
        )
    )


def map_piq(D0, D1=None):
    """
    Compute the probability in queue (piq) for a MAP.

    Args:
        D0: Generator matrix for non-arrival transitions, or MAP container
        D1: Generator matrix for arrival transitions (optional if D0 is container)

    Returns:
        numpy.ndarray: Probability in queue vector
    """
    if D1 is None:
        return jlineMatrixToArray(
            jpype.JPackage('jline').api.mam.Map_piqKt.map_piq(D0)
        )
    else:
        return jlineMatrixToArray(
            jpype.JPackage('jline').api.mam.Map_piqKt.map_piq(
                jlineMatrixFromArray(D0),
                jlineMatrixFromArray(D1)
            )
        )


def map_embedded(D0, D1=None):
    """
    Compute the embedded discrete-time process transition matrix.
    
    Returns the transition matrix P of the discrete-time Markov chain embedded
    at arrival instants. If the MAP is feasible, then P must be an irreducible
    stochastic matrix. P(i,j) gives the probability that the MAP restarts in
    phase j if the last arrival occurred in phase i.

    Args:
        D0: Generator matrix for non-arrival transitions, or MAP container {D0,D1}
        D1: Generator matrix for arrival transitions (optional if D0 is container)

    Returns:
        numpy.ndarray: Transition matrix P = (-D0)^(-1) * D1 of the embedded process
    """
    if D1 is None:
        return jlineMatrixToArray(
            jpype.JPackage('jline').api.mam.Map_embeddedKt.map_embedded(D0)
        )
    else:
        return jlineMatrixToArray(
            jpype.JPackage('jline').api.mam.Map_embeddedKt.map_embedded(
                jlineMatrixFromArray(D0),
                jlineMatrixFromArray(D1)
            )
        )


def map_count_mean(MAP, t):
    """
    Compute the mean of the counting process.
    
    Calculates the expected number of arrivals in time interval (0, t] for
    a Markovian Arrival Process.

    Args:
        MAP: Markovian Arrival Process as container {D0, D1}
        t: The time period considered for each sample of the counting process

    Returns:
        float: Mean number of arrivals in (0, t]
    """
    return jpype.JPackage('jline').api.mam.Map_count_meanKt.map_count_mean(
        MAP, jpype.JDouble(t)
    )


def map_count_var(MAP, t):
    """
    Compute the variance of the counting process.
    
    Calculates the variance of the number of arrivals in time interval (0, t]
    for a Markovian Arrival Process.

    Args:
        MAP: Markovian Arrival Process as container {D0, D1}
        t: The time period considered for each sample of the counting process

    Returns:
        float: Variance of number of arrivals in (0, t]
        
    Reference:
        He and Neuts, "Markov chains with marked transitions", 1998
        Verified special case of MMPP(2) with Andresen and Nielsen, 1998
    """
    return jpype.JPackage('jline').api.mam.Map_count_varKt.map_count_var(
        MAP, jpype.JDouble(t)
    )


def map_varcount(D0, D1, t):
    """
    Compute variance of the counting process for a MAP.
    
    Calculates the variance of the number of arrivals in time interval [0,t]
    for a Markovian Arrival Process. This is an alternative interface for
    map_count_var with direct matrix input.

    Args:
        D0: Generator matrix for non-arrival transitions
        D1: Generator matrix for arrival transitions
        t: Time interval length

    Returns:
        float: Variance of number of arrivals in time interval [0,t]
        
    Note:
        This function provides the same computation as map_count_var but
        accepts the D0, D1 matrices directly rather than a MAP container.
    """
    return jpype.JPackage('jline').api.mam.Map_varcountKt.map_varcount(
        jlineMatrixFromArray(D0),
        jlineMatrixFromArray(D1),
        jpype.JDouble(t)
    )




def map2_fit(e1, e2, e3, g2):
    """
    Fit a second-order MAP from moments and autocorrelation.
    
    Constructs a Markovian Arrival Process of order 2 from given statistical
    moments and lag-1 autocorrelation coefficient using the explicit inverse
    characterization method for acyclic MAPs of second order.

    Args:
        e1: First moment E[X] (mean inter-arrival time)
        e2: Second moment E[X²] 
        e3: Third moment E[X³] (if not provided, automatically selected)
        g2: Lag-1 autocorrelation coefficient

    Returns:
        tuple: (D0, D1) - Generator matrices for the fitted MAP or None if infeasible
        
    Reference:
        A. Heindl, G. Horvath, K. Gross "Explicit inverse characterization of
        acyclic MAPs of second order"
        
    Note:
        The method automatically handles hyperexponential (h2>0) and hypoexponential
        (h2<0) cases, where h2 = (r2-r1²)/r1² with r1=e1, r2=e2/2.
    """
    result = jpype.JPackage('jline').api.mam.Map2_fitKt.map2_fit(
        jpype.JDouble(e1),
        jpype.JDouble(e2),
        jpype.JDouble(e3),
        jpype.JDouble(g2)
    )

    D0 = jlineMatrixToArray(result.get(0))
    D1 = jlineMatrixToArray(result.get(1))

    return D0, D1


def aph_fit(e1, e2, e3, nmax=10):
    """
    Fit an Acyclic Phase-Type distribution to moments.

    Fits an APH distribution to match the first three moments using
    the EC (Expectation-Conditional) algorithm.

    Args:
        e1: First moment (mean)
        e2: Second moment
        e3: Third moment
        nmax: Maximum number of phases (default: 10)

    Returns:
        tuple: (alpha, T, feasible) where:
            - alpha: Initial probability vector
            - T: Sub-generator matrix
            - feasible: Whether the fit was feasible
    """
    result = jpype.JPackage('jline').api.mam.Aph_fitKt.aph_fit(
        jpype.JDouble(e1),
        jpype.JDouble(e2),
        jpype.JDouble(e3),
        jpype.JInt(nmax)
    )

    alpha = jlineMatrixToArray(result.alpha)
    T = jlineMatrixToArray(result.T)
    feasible = result.feasible

    return alpha, T, feasible


def aph2_fit(M1, M2, M3):
    """
    Fit a 2-phase APH distribution to moments.

    Fits a 2-phase Acyclic Phase-Type distribution to match
    the first three moments.

    Args:
        M1: First moment (mean)
        M2: Second moment
        M3: Third moment

    Returns:
        tuple: (alpha, T, feasible) where:
            - alpha: Initial probability vector
            - T: Sub-generator matrix
            - feasible: Whether the fit was feasible
    """
    result = jpype.JPackage('jline').api.mam.Aph2_fitKt.aph2_fit(
        jpype.JDouble(M1),
        jpype.JDouble(M2),
        jpype.JDouble(M3)
    )

    alpha = jlineMatrixToArray(result.alpha)
    T = jlineMatrixToArray(result.T)
    feasible = result.feasible

    return alpha, T, feasible


def aph2_fitall(M1, M2, M3):
    """
    Fit a 2-phase APH distribution with all possible structures.

    Tries all possible 2-phase APH structures to fit the given moments.

    Args:
        M1: First moment (mean)
        M2: Second moment
        M3: Third moment

    Returns:
        Results for all feasible APH structures
    """
    result = jpype.JPackage('jline').api.mam.Aph2_fitallKt.aph2_fitall(
        jpype.JDouble(M1),
        jpype.JDouble(M2),
        jpype.JDouble(M3)
    )

    fits = []
    for i in range(result.size()):
        fit = result.get(i)
        alpha = jlineMatrixToArray(fit.alpha)
        T = jlineMatrixToArray(fit.T)
        fits.append((alpha, T))

    return fits


def aph2_adjust(M1, M2, M3, method="default"):
    """
    Adjust moments to ensure feasible APH fitting.

    Adjusts the input moments to be in the feasible region for
    2-phase APH distribution fitting.

    Args:
        M1: First moment (mean)
        M2: Second moment
        M3: Third moment
        method: Adjustment method (default: "default")

    Returns:
        tuple: (M1, M2, M3) - Adjusted moments
    """
    result = jpype.JPackage('jline').api.mam.Aph2_adjustKt.aph2_adjust(
        jpype.JDouble(M1),
        jpype.JDouble(M2),
        jpype.JDouble(M3),
        jpype.JString(method)
    )

    return result.M1, result.M2, result.M3


def mmpp2_fit(E1, E2, E3, G2):
    """
    Fit a 2-state MMPP (Markov Modulated Poisson Process) to moments.

    Fits a 2-state MMPP to match the first three moments and
    lag-1 autocorrelation of inter-arrival times.

    Args:
        E1: First moment (mean)
        E2: Second moment
        E3: Third moment
        G2: Lag-1 autocorrelation coefficient

    Returns:
        tuple: (D0, D1, feasible) where:
            - D0: Generator matrix for non-arrival transitions
            - D1: Generator matrix for arrival transitions
            - feasible: Whether the fit was successful
    """
    result = jpype.JPackage('jline').api.mam.Mmpp2_fitKt.mmpp2_fit(
        jpype.JDouble(E1),
        jpype.JDouble(E2),
        jpype.JDouble(E3),
        jpype.JDouble(G2)
    )

    D0 = jlineMatrixToArray(result.D0)
    D1 = jlineMatrixToArray(result.D1)
    feasible = result.feasible

    return D0, D1, feasible


def mmpp2_fit1(mean, scv, skew, idc):
    """
    Fit a 2-state MMPP using alternative parameterization.

    Fits a 2-state MMPP using mean, SCV, skewness, and index of dispersion.

    Args:
        mean: Mean inter-arrival time
        scv: Squared coefficient of variation
        skew: Skewness
        idc: Index of dispersion for counts

    Returns:
        tuple: (D0, D1, feasible) where:
            - D0: Generator matrix for non-arrival transitions
            - D1: Generator matrix for arrival transitions
            - feasible: Whether the fit was successful
    """
    result = jpype.JPackage('jline').api.mam.Mmpp2_fit1Kt.mmpp2_fit1(
        jpype.JDouble(mean),
        jpype.JDouble(scv),
        jpype.JDouble(skew),
        jpype.JDouble(idc)
    )

    D0 = jlineMatrixToArray(result.D0)
    D1 = jlineMatrixToArray(result.D1)
    feasible = result.feasible

    return D0, D1, feasible


def mmap_mixture_fit(P2, M1, M2, M3):
    """
    Fit a mixture of MMAPs to given moments.

    Args:
        P2: Mixing probabilities matrix
        M1: First moment matrix
        M2: Second moment matrix
        M3: Third moment matrix

    Returns:
        Fitted MMAP mixture parameters
    """
    result = jpype.JPackage('jline').api.mam.Mmap_mixture_fitKt.mmap_mixture_fit(
        jlineMatrixFromArray(P2),
        jlineMatrixFromArray(M1),
        jlineMatrixFromArray(M2),
        jlineMatrixFromArray(M3)
    )

    D0 = jlineMatrixToArray(result.D0)
    D1 = jlineMatrixToArray(result.D1)
    feasible = result.feasible

    return D0, D1, feasible


def mmap_mixture_fit_mmap(mmap):
    """
    Fit mixture of MMAPs from MMAP representation.

    Args:
        mmap: MMAP object or structure to fit mixture from

    Returns:
        tuple: (D0, D1, components) where:
            - D0: Generator matrix for non-arrival transitions
            - D1: Generator matrix for arrival transitions  
            - components: List of mixture components
    """
    result = jpype.JPackage('jline').api.mam.Mmap_mixture_fit_mmapKt.mmap_mixture_fit_mmap(mmap)

    D0 = jlineMatrixToArray(result.D0)
    D1 = jlineMatrixToArray(result.D1)

    components = []
    for i in range(result.components.size()):
        component = result.components.get(i)
        components.append(jlineMatrixToArray(component))

    return D0, D1, components


def mamap2m_fit_gamma_fb_mmap(mmap):
    """
    Fit 2-moment MAMAP using gamma forward-backward algorithm from MMAP.

    Args:
        mmap: MMAP structure to fit

    Returns:
        Fitted MAMAP parameters
    """
    result = jpype.JPackage('jline').api.mam.Mamap2m_fit_gamma_fb_mmapKt.mamap2m_fit_gamma_fb_mmap(mmap)

    D0 = jlineMatrixToArray(result.get(0))
    D1 = jlineMatrixToArray(result.get(1))

    D_classes = []
    for i in range(2, result.length()):
        D_classes.append(jlineMatrixToArray(result.get(i)))

    return D0, D1, D_classes


def mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B):
    """
    Fit 2-moment MAMAP using gamma forward-backward algorithm.

    Args:
        M1: First moment
        M2: Second moment
        M3: Third moment
        GAMMA: Gamma parameter
        P: Probability parameters
        F: Forward parameters
        B: Backward parameters

    Returns:
        tuple: (D0, D1, D_classes) - Fitted MAMAP matrices and class matrices
    """
    import numpy as np

    P_java = jpype.JArray(jpype.JDouble)(P)
    F_java = jpype.JArray(jpype.JDouble)(F)
    B_java = jpype.JArray(jpype.JDouble)(B)

    result = jpype.JPackage('jline').api.mam.Mamap2m_fit_gamma_fb_mmapKt.mamap2m_fit_gamma_fb(
        jpype.JDouble(M1), jpype.JDouble(M2), jpype.JDouble(M3), jpype.JDouble(GAMMA),
        P_java, F_java, B_java
    )

    D0 = jlineMatrixToArray(result.get(0))
    D1 = jlineMatrixToArray(result.get(1))

    D_classes = []
    for i in range(2, result.length()):
        D_classes.append(jlineMatrixToArray(result.get(i)))

    return D0, D1, D_classes



def map_exponential(mean):
    """
    Fit a Poisson process as a MAP.
    
    Creates a Markovian Arrival Process representing a Poisson process
    (exponential inter-arrival times) with the specified mean.

    Args:
        mean: Mean inter-arrival time of the process

    Returns:
        tuple: (D0, D1) - MAP matrices representing the Poisson process
        
    Examples:
        - map_exponential(2) returns a Poisson process with rate λ=0.5
    """
    result = jpype.JPackage('jline').api.mam.Map_exponentialKt.map_exponential(
        jpype.JDouble(mean)
    )

    D0 = jlineMatrixToArray(result.get(0))
    D1 = jlineMatrixToArray(result.get(1))

    return D0, D1


def map_erlang(mean, k):
    """
    Fit an Erlang-k process as a MAP.
    
    Creates a Markovian Arrival Process representing an Erlang-k distribution
    with the specified mean and number of phases.

    Args:
        mean: Mean inter-arrival time of the process
        k: Number of phases in the Erlang distribution

    Returns:
        tuple: (D0, D1) - MAP matrices representing the Erlang-k process
        
    Examples:
        - map_erlang(2, 3) creates an Erlang-3 process with mean 2
    """
    result = jpype.JPackage('jline').api.mam.Map_erlangKt.map_erlang(
        jpype.JDouble(mean),
        jpype.JInt(k)
    )

    D0 = jlineMatrixToArray(result.get(0))
    D1 = jlineMatrixToArray(result.get(1))

    return D0, D1


def map_hyperexp(mean, scv, p):
    """
    Fit a two-phase Hyperexponential process as a MAP.
    
    Creates a Markovian Arrival Process representing a two-phase hyperexponential
    distribution with specified mean, squared coefficient of variation, and phase
    selection probability.

    Args:
        mean: Mean inter-arrival time of the process
        scv: Squared coefficient of variation of inter-arrival times
        p: Probability of being served in phase 1 (default: p=0.99)

    Returns:
        tuple: (D0, D1) - MAP matrices representing the hyperexponential process
        
    Examples:
        - map_hyperexp(1, 2, 0.99) creates a two-phase hyperexponential process
          where phase 1 is selected with probability 0.99
        - map_hyperexp(1, 2, 0.2) creates a two-phase hyperexponential process
          where phase 1 is selected with probability 0.2
        
    Note:
        For some parameter combinations, a feasible solution may not exist.
        Use map_isfeasible() to check if the returned MAP is valid.
    """
    result = jpype.JPackage('jline').api.mam.Map_hyperexpKt.map_hyperexp(
        jpype.JDouble(mean),
        jpype.JDouble(scv),
        jlineMatrixFromArray(p)
    )

    D0 = jlineMatrixToArray(result.get(0))
    D1 = jlineMatrixToArray(result.get(1))

    return D0, D1


def map_scale(D0, D1, newMean):
    """
    Rescale mean inter-arrival time of a MAP.
    
    Returns a MAP with the same normalized moments and correlations as the
    input MAP, except for the mean inter-arrival time which is set to the
    specified new value.

    Args:
        D0: Generator matrix for non-arrival transitions
        D1: Generator matrix for arrival transitions
        newMean: New mean inter-arrival time

    Returns:
        tuple: (D0_scaled, D1_scaled) - Scaled and normalized MAP matrices
        
    Examples:
        - map_mean(map_scale(map_exponential(1), 2)) has mean 2
    """
    result = jpype.JPackage('jline').api.mam.Map_scaleKt.map_scale(
        jlineMatrixFromArray(D0),
        jlineMatrixFromArray(D1),
        jpype.JDouble(newMean)
    )

    D0_scaled = jlineMatrixToArray(result.get(0))
    D1_scaled = jlineMatrixToArray(result.get(1))

    return D0_scaled, D1_scaled


def map_normalize(D0, D1):
    """
    Try to make a MAP feasible through normalization.
    
    Normalizes the MAP matrices by: (1) ensuring D0+D1 is an infinitesimal
    generator, (2) setting all negative off-diagonal entries to zero, and
    (3) taking the real part of any complex values.

    Args:
        D0: Generator matrix for non-arrival transitions
        D1: Generator matrix for arrival transitions

    Returns:
        tuple: (D0_norm, D1_norm) - Normalized MAP matrices that form a valid MAP
        
    Examples:
        - map_normalize([[0,0],[0,0]], [[1,2],[3,4]]) produces a valid MAP
    """
    result = jpype.JPackage('jline').api.mam.Map_normalizeKt.map_normalize(
        jlineMatrixFromArray(D0),
        jlineMatrixFromArray(D1)
    )

    D0_norm = jlineMatrixToArray(result.get(0))
    D1_norm = jlineMatrixToArray(result.get(1))

    return D0_norm, D1_norm


def map_timereverse(map_input):
    """
    Compute time-reversed MAP.

    Args:
        map_input: Either tuple (D0, D1) or MAP container

    Returns:
        tuple: (D0_rev, D1_rev) - Time-reversed MAP matrices
    """
    if isinstance(map_input, tuple):
        D0, D1 = map_input
        result = jpype.JPackage('jline').api.mam.Map_timereverseKt.map_timereverse(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Map_timereverseKt.map_timereverse(map_input)

    D0_rev = jlineMatrixToArray(result.get(0))
    D1_rev = jlineMatrixToArray(result.get(1))

    return D0_rev, D1_rev


def map_mark(MAP, prob):
    """
    Mark arrivals from a MAP according to given probabilities.
    
    Takes a MAP with a single arrival type and returns a MMAP (Marked MAP)
    with the same unmarked inter-arrival process but marked arrivals
    specified by prob. prob[k] is the probability that an arrival
    will be marked with type k.

    Args:
        MAP: Either tuple (D0, D1) or MAP container representing the input MAP
        prob: Probability vector where prob[k] is the probability that arrivals
              are marked with type k (probabilities are normalized if they don't sum to 1)

    Returns:
        Marked MAP (MMAP) with multiple arrival types
        
    Note:
        Input marking probabilities are automatically normalized if they don't sum to 1.
    """
    if isinstance(MAP, tuple):
        D0, D1 = MAP
        result = jpype.JPackage('jline').api.mam.Map_markKt.map_mark(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1),
            jlineMatrixFromArray(prob)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Map_markKt.map_mark(
            MAP,
            jlineMatrixFromArray(prob)
        )

    D0 = jlineMatrixToArray(result.get(0))
    D1 = jlineMatrixToArray(result.get(1))

    return D0, D1


def map_infgen(D0, D1):
    """
    Compute infinitesimal generator matrix of MAP.

    Args:
        D0: Generator matrix for non-arrival transitions
        D1: Generator matrix for arrival transitions

    Returns:
        numpy.ndarray: Infinitesimal generator matrix
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Map_infgenKt.map_infgen(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1)
        )
    )




def map_super(MAPa, MAPb):
    """
    Compute superposition of two independent MAPs.

    Creates a MAP that represents the superposition (sum) of two
    independent Markovian Arrival Processes using Kronecker sum operations.
    The resulting MAP has arrivals from both input MAPs combined.

    Args:
        MAPa: First MAP as tuple (D0_a, D1_a) or container {D0_a, D1_a}
        MAPb: Second MAP as tuple (D0_b, D1_b) or container {D0_b, D1_b}

    Returns:
        tuple: (D0_super, D1_super) - Superposition MAP matrices where:
               D0_super = kron(D0_a, I_b) + kron(I_a, D0_b)
               D1_super = kron(D1_a, I_b) + kron(I_a, D1_b)
               
    Note:
        The output is automatically normalized to ensure a valid MAP.
    """
    if isinstance(MAPa, tuple) and isinstance(MAPb, tuple):
        D0a, D1a = MAPa
        D0b, D1b = MAPb
        result = jpype.JPackage('jline').api.mam.Map_superKt.map_super(
            jlineMatrixFromArray(D0a),
            jlineMatrixFromArray(D1a),
            jlineMatrixFromArray(D0b),
            jlineMatrixFromArray(D1b)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Map_superKt.map_super(MAPa, MAPb)

    D0_super = jlineMatrixToArray(result.get(0))
    D1_super = jlineMatrixToArray(result.get(1))

    return D0_super, D1_super


def map_sum(MAP, n):
    """
    Compute sum of n identical independent MAPs.

    Args:
        MAP: MAP as tuple (D0, D1) or container
        n: Number of identical MAPs to sum

    Returns:
        tuple: (D0_sum, D1_sum) - Sum MAP matrices
    """
    if isinstance(MAP, tuple):
        D0, D1 = MAP
        result = jpype.JPackage('jline').api.mam.Map_sumKt.map_sum(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1),
            jpype.JInt(n)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Map_sumKt.map_sum(MAP, jpype.JInt(n))

    D0_sum = jlineMatrixToArray(result.get(0))
    D1_sum = jlineMatrixToArray(result.get(1))

    return D0_sum, D1_sum


def map_sumind(MAPs):
    """
    Compute sum of multiple independent MAPs.

    Args:
        MAPs: List of MAPs, each as tuple (D0, D1) or container

    Returns:
        tuple: (D0_sum, D1_sum) - Sum MAP matrices
    """
    java_maps = jpype.java.util.ArrayList()
    for MAP in MAPs:
        if isinstance(MAP, tuple):
            D0, D1 = MAP
            map_matrices = jpype.java.util.ArrayList()
            map_matrices.add(jlineMatrixFromArray(D0))
            map_matrices.add(jlineMatrixFromArray(D1))
            java_maps.add(map_matrices)
        else:
            java_maps.add(MAP)

    result = jpype.JPackage('jline').api.mam.Map_sumindKt.map_sumind(java_maps)

    D0_sum = jlineMatrixToArray(result.get(0))
    D1_sum = jlineMatrixToArray(result.get(1))

    return D0_sum, D1_sum


def map_checkfeasible(MAP, TOL=1e-10):
    """
    Check if MAP satisfies feasibility conditions and return diagnostics.

    Args:
        MAP: MAP as tuple (D0, D1) or container
        TOL: Numerical tolerance (default: 1e-10)

    Returns:
        Feasibility check results with diagnostic information
    """
    if isinstance(MAP, tuple):
        D0, D1 = MAP
        return jpype.JPackage('jline').api.mam.Map_checkfeasibleKt.map_checkfeasible(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1),
            jpype.JDouble(TOL)
        )
    else:
        return jpype.JPackage('jline').api.mam.Map_checkfeasibleKt.map_checkfeasible(
            MAP, jpype.JDouble(TOL)
        )


def map_isfeasible(MAP, TOL=1e-10):
    """
    Evaluate feasibility of a MAP process.
    
    Checks if a Markovian Arrival Process satisfies all feasibility conditions
    including proper matrix structure, non-negativity constraints, stochasticity,
    and irreducibility requirements.

    Args:
        MAP: MAP as tuple (D0, D1) or container {D0, D1}
        TOL: Numerical tolerance for feasibility checks (default: 1e-10)

    Returns:
        bool: True if MAP is feasible, False if infeasible
        
    Examples:
        - map_isfeasible([[0,0],[0,0]], [[1,2],[3,4]]) returns False (infeasible MAP)
        
    Note:
        Numerical tolerance is based on the standard toolbox value in map_feastol().
    """
    if isinstance(MAP, tuple):
        D0, D1 = MAP
        return jpype.JPackage('jline').api.mam.Map_isfeasibleKt.map_isfeasible(
            jlineMatrixFromArray(D0),
            jlineMatrixFromArray(D1),
            jpype.JDouble(TOL)
        )
    else:
        return jpype.JPackage('jline').api.mam.Map_isfeasibleKt.map_isfeasible(
            MAP, jpype.JDouble(TOL)
        )


def map_feastol():
    """
    Get default feasibility tolerance for MAP validation.

    Returns:
        float: Default tolerance value
    """
    return jpype.JPackage('jline').api.mam.Map_feastolKt.map_feastol()


def map_largemap():
    """
    Get threshold for determining if a MAP is considered large.

    Returns:
        int: Size threshold for large MAPs
    """
    return jpype.JPackage('jline').api.mam.Map_largemapKt.map_largemap()


def aph2_assemble(l1, l2, p1):
    """
    Assemble a 2-phase APH distribution from parameters.

    Args:
        l1: Rate parameter for first phase
        l2: Rate parameter for second phase
        p1: Initial probability for first phase

    Returns:
        tuple: (alpha, T) - Initial vector and generator matrix
    """
    result = jpype.JPackage('jline').api.mam.Aph2_assembleKt.aph2_assemble(
        jpype.JDouble(l1),
        jpype.JDouble(l2),
        jpype.JDouble(p1)
    )

    alpha = jlineMatrixToArray(result.alpha)
    T = jlineMatrixToArray(result.T)

    return alpha, T


def ph_reindex(sn):
    """
    Reindex phase-type distributions using NetworkStruct.

    .. deprecated::
        This is an internal MAM implementation detail. Use solvers instead.

    Args:
        sn: NetworkStruct object containing phase-type distributions (sn.proc)

    Returns:
        Reindexed phase-type distributions as a map indexed by station and class indices
    """
    result = jpype.JPackage('jline').api.mam.Ph_reindexKt.ph_reindex(sn)
    return result


def map_rand(K):
    """
    Generate a random MAP with K states.

    Args:
        K: Number of states in the MAP

    Returns:
        tuple: (D0, D1) - Randomly generated MAP matrices
    """
    result = jpype.JPackage('jline').api.mam.Map_randKt.map_rand(jpype.JInt(K))

    D0 = jlineMatrixToArray(result.get(0))
    D1 = jlineMatrixToArray(result.get(1))

    return D0, D1


def map_randn(K, mu, sigma):
    """
    Generate a random MAP with normally distributed matrix elements.
    
    Creates a random Markovian Arrival Process with K states where the matrix
    elements are drawn from normal distributions with specified mean and standard
    deviation parameters. The resulting matrices are normalized to ensure feasibility.

    Args:
        K: Number of states in the MAP
        mu: Mean parameter for the normal distribution of matrix elements
        sigma: Standard deviation parameter for the normal distribution

    Returns:
        tuple: (D0, D1) - Random MAP matrices with normalized feasible structure
        
    Note:
        The generated matrices undergo normalization via map_normalize() to ensure
        they form a valid MAP structure with proper generator properties.
    """
    result = jpype.JPackage('jline').api.mam.Map_randnKt.map_randn(
        jpype.JInt(K),
        jpype.JDouble(mu),
        jpype.JDouble(sigma)
    )

    D0 = jlineMatrixToArray(result.get(0))
    D1 = jlineMatrixToArray(result.get(1))

    return D0, D1




def mmap_lambda(MMAP):
    """
    Calculate arrival rate for each class in a Marked MAP.
    
    Computes the arrival rate (lambda) for each job class in a Marked
    Markovian Arrival Process (MMAP). For a MMAP with C classes, returns
    a vector of C arrival rates.

    Args:
        MMAP: Marked MAP structure as container {D0, D1, D2, ..., DC+1}
              where D0 is the background generator and D1, D2, ..., DC+1
              are the arrival generators for each class

    Returns:
        numpy.ndarray: Vector of arrival rates for each job class
        
    Note:
        This function is equivalent to mmap_count_lambda(MMAP).
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Mmap_lambdaKt.mmap_lambda(MMAP)
    )


def mmap_count_mean(MMAP, t):
    """
    Compute the mean of the counting process for each class in a Marked MAP.
    
    Calculates the expected number of arrivals in time interval (0, t] for
    each job class in a Marked Markovian Arrival Process (MMAP).

    Args:
        MMAP: Marked MAP structure as container {D0, D1, D2, ..., DC+1}
        t: The time period considered for each sample of the counting process

    Returns:
        numpy.ndarray: Vector with the mean number of arrivals in (0, t] for each job class
        
    Note:
        For a MMAP with C classes, returns a C-dimensional vector where
        mk[k] is the mean arrival count for class k in interval (0, t].
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Mmap_count_meanKt.mmap_count_mean(
            MMAP, jpype.JDouble(t)
        )
    )


def mmap_count_var(MMAP, t):
    """
    Compute variance of number of arrivals per class in time interval t for an MMAP.

    Args:
        MMAP: Markovian Arrival Process with Multiple types
        t: Time interval

    Returns:
        numpy.ndarray: Variance of number of arrivals per class in time t
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Mmap_count_varKt.mmap_count_var(
            MMAP, jpype.JDouble(t)
        )
    )


def mmap_count_idc(MMAP, t):
    """
    Compute index of dispersion for counts per class in time interval t for an MMAP.

    Args:
        MMAP: Markovian Arrival Process with Multiple types
        t: Time interval

    Returns:
        numpy.ndarray: Index of dispersion for counts per class
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Mmap_count_idcKt.mmap_count_idc(
            MMAP, jpype.JDouble(t)
        )
    )


def mmap_idc(MMAP):
    """
    Compute index of dispersion for counts of an MMAP.

    Args:
        MMAP: Markovian Arrival Process with Multiple types

    Returns:
        float: Index of dispersion for counts
    """
    return jpype.JPackage('jline').api.mam.Mmap_idcKt.mmap_idc(MMAP)


def mmap_sigma2(mmap):
    """
    Compute second-order statistics (variance) for an MMAP.

    Args:
        mmap: Markovian Arrival Process with Multiple types

    Returns:
        numpy.ndarray: Second-order statistics matrix
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Mmap_sigma2Kt.mmap_sigma2(mmap)
    )


def mmap_exponential(lambda_rates, n):
    """
    Create an exponential MMAP with given arrival rates.

    Args:
        lambda_rates: Array of arrival rates for each class
        n: Number of states

    Returns:
        MMAP: Exponential MMAP with specified rates
    """
    return jpype.JPackage('jline').api.mam.Mmap_exponentialKt.mmap_exponential(
        jlineMatrixFromArray(lambda_rates),
        jpype.JInt(n)
    )


def mmap_mixture(alpha, MAPs):
    """
    Create a mixture of MMAPs with given mixing probabilities.

    Args:
        alpha: Mixing probability vector
        MAPs: List of MMAP components

    Returns:
        MMAP: Mixture MMAP
    """
    java_maps = jpype.java.util.ArrayList()
    for MAP in MAPs:
        java_maps.add(MAP)

    return jpype.JPackage('jline').api.mam.Mmap_mixtureKt.mmap_mixture(
        jlineMatrixFromArray(alpha),
        java_maps
    )


def mmap_super(MMAPa, MMAPb, opt="default"):
    """
    Compute superposition of two MMAPs.

    Args:
        MMAPa: First MMAP
        MMAPb: Second MMAP
        opt: Superposition method (default: "default")

    Returns:
        MMAP: Superposed MMAP
    """
    return jpype.JPackage('jline').api.mam.Mmap_superKt.mmap_super(
        MMAPa, MMAPb, jpype.JString(opt)
    )


def mmap_super_safe(MMAPS, maxorder=10, method="default"):
    """
    Safe superposition of multiple MMAPs with error checking.

    Performs superposition of multiple Marked Markovian Arrival Processes
    with additional validation and error handling.

    Args:
        MMAPS: List of MMAP matrices to superpose
        maxorder: Maximum order for approximation
        method: Superposition method

    Returns:
        Superposed MMAP matrices
    """
    java_mmaps = jpype.java.util.ArrayList()
    for MMAP in MMAPS:
        java_mmaps.add(MMAP)

    return jpype.JPackage('jline').api.mam.Mmap_super_safeKt.mmap_super_safe(
        java_mmaps,
        jpype.JInt(maxorder),
        jpype.JString(method)
    )



def mmap_compress(MMAP, config="default"):
    """
    Compress an MMAP using specified configuration.

    Reduces the number of states in an MMAP while preserving
    key statistical properties using the specified configuration.

    Args:
        MMAP: MMAP matrices to compress
        config: Compression configuration

    Returns:
        Compressed MMAP matrices
    """
    return jpype.JPackage('jline').api.mam.Mmap_compressKt.mmap_compress(
        MMAP, jpype.JString(config)
    )


def mmap_normalize(MMAP):
    """
    Fix MMAP feasibility by normalization.
    
    Normalizes a Marked MAP to ensure feasibility by setting negative values
    to zero and enforcing the proper conditions for a valid MMAP structure.
    This includes ensuring the background matrix D0 has negative diagonal
    elements and the arrival matrices are non-negative.

    Args:
        MMAP: Marked MAP structure as container {D0, D1, D2, ..., DC+1}

    Returns:
        Normalized MMAP structure that satisfies feasibility conditions
        
    Note:
        The normalization process:
        1. Sets negative off-diagonal elements in D0 to zero
        2. Sets negative elements in arrival matrices to zero  
        3. Adjusts diagonal elements of D0 to maintain generator property
    """
    return jpype.JPackage('jline').api.mam.Mmap_normalizeKt.mmap_normalize(MMAP)


def mmap_scale(MMAP, M, maxIter=100):
    """
    Scale an MMAP using iterative algorithm.

    Adjusts the MMAP parameters using matrix M through
    an iterative scaling procedure.

    Args:
        MMAP: MMAP matrices to scale
        M: Scaling matrix
        maxIter: Maximum number of iterations

    Returns:
        Scaled MMAP matrices
    """
    return jpype.JPackage('jline').api.mam.Mmap_scaleKt.mmap_scale(
        MMAP, jlineMatrixFromArray(M), jpype.JInt(maxIter)
    )


def mmap_timereverse(mmap):
    """
    Compute time-reversed MMAP.

    Creates the time-reversed version of a Marked Markovian
    Arrival Process, reversing the direction of time.

    Args:
        mmap: MMAP matrices to time-reverse

    Returns:
        Time-reversed MMAP matrices
    """
    return jpype.JPackage('jline').api.mam.Mmap_timereverseKt.mmap_timereverse(mmap)


def mmap_hide(MMAP, types):
    """
    Hide specific arrival types from an MMAP.

    Removes specified arrival types from the MMAP,
    effectively filtering out those arrival classes.

    Args:
        MMAP: MMAP matrices
        types: Matrix specifying types to hide

    Returns:
        MMAP with specified types hidden
    """
    return jpype.JPackage('jline').api.mam.Mmap_hideKt.mmap_hide(
        MMAP, jlineMatrixFromArray(types)
    )


def mmap_shorten(mmap):
    """
    Shorten an MMAP by removing redundant elements.

    Reduces the size of the MMAP representation by
    eliminating unnecessary components.

    Args:
        mmap: MMAP matrices to shorten

    Returns:
        Shortened MMAP matrices
    """
    return jpype.JPackage('jline').api.mam.Mmap_shortenKt.mmap_shorten(mmap)


def mmap_maps(MMAP):
    """
    Extract individual MAP components from an MMAP.

    Decomposes a Marked Markovian Arrival Process into
    its constituent MAP components.

    Args:
        MMAP: MMAP matrices to decompose

    Returns:
        list: List of individual MAP matrices
    """
    result = jpype.JPackage('jline').api.mam.Mmap_mapsKt.mmap_maps(MMAP)

    maps = []
    for i in range(result.size()):
        maps.append(result.get(i))

    return maps


def mmap_pc(MMAP):
    """
    Compute probability matrix for MMAP.

    Calculates the probability matrix associated with
    the Marked Markovian Arrival Process.

    Args:
        MMAP: MMAP matrices

    Returns:
        numpy.ndarray: Probability matrix
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Mmap_pcKt.mmap_pc(MMAP)
    )


def mmap_forward_moment(MMAP, ORDERS, NORM=True):
    """
    Compute forward moments for MMAP.

    Calculates forward moments of the inter-arrival times
    for a Marked Markovian Arrival Process.

    Args:
        MMAP: MMAP matrices
        ORDERS: Array of moment orders to compute
        NORM: Whether to normalize the moments

    Returns:
        numpy.ndarray: Forward moment matrix
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Mmap_forward_momentKt.mmap_forward_moment(
            MMAP, jlineMatrixFromArray(ORDERS), jpype.JBoolean(NORM)
        )
    )


def mmap_backward_moment(MMAP, ORDERS, NORM=True):
    """
    Compute backward moments for MMAP.

    Calculates backward moments of the inter-arrival times
    for a Marked Markovian Arrival Process.

    Args:
        MMAP: MMAP matrices
        ORDERS: Array of moment orders to compute
        NORM: Whether to normalize the moments

    Returns:
        numpy.ndarray: Backward moment matrix
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Mmap_backward_momentKt.mmap_backward_moment(
            MMAP, jlineMatrixFromArray(ORDERS), jpype.JBoolean(NORM)
        )
    )


def mmap_cross_moment(mmap, k):
    """
    Compute cross moments for MMAP.

    Calculates cross moments between different arrival types
    in a Marked Markovian Arrival Process.

    Args:
        mmap: MMAP matrices
        k: Cross moment order

    Returns:
        numpy.ndarray: Cross moment matrix
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Mmap_cross_momentKt.mmap_cross_moment(
            mmap, jpype.JInt(k)
        )
    )


def mmap_sample(MMAP, n, random=None):
    """
    Generate random samples from an MMAP.

    Produces random arrival times and classes according to
    the specified Marked Markovian Arrival Process.

    Args:
        MMAP: MMAP matrices to sample from
        n: Number of samples to generate
        random: Random number generator (optional)

    Returns:
        tuple: (times, classes) - Arrays of arrival times and class indices
    """
    if random is None:
        result = jpype.JPackage('jline').api.mam.Mmap_sampleKt.mmap_sample(
            MMAP, jpype.JInt(n)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Mmap_sampleKt.mmap_sample(
            MMAP, jpype.JInt(n), random
        )

    times = jlineMatrixToArray(result.times)
    classes = jlineMatrixToArray(result.classes)

    return times, classes


def mmap_rand(order, classes):
    """
    Generate a random MMAP with specified parameters.

    Creates a random Marked Markovian Arrival Process with
    the given order (number of states) and number of classes.

    Args:
        order: Number of states in the MMAP
        classes: Number of arrival classes

    Returns:
        Random MMAP matrices
    """
    return jpype.JPackage('jline').api.mam.Mmap_randKt.mmap_rand(
        jpype.JInt(order), jpype.JInt(classes)
    )




def map_sample(MAP, n, random=None):
    """
    Generate a random sample of inter-arrival times.
    
    Produces random inter-arrival time samples according to the specified
    Markovian Arrival Process using Monte Carlo simulation.

    Args:
        MAP: MAP as tuple (D0, D1) or MAP container {D0, D1}
        n: Number of samples to be generated
        random: Random number generator or seed (optional)

    Returns:
        numpy.ndarray: Set of n inter-arrival time samples
        
    Examples:
        - map_sample(MAP, 10) generates a sample of 10 inter-arrival times
        - map_sample(MAP, 1000, seed=42) generates 1000 samples with fixed seed
        
    Warning:
        The function can be memory consuming and quite slow for sample sizes
        greater than n=10000.
    """
    if isinstance(MAP, tuple):
        D0, D1 = MAP
        if random is None:
            return jlineMatrixToArray(
                jpype.JPackage('jline').api.mam.Map_sampleKt.map_sample(
                    jlineMatrixFromArray(D0),
                    jlineMatrixFromArray(D1),
                    jpype.JInt(n)
                )
            )
        else:
            return jlineMatrixToArray(
                jpype.JPackage('jline').api.mam.Map_sampleKt.map_sample(
                    jlineMatrixFromArray(D0),
                    jlineMatrixFromArray(D1),
                    jpype.JInt(n),
                    random
                )
            )
    else:
        if random is None:
            return jlineMatrixToArray(
                jpype.JPackage('jline').api.mam.Map_sampleKt.map_sample(
                    MAP, jpype.JInt(n)
                )
            )
        else:
            return jlineMatrixToArray(
                jpype.JPackage('jline').api.mam.Map_sampleKt.map_sample(
                    MAP, jpype.JInt(n), random
                )
            )



def mmap_count_lambda(mmap):
    """
    Compute arrival rate lambda for MMAP counting process.

    Calculates the lambda parameters for the counting process
    associated with the Marked Markovian Arrival Process.

    Args:
        mmap: MMAP matrices or MMAP object

    Returns:
        numpy.ndarray: Lambda parameters for counting process
    """
    if isinstance(mmap, list):
        java_mmap = jpype.JPackage('jline').util.matrix.MatrixCell(len(mmap))
        for matrix in mmap:
            java_mmap.add(jlineMatrixFromArray(matrix))
    else:
        java_mmap = mmap

    result = jpype.JPackage('jline').api.mam.Mmap_count_lambdaKt.mmap_count_lambda(java_mmap)
    return jlineMatrixToArray(result)


def mmap_isfeasible(mmap):
    """
    Check if an MMAP is feasible (valid).

    Determines whether the given Marked Markovian Arrival Process
    satisfies all necessary conditions for validity.

    Args:
        mmap: MMAP matrices or MMAP object

    Returns:
        bool: True if MMAP is feasible, False otherwise
    """
    if isinstance(mmap, list):
        java_mmap = jpype.JPackage('jline').util.matrix.MatrixCell(len(mmap))
        for matrix in mmap:
            java_mmap.add(jlineMatrixFromArray(matrix))
    else:
        java_mmap = mmap

    return jpype.JPackage('jline').api.mam.Mmap_isfeasibleKt.mmap_isfeasible(java_mmap)


def mmap_mark(mmap, prob):
    """
    Mark an MMAP with specified probabilities.

    Applies probability markings to a Markovian Arrival Process
    to create a Marked MAP with specified marking probabilities.

    Args:
        mmap: MMAP matrices or MMAP object
        prob: Marking probabilities matrix

    Returns:
        Marked MMAP with applied probabilities
    """
    if isinstance(mmap, list):
        java_mmap = jpype.JPackage('jline').util.matrix.MatrixCell(len(mmap))
        for matrix in mmap:
            java_mmap.add(jlineMatrixFromArray(matrix))
    else:
        java_mmap = mmap

    result = jpype.JPackage('jline').api.mam.Mmap_markKt.mmap_mark(
        java_mmap,
        jlineMatrixFromArray(prob)
    )

    new_mmap = []
    for i in range(result.length()):
        new_mmap.append(jlineMatrixToArray(result.get(i)))

    return new_mmap


def aph_bernstein(f, order):
    """
    Compute Bernstein approximation of function using Acyclic PH.

    Approximates a function using Bernstein polynomials with
    an Acyclic Phase-type distribution representation.

    Args:
        f: Function to approximate
        order: Order of Bernstein approximation

    Returns:
        Acyclic PH approximation result
    """
    class FunctionWrapper:
        def __init__(self, python_func):
            self.python_func = python_func

        def apply(self, x):
            return float(self.python_func(float(x)))

    function_wrapper = FunctionWrapper(f)

    result = jpype.JPackage('jline').api.mam.Aph_bernsteinKt.aph_bernstein(
        function_wrapper.apply,
        jpype.JInt(order)
    )

    D0 = jlineMatrixToArray(result.getFirst())
    D1 = jlineMatrixToArray(result.getSecond())

    return D0, D1


def map_jointpdf_derivative(map_matrices, iset):
    """
    Compute joint PDF derivative for MAP.

    Calculates the derivative of the joint probability density function
    for a Markovian Arrival Process with respect to specified indices.

    Args:
        map_matrices: List of MAP matrices [D0, D1]
        iset: Set of indices for derivative computation

    Returns:
        Joint PDF derivative result
    """
    java_map = jpype.JPackage('jline').util.matrix.MatrixCell(
        jpype.JArray(jpype.JPackage('jline').util.matrix.Matrix)([
            jlineMatrixFromArray(map_matrices[0]),
            jlineMatrixFromArray(map_matrices[1])
        ])
    )

    java_iset = jpype.JArray(jpype.JInt)(len(iset))
    for i, val in enumerate(iset):
        java_iset[i] = int(val)

    return jpype.JPackage('jline').api.mam.Map_jointpdf_derivativeKt.map_jointpdf_derivative(
        java_map, java_iset
    )


def map_ccdf_derivative(map_matrices, i):
    """
    Compute complementary CDF derivative for MAP.

    Calculates the derivative of the complementary cumulative
    distribution function for a Markovian Arrival Process.

    Args:
        map_matrices: List of MAP matrices [D0, D1]
        i: Index parameter for derivative computation

    Returns:
        Complementary CDF derivative result
    """
    java_map = jpype.JPackage('jline').util.matrix.MatrixCell(
        jpype.JArray(jpype.JPackage('jline').util.matrix.Matrix)([
            jlineMatrixFromArray(map_matrices[0]),
            jlineMatrixFromArray(map_matrices[1])
        ])
    )

    return jpype.JPackage('jline').api.mam.Map_ccdf_derivativeKt.map_ccdf_derivative(
        java_map, jpype.JInt(i)
    )


def qbd_R(B, L, F, iter_max=100000):
    """
    Compute the R matrix for a Quasi-Birth-Death process.

    Solves the matrix equation R = F + R*L*R + R^2*B using iterative methods.
    The R matrix is fundamental in QBD analysis.

    Args:
        B: Backward transition matrix
        L: Local transition matrix  
        F: Forward transition matrix
        iter_max: Maximum number of iterations (default: 100000)

    Returns:
        numpy.ndarray: The R matrix
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Qbd_RKt.qbd_R(
            jlineMatrixFromArray(B),
            jlineMatrixFromArray(L),
            jlineMatrixFromArray(F),
            jpype.JInt(iter_max)
        )
    )


def qbd_R_logred(B, L, F, iter_max=100000):
    """
    Compute the R matrix using logarithmic reduction algorithm.

    Alternative method for computing the R matrix that can be more
    numerically stable for certain QBD problems.

    Args:
        B: Backward transition matrix
        L: Local transition matrix
        F: Forward transition matrix
        iter_max: Maximum number of iterations (default: 100000)

    Returns:
        numpy.ndarray: The R matrix
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Qbd_R_logredKt.qbd_R_logred(
            jlineMatrixFromArray(B),
            jlineMatrixFromArray(L),
            jlineMatrixFromArray(F),
            jpype.JInt(iter_max)
        )
    )


def qbd_rg(map_a, map_s, util=None):
    """
    Analyze a MAP/MAP/1 queue using QBD methods.

    Analyzes a single-server queue with MAP arrivals and MAP service times
    using Quasi-Birth-Death process techniques.

    Args:
        map_a: Arrival MAP (list/array with D0, D1 matrices)
        map_s: Service MAP (list/array with D0, D1 matrices)
        util: Utilization constraint (optional)

    Returns:
        QBD analysis results including performance measures
    """
    java_map_a = jpype.JPackage('jline').util.matrix.MatrixCell(
        jpype.JArray(jpype.JPackage('jline').util.matrix.Matrix)([
            jlineMatrixFromArray(map_a[0]),
            jlineMatrixFromArray(map_a[1])
        ])
    )

    java_map_s = jpype.JPackage('jline').util.matrix.MatrixCell(
        jpype.JArray(jpype.JPackage('jline').util.matrix.Matrix)([
            jlineMatrixFromArray(map_s[0]),
            jlineMatrixFromArray(map_s[1])
        ])
    )

    if util is not None:
        result = jpype.JPackage('jline').api.mam.Qbd_rgKt.qbd_rg(
            java_map_a, java_map_s, jpype.JDouble(util)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Qbd_rgKt.qbd_rg(
            java_map_a, java_map_s, None
        )

    return {
        'R': jlineMatrixToArray(result.getR()),
        'G': jlineMatrixToArray(result.getG()),
        'B': jlineMatrixToArray(result.getB()),
        'L': jlineMatrixToArray(result.getL()),
        'F': jlineMatrixToArray(result.getF()),
        'U': jlineMatrixToArray(result.getU())
    }


def map_pdf(MAP, points):
    """
    Compute probability density function of interarrival times.
    
    Evaluates the probability density function (PDF) of the interarrival time
    distribution for a Markovian Arrival Process at the specified time points.
    
    Args:
        MAP: MAP as tuple (D0, D1) representing the Markovian Arrival Process
        points: Array of time points to evaluate the PDF at
        
    Returns:
        numpy.ndarray: PDF values at the specified points, returned in the
                      same order as the points vector
                      
    Examples:
        - map_pdf(MAP, [0.5, 1.0, 2.0]) returns PDF values at t=0.5, 1.0, 2.0
    """
    if isinstance(MAP, (list, tuple)) and len(MAP) == 2:
        D0, D1 = MAP
        map_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1)
        )
    else:
        raise ValueError("MAP must be a tuple/list of (D0, D1) matrices")

    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Map_pdfKt.map_pdf(
            map_cell, jlineMatrixFromArray(points)
        )
    )


def map_prob(MAP, t):
    """
    Compute probability for MAP at time t.
    
    Calculates the probability measure associated with
    the Markovian Arrival Process at the specified time.
    
    Args:
        MAP: MAP as tuple (D0, D1) matrices
        t: Time point for probability computation
        
    Returns:
        float: Probability value at time t
    """
    if isinstance(MAP, (list, tuple)) and len(MAP) == 2:
        D0, D1 = MAP
        map_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1)
        )
    else:
        raise ValueError("MAP must be a tuple/list of (D0, D1) matrices")

    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mam.Map_probKt.map_prob(
            map_cell, jpype.JDouble(t)
        )
    )


def map_joint(MAP1, MAP2):
    """
    Compute joint distribution of two MAPs.
    
    Creates the joint distribution representation of two
    independent Markovian Arrival Processes.
    
    Args:
        MAP1: First MAP as tuple (D0, D1)
        MAP2: Second MAP as tuple (D0, D1)
        
    Returns:
        Joint MAP representation
    """
    if isinstance(MAP1, (list, tuple)) and len(MAP1) == 2:
        D0_1, D1_1 = MAP1
        map_cell_1 = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0_1), jlineMatrixFromArray(D1_1)
        )
    else:
        raise ValueError("MAP1 must be a tuple/list of (D0, D1) matrices")

    if isinstance(MAP2, (list, tuple)) and len(MAP2) == 2:
        D0_2, D1_2 = MAP2
        map_cell_2 = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0_2), jlineMatrixFromArray(D1_2)
        )
    else:
        raise ValueError("MAP2 must be a tuple/list of (D0, D1) matrices")

    result = jpype.JPackage('jline').api.mam.Map_jointKt.map_joint(
        map_cell_1, map_cell_2
    )

    return (jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1)))


def map_mixture(alpha, MAPs):
    """
    Create mixture of multiple MAPs.
    
    Constructs a mixture distribution from multiple Markovian
    Arrival Processes with given mixing probabilities.
    
    Args:
        alpha: Mixing probabilities (weights) for each MAP
        MAPs: List of MAP matrices to mix
        
    Returns:
        Mixture MAP representation
    """
    java_alpha = jpype.JArray(jpype.JDouble)(list(alpha))

    java_maps = jpype.JArray(jpype.JClass("jline.lang.MatrixCell"))(len(MAPs))
    for i, MAP in enumerate(MAPs):
        if isinstance(MAP, (list, tuple)) and len(MAP) == 2:
            D0, D1 = MAP
            java_maps[i] = jpype.JClass("jline.lang.MatrixCell")(
                jlineMatrixFromArray(D0), jlineMatrixFromArray(D1)
            )
        else:
            raise ValueError(f"MAP {i} must be a tuple/list of (D0, D1) matrices")

    result = jpype.JPackage('jline').api.mam.Map_mixtureKt.map_mixture(
        java_alpha, java_maps
    )

    return (jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1)))


def map_max(MAP1, MAP2):
    """
    Compute maximum of two independent MAPs.

    Args:
        MAP1: First MAP as tuple (D0, D1)
        MAP2: Second MAP as tuple (D0, D1)

    Returns:
        MAP representing max(X,Y) where X~MAP1, Y~MAP2
    """
    if isinstance(MAP1, (list, tuple)) and len(MAP1) == 2:
        D0_1, D1_1 = MAP1
        map_cell_1 = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0_1), jlineMatrixFromArray(D1_1)
        )
    else:
        raise ValueError("MAP1 must be a tuple/list of (D0, D1) matrices")

    if isinstance(MAP2, (list, tuple)) and len(MAP2) == 2:
        D0_2, D1_2 = MAP2
        map_cell_2 = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0_2), jlineMatrixFromArray(D1_2)
        )
    else:
        raise ValueError("MAP2 must be a tuple/list of (D0, D1) matrices")

    result = jpype.JPackage('jline').api.mam.Map_maxKt.map_max(
        map_cell_1, map_cell_2
    )

    return (jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1)))


def map_renewal(MAPIN):
    """
    Remove all correlations from a MAP to create a renewal process.
    
    Transforms a Markovian Arrival Process into a renewal MAP process
    with the same cumulative distribution function (CDF) as the input MAP,
    but with no correlations between inter-arrival times.
    
    Args:
        MAPIN: Input MAP as tuple (D0, D1) representing the original MAP
        
    Returns:
        Renewal MAP with same marginal distribution but independent inter-arrival times
        
    Note:
        The resulting process maintains the same inter-arrival time distribution
        but removes all temporal dependencies, making it a renewal process.
    """
    if isinstance(MAPIN, (list, tuple)) and len(MAPIN) == 2:
        D0, D1 = MAPIN
        map_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1)
        )
    else:
        raise ValueError("MAPIN must be a tuple/list of (D0, D1) matrices")

    result = jpype.JPackage('jline').api.mam.Map_renewalKt.map_renewal(map_cell)

    return (jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1)))


def map_stochcomp(MAP):
    """
    Compute stochastic complement of MAP by eliminating states.

    Args:
        MAP: MAP as tuple (D0, D1)

    Returns:
        Reduced MAP with eliminated transient states
    """
    if isinstance(MAP, (list, tuple)) and len(MAP) == 2:
        D0, D1 = MAP
        map_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1)
        )
    else:
        raise ValueError("MAP must be a tuple/list of (D0, D1) matrices")

    result = jpype.JPackage('jline').api.mam.Map_stochcompKt.map_stochcomp(map_cell)

    return (jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1)))


def qbd_mapmap1(MAP_A, MAP_S, mu=None):
    """
    Analyze MAP/MAP/1 queueing system using QBD methods.

    Args:
        MAP_A: Arrival MAP as tuple (D0, D1)
        MAP_S: Service MAP as tuple (D0, D1)
        mu: Service rate parameter (optional)

    Returns:
        dict: QBD analysis results including performance measures
    """
    if isinstance(MAP_A, (list, tuple)) and len(MAP_A) == 2:
        D0_A, D1_A = MAP_A
        map_a_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0_A), jlineMatrixFromArray(D1_A)
        )
    else:
        raise ValueError("MAP_A must be a tuple/list of (D0, D1) matrices")

    if isinstance(MAP_S, (list, tuple)) and len(MAP_S) == 2:
        D0_S, D1_S = MAP_S
        map_s_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0_S), jlineMatrixFromArray(D1_S)
        )
    else:
        raise ValueError("MAP_S must be a tuple/list of (D0, D1) matrices")

    if mu is not None:
        result = jpype.JPackage('jline').api.mam.Qbd_mapmap1Kt.qbd_mapmap1(
            map_a_cell, map_s_cell, jpype.JDouble(mu)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Qbd_mapmap1Kt.qbd_mapmap1(
            map_a_cell, map_s_cell
        )

    return {
        'pi': jlineMatrixToArray(result.getPi()) if hasattr(result, 'getPi') else None,
        'R': jlineMatrixToArray(result.getR()) if hasattr(result, 'getR') else None,
        'utilization': float(result.getUtilization()) if hasattr(result, 'getUtilization') else None,
        'mean_queue_length': float(result.getMeanQueueLength()) if hasattr(result, 'getMeanQueueLength') else None,
        'mean_waiting_time': float(result.getMeanWaitingTime()) if hasattr(result, 'getMeanWaitingTime') else None
    }


def qbd_raprap1(RAP_A, RAP_S, util=None):
    """
    Analyze RAP/RAP/1 queueing system using QBD methods.

    Args:
        RAP_A: Arrival RAP (Rational Arrival Process)
        RAP_S: Service RAP
        util: Utilization parameter (optional)

    Returns:
        dict: QBD analysis results for RAP/RAP/1 queue
    """
    if isinstance(RAP_A, (list, tuple)) and len(RAP_A) == 2:
        D0_A, D1_A = RAP_A
        rap_a_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0_A), jlineMatrixFromArray(D1_A)
        )
    else:
        raise ValueError("RAP_A must be a tuple/list of (D0, D1) matrices")

    if isinstance(RAP_S, (list, tuple)) and len(RAP_S) == 2:
        D0_S, D1_S = RAP_S
        rap_s_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0_S), jlineMatrixFromArray(D1_S)
        )
    else:
        raise ValueError("RAP_S must be a tuple/list of (D0, D1) matrices")

    if util is not None:
        result = jpype.JPackage('jline').api.mam.Qbd_raprap1Kt.qbd_raprap1(
            rap_a_cell, rap_s_cell, jpype.JDouble(util)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Qbd_raprap1Kt.qbd_raprap1(
            rap_a_cell, rap_s_cell
        )

    return {
        'R': jlineMatrixToArray(result.getR()) if hasattr(result, 'getR') else None,
        'G': jlineMatrixToArray(result.getG()) if hasattr(result, 'getG') else None,
        'pi': jlineMatrixToArray(result.getPi()) if hasattr(result, 'getPi') else None,
        'performance_metrics': result.getPerformanceMetrics() if hasattr(result, 'getPerformanceMetrics') else None
    }


def qbd_bmapbmap1(BMAP_A, BMAP_S):
    """
    Analyze BMAP/BMAP/1 queueing system using QBD methods.

    Args:
        BMAP_A: Arrival Batch MAP
        BMAP_S: Service Batch MAP

    Returns:
        dict: QBD analysis results for batch arrival/service system
    """
    if isinstance(BMAP_A, list) and len(BMAP_A) >= 2:
        java_bmap_a = jpype.JClass("jline.lang.MatrixCell")(len(BMAP_A))
        for i, matrix in enumerate(BMAP_A):
            java_bmap_a.set(i, jlineMatrixFromArray(matrix))
    else:
        raise ValueError("BMAP_A must be a list of at least 2 matrices [D0, D1, ...]")

    if isinstance(BMAP_S, list) and len(BMAP_S) >= 2:
        java_bmap_s = jpype.JClass("jline.lang.MatrixCell")(len(BMAP_S))
        for i, matrix in enumerate(BMAP_S):
            java_bmap_s.set(i, jlineMatrixFromArray(matrix))
    else:
        raise ValueError("BMAP_S must be a list of at least 2 matrices [D0, D1, ...]")

    result = jpype.JPackage('jline').api.mam.Qbd_bmapbmap1Kt.qbd_bmapbmap1(
        java_bmap_a, java_bmap_s
    )

    return {
        'R': jlineMatrixToArray(result.getR()) if hasattr(result, 'getR') else None,
        'G': jlineMatrixToArray(result.getG()) if hasattr(result, 'getG') else None,
        'performance_metrics': result.getMetrics() if hasattr(result, 'getMetrics') else None
    }


def qbd_setupdelayoff(MAP_A, MAP_S, setup_time, delay_time, off_time):
    """
    Analyze queue with setup, delay, and off periods using QBD methods.
    
    Models a queueing system with server setup times, processing delays,
    and off periods using Quasi-Birth-Death process analysis.
    
    Args:
        MAP_A: Arrival MAP as tuple (D0, D1).
        MAP_S: Service MAP as tuple (D0, D1).
        setup_time: Server setup time.
        delay_time: Processing delay time.
        off_time: Server off time.
        
    Returns:
        dict: Performance measures for the system with timing constraints.
    """
    if isinstance(MAP_A, (list, tuple)) and len(MAP_A) == 2:
        D0_A, D1_A = MAP_A
        map_a_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0_A), jlineMatrixFromArray(D1_A)
        )
    else:
        raise ValueError("MAP_A must be a tuple/list of (D0, D1) matrices")

    if isinstance(MAP_S, (list, tuple)) and len(MAP_S) == 2:
        D0_S, D1_S = MAP_S
        map_s_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0_S), jlineMatrixFromArray(D1_S)
        )
    else:
        raise ValueError("MAP_S must be a tuple/list of (D0, D1) matrices")

    result = jpype.JPackage('jline').api.mam.Qbd_setupdelayoffKt.qbd_setupdelayoff(
        map_a_cell, map_s_cell,
        jpype.JDouble(setup_time), jpype.JDouble(delay_time), jpype.JDouble(off_time)
    )

    return {
        'steady_state': jlineMatrixToArray(result.getSteadyState()) if hasattr(result, 'getSteadyState') else None,
        'performance_metrics': result.getMetrics() if hasattr(result, 'getMetrics') else None
    }


def aph_simplify(a1, T1, a2, T2, p1, p2, pattern):
    """
    Simplify APH representation using specified pattern.

    Args:
        a1: Initial vector for first APH
        T1: Generator matrix for first APH
        a2: Initial vector for second APH
        T2: Generator matrix for second APH
        p1: Probability parameter for first APH
        p2: Probability parameter for second APH
        pattern: Simplification pattern identifier

    Returns:
        tuple: (alpha, T) simplified APH representation
    """
    result = jpype.JPackage('jline').api.mam.Aph_simplifyKt.aph_simplify(
        jlineMatrixFromArray(a1), jlineMatrixFromArray(T1),
        jlineMatrixFromArray(a2), jlineMatrixFromArray(T2),
        jpype.JDouble(p1), jpype.JDouble(p2), jpype.JInt(pattern)
    )

    return (jlineMatrixToArray(result.getFirst()), jlineMatrixToArray(result.getSecond()))


def randp(P, rows, cols=None):
    """
    Pick random values with relative probability.

    Returns integers in the range from 1 to len(P) with a relative probability,
    so that the value X is present approximately (P[X-1]/sum(P)) times.

    Args:
        P: Probability vector (all values should be >= 0)
        rows: Number of rows for output matrix
        cols: Number of columns for output matrix (defaults to rows if not specified)

    Returns:
        Matrix of random integers with specified probabilities
    """
    if cols is None:
        cols = rows

    # Convert P to Java double array
    if hasattr(P, '__iter__'):
        P_array = jpype.JArray(jpype.JDouble)(len(P))
        for i, val in enumerate(P):
            P_array[i] = float(val)
    else:
        P_array = jpype.JArray(jpype.JDouble)([float(P)])

    result = jpype.JPackage('jline').api.mam.RandpKt.randp(
        P_array, jpype.JInt(rows), jpype.JInt(cols)
    )
    return jlineMatrixToArray(result)


def aph_rand(n):
    """
    Generate random APH distribution with n phases.

    Args:
        n: Number of phases

    Returns:
        tuple: (alpha, T) random APH representation
    """
    result = jpype.JPackage('jline').api.mam.Aph_randKt.aph_rand(jpype.JInt(n))

    return (jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1)))


def amap2_fit_gamma(mean1, var1, mean2, var2, p):
    """
    Fit 2-phase AMAP using gamma matching method.

    Args:
        mean1: Mean for first phase
        var1: Variance for first phase
        mean2: Mean for second phase
        var2: Variance for second phase
        p: Phase probability

    Returns:
        tuple: (D0, D1) fitted AMAP matrices
    """
    result = jpype.JPackage('jline').api.mam.Amap2_fit_gammaKt.amap2_fit_gamma(
        jpype.JDouble(mean1), jpype.JDouble(var1),
        jpype.JDouble(mean2), jpype.JDouble(var2),
        jpype.JDouble(p)
    )

    return (jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1)))


def mamap2m_fit_fb_multiclass(data, classes, options=None):
    """
    Fit multiclass MAMAP using feedback method.

    Args:
        data: Input data for fitting
        classes: Number of classes
        options: Fitting options (optional)

    Returns:
        dict: Fitted AMAP and quality metrics
    """
    java_options = None
    if options is not None:
        pass

    if java_options is not None:
        result = jpype.JPackage('jline').api.mam.Mamap2m_fit_fb_multiclassKt.mamap2m_fit_fb_multiclass(
            jlineMatrixFromArray(data), jpype.JInt(classes), java_options
        )
    else:
        result = jpype.JPackage('jline').api.mam.Mamap2m_fit_fb_multiclassKt.mamap2m_fit_fb_multiclass(
            jlineMatrixFromArray(data), jpype.JInt(classes)
        )

    return {
        'fitted_amap': result.getFittedAmap() if hasattr(result, 'getFittedAmap') else None,
        'quality_metrics': result.getQualityMetrics() if hasattr(result, 'getQualityMetrics') else None
    }


def mmpp_rand(states, lambda_range=(0.1, 5.0)):
    result = jpype.JPackage('jline').api.mam.Mmpp_randKt.mmpp_rand(
        jpype.JInt(states), jpype.JDouble(lambda_range[0]), jpype.JDouble(lambda_range[1])
    )

    return (jlineMatrixToArray(result.getQ()), jlineMatrixToArray(result.getLambda()))


def map_count_moment(MAP, k, lag=0):
    """
    Compute count moments of MAP arrivals.

    Args:
        MAP: MAP as tuple (D0, D1)
        k: Moment order
        lag: Time lag (default: 0)

    Returns:
        float: Count moment value
    """
    if isinstance(MAP, (list, tuple)) and len(MAP) == 2:
        D0, D1 = MAP
        map_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1)
        )
    else:
        raise ValueError("MAP must be a tuple/list of (D0, D1) matrices")

    return float(jpype.JPackage('jline').api.mam.Map_count_momentKt.map_count_moment(
        map_cell, jpype.JInt(k), jpype.JInt(lag)
    ))


def map_kurt(MAP):
    """
    Calculate kurtosis of inter-arrival times for a MAP.
    
    Computes the fourth moment (kurtosis) of the inter-arrival time
    distribution for a Markovian Arrival Process.
    
    Args:
        MAP: MAP as tuple (D0, D1) of generator matrices.
        
    Returns:
        float: Kurtosis of inter-arrival times.
    """
    if isinstance(MAP, (list, tuple)) and len(MAP) == 2:
        D0, D1 = MAP
        map_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1)
        )
    else:
        raise ValueError("MAP must be a tuple/list of (D0, D1) matrices")

    return float(jpype.JPackage('jline').api.mam.Map_kurtKt.map_kurt(map_cell))


def mmap_sigma2_cell(MMAP):
    """
    Compute two-step class transition probabilities for MMAP (cell version).

    Args:
        MMAP: Multi-class MAP matrices as cell structure

    Returns:
        3D matrix of class transition probabilities
    """
    if isinstance(MMAP, list):
        mmap_cell = jpype.JClass("jline.lang.MatrixCell")()
        for i, matrix in enumerate(MMAP):
            mmap_cell.set(i, jlineMatrixFromArray(matrix))
    else:
        mmap_cell = MMAP

    return float(jpype.JPackage('jline').api.mam.Mmap_sigma2_cellKt.mmap_sigma2_cell(mmap_cell))


def amap2_adjust_gamma(mean1, var1, mean2, var2, p, target_mean, target_var):
    """
    Adjust 2-phase AMAP gamma parameters to target moments.

    Args:
        mean1: Mean for first phase
        var1: Variance for first phase
        mean2: Mean for second phase
        var2: Variance for second phase
        p: Phase probability
        target_mean: Target mean value
        target_var: Target variance value

    Returns:
        Adjusted AMAP parameters
    """
    result = jpype.JPackage('jline').api.mam.Amap2_adjust_gammaKt.amap2_adjust_gamma(
        jpype.JDouble(mean1), jpype.JDouble(var1),
        jpype.JDouble(mean2), jpype.JDouble(var2),
        jpype.JDouble(p), jpype.JDouble(target_mean), jpype.JDouble(target_var)
    )

    return {
        'mean1': result.mean1 if hasattr(result, 'mean1') else None,
        'var1': result.var1 if hasattr(result, 'var1') else None,
        'mean2': result.mean2 if hasattr(result, 'mean2') else None,
        'var2': result.var2 if hasattr(result, 'var2') else None,
        'p': result.p if hasattr(result, 'p') else None
    }


def amap2_fitall_gamma(data, options=None):
    """
    Fit all gamma parameters for 2-phase AMAP from data.

    Args:
        data: Input data for fitting
        options: Fitting options (optional)

    Returns:
        dict: Fitted AMAP parameters and quality metrics
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.mam.Amap2_fitall_gammaKt.amap2_fitall_gamma(
            jlineMatrixFromArray(data), java_options
        )
    else:
        result = jpype.JPackage('jline').api.mam.Amap2_fitall_gammaKt.amap2_fitall_gamma(
            jlineMatrixFromArray(data)
        )

    return {
        'parameters': {
            'mean1': result.mean1 if hasattr(result, 'mean1') else None,
            'var1': result.var1 if hasattr(result, 'var1') else None,
            'mean2': result.mean2 if hasattr(result, 'var2') else None,
            'var2': result.var2 if hasattr(result, 'var2') else None,
            'p': result.p if hasattr(result, 'p') else None
        },
        'quality': result.quality if hasattr(result, 'quality') else None,
        'converged': result.converged if hasattr(result, 'converged') else None
    }


def mmpp2_fit_mu00(data, options=None):
    """
    Fit 2-phase MMPP parameter mu00 from data.

    Args:
        data: Input data for parameter fitting
        options: Fitting options (optional)

    Returns:
        Fitted mu00 parameter
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.mam.Mmpp2_fit_mu00Kt.mmpp2_fit_mu00(
            jlineMatrixFromArray(data), java_options
        )
    else:
        result = jpype.JPackage('jline').api.mam.Mmpp2_fit_mu00Kt.mmpp2_fit_mu00(
            jlineMatrixFromArray(data)
        )

    return float(result)


def mmpp2_fit_mu11(data, options=None):
    """
    Fit 2-phase MMPP parameter mu11 from data.

    Args:
        data: Input data for parameter fitting
        options: Fitting options (optional)

    Returns:
        Fitted mu11 parameter
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.mam.Mmpp2_fit_mu11Kt.mmpp2_fit_mu11(
            jlineMatrixFromArray(data), java_options
        )
    else:
        result = jpype.JPackage('jline').api.mam.Mmpp2_fit_mu11Kt.mmpp2_fit_mu11(
            jlineMatrixFromArray(data)
        )

    return float(result)


def mmpp2_fit_q01(data, options=None):
    """
    Fit 2-phase MMPP parameter q01 from data.

    Args:
        data: Input data for parameter fitting
        options: Fitting options (optional)

    Returns:
        Fitted q01 parameter
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.mam.Mmpp2_fit_q01Kt.mmpp2_fit_q01(
            jlineMatrixFromArray(data), java_options
        )
    else:
        result = jpype.JPackage('jline').api.mam.Mmpp2_fit_q01Kt.mmpp2_fit_q01(
            jlineMatrixFromArray(data)
        )

    return float(result)


def mmpp2_fit_q10(data, options=None):
    """
    Fit 2-phase MMPP parameter q10 from data.

    Args:
        data: Input data for parameter fitting
        options: Fitting options (optional)

    Returns:
        Fitted q10 parameter
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.mam.Mmpp2_fit_q10Kt.mmpp2_fit_q10(
            jlineMatrixFromArray(data), java_options
        )
    else:
        result = jpype.JPackage('jline').api.mam.Mmpp2_fit_q10Kt.mmpp2_fit_q10(
            jlineMatrixFromArray(data)
        )

    return float(result)


def assess_compression_quality(original_MAP, compressed_MAP, metrics=['mean', 'var', 'acf']):
    """
    Assess quality of MAP compression by comparing metrics.
    
    Evaluates how well a compressed MAP approximates the original MAP
    by comparing specified statistical metrics.
    
    Args:
        original_MAP: Original MAP as tuple (D0, D1).
        compressed_MAP: Compressed MAP as tuple (D0, D1).
        metrics: List of metrics to compare (default: ['mean', 'var', 'acf']).
        
    Returns:
        dict: Quality assessment results for each metric.
    """
    if isinstance(original_MAP, (list, tuple)) and len(original_MAP) == 2:
        D0_orig, D1_orig = original_MAP
        orig_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0_orig), jlineMatrixFromArray(D1_orig)
        )
    else:
        raise ValueError("original_MAP must be a tuple/list of (D0, D1) matrices")

    if isinstance(compressed_MAP, (list, tuple)) and len(compressed_MAP) == 2:
        D0_comp, D1_comp = compressed_MAP
        comp_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0_comp), jlineMatrixFromArray(D1_comp)
        )
    else:
        raise ValueError("compressed_MAP must be a tuple/list of (D0, D1) matrices")

    java_metrics = jpype.JPackage('java.util').ArrayList()
    for metric in metrics:
        java_metrics.add(jpype.JObject(metric, jpype.JClass("java.lang.String")))

    result = jpype.JPackage('jline').api.mam.Assess_compression_qualityKt.assess_compression_quality(
        orig_cell, comp_cell, java_metrics
    )

    return {
        'error_metrics': jlineMatrixToArray(result.errorMetrics) if hasattr(result, 'errorMetrics') else None,
        'relative_errors': jlineMatrixToArray(result.relativeErrors) if hasattr(result, 'relativeErrors') else None,
        'quality_score': result.qualityScore if hasattr(result, 'qualityScore') else None
    }


def compress_adaptive(MAP, target_order, options=None):
    """
    Compress MAP using adaptive compression algorithm.

    Args:
        MAP: MAP as tuple (D0, D1)
        target_order: Target order for compression
        options: Compression options (optional)

    Returns:
        Compressed MAP with reduced order
    """
    if isinstance(MAP, (list, tuple)) and len(MAP) == 2:
        D0, D1 = MAP
        map_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1)
        )
    else:
        raise ValueError("MAP must be a tuple/list of (D0, D1) matrices")

    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.mam.Compress_adaptiveKt.compress_adaptive(
            map_cell, jpype.JInt(target_order), java_options
        )
    else:
        result = jpype.JPackage('jline').api.mam.Compress_adaptiveKt.compress_adaptive(
            map_cell, jpype.JInt(target_order)
        )

    return (jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1)))


def compress_autocorrelation(MAP, target_order, options=None):
    """
    Compress MAP preserving autocorrelation properties.

    Args:
        MAP: MAP as tuple (D0, D1)
        target_order: Target order for compression
        options: Compression options (optional)

    Returns:
        Compressed MAP preserving autocorrelation structure
    """
    if isinstance(MAP, (list, tuple)) and len(MAP) == 2:
        D0, D1 = MAP
        map_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1)
        )
    else:
        raise ValueError("MAP must be a tuple/list of (D0, D1) matrices")

    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.mam.Compress_autocorrelationKt.compress_autocorrelation(
            map_cell, jpype.JInt(target_order), java_options
        )
    else:
        result = jpype.JPackage('jline').api.mam.Compress_autocorrelationKt.compress_autocorrelation(
            map_cell, jpype.JInt(target_order)
        )

    return (jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1)))


def compress_spectral(MAP, target_order, options=None):
    """
    Compress MAP using spectral decomposition methods.

    Args:
        MAP: MAP as tuple (D0, D1)
        target_order: Target order for compression
        options: Compression options (optional)

    Returns:
        Compressed MAP using spectral techniques
    """
    if isinstance(MAP, (list, tuple)) and len(MAP) == 2:
        D0, D1 = MAP
        map_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1)
        )
    else:
        raise ValueError("MAP must be a tuple/list of (D0, D1) matrices")

    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.mam.Compress_spectralKt.compress_spectral(
            map_cell, jpype.JInt(target_order), java_options
        )
    else:
        result = jpype.JPackage('jline').api.mam.Compress_spectralKt.compress_spectral(
            map_cell, jpype.JInt(target_order)
        )

    return (jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1)))


def compress_with_quality_control(MAP, target_order, quality_threshold=0.95, options=None):
    """
    Compress MAP with quality control constraints.

    Args:
        MAP: MAP as tuple (D0, D1)
        target_order: Target order for compression
        quality_threshold: Quality threshold (default: 0.95)
        options: Compression options (optional)

    Returns:
        dict: Compressed MAP and quality metrics
    """
    if isinstance(MAP, (list, tuple)) and len(MAP) == 2:
        D0, D1 = MAP
        map_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1)
        )
    else:
        raise ValueError("MAP must be a tuple/list of (D0, D1) matrices")

    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.mam.Compress_with_quality_controlKt.compress_with_quality_control(
            map_cell, jpype.JInt(target_order), jpype.JDouble(quality_threshold), java_options
        )
    else:
        result = jpype.JPackage('jline').api.mam.Compress_with_quality_controlKt.compress_with_quality_control(
            map_cell, jpype.JInt(target_order), jpype.JDouble(quality_threshold)
        )

    return {
        'compressed_MAP': (jlineMatrixToArray(result.compressedMAP.get(0)), jlineMatrixToArray(result.compressedMAP.get(1))),
        'quality_achieved': result.qualityAchieved if hasattr(result, 'qualityAchieved') else None,
        'meets_threshold': result.meetsThreshold if hasattr(result, 'meetsThreshold') else None
    }


def qbd_G(A0, A1, A2):
    """
    Compute G matrix for Quasi-Birth-Death (QBD) process.

    Calculates the G matrix which represents the probability
    generating function for the first passage time from level k to k+1.

    Args:
        A0: Transition rate matrix for level down
        A1: Transition rate matrix for same level
        A2: Transition rate matrix for level up

    Returns:
        numpy.ndarray: G matrix
    """
    A0_matrix = jlineMatrixFromArray(A0)
    A1_matrix = jlineMatrixFromArray(A1)
    A2_matrix = jlineMatrixFromArray(A2)

    result = jpype.JPackage('jline').api.mam.Qbd_GKt.qbd_G(
        A0_matrix, A1_matrix, A2_matrix
    )

    return jlineMatrixToArray(result)


def ph_fit(data, max_phases=10):
    """
    Fit phase-type distribution to empirical data.
    
    Estimates the parameters of a phase-type distribution that
    best fits the given data using maximum likelihood estimation.
    
    Args:
        data: Empirical data to fit
        max_phases: Maximum number of phases to use
        
    Returns:
        tuple: (alpha, T) - Initial probability vector and transition rate matrix
    """
    data_array = jpype.JArray(jpype.JDouble)(data)

    result = jpype.JPackage('jline').api.mam.Ph_fitKt.ph_fit(
        data_array, jpype.JInt(max_phases)
    )

    alpha = jlineMatrixToArray(result.getFirst())
    T = jlineMatrixToArray(result.getSecond())

    return alpha, T


def ph_mean(alpha, T):
    """
    Compute mean of phase-type distribution.
    
    Calculates the expected value (first moment) of a phase-type
    distribution with given initial vector and transition matrix.
    
    Args:
        alpha: Initial probability vector
        T: Transition rate matrix
        
    Returns:
        float: Mean of the phase-type distribution
    """
    alpha_matrix = jlineMatrixFromArray(alpha)
    T_matrix = jlineMatrixFromArray(T)

    result = jpype.JPackage('jline').api.mam.Ph_meanKt.ph_mean(
        alpha_matrix, T_matrix
    )

    return float(result)


def ph_var(alpha, T):
    """
    Compute variance of phase-type distribution.
    
    Calculates the variance (second central moment) of a phase-type
    distribution with given parameters.
    
    Args:
        alpha: Initial probability vector
        T: Transition rate matrix
        
    Returns:
        float: Variance of the phase-type distribution
    """
    alpha_matrix = jlineMatrixFromArray(alpha)
    T_matrix = jlineMatrixFromArray(T)

    result = jpype.JPackage('jline').api.mam.Ph_varKt.ph_var(
        alpha_matrix, T_matrix
    )

    return float(result)


def ph_pdf(alpha, T, points):
    """
    Compute PDF of phase-type distribution at given points.
    
    Evaluates the probability density function of a phase-type
    distribution at the specified points.
    
    Args:
        alpha: Initial probability vector
        T: Transition rate matrix
        points: Points at which to evaluate the PDF
        
    Returns:
        numpy.ndarray: PDF values at the specified points
    """
    alpha_matrix = jlineMatrixFromArray(alpha)
    T_matrix = jlineMatrixFromArray(T)
    points_matrix = jlineMatrixFromArray(points)

    result = jpype.JPackage('jline').api.mam.Ph_pdfKt.ph_pdf(
        alpha_matrix, T_matrix, points_matrix
    )

    return jlineMatrixToArray(result)


def ph_cdf(alpha, T, points):
    """
    Compute CDF of phase-type distribution at given points.
    
    Evaluates the cumulative distribution function of a phase-type
    distribution at the specified points.
    
    Args:
        alpha: Initial probability vector
        T: Transition rate matrix
        points: Points at which to evaluate the CDF
        
    Returns:
        numpy.ndarray: CDF values at the specified points
    """
    alpha_matrix = jlineMatrixFromArray(alpha)
    T_matrix = jlineMatrixFromArray(T)
    points_matrix = jlineMatrixFromArray(points)

    result = jpype.JPackage('jline').api.mam.Ph_cdfKt.ph_cdf(
        alpha_matrix, T_matrix, points_matrix
    )

    return jlineMatrixToArray(result)


def qbd_psif(A0, A1, A2, B, options=None):
    """
    Compute finite QBD stationary probability vector.

    Calculates the stationary probability vector for a finite
    Quasi-Birth-Death process with boundary conditions.

    Args:
        A0: Transition rate matrix for level down
        A1: Transition rate matrix for same level
        A2: Transition rate matrix for level up
        B: Boundary condition matrix
        options: Optional solver options

    Returns:
        numpy.ndarray: Stationary probability vector
    """
    A0_matrix = jlineMatrixFromArray(A0)
    A1_matrix = jlineMatrixFromArray(A1)
    A2_matrix = jlineMatrixFromArray(A2)
    B_matrix = jlineMatrixFromArray(B)

    if options is not None:
        result = jpype.JPackage('jline').api.mam.Qbd_psifKt.qbd_psif(
            A0_matrix, A1_matrix, A2_matrix, B_matrix,
            jpype.JObject(options)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Qbd_psifKt.qbd_psif(
            A0_matrix, A1_matrix, A2_matrix, B_matrix
        )

    return jlineMatrixToArray(result)


def qbd_psi(A0, A1, A2, options=None):
    """
    Compute Psi matrix for QBD process using iterative methods.

    Args:
        A0: QBD backward transition matrix
        A1: QBD local transition matrix
        A2: QBD forward transition matrix
        options: Computational options (optional)

    Returns:
        numpy.ndarray: Psi matrix for QBD fundamental matrices
    """
    A0_matrix = jlineMatrixFromArray(A0)
    A1_matrix = jlineMatrixFromArray(A1)
    A2_matrix = jlineMatrixFromArray(A2)

    if options is not None:
        result = jpype.JPackage('jline').api.mam.Qbd_psiKt.qbd_psi(
            A0_matrix, A1_matrix, A2_matrix, jpype.JObject(options)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Qbd_psiKt.qbd_psi(
            A0_matrix, A1_matrix, A2_matrix
        )

    return jlineMatrixToArray(result)


def aph2_check_feasibility(M1, M2, M3):
    """
    Check feasibility of 2-phase APH distribution parameters.

    Args:
        M1: First moment
        M2: Second moment
        M3: Third moment

    Returns:
        bool: True if parameters are feasible for 2-phase APH
    """
    result = jpype.JPackage('jline').api.mam.Aph2_check_feasibilityKt.aph2_check_feasibility(
        jpype.JDouble(M1), jpype.JDouble(M2), jpype.JDouble(M3)
    )

    return bool(result)


def aph2_canonical(a1, t11, a2, t22):
    """
    Convert 2-phase APH to canonical form.

    Args:
        a1: Initial probability for first phase
        t11: Transition rate for first phase
        a2: Initial probability for second phase
        t22: Transition rate for second phase

    Returns:
        tuple: Canonical form APH parameters
    """
    result = jpype.JPackage('jline').api.mam.Aph2_canonicalKt.aph2_canonical(
        jpype.JDouble(a1), jpype.JDouble(t11),
        jpype.JDouble(a2), jpype.JDouble(t22)
    )

    alpha = jlineMatrixToArray(result.getFirst())
    T = jlineMatrixToArray(result.getSecond())

    return alpha, T


def map_cdf_derivative(MAP, x, order=1):
    """
    Compute derivative of MAP cumulative distribution function.

    Args:
        MAP: MAP as tuple (D0, D1)
        x: Point at which to evaluate derivative
        order: Derivative order (default: 1)

    Returns:
        float: CDF derivative value at point x
    """
    if isinstance(MAP, (list, tuple)) and len(MAP) == 2:
        D0, D1 = MAP
        map_cell = jpype.JClass("jline.lang.MatrixCell")(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1)
        )
    else:
        map_cell = MAP

    result = jpype.JPackage('jline').api.mam.Map_cdf_derivativeKt.map_cdf_derivative(
        map_cell, jpype.JDouble(x), jpype.JInt(order)
    )

    return float(result)


def map_rand_moment(K, target_mean=1.0, target_var=1.0):
    """
    Generate random MAP with specified moments.

    Creates a random Markovian Arrival Process with the given
    number of states and target statistical moments.

    Args:
        K: Number of states
        target_mean: Target mean inter-arrival time
        target_var: Target variance of inter-arrival time

    Returns:
        tuple: (D0, D1) - Random MAP matrices with target moments
    """
    result = jpype.JPackage('jline').api.mam.Map_rand_momentKt.map_rand_moment(
        jpype.JInt(K), jpype.JDouble(target_mean), jpype.JDouble(target_var)
    )

    return result


def qbd_solve(A0, A1, A2):
    """
    Solve Quasi-Birth-Death process for stationary distribution.

    Computes the stationary probability distribution and rate matrix
    for an infinite QBD process.

    Args:
        A0: Transition rate matrix for level down
        A1: Transition rate matrix for same level
        A2: Transition rate matrix for level up

    Returns:
        tuple: (pi0, R) - Boundary probabilities and rate matrix
    """
    from .. import jlineMatrixFromArray, jlineMatrixToArray

    java_A0 = jlineMatrixFromArray(A0)
    java_A1 = jlineMatrixFromArray(A1)
    java_A2 = jlineMatrixFromArray(A2)

    java_result = jpype.JPackage('jline').api.mam.Qbd_solveKt.qbd_solve(
        java_A0, java_A1, java_A2
    )

    pi0 = jlineMatrixToArray(java_result.getPi0())
    R = jlineMatrixToArray(java_result.getR())

    return pi0, R



# Additional MAM functions for complete API coverage

def amap2_fit_gamma_map(map_obj):
    """
    Fits AMAP(2) by approximating arbitrary-order MAP with preserved correlation structure.

    Args:
        map_obj: Input MAP object to approximate

    Returns:
        tuple: Fitted AMAP(2) parameters (D0, D1)
    """
    result = jpype.JPackage('jline').api.mam.Amap2_fit_gamma_mapKt.amap2_fit_gamma_map(map_obj)
    return jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1))


def amap2_fit_gamma_trace(trace):
    """
    Fits AMAP(2) from empirical traces while preserving autocorrelation characteristics.

    Args:
        trace: Empirical trace data (array-like)

    Returns:
        tuple: Fitted AMAP(2) parameters (D0, D1)
    """
    trace_matrix = jlineMatrixFromArray(trace)
    result = jpype.JPackage('jline').api.mam.Amap2_fit_gamma_traceKt.amap2_fit_gamma_trace(trace_matrix)
    return jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1))


def aph2_fit_map(map_obj):
    """
    Fits APH(2) distributions by approximating arbitrary-order MAP processes.

    Args:
        map_obj: Input MAP object to approximate with APH(2)

    Returns:
        tuple: Fitted APH(2) parameters (alpha, T)
    """
    result = jpype.JPackage('jline').api.mam.Aph2_fit_mapKt.aph2_fit_map(map_obj)
    return jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1))


def aph2_fit_trace(trace):
    """
    Fits APH(2) distributions from empirical inter-arrival time traces.

    Args:
        trace: Empirical trace data (array-like)

    Returns:
        tuple: Fitted APH(2) parameters (alpha, T)
    """
    trace_matrix = jlineMatrixFromArray(trace)
    result = jpype.JPackage('jline').api.mam.Aph2_fit_traceKt.aph2_fit_trace(trace_matrix)
    return jlineMatrixToArray(result.get(0)), jlineMatrixToArray(result.get(1))


def qbd_r(A0, A1, A2):
    """
    Computes the rate matrix R for QBD processes.
    
    Args:
        A0: Level down transition matrix
        A1: Level transition matrix  
        A2: Level up transition matrix
        
    Returns:
        numpy.ndarray: Rate matrix R
    """
    result = jpype.JPackage('jline').api.mam.Qbd_rKt.qbd_r(
        jlineMatrixFromArray(A0),
        jlineMatrixFromArray(A1),
        jlineMatrixFromArray(A2)
    )
    return jlineMatrixToArray(result)


def map_block(MAP, k):
    """
    Extracts a block from a MAP representation.

    Args:
        MAP: Markovian Arrival Process representation
        k: Block index

    Returns:
        numpy.ndarray: The k-th block of the MAP
    """
    result = jpype.JPackage('jline').api.mam.Map_blockKt.map_block(
        jlineMatrixFromArray(MAP), int(k)
    )
    return jlineMatrixToArray(result)


def map_feasblock(MAP):
    """
    Computes feasibility blocks for MAP representation.

    Args:
        MAP: Markovian Arrival Process representation

    Returns:
        dict: Feasibility block information
    """
    result = jpype.JPackage('jline').api.mam.Map_feasblockKt.map_feasblock(
        jlineMatrixFromArray(MAP)
    )
    return {
        'blocks': jlineMatrixToArray(result.getBlocks()),
        'feasible': bool(result.isFeasible())
    }


def map_kpc(MAP, k):
    """
    Computes k-th order phase-type correlation for MAP.

    Args:
        MAP: Markovian Arrival Process representation
        k: Correlation order

    Returns:
        float: k-th order correlation coefficient
    """
    result = jpype.JPackage('jline').api.mam.Map_kpcKt.map_kpc(
        jlineMatrixFromArray(MAP), int(k)
    )
    return float(result)


def map_pntiter(MAP, n, epsilon=1e-8):
    """
    Point process iteration method for MAP analysis.

    Args:
        MAP: Markovian Arrival Process representation
        n: Number of iterations
        epsilon: Convergence tolerance

    Returns:
        dict: Iteration results
    """
    result = jpype.JPackage('jline').api.mam.Map_pntiterKt.map_pntiter(
        jlineMatrixFromArray(MAP), int(n), float(epsilon)
    )
    return {
        'points': jlineMatrixToArray(result.getPoints()),
        'converged': bool(result.isConverged()),
        'iterations': int(result.getIterations())
    }


def map_pntquad(MAP, n):
    """
    Point process quadrature method for MAP analysis.

    Args:
        MAP: Markovian Arrival Process representation
        n: Number of quadrature points

    Returns:
        dict: Quadrature results
    """
    result = jpype.JPackage('jline').api.mam.Map_pntquadKt.map_pntquad(
        jlineMatrixFromArray(MAP), int(n)
    )
    return {
        'points': jlineMatrixToArray(result.getPoints()),
        'weights': jlineMatrixToArray(result.getWeights())
    }


def map2mmpp(MAP):
    """
    Converts MAP to MMPP representation.

    Args:
        MAP: Markovian Arrival Process representation

    Returns:
        dict: MMPP representation with D0 and D1 matrices
    """
    result = jpype.JPackage('jline').api.mam.Map2mmppKt.map2mmpp(
        jlineMatrixFromArray(MAP)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }


def m3pp_rand(n, params=None):
    """
    Generates random M3PP (third-order Markov-modulated Poisson process).

    Args:
        n: Number of states
        params: Optional parameters for generation

    Returns:
        dict: Random M3PP representation
    """
    if params is not None:
        result = jpype.JPackage('jline').api.mam.M3pp_randKt.m3pp_rand(
            int(n), jlineMatrixFromArray(params)
        )
    else:
        result = jpype.JPackage('jline').api.mam.M3pp_randKt.m3pp_rand(int(n))

    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1()),
        'D2': jlineMatrixToArray(result.getD2())
    }


def m3pp_interleave_fitc_theoretical(M3PP1, M3PP2):
    """
    Fits interleaved M3PP using theoretical approach.

    Args:
        M3PP1: First M3PP process
        M3PP2: Second M3PP process

    Returns:
        dict: Interleaved M3PP representation
    """
    result = jpype.JPackage('jline').api.mam.M3pp_interleave_fitc_theoreticalKt.m3pp_interleave_fitc_theoretical(
        jlineMatrixFromArray(M3PP1), jlineMatrixFromArray(M3PP2)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1()),
        'D2': jlineMatrixToArray(result.getD2())
    }


def m3pp_interleave_fitc_trace(trace1, trace2):
    """
    Fits interleaved M3PP from trace data.

    Args:
        trace1: First trace data
        trace2: Second trace data

    Returns:
        dict: Fitted M3PP representation
    """
    result = jpype.JPackage('jline').api.mam.M3pp_interleave_fitc_traceKt.m3pp_interleave_fitc_trace(
        jlineMatrixFromArray(trace1), jlineMatrixFromArray(trace2)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1()),
        'D2': jlineMatrixToArray(result.getD2()),
        'quality': float(result.getQuality())
    }


def m3pp_superpos_fitc(M3PPs):
    """
    Fits superposition of multiple M3PP processes.

    Args:
        M3PPs: List of M3PP processes

    Returns:
        dict: Superposed M3PP representation
    """
    result = jpype.JPackage('jline').api.mam.M3pp_superpos_fitcKt.m3pp_superpos_fitc(
        jlineMatrixFromArray(M3PPs)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1()),
        'D2': jlineMatrixToArray(result.getD2())
    }


def m3pp_superpos_fitc_theoretical(M3PPs):
    """
    Fits superposition of M3PP processes using theoretical approach.

    Args:
        M3PPs: List of M3PP processes

    Returns:
        dict: Superposed M3PP representation
    """
    result = jpype.JPackage('jline').api.mam.M3pp_superpos_fitc_theoreticalKt.m3pp_superpos_fitc_theoretical(
        jlineMatrixFromArray(M3PPs)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1()),
        'D2': jlineMatrixToArray(result.getD2())
    }


def maph2m_fit(data, m, h):
    """
    Fits PH(2,m) distribution to data.

    Args:
        data: Input data for fitting
        m: Number of phases
        h: Hyper-parameter

    Returns:
        dict: Fitted PH(2,m) representation
    """
    result = jpype.JPackage('jline').api.mam.Maph2m_fitKt.maph2m_fit(
        jlineMatrixFromArray(data), int(m), int(h)
    )
    return {
        'alpha': jlineMatrixToArray(result.getAlpha()),
        'T': jlineMatrixToArray(result.getT()),
        'quality': float(result.getQuality())
    }


def maph2m_fitc_approx(moments, m):
    """
    Approximate fitting of PH(2,m) from moments.

    Args:
        moments: Statistical moments
        m: Number of phases

    Returns:
        dict: Fitted PH(2,m) representation
    """
    result = jpype.JPackage('jline').api.mam.Maph2m_fitc_approxKt.maph2m_fitc_approx(
        jlineMatrixFromArray(moments), int(m)
    )
    return {
        'alpha': jlineMatrixToArray(result.getAlpha()),
        'T': jlineMatrixToArray(result.getT())
    }


def maph2m_fitc_theoretical(params, m):
    """
    Theoretical fitting of PH(2,m) distribution.

    Args:
        params: Theoretical parameters
        m: Number of phases

    Returns:
        dict: Fitted PH(2,m) representation
    """
    result = jpype.JPackage('jline').api.mam.Maph2m_fitc_theoreticalKt.maph2m_fitc_theoretical(
        jlineMatrixFromArray(params), int(m)
    )
    return {
        'alpha': jlineMatrixToArray(result.getAlpha()),
        'T': jlineMatrixToArray(result.getT())
    }


def maph2m_fit_mmap(MMAP, m):
    """
    Fits PH(2,m) from MMAP representation.

    Args:
        MMAP: Marked MAP representation
        m: Number of phases

    Returns:
        dict: Fitted PH(2,m) representation
    """
    result = jpype.JPackage('jline').api.mam.Maph2m_fit_mmapKt.maph2m_fit_mmap(
        jlineMatrixFromArray(MMAP), int(m)
    )
    return {
        'alpha': jlineMatrixToArray(result.getAlpha()),
        'T': jlineMatrixToArray(result.getT())
    }


def maph2m_fit_multiclass(data, classes, m):
    """
    Fits PH(2,m) for multiclass data.

    Args:
        data: Multiclass input data
        classes: Class definitions
        m: Number of phases

    Returns:
        dict: Fitted multiclass PH(2,m) representation
    """
    result = jpype.JPackage('jline').api.mam.Maph2m_fit_multiclassKt.maph2m_fit_multiclass(
        jlineMatrixFromArray(data), jlineMatrixFromArray(classes), int(m)
    )
    return {
        'alpha': jlineMatrixToArray(result.getAlpha()),
        'T': jlineMatrixToArray(result.getT()),
        'class_params': jlineMatrixToArray(result.getClassParams()) if hasattr(result, 'getClassParams') else None
    }


def maph2m_fit_trace(trace, m):
    """
    Fits PH(2,m) from trace data.

    Args:
        trace: Trace data
        m: Number of phases

    Returns:
        dict: Fitted PH(2,m) representation
    """
    result = jpype.JPackage('jline').api.mam.Maph2m_fit_traceKt.maph2m_fit_trace(
        jlineMatrixFromArray(trace), int(m)
    )
    return {
        'alpha': jlineMatrixToArray(result.getAlpha()),
        'T': jlineMatrixToArray(result.getT()),
        'quality': float(result.getQuality())
    }


def mmap_embedded(MMAP):
    """
    Computes embedded process of MMAP.

    Args:
        MMAP: Marked MAP representation

    Returns:
        numpy.ndarray: Embedded process representation
    """
    result = jpype.JPackage('jline').api.mam.Mmap_embeddedKt.mmap_embedded(
        jlineMatrixFromArray(MMAP)
    )
    return jlineMatrixToArray(result)


def mmap_count_mcov(MMAP, lag):
    """
    Computes cross-covariance of counts for MMAP.

    Args:
        MMAP: Marked MAP representation
        lag: Time lag for covariance

    Returns:
        numpy.ndarray: Cross-covariance matrix
    """
    result = jpype.JPackage('jline').api.mam.Mmap_count_mcovKt.mmap_count_mcov(
        jlineMatrixFromArray(MMAP), int(lag)
    )
    return jlineMatrixToArray(result)


def mmap_issym(MMAP):
    """
    Checks if MMAP is symmetric.

    Args:
        MMAP: Marked MAP representation

    Returns:
        bool: True if MMAP is symmetric
    """
    result = jpype.JPackage('jline').api.mam.Mmap_issymKt.mmap_issym(
        jlineMatrixFromArray(MMAP)
    )
    return bool(result)


def mmap_max(MMAPs):
    """
    Computes maximum of multiple MMAPs.

    Args:
        MMAPs: List of MMAP representations

    Returns:
        dict: Maximum MMAP representation
    """
    result = jpype.JPackage('jline').api.mam.Mmap_maxKt.mmap_max(
        jlineMatrixFromArray(MMAPs)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }


def mmap_mixture_order2(MMAPs, weights):
    """
    Creates second-order mixture of MMAPs.

    Args:
        MMAPs: List of MMAP representations
        weights: Mixture weights

    Returns:
        dict: Mixed MMAP representation
    """
    result = jpype.JPackage('jline').api.mam.Mmap_mixture_order2Kt.mmap_mixture_order2(
        jlineMatrixFromArray(MMAPs), jlineMatrixFromArray(weights)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }


def mmap_modulate(MMAP, modulation):
    """
    Modulates MMAP with given modulation function.

    Args:
        MMAP: Marked MAP representation
        modulation: Modulation function or parameters

    Returns:
        dict: Modulated MMAP representation
    """
    result = jpype.JPackage('jline').api.mam.Mmap_modulateKt.mmap_modulate(
        jlineMatrixFromArray(MMAP), jlineMatrixFromArray(modulation)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }


def mmap_pie(MMAP):
    """
    Computes stationary distribution of MMAP.

    Args:
        MMAP: Marked MAP representation

    Returns:
        numpy.ndarray: Stationary distribution
    """
    result = jpype.JPackage('jline').api.mam.Mmap_pieKt.mmap_pie(
        jlineMatrixFromArray(MMAP)
    )
    return jlineMatrixToArray(result)


def mmap_sigma(MMAP):
    """
    Computes sigma parameter for MMAP.

    Args:
        MMAP: Marked MAP representation

    Returns:
        float: Sigma parameter
    """
    result = jpype.JPackage('jline').api.mam.Mmap_sigmaKt.mmap_sigma(
        jlineMatrixFromArray(MMAP)
    )
    return float(result)


def mmap_sum(MMAPs):
    """
    Computes sum of multiple MMAPs.

    Args:
        MMAPs: List of MMAP representations

    Returns:
        dict: Summed MMAP representation
    """
    result = jpype.JPackage('jline').api.mam.Mmap_sumKt.mmap_sum(
        jlineMatrixFromArray(MMAPs)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }


def mmpp2_fitc(data):
    """
    Fits MMPP(2) using correlation fitting.

    Args:
        data: Input data for fitting

    Returns:
        dict: Fitted MMPP(2) representation
    """
    result = jpype.JPackage('jline').api.mam.Mmpp2_fitcKt.mmpp2_fitc(
        jlineMatrixFromArray(data)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1()),
        'quality': float(result.getQuality())
    }


def mmpp2_fitc_approx(moments):
    """
    Approximate fitting of MMPP(2) from moments.

    Args:
        moments: Statistical moments

    Returns:
        dict: Fitted MMPP(2) representation
    """
    result = jpype.JPackage('jline').api.mam.Mmpp2_fitc_approxKt.mmpp2_fitc_approx(
        jlineMatrixFromArray(moments)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }


def qbd_r_logred(A0, A1, A2):
    """
    Computes the rate matrix R using logarithmic reduction.
    
    Args:
        A0: Level down transition matrix
        A1: Level transition matrix
        A2: Level up transition matrix
        
    Returns:
        numpy.ndarray: Rate matrix R computed via logarithmic reduction
    """
    result = jpype.JPackage('jline').api.mam.Qbd_r_logredKt.qbd_r_logred(
        jlineMatrixFromArray(A0),
        jlineMatrixFromArray(A1),
        jlineMatrixFromArray(A2)
    )
    return jlineMatrixToArray(result)



# MAPQN Functions
def mapqn_bnd_lr(MAP, N, Z=None):
    """
    Computes lower and upper response time bounds for MAP queueing networks.

    Args:
        MAP: MAP representation of arrival process
        N: Population vector
        Z: Think times (optional)

    Returns:
        dict: Lower and upper bounds
    """
    if Z is not None:
        result = jpype.JPackage('jline').api.mam.Mapqn_bnd_lrKt.mapqn_bnd_lr(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N), jlineMatrixFromArray(Z)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Mapqn_bnd_lrKt.mapqn_bnd_lr(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N)
        )

    return {
        'lower': jlineMatrixToArray(result.getLower()),
        'upper': jlineMatrixToArray(result.getUpper())
    }


def mapqn_bnd_lr_mva(MAP, N, Z=None):
    """
    Computes MAP queueing network bounds using MVA approach.

    Args:
        MAP: MAP representation of arrival process
        N: Population vector
        Z: Think times (optional)

    Returns:
        dict: MVA-based bounds
    """
    if Z is not None:
        result = jpype.JPackage('jline').api.mam.Mapqn_bnd_lr_mvaKt.mapqn_bnd_lr_mva(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N), jlineMatrixFromArray(Z)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Mapqn_bnd_lr_mvaKt.mapqn_bnd_lr_mva(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N)
        )

    return {
        'lower': jlineMatrixToArray(result.getLower()),
        'upper': jlineMatrixToArray(result.getUpper()),
        'iterations': int(result.getIterations()) if hasattr(result, 'getIterations') else None
    }


def mapqn_bnd_lr_pf(MAP, N, Z=None):
    """
    Computes MAP queueing network bounds using product-form approach.

    Args:
        MAP: MAP representation of arrival process
        N: Population vector
        Z: Think times (optional)

    Returns:
        dict: Product-form based bounds
    """
    if Z is not None:
        result = jpype.JPackage('jline').api.mam.Mapqn_bnd_lr_pfKt.mapqn_bnd_lr_pf(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N), jlineMatrixFromArray(Z)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Mapqn_bnd_lr_pfKt.mapqn_bnd_lr_pf(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N)
        )

    return {
        'lower': jlineMatrixToArray(result.getLower()),
        'upper': jlineMatrixToArray(result.getUpper())
    }


def mapqn_bnd_qr(MAP, N, Z=None):
    """
    Computes queue length and response time bounds for MAP queueing networks.

    Args:
        MAP: MAP representation of arrival process
        N: Population vector
        Z: Think times (optional)

    Returns:
        dict: Queue and response time bounds
    """
    if Z is not None:
        result = jpype.JPackage('jline').api.mam.Mapqn_bnd_qrKt.mapqn_bnd_qr(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N), jlineMatrixFromArray(Z)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Mapqn_bnd_qrKt.mapqn_bnd_qr(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N)
        )

    return {
        'queue_lower': jlineMatrixToArray(result.getQueueLower()),
        'queue_upper': jlineMatrixToArray(result.getQueueUpper()),
        'response_lower': jlineMatrixToArray(result.getResponseLower()),
        'response_upper': jlineMatrixToArray(result.getResponseUpper())
    }


def mapqn_bnd_qr_delay(MAP, N, Z, delay):
    """
    Computes bounds for MAP queueing networks with delay.

    Args:
        MAP: MAP representation of arrival process
        N: Population vector
        Z: Think times
        delay: Delay parameters

    Returns:
        dict: Bounds with delay consideration
    """
    result = jpype.JPackage('jline').api.mam.Mapqn_bnd_qr_delayKt.mapqn_bnd_qr_delay(
        jlineMatrixFromArray(MAP), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jlineMatrixFromArray(delay)
    )

    return {
        'queue_lower': jlineMatrixToArray(result.getQueueLower()),
        'queue_upper': jlineMatrixToArray(result.getQueueUpper()),
        'response_lower': jlineMatrixToArray(result.getResponseLower()),
        'response_upper': jlineMatrixToArray(result.getResponseUpper())
    }


def mapqn_bnd_qr_ld(MAP, N, Z, mu):
    """
    Computes bounds for load-dependent MAP queueing networks.

    Args:
        MAP: MAP representation of arrival process
        N: Population vector
        Z: Think times
        mu: Load-dependent service rates

    Returns:
        dict: Load-dependent bounds
    """
    result = jpype.JPackage('jline').api.mam.Mapqn_bnd_qr_ldKt.mapqn_bnd_qr_ld(
        jlineMatrixFromArray(MAP), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jlineMatrixFromArray(mu)
    )

    return {
        'queue_lower': jlineMatrixToArray(result.getQueueLower()),
        'queue_upper': jlineMatrixToArray(result.getQueueUpper()),
        'response_lower': jlineMatrixToArray(result.getResponseLower()),
        'response_upper': jlineMatrixToArray(result.getResponseUpper())
    }


def mapqn_lpmodel(MAP, N, Z=None):
    """
    Creates linear programming model for MAP queueing network.

    Args:
        MAP: MAP representation of arrival process
        N: Population vector
        Z: Think times (optional)

    Returns:
        dict: LP model components
    """
    if Z is not None:
        result = jpype.JPackage('jline').api.mam.Mapqn_lpmodelKt.mapqn_lpmodel(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N), jlineMatrixFromArray(Z)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Mapqn_lpmodelKt.mapqn_lpmodel(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N)
        )

    return {
        'A': jlineMatrixToArray(result.getA()),
        'b': jlineMatrixToArray(result.getB()),
        'c': jlineMatrixToArray(result.getC()),
        'bounds': jlineMatrixToArray(result.getBounds()) if hasattr(result, 'getBounds') else None
    }


def mapqn_parameters(MAP, N):
    """
    Extracts parameters from MAP queueing network.

    Args:
        MAP: MAP representation of arrival process
        N: Population vector

    Returns:
        dict: Network parameters
    """
    result = jpype.JPackage('jline').api.mam.Mapqn_parametersKt.mapqn_parameters(
        jlineMatrixFromArray(MAP), jlineMatrixFromArray(N)
    )

    return {
        'arrival_rate': float(result.getArrivalRate()),
        'service_rates': jlineMatrixToArray(result.getServiceRates()),
        'routing': jlineMatrixToArray(result.getRouting()),
        'think_times': jlineMatrixToArray(result.getThinkTimes()) if hasattr(result, 'getThinkTimes') else None
    }


def mapqn_parameters_factory(params_dict):
    """
    Factory method to create MAP queueing network parameters.

    Args:
        params_dict: Dictionary of parameters

    Returns:
        dict: Formatted network parameters
    """
    result = jpype.JPackage('jline').api.mam.Mapqn_parameters_factoryKt.mapqn_parameters_factory(
        jpype.JObject(params_dict)
    )

    return {
        'MAP': jlineMatrixToArray(result.getMAP()),
        'N': jlineMatrixToArray(result.getN()),
        'Z': jlineMatrixToArray(result.getZ()) if hasattr(result, 'getZ') else None
    }


def mapqn_qr_bounds_bas(MAP, N, Z=None):
    """
    Computes Balanced Asymptotic System (BAS) bounds for MAP queueing networks.

    Args:
        MAP: MAP representation of arrival process
        N: Population vector
        Z: Think times (optional)

    Returns:
        dict: BAS bounds
    """
    if Z is not None:
        result = jpype.JPackage('jline').api.mam.Mapqn_qr_bounds_basKt.mapqn_qr_bounds_bas(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N), jlineMatrixFromArray(Z)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Mapqn_qr_bounds_basKt.mapqn_qr_bounds_bas(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N)
        )

    return {
        'queue_bounds': jlineMatrixToArray(result.getQueueBounds()),
        'response_bounds': jlineMatrixToArray(result.getResponseBounds())
    }


def mapqn_qr_bounds_rsrd(MAP, N, Z=None):
    """
    Computes Response-time Scaled Routing Delay (RSRD) bounds for MAP queueing networks.

    Args:
        MAP: MAP representation of arrival process
        N: Population vector
        Z: Think times (optional)

    Returns:
        dict: RSRD bounds
    """
    if Z is not None:
        result = jpype.JPackage('jline').api.mam.Mapqn_qr_bounds_rsrdKt.mapqn_qr_bounds_rsrd(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N), jlineMatrixFromArray(Z)
        )
    else:
        result = jpype.JPackage('jline').api.mam.Mapqn_qr_bounds_rsrdKt.mapqn_qr_bounds_rsrd(
            jlineMatrixFromArray(MAP), jlineMatrixFromArray(N)
        )

    return {
        'queue_bounds': jlineMatrixToArray(result.getQueueBounds()),
        'response_bounds': jlineMatrixToArray(result.getResponseBounds())
    }


# M3PP Functions (Markovian Multi-class Point Processes)

def m3pp22_fitc_approx_cov(mean, scv, skew, cov):
    """
    Implements parameter fitting for second-order Marked Markov Modulated Poisson Process.

    Args:
        mean: Mean inter-arrival time
        scv: Squared coefficient of variation
        skew: Skewness
        cov: Covariance matrix

    Returns:
        dict: Fitted M3PP(2,2) parameters
    """
    result = jpype.JPackage('jline').api.mam.m3pp.M3pp22_fitc_approx_covKt.m3pp22_fitc_approx_cov(
        float(mean), float(scv), float(skew), arrayToJLineMatrix(cov)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1()),
        'params': jlineMatrixToArray(result.getParams())
    }


# M3PP and MAMAP Functions
# Merged from mam_extended.py

def mamap2m_coefficients(mean, scv, skew):
    """
    Computes coefficients for MAMAP(2,m) fitting.

    Args:
        mean: Mean value
        scv: Squared coefficient of variation
        skew: Skewness

    Returns:
        numpy.ndarray: Coefficient matrix
    """
    result = jpype.JPackage('jline').api.mam.Mamap2m_coefficientsKt.mamap2m_coefficients(
        float(mean), float(scv), float(skew)
    )
    return jlineMatrixToArray(result)





def m3pp22_fitc_approx_cov_multiclass(means, scvs, skews, covs):
    """
    Implements constrained optimization for fitting M3PP(2,2) parameters given an underlying.

    Args:
        means: Mean inter-arrival times per class
        scvs: Squared coefficients of variation per class
        skews: Skewness values per class
        covs: Covariance matrices

    Returns:
        dict: Fitted M3PP(2,2) parameters for multiple classes
    """
    result = jpype.JPackage('jline').api.mam.m3pp.M3pp22_fitc_approx_cov_multiclassKt.m3pp22_fitc_approx_cov_multiclass(
        arrayToJLineMatrix(means),
        arrayToJLineMatrix(scvs),
        arrayToJLineMatrix(skews),
        arrayToJLineMatrix(covs)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1()),
        'params': jlineMatrixToArray(result.getParams())
    }


def m3pp22_interleave_fitc(processes):
    """
    Implements lumped superposition of multiple M3PP(2,2) processes using interleaved.

    Args:
        processes: List of M3PP(2,2) process parameters

    Returns:
        dict: Interleaved M3PP(2,2) parameters
    """
    java_processes = jpype.JPackage('java.util').ArrayList()
    for proc in processes:
        java_processes.add(proc)

    result = jpype.JPackage('jline').api.mam.m3pp.M3pp22_interleave_fitcKt.m3pp22_interleave_fitc(java_processes)
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }


def m3pp2m_fitc(mean, scv, skew, num_classes):
    """
    Implements exact fitting of second-order Marked Markov Modulated Poisson Process.

    Args:
        mean: Mean inter-arrival time
        scv: Squared coefficient of variation
        skew: Skewness
        num_classes: Number of classes

    Returns:
        dict: Fitted M3PP(2,m) parameters
    """
    result = jpype.JPackage('jline').api.mam.m3pp.M3pp2m_fitcKt.m3pp2m_fitc(
        float(mean), float(scv), float(skew), int(num_classes)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1()),
        'params': jlineMatrixToArray(result.getParams())
    }


def m3pp2m_fitc_approx(mean, scv, skew, num_classes):
    """
    Implements approximation-based fitting for M3PP(2,m) using optimization methods.

    Args:
        mean: Mean inter-arrival time
        scv: Squared coefficient of variation
        skew: Skewness
        num_classes: Number of classes

    Returns:
        dict: Fitted M3PP(2,m) parameters
    """
    result = jpype.JPackage('jline').api.mam.m3pp.M3pp2m_fitc_approxKt.m3pp2m_fitc_approx(
        float(mean), float(scv), float(skew), int(num_classes)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1()),
        'params': jlineMatrixToArray(result.getParams())
    }


def m3pp2m_fitc_approx_ag(mean, scv, skew, num_classes):
    """
    Implements auto-gamma approximation method for M3PP(2,m) parameter fitting.

    Args:
        mean: Mean inter-arrival time
        scv: Squared coefficient of variation
        skew: Skewness
        num_classes: Number of classes

    Returns:
        dict: Fitted M3PP(2,m) parameters using auto-gamma
    """
    result = jpype.JPackage('jline').api.mam.m3pp.M3pp2m_fitc_approx_agKt.m3pp2m_fitc_approx_ag(
        float(mean), float(scv), float(skew), int(num_classes)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1()),
        'params': jlineMatrixToArray(result.getParams())
    }


def m3pp2m_fitc_approx_ag_multiclass(means, scvs, skews, num_classes):
    """
    Implements multiclass auto-gamma fitting for M3PP(2,m) with variance and covariance.

    Args:
        means: Mean inter-arrival times per class
        scvs: Squared coefficients of variation per class
        skews: Skewness values per class
        num_classes: Number of classes

    Returns:
        dict: Fitted M3PP(2,m) parameters for multiple classes
    """
    result = jpype.JPackage('jline').api.mam.m3pp.M3pp2m_fitc_approx_ag_multiclassKt.m3pp2m_fitc_approx_ag_multiclass(
        arrayToJLineMatrix(means),
        arrayToJLineMatrix(scvs),
        arrayToJLineMatrix(skews),
        int(num_classes)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1()),
        'params': jlineMatrixToArray(result.getParams())
    }


def m3pp2m_interleave(processes):
    """
    Implements interleaved superposition of multiple M3PP(2,m) processes to construct.

    Args:
        processes: List of M3PP(2,m) process parameters

    Returns:
        dict: Interleaved M3PP(2,m) parameters
    """
    java_processes = jpype.JPackage('java.util').ArrayList()
    for proc in processes:
        java_processes.add(proc)

    result = jpype.JPackage('jline').api.mam.m3pp.M3pp2m_interleaveKt.m3pp2m_interleave(java_processes)
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }


def m3pp_interleave_fitc(processes):
    """
    Implements fitting and interleaving of k second-order M3PP processes with varying.

    Args:
        processes: List of M3PP process parameters

    Returns:
        dict: Interleaved M3PP parameters
    """
    java_processes = jpype.JPackage('java.util').ArrayList()
    for proc in processes:
        java_processes.add(proc)

    result = jpype.JPackage('jline').api.mam.m3pp.M3pp_interleave_fitcKt.m3pp_interleave_fitc(java_processes)
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }


# Additional MAMAP functions
def mamap22_fit_gamma_fs_trace(trace):
    """
    Fits MAMAP(2,2) from trace data using gamma autocorrelation and forward-sigma characteristics.

    Args:
        trace: Empirical trace data

    Returns:
        dict: Fitted MAMAP(2,2) parameters
    """
    trace_matrix = arrayToJLineMatrix(trace)
    result = jpype.JPackage('jline').api.mam.Mamap22_fit_gamma_fs_traceKt.mamap22_fit_gamma_fs_trace(trace_matrix)
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }


def mamap22_fit_multiclass(class_data):
    """
    Fits MAMAP(2,2) processes for two-class systems with forward moments and sigma characteristics.

    Args:
        class_data: Data for multiple classes

    Returns:
        dict: Fitted MAMAP(2,2) parameters for multiple classes
    """
    result = jpype.JPackage('jline').api.mam.Mamap22_fit_multiclassKt.mamap22_fit_multiclass(class_data)
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }


def mamap2m_fit(mean, scv, skew, num_states):
    """
    Fits MAMAP(2,m) processes matching moments, autocorrelation, and class characteristics.

    Args:
        mean: Mean value
        scv: Squared coefficient of variation
        skew: Skewness
        num_states: Number of states

    Returns:
        dict: Fitted MAMAP(2,m) parameters
    """
    result = jpype.JPackage('jline').api.mam.Mamap2m_fitKt.mamap2m_fit(
        float(mean), float(scv), float(skew), int(num_states)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }


def mamap2m_fit_mmap(mmap_obj):
    """
    Fits MAPH/MAMAP(2,m) by approximating characteristics of input MMAP processes.

    Args:
        mmap_obj: Input MMAP object

    Returns:
        dict: Fitted MAMAP(2,m) parameters
    """
    result = jpype.JPackage('jline').api.mam.Mamap2m_fit_mmapKt.mamap2m_fit_mmap(mmap_obj)
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }


def mamap2m_fit_trace(trace, num_states):
    """
    Fits MAMAP(2,m) processes from empirical trace data with inter-arrival times and class labels.

    Args:
        trace: Empirical trace data
        num_states: Number of states

    Returns:
        dict: Fitted MAMAP(2,m) parameters
    """
    trace_matrix = arrayToJLineMatrix(trace)
    result = jpype.JPackage('jline').api.mam.Mamap2m_fit_traceKt.mamap2m_fit_trace(
        trace_matrix, int(num_states)
    )
    return {
        'D0': jlineMatrixToArray(result.getD0()),
        'D1': jlineMatrixToArray(result.getD1())
    }
