"""
Moment Conversion Algorithms for Discrete Distributions.

Native Python implementations of the one-to-one conversions between the power
(raw), factorial, upward-factorial, binomial, negative-binomial and central
moments of a discrete random variable, together with the combinatorial
triangles (Stirling numbers of the first and second kind, Stirling cycle
numbers, Lah numbers) that underpin them.

References:
    A. Heindl and A. van de Liefvoort. Moment conversions for discrete
    distributions. PMCCS, 2003.

    I. Lah. Eine neue Art von Zahlen, ihre Eigenschaften und Anwendung in der
    mathematischen Statistik. Mitteilungsbl. Math. Statist., 7:203-212, 1955.
"""

from math import factorial

import numpy as np
from scipy.special import comb


def _check_order(n, name):
    """
    Validate that n is a nonnegative integer maximum order.

    Args:
        n: Candidate maximum order.
        name: Name of the calling function, used in the error message.

    Returns:
        The order as a Python int.

    Raises:
        ValueError: If n is not a nonnegative integer scalar.
    """
    if isinstance(n, (bool, np.bool_)):
        raise ValueError('%s: The maximum order n must be a nonnegative integer.' % name)
    arr = np.asarray(n)
    if arr.ndim != 0:
        raise ValueError('%s: The maximum order n must be a nonnegative integer.' % name)
    val = arr.item()
    if not isinstance(val, (int, np.integer)):
        if not float(val).is_integer():
            raise ValueError('%s: The maximum order n must be a nonnegative integer.' % name)
    val = int(val)
    if val < 0:
        raise ValueError('%s: The maximum order n must be a nonnegative integer.' % name)
    return val


def _as_moment_vector(x, name, argname):
    """
    Coerce a moment sequence into a 1-D float array.

    Args:
        x: Array-like of length n+1 holding the moments of order 0,...,n.
        name: Name of the calling function, used in the error message.
        argname: Name of the argument, used in the error message.

    Returns:
        1-D numpy array of floats.

    Raises:
        ValueError: If x is not a nonempty vector.
    """
    arr = np.asarray(x, dtype=float).ravel()
    if arr.size < 1:
        raise ValueError('%s: The input %s must be a nonempty vector of moments.' % (name, argname))
    return arr


def _binom(n, k):
    """
    Exact binomial coefficient nchoosek(n,k) returned as a float.

    Args:
        n: Upper index (nonnegative integer).
        k: Lower index (nonnegative integer).

    Returns:
        The binomial coefficient as a float.
    """
    return float(comb(int(n), int(k), exact=True))


def moment_stirlingcycle(n):
    """
    Triangle of the Stirling cycle numbers.

    The Stirling cycle numbers (unsigned Stirling numbers of the first kind)
    satisfy sigma(i,j) = (-1)^(i-j) * s(i,j) and are obtained from the
    recursion

        sigma(i,j) = (i-1)*sigma(i-1,j) + sigma(i-1,j-1)   for j > 0
        sigma(0,0) = 1,   sigma(i,0) = 0 for i > 0

    These numbers are the coefficients that convert power moments into
    upward-factorial moments.

    Args:
        n: Maximum order (n >= 0).

    Returns:
        (n+1)x(n+1) numpy array with element (i,j) equal to sigma(i,j) in the
        0-based notation of the reference. Entries with j > i are zero.

    Raises:
        ValueError: If n is not a nonnegative integer.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, eq. (12).

    Example:
        sigma = moment_stirlingcycle(3)
    """
    n = _check_order(n, 'moment_stirlingcycle')
    sigma = np.zeros((n + 1, n + 1))
    sigma[0, 0] = 1.0
    for i in range(1, n + 1):
        for j in range(1, i + 1):
            sigma[i, j] = (i - 1) * sigma[i - 1, j] + sigma[i - 1, j - 1]
    return sigma


def moment_stirling1(n):
    """
    Triangle of the signed Stirling numbers of the first kind.

    The signed Stirling numbers of the first kind s(i,j) are defined as the
    coefficients of x^j in the falling factorial

        sum_{j=0}^{i} s(i,j) x^j = x(x-1)(x-2)...(x-i+1)

    These numbers are the coefficients that convert power moments into
    factorial moments. They relate to the Stirling cycle numbers via
    s(i,j) = (-1)^(i-j) * sigma(i,j).

    Args:
        n: Maximum order (n >= 0).

    Returns:
        (n+1)x(n+1) numpy array with element (i,j) equal to s(i,j) in the
        0-based notation of the reference. Entries with j > i are zero.

    Raises:
        ValueError: If n is not a nonnegative integer.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, eq. (10) and eq. (12).

    Example:
        s = moment_stirling1(3)
    """
    n = _check_order(n, 'moment_stirling1')
    sigma = moment_stirlingcycle(n)
    s = np.zeros((n + 1, n + 1))
    for i in range(n + 1):
        for j in range(i + 1):
            s[i, j] = ((-1) ** (i - j)) * sigma[i, j]
    return s


def moment_stirling2(n):
    """
    Triangle of the Stirling numbers of the second kind.

    The Stirling numbers of the second kind S(i,j) are implicitly defined by
    the expansion of a power into falling factorials

        x^i = sum_{j=0}^{i} S(i,j) x(x-1)(x-2)...(x-j+1)

    and computed from the recursion

        S(i,j) = j*S(i-1,j) + S(i-1,j-1)   for j > 0
        S(0,0) = 1,   S(i,0) = 0 for i > 0

    These numbers are the coefficients that convert factorial moments back
    into power moments.

    Args:
        n: Maximum order (n >= 0).

    Returns:
        (n+1)x(n+1) numpy array with element (i,j) equal to S(i,j) in the
        0-based notation of the reference. Entries with j > i are zero.

    Raises:
        ValueError: If n is not a nonnegative integer.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, eq. (11).

    Example:
        S = moment_stirling2(3)
    """
    n = _check_order(n, 'moment_stirling2')
    S = np.zeros((n + 1, n + 1))
    S[0, 0] = 1.0
    for i in range(1, n + 1):
        for j in range(1, i + 1):
            S[i, j] = j * S[i - 1, j] + S[i - 1, j - 1]
    return S


def moment_lah(n):
    """
    Triangle of the Lah numbers.

    The Lah numbers L(i,j) = (i!/j!)*nchoosek(i-1,j-1) link the factorial
    moments to the upward-factorial moments. The triangle is built from the
    equivalent recursion

        L(i,j) = L(i-1,j-1) + (i+j-1)*L(i-1,j)   for j > 0
        L(0,0) = 1,   L(i,0) = 0 for i > 0

    which avoids the overflow of the explicit factorial form for large orders.

    Args:
        n: Maximum order (n >= 0).

    Returns:
        (n+1)x(n+1) numpy array with element (i,j) equal to L(i,j) in the
        0-based notation of the reference. Entries with j > i are zero.

    Raises:
        ValueError: If n is not a nonnegative integer.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, Section 4.

        I. Lah. Eine neue Art von Zahlen, ihre Eigenschaften und Anwendung in
        der mathematischen Statistik. Mitteilungsbl. Math. Statist.,
        7:203-212, 1955.

    Example:
        L = moment_lah(3)
    """
    n = _check_order(n, 'moment_lah')
    L = np.zeros((n + 1, n + 1))
    L[0, 0] = 1.0
    for i in range(1, n + 1):
        for j in range(1, i + 1):
            L[i, j] = L[i - 1, j - 1] + (i + j - 1) * L[i - 1, j]
    return L


def moment_binotrans(x):
    """
    Binomial transform of a sequence.

    The binomial transform maps the sequence x_0,x_1,...,x_n into
    y_0,y_1,...,y_n via

        y_n = sum_{k=0}^{n} (-1)^(n-k) * nchoosek(n,k) * x_k

    Applied to a moment sequence m_i = E[X^i] it returns the moments of the
    unit downshift, y_i = E[(X-1)^i]. It is not an involution: its inverse is
    moment_binotransinv, the unsigned transform.

    Args:
        x: Array-like of length n+1 holding x_0,...,x_n, i.e. element i is the
           element of order i.

    Returns:
        1-D numpy array of length n+1 holding y_0,...,y_n.

    Raises:
        ValueError: If x is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, eq. (8).

    Example:
        y = moment_binotrans([1, 2, 5, 15])
    """
    xcol = _as_moment_vector(x, 'moment_binotrans', 'x')
    n = xcol.size - 1
    y = np.zeros(n + 1)
    for i in range(n + 1):
        for k in range(i + 1):
            y[i] += ((-1) ** (i - k)) * _binom(i, k) * xcol[k]
    return y


def moment_binotransinv(y):
    """
    Inverse binomial transform of a sequence.

    The inverse binomial transform maps the sequence y_0,y_1,...,y_n into
    x_0,x_1,...,x_n via

        x_n = sum_{k=0}^{n} nchoosek(n,k) * y_k

    Args:
        y: Array-like of length n+1 holding y_0,...,y_n, i.e. element i is the
           element of order i.

    Returns:
        1-D numpy array of length n+1 holding x_0,...,x_n.

    Raises:
        ValueError: If y is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, eq. (9).

    Example:
        x = moment_binotransinv(moment_binotrans([1, 2, 5, 15]))
    """
    ycol = _as_moment_vector(y, 'moment_binotransinv', 'y')
    n = ycol.size - 1
    x = np.zeros(n + 1)
    for i in range(n + 1):
        for k in range(i + 1):
            x[i] += _binom(i, k) * ycol[k]
    return x


def moment_factorial_from_raw(m):
    """
    Convert power (raw) moments into factorial moments.

    The power moments m_n = E[N^n] of a discrete random variable N are
    converted into the factorial moments f_n = E[N(N-1)...(N-n+1)] by means of
    the signed Stirling numbers of the first kind,

        f_n = sum_{k=0}^{n} s(n,k) * m_k

    Args:
        m: Array-like of length n+1 holding m_0,...,m_n, i.e. element i is the
           moment of order i and element 0 is m_0 = 1.

    Returns:
        1-D numpy array of length n+1 holding f_0,...,f_n.

    Raises:
        ValueError: If m is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, eq. (13).

    Example:
        f = moment_factorial_from_raw([1, 2, 6, 22])
    """
    mcol = _as_moment_vector(m, 'moment_factorial_from_raw', 'm')
    n = mcol.size - 1
    return moment_stirling1(n) @ mcol


def moment_raw_from_factorial(f):
    """
    Convert factorial moments into power (raw) moments.

    The factorial moments f_n = E[N(N-1)...(N-n+1)] of a discrete random
    variable N are converted into the power moments m_n = E[N^n] by means of
    the Stirling numbers of the second kind,

        m_n = sum_{k=0}^{n} S(n,k) * f_k

    Args:
        f: Array-like of length n+1 holding f_0,...,f_n, i.e. element i is the
           moment of order i and element 0 is f_0 = 1.

    Returns:
        1-D numpy array of length n+1 holding m_0,...,m_n.

    Raises:
        ValueError: If f is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, eq. (13).

    Example:
        m = moment_raw_from_factorial(moment_factorial_from_raw([1, 2, 6, 22]))
    """
    fcol = _as_moment_vector(f, 'moment_raw_from_factorial', 'f')
    n = fcol.size - 1
    return moment_stirling2(n) @ fcol


def moment_upfactorial_from_raw(m):
    """
    Convert power (raw) moments into upward-factorial moments.

    The power moments m_n = E[N^n] of a discrete random variable N are
    converted into the upward-factorial moments f_n^+ = E[N(N+1)...(N+n-1)] by
    means of the Stirling cycle numbers,

        f_n^+ = sum_{k=0}^{n} sigma(n,k) * m_k

    Upward-factorial moments are of use in moment-matching techniques for
    matrix-geometric and discrete phase-type distributions.

    Args:
        m: Array-like of length n+1 holding m_0,...,m_n, i.e. element i is the
           moment of order i and element 0 is m_0 = 1.

    Returns:
        1-D numpy array of length n+1 holding f_0^+,...,f_n^+.

    Raises:
        ValueError: If m is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, Section 4.

    Example:
        fp = moment_upfactorial_from_raw([1, 2, 6, 22])
    """
    mcol = _as_moment_vector(m, 'moment_upfactorial_from_raw', 'm')
    n = mcol.size - 1
    return moment_stirlingcycle(n) @ mcol


def moment_raw_from_upfactorial(fp):
    """
    Convert upward-factorial moments into power (raw) moments.

    The upward-factorial moments f_n^+ = E[N(N+1)...(N+n-1)] of a discrete
    random variable N are converted into the power moments m_n = E[N^n] by
    means of the signed Stirling numbers of the second kind,

        m_n = sum_{k=0}^{n} (-1)^(n-k) * S(n,k) * f_k^+

    Args:
        fp: Array-like of length n+1 holding f_0^+,...,f_n^+, i.e. element i is
            the moment of order i and element 0 is f_0^+ = 1.

    Returns:
        1-D numpy array of length n+1 holding m_0,...,m_n.

    Raises:
        ValueError: If fp is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, Section 4.

    Example:
        m = moment_raw_from_upfactorial(moment_upfactorial_from_raw([1, 2, 6, 22]))
    """
    fpcol = _as_moment_vector(fp, 'moment_raw_from_upfactorial', 'fp')
    n = fpcol.size - 1
    S = moment_stirling2(n)
    T = np.zeros((n + 1, n + 1))
    for i in range(n + 1):
        for j in range(i + 1):
            T[i, j] = ((-1) ** (i - j)) * S[i, j]
    return T @ fpcol


def moment_binomial_from_factorial(f):
    """
    Convert factorial moments into binomial moments.

    The factorial moments f_n = E[N(N-1)...(N-n+1)] of a discrete random
    variable N are converted into the binomial moments b_n = E[nchoosek(N,n)]
    via the one-to-one correspondence

        b_n = f_n / n!

    Args:
        f: Array-like of length n+1 holding f_0,...,f_n, i.e. element i is the
           moment of order i and element 0 is f_0 = 1.

    Returns:
        1-D numpy array of length n+1 holding b_0,...,b_n.

    Raises:
        ValueError: If f is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, Section 4.

    Example:
        b = moment_binomial_from_factorial([1, 2, 4, 8])
    """
    fcol = _as_moment_vector(f, 'moment_binomial_from_factorial', 'f')
    n = fcol.size - 1
    b = np.zeros(n + 1)
    for i in range(n + 1):
        b[i] = fcol[i] / float(factorial(i))
    return b


def moment_factorial_from_binomial(b):
    """
    Convert binomial moments into factorial moments.

    The binomial moments b_n = E[nchoosek(N,n)] of a discrete random variable N
    are converted into the factorial moments f_n = E[N(N-1)...(N-n+1)] via the
    one-to-one correspondence

        f_n = n! * b_n

    Args:
        b: Array-like of length n+1 holding b_0,...,b_n, i.e. element i is the
           moment of order i and element 0 is b_0 = 1.

    Returns:
        1-D numpy array of length n+1 holding f_0,...,f_n.

    Raises:
        ValueError: If b is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, Section 4.

    Example:
        f = moment_factorial_from_binomial([1, 2, 2, 4/3])
    """
    bcol = _as_moment_vector(b, 'moment_factorial_from_binomial', 'b')
    n = bcol.size - 1
    f = np.zeros(n + 1)
    for i in range(n + 1):
        f[i] = float(factorial(i)) * bcol[i]
    return f


def moment_negbinomial_from_upfactorial(fp):
    """
    Convert upward-factorial moments into negative-binomial moments.

    The upward-factorial moments f_n^+ = E[N(N+1)...(N+n-1)] of a discrete
    random variable N are converted into the negative-binomial moments
    b_n^- = E[nchoosek(N+n-1,n)] via the one-to-one correspondence

        b_n^- = f_n^+ / n!

    Args:
        fp: Array-like of length n+1 holding f_0^+,...,f_n^+, i.e. element i is
            the moment of order i and element 0 is f_0^+ = 1.

    Returns:
        1-D numpy array of length n+1 holding b_0^-,...,b_n^-.

    Raises:
        ValueError: If fp is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, eq. (7).

    Example:
        bm = moment_negbinomial_from_upfactorial([1, 2, 8, 44])
    """
    fpcol = _as_moment_vector(fp, 'moment_negbinomial_from_upfactorial', 'fp')
    n = fpcol.size - 1
    bm = np.zeros(n + 1)
    for i in range(n + 1):
        bm[i] = fpcol[i] / float(factorial(i))
    return bm


def moment_upfactorial_from_negbinomial(bm):
    """
    Convert negative-binomial moments into upward-factorial moments.

    The negative-binomial moments b_n^- = E[nchoosek(N+n-1,n)] of a discrete
    random variable N are converted into the upward-factorial moments
    f_n^+ = E[N(N+1)...(N+n-1)] via the one-to-one correspondence

        f_n^+ = n! * b_n^-

    Args:
        bm: Array-like of length n+1 holding b_0^-,...,b_n^-, i.e. element i is
            the moment of order i and element 0 is b_0^- = 1.

    Returns:
        1-D numpy array of length n+1 holding f_0^+,...,f_n^+.

    Raises:
        ValueError: If bm is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, eq. (7).

    Example:
        fp = moment_upfactorial_from_negbinomial([1, 2, 4, 22/3])
    """
    bmcol = _as_moment_vector(bm, 'moment_upfactorial_from_negbinomial', 'bm')
    n = bmcol.size - 1
    fp = np.zeros(n + 1)
    for i in range(n + 1):
        fp[i] = float(factorial(i)) * bmcol[i]
    return fp


def moment_binomial_from_negbinomial(bm):
    """
    Convert negative-binomial moments into binomial moments.

    The negative-binomial moments b_n^- = E[nchoosek(N+n-1,n)] of a discrete
    random variable N are converted into the binomial moments
    b_n = E[nchoosek(N,n)] by means of the shifted binomial transform

        b_n = sum_{k=1}^{n} (-1)^(n-k) * nchoosek(n-1,k-1) * b_k^-   for n >= 1
        b_0 = 1

    Args:
        bm: Array-like of length n+1 holding b_0^-,...,b_n^-, i.e. element i is
            the moment of order i and element 0 is b_0^- = 1.

    Returns:
        1-D numpy array of length n+1 holding b_0,...,b_n.

    Raises:
        ValueError: If bm is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, eq. (14).

    Example:
        b = moment_binomial_from_negbinomial([1, 2, 4, 22/3])
    """
    bmcol = _as_moment_vector(bm, 'moment_binomial_from_negbinomial', 'bm')
    n = bmcol.size - 1
    b = np.zeros(n + 1)
    b[0] = 1.0
    for i in range(1, n + 1):
        for k in range(1, i + 1):
            b[i] += ((-1) ** (i - k)) * _binom(i - 1, k - 1) * bmcol[k]
    return b


def moment_negbinomial_from_binomial(b):
    """
    Convert binomial moments into negative-binomial moments.

    The binomial moments b_n = E[nchoosek(N,n)] of a discrete random variable N
    are converted into the negative-binomial moments
    b_n^- = E[nchoosek(N+n-1,n)] by means of the shifted binomial transform

        b_n^- = sum_{k=1}^{n} nchoosek(n-1,k-1) * b_k   for n >= 1
        b_0^- = 1

    Args:
        b: Array-like of length n+1 holding b_0,...,b_n, i.e. element i is the
           moment of order i and element 0 is b_0 = 1.

    Returns:
        1-D numpy array of length n+1 holding b_0^-,...,b_n^-.

    Raises:
        ValueError: If b is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, eq. (14).

    Example:
        bm = moment_negbinomial_from_binomial([1, 2, 2, 4/3])
    """
    bcol = _as_moment_vector(b, 'moment_negbinomial_from_binomial', 'b')
    n = bcol.size - 1
    bm = np.zeros(n + 1)
    bm[0] = 1.0
    for i in range(1, n + 1):
        for k in range(1, i + 1):
            bm[i] += _binom(i - 1, k - 1) * bcol[k]
    return bm


def moment_factorial_from_upfactorial(fp):
    """
    Convert upward-factorial moments into factorial moments.

    The upward-factorial moments f_n^+ = E[N(N+1)...(N+n-1)] of a discrete
    random variable N are converted into the factorial moments
    f_n = E[N(N-1)...(N-n+1)] by means of the Lah numbers,

        f_n = sum_{k=1}^{n} (-1)^(n-k) * L(n,k) * f_k^+   for n >= 1
        f_0 = 1

    Args:
        fp: Array-like of length n+1 holding f_0^+,...,f_n^+, i.e. element i is
            the moment of order i and element 0 is f_0^+ = 1.

    Returns:
        1-D numpy array of length n+1 holding f_0,...,f_n.

    Raises:
        ValueError: If fp is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, Section 4.

    Example:
        f = moment_factorial_from_upfactorial([1, 2, 8, 44])
    """
    fpcol = _as_moment_vector(fp, 'moment_factorial_from_upfactorial', 'fp')
    n = fpcol.size - 1
    L = moment_lah(n)
    f = np.zeros(n + 1)
    f[0] = 1.0
    for i in range(1, n + 1):
        for k in range(1, i + 1):
            f[i] += ((-1) ** (i - k)) * L[i, k] * fpcol[k]
    return f


def moment_upfactorial_from_factorial(f):
    """
    Convert factorial moments into upward-factorial moments.

    The factorial moments f_n = E[N(N-1)...(N-n+1)] of a discrete random
    variable N are converted into the upward-factorial moments
    f_n^+ = E[N(N+1)...(N+n-1)] by means of the Lah numbers,

        f_n^+ = sum_{k=1}^{n} L(n,k) * f_k   for n >= 1
        f_0^+ = 1

    Args:
        f: Array-like of length n+1 holding f_0,...,f_n, i.e. element i is the
           moment of order i and element 0 is f_0 = 1.

    Returns:
        1-D numpy array of length n+1 holding f_0^+,...,f_n^+.

    Raises:
        ValueError: If f is not a nonempty vector.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, Section 4.

    Example:
        fp = moment_upfactorial_from_factorial([1, 2, 4, 8])
    """
    fcol = _as_moment_vector(f, 'moment_upfactorial_from_factorial', 'f')
    n = fcol.size - 1
    L = moment_lah(n)
    fp = np.zeros(n + 1)
    fp[0] = 1.0
    for i in range(1, n + 1):
        for k in range(1, i + 1):
            fp[i] += L[i, k] * fcol[k]
    return fp


def moment_central_from_raw(m):
    """
    Convert power (raw) moments into central moments.

    The power moments m_n = E[N^n] of a random variable N are converted into
    the central moments m_n^c = E[(N-m_1)^n] by means of the binomial transform
    in the variation that involves the mean m_1,

        m_n^c = sum_{k=0}^{n} (-1)^(n-k) * nchoosek(n,k) * m_k * m_1^(n-k)

    The conversion also holds for continuous random variables.

    Args:
        m: Array-like of length n+1 holding m_0,...,m_n, i.e. element i is the
           moment of order i and element 0 is m_0 = 1. At least the mean m_1
           must be given, hence n >= 1.

    Returns:
        1-D numpy array of length n+1 holding m_0^c,...,m_n^c. By construction
        m_0^c = 1 and m_1^c = 0.

    Raises:
        ValueError: If m has fewer than 2 elements, i.e. the mean m_1 is
            missing.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, Section 4.

    Example:
        mc = moment_central_from_raw([1, 2, 6, 22])
    """
    mcol = _as_moment_vector(m, 'moment_central_from_raw', 'm')
    n = mcol.size - 1
    if n < 1:
        raise ValueError('moment_central_from_raw: The mean m_1 is required for '
                         'this conversion, hence m must have at least 2 elements.')
    m1 = mcol[1]
    mc = np.zeros(n + 1)
    for i in range(n + 1):
        for k in range(i + 1):
            mc[i] += ((-1) ** (i - k)) * _binom(i, k) * mcol[k] * m1 ** (i - k)
    return mc


def moment_raw_from_central(mc, m1):
    """
    Convert central moments into power (raw) moments.

    The central moments m_n^c = E[(N-m_1)^n] of a random variable N are
    converted into the power moments m_n = E[N^n] by means of the inverse
    binomial transform in the variation that involves the mean m_1,

        m_n = sum_{k=0}^{n} nchoosek(n,k) * m_k^c * m_1^(n-k)

    The mean must be supplied separately since m_1^c = 0 carries no information
    on it. The conversion also holds for continuous random variables.

    Args:
        mc: Array-like of length n+1 holding m_0^c,...,m_n^c, i.e. element i is
            the moment of order i and element 0 is m_0^c = 1.
        m1: Mean of N (scalar).

    Returns:
        1-D numpy array of length n+1 holding m_0,...,m_n.

    Raises:
        ValueError: If mc is not a nonempty vector or if m1 is not a scalar.

    References:
        A. Heindl and A. van de Liefvoort. Moment conversions for discrete
        distributions. PMCCS, 2003, Section 4.

    Example:
        m = moment_raw_from_central(moment_central_from_raw([1, 2, 6, 22]), 2)
    """
    mccol = _as_moment_vector(mc, 'moment_raw_from_central', 'mc')
    n = mccol.size - 1
    if np.asarray(m1).ndim != 0:
        raise ValueError('moment_raw_from_central: The mean m1 must be a scalar.')
    m1 = float(m1)
    m = np.zeros(n + 1)
    for i in range(n + 1):
        for k in range(i + 1):
            m[i] += _binom(i, k) * mccol[k] * m1 ** (i - k)
    return m


# ---------------------------------------------------------------------------
# Cumulants and factorial cumulants (univariate)
# ---------------------------------------------------------------------------

def _cumulants_from_moments(mu):
    """
    Cumulant sequence of a moment sequence, from the exponential-formula
    recursion mu_n = sum_{k=1}^{n} nchoosek(n-1,k-1) kappa_k mu_(n-k).

    Args:
        mu: 1-D array of moments of order 0,...,n with mu_0 = 1.

    Returns:
        1-D numpy array of cumulants of order 0,...,n, with element 0 set to 0.
    """
    n = mu.size - 1
    kappa = np.zeros(n + 1)
    for i in range(1, n + 1):
        acc = 0.0
        for k in range(1, i):
            acc += _binom(i - 1, k - 1) * kappa[k] * mu[i - k]
        kappa[i] = mu[i] - acc
    return kappa


def _moments_from_cumulants(kappa):
    """
    Moment sequence of a cumulant sequence, from the same recursion solved in
    the forward direction.

    Args:
        kappa: 1-D array of cumulants of order 0,...,n; element 0 is ignored.

    Returns:
        1-D numpy array of moments of order 0,...,n, with element 0 set to 1.
    """
    n = kappa.size - 1
    mu = np.zeros(n + 1)
    mu[0] = 1.0
    for i in range(1, n + 1):
        acc = 0.0
        for k in range(1, i + 1):
            acc += _binom(i - 1, k - 1) * kappa[k] * mu[i - k]
        mu[i] = acc
    return mu


def moment_cumulant_from_raw(m):
    """
    Convert power (raw) moments into cumulants.

    The cumulants kappa_n of a random variable X are the coefficients of the
    cumulant generating function log E[exp(sX)] = sum_{n>=1} kappa_n s^n / n!.
    They are obtained from the power moments m_n = E[X^n] by the
    exponential-formula recursion

        m_n = sum_{k=1}^{n} nchoosek(n-1,k-1) * kappa_k * m_(n-k)

    solved for kappa_n. Equivalently kappa_n = sum_{pi in P(n)} (|pi|-1)!
    (-1)^(|pi|-1) prod_{B in pi} m_|B| over the set partitions of {1,...,n}.
    The first cumulants are kappa_1 = m_1, kappa_2 = m_2 - m_1^2 (the variance)
    and kappa_3 = m_3 - 3 m_1 m_2 + 2 m_1^3 (the third central moment).

    The conversion is not restricted to discrete random variables.

    Args:
        m: Array-like of length n+1 holding m_0,...,m_n, i.e. element i is the
           moment of order i and element 0 is m_0 = 1.

    Returns:
        1-D numpy array of length n+1 holding kappa_0,...,kappa_n. Element 0 is
        kappa_0 = 0, the value of the cumulant generating function at the
        origin, and not m_0 = 1.

    Raises:
        ValueError: If m is not a nonempty vector.

    References:
        V. P. Leonov and A. N. Shiryaev. On a method of calculation of
        semi-invariants. Theory of Probability and its Applications,
        4(3):319-329, 1959.

    Example:
        kappa = moment_cumulant_from_raw([1, 2, 6, 22])
    """
    mcol = _as_moment_vector(m, 'moment_cumulant_from_raw', 'm')
    return _cumulants_from_moments(mcol)


def moment_raw_from_cumulant(kappa):
    """
    Convert cumulants into power (raw) moments.

    Inverts moment_cumulant_from_raw by running the exponential-formula
    recursion forward,

        m_n = sum_{k=1}^{n} nchoosek(n-1,k-1) * kappa_k * m_(n-k)

    with m_0 = 1. Equivalently m_n = sum_{pi in P(n)} prod_{B in pi}
    kappa_|B| over the set partitions of {1,...,n}.

    Args:
        kappa: Array-like of length n+1 holding kappa_0,...,kappa_n. Element 0
               is ignored, since kappa_0 = 0 carries no information.

    Returns:
        1-D numpy array of length n+1 holding m_0,...,m_n, with m_0 = 1.

    Raises:
        ValueError: If kappa is not a nonempty vector.

    References:
        V. P. Leonov and A. N. Shiryaev. On a method of calculation of
        semi-invariants. Theory of Probability and its Applications,
        4(3):319-329, 1959.

    Example:
        m = moment_raw_from_cumulant(moment_cumulant_from_raw([1, 2, 6, 22]))
    """
    kcol = _as_moment_vector(kappa, 'moment_raw_from_cumulant', 'kappa')
    return _moments_from_cumulants(kcol)


def moment_factcumulant_from_factorial(f):
    """
    Convert factorial moments into factorial cumulants.

    The factorial cumulants kappa_n^[ ] of a discrete random variable N are the
    coefficients of the logarithm of the probability generating function
    expanded about z = 1,

        log E[z^N] = sum_{n>=1} kappa_n^[ ] (z-1)^n / n!

    They stand to the factorial moments f_n = E[N(N-1)...(N-n+1)] exactly as
    the cumulants stand to the power moments, so the same recursion applies,

        f_n = sum_{k=1}^{n} nchoosek(n-1,k-1) * kappa_k^[ ] * f_(n-k)

    For a Poisson variable of rate lambda all factorial cumulants beyond the
    first vanish, which makes them the natural measure of departure from
    Poisson behaviour in the counting process of a MAP.

    Args:
        f: Array-like of length n+1 holding f_0,...,f_n, i.e. element i is the
           moment of order i and element 0 is f_0 = 1.

    Returns:
        1-D numpy array of length n+1 holding kappa_0^[ ],...,kappa_n^[ ], with
        element 0 equal to 0.

    Raises:
        ValueError: If f is not a nonempty vector.

    Example:
        kf = moment_factcumulant_from_factorial([1, 2, 4, 8])
    """
    fcol = _as_moment_vector(f, 'moment_factcumulant_from_factorial', 'f')
    return _cumulants_from_moments(fcol)


def moment_factorial_from_factcumulant(kappa):
    """
    Convert factorial cumulants into factorial moments.

    Inverts moment_factcumulant_from_factorial by running the recursion
    forward,

        f_n = sum_{k=1}^{n} nchoosek(n-1,k-1) * kappa_k^[ ] * f_(n-k)

    with f_0 = 1.

    Args:
        kappa: Array-like of length n+1 holding kappa_0^[ ],...,kappa_n^[ ].
               Element 0 is ignored.

    Returns:
        1-D numpy array of length n+1 holding f_0,...,f_n, with f_0 = 1.

    Raises:
        ValueError: If kappa is not a nonempty vector.

    Example:
        f = moment_factorial_from_factcumulant([0, 2, 0, 0])
    """
    kcol = _as_moment_vector(kappa, 'moment_factorial_from_factcumulant', 'kappa')
    return _moments_from_cumulants(kcol)


# ---------------------------------------------------------------------------
# Multivariate (joint) moment conversions
# ---------------------------------------------------------------------------

def _as_moment_tensor(x, name, argname):
    """
    Coerce a joint moment array into a float ndarray.

    Args:
        x: Array-like of shape (n_1+1,...,n_d+1) holding the joint moments,
           element (i_1,...,i_d) being the moment of multi-order (i_1,...,i_d).
        name: Name of the calling function, used in the error message.
        argname: Name of the argument, used in the error message.

    Returns:
        ndarray of floats with at least one dimension.

    Raises:
        ValueError: If x is empty.
    """
    arr = np.asarray(x, dtype=float)
    if arr.ndim == 0 or arr.size < 1:
        raise ValueError('%s: The input %s must be a nonempty array of joint moments.'
                         % (name, argname))
    return arr


def moment_tensortrans(A, T, mode):
    """
    Apply a conversion matrix along one dimension of a joint moment array.

    This is the mode product of the array with the matrix: every fibre of A
    along the given dimension is replaced by T times that fibre. Applying it
    once per dimension realises the Kronecker product of the univariate
    conversions, which is the structure of every separable edge of the house.

    Args:
        A: Array-like of joint moments, of shape (n_1+1,...,n_d+1).
        T: (n_mode+1)x(n_mode+1) conversion matrix.
        mode: Zero-based dimension to transform.

    Returns:
        ndarray of the same shape as A.

    Raises:
        ValueError: If mode is out of range or T does not match the extent of
            that dimension.

    Example:
        B = moment_tensortrans(A, moment_stirling1(A.shape[0] - 1), 0)
    """
    arr = _as_moment_tensor(A, 'moment_tensortrans', 'A')
    if mode < 0 or mode >= arr.ndim:
        raise ValueError('moment_tensortrans: The mode must be in 0,...,%d.' % (arr.ndim - 1))
    Tm = np.asarray(T, dtype=float)
    if Tm.ndim != 2 or Tm.shape[1] != arr.shape[mode]:
        raise ValueError('moment_tensortrans: The matrix T must have %d columns.'
                         % arr.shape[mode])
    return np.moveaxis(np.tensordot(Tm, arr, axes=([1], [mode])), 0, mode)


def moment_jointtrans(A, edge):
    """
    Apply the conversion matrix of one edge of the house of moments along every
    dimension of a joint moment array.

    This is the separable (Kronecker product) form shared by all the joint
    conversions except the cumulant and the central ones.

    Args:
        A: Array-like of joint moments, of shape (n_1+1,...,n_d+1).
        edge: Edge label accepted by moment_housematrix.

    Returns:
        ndarray of the same shape as A holding the converted joint moments.

    Raises:
        ValueError: If A is empty or edge is not a known label.

    Example:
        f = moment_jointtrans(m, 'factorial_from_raw')
    """
    arr = _as_moment_tensor(A, 'moment_jointtrans', 'A')
    out = arr
    for mode in range(arr.ndim):
        out = moment_tensortrans(out, moment_housematrix(edge, out.shape[mode] - 1), mode)
    return out


def moment_housematrix(edge, n):
    """
    Conversion matrix of one edge of the house of moments.

    The edge is returned as a linear map on the moment subspace {m_0 = 1}. Four
    edges (the Lah pair and the shifted-binomial pair) pin their zeroth output
    to 1 rather than propagating element 0, so as maps of the whole space they
    are affine. Here the offset is folded into column 0, which is empty for
    those edges, making every edge a genuine matrix. On a moment vector, whose
    element 0 is 1 by definition, the two agree. This is also what makes those
    edges usable dimension by dimension in the joint conversions.

    Args:
        edge: One of 'factorial_from_raw', 'raw_from_factorial',
              'binomial_from_tail', 'tail_from_binomial',
              'upfactorial_from_raw', 'raw_from_upfactorial',
              'binomial_from_factorial', 'factorial_from_binomial',
              'negbinomial_from_upfactorial', 'upfactorial_from_negbinomial',
              'factorial_from_upfactorial', 'upfactorial_from_factorial',
              'negbinomial_from_binomial', 'binomial_from_negbinomial'.
        n: Maximum order of the mode (n >= 0).

    Returns:
        (n+1)x(n+1) numpy array.

    Raises:
        ValueError: If edge is not a known label or n is not a nonnegative
            integer.

    Example:
        T = moment_housematrix('factorial_from_raw', 4)
    """
    n = _check_order(n, 'moment_housematrix')
    if edge == 'factorial_from_raw':
        return moment_stirling1(n)
    if edge == 'raw_from_factorial':
        return moment_stirling2(n)
    if edge == 'upfactorial_from_raw':
        return moment_stirlingcycle(n)
    if edge == 'raw_from_upfactorial':
        S = moment_stirling2(n)
        return np.array([[((-1) ** (i - j)) * S[i, j] if j <= i else 0.0
                          for j in range(n + 1)] for i in range(n + 1)])
    if edge in ('binomial_from_factorial', 'negbinomial_from_upfactorial'):
        return np.diag([1.0 / factorial(i) for i in range(n + 1)])
    if edge in ('factorial_from_binomial', 'upfactorial_from_negbinomial'):
        return np.diag([float(factorial(i)) for i in range(n + 1)])
    if edge in ('factorial_from_upfactorial', 'upfactorial_from_factorial'):
        L = moment_lah(n)
        T = np.zeros((n + 1, n + 1))
        T[0, 0] = 1.0
        for i in range(1, n + 1):
            for k in range(1, i + 1):
                T[i, k] = (L[i, k] if edge == 'upfactorial_from_factorial'
                           else ((-1) ** (i - k)) * L[i, k])
        return T
    if edge in ('negbinomial_from_binomial', 'binomial_from_negbinomial'):
        T = np.zeros((n + 1, n + 1))
        T[0, 0] = 1.0
        for i in range(1, n + 1):
            for k in range(1, i + 1):
                c = _binom(i - 1, k - 1)
                T[i, k] = (c if edge == 'negbinomial_from_binomial'
                           else ((-1) ** (i - k)) * c)
        return T
    if edge in ('binomial_from_tail', 'tail_from_binomial'):
        T = np.zeros((n + 1, n + 1))
        T[0, 0] = 1.0
        for i in range(1, n + 1):
            for k in range(i, n + 1):
                c = _binom(k - 1, i - 1)
                T[i, k] = c if edge == 'binomial_from_tail' else ((-1) ** (k - i)) * c
        return T
    raise ValueError('moment_housematrix: Unknown edge %s.' % edge)


def moment_joint_factorial_from_raw(m):
    """
    Convert joint power (raw) moments into joint factorial moments.

    For a random vector (N_1,...,N_d) the joint power moments
    m_(i_1,...,i_d) = E[prod_j N_j^(i_j)] are converted into the joint
    factorial moments f_(i_1,...,i_d) = E[prod_j (N_j)_(i_j)], where
    (N)_i = N(N-1)...(N-i+1), by applying the signed Stirling numbers of the
    first kind separately along every dimension,

        f_(i_1,...,i_d) = sum_(k_1,...,k_d) prod_j s(i_j,k_j) * m_(k_1,...,k_d)

    The joint conversion is the Kronecker product of the univariate ones, which
    is what makes the mode-by-mode evaluation legitimate. Only the cumulant and
    the central conversions are not of this separable form.

    Args:
        m: Array-like of shape (n_1+1,...,n_d+1) holding the joint power
           moments, element (i_1,...,i_d) being the moment of multi-order
           (i_1,...,i_d) and element 0 being 1.

    Returns:
        ndarray of the same shape holding the joint factorial moments.

    Raises:
        ValueError: If m is empty.

    Example:
        f = moment_joint_factorial_from_raw(np.ones((3, 3)))
    """
    arr = _as_moment_tensor(m, 'moment_joint_factorial_from_raw', 'm')
    return moment_jointtrans(arr, 'factorial_from_raw')


def moment_joint_raw_from_factorial(f):
    """
    Convert joint factorial moments into joint power (raw) moments.

    Inverse of moment_joint_factorial_from_raw, driven by the Stirling numbers
    of the second kind along every dimension.

    Args:
        f: Array-like of shape (n_1+1,...,n_d+1) of joint factorial moments.

    Returns:
        ndarray of the same shape holding the joint power moments.

    Raises:
        ValueError: If f is empty.
    """
    arr = _as_moment_tensor(f, 'moment_joint_raw_from_factorial', 'f')
    return moment_jointtrans(arr, 'raw_from_factorial')


def moment_joint_upfactorial_from_raw(m):
    """
    Convert joint power (raw) moments into joint upward-factorial moments.

    Driven by the Stirling cycle numbers along every dimension, giving
    f+_(i_1,...,i_d) = E[prod_j N_j(N_j+1)...(N_j+i_j-1)].

    Args:
        m: Array-like of shape (n_1+1,...,n_d+1) of joint power moments.

    Returns:
        ndarray of the same shape holding the joint upward-factorial moments.

    Raises:
        ValueError: If m is empty.
    """
    arr = _as_moment_tensor(m, 'moment_joint_upfactorial_from_raw', 'm')
    return moment_jointtrans(arr, 'upfactorial_from_raw')


def moment_joint_raw_from_upfactorial(fp):
    """
    Convert joint upward-factorial moments into joint power (raw) moments.

    Inverse of moment_joint_upfactorial_from_raw, driven by the signed Stirling
    numbers of the second kind along every dimension.

    Args:
        fp: Array-like of shape (n_1+1,...,n_d+1) of joint upward-factorial
            moments.

    Returns:
        ndarray of the same shape holding the joint power moments.

    Raises:
        ValueError: If fp is empty.
    """
    arr = _as_moment_tensor(fp, 'moment_joint_raw_from_upfactorial', 'fp')
    return moment_jointtrans(arr, 'raw_from_upfactorial')


def moment_joint_binomial_from_factorial(f):
    """
    Convert joint factorial moments into joint binomial moments.

    b_(i_1,...,i_d) = E[prod_j nchoosek(N_j,i_j)] = f_(i_1,...,i_d) / prod_j
    (i_j!).

    Args:
        f: Array-like of shape (n_1+1,...,n_d+1) of joint factorial moments.

    Returns:
        ndarray of the same shape holding the joint binomial moments.

    Raises:
        ValueError: If f is empty.
    """
    arr = _as_moment_tensor(f, 'moment_joint_binomial_from_factorial', 'f')
    return moment_jointtrans(arr, 'binomial_from_factorial')


def moment_joint_factorial_from_binomial(b):
    """
    Convert joint binomial moments into joint factorial moments.

    Inverse of moment_joint_binomial_from_factorial.

    Args:
        b: Array-like of shape (n_1+1,...,n_d+1) of joint binomial moments.

    Returns:
        ndarray of the same shape holding the joint factorial moments.

    Raises:
        ValueError: If b is empty.
    """
    arr = _as_moment_tensor(b, 'moment_joint_factorial_from_binomial', 'b')
    return moment_jointtrans(arr, 'factorial_from_binomial')


def moment_joint_negbinomial_from_upfactorial(fp):
    """
    Convert joint upward-factorial moments into joint negative-binomial
    moments.

    b-_(i_1,...,i_d) = E[prod_j nchoosek(N_j+i_j-1,i_j)] = f+_(i_1,...,i_d) /
    prod_j (i_j!).

    Args:
        fp: Array-like of shape (n_1+1,...,n_d+1) of joint upward-factorial
            moments.

    Returns:
        ndarray of the same shape holding the joint negative-binomial moments.

    Raises:
        ValueError: If fp is empty.
    """
    arr = _as_moment_tensor(fp, 'moment_joint_negbinomial_from_upfactorial', 'fp')
    return moment_jointtrans(arr, 'negbinomial_from_upfactorial')


def moment_joint_upfactorial_from_negbinomial(bm):
    """
    Convert joint negative-binomial moments into joint upward-factorial
    moments.

    Inverse of moment_joint_negbinomial_from_upfactorial.

    Args:
        bm: Array-like of shape (n_1+1,...,n_d+1) of joint negative-binomial
            moments.

    Returns:
        ndarray of the same shape holding the joint upward-factorial moments.

    Raises:
        ValueError: If bm is empty.
    """
    arr = _as_moment_tensor(bm, 'moment_joint_upfactorial_from_negbinomial', 'bm')
    return moment_jointtrans(arr, 'upfactorial_from_negbinomial')


def moment_joint_factorial_from_upfactorial(fp):
    """
    Convert joint upward-factorial moments into joint factorial moments.

    Driven by the signed Lah numbers along every dimension.

    Args:
        fp: Array-like of shape (n_1+1,...,n_d+1) of joint upward-factorial
            moments.

    Returns:
        ndarray of the same shape holding the joint factorial moments.

    Raises:
        ValueError: If fp is empty.
    """
    arr = _as_moment_tensor(fp, 'moment_joint_factorial_from_upfactorial', 'fp')
    return moment_jointtrans(arr, 'factorial_from_upfactorial')


def moment_joint_upfactorial_from_factorial(f):
    """
    Convert joint factorial moments into joint upward-factorial moments.

    Driven by the Lah numbers along every dimension.

    Args:
        f: Array-like of shape (n_1+1,...,n_d+1) of joint factorial moments.

    Returns:
        ndarray of the same shape holding the joint upward-factorial moments.

    Raises:
        ValueError: If f is empty.
    """
    arr = _as_moment_tensor(f, 'moment_joint_upfactorial_from_factorial', 'f')
    return moment_jointtrans(arr, 'upfactorial_from_factorial')


def moment_joint_negbinomial_from_binomial(b):
    """
    Convert joint binomial moments into joint negative-binomial moments.

    Driven by the shifted binomial transform nchoosek(i-1,k-1) along every
    dimension.

    Args:
        b: Array-like of shape (n_1+1,...,n_d+1) of joint binomial moments.

    Returns:
        ndarray of the same shape holding the joint negative-binomial moments.

    Raises:
        ValueError: If b is empty.
    """
    arr = _as_moment_tensor(b, 'moment_joint_negbinomial_from_binomial', 'b')
    return moment_jointtrans(arr, 'negbinomial_from_binomial')


def moment_joint_binomial_from_negbinomial(bm):
    """
    Convert joint negative-binomial moments into joint binomial moments.

    Inverse of moment_joint_negbinomial_from_binomial, driven by the signed
    shifted binomial transform along every dimension.

    Args:
        bm: Array-like of shape (n_1+1,...,n_d+1) of joint negative-binomial
            moments.

    Returns:
        ndarray of the same shape holding the joint binomial moments.

    Raises:
        ValueError: If bm is empty.
    """
    arr = _as_moment_tensor(bm, 'moment_joint_binomial_from_negbinomial', 'bm')
    return moment_jointtrans(arr, 'binomial_from_negbinomial')


def moment_joint_central_from_raw(m):
    """
    Convert joint power (raw) moments into joint central moments.

    The joint central moments mc_(i_1,...,i_d) = E[prod_j (N_j - E N_j)^(i_j)]
    follow from the multi-index binomial theorem, which is again separable but
    with a different shift per dimension,

        mc_(i) = sum_(k<=i) prod_j (-1)^(i_j-k_j) nchoosek(i_j,k_j)
                 mu_j^(i_j-k_j) * m_(k)

    The means mu_j = m_(e_j) are read off the array itself, so every dimension
    must carry at least the first order. The entry of multi-order e_j+e_l is
    the covariance of N_j and N_l.

    Args:
        m: Array-like of shape (n_1+1,...,n_d+1) of joint power moments, with
           every n_j >= 1.

    Returns:
        ndarray of the same shape holding the joint central moments.

    Raises:
        ValueError: If m is empty or has a dimension of extent 1.

    Example:
        mc = moment_joint_central_from_raw(m)
    """
    arr = _as_moment_tensor(m, 'moment_joint_central_from_raw', 'm')
    if min(arr.shape) < 2:
        raise ValueError('moment_joint_central_from_raw: The means m_(e_j) are '
                         'required for this conversion, hence every dimension of m '
                         'must have at least 2 elements.')
    mu = []
    for j in range(arr.ndim):
        idx = [0] * arr.ndim
        idx[j] = 1
        mu.append(float(arr[tuple(idx)]))
    return moment_joint_central_from_raw_mean(arr, mu)


def moment_joint_central_from_raw_mean(m, mu):
    """
    Convert joint power (raw) moments into joint central moments about a given
    mean vector.

    Same conversion as moment_joint_central_from_raw, with the means supplied
    rather than read off the array, so that it also applies when the array does
    not carry the first-order entries.

    Args:
        m: Array-like of shape (n_1+1,...,n_d+1) of joint power moments.
        mu: Array-like of length d holding the means E[N_1],...,E[N_d].

    Returns:
        ndarray of the same shape holding the joint central moments.

    Raises:
        ValueError: If m is empty or if mu does not have one entry per
            dimension of m.
    """
    arr = _as_moment_tensor(m, 'moment_joint_central_from_raw_mean', 'm')
    muv = np.asarray(mu, dtype=float).ravel()
    if muv.size != arr.ndim:
        raise ValueError('moment_joint_central_from_raw_mean: The mean vector mu '
                         'must have one entry per dimension of m.')
    out = arr
    for mode in range(arr.ndim):
        n = out.shape[mode] - 1
        T = np.array([[_binom(i, k) * ((-muv[mode]) ** (i - k)) if k <= i else 0.0
                       for k in range(n + 1)] for i in range(n + 1)])
        out = moment_tensortrans(out, T, mode)
    return out


def moment_joint_raw_from_central(mc, mu):
    """
    Convert joint central moments into joint power (raw) moments.

    Inverse of moment_joint_central_from_raw, again by the multi-index binomial
    theorem,

        m_(i) = sum_(k<=i) prod_j nchoosek(i_j,k_j) mu_j^(i_j-k_j) * mc_(k)

    The mean vector must be supplied separately, since the first-order central
    moments are zero and carry no information on it.

    Args:
        mc: Array-like of shape (n_1+1,...,n_d+1) of joint central moments.
        mu: Array-like of length d holding the means E[N_1],...,E[N_d].

    Returns:
        ndarray of the same shape holding the joint power moments.

    Raises:
        ValueError: If mc is empty or if mu does not have one entry per
            dimension of mc.
    """
    arr = _as_moment_tensor(mc, 'moment_joint_raw_from_central', 'mc')
    muv = np.asarray(mu, dtype=float).ravel()
    if muv.size != arr.ndim:
        raise ValueError('moment_joint_raw_from_central: The mean vector mu must '
                         'have one entry per dimension of mc.')
    out = arr
    for mode in range(arr.ndim):
        n = out.shape[mode] - 1
        T = np.array([[_binom(i, k) * (muv[mode] ** (i - k)) if k <= i else 0.0
                       for k in range(n + 1)] for i in range(n + 1)])
        out = moment_tensortrans(out, T, mode)
    return out


def _joint_cumulants_from_moments(mu):
    """
    Joint cumulant array of a joint moment array, from the multivariate
    exponential-formula recursion.

    With j the first dimension in which the multi-index a is nonzero,

        mu_a = sum_(0<b<=a) prod_l nchoosek(a_l-[l=j], b_l-[l=j])
               kappa_b mu_(a-b)

    which isolates kappa_a because the b = a term has unit coefficient and
    mu_0 = 1. Multi-indices are swept in lexicographic order, under which every
    b <= a precedes a.

    Args:
        mu: ndarray of joint moments with element 0 equal to 1.

    Returns:
        ndarray of the same shape holding the joint cumulants, element 0 being
        0.
    """
    kappa = np.zeros(mu.shape)
    for a in np.ndindex(*mu.shape):
        if sum(a) == 0:
            continue
        j = next(l for l in range(len(a)) if a[l] > 0)
        acc = 0.0
        for b in np.ndindex(*[ai + 1 for ai in a]):
            if sum(b) == 0 or b == a or b[j] < 1:
                continue
            c = 1.0
            for l in range(len(a)):
                c *= _binom(a[l] - (1 if l == j else 0), b[l] - (1 if l == j else 0))
            if c == 0.0:
                continue
            acc += c * kappa[b] * mu[tuple(np.subtract(a, b))]
        kappa[a] = mu[a] - acc
    return kappa


def _joint_moments_from_cumulants(kappa):
    """
    Joint moment array of a joint cumulant array, from the same recursion run
    forward.

    Args:
        kappa: ndarray of joint cumulants; element 0 is ignored.

    Returns:
        ndarray of the same shape holding the joint moments, element 0 being 1.
    """
    mu = np.zeros(kappa.shape)
    mu[tuple([0] * kappa.ndim)] = 1.0
    for a in np.ndindex(*kappa.shape):
        if sum(a) == 0:
            continue
        j = next(l for l in range(len(a)) if a[l] > 0)
        acc = 0.0
        for b in np.ndindex(*[ai + 1 for ai in a]):
            if sum(b) == 0 or b[j] < 1:
                continue
            c = 1.0
            for l in range(len(a)):
                c *= _binom(a[l] - (1 if l == j else 0), b[l] - (1 if l == j else 0))
            if c == 0.0:
                continue
            acc += c * kappa[b] * mu[tuple(np.subtract(a, b))]
        mu[a] = acc
    return mu


def moment_joint_cumulant_from_raw(m):
    """
    Convert joint power (raw) moments into joint cumulants.

    The joint cumulants of a random vector (N_1,...,N_d) are the coefficients
    of the joint cumulant generating function

        log E[exp(s_1 N_1 + ... + s_d N_d)] = sum_(a != 0) kappa_a prod_j
        s_j^(a_j) / a_j!

    They obey the multivariate exponential formula, equivalently the
    Leonov-Shiryaev partition formula. Unlike every other conversion in the
    house, this one does not factor into a product of univariate transforms:
    the cumulant of multi-order (1,1) is the covariance, which mixes the
    dimensions.

    Args:
        m: Array-like of shape (n_1+1,...,n_d+1) of joint power moments, with
           element 0 equal to 1.

    Returns:
        ndarray of the same shape holding the joint cumulants, element 0 being
        kappa_0 = 0.

    Raises:
        ValueError: If m is empty.

    References:
        V. P. Leonov and A. N. Shiryaev. On a method of calculation of
        semi-invariants. Theory of Probability and its Applications,
        4(3):319-329, 1959.

    Example:
        kappa = moment_joint_cumulant_from_raw(m)  # kappa[1,1] is the covariance
    """
    arr = _as_moment_tensor(m, 'moment_joint_cumulant_from_raw', 'm')
    return _joint_cumulants_from_moments(arr)


def moment_joint_raw_from_cumulant(kappa):
    """
    Convert joint cumulants into joint power (raw) moments.

    Inverse of moment_joint_cumulant_from_raw.

    Args:
        kappa: Array-like of shape (n_1+1,...,n_d+1) of joint cumulants;
               element 0 is ignored.

    Returns:
        ndarray of the same shape holding the joint power moments, element 0
        being 1.

    Raises:
        ValueError: If kappa is empty.
    """
    arr = _as_moment_tensor(kappa, 'moment_joint_raw_from_cumulant', 'kappa')
    return _joint_moments_from_cumulants(arr)


def moment_joint_factcumulant_from_factorial(f):
    """
    Convert joint factorial moments into joint factorial cumulants.

    The joint factorial cumulants are the coefficients of the logarithm of the
    joint probability generating function expanded about z = (1,...,1),

        log E[prod_j z_j^(N_j)] = sum_(a != 0) kappa_a^[ ] prod_j
        (z_j-1)^(a_j) / a_j!

    and stand to the joint factorial moments exactly as the joint cumulants
    stand to the joint power moments, so the same recursion applies. For a
    multivariate Poisson vector with independent components every joint
    factorial cumulant of order two or more vanishes; for the per-class counts
    of a marked MAP they measure the departure from independent Poisson marking.

    Args:
        f: Array-like of shape (n_1+1,...,n_d+1) of joint factorial moments,
           with element 0 equal to 1.

    Returns:
        ndarray of the same shape holding the joint factorial cumulants,
        element 0 being 0.

    Raises:
        ValueError: If f is empty.
    """
    arr = _as_moment_tensor(f, 'moment_joint_factcumulant_from_factorial', 'f')
    return _joint_cumulants_from_moments(arr)


def moment_joint_factorial_from_factcumulant(kappa):
    """
    Convert joint factorial cumulants into joint factorial moments.

    Inverse of moment_joint_factcumulant_from_factorial.

    Args:
        kappa: Array-like of shape (n_1+1,...,n_d+1) of joint factorial
               cumulants; element 0 is ignored.

    Returns:
        ndarray of the same shape holding the joint factorial moments, element
        0 being 1.

    Raises:
        ValueError: If kappa is empty.
    """
    arr = _as_moment_tensor(kappa, 'moment_joint_factorial_from_factcumulant', 'kappa')
    return _joint_moments_from_cumulants(arr)


def moment_joint_marking(f, p, dims):
    """
    Joint factorial moments of the per-class counts under multinomial marking.

    If a count N is marked independently, every event receiving class j with
    probability p_j, then the per-class counts (N_1,...,N_d) have joint
    factorial moments

        E[prod_j (N_j)_(a_j)] = (prod_j p_j^(a_j)) * f_(|a|)

    where f is the factorial moment sequence of the aggregate count N and
    |a| = a_1+...+a_d. This is the counting-process counterpart of the marking
    (class-splitting) formulas of the M3A fitters, and it is exact for the
    per-class counts of a MAP marked in this i.i.d. way, in particular for an
    MMAP whose marking probabilities do not depend on the phase.

    Args:
        f: Array-like of length n+1 holding the factorial moments f_0,...,f_n
           of the aggregate count.
        p: Array-like of length d holding the marking probabilities.
        dims: Array-like of length d holding the maximum order per class. Their
              sum must not exceed n, since an entry of multi-order a consumes
              the aggregate moment of order |a|.

    Returns:
        ndarray of shape (dims_1+1,...,dims_d+1) holding the joint factorial
        moments of the per-class counts.

    Raises:
        ValueError: If f is empty, if p and dims have different lengths, or if
            the requested orders exceed the aggregate moments supplied.

    Example:
        F = moment_joint_marking([1, 2, 4, 8], [0.3, 0.7], [1, 1])
    """
    fcol = _as_moment_vector(f, 'moment_joint_marking', 'f')
    pv = np.asarray(p, dtype=float).ravel()
    dv = np.asarray(dims).ravel().astype(int)
    if pv.size != dv.size:
        raise ValueError('moment_joint_marking: p and dims must have the same length.')
    if dv.min(initial=0) < 0:
        raise ValueError('moment_joint_marking: The maximum orders must be nonnegative.')
    if int(dv.sum()) > fcol.size - 1:
        raise ValueError('moment_joint_marking: The aggregate factorial moments must '
                         'reach order sum(dims) = %d.' % int(dv.sum()))
    out = np.zeros(tuple(dv + 1))
    for a in np.ndindex(*out.shape):
        coef = 1.0
        for j, aj in enumerate(a):
            coef *= pv[j] ** aj
        out[a] = coef * fcol[sum(a)]
    return out


def moment_joint_aggregate(F):
    """
    Factorial moments of the total count from the joint factorial moments of
    its parts.

    For N = N_1+...+N_d the Vandermonde convolution of falling factorials gives

        f_n = sum_(|a|=n) (n! / prod_j a_j!) * F_a

    which holds for ANY joint law of the parts, marked or not, and is the
    inverse direction of moment_joint_marking whenever the marking is
    multinomial. The order reached is limited by the smallest per-class order
    in F, since the term a = n*e_j must be available for every j.

    Args:
        F: Array-like of shape (n_1+1,...,n_d+1) holding the joint factorial
           moments of the parts.

    Returns:
        1-D numpy array of length min_j(n_j)+1 holding f_0,...,f_min_j(n_j),
        the factorial moments of the total.

    Raises:
        ValueError: If F is empty.

    Example:
        f = moment_joint_aggregate(moment_joint_marking([1, 2, 4], [0.3, 0.7], [1, 1]))
    """
    arr = _as_moment_tensor(F, 'moment_joint_aggregate', 'F')
    nmax = min(arr.shape) - 1
    out = np.zeros(nmax + 1)
    for n in range(nmax + 1):
        acc = 0.0
        for a in np.ndindex(*arr.shape):
            if sum(a) != n:
                continue
            coef = float(factorial(n))
            for aj in a:
                coef /= factorial(aj)
            acc += coef * arr[a]
        out[n] = acc
    return out


# ---------------------------------------------------------------------------
# The survival (tail) vertex
# ---------------------------------------------------------------------------

def moment_binomial_from_tail(t):
    """
    Convert survival (tail) probabilities into binomial moments.

    For a nonnegative integer random variable N with survival sequence
    t_m = P(N >= m),

        b_j = E[nchoosek(N,j)] = sum_{m>=j} nchoosek(m-1,j-1) * t_m,  j >= 1

    and b_0 = t_0 = 1. Unlike every other edge of the house this transform is
    UPPER triangular, so it consumes the whole tail: the result is exact only
    if the sequence covers the support, i.e. t_m = 0 beyond the last element
    supplied. This is the natural entry point for a closed queueing network,
    whose queue lengths are bounded by the population, and where the joint
    survival probabilities are ratios of normalizing constants.

    Truncating the tail early yields a strict LOWER bound on every b_j, since
    all the coefficients and all the tail values are nonnegative; the bound is
    monotone in the truncation level. The bound is not inherited by the central
    moments downstream, whose conversion alternates in sign.

    Args:
        t: Array-like of length n+1 holding t_0,...,t_n, i.e. element m is
           P(N >= m) and element 0 is t_0 = 1.

    Returns:
        1-D numpy array of length n+1 holding b_0,...,b_n.

    Raises:
        ValueError: If t is not a nonempty vector.

    Example:
        b = moment_binomial_from_tail([1, 1, 1, 1, 0])   # N = 3 with prob 1
    """
    tcol = _as_moment_vector(t, 'moment_binomial_from_tail', 't')
    n = tcol.size - 1
    return moment_housematrix('binomial_from_tail', n) @ tcol


def moment_tail_from_binomial(b):
    """
    Convert binomial moments into survival (tail) probabilities.

    Inverts moment_binomial_from_tail,

        t_m = sum_{j>=m} (-1)^(j-m) * nchoosek(j-1,m-1) * b_j,  m >= 1

    with t_0 = 1. The inversion is exact on the finite box supplied, the matrix
    being unit upper triangular, but it reconstructs the true tail only if the
    binomial moments were themselves those of a law supported on 0,...,n.

    Args:
        b: Array-like of length n+1 holding b_0,...,b_n.

    Returns:
        1-D numpy array of length n+1 holding t_0,...,t_n.

    Raises:
        ValueError: If b is not a nonempty vector.

    Example:
        t = moment_tail_from_binomial(moment_binomial_from_tail([1, 1, 1, 0]))
    """
    bcol = _as_moment_vector(b, 'moment_tail_from_binomial', 'b')
    n = bcol.size - 1
    return moment_housematrix('tail_from_binomial', n) @ bcol


def moment_joint_binomial_from_tail(t):
    """
    Convert joint survival probabilities into joint binomial moments.

    For a nonnegative integer random vector (N_1,...,N_d) with joint survival
    array t_(m_1,...,m_d) = P(N_1 >= m_1, ..., N_d >= m_d),

        b_(k) = E[prod_j nchoosek(N_j,k_j)]
              = sum_(m>=k) prod_j nchoosek(m_j-1,k_j-1) * t_(m)

    The transform is again the tensor product of the univariate one, which is
    what makes the mode-by-mode evaluation legitimate: an entry of the array
    with k_j = 0 selects m_j = 0, and t_(0,m_2,...) is by construction the
    marginal survival array of the remaining coordinates.

    As in the univariate case the transform is upper triangular, so the array
    must cover the joint support to be exact; truncating gives lower bounds.

    Args:
        t: Array-like of shape (n_1+1,...,n_d+1) of joint survival
           probabilities, element 0 being 1.

    Returns:
        ndarray of the same shape holding the joint binomial moments.

    Raises:
        ValueError: If t is empty.

    Example:
        b = moment_joint_binomial_from_tail(t)
    """
    arr = _as_moment_tensor(t, 'moment_joint_binomial_from_tail', 't')
    return moment_jointtrans(arr, 'binomial_from_tail')


def moment_joint_tail_from_binomial(b):
    """
    Convert joint binomial moments into joint survival probabilities.

    Inverse of moment_joint_binomial_from_tail.

    Args:
        b: Array-like of shape (n_1+1,...,n_d+1) of joint binomial moments.

    Returns:
        ndarray of the same shape holding the joint survival probabilities.

    Raises:
        ValueError: If b is empty.
    """
    arr = _as_moment_tensor(b, 'moment_joint_tail_from_binomial', 'b')
    return moment_jointtrans(arr, 'tail_from_binomial')


def moment_joint_central_from_tail(t):
    """
    Joint central moments of a nonnegative integer random vector from its joint
    survival array.

    Composes the four edges that separate the two vertices, tail ->
    binomial -> factorial -> raw -> central, reading the means off the raw
    array. This is the whole path from a solver that produces survival
    probabilities (a closed queueing network through its normalizing constants,
    a CTMC through its stationary distribution, a simulator through a
    histogram) to the covariances and higher central moments.

    Args:
        t: Array-like of shape (n_1+1,...,n_d+1) of joint survival
           probabilities covering the support, element 0 being 1.

    Returns:
        ndarray of the same shape holding the joint central moments. The entry
        of multi-order e_j+e_l is the covariance of N_j and N_l.

    Raises:
        ValueError: If t is empty or has a dimension of extent 1.

    Example:
        mc = moment_joint_central_from_tail(t)
    """
    arr = _as_moment_tensor(t, 'moment_joint_central_from_tail', 't')
    b = moment_joint_binomial_from_tail(arr)
    f = moment_joint_factorial_from_binomial(b)
    return moment_joint_central_from_raw(moment_joint_raw_from_factorial(f))
