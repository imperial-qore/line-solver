"""
Acyclic Phase-Type (APH) distribution functions for KPC-Toolbox.

Native Python implementations of APH distribution analysis and fitting.
"""

import numpy as np
from enum import Enum
from typing import Tuple, List, Dict
from .mc import dtmc_solve


class ConvolutionPattern(Enum):
    """Convolution patterns for APH simplification."""
    SEQUENCE = 1  # Sequential structure
    PARALLEL = 2  # Parallel structure
    BRANCH = 3    # Branch structure


def aph_simplify(
    a1: np.ndarray,
    T1: np.ndarray,
    a2: np.ndarray,
    T2: np.ndarray,
    p1: float = 1.0,
    p2: float = 1.0,
    pattern: ConvolutionPattern = ConvolutionPattern.SEQUENCE
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simplify/combine two APH distributions using a specified pattern.

    Args:
        a1: Initial probability vector of first distribution
        T1: Rate matrix of first distribution
        a2: Initial probability vector of second distribution
        T2: Rate matrix of second distribution
        p1: Branch probability for first distribution (for BRANCH pattern)
        p2: Branch probability for second distribution (for BRANCH pattern)
        pattern: Convolution pattern to use

    Returns:
        Tuple of (combined alpha, combined T)
    """
    a1 = np.asarray(a1, dtype=np.float64).flatten()
    T1 = np.asarray(T1, dtype=np.float64)
    a2 = np.asarray(a2, dtype=np.float64).flatten()
    T2 = np.asarray(T2, dtype=np.float64)

    order1 = len(a1)
    order2 = len(a2)

    if pattern == ConvolutionPattern.SEQUENCE:
        # Sequential structure: first complete first distribution, then second
        n = order1 + order2

        # Compute sum of a1
        a1_sum = np.sum(a1)

        # alpha = [a1, (1-sum(a1))*a2]
        alpha = np.zeros(n)
        alpha[:order1] = a1
        alpha[order1:] = (1.0 - a1_sum) * a2

        # T = [T1, (-T1*e1)*a2; 0, T2]
        T = np.zeros((n, n))

        # Copy T1
        T[:order1, :order1] = T1

        # Compute exit rates from T1
        exit_rates = -np.sum(T1, axis=1)

        # (-T1*e1)*a2 block
        for i in range(order1):
            for j in range(order2):
                T[i, order1 + j] = exit_rates[i] * a2[j]

        # Copy T2
        T[order1:, order1:] = T2

        return alpha, T

    elif pattern == ConvolutionPattern.PARALLEL:
        # Parallel structure: minimum of two distributions
        n = order1 * order2 + order1 + order2

        # Compute sums
        a1_sum = np.sum(a1)
        a2_sum = np.sum(a2)

        # alpha = [kron(a1,a2), (1-sum(a2))*a1, (1-sum(a1))*a2]
        alpha = np.zeros(n)

        # Kronecker product part
        idx = 0
        for i in range(order1):
            for j in range(order2):
                alpha[idx] = a1[i] * a2[j]
                idx += 1

        # (1-sum(a2))*a1
        for i in range(order1):
            alpha[idx] = (1.0 - a2_sum) * a1[i]
            idx += 1

        # (1-sum(a1))*a2
        for i in range(order2):
            alpha[idx] = (1.0 - a1_sum) * a2[i]
            idx += 1

        # Build T matrix
        T = np.zeros((n, n))

        # Compute exit rates
        exit1 = -np.sum(T1, axis=1)
        exit2 = -np.sum(T2, axis=1)

        # First block: kron(T1, I) + kron(I, T2)
        for i in range(order1):
            for j in range(order2):
                row_idx = i * order2 + j
                for k in range(order1):
                    for l in range(order2):
                        col_idx = k * order2 + l
                        value = 0.0
                        if j == l:
                            value += T1[i, k]
                        if i == k:
                            value += T2[j, l]
                        T[row_idx, col_idx] = value

        # Transition to second block
        for i in range(order1):
            for j in range(order2):
                row_idx = i * order2 + j
                col_idx = order1 * order2 + i
                T[row_idx, col_idx] = exit2[j]

        # Transition to third block
        for i in range(order1):
            for j in range(order2):
                row_idx = i * order2 + j
                col_idx = order1 * order2 + order1 + j
                T[row_idx, col_idx] = exit1[i]

        # Second block: T1
        T[order1 * order2:order1 * order2 + order1,
          order1 * order2:order1 * order2 + order1] = T1

        # Third block: T2
        T[order1 * order2 + order1:,
          order1 * order2 + order1:] = T2

        return alpha, T

    else:  # BRANCH
        # Branch structure: probabilistic choice between distributions
        n = order1 + order2

        # alpha = [p1*a1, p2*a2]
        alpha = np.zeros(n)
        alpha[:order1] = p1 * a1
        alpha[order1:] = p2 * a2

        # T = [T1, 0; 0, T2]
        T = np.zeros((n, n))
        T[:order1, :order1] = T1
        T[order1:, order1:] = T2

        return alpha, T


def aph_convpara(
    distributions: List[Tuple[np.ndarray, np.ndarray]]
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform convolution on parallel structure with any number of elements.

    Args:
        distributions: List of (alpha, T) pairs for each distribution

    Returns:
        Combined (alpha, T)
    """
    if not distributions:
        raise ValueError("Need at least one distribution")

    if len(distributions) == 1:
        return distributions[0]

    alpha, T = aph_simplify(
        distributions[0][0], distributions[0][1],
        distributions[1][0], distributions[1][1],
        pattern=ConvolutionPattern.PARALLEL
    )

    for i in range(2, len(distributions)):
        alpha, T = aph_simplify(
            alpha, T,
            distributions[i][0], distributions[i][1],
            pattern=ConvolutionPattern.SEQUENCE
        )

    return alpha, T


def aph_convseq(
    distributions: List[Tuple[np.ndarray, np.ndarray]]
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform convolution on sequential structure with any number of elements.

    Args:
        distributions: List of (alpha, T) pairs for each distribution

    Returns:
        Combined (alpha, T)
    """
    if not distributions:
        raise ValueError("Need at least one distribution")

    if len(distributions) == 1:
        return distributions[0]

    alpha, T = aph_simplify(
        distributions[0][0], distributions[0][1],
        distributions[1][0], distributions[1][1],
        pattern=ConvolutionPattern.SEQUENCE
    )

    for i in range(2, len(distributions)):
        alpha, T = aph_simplify(
            alpha, T,
            distributions[i][0], distributions[i][1],
            pattern=ConvolutionPattern.SEQUENCE
        )

    return alpha, T


def aph_rand(K: int = 2, seed: int = None) -> Dict[str, np.ndarray]:
    """
    Generate a random APH (acyclic phase-type) distribution as a MAP.

    Args:
        K: Order of the APH distribution (default: 2)
        seed: Random seed (optional)

    Returns:
        MAP representation {'D0': D0, 'D1': D1}
    """
    if seed is not None:
        np.random.seed(seed)

    D0 = np.zeros((K, K))
    D1 = np.zeros((K, K))

    # Generate random upper triangular D0 (acyclic)
    for i in range(K):
        for j in range(i, K):
            D0[i, j] = np.random.rand()

    # Generate random D1
    D1 = np.random.rand(K, K)

    # Normalize to make valid MAP
    MAP = _map_normalize(_map_renewal({'D0': D0, 'D1': D1}))

    return MAP


def _map_renewal(MAP: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Make MAP a renewal process: MAPOUT{2} = MAPIN{2}*ones(n,1)*map_pie(MAPIN)."""
    D0 = MAP['D0']
    D1 = MAP['D1']
    n = D0.shape[0]

    # Compute map_pie: steady-state of embedded DTMC at departure instants
    # map_prob: stationary distribution of underlying CTMC Q = D0 + D1
    Q = D0 + D1
    try:
        from .mc import ctmc_solve as _ctmc_solve
        pi_q = _ctmc_solve(Q)
    except Exception:
        pi_q = np.ones(n) / n

    # map_pie: A = pi_q * D1, PIE = A / (A * ones(n,1))
    if pi_q.ndim == 1:
        pi_q = pi_q.reshape(1, -1)
    A = pi_q @ D1  # row vector
    A_sum = np.sum(A)
    if A_sum > 0:
        pie = A / A_sum
    else:
        pie = np.ones((1, n)) / n

    # D1_new = D1 * ones(n,1) * pie
    ones_col = np.ones((n, 1))
    D1_ones = D1 @ ones_col  # column vector (n,1)
    D1_new = D1_ones @ pie  # (n,1) @ (1,n) = (n,n)

    return {'D0': D0.copy(), 'D1': D1_new}


def _map_normalize(MAP: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Normalize MAP to be valid.

    Following MATLAB map_normalize:
    1. Take real part of all entries
    2. Zero out all negative entries in both D0 and D1
    3. Recompute D0 diagonal so each row of D0+D1 sums to zero
    """
    D0 = MAP['D0'].copy()
    D1 = MAP['D1'].copy()
    n = D0.shape[0]

    # Take real parts (handles complex values from sqrt of negative numbers)
    D0 = np.real(D0)
    D1 = np.real(D1)

    # Zero out negative entries in both matrices
    D0[D0 < 0] = 0
    D1[D1 < 0] = 0

    # Recompute D0 diagonal so each row of D0+D1 sums to zero
    for i in range(n):
        D0[i, i] = 0
        row_sum = np.sum(D0[i, :]) + np.sum(D1[i, :])
        D0[i, i] = -row_sum

    return {'D0': D0, 'D1': D1}


def aph_fit(e1: float, e2: float, e3: float, nmax: int = 10
            ) -> Tuple[Dict[str, np.ndarray], bool]:
    """
    Fit an APH distribution to match first three moments.

    Based on: A.Bobbio, A.Horvath, M.Telek, "Matching three moments
    with minimal acyclic phase type distributions", Stochastic Models
    21:303-326, 2005.

    Args:
        e1: First moment (mean)
        e2: Second moment
        e3: Third moment
        nmax: Maximum order to try (default: 10)

    Returns:
        Tuple of (fitted MAP, isExact flag)
    """
    is_exact = True

    if np.isinf(e2) or np.isinf(e3):
        # Return exponential
        D0 = np.array([[-1.0 / e1]])
        D1 = np.array([[1.0 / e1]])
        return {'D0': D0, 'D1': D1}, True

    n2 = e2 / (e1 * e1)
    n3 = e3 / (e1 * e2)

    # Find suitable order
    n2_feas = False
    n3_ub_feas = False
    n3_lb_feas = False
    n = 1
    un = 0.0
    un_1 = 0.0

    while (not n2_feas or not n3_lb_feas or not n3_ub_feas) and n < nmax:
        n += 1
        pn = ((n + 1) * (n2 - 2) / (3 * n2 * (n - 1))) * \
             (-2 * np.sqrt(n + 1.0) / np.sqrt(4.0 * (n + 1) - 3 * n * n2) - 1)
        an = (n2 - 2) / (pn * (1 - n2) + np.sqrt(pn * pn + pn * n * (n2 - 2) / (n - 1)))
        ln = ((3 + an) * (n - 1) + 2 * an) / ((n - 1) * (1 + an * pn)) - \
             (2 * an * (n + 1)) / (2 * (n - 1) + an * pn * (n * an + 2 * n - 2))
        un_1 = un
        un = (1.0 / (n * n * n2)) * (2 * (n - 2) * (n * n2 - n - 1) *
             np.sqrt(1 + n * (n2 - 2) / (n - 1)) + (n + 2) * (3 * n * n2 - 2 * n - 2))

        if n2 >= (n + 1.0) / n and n2 <= (n + 4.0) / (n + 1):
            n2_feas = True
            if n3 >= ln:
                n3_lb_feas = True
        elif n2 >= (n + 4.0) / (n + 1):
            n2_feas = True
            if n3 >= n2 * (n + 1) / n:
                n3_lb_feas = True

        if n2 >= (n + 1.0) / n and n2 <= n / (n - 1):
            n2_feas = True
            if n3 <= un:
                n3_ub_feas = True
        elif n2 >= n / (n - 1):
            n2_feas = True
            n3_ub_feas = True

    fit_n2 = n2
    fit_n3 = n3

    if not n2_feas or not n3_lb_feas or not n3_ub_feas or n >= nmax:
        # Cannot match exactly, use feasible approximation
        fit_n2 = (n + 1.0) / n
        fit_n3 = 2 * fit_n2 - 1
        is_exact = False

    # Fitting algorithm
    if fit_n2 <= n / (n - 1) or fit_n3 <= 2 * fit_n2 - 1:
        # Case 1 of 2
        b = 2 * (4 - n * (3 * fit_n2 - 4)) / (fit_n2 * (4 + n - n * fit_n3) +
                np.sqrt(n * fit_n2) * np.sqrt(12 * fit_n2 * fit_n2 * (n + 1) +
                16 * fit_n3 * (n + 1) + fit_n2 * (n * (fit_n3 - 15) * (fit_n3 + 1) - 8 * (fit_n3 + 3))))
        a = (b * fit_n2 - 2) * (n - 1) * b / ((b - 1) * n)
        p = (b - 1) / a
        lambda_val = 1.0
        mu = lambda_val * (n - 1) / a

        alpha = np.zeros(n)
        alpha[0] = p
        alpha[n - 1] = 1 - p

        T = np.zeros((n, n))
        for i in range(n):
            T[i, i] = -mu
            if i < n - 1:
                T[i, i + 1] = mu
        T[n - 1, n - 1] = -lambda_val

        # Build D0 = T and D1 = -T*e*alpha
        D0 = T.copy()
        exit_rates = -np.sum(T, axis=1)
        D1 = np.outer(exit_rates, alpha)

    elif fit_n2 > n / (n - 1) and fit_n3 > un_1:
        # Case 2 of 2
        K1 = n - 1.0
        K2 = n - 2.0
        K3 = 3 * fit_n2 - 2 * fit_n3
        K4 = fit_n3 - 3
        K5 = n - fit_n2
        K6 = 1 + fit_n2 - fit_n3
        K7 = n + fit_n2 - n * fit_n2
        K8 = 3 + 3 * fit_n2**2 + fit_n3 - 3 * fit_n2 * fit_n3

        inner_sqrt = -16 * K1**2 * K7**6 + \
            (4 * K1 * K5**3 + K1**2 * K2 * K4**2 * n * fit_n2**2 +
             4 * K2 * n * fit_n2 * (K4 * n**2 - 3 * K6 * fit_n2 + K8 * n))**2
        K9 = 108 * K1**2 * (4 * K2**2 * K3 * n**2 * fit_n2 +
             K1**2 * K2 * K4**2 * n * fit_n2**2 +
             4 * K1 * K5 * (K5**2 - 3 * K2 * K6 * n * fit_n2) +
             np.sqrt(complex(inner_sqrt, 0)).real)

        K10 = K4**2 / (4 * K3**2) - K5 / (K1 * K3 * fit_n2)
        # Use signed cube root for K9^(1/3) to handle negative values
        K9_cbrt = np.abs(K9)**(1.0/3.0) * np.sign(K9) if K9 != 0 else 0.0
        K11 = 2**(1.0/3.0) * (3 * K5**2 + K2 * (K3 + 2 * K4) * n * fit_n2) / \
              (K3 * K9_cbrt * fit_n2) if K9_cbrt != 0 else 0.0
        K12 = K9_cbrt / (3 * 2**(7.0/3.0) * K1**2 * K3 * fit_n2)
        K13 = np.sqrt(K10 + K11 + K12)
        K14 = (6 * K1 * K3 * K4 * K5 + 4 * K2 * K3**2 * n - K1**2 * K4**3 * fit_n2) / \
              (4 * K1**2 * K3**3 * K13 * fit_n2)
        K15 = -K4 / (2 * K3)
        K16 = np.sqrt(2 * K10 - K11 - K12 - K14)
        K17 = np.sqrt(2 * K10 - K11 - K12 + K14)

        inner_sqrt18 = 81 * (4 * K5**3 + 4 * K2 * K4 * K5 * n * fit_n2 +
                       K1 * K2 * K4**2 * n * fit_n2**2)**2 - \
                       48 * (3 * K5**2 + 2 * K2 * K4 * n * fit_n2)**3
        K18 = 36 * K5**3 + 36 * K2 * K4 * K5 * n * fit_n2 + \
              9 * K1 * K2 * K4**2 * n * fit_n2**2 - np.sqrt(complex(inner_sqrt18, 0)).real
        K18_cbrt = np.abs(K18)**(1.0/3.0) * np.sign(K18) if K18 != 0 else 0.0
        K19 = -K5 / (K1 * K4 * fit_n2) - \
              2**(2.0/3.0) * (3 * K5**2 + 2 * K2 * K4 * n * fit_n2) / \
              (3**(1.0/3.0) * K1 * K4 * fit_n2 * K18_cbrt) - \
              K18_cbrt / (6**(2.0/3.0) * K1 * K4 * fit_n2) if K18_cbrt != 0 else 0.0
        K20 = 6 * K1 * K3 * K4 * K5 + 4 * K2 * K3**2 * n - K1**2 * K4**3 * fit_n2
        K21 = K11 + K12 + K5 / (2 * n * K1 * K3)
        K22 = np.sqrt(3 * K4**2 / (4 * K3**2) - 3 * K5 / (K1 * K3 * fit_n2) +
              np.sqrt(4 * K21**2 - n * K2 / (fit_n2 * K1**2 * K3)))

        if fit_n3 > un_1 and fit_n3 < 3 * fit_n2 / 2:
            f = K13 + K15 - K17
        elif fit_n3 == 2 * fit_n2 / 2:
            f = K19
        elif fit_n3 > 3 * fit_n2 / 2 and K20 > 0:
            f = -K13 + K15 + K16
        elif K20 == 0:
            f = K15 + K22
        else:
            # K20 < 0
            f = K13 + K15 + K17

        a = 2 * (f - 1) * (n - 1) / ((n - 1) * (fit_n2 * f**2 - 2 * f + 2) - n)
        p = (f - 1) * a
        lambda_val = 1.0
        mu = lambda_val * (n - 1) / a

        alpha = np.zeros(n)
        alpha[0] = p
        alpha[1] = 1 - p

        T = np.zeros((n, n))
        for i in range(n):
            T[i, i] = -mu
            if i < n - 1:
                T[i, i + 1] = mu
        T[0, 0] = -lambda_val
        T[0, 1] = lambda_val

        # Build D0 = T and D1 = -T*e*alpha
        D0 = T.copy()
        exit_rates = -np.sum(T, axis=1)
        D1 = np.outer(exit_rates, alpha)

    else:
        # Moment set cannot be matched with an APH distribution
        import warnings
        warnings.warn('moment set cannot be matched with an APH distribution')
        is_exact = False
        # Fall back to Erlang approximation
        mu = n / e1

        alpha = np.zeros(n)
        alpha[0] = 1.0

        T = np.zeros((n, n))
        for i in range(n):
            T[i, i] = -mu
            if i < n - 1:
                T[i, i + 1] = mu

        D0 = T.copy()
        exit_rates = -np.sum(T, axis=1)
        D1 = np.outer(exit_rates, alpha)

    # Normalize then scale to match mean (matches MATLAB: map_scale(map_normalize({T,D1}), e1))
    MAP = _map_scale(_map_normalize({'D0': D0, 'D1': D1}), e1)

    return MAP, is_exact


def _map_scale(MAP: Dict[str, np.ndarray], target_mean: float
               ) -> Dict[str, np.ndarray]:
    """Scale MAP to achieve target mean: ratio = map_mean(MAP) / target_mean."""
    D0 = MAP['D0'].copy()
    D1 = MAP['D1'].copy()
    n = D0.shape[0]

    # Compute current mean = 1 / map_lambda
    # map_lambda = pi_q * D1 * ones(n,1) where pi_q = ctmc_solve(D0+D1)
    Q = D0 + D1
    try:
        from .mc import ctmc_solve as _ctmc_solve
        pi_q = _ctmc_solve(Q)
    except Exception:
        pi_q = np.ones(n) / n

    if pi_q.ndim == 1:
        pi_q = pi_q.reshape(1, -1)
    ones_col = np.ones((n, 1))
    lam = float(pi_q @ D1 @ ones_col)

    if lam > 0:
        current_mean = 1.0 / lam
    else:
        current_mean = 1.0

    ratio = current_mean / target_mean

    # Scale and normalize
    D0 *= ratio
    D1 *= ratio
    result = _map_normalize({'D0': D0, 'D1': D1})
    return result


def ph2hyper(PH: Dict[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert a hyper-exponential PH distribution to its rate/probability form.

    Args:
        PH: Phase-type distribution as {'D0': D0, 'D1': D1}

    Returns:
        Tuple of (rates, probabilities)
    """
    D0 = PH['D0']
    n = D0.shape[0]

    # Check if diagonal (hyper-exponential)
    for i in range(n):
        for j in range(n):
            if i != j and abs(D0[i, j]) > 1e-10:
                raise ValueError("The PH distribution is not hyper-exponential")

    # Extract rates (negative diagonal of D0)
    rates = -np.diag(D0)

    # Compute probabilities from stationary distribution of embedded DTMC
    D1 = PH['D1']
    P = np.zeros((n, n))
    for i in range(n):
        if rates[i] > 0:
            P[i, :] = D1[i, :] / rates[i]

    probs = dtmc_solve(P)

    return rates, probs


def hyper_rand(rates: np.ndarray, probs: np.ndarray, n_samples: int,
               seed: int = None) -> np.ndarray:
    """
    Generate random samples from a hyper-exponential distribution.

    Args:
        rates: Exponential rates
        probs: Selection probabilities
        n_samples: Number of samples
        seed: Random seed (optional)

    Returns:
        Array of samples
    """
    rates = np.asarray(rates, dtype=np.float64)
    probs = np.asarray(probs, dtype=np.float64)

    if seed is not None:
        np.random.seed(seed)

    # Normalize probabilities
    probs = probs / np.sum(probs)

    samples = np.zeros(n_samples)
    cum_probs = np.cumsum(probs)

    for s in range(n_samples):
        # Select which exponential to use
        u = np.random.rand()
        selected = np.searchsorted(cum_probs, u)
        selected = min(selected, len(rates) - 1)

        # Generate exponential sample
        samples[s] = -np.log(np.random.rand()) / rates[selected]

    return samples


__all__ = [
    'aph_simplify',
    'aph_convpara',
    'aph_convseq',
    'aph_rand',
    'aph_fit',
    'ph2hyper',
    'hyper_rand',
    'ConvolutionPattern',
]
