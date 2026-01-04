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
    """Make MAP a renewal process (D1 = -D0*e*alpha)."""
    D0 = MAP['D0']
    D1 = MAP['D1']
    n = D0.shape[0]

    # Compute stationary distribution
    P = D1.copy()
    row_sums = np.sum(D1, axis=1)
    for i in range(n):
        if row_sums[i] > 0:
            P[i, :] /= row_sums[i]

    # Simple alpha approximation
    alpha = np.ones(n) / n

    # Exit rates
    exit_rates = -np.sum(D0, axis=1)

    # D1 = exit_rates * alpha
    D1_new = np.outer(exit_rates, alpha)

    return {'D0': D0, 'D1': D1_new}


def _map_normalize(MAP: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Normalize MAP to be valid."""
    D0 = MAP['D0'].copy()
    D1 = MAP['D1'].copy()
    n = D0.shape[0]

    # Ensure D0 has negative diagonal
    for i in range(n):
        off_diag_sum = np.sum(D0[i, :]) - D0[i, i] + np.sum(D1[i, :])
        D0[i, i] = -off_diag_sum

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

    while (not n2_feas or not n3_lb_feas or not n3_ub_feas) and n < nmax:
        n += 1

        # Check n2 feasibility conditions
        if n2 >= (n + 1.0) / n and n2 <= (n + 4.0) / (n + 1):
            n2_feas = True
            # Compute lower bound on n3
            try:
                with np.errstate(divide='ignore', invalid='ignore'):
                    pn = ((n + 1) * (n2 - 2) / (3 * n2 * (n - 1))) * \
                         (-2 * np.sqrt(n + 1.0) / np.sqrt(4.0 * (n + 1) - 3 * n * n2) - 1)
                    an = (n2 - 2) / (pn * (1 - n2) + np.sqrt(pn * pn + pn * n * (n2 - 2) / (n - 1)))
                    ln = ((3 + an) * (n - 1) + 2 * an) / ((n - 1) * (1 + an * pn)) - \
                         (2 * an * (n + 1)) / (2 * (n - 1) + an * pn * (n * an + 2 * n - 2))
                if np.isfinite(ln) and n3 >= ln:
                    n3_lb_feas = True
            except (ValueError, ZeroDivisionError):
                pass
        elif n2 >= (n + 4.0) / (n + 1):
            n2_feas = True
            if n3 >= n2 * (n + 1) / n:
                n3_lb_feas = True

        # Check n3 upper bound
        if n2 >= (n + 1.0) / n and n2 <= n / (n - 1):
            n2_feas = True
            try:
                with np.errstate(divide='ignore', invalid='ignore'):
                    un = (1.0 / (n * n * n2)) * (2 * (n - 2) * (n * n2 - n - 1) *
                         np.sqrt(1 + n * (n2 - 2) / (n - 1)) + (n + 2) * (3 * n * n2 - 2 * n - 2))
                if np.isfinite(un) and n3 <= un:
                    n3_ub_feas = True
            except (ValueError, ZeroDivisionError):
                pass
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
        # Case 1
        try:
            denom = fit_n2 * (4 + n - n * fit_n3) + \
                    np.sqrt(n * fit_n2) * np.sqrt(max(0, 12 * fit_n2 * fit_n2 * (n + 1) +
                    16 * fit_n3 * (n + 1) + fit_n2 * (n * (fit_n3 - 15) * (fit_n3 + 1) - 8 * (fit_n3 + 3))))

            if abs(denom) > 1e-10:
                b = 2 * (4 - n * (3 * fit_n2 - 4)) / denom
            else:
                b = 1.0

            a = (b * fit_n2 - 2) * (n - 1) * b / max(1e-10, (b - 1) * n)
            p = (b - 1) / max(1e-10, a)
            lambda_val = 1.0
            mu = lambda_val * (n - 1) / max(1e-10, a)
        except (ValueError, ZeroDivisionError):
            mu = n / e1
            p = 0.5
            lambda_val = 1.0

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

    else:
        # Case 2 - simplified implementation
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
        is_exact = False

    # Scale to match mean
    MAP = _map_scale({'D0': D0, 'D1': D1}, e1)

    return MAP, is_exact


def _map_scale(MAP: Dict[str, np.ndarray], target_mean: float
               ) -> Dict[str, np.ndarray]:
    """Scale MAP to achieve target mean."""
    D0 = MAP['D0']
    D1 = MAP['D1']
    n = D0.shape[0]

    # Compute current mean
    Q = D0 + D1
    # Stationary distribution of embedded DTMC
    P = np.zeros((n, n))
    for i in range(n):
        row_sum = np.sum(D1[i, :])
        if row_sum > 0:
            P[i, :] = D1[i, :] / row_sum

    # Simple mean calculation
    try:
        pi = dtmc_solve(P)
        current_mean = np.sum(pi / np.maximum(-np.diag(D0), 1e-10))
    except Exception:
        current_mean = 1.0

    if current_mean > 0:
        scale = target_mean / current_mean
    else:
        scale = 1.0

    return {'D0': D0 / scale, 'D1': D1 / scale}


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
