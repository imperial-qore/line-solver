"""
Markov chain and stochastic process utilities.

Low-level functions for:
- Generator matrix validation
- Transition probability checking
- Markov chain solving
- Continuous/discrete renewal processes
"""

import numpy as np
from line_solver import native_to_array


def lib_markov_check_generator(Q):
    """
    Validate continuous-time Markov chain generator matrix.

    Requirements:
    - Diagonal elements non-positive
    - Off-diagonal elements non-negative
    - Row sums equal to zero

    Args:
        Q (array-like): Generator matrix

    Returns:
        bool: True if valid generator matrix
    """
    Q = np.asarray(Q)
    Q_arr = native_to_array(Q)

    try:
            Q_arr
        )
        return bool(result)
    except Exception as e:
        return False


def lib_markov_check_probability_matrix(P):
    """
    Validate discrete-time transition probability matrix.

    Requirements:
    - Elements in [0, 1]
    - Row sums equal to 1

    Args:
        P (array-like): Transition probability matrix

    Returns:
        bool: True if valid probability matrix
    """
    P = np.asarray(P)
    P_arr = native_to_array(P)

    try:
            P_arr
        )
        return bool(result)
    except Exception as e:
        return False


def lib_markov_check_probability_vector(pi):
    """
    Validate probability vector.

    Requirements:
    - All elements in [0, 1]
    - Sum equals 1

    Args:
        pi (array-like): Probability vector

    Returns:
        bool: True if valid probability vector
    """
    pi = np.asarray(pi).flatten()
    pi_arr = native_to_array(pi)

    try:
            pi_arr
        )
        return bool(result)
    except Exception as e:
        return False


def lib_markov_ctmc_solve(Q):
    """
    Solve CTMC for stationary distribution.

    Args:
        Q (array-like): Generator matrix

    Returns:
        ndarray: Stationary probability distribution
    """
    Q = np.asarray(Q)
    Q_arr = native_to_array(Q)

    try:
            Q_arr
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"CTMC solving failed: {str(e)}")


def lib_markov_dtmc_solve(P):
    """
    Solve DTMC for stationary distribution.

    Args:
        P (array-like): Transition probability matrix

    Returns:
        ndarray: Stationary probability distribution
    """
    P = np.asarray(P)
    P_arr = native_to_array(P)

    try:
            P_arr
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"DTMC solving failed: {str(e)}")


def lib_markov_crp_solve(D0, D1):
    """
    Solve continuous renewal process.

    Args:
        D0 (array-like): Generator matrix
        D1 (array-like): Completion matrix

    Returns:
        ndarray: Solution vector
    """
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    D0_arr = native_to_array(D0)
    D1_arr = native_to_array(D1)

    try:
            D0_arr, D1_arr
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"CRP solving failed: {str(e)}")


def lib_markov_drp_solve(D0, D1):
    """
    Solve discrete renewal process.

    Args:
        D0 (array-like): Transition matrix
        D1 (array-like): Completion matrix

    Returns:
        ndarray: Solution vector
    """
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    D0_arr = native_to_array(D0)
    D1_arr = native_to_array(D1)

    try:
            D0_arr, D1_arr
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"DRP solving failed: {str(e)}")


def lib_markov_validate_chain(Q_or_P):
    """
    Validate any Markov chain (continuous or discrete).

    Args:
        Q_or_P (array-like): Generator or transition matrix

    Returns:
        dict: Validation results with details
    """
    matrix = np.asarray(Q_or_P)

    # Check if generator (has negative diagonal and zero row sums)
    is_generator = (
        np.allclose(matrix.sum(axis=1), 0) and
        np.all(np.diag(matrix) <= 0)
    )

    # Check if probability matrix (row sums = 1, all elements in [0,1])
    is_probability = (
        np.allclose(matrix.sum(axis=1), 1) and
        np.all(matrix >= 0) and
        np.all(matrix <= 1)
    )

    return {
        'is_generator': is_generator,
        'is_probability_matrix': is_probability,
        'type': 'CTMC' if is_generator else ('DTMC' if is_probability else 'Invalid'),
        'valid': is_generator or is_probability
    }
