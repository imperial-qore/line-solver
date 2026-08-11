"""
Stationary Moment Computation (SMC) utilities.

Low-level solvers for computing stationary moments in queues:
- QBD (Quasi-Birth-Death) processes
- GIM1 (G/M/1) queue analysis
- MG1 (M/G/1) queue analysis
- Fundamental matrix computation
"""

import numpy as np
from line_solver import native_to_array


def lib_qbd_fundamental_matrix(D0, D1, D2, method='lr'):
    """
    Compute fundamental matrix R for QBD process.

    Args:
        D0 (array-like): Downward transition matrix
        D1 (array-like): Transition within level
        D2 (array-like): Upward transition matrix
        method (str): Solution method - 'lr', 'fi', 'is', 'cr', 'eg', 'ni'

    Returns:
        ndarray: Fundamental matrix R
    """
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    D2 = np.asarray(D2)
    D0_arr = native_to_array(D0)
    D1_arr = native_to_array(D1)
    D2_arr = native_to_array(D2)

    try:
        if method.lower() == 'lr':
                D0_arr, D1_arr, D2_arr
            )
        elif method.lower() == 'fi':
                D0_arr, D1_arr, D2_arr
            )
        elif method.lower() == 'is':
                D0_arr, D1_arr, D2_arr
            )
        elif method.lower() == 'cr':
                D0_arr, D1_arr, D2_arr
            )
        elif method.lower() == 'eg':
                D0_arr, D1_arr, D2_arr
            )
        elif method.lower() == 'ni':
                D0_arr, D1_arr, D2_arr
            )
        else:
            raise ValueError(f"Unknown method: {method}")

        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"QBD fundamental matrix computation failed: {str(e)}")


def lib_qbd_stationary_probability(D0, D1, D2):
    """
    Compute stationary probability vector for QBD process.

    Args:
        D0 (array-like): Downward transition matrix
        D1 (array-like): Transition within level
        D2 (array-like): Upward transition matrix

    Returns:
        ndarray: Stationary probability vector
    """
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    D2 = np.asarray(D2)
    D0_arr = native_to_array(D0)
    D1_arr = native_to_array(D1)
    D2_arr = native_to_array(D2)

    try:
            D0_arr, D1_arr, D2_arr
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"QBD stationary probability failed: {str(e)}")


def lib_gim1_fundamental_matrix(D0, D1):
    """
    Compute fundamental matrix for G/M/1 queue.

    Args:
        D0 (array-like): Downward transitions
        D1 (array-like): Upward transitions

    Returns:
        ndarray: Fundamental matrix
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
        raise RuntimeError(f"GIM1 fundamental matrix failed: {str(e)}")


def lib_gim1_stationary_probability(D0, D1):
    """
    Compute stationary probability for G/M/1 queue.

    Args:
        D0 (array-like): Downward transitions
        D1 (array-like): Upward transitions

    Returns:
        ndarray: Stationary probability vector
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
        raise RuntimeError(f"GIM1 stationary probability failed: {str(e)}")


def lib_mg1_fundamental_matrix(D0, D1):
    """
    Compute fundamental matrix for M/G/1 queue.

    Args:
        D0 (array-like): Transitions
        D1 (array-like): Completion rate matrix

    Returns:
        ndarray: Fundamental matrix
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
        raise RuntimeError(f"MG1 fundamental matrix failed: {str(e)}")


def lib_mg1_performance_moments(D0, D1):
    """
    Compute M/G/1 performance moments.

    Args:
        D0 (array-like): Transitions
        D1 (array-like): Completion rates

    Returns:
        dict: Performance measures including queue length moments
    """
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)

    try:
        R = lib_mg1_fundamental_matrix(D0, D1)
        # Additional computations for moments
        return {'fundamental_matrix': R}
    except Exception as e:
        raise RuntimeError(f"MG1 moments computation failed: {str(e)}")


def lib_qbd_summary(D0, D1, D2):
    """
    Comprehensive summary of QBD process.

    Args:
        D0, D1, D2: QBD transition matrices

    Returns:
        dict: Summary with key metrics
    """
    try:
        pi = lib_qbd_stationary_probability(D0, D1, D2)
        R = lib_qbd_fundamental_matrix(D0, D1, D2)

        return {
            'stationary_prob': pi,
            'fundamental_matrix': R,
            'dimension': R.shape[0]
        }
    except Exception as e:
        raise RuntimeError(f"QBD summary failed: {str(e)}")
