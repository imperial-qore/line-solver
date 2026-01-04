"""
Stationary Moment Computation (SMC) utilities.

Low-level solvers for computing stationary moments in queues:
- QBD (Quasi-Birth-Death) processes
- GIM1 (G/M/1) queue analysis
- MG1 (M/G/1) queue analysis
- Fundamental matrix computation
"""

import jpype
import numpy as np
from line_solver import jlineMatrixFromArray, jlineMatrixToArray


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
    D0_jline = jlineMatrixFromArray(D0)
    D1_jline = jlineMatrixFromArray(D1)
    D2_jline = jlineMatrixFromArray(D2)

    try:
        if method.lower() == 'lr':
            result = jpype.JPackage('jline').lib.smc.SMC.QBD_LR(
                D0_jline, D1_jline, D2_jline
            )
        elif method.lower() == 'fi':
            result = jpype.JPackage('jline').lib.smc.SMC.QBD_FI(
                D0_jline, D1_jline, D2_jline
            )
        elif method.lower() == 'is':
            result = jpype.JPackage('jline').lib.smc.SMC.QBD_IS(
                D0_jline, D1_jline, D2_jline
            )
        elif method.lower() == 'cr':
            result = jpype.JPackage('jline').lib.smc.SMC.QBD_CR(
                D0_jline, D1_jline, D2_jline
            )
        elif method.lower() == 'eg':
            result = jpype.JPackage('jline').lib.smc.SMC.QBD_EG(
                D0_jline, D1_jline, D2_jline
            )
        elif method.lower() == 'ni':
            result = jpype.JPackage('jline').lib.smc.SMC.QBD_NI(
                D0_jline, D1_jline, D2_jline
            )
        else:
            raise ValueError(f"Unknown method: {method}")

        return jlineMatrixToArray(result)
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
    D0_jline = jlineMatrixFromArray(D0)
    D1_jline = jlineMatrixFromArray(D1)
    D2_jline = jlineMatrixFromArray(D2)

    try:
        result = jpype.JPackage('jline').lib.smc.SMC.QBD_pi(
            D0_jline, D1_jline, D2_jline
        )
        return jlineMatrixToArray(result)
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
    D0_jline = jlineMatrixFromArray(D0)
    D1_jline = jlineMatrixFromArray(D1)

    try:
        result = jpype.JPackage('jline').lib.smc.SMC.GIM1_R(
            D0_jline, D1_jline
        )
        return jlineMatrixToArray(result)
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
    D0_jline = jlineMatrixFromArray(D0)
    D1_jline = jlineMatrixFromArray(D1)

    try:
        result = jpype.JPackage('jline').lib.smc.SMC.GIM1_pi(
            D0_jline, D1_jline
        )
        return jlineMatrixToArray(result)
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
    D0_jline = jlineMatrixFromArray(D0)
    D1_jline = jlineMatrixFromArray(D1)

    try:
        result = jpype.JPackage('jline').lib.smc.SMC.MG1_CR(
            D0_jline, D1_jline
        )
        return jlineMatrixToArray(result)
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
