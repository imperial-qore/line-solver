"""
Phase-type distribution manipulation utilities.

Provides direct access to phase-type (PH) and matrix-exponential (ME) operations:
- Creating distributions from moments
- Computing probability measures (CDF, PDF)
- Fitting and optimization
- Canonical forms
- Distribution properties
"""

import numpy as np
from line_solver import native_to_array


def lib_ph_moments_from_matrix_exponential(D0, D1):
    """
    Compute moments from matrix-exponential representation.

    Args:
        D0 (array-like): Generator matrix
        D1 (array-like): Transition matrix

    Returns:
        ndarray: First few moments
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
        raise RuntimeError(f"momentsFromME failed: {str(e)}")


def lib_ph_moments_from_phasetype(alpha, A):
    """
    Compute moments from phase-type representation.

    Args:
        alpha (array-like): Initial distribution (1 x n)
        A (array-like): Generator matrix (n x n)

    Returns:
        ndarray: Moments of the phase-type distribution
    """
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    alpha_arr = native_to_array(alpha)
    A_arr = native_to_array(A)

    try:
            alpha_arr, A_arr
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"momentsFromPH failed: {str(e)}")


def lib_ph_matrix_exponential_from_moments(moments):
    """
    Fit matrix-exponential distribution to moments.

    Args:
        moments (array-like): First n moments

    Returns:
        dict: Dictionary with 'D0' and 'D1' matrices
    """
    moments = np.asarray(moments).flatten()

    try:
            native_to_array(moments)
        )
        D0 = native_to_array(result_obj.first)
        D1 = native_to_array(result_obj.second)

        return {'D0': D0, 'D1': D1}
    except Exception as e:
        raise RuntimeError(f"meFromMoments failed: {str(e)}")


def lib_ph_phasetype2_from_moments(m1, m2, m3):
    """
    Fit 2-phase acyclic phase-type distribution from 3 moments.

    Args:
        m1 (float): First moment
        m2 (float): Second moment
        m3 (float): Third moment

    Returns:
        dict: Dictionary with 'alpha' and 'A' matrices
    """
    try:
            float(m1), float(m2), float(m3)
        )
        alpha = native_to_array(result_obj.first)
        A = native_to_array(result_obj.second)

        return {'alpha': alpha, 'A': A}
    except Exception as e:
        raise RuntimeError(f"ph2From3Moments failed: {str(e)}")


def lib_ph_phasetype3_from_moments(moments):
    """
    Fit 3-phase acyclic phase-type distribution from 5 moments.

    Args:
        moments (array-like): First 5 moments

    Returns:
        dict: Dictionary with 'alpha' and 'A' matrices
    """
    moments = np.asarray(moments).flatten()

    try:
            native_to_array(moments)
        )
        alpha = native_to_array(result_obj.first)
        A = native_to_array(result_obj.second)

        return {'alpha': alpha, 'A': A}
    except Exception as e:
        raise RuntimeError(f"ph3From5Moments failed: {str(e)}")


def lib_ph_cdf(alpha, A, x):
    """
    Compute cumulative distribution function (CDF) of phase-type.

    Args:
        alpha (array-like): Initial distribution
        A (array-like): Generator matrix
        x (float or array-like): Points to evaluate

    Returns:
        float or ndarray: CDF values
    """
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    x = np.atleast_1d(x)
    alpha_arr = native_to_array(alpha)
    A_arr = native_to_array(A)
    x_arr = native_to_array(x)

    try:
            alpha_arr, A_arr, x_arr
        )
        result_array = native_to_array(result)
        return result_array[0] if len(result_array) == 1 else result_array
    except Exception as e:
        raise RuntimeError(f"cdfFromPH failed: {str(e)}")


def lib_ph_pdf(alpha, A, x):
    """
    Compute probability density function (PDF) of phase-type.

    Args:
        alpha (array-like): Initial distribution
        A (array-like): Generator matrix
        x (float or array-like): Points to evaluate

    Returns:
        float or ndarray: PDF values
    """
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    x = np.atleast_1d(x)
    alpha_arr = native_to_array(alpha)
    A_arr = native_to_array(A)
    x_arr = native_to_array(x)

    try:
            alpha_arr, A_arr, x_arr
        )
        result_array = native_to_array(result)
        return result_array[0] if len(result_array) == 1 else result_array
    except Exception as e:
        raise RuntimeError(f"pdfFromPH failed: {str(e)}")


def lib_ph_cdf_matrix_exponential(D0, D1, x):
    """
    Compute CDF of matrix-exponential distribution.

    Args:
        D0 (array-like): Generator matrix
        D1 (array-like): Transition matrix
        x (float or array-like): Points to evaluate

    Returns:
        float or ndarray: CDF values
    """
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    x = np.atleast_1d(x)
    D0_arr = native_to_array(D0)
    D1_arr = native_to_array(D1)
    x_arr = native_to_array(x)

    try:
            D0_arr, D1_arr, x_arr
        )
        result_array = native_to_array(result)
        return result_array[0] if len(result_array) == 1 else result_array
    except Exception as e:
        raise RuntimeError(f"cdfFromME failed: {str(e)}")


def lib_ph_pdf_matrix_exponential(D0, D1, x):
    """
    Compute PDF of matrix-exponential distribution.

    Args:
        D0 (array-like): Generator matrix
        D1 (array-like): Transition matrix
        x (float or array-like): Points to evaluate

    Returns:
        float or ndarray: PDF values
    """
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    x = np.atleast_1d(x)
    D0_arr = native_to_array(D0)
    D1_arr = native_to_array(D1)
    x_arr = native_to_array(x)

    try:
            D0_arr, D1_arr, x_arr
        )
        result_array = native_to_array(result)
        return result_array[0] if len(result_array) == 1 else result_array
    except Exception as e:
        raise RuntimeError(f"pdfFromME failed: {str(e)}")


def lib_ph_canonical_2phase(alpha, A):
    """
    Convert 2-phase to canonical form.

    Args:
        alpha (array-like): Initial distribution
        A (array-like): Generator matrix

    Returns:
        dict: Canonical form with 'alpha' and 'A'
    """
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    alpha_arr = native_to_array(alpha)
    A_arr = native_to_array(A)

    try:
            alpha_arr, A_arr
        )
        alpha_c = native_to_array(result_obj.first)
        A_c = native_to_array(result_obj.second)

        return {'alpha': alpha_c, 'A': A_c}
    except Exception as e:
        raise RuntimeError(f"canonicalFromPH2 failed: {str(e)}")


def lib_ph_canonical_3phase(alpha, A):
    """
    Convert 3-phase to canonical form.

    Args:
        alpha (array-like): Initial distribution
        A (array-like): Generator matrix

    Returns:
        dict: Canonical form with 'alpha' and 'A'
    """
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    alpha_arr = native_to_array(alpha)
    A_arr = native_to_array(A)

    try:
            alpha_arr, A_arr
        )
        alpha_c = native_to_array(result_obj.first)
        A_c = native_to_array(result_obj.second)

        return {'alpha': alpha_c, 'A': A_c}
    except Exception as e:
        raise RuntimeError(f"canonicalFromPH3 failed: {str(e)}")


def lib_ph_check_representation(alpha, A):
    """
    Validate phase-type representation.

    Args:
        alpha (array-like): Initial distribution
        A (array-like): Generator matrix

    Returns:
        bool: True if valid phase-type representation
    """
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    alpha_arr = native_to_array(alpha)
    A_arr = native_to_array(A)

    try:
            alpha_arr, A_arr
        )
        return bool(result)
    except Exception as e:
        return False


def lib_ph_check_matrix_exponential(D0, D1):
    """
    Validate matrix-exponential representation.

    Args:
        D0 (array-like): Generator matrix
        D1 (array-like): Transition matrix

    Returns:
        bool: True if valid matrix-exponential representation
    """
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    D0_arr = native_to_array(D0)
    D1_arr = native_to_array(D1)

    try:
            D0_arr, D1_arr
        )
        return bool(result)
    except Exception as e:
        return False
