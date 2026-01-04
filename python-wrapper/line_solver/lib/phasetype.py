"""
Phase-type distribution manipulation utilities.

Provides direct access to phase-type (PH) and matrix-exponential (ME) operations:
- Creating distributions from moments
- Computing probability measures (CDF, PDF)
- Fitting and optimization
- Canonical forms
- Distribution properties
"""

import jpype
import numpy as np
from line_solver import jlineMatrixFromArray, jlineMatrixToArray


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
    D0_jline = jlineMatrixFromArray(D0)
    D1_jline = jlineMatrixFromArray(D1)

    try:
        result = jpype.JPackage('jline').lib.butools.BuToolsME.momentsFromME(
            D0_jline, D1_jline
        )
        return jlineMatrixToArray(result)
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
    alpha_jline = jlineMatrixFromArray(alpha)
    A_jline = jlineMatrixFromArray(A)

    try:
        result = jpype.JPackage('jline').lib.butools.BuToolsPH.momentsFromPH(
            alpha_jline, A_jline
        )
        return jlineMatrixToArray(result)
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
        result_obj = jpype.JPackage('jline').lib.butools.BuToolsME.meFromMoments(
            jlineMatrixFromArray(moments)
        )
        D0 = jlineMatrixToArray(result_obj.first)
        D1 = jlineMatrixToArray(result_obj.second)

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
        result_obj = jpype.JPackage('jline').lib.butools.BuToolsAPH.ph2From3Moments(
            float(m1), float(m2), float(m3)
        )
        alpha = jlineMatrixToArray(result_obj.first)
        A = jlineMatrixToArray(result_obj.second)

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
        result_obj = jpype.JPackage('jline').lib.butools.BuToolsAPH.ph3From5Moments(
            jlineMatrixFromArray(moments)
        )
        alpha = jlineMatrixToArray(result_obj.first)
        A = jlineMatrixToArray(result_obj.second)

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
    alpha_jline = jlineMatrixFromArray(alpha)
    A_jline = jlineMatrixFromArray(A)
    x_jline = jlineMatrixFromArray(x)

    try:
        result = jpype.JPackage('jline').lib.butools.BuToolsPH.cdfFromPH(
            alpha_jline, A_jline, x_jline
        )
        result_array = jlineMatrixToArray(result)
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
    alpha_jline = jlineMatrixFromArray(alpha)
    A_jline = jlineMatrixFromArray(A)
    x_jline = jlineMatrixFromArray(x)

    try:
        result = jpype.JPackage('jline').lib.butools.BuToolsPH.pdfFromPH(
            alpha_jline, A_jline, x_jline
        )
        result_array = jlineMatrixToArray(result)
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
    D0_jline = jlineMatrixFromArray(D0)
    D1_jline = jlineMatrixFromArray(D1)
    x_jline = jlineMatrixFromArray(x)

    try:
        result = jpype.JPackage('jline').lib.butools.BuToolsME.cdfFromME(
            D0_jline, D1_jline, x_jline
        )
        result_array = jlineMatrixToArray(result)
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
    D0_jline = jlineMatrixFromArray(D0)
    D1_jline = jlineMatrixFromArray(D1)
    x_jline = jlineMatrixFromArray(x)

    try:
        result = jpype.JPackage('jline').lib.butools.BuToolsME.pdfFromME(
            D0_jline, D1_jline, x_jline
        )
        result_array = jlineMatrixToArray(result)
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
    alpha_jline = jlineMatrixFromArray(alpha)
    A_jline = jlineMatrixFromArray(A)

    try:
        result_obj = jpype.JPackage('jline').lib.butools.BuToolsPH.canonicalFromPH2(
            alpha_jline, A_jline
        )
        alpha_c = jlineMatrixToArray(result_obj.first)
        A_c = jlineMatrixToArray(result_obj.second)

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
    alpha_jline = jlineMatrixFromArray(alpha)
    A_jline = jlineMatrixFromArray(A)

    try:
        result_obj = jpype.JPackage('jline').lib.butools.BuToolsPH.canonicalFromPH3(
            alpha_jline, A_jline
        )
        alpha_c = jlineMatrixToArray(result_obj.first)
        A_c = jlineMatrixToArray(result_obj.second)

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
    alpha_jline = jlineMatrixFromArray(alpha)
    A_jline = jlineMatrixFromArray(A)

    try:
        result = jpype.JPackage('jline').lib.butools.BuToolsPH.checkPHRepresentation(
            alpha_jline, A_jline
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
    D0_jline = jlineMatrixFromArray(D0)
    D1_jline = jlineMatrixFromArray(D1)

    try:
        result = jpype.JPackage('jline').lib.butools.BuToolsME.checkMERepresentation(
            D0_jline, D1_jline
        )
        return bool(result)
    except Exception as e:
        return False
