"""
Multivariate Phase-Type (MVPH) Analysis.

Functions for analyzing phase-type distributions with multiple outputs:
- Joint moment computation
- Cross-correlations and covariances
- Multivariate distribution properties
"""

import jpype
import numpy as np
from line_solver import jlineMatrixFromArray, jlineMatrixToArray


def lib_mvph_mean_output1(alpha, A, B):
    """
    Compute mean of first output in multivariate PH.

    Args:
        alpha (array-like): Initial distribution (1 x n)
        A (array-like): Generator matrix (n x n)
        B (array-like): Output matrix (n x 1)

    Returns:
        float: Mean of first output
    """
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    B = np.asarray(B)
    alpha_jline = jlineMatrixFromArray(alpha)
    A_jline = jlineMatrixFromArray(A)
    B_jline = jlineMatrixFromArray(B)

    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.mvph_mean_x(
            alpha_jline, A_jline, B_jline
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"mvph_mean_x failed: {str(e)}")


def lib_mvph_mean_output2(alpha, A, C):
    """
    Compute mean of second output in multivariate PH.

    Args:
        alpha (array-like): Initial distribution (1 x n)
        A (array-like): Generator matrix (n x n)
        C (array-like): Output matrix (n x 1)

    Returns:
        float: Mean of second output
    """
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    C = np.asarray(C)
    alpha_jline = jlineMatrixFromArray(alpha)
    A_jline = jlineMatrixFromArray(A)
    C_jline = jlineMatrixFromArray(C)

    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.mvph_mean_y(
            alpha_jline, A_jline, C_jline
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"mvph_mean_y failed: {str(e)}")


def lib_mvph_covariance(alpha, A, B, C):
    """
    Compute covariance between two outputs in multivariate PH.

    Args:
        alpha (array-like): Initial distribution (1 x n)
        A (array-like): Generator matrix (n x n)
        B (array-like): First output matrix (n x 1)
        C (array-like): Second output matrix (n x 1)

    Returns:
        float: Covariance between outputs
    """
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    B = np.asarray(B)
    C = np.asarray(C)
    alpha_jline = jlineMatrixFromArray(alpha)
    A_jline = jlineMatrixFromArray(A)
    B_jline = jlineMatrixFromArray(B)
    C_jline = jlineMatrixFromArray(C)

    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.mvph_cov(
            alpha_jline, A_jline, B_jline, C_jline
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"mvph_cov failed: {str(e)}")


def lib_mvph_correlation(alpha, A, B, C):
    """
    Compute correlation between two outputs in multivariate PH.

    Args:
        alpha (array-like): Initial distribution (1 x n)
        A (array-like): Generator matrix (n x n)
        B (array-like): First output matrix (n x 1)
        C (array-like): Second output matrix (n x 1)

    Returns:
        float: Correlation coefficient between outputs
    """
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    B = np.asarray(B)
    C = np.asarray(C)
    alpha_jline = jlineMatrixFromArray(alpha)
    A_jline = jlineMatrixFromArray(A)
    B_jline = jlineMatrixFromArray(B)
    C_jline = jlineMatrixFromArray(C)

    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.mvph_corr(
            alpha_jline, A_jline, B_jline, C_jline
        )
        return float(result)
    except Exception as e:
        raise RuntimeError(f"mvph_corr failed: {str(e)}")


def lib_mvph_joint_distribution(alpha, A, B, C, x_vals, y_vals):
    """
    Compute joint distribution of two outputs in multivariate PH.

    Args:
        alpha (array-like): Initial distribution
        A (array-like): Generator matrix
        B (array-like): First output matrix
        C (array-like): Second output matrix
        x_vals (array-like): X values
        y_vals (array-like): Y values

    Returns:
        ndarray: Joint probability values at grid points
    """
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    B = np.asarray(B)
    C = np.asarray(C)
    x_vals = np.asarray(x_vals)
    y_vals = np.asarray(y_vals)

    alpha_jline = jlineMatrixFromArray(alpha)
    A_jline = jlineMatrixFromArray(A)
    B_jline = jlineMatrixFromArray(B)
    C_jline = jlineMatrixFromArray(C)
    x_jline = jlineMatrixFromArray(x_vals)
    y_jline = jlineMatrixFromArray(y_vals)

    try:
        result = jpype.JPackage('jline').lib.kpctoolbox.Kpctoolbox.mvph_joint(
            alpha_jline, A_jline, B_jline, C_jline, x_jline, y_jline
        )
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"mvph_joint failed: {str(e)}")


def lib_mvph_summary(alpha, A, B, C):
    """
    Compute comprehensive summary of multivariate PH distribution.

    Args:
        alpha (array-like): Initial distribution
        A (array-like): Generator matrix
        B (array-like): First output matrix
        C (array-like): Second output matrix

    Returns:
        dict: Summary statistics
    """
    mean_x = lib_mvph_mean_output1(alpha, A, B)
    mean_y = lib_mvph_mean_output2(alpha, A, C)
    cov = lib_mvph_covariance(alpha, A, B, C)
    corr = lib_mvph_correlation(alpha, A, B, C)

    return {
        'mean_x': mean_x,
        'mean_y': mean_y,
        'covariance': cov,
        'correlation': corr
    }
