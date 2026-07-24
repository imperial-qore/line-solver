"""
Multivariate Phase-Type Distribution (MVPH) functions.

Ported from MATLAB: matlab/lib/kpctoolbox/mvph/
JAR reference: jar/src/main/kotlin/jline/lib/kpctoolbox/mvph/MVPH.kt

Provides functions for computing joint moments and derived statistics
(mean, covariance, correlation) of bivariate phase-type distributions.
"""

import numpy as np
from scipy.linalg import inv
from math import factorial, sqrt


def mvph_joint(alpha, S, T, D, n1, n2):
    """
    Compute the joint moment E[X^n1 * Y^n2] of a bivariate phase-type
    distribution.

    The bivariate phase-type distribution has:
    - Initial vector alpha
    - First phase generator matrix S (for X)
    - Transition matrix D between phases (from X completion to Y start)
    - Second phase generator matrix T (for Y)

    Formula:
        E[X^n1 * Y^n2] = n1! * n2! * alpha * inv(-S)^(n1+1) * D
                          * inv(-T)^(n2+1) * (-T) * e

    Parameters
    ----------
    alpha : array_like
        Initial probability row vector, shape (nS,) or (1, nS).
    S : array_like
        First phase generator matrix (sub-generator for X), shape (nS, nS).
    T : array_like
        Second phase generator matrix (sub-generator for Y), shape (nT, nT).
    D : array_like
        Transition matrix between phases, shape (nS, nT).
    n1 : int
        Power for first variable X (non-negative integer).
    n2 : int
        Power for second variable Y (non-negative integer).

    Returns
    -------
    float
        Joint moment E[X^n1 * Y^n2].

    Notes
    -----
    The MATLAB source (mvph_joint.m line 6) has a typo using n1+1 for the
    second inverse power instead of n2+1. The correct formula (as documented
    in the JAR implementation) uses n2+1 for inv(-T), which is implemented
    here.
    """
    alpha = np.atleast_1d(np.asarray(alpha, dtype=float)).ravel()
    S = np.atleast_2d(np.asarray(S, dtype=float))
    T_mat = np.atleast_2d(np.asarray(T, dtype=float))
    D = np.atleast_2d(np.asarray(D, dtype=float))

    nT = T_mat.shape[0]

    # Compute inv(-S) and inv(-T)
    inv_neg_S = inv(-S)
    inv_neg_T = inv(-T_mat)

    # Compute inv(-S)^(n1+1)
    inv_neg_S_pow = np.linalg.matrix_power(inv_neg_S, n1 + 1)

    # Compute inv(-T)^(n2+1)
    inv_neg_T_pow = np.linalg.matrix_power(inv_neg_T, n2 + 1)

    # Chain multiplication: alpha * inv(-S)^(n1+1) * D * inv(-T)^(n2+1) * (-T) * e
    ones_vec = np.ones(nT)

    result = alpha @ inv_neg_S_pow @ D @ inv_neg_T_pow @ (-T_mat) @ ones_vec

    return factorial(n1) * factorial(n2) * float(result)


def mvph_mean_x(alpha, S, T, D):
    """
    Compute the mean of the first variable in a bivariate PH distribution.

    E[X] = mvph_joint(alpha, S, T, D, 1, 0)

    Parameters
    ----------
    alpha : array_like
        Initial probability row vector.
    S : array_like
        First phase generator matrix (for X).
    T : array_like
        Second phase generator matrix (for Y).
    D : array_like
        Transition matrix between phases.

    Returns
    -------
    float
        Mean E[X].
    """
    return mvph_joint(alpha, S, T, D, 1, 0)


def mvph_mean_y(alpha, S, T, D):
    """
    Compute the mean of the second variable in a bivariate PH distribution.

    E[Y] = mvph_joint(alpha, S, T, D, 0, 1)

    Parameters
    ----------
    alpha : array_like
        Initial probability row vector.
    S : array_like
        First phase generator matrix (for X).
    T : array_like
        Second phase generator matrix (for Y).
    D : array_like
        Transition matrix between phases.

    Returns
    -------
    float
        Mean E[Y].
    """
    return mvph_joint(alpha, S, T, D, 0, 1)


def mvph_cov(alpha, S, T, D):
    """
    Compute the covariance of a bivariate PH distribution.

    Cov(X, Y) = E[XY] - E[X] * E[Y]

    Parameters
    ----------
    alpha : array_like
        Initial probability row vector.
    S : array_like
        First phase generator matrix (for X).
    T : array_like
        Second phase generator matrix (for Y).
    D : array_like
        Transition matrix between phases.

    Returns
    -------
    float
        Covariance Cov(X, Y).
    """
    exy = mvph_joint(alpha, S, T, D, 1, 1)
    ex = mvph_joint(alpha, S, T, D, 1, 0)
    ey = mvph_joint(alpha, S, T, D, 0, 1)
    return exy - ex * ey


def mvph_corr(alpha, S, T, D):
    """
    Compute the correlation of a bivariate PH distribution.

    Corr(X, Y) = Cov(X, Y) / (StdDev(X) * StdDev(Y))

    Parameters
    ----------
    alpha : array_like
        Initial probability row vector.
    S : array_like
        First phase generator matrix (for X).
    T : array_like
        Second phase generator matrix (for Y).
    D : array_like
        Transition matrix between phases.

    Returns
    -------
    float
        Correlation Corr(X, Y). Returns 0.0 if either variance is
        non-positive.
    """
    ex = mvph_joint(alpha, S, T, D, 1, 0)
    ey = mvph_joint(alpha, S, T, D, 0, 1)
    ex2 = mvph_joint(alpha, S, T, D, 2, 0)
    ey2 = mvph_joint(alpha, S, T, D, 0, 2)
    exy = mvph_joint(alpha, S, T, D, 1, 1)

    var_x = ex2 - ex * ex
    var_y = ey2 - ey * ey
    cov_xy = exy - ex * ey

    if var_x <= 0 or var_y <= 0:
        return 0.0

    return cov_xy / (sqrt(var_x) * sqrt(var_y))
