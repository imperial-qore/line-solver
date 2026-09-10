# -*- coding: utf-8 -*-
"""
TransformToOnes: similarity transformation mapping a closing vector to ones.
Ported from BUTools-family fluid tools (G. Horvath).
"""
import numpy as np
import numpy.matlib as ml



__all__ = ["TransformToOnes"]

def TransformToOnes(clovec):
    """
    Returns the similarity transformation matrix B such that B*clovec = ones.
    Works even if clovec has zero entries.

    Parameters
    ----------
    clovec : matrix, shape (M,1)
        The original closing (column) vector.

    Returns
    -------
    B : matrix, shape (M,M)
        The matrix for which B*clovec = ones holds.
    """
    clovec = ml.matrix(clovec).reshape(-1, 1)
    m = clovec.shape[0]
    cv = np.asarray(clovec).flatten()

    # permutation moving a non-zero element to the first position by sorting
    # clovec in descending order (i.e. ascending order of -clovec)
    ix = np.argsort(-cv, kind='stable')
    P = ml.zeros((m, m))
    for i in range(m):
        P[i, ix[i]] = 1.0
    cp = np.asarray(P * clovec).flatten()

    B = ml.zeros((m, m))
    for i in range(m):
        B[i, :i + 1] = 1.0 / np.sum(cp[:i + 1])
    B = B * P
    return B
