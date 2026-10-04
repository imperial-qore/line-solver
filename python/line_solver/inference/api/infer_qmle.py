import numpy as np


def infer_qmle(Q, N, Z):
    """Quick MLE closed-form demand estimation.

    D(i,j) = Q(i,j) / (N(j) - sum(Q(:,j))) * Z(j) / (1 + sum(Q(i,:)) - Q(i,j)/N(j))

    Args:
        Q: M x R matrix of mean queue lengths
        N: 1-D array of populations per class (length R)
        Z: 1-D array of think times per class (length R)

    Returns:
        D: M x R matrix of estimated demands

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    Q = np.asarray(Q, dtype=float)
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()

    M, R = Q.shape
    D = np.zeros((M, R))
    for i in range(M):
        for j in range(R):
            denom1 = N[j] - np.sum(Q[:, j])
            denom2 = 1.0 + np.sum(Q[i, :]) - Q[i, j] / N[j]
            if denom1 != 0 and denom2 != 0:
                D[i, j] = Q[i, j] / denom1 * Z[j] / denom2
    return D
