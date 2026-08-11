"""
Continuous-time MAP distance measures.

Reference:
    G. Horvath, "Measuring the distance between MAPs and some
    applications," in Proc. ASMTA 2015, LNCS 9081, pp. 95-109.
    https://link.springer.com/chapter/10.1007/978-3-319-18579-8_8
"""

import numpy as np
from scipy.optimize import minimize


def lyap(A, B, C):
    """Solve the continuous Sylvester equation A X + X B + C = 0.

    Matches MATLAB lyap(A,B,C). scipy's solve_continuous_lyapunov solves only
    the single-matrix form A X + X A^H + Q = 0 and takes two arguments, so it
    cannot be used for the two-matrix Sylvester equation here.
    """
    n = A.shape[0]
    m = B.shape[0]
    lhs = np.kron(np.eye(m), A) + np.kron(B.T, np.eye(n))  # I(x)A + B^T(x)I
    vec_c = np.asarray(C).flatten(order='F')
    vec_x = np.linalg.solve(lhs, -vec_c)
    return vec_x.reshape((n, m), order='F')


def dlyap(A, B, C):
    """Solve the discrete Sylvester equation A X B - X + C = 0 (MATLAB dlyap)."""
    n = A.shape[0]
    m = B.shape[0]
    lhs = np.kron(B.T, A) - np.eye(n * m)
    vec_c = np.asarray(C).flatten(order='F')
    vec_x = np.linalg.solve(lhs, -vec_c)
    return vec_x.reshape((n, m), order='F')


def _map_pie(D0, D1):
    """Stationary vector at arrival epochs: pi such that pi * P = pi, P = (-D0)^{-1} D1."""
    P = np.linalg.solve(-D0, D1)
    n = P.shape[0]
    # Normalization sum(pi)=1 must replace the first ROW (equation 0), not the
    # first column: the row is the equation in M @ pi = rhs.
    M = (P.T - np.eye(n)).copy()
    M[0, :] = 1.0
    rhs = np.zeros(n)
    rhs[0] = 1.0
    return np.linalg.solve(M, rhs)


def map_exp_mul_int(A0, A1, B0, B1, L, alA=None, alB=None):
    """Joint density inner product of two MAPs via recursive Lyapunov equations.

    Args:
        A0, A1: D0, D1 matrices of first MAP.
        B0, B1: D0, D1 matrices of second MAP.
        L: Number of inter-arrival times in the joint density.
        alA, alB: Optional stationary vectors at arrival epochs.

    Returns:
        float: Inner product of the joint densities.
    """
    if alA is None:
        alA = _map_pie(A0, A1)
    if alB is None:
        alB = _map_pie(B0, B1)
    # lyap(A, B, Q) solves A @ X + X @ B + Q = 0
    Z = lyap(B0.T, A0, alB.reshape(-1, 1) @ alA.reshape(1, -1))
    for _ in range(L - 1):
        Z = lyap(B0.T, A0, B1.T @ Z @ A1)
    exitA = (-A0).sum(axis=1)
    exitB = (-B0).sum(axis=1)
    return float(exitB @ Z @ exitA)


def map_dist(A0, A1, B0, B1, L, alA=None, alB=None):
    """Squared L2 distance between lag-L joint densities of two MAPs.

    Args:
        A0, A1: D0, D1 matrices of first MAP.
        B0, B1: D0, D1 matrices of second MAP.
        L: Number of lags (L=1 for lag-1 joint density distance).
        alA, alB: Optional stationary vectors at arrival epochs.

    Returns:
        float: Squared L2 distance.
    """
    if alA is None:
        alA = _map_pie(A0, A1)
    if alB is None:
        alB = _map_pie(B0, B1)
    return (map_exp_mul_int(A0, A1, A0, A1, L + 1, alA, alA)
            - 2 * map_exp_mul_int(A0, A1, B0, B1, L + 1, alA, alB)
            + map_exp_mul_int(B0, B1, B0, B1, L + 1, alB, alB))


def map_dist_lag1(A0, A1, B0, B1, alA=None, alB=None):
    """Lag-1 joint density L2 distance via Kronecker/Lyapunov formulation.

    Equivalent to map_dist(A0, A1, B0, B1, 1) but uses the efficient
    Kronecker product formulation from Theorem 3 of the reference.

    Args:
        A0, A1: D0, D1 matrices of first MAP.
        B0, B1: D0, D1 matrices of second MAP.
        alA, alB: Optional stationary vectors at arrival epochs.

    Returns:
        float: Squared L2 distance of lag-1 joint densities.
    """
    if alA is None:
        alA = _map_pie(A0, A1)
    if alB is None:
        alB = _map_pie(B0, B1)
    a = (-A0).sum(axis=1)
    b = (-B0).sum(axis=1)

    Z_AB = lyap(A0.T, B0, alA.reshape(-1, 1) @ alB.reshape(1, -1))
    Z_AA = lyap(A0.T, A0, alA.reshape(-1, 1) @ alA.reshape(1, -1))
    Z_BB = lyap(B0.T, B0, alB.reshape(-1, 1) @ alB.reshape(1, -1))

    X_AB = lyap(A0, B0.T, a.reshape(-1, 1) @ b.reshape(1, -1))
    X_AA = lyap(A0, A0.T, a.reshape(-1, 1) @ a.reshape(1, -1))
    X_BB = lyap(B0, B0.T, b.reshape(-1, 1) @ b.reshape(1, -1))

    vA1 = A1.flatten(order='F')
    vB1 = B1.flatten(order='F')

    return float(vB1 @ np.kron(X_BB, Z_BB) @ vB1
                 + vA1 @ np.kron(X_AA, Z_AA) @ vA1
                 - 2 * vA1 @ np.kron(X_AB, Z_AB) @ vB1)


def map_geo_mul_sum(A0, A1, B0, B1, alA=None, alB=None):
    """Geometric sum for autocorrelation distance computation.

    Args:
        A0, A1: D0, D1 matrices of first MAP.
        B0, B1: D0, D1 matrices of second MAP.
        alA, alB: Optional stationary vectors at arrival epochs.

    Returns:
        float: Geometric sum value.
    """
    if alA is None:
        alA = _map_pie(A0, A1)
    if alB is None:
        alB = _map_pie(B0, B1)
    A0i = np.linalg.inv(-A0)
    B0i = np.linalg.inv(-B0)
    NA, NB = A0.shape[0], B0.shape[0]
    PAh = A0i @ A1 - np.ones((NA, 1)) @ alA.reshape(1, -1)
    PBh = B0i @ B1 - np.ones((NB, 1)) @ alB.reshape(1, -1)
    M = np.eye(NA * NB) - np.kron(PBh.T, PAh)
    if np.linalg.cond(M) > 1e10:
        return np.linalg.cond(M)
    X = dlyap(PAh, PBh, A0i.sum(axis=1).reshape(-1, 1) @ (alB @ B0i).reshape(1, -1))
    return float((alA @ A0i @ X @ B0i).sum())


def map_dist_acf(A0, A1, B0, B1, alA=None, alB=None):
    """Squared L2 distance between autocorrelation functions of two MAPs.

    Args:
        A0, A1: D0, D1 matrices of first MAP.
        B0, B1: D0, D1 matrices of second MAP.
        alA, alB: Optional stationary vectors at arrival epochs.

    Returns:
        float: Squared L2 distance of autocorrelation functions.
    """
    if alA is None:
        alA = _map_pie(A0, A1)
    if alB is None:
        alB = _map_pie(B0, B1)
    A0i = np.linalg.inv(-A0)
    B0i = np.linalg.inv(-B0)
    e_A = np.ones(A0.shape[0])
    e_B = np.ones(B0.shape[0])
    m1A = float(alA @ A0i @ e_A)
    m2A = float(2 * alA @ A0i @ A0i @ e_A)
    m1B = float(alB @ B0i @ e_B)
    m2B = float(2 * alB @ B0i @ B0i @ e_B)
    varA = m2A - m1A ** 2
    varB = m2B - m1B ** 2
    return ((map_geo_mul_sum(A0, A1, A0, A1, alA, alA) - m2A ** 2 / 4) / varA ** 2
            - 2 * (map_geo_mul_sum(A0, A1, B0, B1, alA, alB) - m2A * m2B / 4) / (varA * varB)
            + (map_geo_mul_sum(B0, B1, B0, B1, alB, alB) - m2B ** 2 / 4) / varB ** 2)


def map_optim_dist(A0, A1, alA, B0, alB, L):
    """Find B1 minimizing the lag-L joint density distance given fixed B0.

    Args:
        A0, A1: D0, D1 matrices of reference MAP.
        alA: Stationary vector at arrivals for reference MAP.
        B0: D0 matrix of approximating MAP (fixed).
        alB: Stationary vector at arrivals for approximating MAP.
        L: Number of lags.

    Returns:
        tuple: (B1, d) - optimal D1 matrix and minimum distance.
    """
    NB = B0.shape[0]
    B0i = np.linalg.inv(-B0)
    b = (-B0).sum(axis=1)
    # Equality constraints: alB * (-B0)^{-1} * B1 = alB, B1 * e = -B0 * e
    Aeq = np.vstack([np.kron(np.eye(NB), (alB @ B0i).reshape(1, -1)),
                     np.kron(np.ones((1, NB)), np.eye(NB))])
    beq = np.concatenate([alB, b])

    def objective(x):
        B1x = x.reshape(NB, NB, order='F')
        return map_dist(A0, A1, B0, B1x, L, alA, alB)

    constraints = {'type': 'eq', 'fun': lambda x: Aeq @ x - beq}
    bounds = [(1e-6, None)] * (NB * NB)
    x0 = np.random.rand(NB * NB)
    result = minimize(objective, x0, method='SLSQP', constraints=constraints, bounds=bounds,
                      options={'disp': False, 'maxiter': 1000})
    B1 = result.x.reshape(NB, NB, order='F')
    d = result.fun
    return B1, d


def map_optim_dist_acf(A0, A1, alA, B0, alB):
    """Find B1 minimizing the autocorrelation distance given fixed B0.

    Args:
        A0, A1: D0, D1 matrices of reference MAP.
        alA: Stationary vector at arrivals for reference MAP.
        B0: D0 matrix of approximating MAP (fixed).
        alB: Stationary vector at arrivals for approximating MAP.

    Returns:
        tuple: (B1, d) - optimal D1 matrix and minimum distance.
    """
    NB = B0.shape[0]
    B0i = np.linalg.inv(-B0)
    b = (-B0).sum(axis=1)
    Aeq = np.vstack([np.kron(np.eye(NB), (alB @ B0i).reshape(1, -1)),
                     np.kron(np.ones((1, NB)), np.eye(NB))])
    beq = np.concatenate([alB, b])

    def objective(x):
        B1x = x.reshape(NB, NB, order='F')
        return map_dist_acf(A0, A1, B0, B1x, alA, alB)

    constraints = {'type': 'eq', 'fun': lambda x: Aeq @ x - beq}
    bounds = [(1e-6, None)] * (NB * NB)
    x0 = np.random.rand(NB * NB)
    result = minimize(objective, x0, method='SLSQP', constraints=constraints, bounds=bounds,
                      options={'disp': False, 'maxiter': 1000})
    B1 = result.x.reshape(NB, NB, order='F')
    d = map_dist_acf(A0, A1, B0, B1, alA, alB)
    return B1, d
