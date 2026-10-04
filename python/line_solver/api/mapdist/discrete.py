"""
Discrete-time MAP (D-MAP) distance measures.

Reference (continuous-time formulation):
    G. Horvath, "Measuring the distance between MAPs and some
    applications," in Proc. ASMTA 2015, LNCS 9081, pp. 95-109.
    https://link.springer.com/chapter/10.1007/978-3-319-18579-8_8

Discrete-time extension by QORE Lab (https://qore.doc.ic.ac.uk/)
"""

import numpy as np
from scipy.optimize import minimize


def dlyap(A, B, C):
    """Solve the discrete Sylvester equation A X B - X + C = 0.

    Matches MATLAB dlyap(A,B,C) and JAR Matrix.dlyap(A,B,C). scipy's
    solve_discrete_lyapunov only solves the single-matrix form A X A^H - X + Q = 0
    and takes two arguments, so it cannot be used here (its 3rd positional is
    the ``method`` string, which raised AttributeError on the passed matrix).
    """
    n = A.shape[0]
    m = B.shape[0]
    lhs = np.kron(B.T, A) - np.eye(n * m)          # (B^T (x) A) - I
    vec_c = np.asarray(C).flatten(order='F')       # column-major vec(C)
    vec_x = np.linalg.solve(lhs, -vec_c)
    return vec_x.reshape((n, m), order='F')


def _dmap_pie(D0, D1):
    """Stationary vector at arrivals: pi*(I-D0)^{-1}*D1 = pi."""
    n = D0.shape[0]
    P = np.linalg.solve(np.eye(n) - D0, D1)
    # Solve (P^T - I) pi = 0 with sum(pi)=1 by replacing the FIRST ROW
    # (equation 0) with the normalization; the row is the equation, so this
    # must be M[0, :] = 1, not M[:, 0] = 1 (which corrupts the solution).
    M = (P.T - np.eye(n)).copy()
    M[0, :] = 1.0
    rhs = np.zeros(n)
    rhs[0] = 1.0
    return np.linalg.solve(M, rhs)


def dmap_geo_mul_sum(D0A, D1A, D0B, D1B, L, alA=None, alB=None):
    """Joint PMF inner product of two D-MAPs via recursive discrete Lyapunov.

    Args:
        D0A, D1A: D0, D1 matrices of first D-MAP.
        D0B, D1B: D0, D1 matrices of second D-MAP.
        L: Number of inter-arrival times in the joint PMF.
        alA, alB: Optional stationary vectors at arrival epochs.

    Returns:
        float: Inner product of joint PMFs.
    """
    NA, NB = D0A.shape[0], D0B.shape[0]
    if alA is None:
        alA = _dmap_pie(D0A, D1A)
    if alB is None:
        alB = _dmap_pie(D0B, D1B)
    Z = dlyap(D0B.T, D0A, alB.reshape(-1, 1) @ alA.reshape(1, -1))
    for _ in range(L - 1):
        Z = dlyap(D0B.T, D0A, D1B.T @ Z @ D1A)
    dA = (np.eye(NA) - D0A).sum(axis=1)
    dB = (np.eye(NB) - D0B).sum(axis=1)
    return float(dB @ Z @ dA)


def dmap_dist(D0A, D1A, D0B, D1B, L, alA=None, alB=None):
    """Squared L2 distance between lag-L joint PMFs of two D-MAPs.

    Args:
        D0A, D1A: D0, D1 matrices of first D-MAP.
        D0B, D1B: D0, D1 matrices of second D-MAP.
        L: Number of lags.
        alA, alB: Optional stationary vectors at arrival epochs.

    Returns:
        float: Squared L2 distance.
    """
    if alA is None:
        alA = _dmap_pie(D0A, D1A)
    if alB is None:
        alB = _dmap_pie(D0B, D1B)
    return (dmap_geo_mul_sum(D0A, D1A, D0A, D1A, L + 1, alA, alA)
            - 2 * dmap_geo_mul_sum(D0A, D1A, D0B, D1B, L + 1, alA, alB)
            + dmap_geo_mul_sum(D0B, D1B, D0B, D1B, L + 1, alB, alB))


def dmap_dist_lag1(D0A, D1A, D0B, D1B, alA=None, alB=None):
    """Lag-1 joint PMF L2 distance via Kronecker/discrete-Lyapunov.

    Args:
        D0A, D1A: D0, D1 matrices of first D-MAP.
        D0B, D1B: D0, D1 matrices of second D-MAP.
        alA, alB: Optional stationary vectors at arrival epochs.

    Returns:
        float: Squared L2 distance.
    """
    NA, NB = D0A.shape[0], D0B.shape[0]
    if alA is None:
        alA = _dmap_pie(D0A, D1A)
    if alB is None:
        alB = _dmap_pie(D0B, D1B)
    dA = (np.eye(NA) - D0A).sum(axis=1)
    dB = (np.eye(NB) - D0B).sum(axis=1)

    Z_AB = dlyap(D0A.T, D0B, alA.reshape(-1, 1) @ alB.reshape(1, -1))
    Z_AA = dlyap(D0A.T, D0A, alA.reshape(-1, 1) @ alA.reshape(1, -1))
    Z_BB = dlyap(D0B.T, D0B, alB.reshape(-1, 1) @ alB.reshape(1, -1))

    X_AB = dlyap(D0A, D0B.T, dA.reshape(-1, 1) @ dB.reshape(1, -1))
    X_AA = dlyap(D0A, D0A.T, dA.reshape(-1, 1) @ dA.reshape(1, -1))
    X_BB = dlyap(D0B, D0B.T, dB.reshape(-1, 1) @ dB.reshape(1, -1))

    vA1 = D1A.flatten(order='F')
    vB1 = D1B.flatten(order='F')

    return float(vB1 @ np.kron(X_BB, Z_BB) @ vB1
                 + vA1 @ np.kron(X_AA, Z_AA) @ vA1
                 - 2 * vA1 @ np.kron(X_AB, Z_AB) @ vB1)


def dmap_geo_mul_sum_acf(D0A, D1A, D0B, D1B, alA=None, alB=None):
    """Geometric sum for discrete autocorrelation distance.

    Args:
        D0A, D1A: D0, D1 matrices of first D-MAP.
        D0B, D1B: D0, D1 matrices of second D-MAP.
        alA, alB: Optional stationary vectors at arrival epochs.

    Returns:
        float: Geometric sum value.
    """
    NA, NB = D0A.shape[0], D0B.shape[0]
    if alA is None:
        alA = _dmap_pie(D0A, D1A)
    if alB is None:
        alB = _dmap_pie(D0B, D1B)
    D0Ai = np.linalg.inv(np.eye(NA) - D0A)
    D0Bi = np.linalg.inv(np.eye(NB) - D0B)
    PAh = D0Ai @ D1A - np.ones((NA, 1)) @ alA.reshape(1, -1)
    PBh = D0Bi @ D1B - np.ones((NB, 1)) @ alB.reshape(1, -1)
    M = np.eye(NA * NB) - np.kron(PBh.T, PAh)
    if np.linalg.cond(M) > 1e10:
        return np.linalg.cond(M)
    X = dlyap(PAh, PBh, D0Ai.sum(axis=1).reshape(-1, 1) @ (alB @ D0Bi).reshape(1, -1))
    return float((alA @ D0Ai @ X @ D0Bi).sum())


def dmap_dist_acf(D0A, D1A, D0B, D1B, alA=None, alB=None):
    """Squared L2 distance between autocorrelation functions of two D-MAPs.

    Args:
        D0A, D1A: D0, D1 matrices of first D-MAP.
        D0B, D1B: D0, D1 matrices of second D-MAP.
        alA, alB: Optional stationary vectors at arrival epochs.

    Returns:
        float: Squared L2 distance.
    """
    NA, NB = D0A.shape[0], D0B.shape[0]
    if alA is None:
        alA = _dmap_pie(D0A, D1A)
    if alB is None:
        alB = _dmap_pie(D0B, D1B)
    D0Ai = np.linalg.inv(np.eye(NA) - D0A)
    D0Bi = np.linalg.inv(np.eye(NB) - D0B)
    eA = np.ones(NA)
    eB = np.ones(NB)
    muA = float(alA @ D0Ai @ eA)
    muB = float(alB @ D0Bi @ eB)
    m2A = float(alA @ (np.eye(NA) + D0A) @ D0Ai @ D0Ai @ eA)
    m2B = float(alB @ (np.eye(NB) + D0B) @ D0Bi @ D0Bi @ eB)
    varA = m2A - muA ** 2
    varB = m2B - muB ** 2
    cA = (m2A + muA) / 2
    cB = (m2B + muB) / 2
    return ((dmap_geo_mul_sum_acf(D0A, D1A, D0A, D1A, alA, alA) - cA ** 2) / varA ** 2
            - 2 * (dmap_geo_mul_sum_acf(D0A, D1A, D0B, D1B, alA, alB) - cA * cB) / (varA * varB)
            + (dmap_geo_mul_sum_acf(D0B, D1B, D0B, D1B, alB, alB) - cB ** 2) / varB ** 2)


def dmap_optim_dist(D0A, D1A, alA, D0B, alB, L):
    """Find D1B minimizing the lag-L joint PMF distance given fixed D0B.

    Args:
        D0A, D1A: D0, D1 matrices of reference D-MAP.
        alA: Stationary vector at arrivals for reference D-MAP.
        D0B: D0 matrix of approximating D-MAP (fixed).
        alB: Stationary vector at arrivals for approximating D-MAP.
        L: Number of lags.

    Returns:
        tuple: (D1B, d) - optimal D1 matrix and minimum distance.
    """
    NB = D0B.shape[0]
    D0Bi = np.linalg.inv(np.eye(NB) - D0B)
    dB = (np.eye(NB) - D0B).sum(axis=1)
    Aeq = np.vstack([np.kron(np.eye(NB), (alB @ D0Bi).reshape(1, -1)),
                     np.kron(np.ones((1, NB)), np.eye(NB))])
    beq = np.concatenate([alB, dB])

    def objective(x):
        D1Bx = x.reshape(NB, NB, order='F')
        return dmap_dist(D0A, D1A, D0B, D1Bx, L, alA, alB)

    constraints = {'type': 'eq', 'fun': lambda x: Aeq @ x - beq}
    bounds = [(1e-6, None)] * (NB * NB)
    x0 = np.random.rand(NB * NB)
    result = minimize(objective, x0, method='SLSQP', constraints=constraints, bounds=bounds,
                      options={'disp': False, 'maxiter': 1000})
    D1B = result.x.reshape(NB, NB, order='F')
    return D1B, result.fun


def dmap_optim_dist_acf(D0A, D1A, alA, D0B, alB):
    """Find D1B minimizing the autocorrelation distance given fixed D0B.

    Args:
        D0A, D1A: D0, D1 matrices of reference D-MAP.
        alA: Stationary vector at arrivals for reference D-MAP.
        D0B: D0 matrix of approximating D-MAP (fixed).
        alB: Stationary vector at arrivals for approximating D-MAP.

    Returns:
        tuple: (D1B, d) - optimal D1 matrix and minimum distance.
    """
    NB = D0B.shape[0]
    D0Bi = np.linalg.inv(np.eye(NB) - D0B)
    dB = (np.eye(NB) - D0B).sum(axis=1)
    Aeq = np.vstack([np.kron(np.eye(NB), (alB @ D0Bi).reshape(1, -1)),
                     np.kron(np.ones((1, NB)), np.eye(NB))])
    beq = np.concatenate([alB, dB])

    def objective(x):
        D1Bx = x.reshape(NB, NB, order='F')
        return dmap_dist_acf(D0A, D1A, D0B, D1Bx, alA, alB)

    constraints = {'type': 'eq', 'fun': lambda x: Aeq @ x - beq}
    bounds = [(1e-6, None)] * (NB * NB)
    x0 = np.random.rand(NB * NB)
    result = minimize(objective, x0, method='SLSQP', constraints=constraints, bounds=bounds,
                      options={'disp': False, 'maxiter': 1000})
    D1B = result.x.reshape(NB, NB, order='F')
    d = dmap_dist_acf(D0A, D1A, D0B, D1B, alA, alB)
    return D1B, d
