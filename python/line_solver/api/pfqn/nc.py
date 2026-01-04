"""
Normalizing Constant methods for Product-Form Queueing Networks.

Native Python implementations of methods for computing normalizing constants:
- Convolution Algorithm (pfqn_ca)
- Related utility functions

References:
    Buzen, J.P. "Computational algorithms for closed queueing networks with
    exponential servers." Communications of the ACM 16.9 (1973): 527-531.
"""

import numpy as np
from math import log, exp, factorial, lgamma
from typing import Tuple, Dict, Any
from functools import lru_cache


def _factln(n: float) -> float:
    """Compute log(n!) using log-gamma function."""
    if n <= 0:
        return 0.0
    return lgamma(n + 1)


def _population_lattice_pprod(n: np.ndarray, N: np.ndarray = None) -> np.ndarray:
    """
    Generate next population vector in lexicographic order.

    Iterates through all population vectors from (0,0,...,0) to N.

    Args:
        n: Current population vector
        N: Maximum population per class (for bounds)

    Returns:
        Next population vector, or (-1,...,-1) when exhausted
    """
    R = len(n)
    n_next = n.copy()

    if N is None:
        # Just increment
        n_next[-1] += 1
        return n_next

    # Find rightmost position that can be incremented
    for i in range(R - 1, -1, -1):
        if n_next[i] < N[i]:
            n_next[i] += 1
            # Reset positions to the right
            for j in range(i + 1, R):
                n_next[j] = 0
            return n_next

    # All exhausted
    return -np.ones(R, dtype=int)


def _hashpop(n: np.ndarray, N: np.ndarray) -> int:
    """
    Compute linear index for population vector.

    Maps population vector n to unique integer index in
    the flattened population lattice [0, prod(N+1)).

    Args:
        n: Population vector
        N: Maximum population per class

    Returns:
        Linear index
    """
    R = len(n)
    idx = 0
    mult = 1
    for i in range(R - 1, -1, -1):
        idx += int(n[i]) * mult
        mult *= int(N[i]) + 1
    return idx


def _pfqn_pff_delay(Z: np.ndarray, n: np.ndarray) -> float:
    """
    Product-form factor for delay stations (think times).

    Computes contribution to normalizing constant from delay stations.

    Args:
        Z: Think times per class
        n: Population vector

    Returns:
        Product-form factor
    """
    R = len(n)
    result = 1.0
    for r in range(R):
        if n[r] > 0:
            # Contribution: Z[r]^n[r] / n[r]!
            result *= (Z[r] ** n[r]) / factorial(int(n[r]))
    return result


def pfqn_ca(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None
            ) -> Tuple[float, float]:
    """
    Convolution Algorithm for normalizing constant computation.

    Computes the normalizing constant G(N) for a closed product-form
    queueing network using Buzen's convolution algorithm.

    Args:
        L: Service demand matrix (M x R) where M is stations, R is classes
        N: Population vector (1 x R or R,) - number of jobs per class
        Z: Think time vector (1 x R or R,) - think time per class (default 0)

    Returns:
        Tuple (G, lG) where:
            - G: Normalizing constant
            - lG: log(G)
    """
    L = np.asarray(L, dtype=np.float64)
    N = np.asarray(N, dtype=np.float64).flatten()
    N = np.ceil(N).astype(int)

    R = len(N)

    if L.ndim == 1:
        L = L.reshape(-1, 1) if R == 1 else L.reshape(1, -1)

    M = L.shape[0]  # Number of stations

    # Handle Z
    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()

    # Special case: no stations (only delay)
    if M == 0:
        # G = prod_r (Z[r]^N[r] / N[r]!)
        lGn = 0.0
        for r in range(R):
            lGn += -_factln(N[r])  # -log(N[r]!)
            if Z[r] > 0 and N[r] > 0:
                lGn += N[r] * log(Z[r])
        Gn = exp(lGn)
        return Gn, lGn

    # Check for negative populations
    if np.any(N < 0):
        return 0.0, float('-inf')

    # Check for zero population
    if N.sum() == 0:
        return 1.0, 0.0

    # Compute total number of population vectors
    product_N_plus_one = int(np.prod(N + 1))

    # G[m, idx] = G_m(n) where idx = hashpop(n, N)
    G = np.ones((M + 1, product_N_plus_one))

    # Iterate through all population vectors
    n = np.zeros(R, dtype=int)  # Start at (0, 0, ..., 0)

    while True:
        # Check if done (n becomes all -1)
        if np.all(n < 0):
            break

        idxn = _hashpop(n, N)

        # Base case: delay station contribution
        G[0, idxn] = _pfqn_pff_delay(Z, n)

        # Convolution recursion: G_m(n) = G_{m-1}(n) + sum_r L[m-1,r] * G_m(n - e_r)
        for m in range(1, M + 1):
            G[m, idxn] = G[m - 1, idxn]
            for r in range(R):
                if n[r] >= 1:
                    n[r] -= 1
                    idxn_1r = _hashpop(n, N)
                    n[r] += 1
                    G[m, idxn] += L[m - 1, r] * G[m, idxn_1r]

        # Next population vector
        n = _population_lattice_pprod(n, N)

    # Final normalizing constant
    Gn = G[M, product_N_plus_one - 1]
    lGn = log(Gn) if Gn > 0 else float('-inf')

    return Gn, lGn


def pfqn_nc(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
            method: str = 'ca') -> Tuple[float, float]:
    """
    Normalizing constant computation dispatcher.

    Selects appropriate algorithm based on method parameter.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,) (default: zeros)
        method: Algorithm to use:
            - 'ca', 'exact': Convolution algorithm
            - 'default': Auto-select based on problem size
            - 'le': Leading eigenvalue asymptotic
            - 'cub': Controllable upper bound
            - 'imci': Importance sampling Monte Carlo integration
            - 'panacea': Hybrid convolution/MVA
            - 'propfair': Proportionally fair allocation
            - 'mmint2': Gauss-Legendre quadrature
            - 'gleint': Gauss-Legendre integration
            - 'sampling': Monte Carlo sampling
            - 'kt': Knessl-Tier expansion
            - 'comom': Conditional moments
            - 'rd': Reduced decomposition
            - 'ls': Linearizer

    Returns:
        Tuple (G, lG) - normalizing constant and its log
    """
    method = method.lower() if method else 'ca'

    if method in ['ca', 'exact']:
        return pfqn_ca(L, N, Z)
    elif method == 'panacea':
        return pfqn_panacea(L, N, Z)
    elif method == 'propfair':
        G, lG, _ = pfqn_propfair(L, N, Z)
        return G, lG
    elif method == 'le':
        from .asymptotic import pfqn_le
        result = pfqn_le(L, N, Z)
        # pfqn_le returns (lG, XN) - extract lG
        lG = result[0] if isinstance(result, tuple) else result
        G = exp(lG) if np.isfinite(lG) else 0.0
        return G, lG
    elif method == 'cub':
        from .asymptotic import pfqn_cub
        result = pfqn_cub(L, N, Z)
        lG = result[0] if isinstance(result, tuple) else result
        G = exp(lG) if np.isfinite(lG) else 0.0
        return G, lG
    elif method == 'imci':
        from .asymptotic import pfqn_mci
        result = pfqn_mci(L, N, Z)
        lG = result[0] if isinstance(result, tuple) else result
        G = exp(lG) if np.isfinite(lG) else 0.0
        return G, lG
    elif method in ['mmint2', 'gleint']:
        from .quadrature import pfqn_mmint2, pfqn_mmint2_gausslegendre
        if method == 'gleint':
            lG, _ = pfqn_mmint2_gausslegendre(L, N, Z)
        else:
            lG, _ = pfqn_mmint2(L, N, Z)
        G = exp(lG) if np.isfinite(lG) else 0.0
        return G, lG
    elif method == 'sampling':
        from .quadrature import pfqn_mmsample2
        lG, _ = pfqn_mmsample2(L, N, Z)
        G = exp(lG) if np.isfinite(lG) else 0.0
        return G, lG
    elif method == 'kt':
        from .kt import pfqn_kt
        lG, _ = pfqn_kt(L, N, Z)
        G = exp(lG) if np.isfinite(lG) else 0.0
        return G, lG
    elif method == 'comom':
        from .comom import pfqn_comom
        result = pfqn_comom(L, N, Z)
        lG = result.lG if hasattr(result, 'lG') else np.nan
        G = exp(lG) if np.isfinite(lG) else 0.0
        return G, lG
    elif method == 'rd':
        from .rd import pfqn_rd
        result = pfqn_rd(L, N, Z)
        lG = result.lG if hasattr(result, 'lG') else np.nan
        G = exp(lG) if np.isfinite(lG) else 0.0
        return G, lG
    elif method == 'ls':
        from .ls import pfqn_ls
        # Logistic sampling approximation
        return pfqn_ls(L, N, Z)
    elif method == 'nrl':
        from .laplace import pfqn_nrl
        lG = pfqn_nrl(L, N, Z)
        G = exp(lG) if np.isfinite(lG) else 0.0
        return G, lG
    elif method == 'nrp':
        from .laplace import pfqn_nrp
        lG = pfqn_nrp(L, N, Z)
        G = exp(lG) if np.isfinite(lG) else 0.0
        return G, lG
    elif method == 'default':
        # Auto-select based on problem size
        L = np.atleast_2d(np.asarray(L, dtype=float))
        N = np.asarray(N, dtype=float).ravel()
        M, R = L.shape
        Ntot = int(np.sum(N))

        # Use CA for small problems, asymptotic for large
        if Ntot * R <= 500:
            return pfqn_ca(L, N, Z)
        else:
            # Use LE for large problems
            from .asymptotic import pfqn_le
            result = pfqn_le(L, N, Z)
            lG = result[0] if isinstance(result, tuple) else result
            G = exp(lG) if np.isfinite(lG) else 0.0
            return G, lG
    else:
        # Default to convolution algorithm
        return pfqn_ca(L, N, Z)


def pfqn_panacea(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None
                 ) -> Tuple[float, float]:
    """
    PANACEA algorithm (hybrid convolution/MVA).

    Currently implemented as wrapper around convolution algorithm.

    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector

    Returns:
        Tuple (G, lG) - normalizing constant and its log
    """
    return pfqn_ca(L, N, Z)


def pfqn_propfair(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None
                  ) -> Tuple[float, float, np.ndarray]:
    """
    Proportionally Fair allocation approximation for normalizing constant.

    Estimates the normalizing constant using a convex optimization program
    that is asymptotically exact in models with single-server PS queues only.

    This method is based on Schweitzer's approach and Walton's proportional
    fairness theory for multi-class networks.

    Args:
        L: Service demand matrix (M x R) where M is stations, R is classes
        N: Population vector (1 x R or R,) - number of jobs per class
        Z: Think time vector (1 x R or R,) - think time per class (default 0)

    Returns:
        Tuple (G, lG, X) where:
            - G: Estimated normalizing constant
            - lG: log(G)
            - X: Asymptotic throughputs per class (1 x R)

    References:
        Schweitzer, P. J. (1979). Approximate analysis of multiclass closed networks
        of queues. In Proceedings of the International Conference on Stochastic
        Control and Optimization.

        Walton, N. (2009). Proportional fairness and its relationship with
        multi-class queueing networks.
    """
    from scipy.optimize import minimize

    L = np.asarray(L, dtype=np.float64)
    N = np.asarray(N, dtype=np.float64).flatten()

    R = len(N)

    if L.ndim == 1:
        L = L.reshape(-1, 1) if R == 1 else L.reshape(1, -1)

    M = L.shape[0]  # Number of stations

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()

    FineTol = 1e-12

    # Objective function: maximize sum_r (N[r] - x[r]*Z[r]) * log(x[r])
    # We minimize the negative
    def objective(x):
        obj = 0.0
        for r in range(R):
            obj += (N[r] - x[r] * Z[r]) * log(abs(x[r]) + FineTol)
        return -obj  # Minimize negative

    # Constraints: sum_r L[i,r] * x[r] <= 1 for all stations i
    # And x[r] >= 0 for all r
    constraints = []

    # Capacity constraints
    for i in range(M):
        constraints.append({
            'type': 'ineq',
            'fun': lambda x, i=i: 1.0 - sum(L[i, r] * x[r] for r in range(R))
        })

    # Non-negativity constraints
    for r in range(R):
        constraints.append({
            'type': 'ineq',
            'fun': lambda x, r=r: x[r]
        })

    # Initial guess - balanced throughput
    x0 = np.zeros(R)
    for r in range(R):
        D_max = L[:, r].max() if M > 0 else 0
        if D_max > 0:
            x0[r] = min(N[r] / (Z[r] + 1), 1.0 / D_max)
        elif Z[r] > 0:
            x0[r] = N[r] / Z[r]
        else:
            x0[r] = N[r]

    # Ensure positive initial guess
    x0 = np.maximum(x0, FineTol)

    # Run optimization with COBYLA
    result = minimize(
        objective,
        x0,
        method='COBYLA',
        constraints=constraints,
        options={'maxiter': 10000, 'rhobeg': 1.0}
    )

    Xasy = result.x

    # Compute lG
    lG = 0.0
    for r in range(R):
        x = Xasy[r]
        if x > FineTol:
            lG += (N[r] - x * Z[r]) * log(1.0 / (x + FineTol))

    # Factorial correction for think times
    for r in range(R):
        thinking = Xasy[r] * Z[r]
        if thinking > 0:
            lG -= _factln(thinking)

    G = exp(lG) if lG > -700 else 0.0  # Avoid underflow

    Xa = Xasy.reshape(1, -1)

    return G, lG, Xa


def pfqn_ls(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
            I: int = 100000) -> Tuple[float, float]:
    """
    Logistic sampling approximation for normalizing constant.

    Approximates the normalizing constant using importance sampling from
    a multivariate normal distribution fitted at the leading eigenvalue mode.

    This method is particularly effective for large networks where
    convolution becomes computationally expensive.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,) (default: zeros)
        I: Number of samples for Monte Carlo integration (default: 100000)

    Returns:
        Tuple (G, lG) where:
            G: Estimated normalizing constant
            lG: log(G)

    Reference:
        G. Casale. "Accelerating performance inference over closed systems by
        asymptotic methods." ACM SIGMETRICS 2017.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()

    # Filter out zero-demand stations
    Lsum = np.sum(L, axis=1)
    L = L[Lsum > 1e-4, :]

    M, R = L.shape

    # Handle empty network
    if L.size == 0 or np.sum(L) < 1e-4 or N.size == 0 or np.sum(N) == 0:
        lGn = -np.sum([_factln(n) for n in N]) + np.sum(N * np.log(np.maximum(np.sum(Z) if Z is not None else 1e-300, 1e-300)))
        return np.exp(lGn), lGn

    if Z is None or len(Z) == 0:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float).flatten()

    # Find the mode using fixed-point iteration
    u, converged = _pfqn_le_fpi(L, N, Z)
    if not converged:
        # Fall back to convolution if mode-finding fails
        return pfqn_ca(L, N, Z)

    Ntot = np.sum(N)

    if np.sum(Z) <= 0:
        # Case without think times
        # Compute Hessian at the mode
        A = _pfqn_le_hessian(L, N, u)
        A = (A + A.T) / 2  # Ensure symmetry

        try:
            iA = np.linalg.inv(A)
        except np.linalg.LinAlgError:
            return pfqn_ca(L, N, Z)

        x0 = np.log(u[:M-1] / u[M-1])

        # Sample from multivariate normal
        samples = np.random.multivariate_normal(x0, iA, I)

        # Evaluate function at samples
        T = np.zeros(I)
        for i in range(I):
            T[i] = _simplex_fun(samples[i, :], L, N)

        # Evaluate PDF at samples
        dpdf = np.zeros(I)
        for i in range(I):
            diff = samples[i, :] - x0
            dpdf[i] = np.exp(-0.5 * diff @ np.linalg.inv(iA) @ diff) / np.sqrt((2 * np.pi) ** (M-1) * np.linalg.det(iA))

        # Compute normalizing constant
        valid = dpdf > 0
        if np.sum(valid) == 0:
            return pfqn_ca(L, N, Z)

        lGn = _multinomialln(np.append(N, M-1)) + _factln(M-1) + np.log(np.mean(T[valid] / dpdf[valid]))
        Gn = np.exp(lGn)

    else:
        # Case with think times Z > 0
        u, v, converged = _pfqn_le_fpiZ(L, N, Z)
        if not converged:
            return pfqn_ca(L, N, Z)

        # Compute Hessian
        A = _pfqn_le_hessianZ(L, N, Z, u, v)
        A = (A + A.T) / 2

        try:
            iA = np.linalg.inv(A)
        except np.linalg.LinAlgError:
            return pfqn_ca(L, N, Z)

        x0 = np.append(np.log(u[:M-1] / u[M-1]), np.log(v))

        # Sample from multivariate normal
        samples = np.random.multivariate_normal(x0, iA, I)

        # Evaluate function at samples
        epsilon = 1e-10
        eN = epsilon * np.sum(N)
        eta = np.sum(N) + M * (1 + eN)
        K = M

        T = np.zeros(I)
        for i in range(I):
            x = samples[i, :]
            term1 = -np.exp(x[K-1]) + K * (1 + eN) * x[M-1]
            term2 = 0.0
            for r in range(R):
                inner = L[K-1, r] * np.exp(x[K-1]) + Z[r]
                for k in range(K-1):
                    inner += np.exp(x[k]) * (L[k, r] * np.exp(x[K-1]) + Z[r])
                term2 += N[r] * np.log(np.maximum(inner, 1e-300))
            term3 = np.sum(x[:K-1])
            term4 = -eta * np.log(1 + np.sum(np.exp(x[:K-1])))
            T[i] = np.exp(term1 + term2 + term3 + term4)

        # Evaluate PDF at samples
        dpdf = np.zeros(I)
        for i in range(I):
            diff = samples[i, :] - x0
            try:
                dpdf[i] = np.exp(-0.5 * diff @ np.linalg.inv(iA) @ diff) / np.sqrt((2 * np.pi) ** len(x0) * np.linalg.det(iA))
            except:
                dpdf[i] = 0

        valid = dpdf > 0
        if np.sum(valid) == 0:
            return pfqn_ca(L, N, Z)

        Gn = np.exp(-np.sum([lgamma(1 + n) for n in N])) * np.mean(T[valid] / dpdf[valid])
        lGn = np.log(Gn) if Gn > 0 else float('-inf')

    return Gn, lGn


def _pfqn_le_fpi(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None):
    """
    Fixed-point iteration to find mode of Gaussian approximation.

    Returns:
        (u, converged): Mode vector and convergence flag
    """
    M, R = L.shape
    Ntot = np.sum(N)
    u = np.ones(M) / M
    u_prev = np.ones(M) * np.inf

    for iteration in range(1000):
        u_prev = u.copy()
        for i in range(M):
            u[i] = 1 / (Ntot + M)
            for r in range(R):
                denom = np.dot(u_prev, L[:, r])
                if denom > 0:
                    u[i] += N[r] / (Ntot + M) * L[i, r] * u_prev[i] / denom

        if np.linalg.norm(u - u_prev, 1) < 1e-10:
            return u, True

    return u, False


def _pfqn_le_fpiZ(L: np.ndarray, N: np.ndarray, Z: np.ndarray):
    """
    Fixed-point iteration with think times.

    Returns:
        (u, v, converged): Mode vector, scale factor, and convergence flag
    """
    M, R = L.shape
    eta = np.sum(N) + M
    u = np.ones(M) / M
    v = eta + 1

    for iteration in range(1000):
        u_prev = u.copy()
        v_prev = v

        for ist in range(M):
            u[ist] = 1 / eta
            for r in range(R):
                denom = Z[r] + v * np.dot(u_prev, L[:, r])
                if denom > 0:
                    u[ist] += (N[r] / eta) * (Z[r] + v * L[ist, r]) * u_prev[ist] / denom

        xi = np.zeros(R)
        for r in range(R):
            denom = Z[r] + v * np.dot(u_prev, L[:, r])
            if denom > 0:
                xi[r] = N[r] / denom

        v = eta + 1 - np.dot(xi, Z)

        if np.linalg.norm(u - u_prev, 1) + abs(v - v_prev) < 1e-10:
            return u, v, True

    return u, v, False


def _pfqn_le_hessian(L: np.ndarray, N: np.ndarray, u: np.ndarray):
    """
    Compute Hessian of Gaussian approximation (without think times).
    """
    M, R = L.shape
    Ntot = np.sum(N)
    hu = np.zeros((M-1, M-1))

    for i in range(M-1):
        for j in range(M-1):
            if i != j:
                hu[i, j] = -(Ntot + M) * u[i] * u[j]
                for r in range(R):
                    denom = np.dot(u, L[:, r]) ** 2
                    if denom > 0:
                        hu[i, j] += N[r] * L[i, r] * L[j, r] * u[i] * u[j] / denom
            else:
                sum_other = np.sum(u) - u[i]
                hu[i, j] = (Ntot + M) * u[i] * sum_other
                for r in range(R):
                    denom = np.dot(u, L[:, r]) ** 2
                    L_other = np.sum(L[:, r]) - L[i, r]
                    if denom > 0:
                        hu[i, j] -= N[r] * L[i, r] * u[i] * (sum_other * L_other) / denom

    return hu


def _pfqn_le_hessianZ(L: np.ndarray, N: np.ndarray, Z: np.ndarray, u: np.ndarray, v: float):
    """
    Compute Hessian of Gaussian approximation (with think times).
    """
    K, R = L.shape
    Ntot = np.sum(N)
    A = np.zeros((K, K))

    csi = np.zeros(R)
    for r in range(R):
        denom = Z[r] + v * np.dot(u, L[:, r])
        if denom > 0:
            csi[r] = N[r] / denom

    Lhat = np.zeros((K, R))
    for k in range(K):
        for r in range(R):
            Lhat[k, r] = Z[r] + v * L[k, r]

    eta = Ntot + K

    for i in range(K):
        for j in range(K):
            if i != j:
                A[i, j] = -eta * u[i] * u[j]
                for r in range(R):
                    if N[r] > 0:
                        A[i, j] += csi[r]**2 * Lhat[i, r] * Lhat[j, r] * u[i] * u[j] / N[r]

    for i in range(K):
        A[i, i] = -np.sum(A[i, :]) + A[i, i]

    # Reduce to (K-1) x (K-1) and add v column
    A_reduced = A[:K-1, :K-1]
    A_result = np.zeros((K, K))
    A_result[:K-1, :K-1] = A_reduced

    A_result[K-1, K-1] = 1
    for r in range(R):
        if N[r] > 0:
            A_result[K-1, K-1] -= (csi[r]**2 / N[r]) * Z[r] * np.dot(u, L[:, r])
    A_result[K-1, K-1] *= v

    for i in range(K-1):
        A_result[i, K-1] = 0
        for r in range(R):
            if N[r] > 0:
                A_result[i, K-1] += v * u[i] * ((csi[r]**2 / N[r]) * Lhat[i, r] * np.dot(u, L[:, r]) - csi[r] * L[i, r])
        A_result[K-1, i] = A_result[i, K-1]

    return A_result


def _simplex_fun(x: np.ndarray, L: np.ndarray, N: np.ndarray) -> float:
    """
    Evaluate simplex function for LS algorithm.
    """
    M = len(x) + 1
    v = np.zeros(M)
    for i in range(M-1):
        v[i] = np.exp(x[i])
    v[M-1] = 1

    term1 = np.sum(N * np.log(np.dot(v, L)))
    term2 = np.sum(x)
    term3 = -(np.sum(N) + M) * np.log(np.sum(v))

    return np.exp(term1 + term2 + term3)


def _multinomialln(n: np.ndarray) -> float:
    """
    Compute log of multinomial coefficient.
    """
    return _factln(np.sum(n)) - np.sum([_factln(ni) for ni in n])


__all__ = [
    'pfqn_ca',
    'pfqn_nc',
    'pfqn_panacea',
    'pfqn_propfair',
    'pfqn_ls',
]
