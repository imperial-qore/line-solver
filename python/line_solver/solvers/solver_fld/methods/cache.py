"""
Cache method for FLD solver - Refined Mean Field analysis of multi-list caches.

Implements mean field and refined mean field (1/N correction) analysis for
multi-list caches with RANDOM(m) replacement policy, based on the DDPP
(Density-Dependent Population Process) framework.

The cache is modeled as a population process where each item can be in one of
h+1 states: list 1, list 2, ..., list h, or outside the cache (list 0).
Transitions occur when an item is requested: the requested item moves up
(toward list 1) while a displaced item moves down (toward list 0).

Reference:
    N. Gast, "Expected Values Estimated via Mean-Field Approximation are
    1/N-Accurate", Proc. ACM Meas. Anal. Comput. Syst., 2017.

Uses rmf_tool library (MIT License) by Nicolas Gast for refined mean field
computations.
"""

import numpy as np
import time as ti
from typing import Optional, Tuple, Dict, List, Any
from scipy.integrate import solve_ivp

from ..options import SolverFLDOptions, FLDResult


class CacheMethod:
    """Refined Mean Field analysis of multi-list caches with RANDOM(m) replacement.

    Models a cache with h lists of capacities m = [m_1, ..., m_h] and n items
    with popularity distribution p = [p_1, ..., p_n]. The state space tracks
    in which list each item resides (or if it is outside the cache).

    State representation:
        x[i + k*n] = probability that item i is in list k
        where k=0 means "outside cache", k=1..h are the cache lists

    Transition dynamics (RANDOM(m) replacement):
        When item i is requested (rate p_i) and item i is in list k:
        - A uniformly random item j from list k+1 is displaced to list k
        - Item i moves from list k to list k+1
        - Net effect: x[i,k] -= 1, x[i,k+1] += 1, x[j,k+1] -= 1, x[j,k] += 1
    """

    def __init__(self, sn, options: SolverFLDOptions):
        self.sn = sn
        self.options = options
        self.iterations = 0
        self.runtime = 0.0

    def solve(self) -> FLDResult:
        """Solve cache model using refined mean field approximation.

        Locates cache nodes in the network, builds a DDPP model for each,
        computes steady-state hit/miss probabilities (with optional 1/N
        refinement), and returns results in the standard FLDResult format.
        """
        start_time = ti.time()

        M = self.sn.nstations
        K = self.sn.nclasses

        # Initialize output arrays
        QN = np.full((M, K), np.nan)
        UN = np.full((M, K), np.nan)
        RN = np.full((M, K), np.nan)
        TN = np.full((M, K), np.nan)
        CN = np.full((1, K), 0.0)
        XN = np.full((1, K), 0.0)

        # Find cache nodes
        cache_results = {}
        for ind in range(self.sn.nnodes):
            if not hasattr(self.sn, 'nodeparam') or self.sn.nodeparam is None:
                continue
            if ind not in self.sn.nodeparam:
                continue

            param = self.sn.nodeparam[ind]
            if not hasattr(param, 'nitems') or param.nitems == 0:
                continue

            # Extract cache parameters from NetworkStruct
            n_items = param.nitems
            itemcap = param.itemcap
            if isinstance(itemcap, (int, float)):
                itemcap = np.array([int(itemcap)])
            else:
                itemcap = np.asarray(itemcap, dtype=int)
            h = len(itemcap)  # number of lists
            m = [int(itemcap[k]) for k in range(h)]

            # Extract per-class popularity distributions
            pread = param.pread  # list of arrays, one per class
            # Build per-class results
            for r in range(K):
                if pread is None or r >= len(pread) or pread[r] is None:
                    continue

                p = np.asarray(pread[r], dtype=float)
                if len(p) != n_items:
                    continue
                # Normalize popularity
                p_sum = np.sum(p)
                if p_sum <= 0:
                    continue
                p = p / p_sum

                # Build and solve the DDPP cache model
                model = CacheRMF(p, m)

                # Compute mean field fixed point
                pi = model.fixed_point()

                # see _kb/06-solver-catalog.md (Fluid: "Fallback/utility
                # numeric methods") for the per-list hit-rate formula
                hit_prob = 0.0
                for k in range(1, h + 1):
                    hit_prob += model.hit_rate(pi, k)
                miss_prob = model.hit_rate(pi, 0)

                # Try refined mean field (1/N correction) if rmf_tool available
                V = None
                try:
                    pi_mf, V_corr, _ = model.meanFieldExpansionSteadyState(order=1)
                    # Refined hit probability: pi + V/N where N = n_items
                    pi_refined = pi_mf + V_corr / n_items
                    # Recompute hit rates with refinement
                    hit_prob_refined = 0.0
                    for k in range(1, h + 1):
                        hit_prob_refined += model.hit_rate(pi_refined, k)
                    miss_prob_refined = model.hit_rate(pi_refined, 0)
                    hit_prob = hit_prob_refined
                    miss_prob = miss_prob_refined
                    V = V_corr
                except Exception:
                    pass  # Fall back to plain mean field

                # Store actual hit/miss probabilities in the cache node param
                if not hasattr(param, 'actualhitprob') or param.actualhitprob is None:
                    param.actualhitprob = np.zeros(K)
                if not hasattr(param, 'actualmissprob') or param.actualmissprob is None:
                    param.actualmissprob = np.zeros(K)
                param.actualhitprob[r] = np.clip(hit_prob, 0.0, 1.0)
                param.actualmissprob[r] = np.clip(miss_prob, 0.0, 1.0)

                cache_results[(ind, r)] = {
                    'pi': pi,
                    'V': V,
                    'hit_prob': hit_prob,
                    'miss_prob': miss_prob,
                    'n_items': n_items,
                    'n_lists': h,
                    'list_capacities': m,
                    'popularity': p,
                }

        elapsed = ti.time() - start_time

        result = FLDResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            CN=CN,
            XN=XN,
            t=None,
            QNt={},
            UNt={},
            TNt={},
            xvec=None,
            iterations=1,
            runtime=elapsed,
            method='rmf'
        )
        # Attach cache-specific results for downstream use
        result.cache_results = cache_results

        return result


def solve_cache(sn, options: Optional[SolverFLDOptions] = None) -> FLDResult:
    """Convenience function to solve cache model.

    Args:
        sn: Compiled NetworkStruct
        options: SolverFLDOptions (uses defaults if None)

    Returns:
        FLDResult with cache hit/miss probabilities
    """
    if options is None:
        options = SolverFLDOptions(method='rmf')
    method = CacheMethod(sn, options)
    return method.solve()


class CacheRMF:
    """Multi-list cache with RANDOM(m) replacement as a DDPP.

    Extends the rmf_tool DDPP framework to model a cache with h lists.
    Each item i can be in list 0 (outside), 1, ..., h (innermost).

    State vector layout:
        x[i + k * n_items] = density of item i in list k
        for i in 0..n_items-1, k in 0..n_lists

    Parameters
    ----------
    p : np.ndarray, shape (n_items,)
        Item request probabilities (popularity distribution), normalized.
    m : list of int, length h
        Capacity of each cache list. m[0] is the outermost list,
        m[h-1] is the innermost.

    Attributes
    ----------
    number_of_items : int
        Number of items in the catalog (n).
    number_of_lists : int
        Number of cache lists (h).
    model_dimension : int
        Total state dimension: n * (h + 1).
    """

    def __init__(self, p: np.ndarray, m: list):
        self.number_of_items = len(p)
        self.number_of_lists = len(m)
        self.model_dimension = self.number_of_items * (self.number_of_lists + 1)
        self.m = m
        self.p = p
        self._ddpp = None  # Lazy-initialized rmf_tool DDPP

        n = self.number_of_items
        h = self.number_of_lists

        # Build initial state: first m[0] items in list 1, next m[1] in list 2, etc.
        # Remaining items outside cache (list 0)
        self._x0 = np.zeros(self.model_dimension)
        obj_idx = 0
        for k in range(h):
            for _ in range(m[k]):
                if obj_idx < n:
                    self._x0[self.index(obj_idx, k + 1)] = 1.0
                    obj_idx += 1
        for i in range(obj_idx, n):
            self._x0[self.index(i, 0)] = 1.0

    def index(self, i: int, k: int) -> int:
        """Map (item i, list k) to flat state index."""
        return i + k * self.number_of_items

    def hit_rate(self, x: np.ndarray, list_number: int) -> float:
        """Compute hit rate contribution from a specific list.

        Args:
            x: State vector of dimension model_dimension.
            list_number: List index (0 = outside cache, 1..h = cache lists).

        Returns:
            Sum of p_i * x[i, list_number] over all items i.
        """
        n = self.number_of_items
        return np.sum([
            self.p[i] * x[self.index(i, list_number)]
            for i in range(n)
        ])

    def drift(self, x: np.ndarray) -> np.ndarray:
        """Compute the mean field drift F(x).

        The drift for RANDOM(m) replacement is:
            dx[i,k]/dt = -p_i * x[i,k] + hitRate[k] * x[i,k+1] / m[k]
            dx[i,k+1]/dt = p_i * x[i,k] - hitRate[k] * x[i,k+1] / m[k]

        where hitRate[k] = sum_j p_j * x[j,k] is the aggregate request rate
        for items currently in list k.

        Args:
            x: State vector of dimension model_dimension.

        Returns:
            Drift vector dx/dt of same dimension.
        """
        n = self.number_of_items
        h = self.number_of_lists
        hit_rates = [self.hit_rate(x, k) for k in range(h + 1)]
        dX = np.zeros(self.model_dimension)
        for i in range(n):
            for k in range(h):
                flow = self.p[i] * x[self.index(i, k)] - hit_rates[k] * x[self.index(i, k + 1)] / self.m[k]
                dX[self.index(i, k)] -= flow
                dX[self.index(i, k + 1)] += flow
        return dX

    def jacobian(self, x: np.ndarray) -> np.ndarray:
        """Compute Jacobian dF/dx at state x.

        Args:
            x: State vector of dimension model_dimension.

        Returns:
            Jacobian matrix of shape (model_dimension, model_dimension).
        """
        n = self.number_of_items
        h = self.number_of_lists
        dim = self.model_dimension
        hit_rates = [self.hit_rate(x, k) for k in range(h + 1)]
        Fp = np.zeros((dim, dim))

        for i in range(n):
            for k in range(h):
                ik = self.index(i, k)
                ik1 = self.index(i, k + 1)

                # Direct rate terms
                Fp[ik, ik] -= self.p[i]
                Fp[ik1, ik] += self.p[i]
                Fp[ik, ik1] += hit_rates[k] / self.m[k]
                Fp[ik1, ik1] -= hit_rates[k] / self.m[k]

                # Indirect terms via hit rate dependence on x[j,k]
                for j in range(n):
                    jk = self.index(j, k)
                    jk1 = self.index(j, k + 1)
                    # d(hitRate[k])/d(x[j,k]) = p[j], affects x[i,k+1] terms
                    Fp[ik, jk1] -= self.p[i] * x[ik] / self.m[k]
                    Fp[ik1, jk1] += self.p[i] * x[ik] / self.m[k]
                    Fp[ik, jk] += self.p[j] * x[ik1] / self.m[k]
                    Fp[ik1, jk] -= self.p[j] * x[ik1] / self.m[k]
        return Fp

    def hessian(self, x: np.ndarray) -> np.ndarray:
        """Compute Hessian d^2F/dx^2 at state x.

        The Hessian is constant (drift is quadratic in x), so the x argument
        is unused but kept for interface consistency.

        Args:
            x: State vector (unused — Hessian is constant for this model).

        Returns:
            Hessian tensor of shape (model_dimension, model_dimension, model_dimension).
        """
        n = self.number_of_items
        h = self.number_of_lists
        dim = self.model_dimension
        Fpp = np.zeros((dim, dim, dim))

        for i in range(n):
            for k in range(h):
                ik = self.index(i, k)
                ik1 = self.index(i, k + 1)
                for j in range(n):
                    if j != i:
                        jk = self.index(j, k)
                        jk1 = self.index(j, k + 1)
                        # d^2 F[ik] / (d x[jk] d x[ik1]) = p[j]/m[k]
                        Fpp[ik, jk, ik1] += self.p[j] / self.m[k]
                        Fpp[ik, ik1, jk] += self.p[j] / self.m[k]
                        # d^2 F[ik] / (d x[jk1] d x[ik]) = -p[i]/m[k]
                        Fpp[ik, jk1, ik] += -self.p[i] / self.m[k]
                        Fpp[ik, ik, jk1] += -self.p[i] / self.m[k]
                        # Symmetric for ik1
                        Fpp[ik1, jk, ik1] -= self.p[j] / self.m[k]
                        Fpp[ik1, ik1, jk] -= self.p[j] / self.m[k]
                        Fpp[ik1, jk1, ik] -= -self.p[i] / self.m[k]
                        Fpp[ik1, ik, jk1] -= -self.p[i] / self.m[k]
        return Fpp

    def noise_matrix(self, x: np.ndarray) -> np.ndarray:
        """Compute noise intensity matrix Q(x) for the DDPP.

        Q[a,b] = sum_ell ell[a] * ell[b] * beta_ell(x)

        where each transition ell is a swap between items i and j across
        lists k and k+1, with rate p_i * x[i,k] * x[j,k+1] / m[k].

        Args:
            x: State vector of dimension model_dimension.

        Returns:
            Noise matrix of shape (model_dimension, model_dimension).
        """
        n = self.number_of_items
        h = self.number_of_lists
        dim = self.model_dimension
        Q = np.zeros((dim, dim))

        for i in range(n):
            for k in range(h):
                for j in range(n):
                    rate = self.p[i] * x[self.index(i, k)] * x[self.index(j, k + 1)] / self.m[k]
                    indices = [
                        self.index(i, k),
                        self.index(j, k),
                        self.index(i, k + 1),
                        self.index(j, k + 1)
                    ]
                    signs = [-1, 1, 1, -1]
                    for ia, a in enumerate(indices):
                        for ib, b in enumerate(indices):
                            Q[a, b] += rate * signs[ia] * signs[ib]
        return Q

    def fixed_point(self, tmax: float = 10000.0) -> np.ndarray:
        """Compute mean field fixed point by ODE integration.

        Integrates dx/dt = F(x) until steady state.

        Args:
            tmax: Maximum integration time (should be large enough for convergence).

        Returns:
            Fixed point state vector of dimension model_dimension.
        """
        sol = solve_ivp(
            lambda t, x: self.drift(x),
            [0, tmax],
            self._x0,
            method='LSODA',
            rtol=1e-8,
            atol=1e-10
        )
        return sol.y[:, -1]

    def _dimension_reduction(self, Fp: np.ndarray):
        """Compute dimension reduction matrices for the singular Jacobian.

        The Jacobian is singular because item populations are conserved
        (sum over lists for each item = 1). This method finds a change of
        basis that separates the rank-deficient directions.

        Args:
            Fp: Jacobian matrix at fixed point.

        Returns:
            Tuple (P, Pinv, rank) where P @ Fp @ Pinv is block-diagonal
            with the top-left (rank x rank) block being non-singular.
        """
        import scipy.linalg as la

        dim = self.model_dimension
        n = self.number_of_items
        h = self.number_of_lists

        rank = np.linalg.matrix_rank(Fp)

        # Build change-of-basis: first rank rows are independent coordinates,
        # remaining rows span the null space of Fp
        C = np.zeros((dim, dim))
        d = 0
        for l_idx in range(h + 1):
            for i in range(n - 1):
                C[d, self.index(i, l_idx)] = 1.0
                d += 1

        U, s, Vh = la.svd(Fp)
        C[rank:, :] = U.T[rank:, :]
        Cinv = la.inv(C)
        return C, Cinv, rank

    def _reduce_FpFppQ(self, Fp, Fpp, Q):
        """Apply dimension reduction to Fp, Fpp, Q.

        Projects the Jacobian, Hessian, and noise matrix onto the
        non-singular subspace identified by dimension reduction.

        Returns:
            Tuple (Fp_r, Fpp_r, Q_r, P, Pinv, rank).
        """
        P, Pinv, rank = self._dimension_reduction(Fp)
        Fp_r = (P @ Fp @ Pinv)[:rank, :rank]
        Fpp_r = np.tensordot(
            np.tensordot(
                np.tensordot(P, Fpp, axes=([1], [0])),
                Pinv, axes=([1], [0])
            ),
            Pinv, axes=([2], [0])
        )[:rank, :rank, :rank]
        Q_r = (P @ Q @ P.T)[:rank, :rank]
        return Fp_r, Fpp_r, Q_r, P, Pinv, rank

    def _expand_VW(self, V_r, W_r, Pinv, rank):
        """Expand reduced V, W back to full dimension."""
        V = Pinv[:, :rank] @ V_r
        W = Pinv[:, :rank] @ W_r @ Pinv.T[:rank, :]
        return V, W

    def meanFieldExpansionSteadyState(self, order: int = 1):
        """Compute refined mean field steady-state expansion.

        Computes the mean field fixed point pi and the 1/N correction V
        using the Lyapunov equation approach with dimension reduction to
        handle the singular Jacobian (due to per-item conservation constraints).

        The refined approximation for a system of N items is:
            E[X] ≈ pi + V/N + O(1/N^2)

        Args:
            order: Expansion order (0 = plain mean field, 1 = with 1/N correction).

        Returns:
            Tuple of (pi, V, (V, W)) where:
            - pi: Mean field fixed point (model_dimension,)
            - V: First-order correction coefficient (model_dimension,)
            - W: Covariance at fixed point (model_dimension, model_dimension)
        """
        from scipy.linalg import solve_lyapunov

        pi = self.fixed_point()

        if order == 0:
            V = np.zeros_like(pi)
            W = np.zeros((self.model_dimension, self.model_dimension))
            return pi, V, (V, W)

        Fp = self.jacobian(pi)
        Fpp = self.hessian(pi)
        Q = self.noise_matrix(pi)

        # Dimension reduction: project onto non-singular subspace
        Fp_r, Fpp_r, Q_r, P, Pinv, rank = self._reduce_FpFppQ(Fp, Fpp, Q)

        # Solve Lyapunov equation in reduced space: Fp_r @ W_r + W_r @ Fp_r.T + Q_r = 0
        W_r = solve_lyapunov(Fp_r, -Q_r)

        # First-order correction in reduced space
        C_r = np.tensordot(Fpp_r, W_r, axes=([1, 2], [0, 1]))
        V_r = -np.linalg.solve(Fp_r, C_r / 2.0)

        # Expand back to full dimension
        V, W = self._expand_VW(V_r, W_r, Pinv, rank)

        return pi, V, (V, W)

    def meanFieldExpansionTransient(self, time: float = 50.0, n_points: int = 200, order: int = 1):
        """Compute refined mean field transient expansion.

        Integrates the coupled ODE system for (X, V, W) where:
        - X(t): mean field trajectory
        - V(t): 1/N correction trajectory
        - W(t): covariance trajectory

        Args:
            time: Maximum integration time.
            n_points: Number of output time points.
            order: Expansion order (0 or 1).

        Returns:
            Tuple of (T, X, V, W) where:
            - T: Time points array (n_points,)
            - X: Mean field trajectory (n_points, model_dimension)
            - V: Correction trajectory (n_points, model_dimension)
            - W: Covariance trajectory (n_points, model_dimension, model_dimension)
        """
        from scipy.linalg import solve_lyapunov

        dim = self.model_dimension
        T = np.linspace(0, time, n_points)

        if order == 0:
            sol = solve_ivp(
                lambda t, x: self.drift(x),
                [0, time],
                self._x0,
                t_eval=T,
                method='LSODA',
                rtol=1e-6
            )
            X = sol.y.T
            V = np.zeros_like(X)
            W = np.zeros((n_points, dim, dim))
            return T, X, V, W

        # Coupled ODE: state = [X (dim), V (dim), W_flat (dim^2)]
        total_dim = dim + dim + dim * dim
        y0 = np.zeros(total_dim)
        y0[:dim] = self._x0

        def coupled_rhs(t, y):
            x = y[:dim]
            v = y[dim:2 * dim]
            w = y[2 * dim:].reshape(dim, dim)

            F = self.drift(x)
            Fp = self.jacobian(x)
            Fpp = self.hessian(x)
            Qmat = self.noise_matrix(x)

            dx = F
            dv = Fp @ v + 0.5 * np.tensordot(Fpp, w, axes=([1, 2], [0, 1]))
            dw = Fp @ w + w @ Fp.T + Qmat

            return np.concatenate([dx, dv, dw.ravel()])

        sol = solve_ivp(
            coupled_rhs,
            [0, time],
            y0,
            t_eval=T,
            method='LSODA',
            rtol=1e-6
        )

        X = sol.y[:dim, :].T
        V = sol.y[dim:2 * dim, :].T
        W_flat = sol.y[2 * dim:, :].T
        W = W_flat.reshape(n_points, dim, dim)

        return T, X, V, W

    def hit_rates_all(self, x: np.ndarray) -> np.ndarray:
        """Compute hit rates for all lists.

        Args:
            x: State vector of dimension model_dimension.

        Returns:
            Array of shape (number_of_lists + 1,) with hit rate per list.
        """
        return np.array([
            self.hit_rate(x, k) for k in range(self.number_of_lists + 1)
        ])

    def convertTo2D(self, x: np.ndarray) -> np.ndarray:
        """Convert flat state vector to 2D (items x lists) matrix.

        Args:
            x: State vector of dimension model_dimension, or
               matrix of shape (T, model_dimension) for trajectories.

        Returns:
            Array of shape (n_items, n_lists+1) or (T, n_items, n_lists+1).
        """
        n = self.number_of_items
        h = self.number_of_lists
        if x.ndim == 1:
            out = np.zeros((n, h + 1))
            for i in range(n):
                for k in range(h + 1):
                    out[i, k] = x[self.index(i, k)]
            return out
        else:
            T = x.shape[0]
            out = np.zeros((T, n, h + 1))
            for t in range(T):
                for i in range(n):
                    for k in range(h + 1):
                        out[t, i, k] = x[t, self.index(i, k)]
            return out
