"""
Core Numba JIT kernels for all LINE solvers.

Provides JIT-compiled numerical kernels used across MVA, MAM, DES, and other solvers.
Includes: M/G/c approximations, Erlang-C calculations, distribution sampling,
and matrix operations optimized for JIT compilation.

Graceful fallback to pure Python if Numba is not available.

License: MIT (same as LINE)
"""

import numpy as np
from typing import Tuple, Optional
from ..solver_fld.utils.numba_support import HAS_NUMBA, njit, prange

if HAS_NUMBA:
    # JIT-compiled numerical kernels

    @njit(fastmath=True, cache=True)
    def erlang_c_jit(rho: float, c: int) -> float:
        """
        JIT-compiled Erlang-C formula for M/M/c queues.

        Computes probability that arriving customer must wait.

        Args:
            rho: Traffic intensity (λ/μ)
            c: Number of servers

        Returns:
            Pw: Erlang-C probability [0, 1]
        """
        if c <= 0 or rho <= 0.0:
            return 0.0
        if rho >= float(c):
            return 1.0

        # Erlang-C: Pw = (c*rho)^c / c! / ((c*rho)^c / c! + (1-rho)*sum_{n=0}^{c-1} (c*rho)^n / n!)
        # Numerically stable via iterative computation
        erlang_b = 1.0
        for n in range(1, c):
            erlang_b = (c * rho * erlang_b) / (n + c * rho * erlang_b)

        erlang_c = erlang_b / (erlang_b + (1.0 - rho) * (1.0 - erlang_b))

        if erlang_c < 0.0:
            return 0.0
        if erlang_c > 1.0:
            return 1.0

        return erlang_c

    @njit(fastmath=True, cache=True)
    def erlang_c_wait_jit(rho: float, c: int) -> float:
        """
        JIT-compiled average waiting time for M/M/c queue (in service time units).

        Wq = Pw / (c * (1 - rho))

        Args:
            rho: Traffic intensity
            c: Number of servers

        Returns:
            Wq: Mean waiting time (fraction of service time)
        """
        if rho >= 1.0 or c <= 0:
            return 1e10  # Return large value for unstable system

        pw = erlang_c_jit(rho, c)
        if 1.0 - rho <= 0.0:
            return 1e10

        return pw / (float(c) * (1.0 - rho))

    @njit(fastmath=True, cache=True, parallel=True)
    def mg1_response_time_jit(
        lambda_k: np.ndarray,
        S: np.ndarray,
        cs2: np.ndarray,
        nservers: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        JIT-compiled M/G/1 response time approximation (vectorized).

        For M/G/c: R = S + Wq where Wq depends on SCV and utilization.

        Args:
            lambda_k: Arrival rates by class (K,)
            S: Service time means (M, K)
            cs2: Service time SCV (M, K)
            nservers: Number of servers per station (M,)

        Returns:
            RN: Response time matrix (M, K)
            QN: Queue length matrix (M, K)
        """
        M = len(nservers)
        K = len(lambda_k)

        RN = np.zeros((M, K))
        QN = np.zeros((M, K))

        for m in prange(M):
            c = int(nservers[m])
            total_arrival = 0.0

            # Compute total arrival rate
            for k in range(K):
                total_arrival += lambda_k[k]

            for k in range(K):
                if total_arrival > 0.0 and S[m, k] > 0.0:
                    rho = lambda_k[k] * S[m, k]

                    if rho >= 1.0:
                        RN[m, k] = 1e10
                        QN[m, k] = 1e10
                    else:
                        if c <= 1:
                            # M/G/1: Wq = (cs2 + 1) * rho * S / (2 * (1 - rho))
                            wq = (cs2[m, k] + 1.0) * rho * S[m, k] / (2.0 * (1.0 - rho))
                        else:
                            # M/G/c approximation
                            wq = erlang_c_wait_jit(rho / float(c), c) * S[m, k]

                        RN[m, k] = S[m, k] + wq
                        QN[m, k] = lambda_k[k] * RN[m, k]

        return RN, QN

    @njit(fastmath=True, cache=True)
    def aggregate_scv_jit(
        lambda_k: np.ndarray,
        scv_k: np.ndarray
    ) -> float:
        """
        JIT-compiled aggregated SCV from per-class values.

        ca2 = sum(lambda_k)^2 / sum(lambda_k^2)

        Args:
            lambda_k: Arrival rates (K,)
            scv_k: SCV values (K,)

        Returns:
            ca2: Aggregated SCV
        """
        K = len(lambda_k)

        sum_lambda = 0.0
        sum_lambda_sq = 0.0

        for k in range(K):
            sum_lambda += lambda_k[k]
            sum_lambda_sq += lambda_k[k] ** 2

        if sum_lambda_sq > 0.0:
            return (sum_lambda ** 2) / sum_lambda_sq
        else:
            return 1.0

    @njit(fastmath=True, cache=True)
    def matrix_vector_multiply_jit(
        A: np.ndarray,
        v: np.ndarray
    ) -> np.ndarray:
        """
        JIT-compiled matrix-vector multiply (optimized for small matrices).

        Args:
            A: Matrix (M, K)
            v: Vector (K,)

        Returns:
            result: A @ v (M,)
        """
        M = A.shape[0]
        result = np.zeros(M)

        for i in range(M):
            for j in range(len(v)):
                result[i] += A[i, j] * v[j]

        return result

else:
    # Fallback: pure Python implementations

    def erlang_c_jit(rho: float, c: int) -> float:
        """Pure Python Erlang-C."""
        if c <= 0 or rho <= 0.0:
            return 0.0
        if rho >= float(c):
            return 1.0

        erlang_b = 1.0
        for n in range(1, c):
            erlang_b = (c * rho * erlang_b) / (n + c * rho * erlang_b)

        erlang_c = erlang_b / (erlang_b + (1.0 - rho) * (1.0 - erlang_b))
        return max(0.0, min(1.0, erlang_c))

    def erlang_c_wait_jit(rho: float, c: int) -> float:
        """Pure Python Erlang-C waiting time."""
        if rho >= 1.0 or c <= 0:
            return 1e10

        pw = erlang_c_jit(rho, c)
        if 1.0 - rho <= 0.0:
            return 1e10

        return pw / (float(c) * (1.0 - rho))

    def mg1_response_time_jit(
        lambda_k: np.ndarray,
        S: np.ndarray,
        cs2: np.ndarray,
        nservers: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Pure Python M/G/1 response time."""
        M = len(nservers)
        K = len(lambda_k)

        RN = np.zeros((M, K))
        QN = np.zeros((M, K))

        for m in range(M):
            c = int(nservers[m])
            total_arrival = np.sum(lambda_k)

            for k in range(K):
                if total_arrival > 0.0 and S[m, k] > 0.0:
                    rho = lambda_k[k] * S[m, k]

                    if rho >= 1.0:
                        RN[m, k] = 1e10
                        QN[m, k] = 1e10
                    else:
                        if c <= 1:
                            wq = (cs2[m, k] + 1.0) * rho * S[m, k] / (2.0 * (1.0 - rho))
                        else:
                            wq = erlang_c_wait_jit(rho / float(c), c) * S[m, k]

                        RN[m, k] = S[m, k] + wq
                        QN[m, k] = lambda_k[k] * RN[m, k]

        return RN, QN

    def aggregate_scv_jit(lambda_k: np.ndarray, scv_k: np.ndarray) -> float:
        """Pure Python aggregated SCV."""
        sum_lambda = np.sum(lambda_k)
        sum_lambda_sq = np.sum(lambda_k ** 2)

        if sum_lambda_sq > 0.0:
            return (sum_lambda ** 2) / sum_lambda_sq
        else:
            return 1.0

    def matrix_vector_multiply_jit(A: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Pure Python matrix-vector multiply."""
        return A @ v


__all__ = [
    'HAS_NUMBA',
    'erlang_c_jit',
    'erlang_c_wait_jit',
    'mg1_response_time_jit',
    'aggregate_scv_jit',
    'matrix_vector_multiply_jit',
]
