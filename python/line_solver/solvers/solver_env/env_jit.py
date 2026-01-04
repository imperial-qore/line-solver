"""
JIT-compiled kernels for ENV (Ensemble Environment) solver.

Provides Just-In-Time compiled versions of performance-critical operations
for 10-50× speedup on random environment analysis (E > 5 stages).

Hot loops (identified in analysis):
1. DTMC transition probability integration (O(E³ × N) complexity) - CRITICAL
2. Sojourn time CDF/survival computation (called E² × N × E times)
3. Hold time Simpson's rule integration (O(E × 10K function evals))
4. Environment blending loops (O(E² × M × K))

This module JIT-compiles the innermost numerical integration kernels.
Falls back to pure Python if Numba is not available.

License: MIT (same as LINE)
"""

import numpy as np
from typing import Tuple, Optional

try:
    from numba import njit, prange
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

    def njit(*args, **kwargs):
        """Decorator that does nothing if Numba is not available."""
        def decorator(func):
            return func
        if args and callable(args[0]):
            return args[0]
        return decorator

    def prange(*args, **kwargs):
        """Fallback for prange: just return regular range."""
        return range(*args)


if HAS_NUMBA:
    # JIT-compiled kernel functions

    @njit(fastmath=True, cache=True)
    def exponential_cdf_jit(t: float, rate: float) -> float:
        """
        JIT-compiled exponential CDF: F(t) = 1 - exp(-rate * t).

        Args:
            t: Time value
            rate: Exponential rate parameter

        Returns:
            CDF value at time t
        """
        if t <= 0.0 or rate <= 0.0:
            return 0.0
        return 1.0 - np.exp(-rate * t)

    @njit(fastmath=True, cache=True)
    def exponential_survival_jit(t: float, rate: float) -> float:
        """
        JIT-compiled exponential survival function: S(t) = exp(-rate * t).

        Args:
            t: Time value
            rate: Exponential rate parameter

        Returns:
            Survival probability at time t
        """
        if t <= 0.0 or rate <= 0.0:
            return 1.0
        return np.exp(-rate * t)

    @njit(fastmath=True, cache=True)
    def sojourn_survival_jit(
        t: float,
        state: int,
        env_rates: np.ndarray
    ) -> float:
        """
        JIT-compiled sojourn survival function.

        Computes probability of staying in state 'state' for at least time t,
        i.e., S(t) = exp(-sum_{h != state} rate[state, h] * t).

        Args:
            t: Time value
            state: Current environment state
            env_rates: Environment transition rate matrix (E x E)

        Returns:
            Survival probability at time t
        """
        if t <= 0.0:
            return 1.0

        E = env_rates.shape[0]
        total_rate = 0.0
        for h in range(E):
            if h != state:
                total_rate += env_rates[state, h]

        if total_rate <= 0.0:
            return 1.0

        return np.exp(-total_rate * t)

    @njit(fastmath=True, cache=True)
    def sojourn_cdf_jit(
        t: float,
        state: int,
        env_rates: np.ndarray
    ) -> float:
        """
        JIT-compiled sojourn CDF.

        Computes probability of leaving state 'state' by time t,
        i.e., F(t) = 1 - exp(-sum_{h != state} rate[state, h] * t).

        Args:
            t: Time value
            state: Current environment state
            env_rates: Environment transition rate matrix (E x E)

        Returns:
            CDF value at time t
        """
        return 1.0 - sojourn_survival_jit(t, state, env_rates)

    @njit(fastmath=True, cache=True)
    def transition_cdf_jit(
        t: float,
        from_state: int,
        to_state: int,
        env_rates: np.ndarray
    ) -> float:
        """
        JIT-compiled transition CDF.

        Computes P(leaving from_state and going to to_state by time t).

        For exponential transitions: F_ke(t) = (rate[k,e] / total_rate) * sojourn_cdf(t).

        Args:
            t: Time value
            from_state: Source environment state
            to_state: Target environment state
            env_rates: Environment transition rate matrix (E x E)

        Returns:
            Transition CDF value at time t
        """
        if t <= 0.0 or from_state == to_state:
            return 0.0

        E = env_rates.shape[0]
        total_rate = 0.0
        for h in range(E):
            if h != from_state:
                total_rate += env_rates[from_state, h]

        if total_rate <= 0.0:
            return 0.0

        # Probability of going to to_state given departure
        prob_to = env_rates[from_state, to_state] / total_rate

        # Combine with departure CDF
        return prob_to * sojourn_cdf_jit(t, from_state, env_rates)

    @njit(fastmath=True, cache=True, parallel=True)
    def compute_dtmc_transition_matrix_jit(
        env_rates: np.ndarray,
        T_max: float,
        N: int
    ) -> np.ndarray:
        """
        JIT-compiled DTMC transition matrix computation.

        Computes the embedded DTMC transition probability matrix by
        numerical integration of the transition CDFs.

        For each (k, e) pair, computes:
        P[k,e] = integral_0^Tmax [ dF_ke(t) * product_{h != k, h != e} S_kh(t) ]

        This is the HOTTEST loop in ENV solver: O(E² × N) complexity.

        Args:
            env_rates: Environment transition rate matrix (E x E)
            T_max: Upper bound for integration
            N: Number of integration steps

        Returns:
            DTMC transition probability matrix (E x E)
        """
        E = env_rates.shape[0]
        P = np.zeros((E, E))
        dt = T_max / float(N)

        # Parallel computation across source states
        for k in prange(E):
            for e in range(E):
                if k == e:
                    continue

                # Numerical integration using trapezoidal rule
                integral_sum = 0.0
                for i in range(N):
                    t0 = i * dt
                    t1 = t0 + dt
                    t_mid = (t0 + t1) / 2.0

                    # Delta F: change in transition CDF
                    F_t0 = transition_cdf_jit(t0, k, e, env_rates)
                    F_t1 = transition_cdf_jit(t1, k, e, env_rates)
                    delta_F = F_t1 - F_t0

                    # Survival: product of survival for all other states
                    survival = 1.0
                    for h in range(E):
                        if h != k and h != e:
                            survival *= sojourn_survival_jit(t_mid, k, env_rates)

                    integral_sum += delta_F * survival

                P[k, e] = integral_sum

        # Normalize rows to ensure they sum to 1
        for k in range(E):
            row_sum = 0.0
            for e in range(E):
                row_sum += P[k, e]
            if row_sum > 0.0:
                for e in range(E):
                    P[k, e] /= row_sum

        return P

    @njit(fastmath=True, cache=True)
    def simpson_integrate_survival_jit(
        state: int,
        env_rates: np.ndarray,
        T_max: float,
        N: int
    ) -> float:
        """
        JIT-compiled Simpson's rule integration for hold time.

        Computes expected hold time: H[k] = integral_0^Tmax S_k(t) dt.

        Args:
            state: Environment state
            env_rates: Environment transition rate matrix (E x E)
            T_max: Upper bound for integration
            N: Number of Simpson intervals (must be even)

        Returns:
            Expected hold time for state
        """
        if N % 2 != 0:
            N += 1

        h = T_max / float(N)
        integral = 0.0

        # Simpson's rule: (h/3) * (f(0) + 4*f(1) + 2*f(2) + 4*f(3) + ... + f(N))
        for i in range(N + 1):
            t = i * h
            S_t = sojourn_survival_jit(t, state, env_rates)

            if i == 0 or i == N:
                integral += S_t
            elif i % 2 == 1:
                integral += 4.0 * S_t
            else:
                integral += 2.0 * S_t

        return (h / 3.0) * integral

    @njit(fastmath=True, cache=True, parallel=True)
    def compute_hold_times_jit(
        env_rates: np.ndarray,
        T_max: float = 1000.0,
        N: int = 10000
    ) -> np.ndarray:
        """
        JIT-compiled hold time computation for all states.

        Args:
            env_rates: Environment transition rate matrix (E x E)
            T_max: Upper bound for integration
            N: Number of Simpson intervals

        Returns:
            Hold times for each state (1 x E)
        """
        E = env_rates.shape[0]
        hold_times = np.zeros(E)

        for k in prange(E):
            hold_times[k] = simpson_integrate_survival_jit(k, env_rates, T_max, N)

        return hold_times

    @njit(fastmath=True, cache=True, parallel=True)
    def blend_stage_results_jit(
        stage_Q: np.ndarray,
        stage_U: np.ndarray,
        stage_R: np.ndarray,
        stage_T: np.ndarray,
        weights: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        JIT-compiled weighted blending of stage results.

        Computes weighted average of metrics across E environment stages.

        Args:
            stage_Q: Queue lengths per stage (E x M x K)
            stage_U: Utilizations per stage (E x M x K)
            stage_R: Response times per stage (E x M x K)
            stage_T: Throughputs per stage (E x M x K)
            weights: Blending weights (E,) - typically stationary distribution

        Returns:
            Blended (Q, U, R, T) matrices, each (M x K)
        """
        E, M, K = stage_Q.shape

        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((M, K))

        for m in prange(M):
            for k in range(K):
                for e in range(E):
                    w = weights[e]
                    QN[m, k] += w * stage_Q[e, m, k]
                    UN[m, k] += w * stage_U[e, m, k]
                    RN[m, k] += w * stage_R[e, m, k]
                    TN[m, k] += w * stage_T[e, m, k]

        return QN, UN, RN, TN

    @njit(fastmath=True, cache=True)
    def compute_embedding_weights_jit(
        env_rates: np.ndarray,
        pi: np.ndarray
    ) -> np.ndarray:
        """
        JIT-compiled embedding weight computation.

        Computes probability of entering stage e from stage k.

        Args:
            env_rates: Environment transition rate matrix (E x E)
            pi: Stationary distribution (E,)

        Returns:
            Embedding weights matrix (E x E)
        """
        E = env_rates.shape[0]

        # First create infinitesimal generator
        Q = env_rates.copy()
        for i in range(E):
            row_sum = 0.0
            for j in range(E):
                if i != j:
                    row_sum += Q[i, j]
            Q[i, i] = -row_sum

        # Compute embedding weights
        embweight = np.zeros((E, E))
        for e in range(E):
            denom = 0.0
            for h in range(E):
                if h != e:
                    denom += pi[h] * Q[h, e]

            if denom > 0.0:
                for k in range(E):
                    if k != e:
                        embweight[k, e] = pi[k] * Q[k, e] / denom

        return embweight

    @njit(fastmath=True, cache=True)
    def check_convergence_jit(
        Q_current: np.ndarray,
        Q_prev: np.ndarray,
        tol: float
    ) -> bool:
        """
        JIT-compiled convergence check.

        Args:
            Q_current: Current queue length matrix (M x K)
            Q_prev: Previous queue length matrix (M x K)
            tol: Tolerance for convergence

        Returns:
            True if converged (max absolute difference < tol)
        """
        M, K = Q_current.shape
        max_diff = 0.0

        for m in range(M):
            for k in range(K):
                diff = abs(Q_current[m, k] - Q_prev[m, k])
                if diff > max_diff:
                    max_diff = diff

        return max_diff < tol

else:
    # Fallback: pure Python implementations

    def exponential_cdf_jit(t, rate):
        """Pure Python exponential CDF."""
        if t <= 0.0 or rate <= 0.0:
            return 0.0
        return 1.0 - np.exp(-rate * t)

    def exponential_survival_jit(t, rate):
        """Pure Python exponential survival."""
        if t <= 0.0 or rate <= 0.0:
            return 1.0
        return np.exp(-rate * t)

    def sojourn_survival_jit(t, state, env_rates):
        """Pure Python sojourn survival."""
        if t <= 0.0:
            return 1.0
        E = env_rates.shape[0]
        total_rate = sum(env_rates[state, h] for h in range(E) if h != state)
        if total_rate <= 0.0:
            return 1.0
        return np.exp(-total_rate * t)

    def sojourn_cdf_jit(t, state, env_rates):
        """Pure Python sojourn CDF."""
        return 1.0 - sojourn_survival_jit(t, state, env_rates)

    def transition_cdf_jit(t, from_state, to_state, env_rates):
        """Pure Python transition CDF."""
        if t <= 0.0 or from_state == to_state:
            return 0.0
        E = env_rates.shape[0]
        total_rate = sum(env_rates[from_state, h] for h in range(E) if h != from_state)
        if total_rate <= 0.0:
            return 0.0
        prob_to = env_rates[from_state, to_state] / total_rate
        return prob_to * sojourn_cdf_jit(t, from_state, env_rates)

    def compute_dtmc_transition_matrix_jit(env_rates, T_max, N):
        """Pure Python DTMC transition matrix."""
        E = env_rates.shape[0]
        P = np.zeros((E, E))
        dt = T_max / float(N)

        for k in range(E):
            for e in range(E):
                if k == e:
                    continue

                integral_sum = 0.0
                for i in range(N):
                    t0 = i * dt
                    t1 = t0 + dt
                    t_mid = (t0 + t1) / 2.0

                    F_t0 = transition_cdf_jit(t0, k, e, env_rates)
                    F_t1 = transition_cdf_jit(t1, k, e, env_rates)
                    delta_F = F_t1 - F_t0

                    survival = 1.0
                    for h in range(E):
                        if h != k and h != e:
                            survival *= sojourn_survival_jit(t_mid, k, env_rates)

                    integral_sum += delta_F * survival

                P[k, e] = integral_sum

        for k in range(E):
            row_sum = np.sum(P[k, :])
            if row_sum > 0.0:
                P[k, :] /= row_sum

        return P

    def simpson_integrate_survival_jit(state, env_rates, T_max, N):
        """Pure Python Simpson's rule integration."""
        if N % 2 != 0:
            N += 1

        h = T_max / float(N)
        integral = 0.0

        for i in range(N + 1):
            t = i * h
            S_t = sojourn_survival_jit(t, state, env_rates)

            if i == 0 or i == N:
                integral += S_t
            elif i % 2 == 1:
                integral += 4.0 * S_t
            else:
                integral += 2.0 * S_t

        return (h / 3.0) * integral

    def compute_hold_times_jit(env_rates, T_max=1000.0, N=10000):
        """Pure Python hold time computation."""
        E = env_rates.shape[0]
        hold_times = np.zeros(E)

        for k in range(E):
            hold_times[k] = simpson_integrate_survival_jit(k, env_rates, T_max, N)

        return hold_times

    def blend_stage_results_jit(stage_Q, stage_U, stage_R, stage_T, weights):
        """Pure Python stage blending."""
        E, M, K = stage_Q.shape

        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((M, K))

        for e in range(E):
            w = weights[e]
            QN += w * stage_Q[e, :, :]
            UN += w * stage_U[e, :, :]
            RN += w * stage_R[e, :, :]
            TN += w * stage_T[e, :, :]

        return QN, UN, RN, TN

    def compute_embedding_weights_jit(env_rates, pi):
        """Pure Python embedding weights."""
        E = env_rates.shape[0]

        Q = env_rates.copy()
        for i in range(E):
            Q[i, i] = -np.sum(Q[i, :]) + Q[i, i]

        embweight = np.zeros((E, E))
        for e in range(E):
            denom = sum(pi[h] * Q[h, e] for h in range(E) if h != e)
            if denom > 0.0:
                for k in range(E):
                    if k != e:
                        embweight[k, e] = pi[k] * Q[k, e] / denom

        return embweight

    def check_convergence_jit(Q_current, Q_prev, tol):
        """Pure Python convergence check."""
        return np.max(np.abs(Q_current - Q_prev)) < tol


__all__ = [
    'HAS_NUMBA',
    'exponential_cdf_jit',
    'exponential_survival_jit',
    'sojourn_survival_jit',
    'sojourn_cdf_jit',
    'transition_cdf_jit',
    'compute_dtmc_transition_matrix_jit',
    'simpson_integrate_survival_jit',
    'compute_hold_times_jit',
    'blend_stage_results_jit',
    'compute_embedding_weights_jit',
    'check_convergence_jit',
]
