"""
Markovian Arrival Process (MAP) and Batch MAP (BMAP) samplers.

This module implements correlated arrival processes:
- MAP: Markov-modulated arrivals (single arrivals)
- BMAP: Batch Markovian Arrival Process (batch arrivals)
- MMPP: Markov-modulated Poisson Process (special case of MAP)
"""

from typing import Optional, Tuple, List
import numpy as np

from .base import DistributionSampler


class MAPSampler(DistributionSampler):
    """
    Markovian Arrival Process sampler.

    A MAP is defined by two matrices:
    - D0: Transition rates without arrivals (n x n)
    - D1: Transition rates with arrivals (n x n)

    The embedded Markov chain has generator Q = D0 + D1.
    Arrivals occur according to D1 transitions.

    The sampler maintains the current phase state between samples.

    Parameters:
        D0: Matrix of transition rates without arrivals
        D1: Matrix of transition rates with arrivals
        initial_phase: Initial phase (0-indexed), or None for stationary
    """

    def __init__(
        self,
        D0: np.ndarray,
        D1: np.ndarray,
        initial_phase: Optional[int] = None,
        rng: Optional[np.random.Generator] = None
    ):
        super().__init__(rng)

        self.D0 = np.asarray(D0, dtype=float)
        self.D1 = np.asarray(D1, dtype=float)

        if self.D0.shape != self.D1.shape:
            raise ValueError("D0 and D1 must have same shape")
        if self.D0.shape[0] != self.D0.shape[1]:
            raise ValueError("D0 must be square")

        self.n_phases = self.D0.shape[0]

        # Compute total rates out of each state
        self.total_rates = -np.diag(self.D0)

        # Precompute transition probabilities for each phase
        self.trans_probs = []
        for i in range(self.n_phases):
            rate = self.total_rates[i]
            if rate > 0:
                probs = np.zeros(2 * self.n_phases)
                # First n_phases: D0 transitions (no arrival)
                for j in range(self.n_phases):
                    if i != j:
                        probs[j] = self.D0[i, j] / rate
                # Next n_phases: D1 transitions (with arrival)
                for j in range(self.n_phases):
                    probs[self.n_phases + j] = self.D1[i, j] / rate
                self.trans_probs.append(np.cumsum(probs))
            else:
                probs = np.zeros(2 * self.n_phases)
                self.trans_probs.append(np.cumsum(probs))

        # Initialize phase
        if initial_phase is not None:
            self.phase = initial_phase
        else:
            # Start from stationary distribution
            self.phase = self._sample_initial_phase()

    def _sample_initial_phase(self) -> int:
        """
        Sample initial phase from stationary distribution.
        """
        Q = self.D0 + self.D1
        # Solve pi * Q = 0 with sum(pi) = 1
        try:
            # Use null space method
            n = self.n_phases
            A = np.vstack([Q.T, np.ones(n)])
            b = np.zeros(n + 1)
            b[-1] = 1.0
            pi, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
            pi = np.maximum(pi, 0)  # Ensure non-negative
            pi = pi / np.sum(pi)  # Normalize

            # Sample from stationary distribution
            cumsum = np.cumsum(pi)
            u = self.rng.random()
            return int(np.searchsorted(cumsum, u))
        except Exception:
            # Fallback to uniform
            return self.rng.integers(0, self.n_phases)

    def sample(self) -> float:
        """
        Sample next inter-arrival time.

        Simulates the MAP until an arrival (D1 transition) occurs.
        Updates internal phase state.
        """
        total_time = 0.0

        while True:
            rate = self.total_rates[self.phase]
            if rate <= 0:
                # Absorbing state - shouldn't happen in valid MAP
                raise RuntimeError(f"Absorbing state {self.phase} in MAP")

            # Time until next transition
            sojourn = self.rng.exponential(1.0 / rate)
            total_time += sojourn

            # Sample transition
            u = self.rng.random()
            next_idx = np.searchsorted(self.trans_probs[self.phase], u)

            if next_idx >= self.n_phases:
                # D1 transition (arrival)
                self.phase = next_idx - self.n_phases
                break
            else:
                # D0 transition (no arrival)
                self.phase = next_idx

        return total_time

    def get_mean(self) -> float:
        """
        Get mean inter-arrival time.

        For MAP: E[X] = pi * (-D0)^{-1} * 1 where pi is stationary
        """
        try:
            D0_inv = np.linalg.inv(-self.D0)
            ones = np.ones(self.n_phases)

            # Get stationary distribution
            Q = self.D0 + self.D1
            n = self.n_phases
            A = np.vstack([Q.T, np.ones(n)])
            b = np.zeros(n + 1)
            b[-1] = 1.0
            pi, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
            pi = np.maximum(pi, 0)
            pi = pi / np.sum(pi)

            return np.dot(pi, np.dot(D0_inv, ones))
        except Exception:
            return float('nan')

    def get_arrival_rate(self) -> float:
        """
        Get the stationary arrival rate (lambda).
        """
        mean = self.get_mean()
        if mean > 0 and not np.isnan(mean):
            return 1.0 / mean
        return float('nan')

    def reset(self, phase: Optional[int] = None):
        """
        Reset the MAP to a specified or random initial phase.
        """
        if phase is not None:
            self.phase = phase
        else:
            self.phase = self._sample_initial_phase()


class MAP2Sampler(DistributionSampler):
    """
    Two-state MAP sampler (optimized).

    Common case of a 2-phase MAP, optimized for efficiency.

    Parameters:
        D0: 2x2 matrix of rates without arrivals
        D1: 2x2 matrix of rates with arrivals
    """

    def __init__(
        self,
        D0: np.ndarray,
        D1: np.ndarray,
        rng: Optional[np.random.Generator] = None
    ):
        super().__init__(rng)

        self.D0 = np.asarray(D0, dtype=float)
        self.D1 = np.asarray(D1, dtype=float)

        if self.D0.shape != (2, 2) or self.D1.shape != (2, 2):
            raise ValueError("D0 and D1 must be 2x2 for MAP2")

        # Extract rates
        self.rate0 = -self.D0[0, 0]  # Total rate out of phase 0
        self.rate1 = -self.D1[1, 1] if self.D1[1, 1] != 0 else -self.D0[1, 1]

        # Compute transition probabilities
        # From phase 0: p00_no_arr (D0[0,1]), p01_no_arr, p0_arr0 (D1[0,0]), p0_arr1 (D1[0,1])
        self.phase = 0

    def sample(self) -> float:
        # Delegate to general MAP implementation
        # For now, use parent class approach but inline for 2x2
        total_time = 0.0
        phase = self.phase

        while True:
            rate = -self.D0[phase, phase]
            sojourn = self.rng.exponential(1.0 / rate)
            total_time += sojourn

            # Compute transition probabilities
            u = self.rng.random() * rate
            cumsum = 0.0

            # Check D0 transitions (no arrival)
            for j in range(2):
                if j != phase:
                    cumsum += self.D0[phase, j]
                    if u < cumsum:
                        phase = j
                        break
            else:
                # Check D1 transitions (arrival)
                for j in range(2):
                    cumsum += self.D1[phase, j]
                    if u < cumsum:
                        self.phase = j
                        return total_time

        return total_time


class BMAPSampler(DistributionSampler):
    """
    Batch Markovian Arrival Process sampler.

    A BMAP generalizes MAP to allow batch arrivals.

    Parameters:
        D0: Matrix of rates without arrivals (n x n)
        D_list: List of matrices D_k for batch size k (k=1,2,...)
        max_batch_size: Maximum batch size to consider
    """

    def __init__(
        self,
        D0: np.ndarray,
        D_list: List[np.ndarray],
        initial_phase: Optional[int] = None,
        rng: Optional[np.random.Generator] = None
    ):
        super().__init__(rng)

        self.D0 = np.asarray(D0, dtype=float)
        self.D_list = [np.asarray(D, dtype=float) for D in D_list]
        self.n_phases = self.D0.shape[0]
        self.max_batch = len(D_list)

        # Validate dimensions
        for k, D in enumerate(self.D_list):
            if D.shape != (self.n_phases, self.n_phases):
                raise ValueError(f"D_{k+1} has wrong shape")

        # Compute total rates
        self.total_rates = -np.diag(self.D0)

        # Precompute transition probabilities
        # Format: (n_phases, n_phases * (1 + max_batch))
        # First n_phases: D0 transitions
        # Next n_phases * k: D_k transitions for batch size k
        self.trans_probs = []
        for i in range(self.n_phases):
            rate = self.total_rates[i]
            n_outcomes = self.n_phases * (1 + self.max_batch)
            probs = np.zeros(n_outcomes)

            if rate > 0:
                # D0 transitions
                for j in range(self.n_phases):
                    if i != j:
                        probs[j] = self.D0[i, j] / rate

                # D_k transitions for k = 1, 2, ...
                offset = self.n_phases
                for k, D_k in enumerate(self.D_list):
                    for j in range(self.n_phases):
                        probs[offset + j] = D_k[i, j] / rate
                    offset += self.n_phases

            self.trans_probs.append(np.cumsum(probs))

        # Initialize phase
        if initial_phase is not None:
            self.phase = initial_phase
        else:
            self.phase = self.rng.integers(0, self.n_phases)

    def sample(self) -> Tuple[float, int]:
        """
        Sample next inter-arrival time and batch size.

        Returns:
            Tuple of (inter-arrival time, batch size)
        """
        total_time = 0.0

        while True:
            rate = self.total_rates[self.phase]
            if rate <= 0:
                raise RuntimeError(f"Absorbing state {self.phase} in BMAP")

            sojourn = self.rng.exponential(1.0 / rate)
            total_time += sojourn

            # Sample transition
            u = self.rng.random()
            next_idx = np.searchsorted(self.trans_probs[self.phase], u)

            if next_idx < self.n_phases:
                # D0 transition (no arrival)
                self.phase = next_idx
            else:
                # D_k transition (batch arrival of size k)
                offset_idx = next_idx - self.n_phases
                batch_size = 1 + (offset_idx // self.n_phases)
                self.phase = offset_idx % self.n_phases
                return (total_time, batch_size)

    def sample_single(self) -> float:
        """
        Sample inter-arrival time (ignoring batch size).

        For compatibility with single-arrival interface.
        """
        time, _ = self.sample()
        return time


class MMPPSampler(DistributionSampler):
    """
    Markov-Modulated Poisson Process sampler.

    MMPP is a special case of MAP where arrivals don't change the phase.
    The arrival rate depends on the current phase.

    Parameters:
        Q: Phase transition rate matrix (n x n)
        lambdas: Arrival rates for each phase (n,)
    """

    def __init__(
        self,
        Q: np.ndarray,
        lambdas: np.ndarray,
        initial_phase: Optional[int] = None,
        rng: Optional[np.random.Generator] = None
    ):
        super().__init__(rng)

        self.Q = np.asarray(Q, dtype=float)
        self.lambdas = np.asarray(lambdas, dtype=float)

        self.n_phases = self.Q.shape[0]
        if len(self.lambdas) != self.n_phases:
            raise ValueError("Number of lambdas must match number of phases")

        # Convert to MAP representation
        # D0[i,j] = Q[i,j] for i != j, D0[i,i] = Q[i,i] - lambda[i]
        # D1[i,i] = lambda[i], D1[i,j] = 0 for i != j
        self.D0 = self.Q.copy()
        for i in range(self.n_phases):
            self.D0[i, i] -= self.lambdas[i]

        self.D1 = np.diag(self.lambdas)

        # Total rates out of each state
        self.total_rates = -np.diag(self.D0)

        # Initialize phase
        if initial_phase is not None:
            self.phase = initial_phase
        else:
            self.phase = self._sample_initial_phase()

    def _sample_initial_phase(self) -> int:
        """Sample from stationary distribution of Q."""
        try:
            n = self.n_phases
            A = np.vstack([self.Q.T, np.ones(n)])
            b = np.zeros(n + 1)
            b[-1] = 1.0
            pi, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
            pi = np.maximum(pi, 0)
            pi = pi / np.sum(pi)
            cumsum = np.cumsum(pi)
            return int(np.searchsorted(cumsum, self.rng.random()))
        except Exception:
            return self.rng.integers(0, self.n_phases)

    def sample(self) -> float:
        """
        Sample next inter-arrival time.
        """
        total_time = 0.0

        while True:
            rate = self.total_rates[self.phase]
            if rate <= 0:
                raise RuntimeError("Invalid MMPP: zero total rate")

            sojourn = self.rng.exponential(1.0 / rate)
            total_time += sojourn

            # Check for arrival
            lambda_i = self.lambdas[self.phase]
            arr_prob = lambda_i / rate

            u = self.rng.random()
            if u < arr_prob:
                # Arrival occurred (phase unchanged for MMPP)
                return total_time
            else:
                # Phase transition (no arrival)
                # Sample next phase from Q transitions
                trans_rate = -self.Q[self.phase, self.phase]
                if trans_rate > 0:
                    u2 = self.rng.random() * trans_rate
                    cumsum = 0.0
                    for j in range(self.n_phases):
                        if j != self.phase:
                            cumsum += self.Q[self.phase, j]
                            if u2 < cumsum:
                                self.phase = j
                                break

    def get_mean(self) -> float:
        """
        Get mean inter-arrival time.

        E[X] = 1 / (pi * lambda) where pi is stationary distribution of Q.
        """
        try:
            n = self.n_phases
            A = np.vstack([self.Q.T, np.ones(n)])
            b = np.zeros(n + 1)
            b[-1] = 1.0
            pi, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
            pi = np.maximum(pi, 0)
            pi = pi / np.sum(pi)
            avg_rate = np.dot(pi, self.lambdas)
            return 1.0 / avg_rate if avg_rate > 0 else float('inf')
        except Exception:
            return float('nan')


def create_map_from_ph(
    ph_alpha: np.ndarray,
    ph_T: np.ndarray,
    rng: Optional[np.random.Generator] = None
) -> MAPSampler:
    """
    Create a MAP from a phase-type renewal process.

    For a PH renewal process, D0 = T and D1 = t * alpha where t = -T*1.

    Parameters:
        ph_alpha: Initial probability vector of PH
        ph_T: Sub-generator matrix of PH

    Returns:
        MAPSampler instance
    """
    alpha = np.asarray(ph_alpha).flatten()
    T = np.asarray(ph_T)
    n = len(alpha)

    # Exit rate vector
    t = -np.sum(T, axis=1)

    # D1 = t * alpha (outer product)
    D1 = np.outer(t, alpha)

    return MAPSampler(T, D1, rng=rng)
