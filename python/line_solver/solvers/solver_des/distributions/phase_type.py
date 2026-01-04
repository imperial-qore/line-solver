"""
Phase-type distribution samplers.

This module implements matrix-exponential distributions:
- PH: General phase-type distribution
- APH: Acyclic phase-type distribution
- Coxian: Coxian distribution (series of exponentials with escape)
- HyperExponential: Mixture of exponentials
"""

from typing import Optional, List, Tuple
import numpy as np

from .base import DistributionSampler


class PHSampler(DistributionSampler):
    """
    General Phase-Type distribution sampler.

    A phase-type distribution is defined by:
    - alpha: Initial probability vector (1 x n)
    - T: Sub-generator matrix (n x n), negative diagonal
    - t: Exit rate vector t = -T * 1 (n x 1)

    The distribution represents the absorption time of a Markov chain.

    Parameters:
        alpha: Initial phase probabilities (row vector)
        T: Sub-generator matrix
    """

    def __init__(
        self,
        alpha: np.ndarray,
        T: np.ndarray,
        rng: Optional[np.random.Generator] = None
    ):
        super().__init__(rng)
        self.alpha = np.asarray(alpha).flatten()
        self.T = np.asarray(T)

        if len(self.alpha) != self.T.shape[0]:
            raise ValueError("Alpha length must match T dimension")
        if self.T.shape[0] != self.T.shape[1]:
            raise ValueError("T must be square")

        self.n_phases = len(self.alpha)
        # Exit rates: t = -T * ones
        self.t = -np.sum(self.T, axis=1)
        # Total rates out of each phase
        self.total_rates = -np.diag(self.T)

        # Precompute cumulative probabilities for initial state selection
        self.alpha_cumsum = np.cumsum(self.alpha)

        # Precompute transition probabilities for each phase
        self.trans_probs = []
        for i in range(self.n_phases):
            if self.total_rates[i] > 0:
                # Probability of transition to each other phase
                probs = np.zeros(self.n_phases + 1)  # +1 for absorption
                for j in range(self.n_phases):
                    if i != j:
                        probs[j] = self.T[i, j] / self.total_rates[i]
                probs[self.n_phases] = self.t[i] / self.total_rates[i]  # Absorption
                self.trans_probs.append(np.cumsum(probs))
            else:
                # Absorbing state
                probs = np.zeros(self.n_phases + 1)
                probs[self.n_phases] = 1.0
                self.trans_probs.append(np.cumsum(probs))

    def sample(self) -> float:
        """
        Sample from the PH distribution using phase simulation.

        Simulates the Markov chain until absorption.
        """
        total_time = 0.0

        # Select initial phase
        u = self.rng.random()
        phase = np.searchsorted(self.alpha_cumsum, u)
        phase = min(phase, self.n_phases - 1)

        # Simulate until absorption
        while True:
            rate = self.total_rates[phase]
            if rate <= 0:
                break

            # Time in current phase
            sojourn = self.rng.exponential(1.0 / rate)
            total_time += sojourn

            # Transition
            u = self.rng.random()
            next_state = np.searchsorted(self.trans_probs[phase], u)

            if next_state >= self.n_phases:
                # Absorbed
                break

            phase = next_state

        return total_time

    def get_mean(self) -> float:
        """
        Compute mean of PH distribution.

        E[X] = -alpha * T^{-1} * 1
        """
        try:
            T_inv = np.linalg.inv(self.T)
            ones = np.ones(self.n_phases)
            return -np.dot(self.alpha, np.dot(T_inv, ones))
        except np.linalg.LinAlgError:
            return float('nan')

    def get_variance(self) -> float:
        """
        Compute variance of PH distribution.

        Var(X) = E[X^2] - E[X]^2
        E[X^2] = 2 * alpha * T^{-2} * 1
        """
        try:
            T_inv = np.linalg.inv(self.T)
            ones = np.ones(self.n_phases)
            T_inv_sq = np.dot(T_inv, T_inv)
            mean = -np.dot(self.alpha, np.dot(T_inv, ones))
            second_moment = 2 * np.dot(self.alpha, np.dot(T_inv_sq, ones))
            return second_moment - mean * mean
        except np.linalg.LinAlgError:
            return float('nan')


class APHSampler(PHSampler):
    """
    Acyclic Phase-Type distribution sampler.

    Same as PHSampler but the T matrix is upper triangular (acyclic).
    This is more efficient to simulate as phases are visited in order.
    """

    def __init__(
        self,
        alpha: np.ndarray,
        T: np.ndarray,
        rng: Optional[np.random.Generator] = None
    ):
        super().__init__(alpha, T, rng)
        # Verify acyclicity (T should be upper triangular-ish)
        # We don't enforce this strictly, but could validate

    def sample(self) -> float:
        """
        Sample from APH using efficient forward simulation.

        For acyclic PH, we can use the fact that phases are visited
        in a specific order (forward only).
        """
        total_time = 0.0

        # Select initial phase
        u = self.rng.random()
        phase = np.searchsorted(self.alpha_cumsum, u)
        phase = min(phase, self.n_phases - 1)

        # Simulate until absorption (forward only)
        while phase < self.n_phases:
            rate = self.total_rates[phase]
            if rate <= 0:
                break

            # Time in current phase
            sojourn = self.rng.exponential(1.0 / rate)
            total_time += sojourn

            # Transition (only forward or absorb in acyclic)
            u = self.rng.random()
            next_state = np.searchsorted(self.trans_probs[phase], u)

            if next_state >= self.n_phases:
                break

            phase = next_state

        return total_time


class CoxianSampler(DistributionSampler):
    """
    Coxian distribution sampler.

    A Coxian distribution is a series of exponential phases where
    at each phase the process can either continue or absorb.

    Parameters:
        rates: Service rates at each phase (mu_1, mu_2, ..., mu_n)
        probs: Continuation probabilities (p_1, p_2, ..., p_{n-1})
               p_i = probability of continuing to phase i+1
    """

    def __init__(
        self,
        rates: List[float],
        probs: List[float],
        rng: Optional[np.random.Generator] = None
    ):
        super().__init__(rng)

        self.rates = np.asarray(rates)
        self.probs = np.asarray(probs) if probs else np.array([])

        self.n_phases = len(rates)
        if len(self.probs) != max(0, self.n_phases - 1):
            raise ValueError("Need n-1 continuation probabilities for n phases")

        # Precompute scales for exponential sampling
        self.scales = 1.0 / self.rates

    def sample(self) -> float:
        """
        Sample from Coxian by phase simulation.
        """
        total_time = 0.0
        phase = 0

        while phase < self.n_phases:
            # Time in current phase
            sojourn = self.rng.exponential(self.scales[phase])
            total_time += sojourn

            # Check if we continue to next phase
            if phase < self.n_phases - 1:
                if self.rng.random() < self.probs[phase]:
                    phase += 1
                else:
                    break
            else:
                break

        return total_time

    def get_mean(self) -> float:
        """
        Compute mean of Coxian distribution.
        """
        mean = 0.0
        survival_prob = 1.0  # Probability of reaching each phase

        for i in range(self.n_phases):
            mean += survival_prob / self.rates[i]
            if i < len(self.probs):
                survival_prob *= self.probs[i]

        return mean

    def get_variance(self) -> float:
        """
        Compute variance of Coxian distribution.

        Uses the formula for variance of a sum of dependent random variables.
        """
        # Build PH representation and use that
        alpha, T = self._to_ph()
        ph = PHSampler(alpha, T, self.rng)
        return ph.get_variance()

    def _to_ph(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Convert Coxian to PH representation.
        """
        n = self.n_phases
        alpha = np.zeros(n)
        alpha[0] = 1.0

        T = np.zeros((n, n))
        for i in range(n):
            T[i, i] = -self.rates[i]
            if i < n - 1:
                T[i, i + 1] = self.rates[i] * self.probs[i]

        return alpha, T


class Coxian2Sampler(DistributionSampler):
    """
    Two-phase Coxian distribution sampler.

    Optimized implementation for the common 2-phase case.

    Parameters:
        mu1: Rate of phase 1
        mu2: Rate of phase 2
        phi: Probability of completing in phase 1 (bypassing phase 2)
    """

    def __init__(
        self,
        mu1: float,
        mu2: float,
        phi: float,
        rng: Optional[np.random.Generator] = None
    ):
        super().__init__(rng)
        if mu1 <= 0 or mu2 <= 0:
            raise ValueError("Rates must be positive")
        if not 0 <= phi <= 1:
            raise ValueError("Phi must be in [0, 1]")

        self.mu1 = mu1
        self.mu2 = mu2
        self.phi = phi
        self.scale1 = 1.0 / mu1
        self.scale2 = 1.0 / mu2

    def sample(self) -> float:
        # Phase 1
        time = self.rng.exponential(self.scale1)

        # With probability (1-phi), continue to phase 2
        if self.rng.random() >= self.phi:
            time += self.rng.exponential(self.scale2)

        return time

    def get_mean(self) -> float:
        return 1.0 / self.mu1 + (1.0 - self.phi) / self.mu2

    def get_variance(self) -> float:
        var1 = 1.0 / (self.mu1 ** 2)
        var2 = (1.0 - self.phi) / (self.mu2 ** 2)
        # Covariance term is 0 since phases are independent given path
        return var1 + var2


class HyperExponentialSampler(DistributionSampler):
    """
    Hyper-exponential (mixture of exponentials) distribution.

    f(x) = sum_i p_i * mu_i * exp(-mu_i * x)

    Parameters:
        probabilities: Mixing probabilities (must sum to 1)
        rates: Exponential rates for each component
    """

    def __init__(
        self,
        probabilities: List[float],
        rates: List[float],
        rng: Optional[np.random.Generator] = None
    ):
        super().__init__(rng)

        self.probs = np.asarray(probabilities)
        self.rates = np.asarray(rates)

        if len(self.probs) != len(self.rates):
            raise ValueError("Number of probabilities must match number of rates")
        if abs(np.sum(self.probs) - 1.0) > 1e-10:
            raise ValueError("Probabilities must sum to 1")
        if np.any(self.rates <= 0):
            raise ValueError("All rates must be positive")

        self.scales = 1.0 / self.rates
        self.cumprobs = np.cumsum(self.probs)

    def sample(self) -> float:
        # Select component
        u = self.rng.random()
        idx = np.searchsorted(self.cumprobs, u)
        idx = min(idx, len(self.rates) - 1)

        # Sample from selected exponential
        return self.rng.exponential(self.scales[idx])

    def get_mean(self) -> float:
        return np.sum(self.probs / self.rates)

    def get_variance(self) -> float:
        mean = self.get_mean()
        second_moment = 2 * np.sum(self.probs / (self.rates ** 2))
        return second_moment - mean ** 2

    def get_scv(self) -> float:
        """
        For hyper-exponential, SCV >= 1 always.
        """
        mean = self.get_mean()
        variance = self.get_variance()
        return variance / (mean ** 2)


class HyperExp2Sampler(DistributionSampler):
    """
    Two-phase hyper-exponential distribution.

    Optimized for the common 2-component case.

    Parameters:
        p: Probability of component 1
        mu1: Rate of component 1
        mu2: Rate of component 2
    """

    def __init__(
        self,
        p: float,
        mu1: float,
        mu2: float,
        rng: Optional[np.random.Generator] = None
    ):
        super().__init__(rng)
        if not 0 <= p <= 1:
            raise ValueError("p must be in [0, 1]")
        if mu1 <= 0 or mu2 <= 0:
            raise ValueError("Rates must be positive")

        self.p = p
        self.mu1 = mu1
        self.mu2 = mu2
        self.scale1 = 1.0 / mu1
        self.scale2 = 1.0 / mu2

    def sample(self) -> float:
        if self.rng.random() < self.p:
            return self.rng.exponential(self.scale1)
        else:
            return self.rng.exponential(self.scale2)

    def get_mean(self) -> float:
        return self.p / self.mu1 + (1.0 - self.p) / self.mu2

    def get_variance(self) -> float:
        mean = self.get_mean()
        second_moment = 2 * (self.p / self.mu1 ** 2 + (1.0 - self.p) / self.mu2 ** 2)
        return second_moment - mean ** 2


def create_ph_from_moments(mean: float, scv: float, rng: Optional[np.random.Generator] = None) -> DistributionSampler:
    """
    Create a phase-type distribution matching given mean and SCV.

    Uses appropriate distribution based on SCV:
    - SCV = 1: Exponential
    - SCV < 1: Erlang (approximation)
    - SCV > 1: Hyper-exponential (approximation)

    Parameters:
        mean: Target mean
        scv: Target squared coefficient of variation

    Returns:
        DistributionSampler instance
    """
    from .continuous import ExponentialSampler, ErlangSampler

    if abs(scv - 1.0) < 0.01:
        # Use exponential
        rate = 1.0 / mean
        return ExponentialSampler(rate, rng)

    elif scv < 1.0:
        # Use Erlang with k = ceil(1/SCV) phases
        k = max(1, int(np.ceil(1.0 / scv)))
        rate = k / mean
        return ErlangSampler(k, rate, rng)

    else:
        # Use 2-phase hyper-exponential
        # Match mean and SCV using standard fitting
        # p * (1/mu1) + (1-p) * (1/mu2) = mean
        # 2 * [p/mu1^2 + (1-p)/mu2^2] / mean^2 = scv + 1

        # Use balanced hyperexponential: p = 0.5
        # Then: 1/(2*mu1) + 1/(2*mu2) = mean
        #       1/mu1^2 + 1/mu2^2 = mean^2 * (scv + 1)

        # Solve for rates
        cv = np.sqrt(scv)
        p = 0.5

        # Using ratio r = mu2/mu1
        # From balanced fitting: r = (cv + 1) / (cv - 1) when cv > 1
        if cv > 1.001:
            r = (cv + 1) / (cv - 1)
            # From mean equation: 1/mu1 + 1/mu2 = 2*mean
            # 1/mu1 * (1 + 1/r) = 2*mean
            mu1 = (1 + 1/r) / (2 * mean)
            mu2 = r * mu1
        else:
            # cv ~ 1, use exponential-like
            mu1 = mu2 = 1.0 / mean

        return HyperExp2Sampler(p, mu1, mu2, rng)
