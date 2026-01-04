"""
Base class for distribution samplers.

This module provides the abstract DistributionSampler base class that
all distribution implementations inherit from.
"""

from abc import ABC, abstractmethod
from typing import Optional
import numpy as np


class DistributionSampler(ABC):
    """
    Abstract base class for distribution samplers.

    Distribution samplers generate random variates for arrival times,
    service times, and other stochastic elements in the simulation.

    Subclasses must implement:
        - sample(): Generate a single random variate

    The random number generator can be provided externally for
    reproducibility (using the same seed produces the same sequence).

    Attributes:
        rng: NumPy random generator instance
    """

    def __init__(self, rng: Optional[np.random.Generator] = None):
        """
        Initialize the distribution sampler.

        Args:
            rng: NumPy random generator. If None, creates a new default generator.
        """
        self.rng = rng if rng is not None else np.random.default_rng()

    @abstractmethod
    def sample(self) -> float:
        """
        Generate a single random sample from this distribution.

        Returns:
            A random variate from the distribution (always non-negative)
        """
        pass

    def sample_n(self, n: int) -> np.ndarray:
        """
        Generate n random samples from this distribution.

        Args:
            n: Number of samples to generate

        Returns:
            Array of n random variates
        """
        return np.array([self.sample() for _ in range(n)])

    def get_mean(self) -> float:
        """
        Get the theoretical mean of the distribution.

        Returns:
            Expected value (mean) of the distribution.
            Returns NaN if not implemented or undefined.
        """
        return float('nan')

    def get_variance(self) -> float:
        """
        Get the theoretical variance of the distribution.

        Returns:
            Variance of the distribution.
            Returns NaN if not implemented or undefined.
        """
        return float('nan')

    def get_scv(self) -> float:
        """
        Get the squared coefficient of variation (SCV).

        SCV = Var(X) / E[X]^2

        Returns:
            Squared coefficient of variation.
            Returns NaN if not implemented or undefined.
        """
        mean = self.get_mean()
        if mean == 0 or np.isnan(mean):
            return float('nan')
        variance = self.get_variance()
        if np.isnan(variance):
            return float('nan')
        return variance / (mean * mean)


class DisabledSampler(DistributionSampler):
    """
    Sampler for disabled distributions (service bypassed).

    Jobs with disabled service immediately pass through the station
    without any service time.

    Raises:
        RuntimeError: If sample() is called (should never happen)
    """

    def sample(self) -> float:
        """
        Should never be called for disabled distributions.

        Raises:
            RuntimeError: Always raises, as disabled should bypass service
        """
        raise RuntimeError("Cannot sample from disabled distribution")

    def get_mean(self) -> float:
        return 0.0

    def get_variance(self) -> float:
        return 0.0


class ImmediateSampler(DistributionSampler):
    """
    Sampler for immediate (zero) service time.

    Jobs receive instant service (zero service time).
    """

    def sample(self) -> float:
        """Return zero service time."""
        return 0.0

    def get_mean(self) -> float:
        return 0.0

    def get_variance(self) -> float:
        return 0.0
