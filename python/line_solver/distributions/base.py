"""
Base distribution classes for LINE (pure Python).

This module provides the base classes for probability distributions
implemented in pure Python.
"""

from abc import ABC, abstractmethod
from typing import Any, Optional, Tuple, Union
import numpy as np


class Distribution(ABC):
    """
    Base class for all probability distributions.

    This class provides the common interface for all probability distributions
    in LINE, including service time distributions and inter-arrival time distributions.
    """

    def __init__(self):
        """Initialize a new distribution."""
        self._name = 'Distribution'

    @property
    def name(self) -> str:
        """Get the distribution name."""
        return self._name

    def getName(self) -> str:
        """Get the name of this distribution."""
        return self._name

    @abstractmethod
    def getMean(self) -> float:
        """Get the mean (expected value) of the distribution."""
        pass

    @abstractmethod
    def getVar(self) -> float:
        """Get the variance of the distribution."""
        pass

    # Abstract methods for snake_case aliases
    def get_mean(self) -> float:
        """Get the mean (expected value) of the distribution."""
        return self.getMean()

    def get_var(self) -> float:
        """Get the variance of the distribution."""
        return self.getVar()

    def getSCV(self) -> float:
        """
        Get the squared coefficient of variation (SCV).

        SCV = Var[X] / E[X]^2
        """
        mean = self.getMean()
        if mean == 0:
            return float('inf')
        return self.getVar() / (mean ** 2)

    def getRate(self) -> float:
        """Get the rate parameter (1/mean)."""
        mean = self.getMean()
        if mean == 0:
            return float('inf')
        return 1.0 / mean

    def getSkew(self) -> float:
        """Get the skewness of the distribution."""
        return 0.0  # Default, override in subclasses

    def getSupport(self) -> Tuple[float, float]:
        """Get the support range [min, max] of the distribution."""
        return (0.0, float('inf'))

    def isContinuous(self) -> bool:
        """Check if this distribution is continuous."""
        return True

    def isDiscrete(self) -> bool:
        """Check if this distribution is discrete."""
        return False

    def isDisabled(self) -> bool:
        """Check if this distribution is disabled."""
        return False

    def isImmediate(self) -> bool:
        """Check if this distribution represents immediate service."""
        return False

    def evalCDF(self, x: float) -> float:
        """Evaluate the cumulative distribution function at point x."""
        raise NotImplementedError("CDF not implemented for this distribution")

    def evalPDF(self, x: float) -> float:
        """Evaluate the probability density function at point x."""
        raise NotImplementedError("PDF not implemented for this distribution")

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """
        Generate random samples from this distribution.

        Args:
            n: Number of samples to generate.
            rng: Optional random number generator.

        Returns:
            Array of n random samples.
        """
        raise NotImplementedError("Sampling not implemented for this distribution")

    # Aliases (note: get_mean and get_var are already implemented as methods above)
    get_name = getName
    get_scv = getSCV
    get_rate = getRate
    get_skew = getSkew
    get_support = getSupport


class ContinuousDistribution(Distribution):
    """Base class for continuous probability distributions."""

    def __init__(self):
        super().__init__()

    def isContinuous(self) -> bool:
        return True

    def isDiscrete(self) -> bool:
        return False


class DiscreteDistribution(Distribution):
    """Base class for discrete probability distributions."""

    def __init__(self):
        super().__init__()

    def isContinuous(self) -> bool:
        return False

    def isDiscrete(self) -> bool:
        return True

    def evalPMF(self, x: int) -> float:
        """Evaluate the probability mass function at point x."""
        raise NotImplementedError("PMF not implemented for this distribution")


class Markovian(Distribution):
    """
    Base class for Markovian (phase-type) distributions.

    Markovian distributions can be represented in terms of
    initial probability vectors and transition rate matrices.
    """

    def __init__(self):
        super().__init__()

    def getD0(self) -> np.ndarray:
        """Get the D0 matrix (for MAP representations)."""
        raise NotImplementedError

    def getD1(self) -> np.ndarray:
        """Get the D1 matrix (for MAP representations)."""
        raise NotImplementedError

    def getMu(self) -> np.ndarray:
        """Get the service rates in each phase."""
        raise NotImplementedError

    def getPhi(self) -> np.ndarray:
        """Get the completion probabilities from each phase."""
        raise NotImplementedError

    def getInitProb(self) -> np.ndarray:
        """Get the initial probability vector."""
        raise NotImplementedError

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        raise NotImplementedError

    # Aliases
    get_d0 = getD0
    get_d1 = getD1
    get_mu = getMu
    get_phi = getPhi
    get_init_prob = getInitProb
    get_number_of_phases = getNumberOfPhases
