"""
Continuous probability distributions for LINE (pure Python).

This module provides continuous distribution implementations including
exponential, deterministic, Erlang, hyperexponential, and other common
service time distributions.
"""

from typing import Optional, Tuple, Union, List
import numpy as np
from scipy import stats, linalg

from .base import ContinuousDistribution, Markovian


class Exp(ContinuousDistribution, Markovian):
    """
    Exponential distribution.

    The exponential distribution is the simplest continuous distribution
    for modeling service times in queueing systems. It has the memoryless
    property and SCV = 1.

    Args:
        rate: The rate parameter (lambda = 1/mean).
    """

    def __init__(self, rate: float):
        super().__init__()
        self._name = 'Exp'
        if rate <= 0:
            raise ValueError("Rate must be positive")
        self._rate = rate

    @classmethod
    def fit_mean(cls, mean: float) -> 'Exp':
        """
        Create an exponential distribution with the given mean.

        Args:
            mean: Target mean.

        Returns:
            Exp distribution with rate = 1/mean.
        """
        if mean <= 0:
            raise ValueError("Mean must be positive")
        return cls(rate=1.0 / mean)

    # CamelCase alias
    fitMean = fit_mean

    @property
    def rate(self) -> float:
        """Get the rate parameter."""
        return self._rate

    @rate.setter
    def rate(self, value: float):
        """Set the rate parameter."""
        if value <= 0:
            raise ValueError("Rate must be positive")
        self._rate = value
        

    def getMean(self) -> float:
        """Get the mean (1/rate)."""
        return 1.0 / self._rate

    def getVar(self) -> float:
        """Get the variance (1/rate^2)."""
        return 1.0 / (self._rate ** 2)

    def getSCV(self) -> float:
        """Get the squared coefficient of variation (always 1 for exponential)."""
        return 1.0

    def getSkew(self) -> float:
        """Get the skewness (always 2 for exponential)."""
        return 2.0

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < 0:
            return 0.0
        return 1.0 - np.exp(-self._rate * x)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < 0:
            return 0.0
        return self._rate * np.exp(-self._rate * x)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.exponential(scale=1.0/self._rate, size=n)

    def getNumberOfPhases(self) -> int:
        """Get the number of phases (1 for exponential)."""
        return 1

    def getD0(self) -> np.ndarray:
        """Get the D0 matrix for MAP representation."""
        return np.array([[-self._rate]])

    def getD1(self) -> np.ndarray:
        """Get the D1 matrix for MAP representation."""
        return np.array([[self._rate]])

    def getMu(self) -> np.ndarray:
        """Get the service rates in each phase."""
        return np.array([self._rate])

    def getPhi(self) -> np.ndarray:
        """Get the completion probabilities from each phase."""
        return np.array([1.0])

    def getInitProb(self) -> np.ndarray:
        """Get the initial probability vector."""
        return np.array([1.0])


    @classmethod
    def fitRate(cls, rate: float) -> 'Exp':
        """
        Create an exponential distribution with the given rate.

        Args:
            rate: The rate parameter (lambda).

        Returns:
            Exp distribution with specified rate.
        """
        return cls(rate=rate)

    # Snake_case alias
    fit_rate = fitRate

    # Aliases
    get_rate = lambda self: self._rate


class Det(ContinuousDistribution):
    """
    Deterministic (constant) distribution.

    All service times are exactly equal to the specified value.
    Has SCV = 0 (no variability).

    Args:
        value: The constant service time value.
    """

    def __init__(self, value: float):
        super().__init__()
        self._name = 'Det'
        if value < 0:
            raise ValueError("Value must be non-negative")
        self._value = value

    @property
    def value(self) -> float:
        """Get the constant value."""
        return self._value

    @value.setter
    def value(self, val: float):
        """Set the constant value."""
        if val < 0:
            raise ValueError("Value must be non-negative")
        self._value = val
        

    def getMean(self) -> float:
        """Get the mean (equals the constant value)."""
        return self._value

    def getVar(self) -> float:
        """Get the variance (always 0)."""
        return 0.0

    def getSCV(self) -> float:
        """Get the SCV (always 0)."""
        return 0.0

    def getSkew(self) -> float:
        """Get the skewness (undefined, return 0)."""
        return 0.0

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        return 1.0 if x >= self._value else 0.0

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x (delta function, return inf at value)."""
        return float('inf') if x == self._value else 0.0

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples (all equal to value)."""
        return np.full(n, self._value)

    def isImmediate(self) -> bool:
        """Check if this is an immediate (zero) service."""
        return self._value == 0.0



class Immediate(Det):
    """
    Immediate (zero delay) distribution.

    Represents instantaneous service with zero delay.
    """

    _instance = None

    def __init__(self):
        super().__init__(0.0)
        self._name = 'Immediate'

    @classmethod
    def getInstance(cls) -> 'Immediate':
        """Get singleton instance of Immediate distribution."""
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    # snake_case alias
    get_instance = getInstance

    def isImmediate(self) -> bool:
        """Check if this is immediate service."""
        return True


class Disabled(ContinuousDistribution):
    """
    Disabled distribution.

    Represents a disabled service (no service at all).
    Used for nodes that don't serve a particular job class.
    """

    _instance = None

    def __init__(self):
        super().__init__()
        self._name = 'Disabled'

    @classmethod
    def getInstance(cls) -> 'Disabled':
        """Get singleton instance of Disabled distribution."""
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    # snake_case alias
    get_instance = getInstance

    def getMean(self) -> float:
        """Get the mean (infinity for disabled)."""
        return float('inf')

    def getVar(self) -> float:
        """Get the variance."""
        return float('inf')

    def isDisabled(self) -> bool:
        """Check if this distribution is disabled."""
        return True



class Erlang(ContinuousDistribution, Markovian):
    """
    Erlang distribution (sum of k exponentials).

    The Erlang distribution is the distribution of the sum of k
    independent exponential random variables with the same rate.
    It has SCV = 1/k.

    Args:
        phase_rate: Rate parameter for each exponential phase (alpha).
        nphases: Number of sequential exponential phases (r).

    Mean = nphases / phase_rate = r / alpha
    """

    def __init__(self, phase_rate: float, nphases: int):
        super().__init__()
        self._name = 'Erlang'
        if phase_rate <= 0:
            raise ValueError("Phase rate must be positive")
        if nphases < 1:
            raise ValueError("Number of phases must be at least 1")
        self._phase_rate = phase_rate
        self._phases = int(round(nphases))
        # Mean = nphases / phase_rate
        self._mean = self._phases / self._phase_rate

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float) -> 'Erlang':
        """
        Create an Erlang distribution from mean and SCV.

        For Erlang, SCV = 1/k where k is number of phases.
        So phases = round(1/SCV), constrained to be >= 1.

        Args:
            mean: Target mean.
            scv: Target squared coefficient of variation.

        Returns:
            Erlang distribution with given mean and closest achievable SCV.
        """
        if scv <= 0:
            raise ValueError("SCV must be positive")
        if mean <= 0:
            raise ValueError("Mean must be positive")
        phases = max(1, round(1.0 / scv))
        phase_rate = phases / mean
        return cls(phase_rate=phase_rate, nphases=phases)

    @classmethod
    def fit_mean_and_order(cls, mean: float, phases: int) -> 'Erlang':
        """
        Create an Erlang distribution from mean and number of phases.

        Args:
            mean: Target mean.
            phases: Number of phases (order).

        Returns:
            Erlang distribution with given mean and phases.
        """
        if mean <= 0:
            raise ValueError("Mean must be positive")
        phase_rate = phases / mean
        return cls(phase_rate=phase_rate, nphases=phases)

    # CamelCase aliases
    fitMeanAndScv = fit_mean_and_scv
    fitMeanAndOrder = fit_mean_and_order

    @property
    def phases(self) -> int:
        """Get the number of phases."""
        return self._phases

    def getMean(self) -> float:
        """Get the mean."""
        return self._mean

    def getVar(self) -> float:
        """Get the variance."""
        return self._mean ** 2 / self._phases

    def getSCV(self) -> float:
        """Get the SCV (1/phases)."""
        return 1.0 / self._phases

    def getSkew(self) -> float:
        """Get the skewness."""
        return 2.0 / np.sqrt(self._phases)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x <= 0:
            return 0.0
        return stats.gamma.cdf(x, a=self._phases, scale=1.0/self._phase_rate)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < 0:
            return 0.0
        return stats.gamma.pdf(x, a=self._phases, scale=1.0/self._phase_rate)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.gamma(shape=self._phases, scale=1.0/self._phase_rate, size=n)

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        return self._phases

    def getD0(self) -> np.ndarray:
        """Get the D0 matrix for MAP representation."""
        k = self._phases
        D0 = np.zeros((k, k))
        for i in range(k):
            D0[i, i] = -self._phase_rate
            if i < k - 1:
                D0[i, i + 1] = self._phase_rate
        return D0

    def getD1(self) -> np.ndarray:
        """Get the D1 matrix for MAP representation."""
        k = self._phases
        D1 = np.zeros((k, k))
        D1[k - 1, 0] = self._phase_rate
        return D1

    def getMu(self) -> np.ndarray:
        """Get the service rates in each phase."""
        return np.full(self._phases, self._phase_rate)

    def getPhi(self) -> np.ndarray:
        """Get the completion probabilities from each phase."""
        phi = np.zeros(self._phases)
        phi[-1] = 1.0
        return phi

    def getInitProb(self) -> np.ndarray:
        """Get the initial probability vector."""
        alpha = np.zeros(self._phases)
        alpha[0] = 1.0
        return alpha



class HyperExp(ContinuousDistribution, Markovian):
    """
    Hyperexponential distribution (mixture of exponentials).

    The hyperexponential distribution is a mixture of exponential
    distributions. It has SCV >= 1.

    Supports two calling conventions (matching MATLAB API):
        - HyperExp(p, rate1, rate2): 2-phase with probability p of rate1,
          probability (1-p) of rate2
        - HyperExp(probs, rates): n-phase with vectors of probabilities and rates

    Args:
        p_or_probs: For 2-phase: probability of first component (scalar).
                    For n-phase: list/array of probabilities.
        rate1_or_rates: For 2-phase: rate of first component (scalar).
                        For n-phase: list/array of rates.
        rate2: For 2-phase only: rate of second component.
    """

    def __init__(self, p_or_probs: Union[float, list, np.ndarray],
                 rate1_or_rates: Union[float, list, np.ndarray],
                 rate2: Optional[float] = None):
        super().__init__()
        self._name = 'HyperExp'

        # Determine if this is 2-phase (3 scalar args) or n-phase (2 vector args)
        if rate2 is not None:
            # 2-phase form: HyperExp(p, rate1, rate2)
            p = float(p_or_probs)
            r1 = float(rate1_or_rates)
            r2 = float(rate2)
            if not (0 <= p <= 1):
                raise ValueError("Probability p must be in [0, 1]")
            if r1 <= 0 or r2 <= 0:
                raise ValueError("All rates must be positive")
            self._probs = np.array([p, 1.0 - p], dtype=float)
            self._rates = np.array([r1, r2], dtype=float)
        else:
            # n-phase form: HyperExp(probs, rates)
            self._probs = np.array(p_or_probs, dtype=float)
            self._rates = np.array(rate1_or_rates, dtype=float)

            if len(self._probs) != len(self._rates):
                raise ValueError("Probabilities and rates must have same length")
            if not np.allclose(np.sum(self._probs), 1.0):
                raise ValueError("Probabilities must sum to 1")
            if np.any(self._rates <= 0):
                raise ValueError("All rates must be positive")
            if np.any(self._probs < 0):
                raise ValueError("All probabilities must be non-negative")

        self._means = 1.0 / self._rates

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float) -> 'HyperExp':
        """
        Create a 2-phase hyperexponential distribution from mean and SCV.

        Uses balanced means representation for SCV >= 1.

        Args:
            mean: Target mean.
            scv: Target squared coefficient of variation (must be >= 1).

        Returns:
            HyperExp distribution with given mean and SCV.
        """
        if scv < 1.0:
            raise ValueError("HyperExp requires SCV >= 1")
        if mean <= 0:
            raise ValueError("Mean must be positive")

        # For 2-phase HyperExp with balanced means:
        # p1 = 0.5 * (1 + sqrt((scv-1)/(scv+1)))
        # p2 = 1 - p1
        # mu1 = 2 * p1 / mean
        # mu2 = 2 * p2 / mean
        cv = np.sqrt(scv)
        p1 = 0.5 * (1.0 + np.sqrt((scv - 1.0) / (scv + 1.0)))
        p2 = 1.0 - p1

        # Rates
        mu1 = 2.0 * p1 / mean
        mu2 = 2.0 * p2 / mean

        # Use the 3-argument form: HyperExp(p, rate1, rate2)
        return cls(p1, mu1, mu2)

    # CamelCase alias
    fitMeanAndScv = fit_mean_and_scv

    @property
    def means(self) -> np.ndarray:
        """Get the means of each component."""
        return self._means

    @property
    def probs(self) -> np.ndarray:
        """Get the probabilities of each component."""
        return self._probs

    @property
    def rates(self) -> np.ndarray:
        """Get the rates of each component."""
        return self._rates

    def getMean(self) -> float:
        """Get the mean."""
        return float(np.sum(self._probs * self._means))

    def getVar(self) -> float:
        """Get the variance."""
        mean = self.getMean()
        second_moment = 2 * np.sum(self._probs * self._means ** 2)
        return second_moment - mean ** 2

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x <= 0:
            return 0.0
        cdf = 0.0
        for p, r in zip(self._probs, self._rates):
            cdf += p * (1.0 - np.exp(-r * x))
        return cdf

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < 0:
            return 0.0
        pdf = 0.0
        for p, r in zip(self._probs, self._rates):
            pdf += p * r * np.exp(-r * x)
        return pdf

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        components = rng.choice(len(self._means), size=n, p=self._probs)
        samples = np.zeros(n)
        for i, c in enumerate(components):
            samples[i] = rng.exponential(scale=self._means[c])
        return samples

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        return len(self._means)

    def getD0(self) -> np.ndarray:
        """Get the D0 matrix for MAP representation."""
        k = len(self._rates)
        return np.diag(-self._rates)

    def getD1(self) -> np.ndarray:
        """Get the D1 matrix for MAP representation."""
        k = len(self._rates)
        D1 = np.zeros((k, k))
        for i in range(k):
            for j in range(k):
                D1[i, j] = self._rates[i] * self._probs[j]
        return D1

    def getMu(self) -> np.ndarray:
        """Get the service rates in each phase."""
        return self._rates.copy()

    def getPhi(self) -> np.ndarray:
        """Get the completion probabilities from each phase."""
        return np.ones(len(self._rates))

    def getInitProb(self) -> np.ndarray:
        """Get the initial probability vector."""
        return self._probs.copy()



class Gamma(ContinuousDistribution):
    """
    Gamma distribution.

    The gamma distribution is a two-parameter continuous distribution
    that generalizes the exponential and Erlang distributions.

    Args:
        shape: Shape parameter (k or alpha).
        scale: Scale parameter (theta).
    """

    def __init__(self, shape: float, scale: float):
        super().__init__()
        self._name = 'Gamma'
        if shape <= 0:
            raise ValueError("Shape must be positive")
        if scale <= 0:
            raise ValueError("Scale must be positive")
        self._shape = shape
        self._scale = scale

    @property
    def shape(self) -> float:
        """Get the shape parameter."""
        return self._shape

    @property
    def scale(self) -> float:
        """Get the scale parameter."""
        return self._scale

    def getMean(self) -> float:
        """Get the mean (shape * scale)."""
        return self._shape * self._scale

    def getVar(self) -> float:
        """Get the variance (shape * scale^2)."""
        return self._shape * self._scale ** 2

    def getSCV(self) -> float:
        """Get the SCV (1/shape)."""
        return 1.0 / self._shape

    def getSkew(self) -> float:
        """Get the skewness."""
        return 2.0 / np.sqrt(self._shape)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x <= 0:
            return 0.0
        return stats.gamma.cdf(x, a=self._shape, scale=self._scale)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < 0:
            return 0.0
        return stats.gamma.pdf(x, a=self._shape, scale=self._scale)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.gamma(shape=self._shape, scale=self._scale, size=n)



class Lognormal(ContinuousDistribution):
    """
    Lognormal distribution.

    A random variable X has a lognormal distribution if log(X) is
    normally distributed.

    Args:
        mu: Mean of the underlying normal distribution.
        sigma: Standard deviation of the underlying normal distribution.
    """

    def __init__(self, mu: float, sigma: float):
        super().__init__()
        self._name = 'Lognormal'
        if sigma <= 0:
            raise ValueError("Sigma must be positive")
        self._mu = mu
        self._sigma = sigma

    @property
    def mu(self) -> float:
        """Get the mu parameter."""
        return self._mu

    @property
    def sigma(self) -> float:
        """Get the sigma parameter."""
        return self._sigma

    def getMean(self) -> float:
        """Get the mean."""
        return np.exp(self._mu + self._sigma ** 2 / 2)

    def getVar(self) -> float:
        """Get the variance."""
        return (np.exp(self._sigma ** 2) - 1) * np.exp(2 * self._mu + self._sigma ** 2)

    def getSkew(self) -> float:
        """Get the skewness."""
        es2 = np.exp(self._sigma ** 2)
        return (es2 + 2) * np.sqrt(es2 - 1)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x <= 0:
            return 0.0
        return stats.lognorm.cdf(x, s=self._sigma, scale=np.exp(self._mu))

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x <= 0:
            return 0.0
        return stats.lognorm.pdf(x, s=self._sigma, scale=np.exp(self._mu))

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.lognormal(mean=self._mu, sigma=self._sigma, size=n)



class Pareto(ContinuousDistribution):
    """
    Pareto distribution.

    The Pareto distribution is a power-law distribution often used
    to model heavy-tailed phenomena.

    Args:
        alpha: Shape parameter (tail index).
        scale: Scale parameter (minimum value).
    """

    def __init__(self, alpha: float, scale: float):
        super().__init__()
        self._name = 'Pareto'
        if alpha <= 0:
            raise ValueError("Alpha must be positive")
        if scale <= 0:
            raise ValueError("Scale must be positive")
        self._alpha = alpha
        self._scale = scale

    @property
    def alpha(self) -> float:
        """Get the alpha parameter."""
        return self._alpha

    @property
    def scale(self) -> float:
        """Get the scale parameter."""
        return self._scale

    def getMean(self) -> float:
        """Get the mean."""
        if self._alpha <= 1:
            return float('inf')
        return self._alpha * self._scale / (self._alpha - 1)

    def getVar(self) -> float:
        """Get the variance."""
        if self._alpha <= 2:
            return float('inf')
        return (self._scale ** 2 * self._alpha) / ((self._alpha - 1) ** 2 * (self._alpha - 2))

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [scale, inf)."""
        return (self._scale, float('inf'))

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < self._scale:
            return 0.0
        return 1.0 - (self._scale / x) ** self._alpha

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < self._scale:
            return 0.0
        return self._alpha * self._scale ** self._alpha / x ** (self._alpha + 1)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return (rng.pareto(a=self._alpha, size=n) + 1) * self._scale



class Uniform(ContinuousDistribution):
    """
    Uniform distribution on [min, max].

    Args:
        min_val: Minimum value.
        max_val: Maximum value.
    """

    def __init__(self, min_val: float, max_val: float):
        super().__init__()
        self._name = 'Uniform'
        if max_val < min_val:
            raise ValueError("max_val must be >= min_val")
        self._min = min_val
        self._max = max_val

    @property
    def min_val(self) -> float:
        """Get the minimum value."""
        return self._min

    @property
    def max_val(self) -> float:
        """Get the maximum value."""
        return self._max

    def getMean(self) -> float:
        """Get the mean."""
        return (self._min + self._max) / 2

    def getVar(self) -> float:
        """Get the variance."""
        return (self._max - self._min) ** 2 / 12

    def getSkew(self) -> float:
        """Get the skewness (always 0 for uniform)."""
        return 0.0

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [min, max]."""
        return (self._min, self._max)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < self._min:
            return 0.0
        if x > self._max:
            return 1.0
        return (x - self._min) / (self._max - self._min)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < self._min or x > self._max:
            return 0.0
        return 1.0 / (self._max - self._min)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.uniform(self._min, self._max, size=n)



class Weibull(ContinuousDistribution):
    """
    Weibull distribution.

    The Weibull distribution is commonly used in reliability engineering
    to model time to failure.

    Args:
        shape: Shape parameter (k).
        scale: Scale parameter (lambda).
    """

    def __init__(self, shape: float, scale: float):
        super().__init__()
        self._name = 'Weibull'
        if shape <= 0:
            raise ValueError("Shape must be positive")
        if scale <= 0:
            raise ValueError("Scale must be positive")
        self._shape = shape
        self._scale = scale

    @property
    def shape(self) -> float:
        """Get the shape parameter."""
        return self._shape

    @property
    def scale(self) -> float:
        """Get the scale parameter."""
        return self._scale

    def getMean(self) -> float:
        """Get the mean."""
        from scipy.special import gamma
        return self._scale * gamma(1 + 1 / self._shape)

    def getVar(self) -> float:
        """Get the variance."""
        from scipy.special import gamma
        g1 = gamma(1 + 1 / self._shape)
        g2 = gamma(1 + 2 / self._shape)
        return self._scale ** 2 * (g2 - g1 ** 2)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < 0:
            return 0.0
        return 1.0 - np.exp(-(x / self._scale) ** self._shape)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < 0:
            return 0.0
        return (self._shape / self._scale) * (x / self._scale) ** (self._shape - 1) * \
               np.exp(-(x / self._scale) ** self._shape)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return self._scale * rng.weibull(a=self._shape, size=n)



class Normal(ContinuousDistribution):
    """
    Normal (Gaussian) distribution.

    Note: For queueing applications, a truncated or shifted version
    may be needed since normal distributions can take negative values.

    Args:
        mean: Mean of the distribution.
        std: Standard deviation.
    """

    def __init__(self, mean: float, std: float):
        super().__init__()
        self._name = 'Normal'
        if std <= 0:
            raise ValueError("Standard deviation must be positive")
        self._mean_val = mean
        self._std = std

    @property
    def std(self) -> float:
        """Get the standard deviation."""
        return self._std

    def getMean(self) -> float:
        """Get the mean."""
        return self._mean_val

    def getVar(self) -> float:
        """Get the variance."""
        return self._std ** 2

    def getSCV(self) -> float:
        """Get the squared coefficient of variation (var/mean^2)."""
        if self._mean_val == 0:
            return float('inf')
        return (self._std ** 2) / (self._mean_val ** 2)

    def getStd(self) -> float:
        """Get the standard deviation."""
        return self._std

    def getSkew(self) -> float:
        """Get the skewness (always 0 for normal)."""
        return 0.0

    def getSupport(self) -> Tuple[float, float]:
        """Get the support (-inf, inf)."""
        return (float('-inf'), float('inf'))

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        return stats.norm.cdf(x, loc=self._mean_val, scale=self._std)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        return stats.norm.pdf(x, loc=self._mean_val, scale=self._std)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.normal(self._mean_val, self._std, size=n)

    @classmethod
    def fitMean(cls, mean: float, std: float = 1.0) -> 'Normal':
        """Create a Normal distribution with given mean and std."""
        return cls(mean, std)

    @classmethod
    def fitMeanAndStd(cls, mean: float, std: float) -> 'Normal':
        """Create a Normal distribution with given mean and std."""
        return cls(mean, std)

    @classmethod
    def fitMeanAndVar(cls, mean: float, var: float) -> 'Normal':
        """Create a Normal distribution with given mean and variance."""
        return cls(mean, np.sqrt(var))

    def getSkewness(self) -> float:
        """Get the skewness (always 0 for normal)."""
        return 0.0

    # Snake_case aliases
    fit_mean = fitMean
    fit_mean_and_std = fitMeanAndStd
    fit_mean_and_var = fitMeanAndVar
    get_mean = lambda self: self._mean_val
    get_std = getStd
    get_var = getVar
    get_scv = getSCV
    get_skewness = getSkewness
    eval_cdf = evalCDF
    eval_pdf = evalPDF


class MultivariateNormal(ContinuousDistribution):
    """
    Multivariate Normal (Gaussian) distribution.

    Represents a d-dimensional normal distribution with mean vector mu
    and covariance matrix Sigma.

    Args:
        mu: d-dimensional mean vector.
        Sigma: d x d positive definite covariance matrix.
    """

    def __init__(self, mu: Union[list, np.ndarray], Sigma: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'MultivariateNormal'

        self._mu = np.atleast_1d(np.array(mu, dtype=float)).flatten()
        self._Sigma = np.atleast_2d(np.array(Sigma, dtype=float))

        d = len(self._mu)
        if self._Sigma.shape != (d, d):
            raise ValueError(f"Sigma must be {d}x{d} to match mu of length {d}")

        # Check positive definite via Cholesky
        try:
            self._L = linalg.cholesky(self._Sigma, lower=True)
        except linalg.LinAlgError:
            raise ValueError("Sigma must be positive definite")

        self._dimension = d

    @property
    def dimension(self) -> int:
        """Get the dimensionality."""
        return self._dimension

    def getMeanVector(self) -> np.ndarray:
        """Get the mean vector."""
        return self._mu.copy()

    def getCovariance(self) -> np.ndarray:
        """Get the covariance matrix."""
        return self._Sigma.copy()

    def getCorrelation(self) -> np.ndarray:
        """Get the correlation matrix."""
        d = self._dimension
        R = np.zeros((d, d))
        for i in range(d):
            for j in range(d):
                std_i = np.sqrt(self._Sigma[i, i])
                std_j = np.sqrt(self._Sigma[j, j])
                if std_i > 1e-10 and std_j > 1e-10:
                    R[i, j] = self._Sigma[i, j] / (std_i * std_j)
                else:
                    R[i, j] = float(i == j)
        return R

    def getMean(self) -> float:
        """Get the mean of the first component (for compatibility)."""
        return float(self._mu[0])

    def getVar(self) -> float:
        """Get the variance of the first component."""
        return float(self._Sigma[0, 0])

    def getSkew(self) -> float:
        """Get skewness (0 for normal)."""
        return 0.0

    def evalPDF(self, x: Union[list, np.ndarray]) -> Union[float, np.ndarray]:
        """Evaluate the multivariate normal PDF at point(s) x.

        Args:
            x: Single point (1D array of length d) or multiple points (2D array of shape n x d)

        Returns:
            Single float for one point, or numpy array for multiple points.
        """
        x_arr = np.atleast_1d(np.array(x, dtype=float))

        # Handle 2D array (multiple points)
        if x_arr.ndim == 2:
            n_points = x_arr.shape[0]
            if x_arr.shape[1] != self._dimension:
                raise ValueError(f"Each point must have dimension {self._dimension}")

            inv_Sigma = linalg.inv(self._Sigma)
            det_Sigma = linalg.det(self._Sigma)
            norm_const = 1.0 / np.sqrt((2 * np.pi) ** self._dimension * det_Sigma)

            results = np.zeros(n_points)
            for i in range(n_points):
                diff = x_arr[i] - self._mu
                exponent = -0.5 * diff @ inv_Sigma @ diff
                results[i] = norm_const * np.exp(exponent)
            return results

        # Handle 1D array (single point)
        x_arr = x_arr.flatten()
        if len(x_arr) != self._dimension:
            raise ValueError(f"x must have dimension {self._dimension}")

        diff = x_arr - self._mu
        inv_Sigma = linalg.inv(self._Sigma)
        det_Sigma = linalg.det(self._Sigma)

        norm_const = 1.0 / np.sqrt((2 * np.pi) ** self._dimension * det_Sigma)
        exponent = -0.5 * diff @ inv_Sigma @ diff

        return float(norm_const * np.exp(exponent))

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """
        Generate n samples from the multivariate normal.

        Returns:
            n x d matrix of samples.
        """
        if rng is None:
            rng = np.random.default_rng()

        # X = mu + L @ Z where L = chol(Sigma), Z ~ N(0, I)
        Z = rng.standard_normal((self._dimension, n))
        X = (self._mu[:, np.newaxis] + self._L @ Z).T

        return X

    def getMarginal(self, indices: Union[list, np.ndarray]) -> 'MultivariateNormal':
        """Extract a marginal distribution for a subset of dimensions."""
        indices = np.atleast_1d(np.array(indices, dtype=int))
        mu_marg = self._mu[indices]
        Sigma_marg = self._Sigma[np.ix_(indices, indices)]
        return MultivariateNormal(mu_marg, Sigma_marg)

    def getMarginalUniv(self, index: int) -> 'Normal':
        """Extract a univariate marginal distribution."""
        mean_marg = self._mu[index]
        std_marg = np.sqrt(self._Sigma[index, index])
        return Normal(mean_marg, std_marg)

    def getDimension(self) -> int:
        """Get the dimensionality."""
        return self._dimension

    # Snake_case aliases
    def get_dimension(self) -> int:
        """Get the dimensionality."""
        return self._dimension

    def get_mean_vector(self) -> np.ndarray:
        """Get the mean vector."""
        return self._mu.copy()

    def get_covariance(self) -> np.ndarray:
        """Get the covariance matrix."""
        return self._Sigma.copy()

    def get_correlation(self) -> np.ndarray:
        """Get the correlation matrix."""
        return self.getCorrelation()

    def get_marginal(self, indices: Union[list, np.ndarray]) -> 'MultivariateNormal':
        """Extract a marginal distribution for a subset of dimensions."""
        return self.getMarginal(indices)

    def get_marginal_univ(self, index: int) -> 'Normal':
        """Extract a univariate marginal distribution."""
        return self.getMarginalUniv(index)

    def eval_pdf(self, x: Union[list, np.ndarray]) -> Union[float, np.ndarray]:
        """Evaluate the multivariate normal PDF at point(s) x."""
        x_arr = np.atleast_1d(np.array(x, dtype=float))

        # Handle multiple points: if x is 2D (n_points x d)
        if x_arr.ndim == 2:
            results = np.zeros(x_arr.shape[0])
            for i in range(x_arr.shape[0]):
                results[i] = self.evalPDF(x_arr[i])
            return results
        else:
            return self.evalPDF(x_arr)

    @classmethod
    def fitMeanAndCovariance(cls, mu: Union[list, np.ndarray],
                              Sigma: Union[list, np.ndarray]) -> 'MultivariateNormal':
        """Create a MultivariateNormal distribution with given mean and covariance."""
        return cls(mu, Sigma)

    # CamelCase to snake_case alias for fit method
    fit_mean_and_covariance = fitMeanAndCovariance


class Prior(ContinuousDistribution):
    """
    Discrete prior distribution over alternative distributions.

    Prior represents parameter uncertainty by specifying a discrete set of
    alternative distributions with associated probabilities. Used with the
    Posterior solver for Bayesian-style analysis.

    This is NOT a mixture distribution - each alternative represents a
    separate model realization.

    Args:
        distributions: List of Distribution objects.
        probabilities: List of probabilities (must sum to 1).
    """

    def __init__(self, distributions: list, probabilities: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'Prior'

        if not isinstance(distributions, list) or len(distributions) == 0:
            raise ValueError("distributions must be a non-empty list")

        self._distributions = distributions
        self._probabilities = np.array(probabilities, dtype=float)

        if len(self._distributions) != len(self._probabilities):
            raise ValueError("Number of distributions must match number of probabilities")

        if np.any(self._probabilities < 0):
            raise ValueError("Probabilities must be non-negative")

        if not np.isclose(self._probabilities.sum(), 1.0, atol=1e-6):
            raise ValueError(f"Probabilities must sum to 1 (got {self._probabilities.sum()})")

    @property
    def distributions(self) -> list:
        """Get the alternative distributions."""
        return self._distributions

    @property
    def probabilities(self) -> np.ndarray:
        """Get the probabilities."""
        return self._probabilities.copy()

    def getNumAlternatives(self) -> int:
        """Get the number of alternative distributions."""
        return len(self._distributions)

    def getAlternative(self, idx: int):
        """Get the distribution at index idx."""
        if idx < 0 or idx >= len(self._distributions):
            raise IndexError("Index out of bounds")
        return self._distributions[idx]

    def getProbability(self, idx: int) -> float:
        """Get the probability of alternative idx."""
        if idx < 0 or idx >= len(self._probabilities):
            raise IndexError("Index out of bounds")
        return float(self._probabilities[idx])

    def getMean(self) -> float:
        """Get prior-weighted mean (expected mean over alternatives)."""
        mean = 0.0
        for i in range(len(self._distributions)):
            mean += self._probabilities[i] * self._distributions[i].getMean()
        return mean

    def getVar(self) -> float:
        """Get prior-weighted variance using law of total variance."""
        E_mean = 0.0       # E[E[X|D]]
        E_var = 0.0        # E[Var(X|D)]
        E_mean_sq = 0.0    # E[E[X|D]^2]

        for i in range(len(self._distributions)):
            m = self._distributions[i].getMean()
            v = self._distributions[i].getVar()
            E_mean += self._probabilities[i] * m
            E_var += self._probabilities[i] * v
            E_mean_sq += self._probabilities[i] * m ** 2

        # Var(X) = E[Var(X|D)] + Var(E[X|D])
        return E_var + (E_mean_sq - E_mean ** 2)

    def getSCV(self) -> float:
        """Get prior-weighted SCV."""
        mean = self.getMean()
        var = self.getVar()
        return var / mean ** 2 if mean > 0 else 0.0

    def evalCDF(self, t: float) -> float:
        """Evaluate mixture CDF at t."""
        cdf = 0.0
        for i in range(len(self._distributions)):
            cdf += self._probabilities[i] * self._distributions[i].evalCDF(t)
        return cdf

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Sample from prior (mixture sampling)."""
        if rng is None:
            rng = np.random.default_rng()

        samples = np.zeros(n)
        cumprob = np.cumsum(self._probabilities)

        for i in range(n):
            r = rng.random()
            idx = np.searchsorted(cumprob, r)
            idx = min(idx, len(self._distributions) - 1)
            samples[i] = self._distributions[idx].sample(1, rng)[0]

        return samples

    def isPrior(self) -> bool:
        """Return True (used for detection by Posterior solver)."""
        return True

    def isPriorDistribution(self) -> bool:
        """Alias for isPrior (used for detection by Posterior solver)."""
        return True

    def getProbabilities(self) -> np.ndarray:
        """Get all probabilities."""
        return self._probabilities.copy()

    # Snake_case aliases
    get_num_alternatives = getNumAlternatives
    get_alternative = getAlternative
    get_probability = getProbability
    get_probabilities = getProbabilities
    is_prior = isPrior
    is_prior_distribution = isPriorDistribution

