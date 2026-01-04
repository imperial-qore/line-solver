"""
Discrete probability distributions for LINE (pure Python).

This module provides discrete distribution implementations including
Poisson, geometric, binomial, and other common distributions.
"""

from typing import Optional, Tuple, Union
import numpy as np
from scipy import stats
from scipy.special import factorial

from .base import DiscreteDistribution


class Poisson(DiscreteDistribution):
    """
    Poisson distribution.

    The Poisson distribution models the number of events occurring
    in a fixed interval of time or space.

    Args:
        lambda_: Rate parameter (mean = variance = lambda).
    """

    def __init__(self, lambda_: float):
        super().__init__()
        self._name = 'Poisson'
        if lambda_ <= 0:
            raise ValueError("Lambda must be positive")
        self._lambda = lambda_

    @property
    def lambda_(self) -> float:
        """Get the rate parameter."""
        return self._lambda

    def getMean(self) -> float:
        """Get the mean (equals lambda)."""
        return self._lambda

    def getVar(self) -> float:
        """Get the variance (equals lambda)."""
        return self._lambda

    def getSCV(self) -> float:
        """Get the SCV (1/lambda)."""
        return 1.0 / self._lambda

    def getSkew(self) -> float:
        """Get the skewness."""
        return 1.0 / np.sqrt(self._lambda)

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [0, inf)."""
        return (0.0, float('inf'))

    def evalPMF(self, x: int) -> float:
        """Evaluate the probability mass function at point x."""
        if x < 0:
            return 0.0
        return stats.poisson.pmf(x, self._lambda)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < 0:
            return 0.0
        return stats.poisson.cdf(int(x), self._lambda)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.poisson(lam=self._lambda, size=n)



class Geometric(DiscreteDistribution):
    """
    Geometric distribution.

    The geometric distribution models the number of failures before
    the first success in a sequence of Bernoulli trials.

    Args:
        p: Probability of success on each trial (0 < p <= 1).
    """

    def __init__(self, p: float):
        super().__init__()
        self._name = 'Geometric'
        if p <= 0 or p > 1:
            raise ValueError("Probability p must be in (0, 1]")
        self._p = p

    @property
    def p(self) -> float:
        """Get the probability of success."""
        return self._p

    def getMean(self) -> float:
        """Get the mean (1/p for number of trials)."""
        return 1.0 / self._p

    def getVar(self) -> float:
        """Get the variance."""
        return (1 - self._p) / (self._p ** 2)

    def getSkew(self) -> float:
        """Get the skewness."""
        return (2 - self._p) / np.sqrt(1 - self._p)

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [1, inf)."""
        return (1.0, float('inf'))

    def evalPMF(self, x: int) -> float:
        """Evaluate the PMF at point x (number of trials until first success)."""
        if x < 1:
            return 0.0
        return self._p * (1 - self._p) ** (x - 1)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < 1:
            return 0.0
        return 1 - (1 - self._p) ** int(x)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.geometric(p=self._p, size=n)



class Binomial(DiscreteDistribution):
    """
    Binomial distribution.

    The binomial distribution models the number of successes in
    n independent Bernoulli trials.

    Args:
        n: Number of trials.
        p: Probability of success on each trial.
    """

    def __init__(self, n: int, p: float):
        super().__init__()
        self._name = 'Binomial'
        if n < 0:
            raise ValueError("Number of trials must be non-negative")
        if p < 0 or p > 1:
            raise ValueError("Probability p must be in [0, 1]")
        self._n = n
        self._p = p

    @property
    def n(self) -> int:
        """Get the number of trials."""
        return self._n

    @property
    def p(self) -> float:
        """Get the probability of success."""
        return self._p

    def getMean(self) -> float:
        """Get the mean (n * p)."""
        return self._n * self._p

    def getVar(self) -> float:
        """Get the variance (n * p * (1 - p))."""
        return self._n * self._p * (1 - self._p)

    def getSkew(self) -> float:
        """Get the skewness."""
        if self._n == 0 or self._p == 0 or self._p == 1:
            return 0.0
        return (1 - 2 * self._p) / np.sqrt(self._n * self._p * (1 - self._p))

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [0, n]."""
        return (0.0, float(self._n))

    def evalPMF(self, x: int) -> float:
        """Evaluate the PMF at point x."""
        if x < 0 or x > self._n:
            return 0.0
        return stats.binom.pmf(x, self._n, self._p)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < 0:
            return 0.0
        if x >= self._n:
            return 1.0
        return stats.binom.cdf(int(x), self._n, self._p)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.binomial(n=self._n, p=self._p, size=n)



class NegBinomial(DiscreteDistribution):
    """
    Negative Binomial distribution.

    Models the number of failures before r successes in
    a sequence of Bernoulli trials.

    Args:
        r: Number of successes required.
        p: Probability of success on each trial.
    """

    def __init__(self, r: int, p: float):
        super().__init__()
        self._name = 'NegBinomial'
        if r < 1:
            raise ValueError("Number of successes r must be at least 1")
        if p <= 0 or p > 1:
            raise ValueError("Probability p must be in (0, 1]")
        self._r = r
        self._p = p

    @property
    def r(self) -> int:
        """Get the number of successes required."""
        return self._r

    @property
    def p(self) -> float:
        """Get the probability of success."""
        return self._p

    def getMean(self) -> float:
        """Get the mean (r * (1 - p) / p)."""
        return self._r * (1 - self._p) / self._p

    def getVar(self) -> float:
        """Get the variance."""
        return self._r * (1 - self._p) / (self._p ** 2)

    def getSkew(self) -> float:
        """Get the skewness."""
        return (2 - self._p) / np.sqrt(self._r * (1 - self._p))

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [0, inf)."""
        return (0.0, float('inf'))

    def evalPMF(self, x: int) -> float:
        """Evaluate the PMF at point x."""
        if x < 0:
            return 0.0
        return stats.nbinom.pmf(x, self._r, self._p)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < 0:
            return 0.0
        return stats.nbinom.cdf(int(x), self._r, self._p)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.negative_binomial(n=self._r, p=self._p, size=n)



class Zipf(DiscreteDistribution):
    """
    Zipf distribution.

    The Zipf distribution is a power-law distribution often used
    to model rank-frequency relationships (e.g., word frequencies).

    Args:
        s: Shape parameter (s > 0).
        n: Upper bound on support (optional, defaults to large value).
    """

    def __init__(self, s: float, n: int = 10000):
        super().__init__()
        self._name = 'Zipf'
        if s <= 0:
            raise ValueError("Shape parameter s must be > 0")
        if n < 1:
            raise ValueError("Upper bound n must be at least 1")
        self._s = s
        self._n = n
        # Compute normalization constant
        self._H = np.sum(1.0 / np.arange(1, n + 1) ** s)

    @property
    def s(self) -> float:
        """Get the shape parameter."""
        return self._s

    @property
    def n(self) -> int:
        """Get the upper bound."""
        return self._n

    def getMean(self) -> float:
        """Get the mean."""
        k = np.arange(1, self._n + 1)
        return np.sum(k * (1.0 / k ** self._s) / self._H)

    def getVar(self) -> float:
        """Get the variance."""
        k = np.arange(1, self._n + 1)
        probs = (1.0 / k ** self._s) / self._H
        mean = np.sum(k * probs)
        return np.sum(k ** 2 * probs) - mean ** 2

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [1, n]."""
        return (1.0, float(self._n))

    def evalPMF(self, x: int) -> float:
        """Evaluate the PMF at point x."""
        if x < 1 or x > self._n:
            return 0.0
        return (1.0 / x ** self._s) / self._H

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < 1:
            return 0.0
        if x >= self._n:
            return 1.0
        k = int(x)
        return np.sum(1.0 / np.arange(1, k + 1) ** self._s) / self._H

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples using inverse transform."""
        if rng is None:
            rng = np.random.default_rng()
        # Build CDF
        k = np.arange(1, self._n + 1)
        pmf = (1.0 / k ** self._s) / self._H
        cdf = np.cumsum(pmf)
        # Inverse transform sampling
        u = rng.random(size=n)
        return np.searchsorted(cdf, u) + 1



class Empirical(DiscreteDistribution):
    """
    Empirical discrete distribution from data.

    Creates a distribution from observed data by computing
    the empirical probability mass function.

    Args:
        data: Array of observed values.
    """

    def __init__(self, data: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'Empirical'
        self._data = np.array(data)
        # Compute unique values and their probabilities
        unique, counts = np.unique(self._data, return_counts=True)
        self._values = unique
        self._probs = counts / len(self._data)
        self._cdf = np.cumsum(self._probs)

    @property
    def values(self) -> np.ndarray:
        """Get the unique values."""
        return self._values

    @property
    def probs(self) -> np.ndarray:
        """Get the probabilities for each value."""
        return self._probs

    def getMean(self) -> float:
        """Get the mean."""
        return float(np.sum(self._values * self._probs))

    def getVar(self) -> float:
        """Get the variance."""
        mean = self.getMean()
        return float(np.sum(self._probs * (self._values - mean) ** 2))

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [min, max]."""
        return (float(np.min(self._values)), float(np.max(self._values)))

    def evalPMF(self, x: int) -> float:
        """Evaluate the PMF at point x."""
        idx = np.where(self._values == x)[0]
        if len(idx) == 0:
            return 0.0
        return self._probs[idx[0]]

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        idx = np.searchsorted(self._values, x, side='right') - 1
        if idx < 0:
            return 0.0
        return self._cdf[min(idx, len(self._cdf) - 1)]

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        indices = rng.choice(len(self._values), size=n, p=self._probs)
        return self._values[indices]



class Bernoulli(DiscreteDistribution):
    """
    Bernoulli distribution.

    The simplest discrete distribution with two outcomes: success (1)
    with probability p and failure (0) with probability 1-p.

    Args:
        p: Probability of success (0 <= p <= 1).
    """

    def __init__(self, p: float):
        super().__init__()
        self._name = 'Bernoulli'
        if p < 0 or p > 1:
            raise ValueError("Probability p must be in [0, 1]")
        self._p = p

    @property
    def p(self) -> float:
        """Get the probability of success."""
        return self._p

    def getMean(self) -> float:
        """Get the mean (equals p)."""
        return self._p

    def getVar(self) -> float:
        """Get the variance (p * (1 - p))."""
        return self._p * (1 - self._p)

    def getSkew(self) -> float:
        """Get the skewness."""
        if self._p == 0 or self._p == 1:
            return 0.0
        return (1 - 2 * self._p) / np.sqrt(self._p * (1 - self._p))

    def getSupport(self) -> Tuple[float, float]:
        """Get the support {0, 1}."""
        return (0.0, 1.0)

    def evalPMF(self, x: int) -> float:
        """Evaluate the PMF at point x."""
        if x == 0:
            return 1 - self._p
        elif x == 1:
            return self._p
        else:
            return 0.0

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < 0:
            return 0.0
        elif x < 1:
            return 1 - self._p
        else:
            return 1.0

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return (rng.random(size=n) < self._p).astype(int)



class DiscreteUniform(DiscreteDistribution):
    """
    Discrete Uniform distribution.

    Assigns equal probability to all integers in [a, b].

    Args:
        a: Lower bound (inclusive).
        b: Upper bound (inclusive).
    """

    def __init__(self, a: int, b: int):
        super().__init__()
        self._name = 'DiscreteUniform'
        if a > b:
            raise ValueError("Lower bound a must be <= upper bound b")
        self._a = int(a)
        self._b = int(b)
        self._n = self._b - self._a + 1  # Number of values

    @property
    def a(self) -> int:
        """Get the lower bound."""
        return self._a

    @property
    def b(self) -> int:
        """Get the upper bound."""
        return self._b

    def getMean(self) -> float:
        """Get the mean ((a + b) / 2)."""
        return (self._a + self._b) / 2.0

    def getVar(self) -> float:
        """Get the variance ((n^2 - 1) / 12)."""
        return (self._n ** 2 - 1) / 12.0

    def getSkew(self) -> float:
        """Get the skewness (0 for symmetric distribution)."""
        return 0.0

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [a, b]."""
        return (float(self._a), float(self._b))

    def evalPMF(self, x: int) -> float:
        """Evaluate the PMF at point x."""
        if x < self._a or x > self._b:
            return 0.0
        return 1.0 / self._n

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < self._a:
            return 0.0
        if x >= self._b:
            return 1.0
        return (int(x) - self._a + 1) / self._n

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.integers(low=self._a, high=self._b + 1, size=n)



class DiscreteSampler(DiscreteDistribution):
    """
    General discrete distribution from values and probabilities.

    Allows specifying an arbitrary discrete distribution by providing
    the possible values and their probabilities.

    Supports two calling conventions (matching MATLAB API):
        - DiscreteSampler(probs): Values are implicitly 1, 2, ..., n (MATLAB-style)
        - DiscreteSampler(values, probs): Explicit values and probabilities

    Args:
        probs_or_values: If single argument: array of probabilities.
                         If two arguments: array of possible values.
        probs: Array of probabilities (only when first arg is values).
    """

    def __init__(self, probs_or_values: Union[list, np.ndarray],
                 probs: Optional[Union[list, np.ndarray]] = None):
        super().__init__()
        self._name = 'DiscreteSampler'

        if probs is None:
            # MATLAB-style: DiscreteSampler(probs) - values are 1, 2, ..., n
            self._probs = np.array(probs_or_values, dtype=float)
            n = len(self._probs)
            self._values = np.arange(1, n + 1, dtype=float)  # 1-indexed like MATLAB
        else:
            # Full form: DiscreteSampler(values, probs)
            self._values = np.array(probs_or_values, dtype=float)
            self._probs = np.array(probs, dtype=float)

        if len(self._values) != len(self._probs):
            raise ValueError("values and probs must have the same length")

        if len(self._values) == 0:
            raise ValueError("Must have at least one value")

        if np.any(self._probs < 0):
            raise ValueError("Probabilities must be non-negative")

        # Normalize probabilities
        total = self._probs.sum()
        if total <= 0:
            raise ValueError("Probabilities must sum to a positive value")
        self._probs = self._probs / total

        # Sort by values for CDF computation
        sort_idx = np.argsort(self._values)
        self._values = self._values[sort_idx]
        self._probs = self._probs[sort_idx]
        self._cdf = np.cumsum(self._probs)

    @property
    def values(self) -> np.ndarray:
        """Get the possible values."""
        return self._values.copy()

    @property
    def probs(self) -> np.ndarray:
        """Get the probabilities."""
        return self._probs.copy()

    def getMean(self) -> float:
        """Get the mean."""
        return float(np.sum(self._values * self._probs))

    def getVar(self) -> float:
        """Get the variance."""
        mean = self.getMean()
        return float(np.sum(self._probs * (self._values - mean) ** 2))

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [min, max]."""
        return (float(self._values[0]), float(self._values[-1]))

    def evalPMF(self, x: float) -> float:
        """Evaluate the PMF at point x."""
        idx = np.where(np.isclose(self._values, x))[0]
        if len(idx) == 0:
            return 0.0
        return float(self._probs[idx[0]])

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < self._values[0]:
            return 0.0
        idx = np.searchsorted(self._values, x, side='right') - 1
        return float(self._cdf[min(idx, len(self._cdf) - 1)])

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        indices = rng.choice(len(self._values), size=n, p=self._probs)
        return self._values[indices]



class EmpiricalCDF(DiscreteDistribution):
    """
    Empirical distribution from CDF data.

    Creates a distribution from specified CDF points (values and
    cumulative probabilities).

    Args:
        values: Sorted array of values.
        cdf: Cumulative probabilities at each value (must end with 1.0).
    """

    def __init__(self, values: Union[list, np.ndarray],
                 cdf: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'EmpiricalCDF'

        self._values = np.array(values, dtype=float)
        self._cdf = np.array(cdf, dtype=float)

        if len(self._values) != len(self._cdf):
            raise ValueError("values and cdf must have the same length")

        # Verify CDF is monotonically increasing
        if not np.all(np.diff(self._cdf) >= -1e-10):
            raise ValueError("CDF must be monotonically non-decreasing")

        # Normalize CDF to end at 1.0
        if self._cdf[-1] > 0:
            self._cdf = self._cdf / self._cdf[-1]

        # Compute PMF from CDF
        self._probs = np.diff(np.concatenate([[0], self._cdf]))

    @property
    def values(self) -> np.ndarray:
        """Get the values."""
        return self._values.copy()

    @property
    def cdf_values(self) -> np.ndarray:
        """Get the CDF values."""
        return self._cdf.copy()

    def getMean(self) -> float:
        """Get the mean."""
        return float(np.sum(self._values * self._probs))

    def getVar(self) -> float:
        """Get the variance."""
        mean = self.getMean()
        return float(np.sum(self._probs * (self._values - mean) ** 2))

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [min, max]."""
        return (float(self._values[0]), float(self._values[-1]))

    def evalPMF(self, x: float) -> float:
        """Evaluate the PMF at point x."""
        idx = np.where(np.isclose(self._values, x))[0]
        if len(idx) == 0:
            return 0.0
        return float(self._probs[idx[0]])

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < self._values[0]:
            return 0.0
        idx = np.searchsorted(self._values, x, side='right') - 1
        return float(self._cdf[min(idx, len(self._cdf) - 1)])

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples using inverse transform."""
        if rng is None:
            rng = np.random.default_rng()
        u = rng.random(size=n)
        indices = np.searchsorted(self._cdf, u)
        return self._values[np.minimum(indices, len(self._values) - 1)]



class Replayer(DiscreteDistribution):
    """
    Trace-based distribution that replays recorded values.

    Used for simulation where inter-arrival or service times are
    read from a trace file or data array.

    Args:
        trace: Array of values to replay.
        loop: Whether to loop when trace is exhausted (default: True).
    """

    def __init__(self, trace: Union[str, list, np.ndarray], loop: bool = True):
        super().__init__()
        self._name = 'Replayer'

        # Handle file path string
        if isinstance(trace, str):
            # Read trace from file
            self._trace = np.loadtxt(trace, dtype=float)
        else:
            self._trace = np.array(trace, dtype=float)

        if len(self._trace) == 0:
            raise ValueError("Trace must have at least one value")

        self._loop = loop
        self._index = 0

    @property
    def trace(self) -> np.ndarray:
        """Get the trace data."""
        return self._trace.copy()

    @property
    def loop(self) -> bool:
        """Check if looping is enabled."""
        return self._loop

    def reset(self) -> None:
        """Reset the replay index to the beginning."""
        self._index = 0

    def getMean(self) -> float:
        """Get the mean of the trace."""
        return float(np.mean(self._trace))

    def getVar(self) -> float:
        """Get the variance of the trace."""
        return float(np.var(self._trace))

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [min, max] of trace values."""
        return (float(np.min(self._trace)), float(np.max(self._trace)))

    def evalPMF(self, x: float) -> float:
        """Evaluate PMF based on trace frequency."""
        count = np.sum(np.isclose(self._trace, x))
        return count / len(self._trace)

    def evalCDF(self, x: float) -> float:
        """Evaluate CDF based on trace."""
        return np.sum(self._trace <= x) / len(self._trace)

    def next_value(self) -> float:
        """Get the next value in the trace."""
        if self._index >= len(self._trace):
            if self._loop:
                self._index = 0
            else:
                raise StopIteration("Trace exhausted and loop is disabled")

        value = self._trace[self._index]
        self._index += 1
        return float(value)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Get the next n values from the trace."""
        samples = np.zeros(n)
        for i in range(n):
            samples[i] = self.next_value()
        return samples

    def fit_aph(self, order: int = 2):
        """
        Fit an acyclic phase-type (APH) distribution to the trace data.

        Args:
            order: Number of phases for the APH distribution (default: 2).

        Returns:
            APH distribution fitted to the trace data.
        """
        from .markovian import APH

        mean = self.getMean()
        var = self.getVar()
        scv = var / (mean ** 2) if mean > 0 else 1.0

        # APH uses moment matching - just pass mean and scv
        return APH(mean=mean, scv=scv)


class Trace(Replayer):
    """
    Empirical time series from a trace file.

    Alias for Replayer with additional moment computation for
    histogram-style trace data (value, count pairs).

    Args:
        data: Array of values to replay, or 2-column array of (x, cdf) pairs.
        loop: Whether to loop when trace is exhausted (default: True).
    """

    def __init__(self, data: Union[list, np.ndarray], loop: bool = True):
        data = np.atleast_2d(np.array(data, dtype=float))

        # If 2-column data (x, cdf format), extract values
        if data.shape[1] == 2:
            self._histogram_data = data.copy()
            # Convert CDF to samples (approximate)
            values = data[:, 0]
            super().__init__(values, loop)
        else:
            self._histogram_data = None
            super().__init__(data.flatten(), loop)

        self._name = 'Trace'

    def getMoments(self) -> Tuple[float, float, float, float, float]:
        """
        Compute moments from histogram-style trace data.

        Returns:
            Tuple of (m1, m2, m3, scv, skew) - first three moments, SCV, and skewness.
        """
        if self._histogram_data is not None:
            data = self._histogram_data
            m1 = 0.0  # First moment
            m2 = 0.0  # Second moment
            m3 = 0.0  # Third moment

            for i in range(len(data) - 1):
                # Trapezoidal integration
                mid_val = (data[i + 1, 1] - data[i, 1]) / 2 + data[i, 1]
                width = data[i + 1, 0] - data[i, 0]
                m1 += mid_val * width
                m2 += mid_val ** 2 * width
                m3 += mid_val ** 3 * width

            if m1 > 0:
                scv = (m2 / m1 ** 2) - 1
                var = m2 - m1 ** 2
                if var > 0:
                    skew = (m3 - 3 * m1 * var - m1 ** 3) / (var ** 1.5)
                else:
                    skew = 0.0
            else:
                scv = 0.0
                skew = 0.0

            return (m1, m2, m3, scv, skew)
        else:
            # Compute from raw trace data
            m1 = float(np.mean(self._trace))
            m2 = float(np.mean(self._trace ** 2))
            m3 = float(np.mean(self._trace ** 3))
            var = m2 - m1 ** 2
            scv = var / m1 ** 2 if m1 > 0 else 0.0
            skew = (m3 - 3 * m1 * var - m1 ** 3) / (var ** 1.5) if var > 0 else 0.0
            return (m1, m2, m3, scv, skew)

