"""
Discrete distribution samplers for DES simulation.

This module provides samplers for discrete probability distributions
used in queueing network simulations, particularly for:
- Batch sizes (Poisson, Binomial, Geometric)
- Discrete arrival patterns
- Integer-valued random variables
"""

import numpy as np
from typing import Optional, List, Union
from .base import DistributionSampler


class PoissonSampler(DistributionSampler):
    """
    Poisson distribution sampler.

    Used for modeling:
    - Number of arrivals in a fixed time interval
    - Batch sizes with Poisson-distributed counts
    - Discrete count data

    Parameters:
        lambda_: Rate parameter (mean = variance = lambda)

    Properties:
        - P(X = k) = (lambda^k * e^(-lambda)) / k!
        - E[X] = Var[X] = lambda
        - Support: {0, 1, 2, ...}
    """

    def __init__(self, lambda_: float, rng: Optional[np.random.Generator] = None):
        """
        Initialize Poisson sampler.

        Args:
            lambda_: Rate parameter (must be > 0)
            rng: Random number generator
        """
        super().__init__(rng)
        if lambda_ <= 0:
            raise ValueError("Poisson lambda must be positive")
        self.lambda_ = lambda_

    def sample(self) -> float:
        """Generate a Poisson random variate."""
        return float(self.rng.poisson(self.lambda_))

    def sample_n(self, n: int) -> np.ndarray:
        """Generate n Poisson random variates."""
        return self.rng.poisson(self.lambda_, size=n).astype(float)

    @property
    def mean(self) -> float:
        """Return the mean of the distribution."""
        return self.lambda_

    @property
    def variance(self) -> float:
        """Return the variance of the distribution."""
        return self.lambda_

    def pmf(self, k: int) -> float:
        """Probability mass function P(X = k)."""
        if k < 0:
            return 0.0
        from scipy.stats import poisson
        return poisson.pmf(k, self.lambda_)


class BinomialSampler(DistributionSampler):
    """
    Binomial distribution sampler.

    Used for modeling:
    - Number of successes in n independent trials
    - Batch sizes with fixed maximum
    - Discrete counts with upper bound

    Parameters:
        n: Number of trials
        p: Success probability per trial

    Properties:
        - P(X = k) = C(n,k) * p^k * (1-p)^(n-k)
        - E[X] = n * p
        - Var[X] = n * p * (1-p)
        - Support: {0, 1, 2, ..., n}
    """

    def __init__(self, n: int, p: float, rng: Optional[np.random.Generator] = None):
        """
        Initialize Binomial sampler.

        Args:
            n: Number of trials (must be >= 1)
            p: Success probability (must be in [0, 1])
            rng: Random number generator
        """
        super().__init__(rng)
        if n < 1:
            raise ValueError("Binomial n must be >= 1")
        if not 0 <= p <= 1:
            raise ValueError("Binomial p must be in [0, 1]")
        self.n = n
        self.p = p

    def sample(self) -> float:
        """Generate a Binomial random variate."""
        return float(self.rng.binomial(self.n, self.p))

    def sample_n(self, count: int) -> np.ndarray:
        """Generate count Binomial random variates."""
        return self.rng.binomial(self.n, self.p, size=count).astype(float)

    @property
    def mean(self) -> float:
        """Return the mean of the distribution."""
        return self.n * self.p

    @property
    def variance(self) -> float:
        """Return the variance of the distribution."""
        return self.n * self.p * (1 - self.p)

    def pmf(self, k: int) -> float:
        """Probability mass function P(X = k)."""
        if k < 0 or k > self.n:
            return 0.0
        from scipy.stats import binom
        return binom.pmf(k, self.n, self.p)


class GeometricSampler(DistributionSampler):
    """
    Geometric distribution sampler.

    Used for modeling:
    - Number of trials until first success
    - Number of failures before first success
    - Memoryless discrete waiting times

    Parameters:
        p: Success probability per trial

    Properties:
        - P(X = k) = (1-p)^(k-1) * p  (for k = 1, 2, 3, ...)
        - E[X] = 1/p
        - Var[X] = (1-p)/p^2
        - Support: {1, 2, 3, ...} (number of trials until success)

    Note:
        This uses the "shifted" convention where X >= 1.
        For the "unshifted" version (X >= 0, counting failures), use shift=False.
    """

    def __init__(self, p: float, shift: bool = True,
                 rng: Optional[np.random.Generator] = None):
        """
        Initialize Geometric sampler.

        Args:
            p: Success probability (must be in (0, 1])
            shift: If True, returns trials until success (X >= 1).
                   If False, returns failures before success (X >= 0).
            rng: Random number generator
        """
        super().__init__(rng)
        if not 0 < p <= 1:
            raise ValueError("Geometric p must be in (0, 1]")
        self.p = p
        self.shift = shift

    def sample(self) -> float:
        """Generate a Geometric random variate."""
        # NumPy's geometric returns number of failures (X >= 0)
        val = self.rng.geometric(self.p)
        if self.shift:
            return float(val)  # Already 1-indexed in NumPy
        return float(val - 1)  # Convert to 0-indexed

    def sample_n(self, n: int) -> np.ndarray:
        """Generate n Geometric random variates."""
        vals = self.rng.geometric(self.p, size=n)
        if self.shift:
            return vals.astype(float)
        return (vals - 1).astype(float)

    @property
    def mean(self) -> float:
        """Return the mean of the distribution."""
        if self.shift:
            return 1.0 / self.p
        return (1.0 - self.p) / self.p

    @property
    def variance(self) -> float:
        """Return the variance of the distribution."""
        return (1.0 - self.p) / (self.p ** 2)

    def pmf(self, k: int) -> float:
        """Probability mass function P(X = k)."""
        if self.shift:
            if k < 1:
                return 0.0
            return ((1 - self.p) ** (k - 1)) * self.p
        else:
            if k < 0:
                return 0.0
            return ((1 - self.p) ** k) * self.p


class NegativeBinomialSampler(DistributionSampler):
    """
    Negative Binomial distribution sampler.

    Used for modeling:
    - Number of failures before r successes
    - Over-dispersed count data (variance > mean)
    - Compound Poisson processes

    Parameters:
        r: Number of successes (can be non-integer for generalized form)
        p: Success probability per trial

    Properties:
        - E[X] = r(1-p)/p
        - Var[X] = r(1-p)/p^2
        - Support: {0, 1, 2, ...}
    """

    def __init__(self, r: float, p: float,
                 rng: Optional[np.random.Generator] = None):
        """
        Initialize Negative Binomial sampler.

        Args:
            r: Number of successes (must be > 0)
            p: Success probability (must be in (0, 1])
            rng: Random number generator
        """
        super().__init__(rng)
        if r <= 0:
            raise ValueError("Negative Binomial r must be > 0")
        if not 0 < p <= 1:
            raise ValueError("Negative Binomial p must be in (0, 1]")
        self.r = r
        self.p = p

    def sample(self) -> float:
        """Generate a Negative Binomial random variate."""
        return float(self.rng.negative_binomial(self.r, self.p))

    def sample_n(self, n: int) -> np.ndarray:
        """Generate n Negative Binomial random variates."""
        return self.rng.negative_binomial(self.r, self.p, size=n).astype(float)

    @property
    def mean(self) -> float:
        """Return the mean of the distribution."""
        return self.r * (1 - self.p) / self.p

    @property
    def variance(self) -> float:
        """Return the variance of the distribution."""
        return self.r * (1 - self.p) / (self.p ** 2)


class BernoulliSampler(DistributionSampler):
    """
    Bernoulli distribution sampler.

    Used for modeling:
    - Binary outcomes (success/failure)
    - Probabilistic routing decisions
    - Random choices with two outcomes

    Parameters:
        p: Success probability

    Properties:
        - P(X = 1) = p
        - P(X = 0) = 1 - p
        - E[X] = p
        - Var[X] = p(1-p)
        - Support: {0, 1}
    """

    def __init__(self, p: float, rng: Optional[np.random.Generator] = None):
        """
        Initialize Bernoulli sampler.

        Args:
            p: Success probability (must be in [0, 1])
            rng: Random number generator
        """
        super().__init__(rng)
        if not 0 <= p <= 1:
            raise ValueError("Bernoulli p must be in [0, 1]")
        self.p = p

    def sample(self) -> float:
        """Generate a Bernoulli random variate (0 or 1)."""
        return float(self.rng.random() < self.p)

    def sample_n(self, n: int) -> np.ndarray:
        """Generate n Bernoulli random variates."""
        return (self.rng.random(n) < self.p).astype(float)

    @property
    def mean(self) -> float:
        """Return the mean of the distribution."""
        return self.p

    @property
    def variance(self) -> float:
        """Return the variance of the distribution."""
        return self.p * (1 - self.p)


class DiscreteUniformSampler(DistributionSampler):
    """
    Discrete Uniform distribution sampler.

    Used for modeling:
    - Equally likely discrete outcomes
    - Random selection from a finite set
    - Uniform batch sizes within a range

    Parameters:
        low: Minimum value (inclusive)
        high: Maximum value (inclusive)

    Properties:
        - P(X = k) = 1/(high - low + 1) for k in {low, ..., high}
        - E[X] = (low + high) / 2
        - Var[X] = ((high - low + 1)^2 - 1) / 12
        - Support: {low, low+1, ..., high}
    """

    def __init__(self, low: int, high: int,
                 rng: Optional[np.random.Generator] = None):
        """
        Initialize Discrete Uniform sampler.

        Args:
            low: Minimum value (inclusive)
            high: Maximum value (inclusive, must be >= low)
            rng: Random number generator
        """
        super().__init__(rng)
        if high < low:
            raise ValueError("DiscreteUniform high must be >= low")
        self.low = low
        self.high = high

    def sample(self) -> float:
        """Generate a Discrete Uniform random variate."""
        return float(self.rng.integers(self.low, self.high + 1))

    def sample_n(self, n: int) -> np.ndarray:
        """Generate n Discrete Uniform random variates."""
        return self.rng.integers(self.low, self.high + 1, size=n).astype(float)

    @property
    def mean(self) -> float:
        """Return the mean of the distribution."""
        return (self.low + self.high) / 2.0

    @property
    def variance(self) -> float:
        """Return the variance of the distribution."""
        n = self.high - self.low + 1
        return (n ** 2 - 1) / 12.0


class ZipfSampler(DistributionSampler):
    """
    Zipf (Zeta) distribution sampler.

    Used for modeling:
    - Word frequency distributions
    - Popularity rankings (websites, content)
    - Power-law phenomena

    Parameters:
        s: Shape parameter (exponent, must be > 0; > 1 for untruncated)
        n: Upper bound (optional, for truncated Zipf)

    Properties:
        - P(X = k) ∝ 1/k^s for k = 1, 2, 3, ...
        - Heavy-tailed distribution
        - Support: {1, 2, 3, ...} or {1, 2, ..., n} if truncated
    """

    def __init__(self, s: float, n: Optional[int] = None,
                 rng: Optional[np.random.Generator] = None):
        """
        Initialize Zipf sampler.

        Args:
            s: Shape parameter (must be > 0; > 1 for untruncated Zipf)
            n: Optional upper bound for truncated Zipf
            rng: Random number generator
        """
        super().__init__(rng)
        if n is None and s <= 1:
            raise ValueError("Zipf s must be > 1 for untruncated Zipf")
        if s <= 0:
            raise ValueError("Zipf s must be > 0")
        self.s = s
        self.n = n

        # Pre-compute CDF for truncated version
        if n is not None:
            self._precompute_cdf()

    def _precompute_cdf(self):
        """Pre-compute CDF for truncated Zipf sampling."""
        k = np.arange(1, self.n + 1)
        pmf = 1.0 / (k ** self.s)
        pmf /= pmf.sum()
        self.cdf = np.cumsum(pmf)
        self.pmf_values = pmf

    def sample(self) -> float:
        """Generate a Zipf random variate."""
        if self.n is not None:
            # Truncated Zipf: use inverse CDF
            u = self.rng.random()
            return float(np.searchsorted(self.cdf, u) + 1)
        else:
            # Standard Zipf (uses NumPy's implementation)
            return float(self.rng.zipf(self.s))

    def sample_n(self, count: int) -> np.ndarray:
        """Generate count Zipf random variates."""
        if self.n is not None:
            u = self.rng.random(count)
            return (np.searchsorted(self.cdf, u) + 1).astype(float)
        else:
            return self.rng.zipf(self.s, size=count).astype(float)

    @property
    def mean(self) -> float:
        """Return the mean of the distribution (if it exists)."""
        if self.n is not None:
            k = np.arange(1, self.n + 1)
            return float(np.sum(k * self.pmf_values))
        elif self.s > 2:
            from scipy.special import zeta
            return zeta(self.s - 1) / zeta(self.s)
        else:
            return float('inf')


class CategoricalSampler(DistributionSampler):
    """
    Categorical (Generalized Bernoulli) distribution sampler.

    Used for modeling:
    - Multi-class routing decisions
    - Discrete outcomes with arbitrary probabilities
    - Class switching probabilities

    Parameters:
        probabilities: List of probabilities (must sum to 1)
        values: Optional list of values to return (default: 0, 1, 2, ...)

    Properties:
        - P(X = values[i]) = probabilities[i]
        - Support: {values[0], values[1], ...}
    """

    def __init__(self, probabilities: List[float],
                 values: Optional[List[Union[int, float]]] = None,
                 rng: Optional[np.random.Generator] = None):
        """
        Initialize Categorical sampler.

        Args:
            probabilities: List of probabilities (must sum to ~1)
            values: Optional values to return (default: 0, 1, 2, ...)
            rng: Random number generator
        """
        super().__init__(rng)
        probs = np.array(probabilities, dtype=float)
        probs /= probs.sum()  # Normalize

        self.probabilities = probs
        self.cdf = np.cumsum(probs)

        if values is not None:
            if len(values) != len(probabilities):
                raise ValueError("values and probabilities must have same length")
            self.values = np.array(values, dtype=float)
        else:
            self.values = np.arange(len(probabilities), dtype=float)

    def sample(self) -> float:
        """Generate a Categorical random variate."""
        u = self.rng.random()
        idx = np.searchsorted(self.cdf, u)
        return float(self.values[min(idx, len(self.values) - 1)])

    def sample_n(self, n: int) -> np.ndarray:
        """Generate n Categorical random variates."""
        u = self.rng.random(n)
        indices = np.searchsorted(self.cdf, u)
        indices = np.clip(indices, 0, len(self.values) - 1)
        return self.values[indices]

    @property
    def mean(self) -> float:
        """Return the mean of the distribution."""
        return float(np.sum(self.values * self.probabilities))

    @property
    def variance(self) -> float:
        """Return the variance of the distribution."""
        mu = self.mean
        return float(np.sum(self.probabilities * (self.values - mu) ** 2))


class EmpiricalDiscreteSampler(DistributionSampler):
    """
    Empirical discrete distribution sampler.

    Used for modeling:
    - Distributions fitted from observed data
    - Histogram-based sampling
    - Custom discrete distributions

    Parameters:
        values: Observed values
        counts: Frequency counts for each value (optional)
    """

    def __init__(self, values: List[Union[int, float]],
                 counts: Optional[List[int]] = None,
                 rng: Optional[np.random.Generator] = None):
        """
        Initialize Empirical Discrete sampler.

        Args:
            values: List of distinct values
            counts: Frequency counts (optional, default: all equal)
            rng: Random number generator
        """
        super().__init__(rng)
        self.unique_values = np.array(values, dtype=float)

        if counts is not None:
            counts_arr = np.array(counts, dtype=float)
        else:
            counts_arr = np.ones(len(values))

        self.probabilities = counts_arr / counts_arr.sum()
        self.cdf = np.cumsum(self.probabilities)

    def sample(self) -> float:
        """Generate a sample from the empirical distribution."""
        u = self.rng.random()
        idx = np.searchsorted(self.cdf, u)
        return float(self.unique_values[min(idx, len(self.unique_values) - 1)])

    def sample_n(self, n: int) -> np.ndarray:
        """Generate n samples from the empirical distribution."""
        u = self.rng.random(n)
        indices = np.searchsorted(self.cdf, u)
        indices = np.clip(indices, 0, len(self.unique_values) - 1)
        return self.unique_values[indices]

    @property
    def mean(self) -> float:
        """Return the mean of the distribution."""
        return float(np.sum(self.unique_values * self.probabilities))

    @property
    def variance(self) -> float:
        """Return the variance of the distribution."""
        mu = self.mean
        return float(np.sum(self.probabilities * (self.unique_values - mu) ** 2))


# Export all discrete samplers
__all__ = [
    'PoissonSampler',
    'BinomialSampler',
    'GeometricSampler',
    'NegativeBinomialSampler',
    'BernoulliSampler',
    'DiscreteUniformSampler',
    'ZipfSampler',
    'CategoricalSampler',
    'EmpiricalDiscreteSampler',
]
