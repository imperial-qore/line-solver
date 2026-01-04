"""
Continuous distribution samplers.

This module implements standard continuous distributions:
- Exponential
- Erlang
- Gamma
- Uniform
- Deterministic
- Pareto
- Weibull
- Lognormal
- Normal (truncated)
"""

from typing import Optional
import numpy as np
from scipy import stats

from .base import DistributionSampler


class ExponentialSampler(DistributionSampler):
    """
    Exponential distribution sampler.

    f(x) = rate * exp(-rate * x) for x >= 0

    Parameters:
        rate: Rate parameter (lambda = 1/mean)
    """

    def __init__(self, rate: float, rng: Optional[np.random.Generator] = None):
        super().__init__(rng)
        if rate <= 0:
            raise ValueError("Rate must be positive")
        self.rate = rate
        self.scale = 1.0 / rate  # numpy uses scale = 1/rate

    def sample(self) -> float:
        return self.rng.exponential(self.scale)

    def get_mean(self) -> float:
        return self.scale

    def get_variance(self) -> float:
        return self.scale * self.scale

    def get_scv(self) -> float:
        return 1.0  # SCV of exponential is always 1


class ErlangSampler(DistributionSampler):
    """
    Erlang distribution sampler (Gamma with integer shape).

    Sum of k independent exponential distributions with rate lambda.

    Parameters:
        shape: Number of stages (k, must be positive integer)
        rate: Rate parameter per stage
    """

    def __init__(self, shape: int, rate: float, rng: Optional[np.random.Generator] = None):
        super().__init__(rng)
        if shape < 1:
            raise ValueError("Shape must be a positive integer")
        if rate <= 0:
            raise ValueError("Rate must be positive")
        self.shape = int(shape)
        self.rate = rate
        self.scale = 1.0 / rate

    def sample(self) -> float:
        return self.rng.gamma(self.shape, self.scale)

    def get_mean(self) -> float:
        return self.shape * self.scale

    def get_variance(self) -> float:
        return self.shape * self.scale * self.scale

    def get_scv(self) -> float:
        return 1.0 / self.shape


class GammaSampler(DistributionSampler):
    """
    Gamma distribution sampler.

    Parameters:
        shape: Shape parameter (alpha > 0)
        rate: Rate parameter (beta = 1/scale)
    """

    def __init__(self, shape: float, rate: float, rng: Optional[np.random.Generator] = None):
        super().__init__(rng)
        if shape <= 0:
            raise ValueError("Shape must be positive")
        if rate <= 0:
            raise ValueError("Rate must be positive")
        self.shape = shape
        self.rate = rate
        self.scale = 1.0 / rate

    def sample(self) -> float:
        return self.rng.gamma(self.shape, self.scale)

    def get_mean(self) -> float:
        return self.shape * self.scale

    def get_variance(self) -> float:
        return self.shape * self.scale * self.scale

    def get_scv(self) -> float:
        return 1.0 / self.shape


class DeterministicSampler(DistributionSampler):
    """
    Deterministic (constant) sampler.

    Always returns the same value.

    Parameters:
        value: The constant value to return
    """

    def __init__(self, value: float, rng: Optional[np.random.Generator] = None):
        super().__init__(rng)
        if value < 0:
            raise ValueError("Value must be non-negative")
        self.value = value

    def sample(self) -> float:
        return self.value

    def get_mean(self) -> float:
        return self.value

    def get_variance(self) -> float:
        return 0.0

    def get_scv(self) -> float:
        return 0.0


class UniformSampler(DistributionSampler):
    """
    Uniform distribution sampler on [a, b].

    Parameters:
        a: Lower bound
        b: Upper bound
    """

    def __init__(self, a: float, b: float, rng: Optional[np.random.Generator] = None):
        super().__init__(rng)
        if a >= b:
            raise ValueError("Lower bound must be less than upper bound")
        self.a = a
        self.b = b

    def sample(self) -> float:
        return self.rng.uniform(self.a, self.b)

    def get_mean(self) -> float:
        return (self.a + self.b) / 2.0

    def get_variance(self) -> float:
        return (self.b - self.a) ** 2 / 12.0


class ParetoSampler(DistributionSampler):
    """
    Pareto distribution sampler.

    f(x) = alpha * x_m^alpha / x^(alpha+1) for x >= x_m

    Parameters:
        shape: Shape parameter (alpha > 0)
        scale: Scale parameter (x_m > 0, minimum value)
    """

    def __init__(self, shape: float, scale: float, rng: Optional[np.random.Generator] = None):
        super().__init__(rng)
        if shape <= 0:
            raise ValueError("Shape must be positive")
        if scale <= 0:
            raise ValueError("Scale must be positive")
        self.shape = shape
        self.scale = scale

    def sample(self) -> float:
        # Pareto type I: x = x_m * (1-U)^(-1/alpha)
        return self.scale * self.rng.pareto(self.shape) + self.scale

    def get_mean(self) -> float:
        if self.shape <= 1:
            return float('inf')
        return self.shape * self.scale / (self.shape - 1)

    def get_variance(self) -> float:
        if self.shape <= 2:
            return float('inf')
        return (self.scale ** 2 * self.shape) / ((self.shape - 1) ** 2 * (self.shape - 2))


class WeibullSampler(DistributionSampler):
    """
    Weibull distribution sampler.

    f(x) = (k/lambda) * (x/lambda)^(k-1) * exp(-(x/lambda)^k)

    Parameters:
        shape: Shape parameter (k > 0)
        scale: Scale parameter (lambda > 0)
    """

    def __init__(self, shape: float, scale: float, rng: Optional[np.random.Generator] = None):
        super().__init__(rng)
        if shape <= 0:
            raise ValueError("Shape must be positive")
        if scale <= 0:
            raise ValueError("Scale must be positive")
        self.shape = shape
        self.scale = scale

    def sample(self) -> float:
        return self.scale * self.rng.weibull(self.shape)

    def get_mean(self) -> float:
        from scipy.special import gamma as gamma_func
        return self.scale * gamma_func(1 + 1 / self.shape)

    def get_variance(self) -> float:
        from scipy.special import gamma as gamma_func
        g1 = gamma_func(1 + 1 / self.shape)
        g2 = gamma_func(1 + 2 / self.shape)
        return self.scale ** 2 * (g2 - g1 ** 2)


class LognormalSampler(DistributionSampler):
    """
    Lognormal distribution sampler.

    log(X) ~ Normal(mu, sigma^2)

    Parameters:
        mu: Mean of the underlying normal distribution
        sigma: Standard deviation of the underlying normal distribution
    """

    def __init__(self, mu: float, sigma: float, rng: Optional[np.random.Generator] = None):
        super().__init__(rng)
        if sigma <= 0:
            raise ValueError("Sigma must be positive")
        self.mu = mu
        self.sigma = sigma

    def sample(self) -> float:
        return self.rng.lognormal(self.mu, self.sigma)

    def get_mean(self) -> float:
        return np.exp(self.mu + self.sigma ** 2 / 2)

    def get_variance(self) -> float:
        return (np.exp(self.sigma ** 2) - 1) * np.exp(2 * self.mu + self.sigma ** 2)


class NormalSampler(DistributionSampler):
    """
    Normal distribution sampler (truncated at 0).

    For queueing applications, negative values are truncated to 0.

    Parameters:
        mean: Mean of the distribution
        std: Standard deviation
    """

    def __init__(self, mean: float, std: float, rng: Optional[np.random.Generator] = None):
        super().__init__(rng)
        if std <= 0:
            raise ValueError("Standard deviation must be positive")
        self.mean = mean
        self.std = std

    def sample(self) -> float:
        value = self.rng.normal(self.mean, self.std)
        return max(0.0, value)  # Truncate at 0

    def get_mean(self) -> float:
        # This is approximate (ignores truncation effect)
        return self.mean

    def get_variance(self) -> float:
        # This is approximate (ignores truncation effect)
        return self.std ** 2


class BetaSampler(DistributionSampler):
    """
    Beta distribution sampler on [0, 1].

    Parameters:
        alpha: Shape parameter alpha > 0
        beta: Shape parameter beta > 0
    """

    def __init__(self, alpha: float, beta: float, rng: Optional[np.random.Generator] = None):
        super().__init__(rng)
        if alpha <= 0 or beta <= 0:
            raise ValueError("Alpha and beta must be positive")
        self.alpha = alpha
        self.beta = beta

    def sample(self) -> float:
        return self.rng.beta(self.alpha, self.beta)

    def get_mean(self) -> float:
        return self.alpha / (self.alpha + self.beta)

    def get_variance(self) -> float:
        ab = self.alpha + self.beta
        return (self.alpha * self.beta) / (ab ** 2 * (ab + 1))


class ReplayTraceSampler(DistributionSampler):
    """
    Replay a trace of values.

    Cycles through a list of pre-recorded values.

    Parameters:
        values: List of values to replay
        loop: Whether to loop back to start after exhausting values
    """

    def __init__(self, values: list, loop: bool = True, rng: Optional[np.random.Generator] = None):
        super().__init__(rng)
        if not values:
            raise ValueError("Values list cannot be empty")
        self.values = list(values)
        self.loop = loop
        self.index = 0

    def sample(self) -> float:
        if self.index >= len(self.values):
            if self.loop:
                self.index = 0
            else:
                return self.values[-1]  # Return last value if not looping

        value = self.values[self.index]
        self.index += 1
        return value

    def get_mean(self) -> float:
        return np.mean(self.values)

    def get_variance(self) -> float:
        return np.var(self.values)


class MixtureSampler(DistributionSampler):
    """
    Mixture of multiple distributions.

    Samples from one of several distributions according to probabilities.

    Parameters:
        samplers: List of DistributionSampler instances
        probabilities: Mixing probabilities (must sum to 1)
    """

    def __init__(
        self,
        samplers: list,
        probabilities: list,
        rng: Optional[np.random.Generator] = None
    ):
        super().__init__(rng)
        if len(samplers) != len(probabilities):
            raise ValueError("Number of samplers must match number of probabilities")
        if abs(sum(probabilities) - 1.0) > 1e-10:
            raise ValueError("Probabilities must sum to 1")

        self.samplers = samplers
        self.probabilities = np.array(probabilities)
        self.cumulative = np.cumsum(self.probabilities)

    def sample(self) -> float:
        u = self.rng.random()
        idx = np.searchsorted(self.cumulative, u)
        idx = min(idx, len(self.samplers) - 1)
        return self.samplers[idx].sample()

    def get_mean(self) -> float:
        means = [s.get_mean() for s in self.samplers]
        return np.dot(self.probabilities, means)

    def get_variance(self) -> float:
        means = np.array([s.get_mean() for s in self.samplers])
        variances = np.array([s.get_variance() for s in self.samplers])
        # Var(X) = E[Var(X|I)] + Var(E[X|I])
        overall_mean = np.dot(self.probabilities, means)
        within_variance = np.dot(self.probabilities, variances)
        between_variance = np.dot(self.probabilities, (means - overall_mean) ** 2)
        return within_variance + between_variance
