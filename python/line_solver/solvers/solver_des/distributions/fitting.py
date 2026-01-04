"""
Parameter fitting and estimation for distributions.

This module provides methods for:
- Moment matching (fitting distributions to match mean, variance, etc.)
- Maximum Likelihood Estimation (MLE)
- Method of moments estimation
- Distribution fitting from data samples
"""

from typing import List, Optional, Tuple, Dict, Any, Union
from dataclasses import dataclass
import numpy as np
from scipy import optimize, stats

from .base import DistributionSampler
from .continuous import (
    ExponentialSampler,
    ErlangSampler,
    GammaSampler,
    UniformSampler,
    ParetoSampler,
    WeibullSampler,
    LognormalSampler,
    NormalSampler,
)
from .phase_type import (
    PHSampler,
    APHSampler,
    HyperExponentialSampler,
    HyperExp2Sampler,
    CoxianSampler,
    Coxian2Sampler,
)


@dataclass
class FitResult:
    """
    Result of a distribution fitting operation.

    Attributes:
        distribution: Name of the fitted distribution
        parameters: Dictionary of fitted parameters
        sampler: Configured DistributionSampler instance
        log_likelihood: Log-likelihood of the fit (if computed)
        aic: Akaike Information Criterion (if computed)
        bic: Bayesian Information Criterion (if computed)
        ks_statistic: Kolmogorov-Smirnov test statistic
        ks_pvalue: Kolmogorov-Smirnov p-value
        goodness_of_fit: Summary goodness-of-fit metric
    """
    distribution: str
    parameters: Dict[str, Any]
    sampler: DistributionSampler
    log_likelihood: Optional[float] = None
    aic: Optional[float] = None
    bic: Optional[float] = None
    ks_statistic: Optional[float] = None
    ks_pvalue: Optional[float] = None
    goodness_of_fit: Optional[float] = None


class DistributionFitter:
    """
    Fits distributions to data samples.

    Supports multiple fitting methods:
    - Moment matching (fast, approximate)
    - Maximum Likelihood Estimation (accurate, slower)
    - Automatic distribution selection (fits multiple, picks best)
    """

    def __init__(self, rng: Optional[np.random.Generator] = None):
        """
        Initialize the fitter.

        Args:
            rng: Random number generator for created samplers
        """
        self.rng = rng if rng is not None else np.random.default_rng()

    def fit_exponential(
        self,
        data: Union[List[float], np.ndarray],
        method: str = 'mle'
    ) -> FitResult:
        """
        Fit exponential distribution to data.

        Args:
            data: Sample data
            method: 'mle' or 'mom' (method of moments)

        Returns:
            FitResult with fitted exponential distribution
        """
        data = np.array(data)
        n = len(data)

        # Rate = 1 / mean
        mean = data.mean()
        rate = 1.0 / mean

        sampler = ExponentialSampler(rate=rate, rng=self.rng)

        # Compute log-likelihood
        log_lik = n * np.log(rate) - rate * data.sum()

        # Information criteria
        k = 1  # number of parameters
        aic = 2 * k - 2 * log_lik
        bic = k * np.log(n) - 2 * log_lik

        # KS test
        ks_stat, ks_pval = stats.kstest(data, 'expon', args=(0, 1/rate))

        return FitResult(
            distribution='exponential',
            parameters={'rate': rate, 'mean': mean},
            sampler=sampler,
            log_likelihood=log_lik,
            aic=aic,
            bic=bic,
            ks_statistic=ks_stat,
            ks_pvalue=ks_pval,
            goodness_of_fit=1 - ks_stat
        )

    def fit_erlang(
        self,
        data: Union[List[float], np.ndarray],
        shape: Optional[int] = None
    ) -> FitResult:
        """
        Fit Erlang distribution to data.

        If shape is not specified, estimates it from the data.

        Args:
            data: Sample data
            shape: Fixed shape parameter (optional)

        Returns:
            FitResult with fitted Erlang distribution
        """
        data = np.array(data)
        n = len(data)

        mean = data.mean()
        var = data.var(ddof=1)
        scv = var / (mean ** 2) if mean > 0 else 1.0

        if shape is None:
            # Estimate shape from squared coefficient of variation
            # For Erlang: scv = 1/k, so k = 1/scv
            shape = max(1, int(round(1 / scv)))

        rate = shape / mean

        sampler = ErlangSampler(shape=shape, rate=rate, rng=self.rng)

        # KS test
        ks_stat, ks_pval = stats.kstest(data, 'gamma', args=(shape, 0, 1/rate))

        return FitResult(
            distribution='erlang',
            parameters={'shape': shape, 'rate': rate, 'mean': mean, 'scv': scv},
            sampler=sampler,
            ks_statistic=ks_stat,
            ks_pvalue=ks_pval,
            goodness_of_fit=1 - ks_stat
        )

    def fit_gamma(
        self,
        data: Union[List[float], np.ndarray],
        method: str = 'mle'
    ) -> FitResult:
        """
        Fit gamma distribution to data.

        Args:
            data: Sample data
            method: 'mle' or 'mom'

        Returns:
            FitResult with fitted gamma distribution
        """
        data = np.array(data)
        n = len(data)

        if method == 'mom':
            mean = data.mean()
            var = data.var(ddof=1)
            shape = (mean ** 2) / var
            rate = mean / var
        else:
            # MLE using scipy
            shape, loc, scale = stats.gamma.fit(data, floc=0)
            rate = 1.0 / scale

        sampler = GammaSampler(shape=shape, rate=rate, rng=self.rng)

        # KS test
        ks_stat, ks_pval = stats.kstest(data, 'gamma', args=(shape, 0, 1/rate))

        return FitResult(
            distribution='gamma',
            parameters={'shape': shape, 'rate': rate, 'mean': shape/rate},
            sampler=sampler,
            ks_statistic=ks_stat,
            ks_pvalue=ks_pval,
            goodness_of_fit=1 - ks_stat
        )

    def fit_lognormal(
        self,
        data: Union[List[float], np.ndarray],
        method: str = 'mle'
    ) -> FitResult:
        """
        Fit lognormal distribution to data.

        Args:
            data: Sample data
            method: 'mle' or 'mom'

        Returns:
            FitResult with fitted lognormal distribution
        """
        data = np.array(data)
        data = data[data > 0]  # Lognormal requires positive data

        if method == 'mom':
            log_data = np.log(data)
            mu = log_data.mean()
            sigma = log_data.std(ddof=1)
        else:
            # MLE
            sigma, loc, scale = stats.lognorm.fit(data, floc=0)
            mu = np.log(scale)

        sampler = LognormalSampler(mu=mu, sigma=sigma, rng=self.rng)

        # KS test
        ks_stat, ks_pval = stats.kstest(data, 'lognorm', args=(sigma, 0, np.exp(mu)))

        return FitResult(
            distribution='lognormal',
            parameters={'mu': mu, 'sigma': sigma},
            sampler=sampler,
            ks_statistic=ks_stat,
            ks_pvalue=ks_pval,
            goodness_of_fit=1 - ks_stat
        )

    def fit_weibull(
        self,
        data: Union[List[float], np.ndarray],
        method: str = 'mle'
    ) -> FitResult:
        """
        Fit Weibull distribution to data.

        Args:
            data: Sample data
            method: 'mle' or 'mom'

        Returns:
            FitResult with fitted Weibull distribution
        """
        data = np.array(data)
        data = data[data > 0]

        # MLE using scipy
        c, loc, scale = stats.weibull_min.fit(data, floc=0)

        sampler = WeibullSampler(shape=c, scale=scale, rng=self.rng)

        # KS test
        ks_stat, ks_pval = stats.kstest(data, 'weibull_min', args=(c, 0, scale))

        return FitResult(
            distribution='weibull',
            parameters={'shape': c, 'scale': scale},
            sampler=sampler,
            ks_statistic=ks_stat,
            ks_pvalue=ks_pval,
            goodness_of_fit=1 - ks_stat
        )

    def fit_pareto(
        self,
        data: Union[List[float], np.ndarray]
    ) -> FitResult:
        """
        Fit Pareto distribution to data.

        Args:
            data: Sample data

        Returns:
            FitResult with fitted Pareto distribution
        """
        data = np.array(data)
        data = data[data > 0]

        # MLE for Pareto
        xm = data.min()  # Scale (minimum)
        alpha = len(data) / np.sum(np.log(data / xm))  # Shape

        sampler = ParetoSampler(scale=xm, shape=alpha, rng=self.rng)

        # KS test
        ks_stat, ks_pval = stats.kstest(data, 'pareto', args=(alpha, 0, xm))

        return FitResult(
            distribution='pareto',
            parameters={'scale': xm, 'shape': alpha},
            sampler=sampler,
            ks_statistic=ks_stat,
            ks_pvalue=ks_pval,
            goodness_of_fit=1 - ks_stat
        )

    def fit_hyperexp2(
        self,
        data: Union[List[float], np.ndarray]
    ) -> FitResult:
        """
        Fit 2-phase hyperexponential (HyperExp2) to data.

        Uses moment matching (mean and SCV).

        Args:
            data: Sample data

        Returns:
            FitResult with fitted HyperExp2 distribution
        """
        data = np.array(data)

        mean = data.mean()
        var = data.var(ddof=1)
        scv = var / (mean ** 2)

        if scv <= 1:
            # SCV <= 1: use Erlang or Coxian instead
            return self.fit_erlang(data)

        sampler = HyperExp2Sampler(mean=mean, scv=scv, rng=self.rng)

        return FitResult(
            distribution='hyperexp2',
            parameters={'mean': mean, 'scv': scv, 'variance': var},
            sampler=sampler,
            goodness_of_fit=0.5  # No KS test for moment-matched
        )

    def fit_coxian2(
        self,
        data: Union[List[float], np.ndarray]
    ) -> FitResult:
        """
        Fit 2-phase Coxian distribution to data.

        Uses moment matching (mean and SCV).

        Args:
            data: Sample data

        Returns:
            FitResult with fitted Coxian2 distribution
        """
        data = np.array(data)

        mean = data.mean()
        var = data.var(ddof=1)
        scv = var / (mean ** 2)

        if scv >= 1:
            # SCV >= 1: use HyperExp instead
            return self.fit_hyperexp2(data)

        sampler = Coxian2Sampler(mean=mean, scv=scv, rng=self.rng)

        return FitResult(
            distribution='coxian2',
            parameters={'mean': mean, 'scv': scv, 'variance': var},
            sampler=sampler,
            goodness_of_fit=0.5
        )

    def fit_phase_type(
        self,
        data: Union[List[float], np.ndarray],
        num_moments: int = 3
    ) -> FitResult:
        """
        Fit general phase-type distribution using moment matching.

        Args:
            data: Sample data
            num_moments: Number of moments to match (2 or 3)

        Returns:
            FitResult with fitted PH distribution
        """
        data = np.array(data)

        # Compute moments
        m1 = data.mean()
        m2 = (data ** 2).mean()
        m3 = (data ** 3).mean() if num_moments >= 3 else None

        moments = [m1, m2]
        if m3 is not None:
            moments.append(m3)

        sampler = APHSampler.from_moments(moments=moments, rng=self.rng)

        return FitResult(
            distribution='aph',
            parameters={'moments': moments},
            sampler=sampler,
            goodness_of_fit=0.5
        )

    def auto_fit(
        self,
        data: Union[List[float], np.ndarray],
        candidates: Optional[List[str]] = None
    ) -> FitResult:
        """
        Automatically fit best distribution to data.

        Tries multiple distributions and returns the one with
        best goodness-of-fit (based on KS test or AIC).

        Args:
            data: Sample data
            candidates: List of distribution names to try
                       (default: common service time distributions)

        Returns:
            FitResult with best fitting distribution
        """
        if candidates is None:
            candidates = [
                'exponential', 'erlang', 'gamma',
                'lognormal', 'weibull', 'hyperexp2'
            ]

        data = np.array(data)

        # Compute statistics to help select
        mean = data.mean()
        var = data.var(ddof=1)
        scv = var / (mean ** 2) if mean > 0 else 1.0

        best_fit = None
        best_score = -float('inf')

        for dist in candidates:
            try:
                if dist == 'exponential':
                    fit = self.fit_exponential(data)
                elif dist == 'erlang':
                    fit = self.fit_erlang(data)
                elif dist == 'gamma':
                    fit = self.fit_gamma(data)
                elif dist == 'lognormal':
                    fit = self.fit_lognormal(data)
                elif dist == 'weibull':
                    fit = self.fit_weibull(data)
                elif dist == 'pareto':
                    fit = self.fit_pareto(data)
                elif dist == 'hyperexp2':
                    if scv > 1:
                        fit = self.fit_hyperexp2(data)
                    else:
                        continue
                elif dist == 'coxian2':
                    if scv < 1:
                        fit = self.fit_coxian2(data)
                    else:
                        continue
                else:
                    continue

                score = fit.goodness_of_fit or 0.0
                if score > best_score:
                    best_score = score
                    best_fit = fit

            except Exception:
                continue

        if best_fit is None:
            # Fallback to exponential
            best_fit = self.fit_exponential(data)

        return best_fit

    def fit_from_moments(
        self,
        mean: float,
        scv: float,
        target: str = 'auto'
    ) -> DistributionSampler:
        """
        Create sampler from specified moments (mean and SCV).

        Args:
            mean: Target mean
            scv: Target squared coefficient of variation (variance/mean^2)
            target: Distribution type: 'auto', 'exponential', 'erlang',
                   'hyperexp2', 'coxian2', 'ph'

        Returns:
            Configured DistributionSampler
        """
        if target == 'auto':
            if abs(scv - 1.0) < 0.01:
                target = 'exponential'
            elif scv < 1:
                target = 'erlang'
            else:
                target = 'hyperexp2'

        if target == 'exponential':
            return ExponentialSampler(rate=1.0/mean, rng=self.rng)

        elif target == 'erlang':
            shape = max(1, int(round(1 / scv)))
            rate = shape / mean
            return ErlangSampler(shape=shape, rate=rate, rng=self.rng)

        elif target == 'hyperexp2':
            return HyperExp2Sampler(mean=mean, scv=scv, rng=self.rng)

        elif target == 'coxian2':
            return Coxian2Sampler(mean=mean, scv=scv, rng=self.rng)

        elif target == 'ph':
            variance = scv * (mean ** 2)
            m2 = variance + mean ** 2
            return APHSampler.from_moments([mean, m2], rng=self.rng)

        else:
            raise ValueError(f"Unknown target distribution: {target}")


def fit_distribution(
    data: Union[List[float], np.ndarray],
    distribution: str = 'auto',
    rng: Optional[np.random.Generator] = None,
    **kwargs
) -> FitResult:
    """
    Convenience function to fit a distribution to data.

    Args:
        data: Sample data
        distribution: Distribution type or 'auto'
        rng: Random number generator
        **kwargs: Additional arguments for specific fitters

    Returns:
        FitResult with fitted distribution
    """
    fitter = DistributionFitter(rng)

    if distribution == 'auto':
        return fitter.auto_fit(data, **kwargs)
    elif distribution == 'exponential':
        return fitter.fit_exponential(data, **kwargs)
    elif distribution == 'erlang':
        return fitter.fit_erlang(data, **kwargs)
    elif distribution == 'gamma':
        return fitter.fit_gamma(data, **kwargs)
    elif distribution == 'lognormal':
        return fitter.fit_lognormal(data, **kwargs)
    elif distribution == 'weibull':
        return fitter.fit_weibull(data, **kwargs)
    elif distribution == 'pareto':
        return fitter.fit_pareto(data, **kwargs)
    elif distribution in ['hyperexp', 'hyperexp2']:
        return fitter.fit_hyperexp2(data, **kwargs)
    elif distribution in ['coxian', 'coxian2']:
        return fitter.fit_coxian2(data, **kwargs)
    elif distribution in ['ph', 'aph']:
        return fitter.fit_phase_type(data, **kwargs)
    else:
        raise ValueError(f"Unknown distribution: {distribution}")


def moment_match(
    mean: float,
    scv: float,
    rng: Optional[np.random.Generator] = None
) -> DistributionSampler:
    """
    Create sampler matching specified mean and SCV.

    Args:
        mean: Target mean
        scv: Target squared coefficient of variation
        rng: Random number generator

    Returns:
        DistributionSampler matching the moments
    """
    fitter = DistributionFitter(rng)
    return fitter.fit_from_moments(mean, scv)


def compute_scv(data: Union[List[float], np.ndarray]) -> float:
    """
    Compute squared coefficient of variation from data.

    Args:
        data: Sample data

    Returns:
        SCV = variance / mean^2
    """
    data = np.array(data)
    mean = data.mean()
    if mean <= 0:
        return 1.0
    return data.var(ddof=1) / (mean ** 2)


def compute_moments(
    data: Union[List[float], np.ndarray],
    n: int = 3
) -> List[float]:
    """
    Compute first n moments from data.

    Args:
        data: Sample data
        n: Number of moments

    Returns:
        List of moments [E[X], E[X^2], E[X^3], ...]
    """
    data = np.array(data)
    return [float((data ** (i+1)).mean()) for i in range(n)]
