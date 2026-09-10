"""
Pure Python probability distributions for LINE queueing network models.

This package provides all probability distributions without requiring
pure Python without Java dependencies
and converted to Java only when needed (e.g., when running a solver


Usage:
    from line_solver.distributions import Exp, Erlang, HyperExp
    from line_solver.distributions import Poisson, Geometric
    from line_solver.distributions import PH, MAP, Coxian

    # Create distributions in pure Python
    service = Exp(rate=2.0)
    arrival = Erlang(mean=1.0, phases=3)

    # Distributions can be converted to Java when needed
    """

# Base classes
from .base import (
    Distribution,
    ContinuousDistribution,
    DiscreteDistribution,
    Markovian,
)

# Continuous distributions
from .continuous import (
    Exp,
    Det,
    Immediate,
    Disabled,
    Erlang,
    HyperExp,
    Gamma,
    Lognormal,
    Pareto,
    Uniform,
    Weibull,
    Normal,
    MultivariateNormal,
    Prior,
    Expolynomial,
    NHPP,
    MAPt,
    PHt,
)

# Discrete distributions
from .discrete import (
    Poisson,
    Geometric,
    Binomial,
    NegBinomial,
    Zipf,
    Empirical,
    Bernoulli,
    DiscreteUniform,
    DiscreteSampler,
    EmpiricalCdf,
    EmpiricalCDF,
    Replayer,
    Trace,
)

# Markovian distributions
from .markovian import (
    PH,
    APH,
    Coxian,
    Cox2,
    MAP,
    MMPP,
    MMPP2,
    DMAP,
    ME,
    CME,
    fit_me_mean_scv,
    MarkedMAP,
    MarkedMMPP,
    BMAP,
    RAP,
    MMDP,
    MMDP2,
)

# Rate scaling
from .scaling import dist_scale_rate, distScaleRate

__all__ = [
    # Base
    'Distribution',
    'ContinuousDistribution',
    'DiscreteDistribution',
    'Markovian',
    # Continuous
    'Exp',
    'Det',
    'Immediate',
    'Disabled',
    'Erlang',
    'HyperExp',
    'Gamma',
    'Lognormal',
    'Pareto',
    'Uniform',
    'Weibull',
    'Normal',
    'MultivariateNormal',
    'Prior',
    'Expolynomial',
    'NHPP',
    'MAPt',
    'PHt',
    # Discrete
    'Poisson',
    'Geometric',
    'Binomial',
    'NegBinomial',
    'Zipf',
    'Empirical',
    'Bernoulli',
    'DiscreteUniform',
    'DiscreteSampler',
    'EmpiricalCdf',
    'EmpiricalCDF',
    'Replayer',
    'Trace',
    # Markovian
    'PH',
    'APH',
    'Coxian',
    'Cox2',
    'MAP',
    'MMPP',
    'MMPP2',
    'DMAP',
    'ME',
    'CME',
    'fit_me_mean_scv',
    'MarkedMAP',
    'MarkedMMPP',
    'BMAP',
    'RAP',
    'MMDP',
    'MMDP2',
    # Scaling
    'dist_scale_rate',
    'distScaleRate',
]
