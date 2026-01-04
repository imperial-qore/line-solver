"""
Distribution sampler implementations for DES solver.

This package contains distribution samplers for generating random variates
used in arrival and service time generation.
"""

from .base import DistributionSampler, DisabledSampler, ImmediateSampler

from .continuous import (
    ExponentialSampler,
    ErlangSampler,
    GammaSampler,
    DeterministicSampler,
    UniformSampler,
    ParetoSampler,
    WeibullSampler,
    LognormalSampler,
    NormalSampler,
    BetaSampler,
    ReplayTraceSampler,
    MixtureSampler,
)

from .phase_type import (
    PHSampler,
    APHSampler,
    CoxianSampler,
    Coxian2Sampler,
    HyperExponentialSampler,
    HyperExp2Sampler,
    create_ph_from_moments,
)

from .map_bmap import (
    MAPSampler,
    MAP2Sampler,
    BMAPSampler,
    MMPPSampler,
    create_map_from_ph,
)

from .discrete import (
    PoissonSampler,
    BinomialSampler,
    GeometricSampler,
    NegativeBinomialSampler,
    BernoulliSampler,
    DiscreteUniformSampler,
    ZipfSampler,
    CategoricalSampler,
    EmpiricalDiscreteSampler,
)

from .trace import (
    FileTraceSampler,
    MultiColumnTraceSampler,
    OnlineTraceSampler,
    load_trace,
)

from .fitting import (
    FitResult,
    DistributionFitter,
    fit_distribution,
    moment_match,
    compute_scv,
    compute_moments,
)

__all__ = [
    # Base
    'DistributionSampler',
    'DisabledSampler',
    'ImmediateSampler',
    # Continuous
    'ExponentialSampler',
    'ErlangSampler',
    'GammaSampler',
    'DeterministicSampler',
    'UniformSampler',
    'ParetoSampler',
    'WeibullSampler',
    'LognormalSampler',
    'NormalSampler',
    'BetaSampler',
    'ReplayTraceSampler',
    'MixtureSampler',
    # Phase-type
    'PHSampler',
    'APHSampler',
    'CoxianSampler',
    'Coxian2Sampler',
    'HyperExponentialSampler',
    'HyperExp2Sampler',
    'create_ph_from_moments',
    # MAP/BMAP
    'MAPSampler',
    'MAP2Sampler',
    'BMAPSampler',
    'MMPPSampler',
    'create_map_from_ph',
    # Discrete
    'PoissonSampler',
    'BinomialSampler',
    'GeometricSampler',
    'NegativeBinomialSampler',
    'BernoulliSampler',
    'DiscreteUniformSampler',
    'ZipfSampler',
    'CategoricalSampler',
    'EmpiricalDiscreteSampler',
    # Trace replay
    'FileTraceSampler',
    'MultiColumnTraceSampler',
    'OnlineTraceSampler',
    'load_trace',
    # Fitting
    'FitResult',
    'DistributionFitter',
    'fit_distribution',
    'moment_match',
    'compute_scv',
    'compute_moments',
]
