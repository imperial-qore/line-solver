"""
Statistics collection and analysis for DES solver.

This package contains classes for collecting simulation statistics,
warmup detection (MSER-5), and confidence interval computation (OBM).
"""

from .metrics import (
    TimeWeightedAccumulator,
    ResponseTimeTally,
    ClassMetrics,
    StationMetrics,
    SimulationMetrics,
)

from .warmup import (
    MSER5TransientDetector,
    WindowedMSER5Detector,
    MultivariateMSER5,
    schruben_rule,
    welch_method,
)

from .confidence import (
    ConfidenceInterval,
    OBMConfidenceInterval,
    NonOverlappingBatchMeans,
    SpectralConfidenceInterval,
    compute_confidence_interval,
    relative_precision_stopping,
)

from .transient import (
    TransientSample,
    TransientCollector,
    SteadyStateDetector,
    compute_transient_statistics,
)

__all__ = [
    # Metrics
    'TimeWeightedAccumulator',
    'ResponseTimeTally',
    'ClassMetrics',
    'StationMetrics',
    'SimulationMetrics',
    # Warmup detection
    'MSER5TransientDetector',
    'WindowedMSER5Detector',
    'MultivariateMSER5',
    'schruben_rule',
    'welch_method',
    # Confidence intervals
    'ConfidenceInterval',
    'OBMConfidenceInterval',
    'NonOverlappingBatchMeans',
    'SpectralConfidenceInterval',
    'compute_confidence_interval',
    'relative_precision_stopping',
    # Transient analysis
    'TransientSample',
    'TransientCollector',
    'SteadyStateDetector',
    'compute_transient_statistics',
]
