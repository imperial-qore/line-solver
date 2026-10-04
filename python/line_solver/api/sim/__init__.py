"""
Native Python implementations for simulation output analysis.

This module turns a simulation sample path into an interval estimate. Its subject
is the statistics of the output process, not the queueing model that produced it,
so nothing here takes a NetworkStruct: the input is a sequence of observations
such as successive waiting times exported from a simulation run.

The steady-state mean is already covered elsewhere: the LDES engine forms batch
means and reports confidence-interval half-widths for every average metric, with
MSER-5 warmup detection. What this module adds is the steady-state quantile, which
is the quantity a tail service-level agreement is written against and for which no
interval was previously available anywhere in LINE.

Key algorithms:
    Quantile intervals: fquest (one sample path), firquest (replications)
    STS machinery: sts_quantile_areas
    Hypothesis tests: vonneumann (randomness), shapirowilk (normality)
    Distributions: normcdf, norminv, tinv
"""

from .dist import normcdf, norminv, tinv
from .fquest import firquest, firquest_batchcounts, fquest, quest_options
from .sts import DEFAULT_WEIGHT, sts_quantile_areas
from .tests import shapirowilk, shapirowilk_weights, vonneumann

from .runlength import sim_asymvar_mm1, sim_asymvar_ctmc, sim_runlength, sim_runlength_plan

__all__ = [
    'normcdf',
    'norminv',
    'tinv',
    'vonneumann',
    'shapirowilk',
    'shapirowilk_weights',
    'sts_quantile_areas',
    'DEFAULT_WEIGHT',
    'fquest',
    'firquest',
    'quest_options',
    'firquest_batchcounts',
    'sim_asymvar_mm1',
    'sim_asymvar_ctmc',
    'sim_runlength',
    'sim_runlength_plan',
]
