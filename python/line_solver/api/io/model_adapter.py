"""
Re-export of the model-transformation primitives.

Every primitive lives in a single ModelAdapter, in line_solver.io.model_adapter,
mirroring MATLAB's matlab/src/io/@ModelAdapter/ directory and the JAR's
jline.lang.ModelAdapter. This module exists only so that the historical
line_solver.api.io import path keeps working.

Do not add a transformation here. A second copy of mmt, ht, sort_forks and the
four path enumerators used to live in this file while every production caller
used the line_solver.io copy, so the two drifted: the stale copy's MMTResult
had lost the routing_matrix, service_src, immediate_slots, arrival_src,
base_model and fork_lambda_init fields that the fork-join driver relies on.
"""

from ...io.model_adapter import (
    ModelAdapter,
    DeaggregationInfo,
    TaggedModelResult,
    FESAggregationInfo,
    MMTResult,
    HTResult,
    tag_chain,
    aggregate_chains,
    remove_class,
    sort_forks,
    mmt,
    ht,
    paths,
    paths_cs,
    find_paths,
    find_paths_cs,
    fjtag,
    aggregate_fes,
)

__all__ = [
    'ModelAdapter',
    'DeaggregationInfo',
    'TaggedModelResult',
    'FESAggregationInfo',
    'MMTResult',
    'HTResult',
    'tag_chain',
    'aggregate_chains',
    'remove_class',
    'sort_forks',
    'mmt',
    'ht',
    'paths',
    'paths_cs',
    'find_paths',
    'find_paths_cs',
    'fjtag',
    'aggregate_fes',
]
