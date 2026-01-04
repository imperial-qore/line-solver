"""
NC Solver handlers for native Python implementation.

Provides the core analysis handlers for the Normalizing Constant (NC) solver,
including support for standard product-form networks, load-dependent service,
probability computation, and analyzers.

Key handlers:
    solver_nc: Main NC solver handler
    solver_ncld: Load-dependent NC handler
    solver_nc_margaggr: Aggregated marginal probability handler
    solver_nc_jointaggr: Aggregated joint probability handler

Analyzers:
    solver_nc_analyzer: Main analyzer with interpolation
    solver_ncld_analyzer: Load-dependent analyzer with interpolation
"""

from .handler import (
    solver_nc,
    solver_ncld,
    SolverNCReturn,
    SolverNCLDReturn,
    SolverOptions,
)

from .prob import (
    solver_nc_margaggr,
    solver_nc_jointaggr,
    to_marginal_aggr,
    StateMarginalStatistics,
    SolverNCMargReturn,
    SolverNCJointReturn,
)

from .analyzers import (
    solver_nc_analyzer,
    solver_ncld_analyzer,
    NCResult,
    NCResultProb,
)

__all__ = [
    'solver_nc',
    'solver_ncld',
    'SolverNCReturn',
    'SolverNCLDReturn',
    'SolverOptions',
    'solver_nc_margaggr',
    'solver_nc_jointaggr',
    'to_marginal_aggr',
    'StateMarginalStatistics',
    'SolverNCMargReturn',
    'SolverNCJointReturn',
    'solver_nc_analyzer',
    'solver_ncld_analyzer',
    'NCResult',
    'NCResultProb',
]
