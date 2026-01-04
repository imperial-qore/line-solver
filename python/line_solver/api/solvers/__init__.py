"""
Native Python solver implementations.

Provides solver handlers for various queueing network analysis methods.

Available solvers:
    nc: Normalizing Constant solver for product-form networks
    mva: Mean Value Analysis solver for product-form networks
    fld: Fluid/Mean-Field Approximation solver
    mam: Matrix-Analytic Methods solver
    env: Environment/Ensemble solver for CTMC-modulated systems
    ctmc: Continuous-Time Markov Chain solver for exact analysis
    ssa: Stochastic Simulation Algorithm solver
"""

from .nc import (
    solver_nc,
    solver_ncld,
    SolverNCReturn,
    SolverNCLDReturn,
    SolverOptions,
    solver_nc_margaggr,
    solver_nc_jointaggr,
    to_marginal_aggr,
    StateMarginalStatistics,
    SolverNCMargReturn,
    SolverNCJointReturn,
    solver_nc_analyzer,
    solver_ncld_analyzer,
    NCResult,
    NCResultProb,
)

from .mva import (
    solver_mva,
    SolverMVAReturn,
    SolverMVAOptions,
    solver_mva_analyzer,
    solver_amva,
    solver_qna,
    MVAResult,
)

from .fld import (
    solver_fld,
    SolverFLDReturn,
    SolverFLDOptions,
    solver_fld_analyzer,
    FLDResult,
    matrix_method_analyzer,
)

from .mam import (
    solver_mam,
    solver_mam_basic,
    SolverMAMReturn,
    SolverMAMOptions,
    solver_mam_analyzer,
    MAMResult,
)

from .env import (
    solver_env,
    solver_env_basic,
    SolverENVReturn,
    SolverENVOptions,
    solver_env_analyzer,
    ENVResult,
)

from .ctmc import (
    solver_ctmc,
    solver_ctmc_basic,
    SolverCTMCReturn,
    SolverCTMCOptions,
    solver_ctmc_analyzer,
    CTMCResult,
)

from .ssa import (
    solver_ssa,
    solver_ssa_basic,
    SolverSSAReturn,
    SolverSSAOptions,
    solver_ssa_analyzer,
    SSAResult,
)

from .qns import (
    solver_qns,
    SolverQNSOptions,
    SolverQNSReturn,
    is_qns_available,
)

from .jmt import (
    solver_jmt,
    SolverJMTOptions,
    SolverJMTReturn,
    is_jmt_available,
)

from .lqns import (
    solver_lqns,
    solver_lqns_from_model,
    SolverLQNSOptions,
    SolverLQNSReturn,
    is_lqns_available,
)

__all__ = [
    # NC solver
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
    # MVA solver
    'solver_mva',
    'SolverMVAReturn',
    'SolverMVAOptions',
    'solver_mva_analyzer',
    'solver_amva',
    'solver_qna',
    'MVAResult',
    # FLD solver
    'solver_fld',
    'SolverFLDReturn',
    'SolverFLDOptions',
    'solver_fld_analyzer',
    'FLDResult',
    'matrix_method_analyzer',
    # MAM solver
    'solver_mam',
    'solver_mam_basic',
    'SolverMAMReturn',
    'SolverMAMOptions',
    'solver_mam_analyzer',
    'MAMResult',
    # ENV solver
    'solver_env',
    'solver_env_basic',
    'SolverENVReturn',
    'SolverENVOptions',
    'solver_env_analyzer',
    'ENVResult',
    # CTMC solver
    'solver_ctmc',
    'solver_ctmc_basic',
    'SolverCTMCReturn',
    'SolverCTMCOptions',
    'solver_ctmc_analyzer',
    'CTMCResult',
    # SSA solver
    'solver_ssa',
    'solver_ssa_basic',
    'SolverSSAReturn',
    'SolverSSAOptions',
    'solver_ssa_analyzer',
    'SSAResult',
    # QNS solver (external tool)
    'solver_qns',
    'SolverQNSOptions',
    'SolverQNSReturn',
    'is_qns_available',
    # JMT solver (external tool)
    'solver_jmt',
    'SolverJMTOptions',
    'SolverJMTReturn',
    'is_jmt_available',
    # LQNS solver (external tool)
    'solver_lqns',
    'solver_lqns_from_model',
    'SolverLQNSOptions',
    'SolverLQNSReturn',
    'is_lqns_available',
]
