"""
Native Python implementations of LINE solvers.

These implementations use pure Python/NumPy algorithms instead of
the Java/JPype wrapper, providing faster execution and no .
"""

# Import base Solver and SolverOptions from parent module
import sys
import os
_parent_module = os.path.join(os.path.dirname(__file__), '..', 'solvers.py')
import importlib.util
_spec = importlib.util.spec_from_file_location("_solvers_base", _parent_module)
_solvers_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_solvers_base)
Solver = _solvers_base.Solver
SolverOptions = _solvers_base.SolverOptions

from .solver_mva import SolverMVA
from .solver_ln import SolverLN
from .solver_auto import SolverAuto, SolverAutoOptions, ModelAnalyzer
from .solver_qns import SolverQNS, QNSOptions, QNSResult
from .solver_lqns import SolverLQNS, LQNSOptions, LQNSResult
from .solver_des.solver_des import SolverDES
from .solver_des.des_options import DESOptions, DESResult
from .solver_mam import SolverMAM, SolverMAMOptions
from .solver_fld import SolverFLD
from .solver_fld.options import SolverFLDOptions, FLDResult
from .solver_ctmc import SolverCTMC, SolverCTMCOptions
from .solver_ssa import SolverSSA, SolverSSAOptions
from .solver_nc import SolverNC, SolverNCOptions
from .solver_jmt import SolverJMT, SolverJMTOptions
from .solver_posterior import SolverPosterior, PosteriorOptions, PosteriorResult
from .convert import wrapper_sn_to_native, native_model_to_struct, get_native_sched_strategy

# Short aliases for solver classes (MATLAB-style)
MVA = SolverMVA
NC = SolverNC
CTMC = SolverCTMC
SSA = SolverSSA
FLD = SolverFLD
MAM = SolverMAM
JMT = SolverJMT
DES = SolverDES
AUTO = SolverAuto
LINE = SolverAuto
LN = SolverLN
QNS = SolverQNS
LQNS = SolverLQNS

__all__ = [
    'Solver',
    'SolverOptions',
    'SolverMVA',
    'SolverLN',
    'SolverAuto',
    'SolverAutoOptions',
    'ModelAnalyzer',
    'SolverQNS',
    'QNSOptions',
    'QNSResult',
    'SolverLQNS',
    'LQNSOptions',
    'LQNSResult',
    'SolverDES',
    'DESOptions',
    'DESResult',
    'SolverMAM',
    'SolverMAMOptions',
    'SolverFLD',
    'SolverFLDOptions',
    'FLDResult',
    'SolverCTMC',
    'SolverCTMCOptions',
    'SolverSSA',
    'SolverSSAOptions',
    'SolverNC',
    'SolverNCOptions',
    'SolverJMT',
    'SolverJMTOptions',
    'SolverPosterior',
    'PosteriorOptions',
    'PosteriorResult',
    'wrapper_sn_to_native',
    'native_model_to_struct',
    'get_native_sched_strategy',
    # Short aliases
    'MVA',
    'NC',
    'CTMC',
    'SSA',
    'FLD',
    'MAM',
    'JMT',
    'DES',
    'AUTO',
    'LINE',
    'LN',
    'QNS',
    'LQNS',
]
