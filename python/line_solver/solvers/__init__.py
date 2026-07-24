"""
Native Python implementations of LINE solvers.

These implementations use pure Python/NumPy algorithms.
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
from .solver_ba import SolverBA
from .solver_ln import SolverLN, SolverLNOptions
from .solver_auto import SolverAuto, SolverAutoOptions, ModelAnalyzer
from .wrappers.solver_qns import SolverQNS, QNSOptions, QNSResult
from .wrappers.solver_lqns import SolverLQNS, LQNSOptions, LQNSResult
try:
    from .wrappers.solver_ldes.solver_ldes import SolverLDES
except ImportError:
    SolverLDES = None
from .wrappers.solver_ldes.ldes_options import LDESOptions, LDESResult
from .solver_mam import SolverMAM, SolverMAMOptions
from .solver_fld import SolverFLD
from .solver_fld.options import SolverFLDOptions, FLDResult
from .solver_ctmc import SolverCTMC, SolverCTMCOptions
from .solver_ssa import SolverSSA, SolverSSAOptions, SamplePath, SampleEvent
from .solver_nc import SolverNC, SolverNCOptions
from .wrappers.solver_jmt import SolverJMT, SolverJMTOptions
from .solver_uq import SolverUQ, UQOptions, UQResult, EmpiricalCDF

# Short aliases for solver classes (MATLAB-style)
MVA = SolverMVA
BA = SolverBA
NC = SolverNC
CTMC = SolverCTMC
SSA = SolverSSA
FLD = SolverFLD
SolverFluid = SolverFLD  # MATLAB/wrapper class-name alias for portability
Fluid = SolverFLD
MAM = SolverMAM
JMT = SolverJMT
LDES = SolverLDES
AUTO = SolverAuto
SolverAUTO = SolverAuto  # MATLAB/JAR class-name alias for portability
LINE = SolverAuto
LN = SolverLN
QNS = SolverQNS
LQNS = SolverLQNS
UQ = SolverUQ

__all__ = [
    'Solver',
    'SolverOptions',
    'SolverMVA',
    'SolverBA',
    'BA',
    'SolverLN',
    'SolverLNOptions',
    'SolverAuto',
    'SolverAUTO',
    'SolverAutoOptions',
    'ModelAnalyzer',
    'SolverQNS',
    'QNSOptions',
    'QNSResult',
    'SolverLQNS',
    'LQNSOptions',
    'LQNSResult',
    'SolverLDES',
    'LDESOptions',
    'LDESResult',
    'SolverMAM',
    'SolverMAMOptions',
    'SolverFLD',
    'SolverFluid',
    'Fluid',
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
    'SolverUQ',
    'UQOptions',
    'UQResult',
    'EmpiricalCDF',
    # Short aliases
    'MVA',
    'NC',
    'CTMC',
    'SSA',
    'FLD',
    'MAM',
    'JMT',
    'LDES',
    'AUTO',
    'LINE',
    'LN',
    'QNS',
    'LQNS',
    'UQ',
]
