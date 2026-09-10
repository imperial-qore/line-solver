"""
LINE Solver for Python - Queueing Network Analysis

LINE (Library for INteractive Evaluation) is a library for analyzing queueing
networks via analytical methods and simulation. This Python package provides
native Python implementations of the LINE solver algorithms.

Key Features:
- Analytical solvers (MVA, Fluid, NC, CTMC, SSA)
- Support for open, closed, and mixed networks
- Layered queueing networks (LQN) for software models
- Rich set of probability distributions
- Performance metrics and statistical analysis

Basic Usage:
    >>> from line_solver import *
    >>> model = Network('MyModel')
    >>> source = Source(model, 'Source')
    >>> queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    >>> sink = Sink(model, 'Sink')
    >>>
    >>> jobclass = OpenClass(model, 'Class1')
    >>> source.setArrival(jobclass, Exp(1.0))
    >>> queue.setService(jobclass, Exp(2.0))
    >>>
    >>> model.link(Network.serial_routing([source, queue, sink]))
    >>> solver = SolverMVA(model)
    >>> results = solver.avg_table()
    >>> print(results)

For more information, see https://line-solver.sf.net
"""

import pandas as pd
import numpy as np
import os
import sys

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(1, dir_path)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.precision', 5)


class GlobalImport:
    def __enter__(self):
        return self

    def __call__(self):
        import inspect
        self.collector = inspect.getargvalues(inspect.getouterframes(inspect.currentframe())[1].frame).locals

    def __exit__(self, *args):
        try:
            globals().update(self.collector)
        except:
            pass


def lineRootFolder():
    """
    Get the root folder path of the LINE solver installation.

    Returns:
        str: Absolute path to the LINE solver root directory.
    """
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def lineDefaults():
    """
    Get default solver options.

    MATLAB-style function to return default solver configuration options.

    Returns:
        OptionsDict: Dictionary-like object with default solver options.
    """
    from .solvers import Solver
    return Solver.default_options()


def is_interactive():
    import __main__ as main
    return not hasattr(main, '__file__')


def native_to_array(data):
    """
    Convert input data to a numpy array.

    Args:
        data: Input data (numpy array, list, or matrix-like object)

    Returns:
        numpy array representation
    """
    return np.asarray(data)


def jlineMatrixFromArray(data):
    """
    Convert input data to a numpy array (Native compatibility alias).
    
    In the wrapper version, this converts to a Java Matrix.
    In the native version, this ensures we have a numpy array.
    """
    return np.asarray(data)


def jlineMatrixToArray(data):
    """
    Convert input data to a numpy array (Native compatibility alias).
    
    In the wrapper version, this converts from a Java Matrix.
    In the native version, this ensures we have a numpy array.
    """
    return np.asarray(data)


# Import from standard modules
from .api import *
from .constants import *
from .utils import *
from .solvers import *

# Import native Python implementations
from .lang import *
from .distributions import *

# Import from layered
from .layered import *

# Import workflow classes
from .lang import Workflow

# Import reward classes
from .lang import Reward, RewardState, RewardStateView

# Import environment (random environment models)
from .environment import Environment, SolverENV, SolverEnv, ENV

# Import I/O functions (native only)
from .api.io import qn2jsimg, lqn2qn
# The exception every deliberate LINE refusal raises, so a caller can
# `except LineError` instead of catching every RuntimeError in sight. It is
# re-homed onto the package the way numpy re-homes its own public exceptions:
# a traceback then ends in `line_solver.LineError`, not in the private
# `line_solver.api.io.logging.LineError` the reader has no business knowing.
# Unpickling still resolves, because the name is bound right here.
from .api.io import LineError
LineError.__module__ = __name__
from .io.linemodel_io import save_model, load_model
from .io.pnml_io import save_pnml, load_pnml

# Environment check (mirror of MATLAB lineInstall)
from .install import line_install

# Import gallery after lang and distributions to avoid circular import
from .gallery import *

# Server-farm builder (mirror of jline.gen.Cluster)
from .gen import Cluster

# Native implementations
from . import distributions

# see _kb/05-solvers-overview.md (LineOpt) / _kb/11-conventions-and-gotchas.md (long-tail gotchas) for rationale
from .opt import (
    OptimizationProblem,
    DecisionVariable,
    ServerAllocation,
    StationReplicas,
    RoutingProbabilities,
    ClassServiceMapping,
    ServiceRate,
    JobPopulation,
    ClassPriority,
    HostDemand,
    ActivityThinkTime,
    TaskThinkTime,
    TaskMultiplicity,
    TaskReplication,
    ProcessorMultiplicity,
    Objective,
    MinimizeCost,
    MaximizePerformance,
    MinimizeSystemResponseTime,
    Constraint,
    ResponseTimeConstraint,
    SystemResponseTimeConstraint,
    ThroughputConstraint,
    UtilizationConstraint,
    BudgetConstraint,
    LineOptSolver,
    LineOptSolverOptions,
    BisectionSolver,
    ParetoPoint,
    ParetoSweep,
    SubProblem,
    DecompositionWorkflow,
    EvaluationResult,
    OptimizationResult,
    SubProblemResult,
    WorkflowResult,
    LineEvaluator,
)

# see _kb/11-conventions-and-gotchas.md (Python long-tail low-hit gotchas) for rationale
import types as _types
__all__ = [
    _name for _name in dir()
    if not _name.startswith('_')
    and not isinstance(globals().get(_name), _types.ModuleType)
]
del _types


def lineStart(verbose=True):
    """Print the LINE startup banner, as MATLAB's lineStart does.

    Importing line_solver stays silent, so this is the explicit entry point for
    a session banner. It also names where to obtain the algorithm attribution:
    LINE follows a pull-based model (as in Sage's sage.misc.citation), so
    nothing is printed during a solve.

    Args:
        verbose: print the banner when True

    Returns:
        the LINE version string
    """
    from .constants import GlobalConstants
    if verbose:
        print('Starting LINE version %s: StdOut=console, VerboseLevel=STD, '
              'DoChecks=true, CoarseTol=%.1e, FineTol=%.1e, Zero=%.1e, MaxInt=%d'
              % (GlobalConstants.Version, GlobalConstants.CoarseTol,
                 GlobalConstants.FineTol, GlobalConstants.Zero, GlobalConstants.MaxInt))
        print('Type model.help() for the solvers that support a model, '
              'solver.libraries() for third-party dependencies, '
              'solver.citations() for references.')
    return GlobalConstants.Version


line_start = lineStart


# RESULT RECORDING, off unless asked for. `LINE_RECORD_RESULTS=1` makes every
# result-table getter append what it returns to a process-wide buffer, together
# with the solver and method that produced it, so a consumer reads the values a
# script COMPUTED rather than the text it printed. That is how the parity suites
# assert against the shared goldens (see line_solver/result_recorder.py and
# python/tests/parity/). Imported last, because installing the hooks touches
# every solver module and they must already be defined; and imported at all only
# under the flag, so an ordinary run pays nothing.
if os.environ.get('LINE_RECORD_RESULTS'):
    from . import result_recorder as _result_recorder  # noqa: F401,E402
