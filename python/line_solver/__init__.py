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
pd.set_option('display.precision', 4)


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
from .environment import Environment, SolverENV, ENV

# Import I/O functions (native only)
from .api.io import qn2jsimg, lqn2qn

# Import gallery after lang and distributions to avoid circular import
from .gallery import *

# Native implementations
from . import distributions
