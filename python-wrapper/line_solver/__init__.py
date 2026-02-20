
"""
LINE Solver for Python - Queueing Network Analysis

LINE (Library for INteractive Evaluation) is a library for analyzing queueing
networks via analytical methods and simulation. This Python package provides
a wrapper around the core LINE JAR implementation.

Key Features:
- Analytical solvers (MVA, Fluid, NC, CTMC, SSA)
- Simulation solvers (JMT integration)
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

For more information, see https://line-solver.org
"""

import jpype
import pandas as pd

from urllib.request import urlretrieve
import jpype.imports
from jpype import startJVM, shutdownJVM, java
import numpy as np
import os, sys

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


def jlineStart():
    """
    Initialize the Java LINE solver backend.

    Downloads the Java JAR file if not present, starts the JVM,
    and imports necessary Java classes for LINE solver functionality.
    This function is called automatically when importing line_solver.
    """
    with GlobalImport() as gi:
        package_dir = os.path.dirname(os.path.realpath(__file__))

        # Look for JAR in multiple locations (in order of preference):
        # 1. Inside the package (pip install)
        # 2. In common/ directory relative to repo root (dev environment)
        # 3. Download from GitHub releases
        jar_locations = [
            os.path.join(package_dir, "jline.jar"),  # Inside package (pip install)
            os.path.join(os.path.dirname(os.path.dirname(package_dir)), 'common', "jline.jar"),  # Dev: ../../common/
        ]

        jar_file_path = None
        for loc in jar_locations:
            if os.path.isfile(loc):
                jar_file_path = loc
                break

        # If not found, download to package directory
        if jar_file_path is None:
            jar_file_path = os.path.join(package_dir, "jline.jar")
            print("Downloading LINE solver JAR, please wait... ", end='')
            urlretrieve("https://sourceforge.net/projects/line-solver/files/jline.jar/download",
                        jar_file_path)
            print("done.")

        # Look for JMT.jar (contains dependencies like EJML)
        jmt_jar_path = None
        jmt_locations = [
            os.path.join(package_dir, "JMT.jar"),  # Inside package (pip install)
            os.path.join(os.path.dirname(os.path.dirname(package_dir)), 'common', "JMT.jar"),  # Dev: ../../common/
        ]
        for loc in jmt_locations:
            if os.path.isfile(loc):
                jmt_jar_path = loc
                break

        # Build classpath with both JARs if JMT.jar exists
        classpath = [jar_file_path]
        if jmt_jar_path:
            classpath.append(jmt_jar_path)

        # Only start JVM if not already running
        if not jpype.isJVMStarted():
            # Pass classpath as argument to startJVM for reliable loading
            jpype.startJVM(classpath=classpath)
            #jpype.startJVM(jpype.getDefaultJVMPath(), classpath=[jar_file_path],
            #               "-agentlib:jdwp=transport=dt_socket,server=n,suspend=y,address=5005")
        else:
            # JVM already running, add JAR to classpath dynamically
            jpype.addClassPath(jar_file_path)
        # Use JPackage to access Java classes (avoids conflict with line_solver/jline.py)
        jline_pkg = jpype.JPackage('jline')
        GlobalConstants = jline_pkg.GlobalConstants
        Chain = jline_pkg.lang.Chain
        Element = jline_pkg.lang.Element
        Ensemble = jline_pkg.lang.Ensemble
        Metric = jline_pkg.lang.Metric
        FeatureSet = jline_pkg.lang.FeatureSet
        Region = jline_pkg.lang.Region
        Model = jline_pkg.lang.Model
        NetworkAttribute = jline_pkg.lang.NetworkAttribute
        NetworkElement = jline_pkg.lang.NetworkElement
        Event = jline_pkg.lang.Event
        ItemSet = jline_pkg.lang.ItemSet
        NodeAttribute = jline_pkg.lang.NodeAttribute
        OutputStrategy = jline_pkg.lang.OutputStrategy
        ServiceBinding = jline_pkg.lang.ServiceBinding
        # Store Java ActivityPrecedence for wrapper initialization (not exported directly)
        _JavaActivityPrecedence = jline_pkg.lang.layered.ActivityPrecedence
        # ActivityPrecedence will be exported as the wrapper
        CacheTask = jline_pkg.lang.layered.CacheTask
        FunctionTask = jline_pkg.lang.layered.FunctionTask
        LayeredNetworkElement = jline_pkg.lang.layered.LayeredNetworkElement
        LayeredNetworkStruct = jline_pkg.lang.layered.LayeredNetworkStruct
        ItemEntry = jline_pkg.lang.layered.ItemEntry
        Host = jline_pkg.lang.layered.Host
        ContinuousDistribution = jline_pkg.lang.processes.ContinuousDistribution
        Coxian = jline_pkg.lang.processes.Coxian
        DiscreteDistribution = jline_pkg.lang.processes.DiscreteDistribution
        DiscreteSampler = jline_pkg.lang.processes.DiscreteSampler
        Distribution = jline_pkg.lang.processes.Distribution
        Markovian = jline_pkg.lang.processes.Markovian
        Logger = jline_pkg.lang.nodes.Logger
        Place = jline_pkg.lang.nodes.Place
        StatefulNode = jline_pkg.lang.nodes.StatefulNode
        Station = jline_pkg.lang.nodes.Station
        Transition = jline_pkg.lang.nodes.Transition
        MarkedMAP = jline_pkg.lang.processes.MarkedMAP
        MarkedMMPP = jline_pkg.lang.processes.MarkedMMPP
        Buffer = jline_pkg.lang.sections.Buffer
        CacheClassSwitcher = jline_pkg.lang.sections.CacheClassSwitcher
        ClassSwitcher = jline_pkg.lang.sections.ClassSwitcher
        Dispatcher = jline_pkg.lang.sections.Dispatcher
        Forker = jline_pkg.lang.sections.Forker
        InfiniteServer = jline_pkg.lang.sections.InfiniteServer
        InputSection = jline_pkg.lang.sections.InputSection
        Joiner = jline_pkg.lang.sections.Joiner
        OutputSection = jline_pkg.lang.sections.OutputSection
        PreemptiveServer = jline_pkg.lang.sections.PreemptiveServer
        RandomSource = jline_pkg.lang.sections.RandomSource
        Section = jline_pkg.lang.sections.Section
        Server = jline_pkg.lang.sections.Server
        ServiceSection = jline_pkg.lang.sections.ServiceSection
        ServiceTunnel = jline_pkg.lang.sections.ServiceTunnel
        SharedServer = jline_pkg.lang.sections.SharedServer
        StatefulClassSwitcher = jline_pkg.lang.sections.StatefulClassSwitcher
        StatelessClassSwitcher = jline_pkg.lang.sections.StatelessClassSwitcher
        State = jline_pkg.lang.state.State
        EnsembleSolver = jline_pkg.solvers.EnsembleSolver
        NetworkAvgTable = jline_pkg.solvers.NetworkAvgTable
        NetworkSolver = jline_pkg.solvers.NetworkSolver
        SolverAvgHandles = jline_pkg.solvers.SolverAvgHandles
        SolverTranHandles = jline_pkg.solvers.SolverTranHandles
        gi()
        jline_pkg.util.Maths.setMatlabRandomSeed(True)


def jlineMapMatrixToArray(mapmatrix):
    d = dict(mapmatrix)
    for i in range(len(d)):
        d[i] = jlineMatrixToArray(d[i])
    return d

def jlineMatrixCellToArray(matrixcell):
    d = {}
    for i in range(matrixcell.size()):
        matrix = matrixcell.get(i)
        if matrix is not None:
            d[i] = jlineMatrixToArray(matrix)
        else:
            d[i] = None
    return d


def jlineFromDistribution(distrib):
    python_distrib = None
    if distrib is not None:
        # Handle both Python wrappers (get_name) and raw Java objects (getName)
        if hasattr(distrib, 'get_name'):
            distrib_name = distrib.get_name()
        else:
            distrib_name = distrib.getName()
        if distrib_name == 'APH':
            python_distrib = APH(distrib)
        elif distrib_name == 'Cox2':
            python_distrib = Cox2(distrib)
        elif distrib_name == 'Det':
            python_distrib = Det(distrib)
        elif distrib_name == 'Disabled':
            python_distrib = Disabled()
        elif distrib_name == 'Erlang':
            python_distrib = Erlang(distrib)
        elif distrib_name == 'Exp':
            python_distrib = Exp(distrib)
        elif distrib_name == 'Gamma':
            python_distrib = Gamma(distrib)
        elif distrib_name == 'HyperExp':
            python_distrib = HyperExp(distrib)
        elif distrib_name == 'Immediate':
            python_distrib = Immediate()
        elif distrib_name == 'Lognormal':
            python_distrib = Lognormal(distrib)
        elif distrib_name == 'MAP':
            python_distrib = MAP(distrib)
        elif distrib_name == 'Pareto':
            python_distrib = Pareto(distrib)
        elif distrib_name == 'PH':
            python_distrib = PH(distrib)
        elif distrib_name == 'Replayer':
            python_distrib = Replayer(distrib)
        elif distrib_name == 'Uniform':
            python_distrib = Uniform(distrib)
        elif distrib_name == 'Weibull':
            python_distrib = Weibull(distrib)
        elif distrib_name == 'Binomial':
            python_distrib = Binomial(distrib)
        elif distrib_name == 'DiscreteSampler':
            python_distrib = DiscreteSampler(distrib)
        elif distrib_name == 'Geometric':
            python_distrib = Geometric(distrib)
        elif distrib_name == 'Poisson':
            python_distrib = Poisson(distrib)
        elif distrib_name == 'Zipf':
            python_distrib = Zipf(distrib)
    return python_distrib


def jlineMatrixToArray(matrix):
    """
    Convert Java LINE Matrix to numpy array.
    
    Args:
        matrix: Java Matrix object from LINE solver.
        
    Returns:
        ndarray: NumPy array representation, or None if input is None.
    """
    if matrix is None:
        return None
    else:
        return np.array(list(matrix.toArray2D()))


def jlineMatrixFromArray(array):
    """
    Convert numpy array or list to Java LINE Matrix.

    Args:
        array (array_like): NumPy array or Python list to convert.

    Returns:
        Java Matrix object compatible with LINE solver.
    """
    # Handle scalar inputs (int, float, numpy scalar) by wrapping as 1x1 2D array
    if np.isscalar(array):
        array = np.array([[array]], dtype=float)
    elif isinstance(array, list):
        array = np.array(array)

    # Use Matrix(int, int) constructor and set values individually
    if len(np.shape(array)) > 1:
        rows, cols = array.shape[0], array.shape[1]
        ret = jpype.JPackage('jline').util.matrix.Matrix(rows, cols)
        for i in range(rows):
            for j in range(cols):
                # Convert to Python float first to avoid numpy type issues
                ret.set(i, j, float(array[i, j]))
    else:
        # For 1D array, create a row vector
        length = array.shape[0]
        ret = jpype.JPackage('jline').util.matrix.Matrix(1, length)
        for i in range(length):
            ret.set(0, i, float(array[i]))
    return ret


def jlineMatrixZeros(rows, cols):
    return jlineMatrixFromArray([[0.0] * cols for _ in range(rows)])


def jlineMatrixSingleton(value):
    return jlineMatrixFromArray([value])


def is_interactive():
    import __main__ as main
    return not hasattr(main, '__file__')


def lineDefaults(solverName='Solver'):
    from .solvers import SolverOptions
    from .constants import SolverType

    if solverName == 'Solver':
        return SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.MVA)
    else:
        solver_type_map = {
            'AG': jpype.JPackage('jline').lang.constant.SolverType.AG,
            'MVA': jpype.JPackage('jline').lang.constant.SolverType.MVA,
            'CTMC': jpype.JPackage('jline').lang.constant.SolverType.CTMC,
            'JMT': jpype.JPackage('jline').lang.constant.SolverType.JMT,
            'SSA': jpype.JPackage('jline').lang.constant.SolverType.SSA,
            'MAM': jpype.JPackage('jline').lang.constant.SolverType.MAM,
            'FLUID': jpype.JPackage('jline').lang.constant.SolverType.FLUID,
            'NC': jpype.JPackage('jline').lang.constant.SolverType.NC,
            'AUTO': jpype.JPackage('jline').lang.constant.SolverType.AUTO,
            'ENV': jpype.JPackage('jline').lang.constant.SolverType.ENV,
            'LQNS': jpype.JPackage('jline').lang.constant.SolverType.LQNS,
            'LN': jpype.JPackage('jline').lang.constant.SolverType.LN,
            'QNS': jpype.JPackage('jline').lang.constant.SolverType.QNS
        }

        solver_type = solver_type_map.get(solverName.upper(), jpype.JPackage('jline').lang.constant.SolverType.MVA)
        return SolverOptions(solver_type)


# Helper to add snake_case support to Java ActivityPrecedence
_ActivityPrecedenceJavaClass = None

def _add_activity_precedence_snake_case(java_class):
    """Add snake_case aliases to Java ActivityPrecedence class via a wrapper module."""
    class ActivityPrecedenceProxy:
        """Proxy that provides both PascalCase and snake_case methods."""
        def __getattr__(self, name):
            return getattr(java_class, name)

        @staticmethod
        def Serial(act0, act1):
            return java_class.Serial(act0, act1)

        @staticmethod
        def serial(acts):
            """Python snake_case alias for Serial(). Takes a list of activities.

            Returns an array of ActivityPrecedence objects representing the serial chain.
            For N activities, returns N-1 precedence relationships.
            """
            if not acts:
                raise ValueError("serial() requires at least one activity")
            if len(acts) == 1:
                return java_class.Serial(acts[0], acts[0])
            # Build a Java list and call Serial(List) which returns ActivityPrecedence[]
            java_list = jpype.java.util.ArrayList()
            for act in acts:
                java_list.add(act)
            return java_class.Serial(java_list)

        @staticmethod
        def AndFork(preAct, postActs):
            return java_class.AndFork(preAct, postActs)

        @staticmethod
        def and_fork(preAct, postActs):
            return ActivityPrecedenceProxy.AndFork(preAct, postActs)

        @staticmethod
        def AndJoin(preActs, postAct):
            return java_class.AndJoin(preActs, postAct)

        @staticmethod
        def and_join(preActs, postAct):
            return ActivityPrecedenceProxy.AndJoin(preActs, postAct)

        @staticmethod
        def OrFork(preAct, postActs, probs):
            return java_class.OrFork(preAct, postActs, probs)

        @staticmethod
        def or_fork(preAct, postActs, probs):
            return ActivityPrecedenceProxy.OrFork(preAct, postActs, probs)

        @staticmethod
        def OrJoin(preActs, postAct):
            return java_class.OrJoin(preActs, postAct)

        @staticmethod
        def or_join(preActs, postAct):
            return ActivityPrecedenceProxy.OrJoin(preActs, postAct)

        @staticmethod
        def Loop(preAct, postActs, counts):
            return java_class.Loop(preAct, postActs, counts)

        @staticmethod
        def loop(preAct, postActs, counts):
            return ActivityPrecedenceProxy.Loop(preAct, postActs, counts)

        @staticmethod
        def CacheAccess(*args):
            return java_class.CacheAccess(*args)

        @staticmethod
        def cache_access(*args):
            return ActivityPrecedenceProxy.CacheAccess(*args)

        # Forward all other attributes to java_class
        Xor = staticmethod(lambda *args, **kw: java_class.Xor(*args, **kw))
        xor = staticmethod(lambda *args, **kw: java_class.Xor(*args, **kw))
        fromActivities = staticmethod(lambda *args, **kw: java_class.fromActivities(*args, **kw))
        getPrecedenceId = staticmethod(lambda *args, **kw: java_class.getPrecedenceId(*args, **kw))

    return ActivityPrecedenceProxy


# Start JVM on import
jlineStart()
# Create ActivityPrecedence proxy with snake_case support
_ActivityPrecedenceProxy = _add_activity_precedence_snake_case(jpype.JPackage('jline').lang.layered.ActivityPrecedence)

# Import from standard modules
from .api import *
from .constants import *
from .utils import *
from .solvers import *

# Import from lang and distributions (wrapper classes)
from .lang import *
from .distributions import *

# Import from layered
from .layered import *

from .io import QN2JSIMG, qn2jsimg, JLINE, LQN2QN, lqn2qn

# Import gallery after lang and distributions to avoid circular import
from .gallery import *

