# In __init__.py
import jpype

from urllib.request import urlretrieve
import jpype.imports
from jpype import startJVM, shutdownJVM, java
import numpy as np
import os, sys

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(1, dir_path)


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

    # is called before the end of this block


def jlineStart():
    with GlobalImport() as gi:
        package_dir = os.path.dirname(os.path.realpath(__file__))
        jar_file_path = os.path.join(package_dir, f"jline.jar")
        if not os.path.isfile(jar_file_path):
            print("Downloading LINE solver JAR, please wait... ", end='')
            urlretrieve("https://github.com/imperial-qore/line-solver/raw/main/python/line_solver/jline.jar",
                        jar_file_path)
            print("done.")
        jpype.startJVM()
        #jpype.startJVM("-Xint", "-Xdebug", "-Xnoagent","-agentlib:jdwp=transport=dt_socket,server=y,suspend=n,address=*:5005")
        jpype.addClassPath('jline.jar')
        from jline.lang.constant import GlobalConstants
        from jline.lang import Chain, Element, Ensemble, Metric
        from jline.lang import FeatureSet, FiniteCapacityRegion, InputBinding
        from jline.lang import Model, NetworkAttribute, NetworkElement, NetworkEvent
        from jline.lang import ItemSet, NodeAttribute, OutputStrategy, ServiceBinding
        from jline.lang.layered import ActivityPrecedence, CacheTask, LayeredNetworkElement
        from jline.lang.layered import LayeredNetworkStruct, ItemEntry, Host
        from jline.lang.distributions import ContinuousDistribution, Coxian, CumulativeDistribution
        from jline.lang.distributions import DiscreteDistribution, DiscreteSampler, Distribution
        from jline.lang.distributions import MarkovianDistribution
        from jline.lang.nodes import Logger, Place
        from jline.lang.nodes import StatefulNode, Station, Transition
        from jline.lang.processes import MMAP, Process
        from jline.lang.sections import Buffer, CacheClassSwitcher, ClassSwitcher, ClassSwitchOutputSection, Dispatcher
        from jline.lang.sections import Forker, InfiniteServer, InputSection, Joiner, OutputSection, PreemptiveServer
        from jline.lang.sections import RandomSource, Section, Server, ServiceSection, ServiceTunnel, SharedServer
        from jline.lang.sections import StatefulClassSwitcher, StatelessClassSwitcher
        from jline.lang.state import State
        from jline.solvers import EnsembleSolver, NetworkAvgTable, NetworkSolver, SolverHandles

        gi()


def jlineMapMatrixToArray(mapmatrix):
    d = dict(mapmatrix)
    for i in range(len(d)):
        d[i] = jlineMatrixToArray(d[i])
    return d


def jlineFromDistribution(distrib):
    python_distrib = None
    if distrib is not None:
        distrib_name = distrib.getName()
        match distrib_name:
            case 'APH':
                python_distrib = APH(distrib)
            case 'Cox2':
                python_distrib = Cox2(distrib)
            case 'Det':
                python_distrib = Det(distrib)
            case 'Disabled':
                python_distrib = Disabled()
            case 'Erlang':
                python_distrib = Erlang(distrib)
            case 'Exp':
                python_distrib = Exp(distrib)
            case 'Gamma':
                python_distrib = Gamma(distrib)
            case 'HyperExp':
                python_distrib = HyperExp(distrib)
            case 'Immediate':
                python_distrib = Immediate()
            case 'Lognormal':
                python_distrib = Lognormal(distrib)
            case 'MAP':
                python_distrib = MAP(distrib)
            case 'Pareto':
                python_distrib = Pareto(distrib)
            case 'PH':
                python_distrib = PH(distrib)
            case 'Replayer':
                python_distrib = Replayer(distrib)
            case 'Uniform':
                python_distrib = Uniform(distrib)
            case 'Weibull':
                python_distrib = Weibull(distrib)
            case 'Binomial':
                python_distrib = Binomial(distrib)
            case 'DiscreteSampler':
                python_distrib = DiscreteSampler(distrib)
            case 'Geometric':
                python_distrib = Geometric(distrib)
            case 'Poisson':
                python_distrib = Poisson(distrib)
            case 'Zipf':
                python_distrib = Zipf(distrib)
    return python_distrib

def jlineMatrixToArray(matrix):
    if matrix is None:
        return None
    else:
        return np.array(list(matrix.toArray2D()))


def jlineMatrixFromArray(array):
    if isinstance(array, list):
        array = np.array(array)
    if len(np.shape(array)) > 1:
        ret = jpype.JPackage('jline').util.Matrix(np.size(array, 0), np.size(array, 1), array.size)
        for i in range(np.size(array, 0)):
            for j in range(np.size(array, 1)):
                ret.set(i, j, array[i][j])
    else:
        ret = jpype.JPackage('jline').util.Matrix(1, np.size(array, 0), array.size)
        for i in range(np.size(array, 0)):
            ret.set(0, i, array[i])
    return ret


def is_interactive():
    import __main__ as main
    return not hasattr(main, '__file__')


jlineStart()
from .api import *
from .constants import *
from .lang import *
from .utils import *
from .solvers import *
from .distributions import *
from .layered import *
from .lib import *
from .gallery import *
