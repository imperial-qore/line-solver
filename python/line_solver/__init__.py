# In __init__.py
import jpype
import jpype.imports
from jpype import startJVM, shutdownJVM, java
import numpy as np


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
        jpype.startJVM()
        jpype.addClassPath('jline.jar')
        from jline.lang import Chain, Element, Ensemble
        from jline.lang import Env, FeatureSet, FiniteCapacityRegion, InputBinding
        from jline.lang import Model, NetworkAttribute, NetworkElement, NetworkEvent, NetworkStruct
        from jline.lang import ItemSet, JobClass, NodeAttribute, OutputStrategy, ServiceBinding
        from jline.lang.layerednetworks import ActivityPrecedence, CacheTask, LayeredNetworkElement
        from jline.lang.layerednetworks import LayeredNetworkStruct, ItemEntry, Host
        from jline.lang.constant import ActivityPrecedenceType, CallType, DropStrategy, EventType, GlobalConstants
        from jline.lang.constant import JobClassType, JoinStrategy, Metric, MetricType, NodeType, ProcessType
        from jline.lang.constant import ReplacementStrategy, RoutingStrategy, SchedStrategyType
        from jline.lang.constant import ServiceStrategy, SolverType, TimingStrategy, VerboseLevel
        from jline.lang.distributions import APH, Binomial, ContinuousDistribution, Coxian, CumulativeDistribution
        from jline.lang.distributions import DiscreteDistribution, DiscreteSampler, Distribution
        from jline.lang.distributions import Gamma, Geometric, LogNormal, MarkovianDistribution
        from jline.lang.distributions import Pareto, PH, Poisson, Uniform, Weibull
        from jline.lang.nodes import ClassSwitch, Fork, Join, Logger, Node, Place
        from jline.lang.nodes import StatefulNode, Station, Transition
        from jline.lang.processes import MAP, Process
        from jline.lang.sections import Buffer, CacheClassSwitcher, ClassSwitcher, ClassSwitchOutputSection, Dispatcher
        from jline.lang.sections import Forker, InfiniteServer, InputSection, Joiner, OutputSection, PreemptiveServer
        from jline.lang.sections import RandomSource, Section, Server, ServiceSection, ServiceTunnel, SharedServer
        from jline.lang.sections import StatefulClassSwitcher, StatelessClassSwitcher
        from jline.lang.state import State
        from jline.solvers import EnsembleSolver, NetworkAvgTable, NetworkSolver, Solver, SolverHandles, SolverMetrics

        gi()


def jlineToArray(matrix):
    return np.array(list(matrix.toArray2D()))

jlineStart()
from .constants import *
from .lang import *
from .utils import *
from .solvers import *
from .distributions import *
from .layerednetworks import *
