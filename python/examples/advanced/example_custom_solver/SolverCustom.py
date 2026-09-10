"""
SolverCustom  Template for a user-written LINE solver.

The MATLAB twin is a `@SolverCustom` class folder holding SolverCustom.m and
runAnalyzer.m; Python has no class-folder idiom, so both live here, and the two
free functions the analyzer calls stay in solver_custom.py and
solver_custom_analyzer.py exactly as they do there.

To write a real solver: fill in solver_custom(), widen getFeatureSet() to the
features the algorithm actually supports, and list the method names it accepts
in listValidMethods().
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from line_solver import lineDefaults                                   # noqa: E402
from line_solver.solvers.base import NetworkSolver                     # noqa: E402
from solver_custom_analyzer import solver_custom_analyzer              # noqa: E402


class SolverCustom(NetworkSolver):
    """Example of custom solver."""

    def __init__(self, model, *args, **kwargs):
        # NetworkSolver declares no __init__ here: each concrete solver sets
        # model / sn / options itself, which is what MATLAB's
        # self@NetworkSolver(model, mfilename) does for the class-folder twin.
        self.model = model
        self.network = model
        self.name = 'SolverCustom'
        opts = self.defaultOptions()
        if args and isinstance(args[0], str):
            opts.method = args[0]
        for key, value in kwargs.items():
            opts[key] = value
        self.options = opts
        self._sn = None
        self._result = None

    def getStruct(self):
        """The data structure summarizing the model (no initial state needed)."""
        return self.model.getStruct()

    def listValidMethods(self):
        """The method names this solver accepts."""
        return ['default']

    @staticmethod
    def getFeatureSet():
        """The model features this solver claims to support.

        Kept identical to the MATLAB template's list, which is deliberately
        broad: narrow it to what the algorithm can really answer, or `supports`
        will accept a model that solver_custom then gets wrong.
        """
        return {
            'Sink', 'Source', 'Router', 'ClassSwitch', 'DelayStation', 'Queue',
            'Fork', 'Join', 'Forker', 'Joiner', 'Logger',
            'Coxian', 'Cox2', 'APH', 'Erlang', 'Exp', 'HyperExp', 'Det',
            'Gamma', 'Lognormal', 'MAP', 'MMPP2', 'Normal', 'PH', 'Pareto',
            'Weibull', 'Replayer', 'Uniform',
            'StatelessClassSwitcher', 'InfiniteServer', 'SharedServer',
            'Buffer', 'Dispatcher', 'Server', 'JobSink', 'RandomSource',
            'ServiceTunnel', 'LogTunnel', 'Linkage',
            'Enabling', 'Timing', 'Firing', 'Storage', 'Place', 'Transition',
            'SchedStrategy_INF', 'SchedStrategy_PS', 'SchedStrategy_DPS',
            'SchedStrategy_FCFS', 'SchedStrategy_GPS', 'SchedStrategy_SIRO',
            'SchedStrategy_HOL', 'SchedStrategy_LCFS', 'SchedStrategy_LCFSPR',
            'SchedStrategy_SEPT', 'SchedStrategy_LEPT', 'SchedStrategy_SJF',
            'SchedStrategy_LJF',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND',
            'RoutingStrategy_RROBIN', 'RoutingStrategy_WRROBIN',
            'RoutingStrategy_SQ', 'SchedStrategy_EXT',
            'ClosedClass', 'OpenClass',
        }

    @staticmethod
    def defaultOptions():
        options = lineDefaults()
        options['timespan'] = [float('inf'), float('inf')]
        return options

    def runAnalyzer(self):
        """Run the solver: check the model, call the analyzer, publish results."""
        t0 = time.time()
        sn = self.getStruct()          # doesn't need initial state

        QN, UN, RN, TN, CN, XN, runtime = solver_custom_analyzer(sn, self.options)

        self._result = {
            'QN': QN, 'UN': UN, 'RN': RN, 'TN': TN,
            'CN': CN, 'XN': XN, 'AN': TN, 'WN': RN,
            'runtime': runtime,
            'method': str(getattr(self.options, 'method', 'default')),
        }
        self.runtime = time.time() - t0
        return self._result


if __name__ == '__main__':
    from line_solver import (ClosedClass, Delay, Exp, GlobalConstants, Network,
                             Queue, SchedStrategy, VerboseLevel)

    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('customExample')
    delay = Delay(model, 'Think')
    queue = Queue(model, 'Queue1', SchedStrategy.PS)
    jobclass = ClosedClass(model, 'Class1', 2, delay)
    delay.setService(jobclass, Exp(1.0))
    queue.setService(jobclass, Exp(2.0))
    model.link(Network.serialRouting(delay, queue))

    print(SolverCustom(model).runAnalyzer())
