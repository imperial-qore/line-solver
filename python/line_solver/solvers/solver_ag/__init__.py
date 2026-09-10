"""Agent-based (RCAT) solver.

Solves a network by the Reversed Compound Agent Theorem: every (station, class)
pair becomes an isolated agent, and the agents are coupled ONLY through the
reversed rates of the synchronizing actions. Agent k carries

    Q_k(x) = L_k + sum_{c passive at k} x_c Pb_c

and publishes, for every action it is active on, a scalar read off its own
stationary vector. The fixed point over that scalar vector is the whole
analysis, which is why it parallelises exactly rather than approximately -- see
:mod:`line_solver.solvers.solver_ag.exec_backend`.

Methods: ``inap``, ``inapplus``, ``inapinf``, and the vestigial ``exact`` alias
which warns and falls back to ``inap``. Per Marin, Rota Bulo and Balsamo, "A
Numerical Algorithm for the Decomposition of Cooperating Structured Markov
Processes", MASCOTS 2012.

These methods used to live in SolverMAM. They are the only algorithms in LINE
that read ``sn.issignal``, so the G-network feature names belong here and to no
other solver.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import os
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .algorithms import (
    AGResult,
    INAPAlgorithm,
    INAPPlusAlgorithm,
    INAPInfAlgorithm,
)
from ..avg_results import AvgResultsMixin
from ..base import NetworkSolver
from ...constants import default_verbose
from ...api.io.logging import LineError
from ...api.sn.getters import sn_get_arvr_from_tput
from ...api.sn.transforms import sn_get_residt_from_respt

# The RCAT family: every method this solver serves.
_AG_METHODS = ('default', 'inap', 'inapplus', 'inapinf', 'exact')

# Process types the RCAT phase construction can represent. An ALLOW-list on
# purpose: a type nobody has checked against this construction must be refused,
# not answered. Compared BY NAME, never by the raw enum value, because the
# numbering differs across codebases.
_RCAT_PROCESS_NAMES = frozenset((
    'EXP', 'ERLANG', 'HYPEREXP', 'PH', 'APH', 'COXIAN', 'COX2', 'MAP', 'MMPP2',
    'DET', 'UNIFORM', 'GAMMA', 'PARETO', 'WEIBULL', 'LOGNORMAL', 'REPLAYER',
    'IMMEDIATE', 'DISABLED',
))


def _proc_name(procid) -> Optional[str]:
    """The process type's NAME, whatever representation sn.procid holds."""
    if procid is None:
        return None
    if hasattr(procid, 'name'):
        return procid.name
    try:
        if np.isnan(float(procid)):
            return None
    except (TypeError, ValueError):
        return str(procid)
    from ...constants import ProcessType
    try:
        return ProcessType(int(procid)).name
    except (ValueError, TypeError):
        return None


def _has_signal_class(sn) -> bool:
    """True when the model declares at least one G-network signal class."""
    issignal = getattr(sn, 'issignal', None)
    if issignal is None:
        return False
    issignal = np.asarray(issignal, dtype=float)
    return issignal.size > 0 and bool(np.any(issignal > 0))


@dataclass
class SolverAGOptions:
    """Options for SolverAG.

    Attributes:
        method: 'default', 'inap', 'inapplus', 'inapinf' or 'exact'
        tol: Convergence tolerance of the reversed-rate fixed point. 1e-4 and
            100 sweeps are the MATLAB defaults (SolverOptions('AG') sets
            iter_max, iter_tol falls through to the global 1e-4); the fixed
            point is only linearly convergent, so a looser tolerance stops it
            ~1e-4 short of the limit and the two codebases must agree on where.
        max_iter: Maximum sweeps
        config: Per-run knobs. 'maxStates' truncates an open agent's
            queue-length dimension (default 100; 'inapinf' ignores it). 'exec'
            selects the execution backend -- 'serial' (default), 'parallel' (alias 'para') or
            'cluster' -- with 'nworkers' sizing the thread pool, 'endpoints'
            listing the ag-worker addresses and 'worker_timeout' bounding the
            wait on one. See exec_backend for why every backend produces the
            same iterates.
        verbose: Print debug information
    """
    method: str = 'default'
    tol: float = 1e-4
    max_iter: int = 100
    config: Dict[str, Any] = field(default_factory=dict)
    verbose: bool = field(default_factory=default_verbose)
    timeout: float = float('inf')
    lang: str = field(default_factory=lambda: os.environ.get('LINE_SOLVER_LANG', 'python'))
    # Arithmetic backend, lang='cpp' ONLY: 'double' (default), 'exact' or
    # 'real:<digits>'. Meaningless for the other langs, which are IEEE double
    # throughout, so line-cli is invoked without --arith unless the caller sets
    # it. `-s ag` carries every backend the model-solve path has, the RCAT fixed
    # point being linear algebra over the agent generators.
    arith: Optional[str] = None


class SolverAG(AvgResultsMixin, NetworkSolver):
    """Native Python agent-based (RCAT) solver."""

    ALGORITHMS = {
        'inap': INAPAlgorithm,
        'inapplus': INAPPlusAlgorithm,
        'inapinf': INAPInfAlgorithm,
    }

    def __init__(self, network, method: str = 'default',
                 options: Optional[SolverAGOptions] = None, **kwargs):
        self.network = network
        # The shared feature gate in NetworkSolver, and the solver console, both
        # read self.model, as the other native solvers do.
        self.model = network
        self.sn = self._get_network_struct(network)
        self._seed = kwargs.get('seed', None)

        if options is None:
            options = SolverAGOptions(method=method)
        elif method != 'default':
            options.method = method
        for key, value in kwargs.items():
            if hasattr(options, key):
                setattr(options, key, value)
        self.options = options
        self.result = None
        self.enableChecks = kwargs.get('enableChecks', True)

    def _get_network_struct(self, network):
        if hasattr(network, 'getStruct'):
            return network.getStruct()
        if hasattr(network, 'get_struct'):
            return network.get_struct()
        return network

    def reset(self):
        """Drop the cached result so the next getter re-solves."""
        self.result = None
        return self

    def getName(self) -> str:
        return 'SolverAG'

    @staticmethod
    def listValidMethods() -> List[str]:
        return list(_AG_METHODS)

    def resolveMethod(self, options=None):
        method = getattr(options or self.options, 'method', 'default')
        return 'inap' if method in ('default', 'exact') else method

    @staticmethod
    def getFeatureSet() -> set:
        """The RCAT feature envelope.

        The G-network names live here because the RCAT builder is the only code
        in LINE that reads sn.issignal.
        """
        return {
            'Sink', 'Source',
            'Fork', 'Join', 'Forker', 'Joiner',
            'Delay', 'DelayStation', 'Queue',
            'APH', 'Coxian', 'Erlang', 'Exp', 'HyperExp', 'MAP', 'MMPP2',
            'Det', 'Gamma', 'Lognormal', 'Pareto', 'Uniform', 'Weibull',
            'StatelessClassSwitcher', 'InfiniteServer',
            'SharedServer', 'Buffer', 'Dispatcher',
            'Server', 'JobSink', 'RandomSource', 'ServiceTunnel',
            'SchedStrategy_INF', 'SchedStrategy_PS',
            'SchedStrategy_FCFS',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND',
            'ClosedClass',
            # a self-looping class is its own single-station component
            # (solver_ag reads sn.isslc). MultiServer and FiniteCapacity are
            # deliberately absent: supportsModelMethod words both refusals
            # (RCAT drives rho = lambda/mu, and there is no buffer).
            'SelfLoopingClass',
            'OpenClass',
            'OpenSignal', 'ClosedSignal',
            'SignalType_NEGATIVE', 'SignalType_CATASTROPHE',
            'SignalBatchRemoval',
        }

    def getMethodFeatureSet(self, method):
        """Every AG method is the same decomposition, so they share one envelope.

        The genuine restrictions -- process type, single server -- are structural
        and are applied in supportsModelMethod.
        """
        return SolverAG.getFeatureSet()

    def supportsModelMethod(self, method) -> Tuple[bool, Optional[str]]:
        """Structural gate for the RCAT decomposition."""
        sn = self.sn
        procid = getattr(sn, 'procid', None)
        if procid is None:
            return True, None

        procid = np.asarray(procid, dtype=object)
        rates = np.asarray(getattr(sn, 'rates', np.zeros(procid.shape)), dtype=float)
        issignal = np.asarray(getattr(sn, 'issignal', np.zeros(procid.shape[1])),
                              dtype=float).ravel()

        nodetype = getattr(sn, 'nodetype', None)
        station_to_node = getattr(sn, 'stationToNode', None)

        def is_source(ist):
            if nodetype is None or station_to_node is None:
                return False
            from ...constants import NodeType
            node = int(np.asarray(station_to_node).ravel()[ist])
            nt = np.asarray(nodetype).ravel()[node]
            # NodeType is a plain Enum, not an IntEnum, so int() on a member
            # raises. Compare BY NAME when one is in hand and by identity
            # otherwise; never by an integer, which is also the cross-codebase
            # rule for every enum in sn.
            name = getattr(nt, 'name', None)
            if name is not None:
                return name == 'Source'
            return nt == NodeType.Source

        for ist in range(procid.shape[0]):
            for r in range(procid.shape[1]):
                rate = rates[ist, r] if ist < rates.shape[0] and r < rates.shape[1] else 0.0
                # Only a process that is actually in use can mis-answer.
                if not np.isfinite(rate) or rate <= 0:
                    continue
                signal = r < issignal.size and issignal[r] > 0
                # A signal is a trigger with no service, so its service entry is
                # never read; only its Source arrival rate is, and that one must
                # stay exponential because the removal is folded into the agent
                # as a scalar rate.
                if signal and not is_source(ist):
                    continue
                name = _proc_name(procid[ist, r])
                if signal:
                    if name is None or name == 'EXP':
                        continue
                    return False, (
                        'The %s method needs an exponential signal arrival process (a '
                        'removal signal is folded into the agent as a scalar rate), but '
                        'station %d class %d is %s. Use SolverMAM (method \'dec.source\') '
                        'for such models.' % (method, ist + 1, r + 1, name))
                if name is None or name in _RCAT_PROCESS_NAMES:
                    continue
                return False, (
                    'The %s method supports processes with a Markovian (D0,D1) '
                    'representation only (RCAT builds a phase dimension per agent out of '
                    'it), but station %d class %d is %s. Use SolverMAM (method '
                    '\'dec.source\') for such models.' % (method, ist + 1, r + 1, name))

        # RCAT models every station single-server; see _kb/06-solver-catalog.md
        nservers = np.asarray(getattr(sn, 'nservers', []), dtype=float).ravel()
        for ist in range(nservers.size):
            c = nservers[ist]
            if np.isfinite(c) and c > 1:
                return False, (
                    'The %s method supports single-server stations only (RCAT does not '
                    'model sn.nservers, so a multiserver station is driven at rho = '
                    'lambda/mu instead of lambda/(c*mu)), but station %d has %d servers. '
                    'Use SolverMAM (method \'dec.source\') for multiserver models.'
                    % (method, ist + 1, int(c)))

        # A FINITE BUFFER IS NOT SOMETHING RCAT CAN CARRY, and it was being
        # answered rather than refused: nothing under solver_ag reads sn.cap or
        # sn.classcap, so a capped station was decomposed as an unbounded one and
        # the table reported the UNCONSTRAINED figures (a closed 2-job tandem with
        # cap 1 on the second queue returned the same numbers with and without the
        # cap, 1.09 jobs in a buffer of 1). Lowering the component's level bound
        # (ag_inap: nlev = njobs[r]+1) to the buffer would not fix it: the
        # top-level boundary is a self-loop, which LOSES the arrival, whereas a
        # closed job refused at a full buffer must BLOCK the upstream departure,
        # and that coupling is exactly the independence RCAT assumes.
        # see _kb/06-solver-catalog.md
        cap_ok, cap_why = NetworkSolver.checkBindingCapacity(self.model, 'SolverAG')
        if not cap_ok:
            return False, cap_why

        return True, None

    @staticmethod
    def defaultOptions() -> SolverAGOptions:
        return SolverAGOptions()

    def runAnalyzer(self) -> 'SolverAG':
        """Run the reversed-rate fixed point and store the mean measures."""
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import populate_java_result
            populate_java_result(self)
            return self
        # see _kb/06-solver-catalog.md ("SolverAG under lang='cpp'"); an absent
        # binary is the only automatic fallback, a C++ refusal propagates.
        # `_CPP_SOLVER_TOKENS` carried 'AG' from the start, so the method name
        # resolved and the arm still could not be reached: without this branch
        # lang='cpp' ran the NATIVE solver and reported it as C++.
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import LineCliNotAvailable, populate_cpp_result
            try:
                populate_cpp_result(self)
                return self
            except LineCliNotAvailable as e:
                import warnings
                warnings.warn("SolverAG: lang='cpp' requested but the C++ solver is "
                              "unavailable (%s); falling back to lang='python'." % e)

        method = self.resolveMethod(self.options)

        if getattr(self, 'enableChecks', True):
            ok, reason = self.supportsModelMethod(method)
            if not ok:
                raise LineError(
                    'This model contains features not supported by the solver. ' + reason)

        start_time = time.time()

        algo_class = self.ALGORITHMS.get(method)
        if algo_class is None:
            raise ValueError(
                "Unknown AG method '%s'. Valid methods: %s"
                % (method, ', '.join(_AG_METHODS)))

        self.result = algo_class().solve(self.sn, self.options)

        # 'default' resolves to inap and 'exact' falls back to it with a warning
        # (AutoCAT is unreachable). Report what actually ran: 'exact' is
        # classified globally as an exact method, so leaving the name in place
        # would banner an iterative approximation as exact.
        self.result.method = method

        # A Source station's throughput is its arrival rate, not a solved value.
        if self.result.TN is not None and hasattr(self.sn, 'sched'):
            rates = getattr(self.sn, 'rates', None)
            if rates is not None:
                rates = np.asarray(rates, dtype=float)
                for ist in range(min(self.result.TN.shape[0], rates.shape[0])):
                    sched = self.sn.sched.get(ist, None) if hasattr(self.sn.sched, 'get') else None
                    name = getattr(sched, 'name', None)
                    if name == 'EXT':
                        self.result.TN[ist, :] = rates[ist, :]

        # THE ARRIVAL RATE IS A ROUTING TRANSFORM OF THE THROUGHPUTS, not a copy of
        # them. Leaving it unset made the table fall back to AN = TN, which is only
        # right when flow balances -- and an RCAT fixed point does NOT balance flow
        # on a closed network (INAP is stated for open ones), so the column read as
        # the station's own throughput while MATLAB and the JAR, which both apply
        # this transform, reported the upstream one. On the closed 2-queue cycle
        # that is 0.623 against MATLAB's 0.570 in the same cell.
        if self.result.TN is not None:
            self.result.AN = sn_get_arvr_from_tput(self.sn, self.result.TN)

        if self.result.RN is not None:
            self.result.WN = sn_get_residt_from_respt(self.sn, self.result.RN, None)

        self.result.runtime = time.time() - start_time
        return self

    def _ensureAvgResults(self):
        if self.result is None:
            self.runAnalyzer()

    def getAGResult(self) -> AGResult:
        """The result container, solving first if needed."""
        self._ensureAvgResults()
        return self.result


__all__ = ['SolverAG', 'SolverAGOptions', 'AGResult']
