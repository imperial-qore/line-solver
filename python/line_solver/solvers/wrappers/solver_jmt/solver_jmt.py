"""
JMT solver integration.

This implementation calls JMT via subprocess (command line), matching how
MATLAB's SolverJMT works. No JVM integration in Python itself.

The solver:
1. Writes the model to JSIM/JMVA XML format
2. Calls JMT via command line
3. Parses the result XML file

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np
import pandas as pd
import os
import re
import sys
import tempfile
import shutil
from typing import Optional, Dict, Any, List, Tuple, Set

from ....api.sn.transforms import sn_get_residt_from_respt
from ....api.io.logging import line_debug, line_warning, line_ack


class OptionsDict(dict):
    """A dict that supports attribute-style access."""
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(f"'OptionsDict' object has no attribute '{name}'")

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        try:
            del self[name]
        except KeyError:
            raise AttributeError(f"'OptionsDict' object has no attribute '{name}'")
from dataclasses import dataclass, field
from ....constants import default_verbose

from ....api.solvers.jmt.handler import (
    solver_jmt,
    SolverJMTOptions as _SolverJMTOptions,
    SolverJMTReturn,
    is_jmt_available,
    _get_jmt_jar_path,
)
from ...base import NetworkSolver, method_label, method_type


#: The JMVA algorithms that solve a CLOSED product-form network only.
#:
#: RECAL, CoMoM, Chow, Bard-Schweitzer (both spellings), AQL, Linearizer and De
#: Souza-Muntz Linearizer. Measured against JMT 1.2.x: each answers an open or a
#: mixed model with ``jmt.common.exception.UnsupportedModelException: The
#: selected solver cannot handle open classes, please choose another.`` and a
#: load-dependent one with the same exception naming load-dependent stations,
#: while the exact MVA engine behind 'jmva' and 'jmva.mva' serves both.
_JMVA_CLOSED_ONLY = frozenset((
    'jmva.amva', 'jmva.recal', 'jmva.comom', 'jmva.chow',
    'jmva.bs', 'jmva.aql', 'jmva.lin', 'jmva.dmlin',
))


def jmva_is_closed_only(method):
    """True for a JMVA algorithm restricted to closed single-server networks."""
    return str(method or '').lower() in _JMVA_CLOSED_ONLY


def jmt_method_refusal(sn, method, options=None):
    """The structural half of SolverJMT's method gate; '' when admissible.

    Three rules: a finite timespan for 'replication', single-server stations for
    the eight closed-form JMVA algorithms, and a binding finite buffer, which
    neither engine carries. None of them has a registry feature name.
    ONE PREDICATE, TWO CALLERS:
    ``supportsModelMethod`` asks it so findSolver and SolverAUTO never offer a
    pair that would die at run time, and the analyzer asks it again so a caller
    naming the method by hand gets the same sentence rather than a JMT stack
    trace. A second copy of either rule is how the gate and the run drift into
    two different answers.

    Everything a feature name CAN state lives in ``getMethodFeatureSet``
    instead: the JMVA envelope is narrower than the JSIM one, and the
    closed-only algorithms additionally drop OpenClass and LoadDependence.
    """
    method = str(method or '')
    if method.lower() == 'replication':
        # The transient arm integrates over [0, T]: a mean at an unstated
        # horizon is not a quantity. The horizon is an option, not a model
        # feature, so a feature set cannot see it.
        horizon = float('inf')
        if options is not None:
            horizon = getattr(options, 'max_simulated_time', float('inf'))
        try:
            horizon = float(horizon)
        except (TypeError, ValueError):
            horizon = float('inf')
        if not np.isfinite(horizon):
            return ("The replication method needs a finite timespan, e.g. "
                    "SolverJMT(model, timespan=[0, 10]).")
        return ''
    if jmva_is_closed_only(method):
        # JMVA implements these seven algorithms for SINGLE-SERVER stations
        # only, which is why the JMVA writer refuses the model rather than
        # emitting an <ldstation> the algorithm cannot read. A server count is
        # not a declared feature, so it cannot ride in the feature set the way
        # the load-dependent scaling of the same restriction does.
        nservers = getattr(sn, 'nservers', None)
        if nservers is not None:
            ns = np.asarray(nservers, dtype=float).ravel()
            ns = ns[np.isfinite(ns)]
            if ns.size and float(np.max(ns)) > 1.0:
                return '%s does not support multi-server stations.' % method
    return jmt_buffer_capacity_refusal(sn, method)


def jmt_buffer_capacity_refusal(sn, method):
    """A binding finite buffer, which NEITHER engine can carry; '' otherwise.

    The two engines fail it for opposite reasons, so the binding TEST is shared
    and the verdict is not.

    What makes a buffer BIND is not that ``sn.cap`` is finite: ``refresh_capacity``
    DERIVES a finite cap for every station nobody capped. It is that the cap is
    strictly below the population that can REACH the station, which is the JSIM
    writer's own test, and an infinite-server station has no buffer at all. Both
    are the writer's (``_jmt_reachable_population``), so the gate binds exactly
    where the writer binds.

    JSIM exports the buffer, but only for the rules JMT can read, and
    ``_jmt_station_cap_assert`` -- the writer's own predicate, asked here without
    letting it raise -- is what decides which. An open loss buffer and a declared
    BAS one stay runnable; only the cases JMT would answer unconstrained go.

    JMVA is refused OUTRIGHT: ``write_jmva`` emits a station type, a per-chain
    service demand and a per-chain visit count and nothing else, so the document
    has no capacity element for the buffer to ride in. Measured on a closed
    Delay+FCFS model, N=4, cap 2: every jmva method reported 2.19 jobs at a
    station that can hold 2, against the exact 1.33. This is the rule SolverMVA,
    SolverNC and SolverQNS already apply -- and SolverQNS writes THIS SAME
    DOCUMENT, so the jmva arm was the one hole in it.

    WHICH ENGINE IS ASKED ABOUT matters, and this port derives it from the method
    name because nothing else reaches the predicate: unlike MATLAB and the JAR,
    SolverQNS here writes its own JMVA document rather than borrowing
    SolverJMT.writeJMVA. A name SolverJMT does not implement therefore gets no
    verdict rather than JSIM's.
    """
    from ....api.solvers.jmt.handler import (
        _jmt_reachable_population, _jmt_station_cap_assert)
    name = str(method or '').lower()
    is_jmva = name.startswith('jmva')
    if not is_jmva and name not in ('default', 'jsim', 'replication'):
        return ''
    cap = getattr(sn, 'cap', None)
    if cap is None:
        return ''
    cap = np.asarray(cap, dtype=float).ravel()
    nservers = np.asarray(getattr(sn, 'nservers', []), dtype=float).ravel()
    from ....api.sn import NodeType
    nodetype = getattr(sn, 'nodetype', None)
    for ist in range(int(getattr(sn, 'nstations', 0))):
        # A SOURCE AND A SINK HAVE NO BUFFER THAT CAN BIND. The Source IS the
        # external world and the Sink absorbs, so neither ever holds a job a
        # capacity could refuse, yet refresh_capacity writes them a row like any
        # other station. Excluded on NODE TYPE, as the shared binding-capacity
        # gate of SolverMVA/SolverNC excludes them, and not by name.
        if nodetype is not None:
            ntype = nodetype[int(sn.stationToNode[ist])]
            if ntype in (NodeType.SOURCE, NodeType.SINK):
                continue
        # UNBOUNDED IS inf HERE. The JAR cannot carry inf on sn.cap -- its
        # Station.cap is an int whose "no bound" value is Integer.MAX_VALUE, and
        # refreshCapacity SUMS that sentinel across the classes served, so a
        # mixed station comes out as 2147483647 + N there and needs
        # SaveHandlers.jmtCapIsUnbounded. MATLAB, this port and C++ all default
        # station.cap to inf, so isfinite is the whole test.
        if ist >= cap.size or not np.isfinite(cap[ist]):
            continue
        if cap[ist] >= _jmt_reachable_population(sn, ist):
            continue
        if ist < nservers.size and np.isinf(nservers[ist]):
            continue
        if is_jmva:
            return ("Station %s carries a finite capacity %d that binds. The JMVA document "
                    "has no capacity element at all, so the analytical engine would solve the model "
                    "as if the buffer were unbounded and report that as the answer. Use the "
                    "'jsim' method, which exports the buffer with its drop rule when JMT can "
                    "express it, or SolverCTMC, SolverSSA or SolverLDES."
                    % (sn.nodenames[int(sn.stationToNode[ist])], int(cap[ist])))
        try:
            _jmt_station_cap_assert(sn, ist)
        except ValueError as exc:
            return str(exc)
    return ''


@dataclass
class SolverJMTOptions:
    """Options for the native JMT solver."""
    method: str = 'jsim'
    samples: int = 10000
    seed: int = 23000
    max_simulated_time: float = float('inf')
    conf_int: float = 0.99
    max_rel_err: float = 0.03
    verbose: bool = field(default_factory=default_verbose)
    keep: bool = False  # Keep temp files after execution
    timeout: float = float('inf')
    # Backend selection, see api/solvers/jmt/runner.py. rest_url points at a
    # JMT REST server (imperialqore/jmt-rest); container overrides the Docker
    # image used when no local JVM exists. Both empty means the local JVM.
    rest_url: Optional[str] = None
    container: Optional[str] = None
    # THIS FIELD WAS MISSING, and its absence is why `SolverJMT(model,
    # lang='cpp')` looked wired and was not: the kwarg went into **kwargs, was
    # never read, and every `options.lang` test in this file saw None. The
    # solve itself stays JSIM either way -- what lang='cpp' selects is WHICH
    # wrapper drives it, this one or `line-cli -s jmt`, so the two can be
    # compared. Same default resolution as every other solver.
    lang: str = field(default_factory=lambda: os.environ.get('LINE_SOLVER_LANG', 'python'))


class SolverJMT(NetworkSolver):
    """
    JMT solver integration.

    This solver provides discrete-event simulation and analytical methods
    via command line
    is launched as an external process, exactly like MATLAB's SolverJMT.

    Supported methods:
        - 'jsim' / 'default': Discrete event simulation
        - 'jmva' / 'jmva.mva': Mean Value Analysis
        - 'jmva.amva': Approximate MVA
        - 'jmva.recal': RECALsimulation
        - 'jmva.comom': CoMoM algorithm
        - 'jmva.chow': Chow algorithm
        - 'jmva.bs': Bard-Schweitzer
        - 'jmva.aql': AQL algorithm
        - 'jmva.lin': Linearizer
        - 'jmva.dmlin': De Souza-Muntz Linearizer

    Args:
        model: Network model (Python wrapper or native structure)
        method: Solution method (default: 'jsim')
        **kwargs: Additional solver options (samples, seed, etc.)

    Example:
        >>> solver = SolverJMT(model, samples=10000, seed=42)
        >>> solver.runAnalyzer()
        >>> table = solver.getAvgTable()
    """

    def __init__(self, model, method_or_options=None, **kwargs):
        self.model = model

        # see _kb/09-ldes-and-cache.md (Warm start) for initFromSolver contract
        init_solver = None
        if method_or_options is not None and hasattr(method_or_options, 'getAvgQLen'):
            init_solver = method_or_options
            method_or_options = None

        # Handle options passed as second argument (MATLAB-style)
        if method_or_options is None:
            # honor method= keyword (consistent with SolverNC/SolverMVA)
            self.method = str(kwargs.pop('method', 'jsim')).lower()
        elif isinstance(method_or_options, str):
            self.method = method_or_options.lower()
        elif hasattr(method_or_options, 'get'):
            # Dict-like options object
            self.method = method_or_options.get('method', 'jsim')
            if 'samples' in method_or_options:
                kwargs.setdefault('samples', method_or_options['samples'])
            if 'seed' in method_or_options:
                kwargs.setdefault('seed', method_or_options['seed'])
            if 'keep' in method_or_options:
                kwargs.setdefault('keep', method_or_options['keep'])
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
        elif hasattr(method_or_options, 'method'):
            # SolverOptions-like object
            self.method = getattr(method_or_options, 'method', 'jsim')
            if hasattr(method_or_options, 'samples'):
                kwargs.setdefault('samples', method_or_options.samples)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
        else:
            self.method = 'jsim'

        # Parse options
        samples = kwargs.get('samples', 10000)
        seed = kwargs.get('seed', 23000)
        verbose = kwargs.get('verbose', default_verbose())
        keep = kwargs.get('keep', False)
        conf_int = kwargs.get('conf_int', kwargs.get('confint', 0.99))
        max_rel_err = kwargs.get('max_rel_err', 0.03)
        max_simulated_time = kwargs.get('max_simulated_time',
                                        kwargs.get('timespan', [0, float('inf')])[1]
                                        if isinstance(kwargs.get('timespan'), list) else float('inf'))
        timeout = kwargs.get('timeout', float('inf'))
        rest_url = kwargs.get('rest_url', None)
        container = kwargs.get('container', None)
        lang = kwargs.get('lang', os.environ.get('LINE_SOLVER_LANG', 'python'))

        self.options = SolverJMTOptions(
            lang=lang,
            method=self.method,
            samples=samples,
            seed=seed,
            max_simulated_time=max_simulated_time,
            conf_int=conf_int,
            max_rel_err=max_rel_err,
            verbose=verbose,
            keep=keep,
            timeout=timeout,
            rest_url=rest_url,
            container=container
        )

        self._result: Optional[SolverJMTReturn] = None
        self._sn = None

        # Extract network structure
        self._extract_network_params()

        if init_solver is not None:
            self.initFromSolver(init_solver)

    def getName(self) -> str:
        """Get the name of this solver."""
        return "JMT"

    get_name = getName

    def _extract_network_params(self):
        """Extract parameters from the model."""
        model = self.model

        # see _kb/06-solver-catalog.md (Wrappers: "JMT python wrapper: export/import internals")
        if hasattr(model, 'refresh_struct'):
            model.refresh_struct()
            if hasattr(model, '_sn') and model._sn is not None:
                self._sn = model._sn
                return

        # Fallback: Use existing _sn if refresh_struct is not available
        if hasattr(model, '_sn') and model._sn is not None:
            self._sn = model._sn
            return

        # Native model (snake-case get_struct()); no JAR-wrapper bridge here.
        if hasattr(model, 'get_struct'):
            self._sn = model.get_struct()
            if self._sn is not None:
                return

        # Already a native NetworkStruct
        if hasattr(model, 'nclasses') and hasattr(model, 'nstations'):
            self._sn = model
            return

        raise ValueError(
            "Cannot extract a native NetworkStruct from model. Native solvers "
            "accept only native Network / NetworkStruct inputs (no JAR wrapper).")

    def supportsTransientAnalysis(self):
        """Transient averages are available (simulation restricted to options.timespan)."""
        return True

    supports_transient_analysis = supportsTransientAnalysis

    def runAnalyzer(self) -> 'SolverJMT':
        """
        Run the JMT analyzer.

        Calls JMT via command line and stores the results.

        Returns:
            self for method chaining
        """
        line_ack('JMT', self.options.verbose)
        line_debug("JMT analyzer starting: lang=python, samples=%d, seed=%d",
                   self.options.samples, self.options.seed, options=self.options)

        if self._sn is None:
            raise RuntimeError("Network structure not available")

        # see _kb/06-solver-catalog.md (Wrappers: "JMT python wrapper: export/import internals")
        model = getattr(self, 'model', None)
        if model is not None and hasattr(model, 'get_used_lang_features'):
            self.runAnalyzerChecks(self.options)

        # The structural half of the gate, asked again here so a caller who
        # reaches the analyzer with the checks disabled still gets the gate's
        # own sentence rather than a JMT stack trace.
        structural = jmt_method_refusal(self._sn, self.options.method, self.options)
        if structural:
            raise RuntimeError(structural)

        if getattr(self._sn, 'immfeed', None) is not None and np.any(self._sn.immfeed):
            line_warning("SolverJMT", "SolverJMT does not support immediate feedback (immfeed); no solution returned.")
            return self

        method = self.options.method
        if method == 'replication':
            # TRANSIENT AVERAGES BY INDEPENDENT REPLICATION, the reference's own
            # transient route for JMT: a single sample path is not the transient
            # mean E[N](t), there being no time-ergodicity at a fixed t, so
            # iter_max seeded replications are sampled and averaged onto a common
            # time grid. Port of the 'replication' arm of
            # @SolverJMT/runAnalyzer.m, which python did not carry at all.
            return self._runReplication()
        if method in ('jsim', 'default'):
            line_debug("JMT: using JSIM method (discrete-event simulation), samples=%d, seed=%d",
                       self.options.samples, self.options.seed, options=self.options)
        elif method.startswith('jmva'):
            line_debug("JMT: using JMVA method: %s", method, options=self.options)
        else:
            line_debug("JMT: using method: %s", method, options=self.options)

        if self.options.samples < 5000:
            line_debug("JMT: sample size adjusted to minimum 5000", options=self.options)

        # Convert options to handler format
        handler_options = _SolverJMTOptions(
            method=self.options.method,
            samples=self.options.samples,
            seed=self.options.seed,
            max_simulated_time=self.options.max_simulated_time,
            conf_int=self.options.conf_int,
            max_rel_err=self.options.max_rel_err,
            verbose=self.options.verbose,
            keep=self.options.keep,
            timeout=getattr(self.options, 'timeout', float('inf')),
            rest_url=getattr(self.options, 'rest_url', None),
            container=getattr(self.options, 'container', None)
        )

        # Call the handler (pass model for FCR region support)
        self._result = solver_jmt(self._sn, handler_options, self.model)

        # Print completion message (matches MATLAB verbose guard)
        if self.options.verbose:
            py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
            runtime = self._result.runtime if hasattr(self._result, 'runtime') else 0.0
            method = self._result.method if hasattr(self._result, 'method') else self.options.method
            from line_solver.solvers.base import print_solver_banner
            print_solver_banner(f"JMT analysis [method: {method_label(self.options.method, method)}; type: {method_type('JMT', method_label(self.options.method, method))}; lang: python; env: {py_version}] completed in {runtime:.6f}s.")

        return self

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get average performance metrics as a DataFrame.

        Returns:
            DataFrame with columns: Station, Class, QLen, Util, RespT, Tput, ArvR
        """
        if self._result is None:
            self._ensureAvgResults()

        M = self._sn.nstations
        K = self._sn.nclasses

        nodenames = self._sn.nodenames if self._sn.nodenames else [f'Station{i}' for i in range(M)]
        classnames = self._sn.classnames if self._sn.classnames else [f'Class{r}' for r in range(K)]

        # Compute ResidT using proper visit ratios from network structure
        # This uses the correct formula: WN[ist,k] = RN[ist,k] * V[ist,k] / V[refstat,refclass]
        if self._sn is not None and self._sn.visits and self._result.R is not None:
            WN = sn_get_residt_from_respt(self._sn, self._result.R, None)
        else:
            # Fallback: ResidT = RespT (no visit information available)
            WN = self._result.R.copy() if self._result.R is not None else np.zeros((M, K))

        # Get station names and identify source stations
        station_names = []
        source_stations = set()
        nodetype = self._sn.nodetype if hasattr(self._sn, 'nodetype') else None

        for i in range(M):
            node_idx = int(self._sn.stationToNode[i]) if self._sn.stationToNode is not None else i
            if node_idx < len(nodenames):
                station_names.append(nodenames[node_idx])
            else:
                station_names.append(f'Station{i}')

            # Check if source station
            if nodetype is not None and node_idx < len(nodetype):
                if int(nodetype[node_idx]) == 0:  # SOURCE = 0
                    source_stations.add(i)

        # Get arrival rates from rates matrix for source stations
        rates = np.asarray(self._sn.rates) if hasattr(self._sn, 'rates') and self._sn.rates is not None else None

        rows = []
        for i in range(M):
            for r in range(K):
                is_source = i in source_stations

                # Get values with NaN handling
                qlen = self._result.Q[i, r] if self._result.Q is not None else np.nan
                util = self._result.U[i, r] if self._result.U is not None else np.nan
                respt = self._result.R[i, r] if self._result.R is not None else np.nan
                residt = WN[i, r] if i < WN.shape[0] and r < WN.shape[1] else respt
                arvr = self._result.A[i, r] if self._result.A is not None else np.nan
                tput = self._result.T[i, r] if self._result.T is not None else np.nan

                # For source stations, replace NaN with 0 and set Tput to arrival rate
                if is_source:
                    qlen = 0.0 if np.isnan(qlen) else qlen
                    util = 0.0 if np.isnan(util) else util
                    respt = 0.0 if np.isnan(respt) else respt
                    residt = 0.0 if np.isnan(residt) else residt
                    arvr = 0.0  # Source has no arrivals to itself

                    # Set Tput from arrival rate
                    if np.isnan(tput) and rates is not None:
                        stationToNode = np.asarray(self._sn.stationToNode).flatten()
                        node_idx = int(stationToNode[i])
                        if node_idx < rates.shape[0] and r < rates.shape[1]:
                            tput = rates[node_idx, r]

                # Filter out rows where all metrics are zero or NaN (matching MATLAB behavior)
                # Only include row if at least one metric is non-zero and not NaN
                metrics = [qlen, util, respt, residt, arvr, tput]
                has_significant_value = any(
                    (not np.isnan(v) and v > 0) for v in metrics
                )
                if not has_significant_value:
                    continue

                rows.append({
                    'Station': station_names[i],
                    'JobClass': classnames[r],
                    'QLen': qlen,
                    'Util': util,
                    'RespT': respt,
                    'ResidT': residt,
                    'ArvR': arvr,
                    'Tput': tput,
                })

        df = pd.DataFrame(rows)

        if not getattr(self, '_table_silent', False):
            print(df.to_string(index=False))

        return df

    def getAvgQLen(self) -> np.ndarray:
        """Get average queue lengths (M x K matrix)."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.Q if self._result.Q is not None else np.array([])

    def getAvgUtil(self) -> np.ndarray:
        """Get average utilizations (M x K matrix)."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.U if self._result.U is not None else np.array([])

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times (M x K matrix)."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.R if self._result.R is not None else np.array([])

    def getAvgTput(self) -> np.ndarray:
        """Get average throughputs (M x K matrix)."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.T if self._result.T is not None else np.array([])

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates (M x K matrix)."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.A if self._result.A is not None else np.array([])

    def getAvgFcr(self) -> np.ndarray:
        """The finite-capacity-region rows, ((nregions*6) x nclasses) or empty.

        WHY A SEPARATE ACCESSOR. `getAvg` returns the STATION metrics alone, and
        a region is not a station: its rows live past the last one, which is
        where `getAvgNodeTable` reads them from to fill the FCR pseudo-node. A
        host bridging through `getAvg` (MATLAB `lang='python'`) therefore saw no
        region at all and dropped the FCR row from its node table
        (fcr_mm1waitq[M2P], "row FCR1 missing").

        The six blocks are stacked in the order Q, U, R, W, A, T, each
        (nregions x nclasses), so one marshalled matrix carries all of them.
        Util and ArvR are NaN: JMT reports neither for a region.
        """
        if self._result is None:
            self._ensureAvgResults()
        F = int(getattr(self._sn, 'nregions', 0) or 0)
        Qfcr = getattr(self._result, 'Qfcr', None)
        if F <= 0 or Qfcr is None:
            return np.array([])
        K = int(self._sn.nclasses)
        nan = np.full((F, K), np.nan)
        zero = np.zeros((F, K))

        def blk(name, default):
            v = getattr(self._result, name, None)
            return default if v is None else np.asarray(v, dtype=float).reshape(F, K)

        return np.vstack([blk('Qfcr', zero), nan, blk('Rfcr', zero),
                          blk('Wfcr', zero), nan, blk('Tfcr', zero)])

    def getAvgChainTable(self) -> pd.DataFrame:
        """
        Get average performance metrics aggregated by chain.

        Returns:
            DataFrame with columns: Chain, QLen, Util, RespT, Tput
        """
        if self._result is None:
            self._ensureAvgResults()

        # Get chain information from model structure
        nchains = self._sn.nchains if hasattr(self._sn, 'nchains') else self._sn.nclasses
        inchain = self._sn.inchain if hasattr(self._sn, 'inchain') else None

        rows = []
        for c in range(nchains):
            chain_name = f'Chain{c+1}'

            # Get classes in this chain
            if inchain is not None and c in inchain:
                chain_classes = inchain[c].flatten().astype(int)
            else:
                chain_classes = [c]  # Single class per chain

            # Aggregate metrics across stations and classes in chain
            total_qlen = 0.0
            total_util = 0.0
            total_respt = 0.0
            total_tput = 0.0

            M = self._sn.nstations
            for i in range(M):
                for k in chain_classes:
                    if k < self._result.Q.shape[1]:
                        total_qlen += self._result.Q[i, k] if not np.isnan(self._result.Q[i, k]) else 0.0
                        total_util += self._result.U[i, k] if not np.isnan(self._result.U[i, k]) else 0.0
                        total_respt += self._result.R[i, k] if not np.isnan(self._result.R[i, k]) else 0.0
                        total_tput = max(total_tput, self._result.T[i, k] if not np.isnan(self._result.T[i, k]) else 0.0)

            rows.append({
                'Chain': chain_name,
                'QLen': total_qlen,
                'Util': total_util,
                'RespT': total_respt,
                'Tput': total_tput,
            })

        # five SIGNIFICANT digits like MATLAB's table, not pandas' five decimals
        from line_solver.indexed_table import IndexedTable
        return IndexedTable(pd.DataFrame(rows))

    def getAvgSysTable(self) -> pd.DataFrame:
        """
        Get system-level average performance metrics.

        Returns:
            DataFrame with columns: Chain, SysRespT, SysTput
        """
        if self._result is None:
            self._ensureAvgResults()

        chain_table = self.getAvgChainTable()
        CN = []
        XN = []
        for _, row in chain_table.iterrows():
            CN.append(row['RespT'])
            XN.append(row['Tput'])
        return self._make_sys_table(CN, XN)

    def getAvgSysRespT(self) -> np.ndarray:
        """Get system response times (1 x K)."""
        if self._result is None:
            self._ensureAvgResults()
        # Sum response times across all stations for each class
        if self._result.R is not None:
            return np.nansum(self._result.R, axis=0, keepdims=True)
        return np.array([[]])

    def getAvgSysTput(self) -> np.ndarray:
        """Get system throughputs (1 x K)."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.X if self._result.X is not None else np.array([[]])

    def sampleSysAggr(self, num_events: Optional[int] = None) -> Optional[Dict[str, Any]]:
        """Sample system-wide aggregated state trajectories via JMT logging.

        Faithful port of MATLAB SolverJMT.sampleSysAggr: all non-Source stations
        are logged in a temporary model copy, simulated, and their per-class
        queue-length trajectories reconstructed and interpolated (previous/step)
        onto a common timeline (the union of all stations' event times).

        Args:
            num_events: Number of events to sample (default: uses options.samples)

        Returns:
            Dict with 't' (common timeline), 'state' (list of per-station
            len(t) x nclasses matrices; Source stations = [[inf]]), 'handle',
            'event' (chronological event list) and 'isaggregate'=True.
        """
        import os
        from ....api.solvers.jmt.handler import parse_tran_state
        from ....lang.base import NodeType

        sn = self._sn if self._sn is not None else self.model.get_struct()
        if num_events is None:
            num_events = getattr(self.options, 'samples', 10000)
        K = sn.nclasses
        M = sn.nstations

        def _is_source(ind):
            nt = sn.nodetype[ind]
            nt_val = nt.value if hasattr(nt, 'value') else int(nt)
            src_val = NodeType.SOURCE.value if hasattr(NodeType.SOURCE, 'value') else int(NodeType.SOURCE)
            return nt_val == src_val

        nnodes = self.model.get_number_of_nodes()
        is_node_logged = [False] * nnodes
        station_node = [int(sn.stationToNode[ist]) for ist in range(M)]
        for ist in range(M):
            ind = station_node[ist]
            if not _is_source(ind):
                is_node_logged[ind] = True

        model_copy, log_path = self._run_logged_copy(is_node_logged, num_events)
        class_names = [sn.classnames[r] for r in range(K)] if sn.classnames is not None else None

        stat_t = [None] * M
        stat_q = [None] * M
        all_events = []
        for ist in range(M):
            ind = station_node[ist]
            if not is_node_logged[ind]:
                continue
            name = sn.nodenames[ind]
            preload = self._node_preload(ind, K)
            arv = os.path.join(log_path, f"{name}-Arv.csv")
            dep = os.path.join(log_path, f"{name}-Dep.csv")
            state, evtype, evclass, evjob = parse_tran_state(arv, dep, preload, class_names)
            _, uniq = np.unique(state[:, 0], return_index=True)
            uniq = np.sort(uniq)
            stat_t[ist] = state[uniq, 0]
            stat_q[ist] = state[uniq, 1:1 + K]
            for e in range(len(evtype)):
                if not np.isnan(evjob[e]):
                    all_events.append({
                        'event': evtype[e], 'node': ind,
                        'class': int(evclass[e]) if not np.isnan(evclass[e]) else None,
                        't': float(state[e, 0]),
                        'job': int(evjob[e]) if not np.isnan(evjob[e]) else None,
                    })

        # Common timeline = union of all stations' timestamps, capped at the
        # earliest station end so previous-interpolation never extrapolates.
        t_union = np.array([])
        maxes = []
        for ist in range(M):
            if stat_t[ist] is not None and len(stat_t[ist]) > 0:
                t_union = np.union1d(t_union, stat_t[ist])
                maxes.append(np.max(stat_t[ist]))
        if maxes:
            t_union = t_union[t_union <= min(maxes)]

        def _prev_interp(ts, ys, tq):
            if len(ts) == 0:
                return np.full(len(tq), np.nan)
            idx = np.searchsorted(ts, tq, side='right') - 1
            idx = np.clip(idx, 0, len(ys) - 1)
            return ys[idx]

        state_list = []
        for ist in range(M):
            ind = station_node[ist]
            if _is_source(ind) or stat_t[ist] is None:
                state_list.append(np.array([[np.inf]]))
                continue
            cols = np.zeros((len(t_union), K))
            for r in range(K):
                cols[:, r] = _prev_interp(stat_t[ist], stat_q[ist][:, r], t_union)
            state_list.append(cols)

        all_events.sort(key=lambda ev: ev['t'])

        handles = [self.model.get_stations()[ist] if hasattr(self.model, 'get_stations')
                   else ist for ist in range(M)]
        self._cleanup_log_dir(log_path)
        return {
            'handle': handles,
            't': t_union,
            'state': state_list,
            'event': all_events,
            'isaggregate': True,
        }

    def getProbSysAggr(self) -> float:
        """Get system state probability via simulation sampling.

        Uses JMT simulation with logging to estimate the probability of the
        system being in the current aggregated state (as set via setState).

        The probability is computed as the fraction of time the system spends
        in the target state during simulation.

        Note: This method requires Logger nodes to be present in the model for
        full functionality. If log files are not available, it attempts to
        estimate probabilities from the simulation's average metrics.

        Returns:
            float: Estimated probability of the current system state.
                   Returns 0.0 if the state was not observed during simulation.
        """
        if getattr(self.options, 'lang', 'python') == 'cpp':
            import warnings
            from ...cpp_dispatch import jmt_prob_aggr_via_cpp
            r = jmt_prob_aggr_via_cpp(self)
            if not r['sysStateSeen']:
                warnings.warn("the system state was not seen in the simulation, "
                              "so its probability is reported as 0")
            return r['probSysAggr']

        if self._sn is None:
            self._extract_network_params()

        sn = self._sn

        # Try to get state samples from simulation with logging
        sample_result = self.sampleSysAggr()

        if sample_result is not None and 't' in sample_result and 'state' in sample_result:
            timestamps = sample_result['t']
            states = sample_result['state']

            if len(timestamps) >= 2:
                # Get target state from the model's current state
                nstations = sn.nstations
                nclasses = sn.nclasses
                target_state = np.zeros((nstations, nclasses))

                # Get state from sn.state dict
                if hasattr(sn, 'state') and sn.state is not None:
                    for ist in range(nstations):
                        node_idx = int(sn.stationToNode[ist]) if sn.stationToNode is not None else ist
                        stateful_idx = int(sn.nodeToStateful[node_idx]) if sn.nodeToStateful is not None and node_idx < len(sn.nodeToStateful) else -1
                        if stateful_idx >= 0 and isinstance(sn.state, dict):
                            for key, state_vec in sn.state.items():
                                key_idx = -1
                                if hasattr(key, 'statefulIndex'):
                                    key_idx = key.statefulIndex
                                elif hasattr(key, 'getStatefulIndex'):
                                    key_idx = key.getStatefulIndex()
                                if key_idx == stateful_idx:
                                    if state_vec is not None:
                                        state_arr = np.asarray(state_vec).flatten()
                                        for r in range(min(nclasses, len(state_arr))):
                                            target_state[ist, r] = state_arr[r]
                                    break

                # see _kb/06-solver-catalog.md (Wrappers: "JMT python wrapper:
                # export/import internals") for the vectorized snapshot matching
                time_diffs = np.diff(timestamps)
                total_time = np.sum(time_diffs)

                if total_time > 0:
                    n_steps = len(time_diffs)
                    tol = 0.5 + 1e-5 * np.abs(target_state)
                    match = np.ones(n_steps, dtype=bool)
                    for ist in range(nstations):
                        si = states[ist] if ist < len(states) else None
                        row = np.zeros((n_steps, nclasses))
                        if si is not None and np.all(np.isfinite(si)):
                            # Source station ([inf]) keeps the zero row.
                            rows_avail = min(n_steps, si.shape[0])
                            ncol = min(nclasses, si.shape[1])
                            row[:rows_avail, :ncol] = si[:rows_avail, :ncol]
                        match &= np.all(
                            np.abs(row - target_state[ist]) <= tol[ist], axis=1)
                    matching_time = float(np.sum(time_diffs[match]))

                    if matching_time > 0:
                        return matching_time / total_time

        return 0.0

    def getProbAggr(self, station: int) -> float:
        """Get the aggregated state probability at a station.

        Under lang='cpp' this is `-s jmt -a prob`: one logged JSIM run,
        dwell-weighted over the time the station holds its declared state. The
        native path here has no equivalent -- JSIM reports means and not state
        occupancies to this wrapper -- and returns an empty array, which is what
        it has always done.
        """
        if getattr(self.options, 'lang', 'python') == 'cpp':
            import warnings
            from ...cpp_dispatch import jmt_prob_aggr_via_cpp
            r = jmt_prob_aggr_via_cpp(self)
            ist = int(station)
            if not (0 <= ist < r['probAggr'].size):
                raise ValueError("station index %r is outside 0..%d"
                                 % (station, r['probAggr'].size - 1))
            if ist < len(r['stateSeen']) and not r['stateSeen'][ist]:
                warnings.warn("station %d's state was not seen in the simulation, "
                              "so its probability is reported as 0" % ist)
            return float(r['probAggr'][ist])
        if self._result is None:
            self._ensureAvgResults()
        # Not supported in simulation - return empty
        return np.array([])

    def getRuntime(self) -> float:
        """Get solver runtime in seconds."""
        if self._result is None:
            return 0.0
        return self._result.runtime

    def getMethod(self) -> str:
        """Get the method used."""
        if self._result is None:
            return self.method
        return self._result.method

    def listValidMethods(self) -> List[str]:
        """List valid methods for this solver."""
        return [
            'default', 'jsim',
            # TRANSIENT AVERAGES BY INDEPENDENT REPLICATION, the reference's own
            # transient route for JMT; ported and dispatched in runAnalyzer.
            'replication',
            'jmva', 'jmva.mva', 'jmva.amva', 'jmva.recal',
            'jmva.comom', 'jmva.chow', 'jmva.bs', 'jmva.aql',
            'jmva.lin', 'jmva.dmlin'
        ]

    def isStochasticMethod(self, method):
        """Simulation-based methods (default, jsim, replication) return
        stochastic estimates. The analytical JMVA methods do not, except
        for the sampling-based variants (e.g. jmva.ls).
        """
        if not method:
            return True  # default resolves to jsim simulation
        method_names = re.split(r'[./]', str(method).lower())
        if 'jmva' in method_names:
            return any(tok in ('ls', 'mci', 'imci', 'sampling') for tok in method_names)
        return True

    is_stochastic_method = isStochasticMethod

    @staticmethod
    def isAvailable() -> bool:
        """Check if JMT solver is available."""
        return is_jmt_available()

    @staticmethod
    def getFeatureSet() -> Set[str]:
        """Get the set of features supported by this solver."""
        return {
            'Sink', 'Source', 'Router', 'ClassSwitch',
            'Delay', 'DelayStation', 'Queue',
            'Fork', 'Join', 'Forker', 'Joiner', 'Logger',
            'JoinPartial',  # quorum join, written out as a jmt PartialJoin
            # A variable forking level: saveForkStrategy turns isSimplifiedFork
            # off and writes the per-branch entries, so jmt reads the counts,
            # the probabilities and the degree distribution rather than sending
            # one job down every link.
            'ForkFanoutVector', 'ForkFanoutRandom', 'ForkBranchProbability',
            'Coxian', 'Cox2', 'APH', 'Erlang', 'Exp', 'HyperExp',
            'Det', 'Gamma', 'Lognormal', 'MAP', 'MMPP2',
            'Normal', 'PH', 'Pareto', 'Weibull', 'Replayer', 'Uniform',
            'StatelessClassSwitcher', 'InfiniteServer', 'SharedServer',
            'Buffer', 'Dispatcher', 'Server', 'JobSink', 'RandomSource',
            'ServiceTunnel', 'LogTunnel', 'Linkage',
            'Enabling', 'Inhibiting', 'Timing', 'Firing', 'Storage', 'Place', 'Transition',
            'SchedStrategy_INF', 'SchedStrategy_PS', 'SchedStrategy_DPS',
            'SchedStrategy_FCFS', 'SchedStrategy_GPS', 'SchedStrategy_SIRO',
            'SchedStrategy_HOL', 'SchedStrategy_PSPRIO', 'SchedStrategy_DPSPRIO',
            'SchedStrategy_GPSPRIO', 'SchedStrategy_LCFS', 'SchedStrategy_LCFSPR',
            'SchedStrategy_LCFSPRIO', 'SchedStrategy_LCFSPRPRIO',
            # LCFSPI is emitted by the writer (QueuePutStrategies.LCFSPIStrategy)
            # and was not declared. The FCFS preemptive family is NOT added: this
            # writer has no FCFSPRStrategy/FCFSPIStrategy branch, unlike MATLAB,
            # the JAR and C++, so declaring it would promise an export that then
            # falls through to the non-preemptive tail.
            'SchedStrategy_LCFSPI',
            'SchedStrategy_SEPT', 'SchedStrategy_SRPT', 'SchedStrategy_SRPTPRIO', 'SchedStrategy_LEPT',
            'SchedStrategy_SJF', 'SchedStrategy_LJF', 'SchedStrategy_LPS',
            'SchedStrategy_POLLING', 'SchedStrategy_EXT',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND',
            'RoutingStrategy_RROBIN', 'RoutingStrategy_WRROBIN',
            'RoutingStrategy_JSQ',
            'RoutingStrategy_SQ',
            'ClosedClass', 'SelfLoopingClass', 'OpenClass',
            # Cache, CacheClassSwitcher and the four replacement strategies are
            # NOT declared here, and this port is the only one that withholds
            # them. MATLAB (saveCacheStrategy.m), the JAR (SaveHandlers) and C++
            # (jmt_writer.h) all serialize a Cache node; this writer has no
            # cache branch at all -- its module header says so -- so it emits
            # <node name="Cache"/> with no sections and jsim dies inside JMT
            # with "Cannot invoke NodeSection.updateVisitPath ... inputSection
            # is null". Declaring a name the writer cannot emit promises an
            # export that is not there, which is the same rule that keeps the
            # FCFS preemptive family out above.
            'Region',
            # Limited load dependence reaches JMT only as a SERVER COUNT: the
            # JSIM writer exports max(nservers, max(alpha)) and the JMVA writer
            # the matching <ldstation>. That is exact for alpha(n) = min(n,c)
            # and for nothing else, so supportsModelMethod refuses any other
            # scaling by name.
            'LoadDependence',
            # Exported as delayOffTime/setUpTime (_write_delayoff_strategy).
            'SetupDelayOff',
            # Exported as classParallelism (_write_class_parallelism).
            'ServerParallelism',
            # Heterogeneous server pools: the type names, the servers per type
            # and the compatibility matrix are exported as serverTypesNames /
            # serverTypesNumOfServers / serverTypesCompatibilities, so jsim
            # simulates the pools rather than a station of the same total size.
            'HeteroServers',
            # Exported as Impatience/Reneging and Impatience/Balking strategies.
            'Reneging', 'Balking',
            # c-server stations (the writer emits numberOfServers) and finite
            # buffers with their drop rule (a capacity plus the dropStrategy
            # text): jsim exports both, and the structural capacity gate keeps
            # refusing the buffers JMT would answer unconstrained (closed WAITQ,
            # BBS, RSRD, retrial-with-limit); the JMVA set withdraws the buffer.
            # 'Retrial' is NOT declared, and this port is the only one that
            # withholds it: MATLAB, the JAR and C++ pick the retrial Queue
            # constructor and write the per-class orbit delay, while this writer
            # has no such branch, so the orbit would simply not be exported.
            'MultiServer', 'FiniteCapacity',
        }

    @staticmethod
    def getJMVAFeatureSet() -> Set[str]:
        """What the JMVA ANALYTICAL engine accepts, much less than JSIM.

        Derived from the writer rather than guessed: ``write_jmva`` emits, per
        station, a ``<delaystation>``, a ``<listation>`` or an ``<ldstation>``,
        a per-chain ``<servicetime>`` and a per-chain ``<visit>``, and at model
        level the closed populations, the open arrival rates and the reference
        station. NOTHING ELSE IN THE MODEL REACHES JMVA, so a construct whose
        whole effect is not carried by (station type, demand, visits,
        population) would be solved away silently.

        Dropped from the JSIM set, and why:
          * Fork/Join and the fan-out names -- no fork element exists, and a
            visit ratio cannot express the join synchronization.
          * Place/Transition and the Petri-net sections -- no counterpart.
          * Region -- JMVA has no finite capacity region.
          * Reneging/Balking -- no impatience element; the abandonment would
            simply not happen.
          * SetupDelayOff, ServerParallelism, HeteroServers -- each a
            server-side attribute the JMVA document has no slot for.
          * the non-BCMP disciplines -- the writer emits NO discipline at all,
            so a priority, weighted, size-based or limited-sharing station
            would be solved as an ordinary load-independent one. Only the four
            BCMP station types survive the encoding, the same line SolverNC and
            SolverMVA draw.
          * the state-dependent routings (RROBIN, WRROBIN, JSQ, SQ) -- the
            document carries mean visit counts, which is not what makes a
            join-the-shortest-queue model behave as it does.

        The DISTRIBUTIONS are deliberately kept: JMVA consumes a mean service
        demand, so any renewal law with a finite mean is admissible, exactly as
        it is for SolverMVA and SolverNC. Cache is absent from this port's JSIM
        set already and so does not appear here either.
        """
        return SolverJMT.getFeatureSet() - {
            'Fork', 'Join', 'Forker', 'Joiner', 'JoinPartial',
            'ForkFanoutVector', 'ForkFanoutRandom', 'ForkBranchProbability',
            'Place', 'Transition', 'Enabling', 'Inhibiting', 'Timing',
            'Firing', 'Storage',
            'Region',
            'Reneging', 'Balking',
            'SetupDelayOff', 'ServerParallelism', 'HeteroServers',
            'SchedStrategy_DPS', 'SchedStrategy_GPS', 'SchedStrategy_HOL',
            'SchedStrategy_PSPRIO', 'SchedStrategy_DPSPRIO', 'SchedStrategy_GPSPRIO',
            'SchedStrategy_LCFSPI', 'SchedStrategy_LCFSPIPRIO',
            'SchedStrategy_LCFSPRIO', 'SchedStrategy_LCFSPRPRIO',
            'SchedStrategy_FCFSPR', 'SchedStrategy_FCFSPI',
            'SchedStrategy_FCFSPRPRIO', 'SchedStrategy_FCFSPIPRIO',
            'SchedStrategy_SEPT', 'SchedStrategy_LEPT',
            'SchedStrategy_SJF', 'SchedStrategy_LJF',
            'SchedStrategy_SRPT', 'SchedStrategy_SRPTPRIO',
            'SchedStrategy_LPS', 'SchedStrategy_POLLING',
            'RoutingStrategy_RROBIN', 'RoutingStrategy_WRROBIN',
            'RoutingStrategy_JSQ', 'RoutingStrategy_SQ',
            # the JMVA document has no capacity element at all
            'FiniteCapacity',
        }

    get_jmva_feature_set = getJMVAFeatureSet

    def getMethodFeatureSet(self, method):
        """SolverJMT drives TWO ENGINES, and they accept different models.

        'default', 'jsim' and 'replication' run the JSIM SIMULATOR, whose
        envelope is ``getFeatureSet``. The 'jmva.*' names run the JMVA
        ANALYTICAL engine, which reads a document carrying only a station type,
        a per-chain demand, a per-chain visit count, the populations or arrival
        rates and a reference station -- so declaring the JSIM envelope for
        jmva was a promise the writer could not keep.

        Defining this is also what lets the base runAnalyzerChecks gate name the
        offending features (mirrors MATLAB SolverJMT.getMethodFeatureSet):
        without it the coarse supports(model) is used, which accepts every
        model. A non-Network model (e.g. a LayeredNetwork) keeps the coarse
        path and any structural checks that operate on such models."""
        from ....lang.network import Network
        model = getattr(self, 'model', None)
        if not isinstance(model, Network):
            return None
        if str(method or '').lower().startswith('jmva'):
            feats = SolverJMT.getJMVAFeatureSet()
            if jmva_is_closed_only(method):
                # RECAL, CoMoM, Chow, Bard-Schweitzer, AQL, Linearizer and De
                # Souza-Muntz Linearizer are closed-network algorithms: JMT
                # answers an open or a mixed model with "The selected solver
                # cannot handle open classes" and a load-dependent one with the
                # matching refusal. Exact MVA, which 'jmva' and 'jmva.mva'
                # select, serves both.
                # the eight closed-form algorithms are single-server ones
                # (jmt_method_refusal words it); exact MVA carries the count
                feats = feats - {'OpenClass', 'LoadDependence', 'MultiServer'}
            return feats
        return SolverJMT.getFeatureSet()

    get_method_feature_set = getMethodFeatureSet

    def supportsModelMethod(self, method):
        """Structural gate for what no registry name can state.

        Three rules: the finite timespan the 'replication' arm integrates over,
        the single-server restriction of the closed-form JMVA algorithms (a
        server count is not a declared feature), and the one feature JMT admits
        in a RE-ENCODED form only. Limited load dependence has no
        representation of its own in either JMT document: the JSIM writer turns
        it into a server count and the JMVA writer into the matching
        <ldstation>, so alpha(n) = min(n,c) with an integer c is written
        exactly and any other scaling would be solved at a service rate JMT
        never saw. The first two come from ``jmt_method_refusal``, which the
        analyzer asks as well."""
        from ....constants import GlobalConstants
        from ....lang.network import Network
        model = getattr(self, 'model', None)
        if isinstance(model, Network):
            # The same predicate the JMVA arm and the replication arm ask, so
            # this gate and those runs cannot answer differently.
            structural = jmt_method_refusal(model.getStruct(), method, self.options)
            if structural:
                return False, structural
        lld = model.getStruct().lldscaling if isinstance(model, Network) else None
        if lld is not None:
            lld = np.atleast_2d(np.asarray(lld, dtype=float))
            for ist in range(lld.shape[0] if lld.ndim == 2 else 0):
                alpha = lld[ist, :]
                if alpha.size == 0 or np.all(alpha == 1.0):
                    continue
                c = float(np.max(alpha))
                shape = np.minimum(np.arange(1, alpha.size + 1, dtype=float), c)
                if c != round(c) or c < 1 or np.any(np.abs(alpha - shape) > GlobalConstants.Zero):
                    return False, (
                        'Station %d uses a load-dependent scaling that is not the multiserver '
                        'encoding alpha(n) = min(n,c): JMT has no representation for it, since '
                        'both the JSIM and the JMVA writer carry the scaling as a server count, '
                        'and the model would be solved at the nominal service rate. Use '
                        'SolverCTMC, SolverNC, SolverMVA or SolverSSA, which read sn.lldscaling '
                        'directly.' % (ist + 1))
        return super().supportsModelMethod(method)

    supports_model_method = supportsModelMethod

    @staticmethod
    def supports(model) -> bool:
        """Check if this solver supports the given model.

        Mirrors MATLAB SolverJMT.supports. This previously returned True
        unconditionally, so it accepted models built on features the JSIM writer
        cannot represent (e.g. FCFSPR).
        """
        from ...base import supports_via_featureset
        return supports_via_featureset(SolverJMT, model)

    @staticmethod
    def defaultOptions() -> OptionsDict:
        """Get default solver options."""
        return OptionsDict({
            'method': 'jsim',
            'samples': 10000,
            'seed': 23000,
            'verbose': default_verbose(),
            'keep': False,
            'conf_int': 0.99,
            'max_rel_err': 0.03,
        })

    # =========================================================================
    # File Management Methods (Gap 3a)
    # =========================================================================

    def _ensureTempDir(self, solvername: str) -> str:
        """Allocate, once, a private temporary directory for this solver instance.

        Mirrors MATLAB getJSIMTempPath/getJMVATempPath, which set self.filePath =
        lineTempName(solvername) on first use and reuse it thereafter. The
        allocation must be cached: callers write a model with writeJSIM() and then
        read the path back, so handing out a fresh directory per call would point
        them at an empty one. The directory must also be private to the instance,
        because the previous fixed fallback (/tmp/line_jmt/model.jsimg) made two
        concurrent writeJSIM() calls collide on one file.

        Whichever accessor runs first fixes the directory for both, exactly as in
        MATLAB, where the second call finds self.filePath already set.
        """
        if not getattr(self, '_temp_dir', None):
            base = os.path.join(tempfile.gettempdir(), 'line_workspace', solvername)
            os.makedirs(base, exist_ok=True)
            self._temp_dir = tempfile.mkdtemp(prefix='tmp_', dir=base)
        return self._temp_dir

    def getFileName(self) -> str:
        """Get the model file name, WITHOUT directory and WITHOUT extension.

        Matches MATLAB getFileName.m, whose callers build the name as
        [fileName '.jsim'], and the JAR, which does fileName + ".jsim". This
        previously returned 'model.jsimg', i.e. it included the extension, so it
        did not compose the way the other two codebases' callers expect.
        """
        return 'model'

    def getFilePath(self) -> str:
        """Get the directory holding the model file.

        Returns the DIRECTORY, matching MATLAB getFilePath (out = self.filePath);
        the file name is getFileName() and the joined path is getJSIMTempPath().
        """
        return self._ensureTempDir('jsim')

    @staticmethod
    def getJMTJarPath() -> str:
        """Get path to JMT.jar."""
        from ....api.solvers.jmt.handler import _get_jmt_jar_path
        return _get_jmt_jar_path()

    def getJMVATempPath(self) -> str:
        """Get path to the temporary JMVA model file."""
        return os.path.join(self._ensureTempDir('jmva'), 'model.jmva')

    def getJSIMTempPath(self) -> str:
        """Get path to the temporary JSIM model file."""
        return os.path.join(self._ensureTempDir('jsim'), 'model.jsim')

    # =========================================================================
    # Export Methods (Gap 3b)
    # =========================================================================

    def writeJMVA(self, outputFileName: str = None) -> str:
        """Write model to JMVA format.

        Args:
            outputFileName: Output file path. If None, writes to temp directory.

        Returns:
            Path to the written file.
        """
        from ..solver_qns.jmva_writer import write_jmva
        sn = self._sn
        if sn is None:
            self._extract_network_params()
            sn = self._sn

        if outputFileName is None:
            outputFileName = self.getJMVATempPath()

        os.makedirs(os.path.dirname(outputFileName), exist_ok=True)
        write_jmva(sn, outputFileName, {
            'method': self.method,
            'samples': self.options.samples,
        })
        return outputFileName

    def writeJSIM(self, outputFileName: str = None) -> str:
        """Write model to JSIM XML format.

        Args:
            outputFileName: Output file path. If None, writes to temp directory.

        Returns:
            Path to the written file.
        """
        # handler exports no underscored SolverJMTOptions name; import the real one.
        from ....api.solvers.jmt.handler import _write_jsim_file
        sn = self._sn
        if sn is None:
            self._extract_network_params()
            sn = self._sn

        if outputFileName is None:
            outputFileName = self.getJSIMTempPath()

        os.makedirs(os.path.dirname(outputFileName), exist_ok=True)
        handler_options = _SolverJMTOptions(
            method=self.method,
            samples=self.options.samples,
            seed=self.options.seed,
            max_simulated_time=self.options.max_simulated_time,
            conf_int=self.options.conf_int,
            max_rel_err=self.options.max_rel_err,
        )
        _write_jsim_file(sn, outputFileName, handler_options, model=self.model)
        return outputFileName

    def QN2JSIMG(self, outputFileName: str = None) -> str:
        """Convert queueing network to JSIMG format. Wrapper for writeJSIM.

        Args:
            outputFileName: Output file path.

        Returns:
            Path to the written file.
        """
        return self.writeJSIM(outputFileName)

    # =========================================================================
    # Transient Methods (Gap 3c)
    # =========================================================================

    def getTranCdfRespT(self, R=None):
        """Get transient CDF of response times.

        The same logged pipeline as getCdfRespT WITHOUT the steady-state seed:
        the reference @SolverJMT/getTranCdfRespT.m starts the logged run from
        the model's default initial state, so the collected samples cover the
        transient, where getCdfRespT preloads the rounded steady-state queue
        lengths to shorten the warmup.
        """
        return self._cdfRespTPipeline(R, init_from_steady=False)

    def getTranCdfPassT(self, R=None):
        """Get transient CDF of passage times. Delegates to getTranCdfRespT,
        its own name in the reference's sibling file."""
        return self.getTranCdfRespT(R)

    def getTranProbAggr(self, node=None):
        """Get transient aggregated state probabilities from simulation.

        Runs simulation and computes time-windowed probability from trajectory.

        Args:
            node: Node index or node object. If None, returns for all nodes.

        Returns:
            Dict with 't' (time vector) and 'prob' (probability trajectory).
        """
        # Run simulation to get trajectory
        result = self.sampleSysAggr(num_events=self.options.samples)
        if result is None:
            return None

        if node is not None:
            node_idx = node if isinstance(node, int) else getattr(node, '_station_index', 0)
            if 'states' in result and node_idx < len(result['states']):
                return {
                    't': result.get('t', np.array([])),
                    'state': result['states'][node_idx],
                }
        return result

    def _runReplication(self):
        """Transient averages by independent replication.

        Samples ``options.iter_max`` seeded system trajectories, interpolates each
        onto the union of their time grids (previous-neighbour, capped at the
        MINIMUM of their maxima so the state predictor never runs past the
        constraints the shortest replication established) and averages them.

        Utilization is read as ``min(n, c)/c`` at a finite server and as the raw
        queue length at a delay; throughput follows it as ``U*c*mu`` and ``U*mu``.
        Mirrors the 'replication' arm of MATLAB @SolverJMT/runAnalyzer.m.
        """
        import time as _time
        sn = self._sn
        M, K = sn.nstations, sn.nclasses
        # The predicate supportsModelMethod asks, so the gate that decides
        # whether to OFFER 'replication' and this run cannot drift apart.
        structural = jmt_method_refusal(sn, 'replication', self.options)
        if structural:
            raise RuntimeError(structural)
        t0 = _time.time()
        init_seed = self.options.seed
        reps = max(1, int(getattr(self.options, 'iter_max', 10) or 10))

        paths = []
        tumax = float('inf')
        grid = set()
        for it in range(reps):
            self.options.seed = init_seed + it
            try:
                path = self.sampleSysAggr()
            except Exception as exc:
                line_warning("SolverJMT", "Replication %d failed (%s), skipping.", it + 1, exc)
                continue
            if not path or path.get('t') is None or len(path['t']) == 0:
                line_warning("SolverJMT",
                             "Replication %d produced empty/invalid time series, skipping.", it + 1)
                continue
            tv = np.asarray(path['t'], dtype=float).ravel()
            paths.append((tv, path['state']))
            tumax = min(tumax, float(np.max(tv)))
            grid.update(tv.tolist())
        self.options.seed = init_seed
        if not paths:
            raise RuntimeError("No valid replications produced. Cannot compute transient averages.")

        tu = np.array(sorted(v for v in grid if v <= tumax), dtype=float)
        nvalid = len(paths)

        QNt, UNt, TNt = {}, {}, {}
        nservers = np.asarray(sn.nservers, dtype=float).ravel()
        for i in range(M):
            c = nservers[i] if i < len(nservers) else 1.0
            for k in range(K):
                q = np.zeros(len(tu))
                u = np.zeros(len(tu))
                for tv, states in paths:
                    st = states[i] if i < len(states) else None
                    if st is None:
                        continue
                    st = np.asarray(st, dtype=float)
                    if st.ndim != 2 or k >= st.shape[1] or not np.isfinite(st).any():
                        continue
                    col = st[:, k]
                    # previous-neighbour interpolation; a grid point before the
                    # first sample has no predecessor and reads 0, which is what
                    # the reference's own NaN-to-zero step leaves
                    idx = np.searchsorted(tv, tu, side='right') - 1
                    valid = idx >= 0
                    qv = np.zeros(len(tu))
                    qv[valid] = np.nan_to_num(col[idx[valid]])
                    q += qv / nvalid
                    if np.isfinite(c) and c > 0:
                        uv = np.zeros(len(tu))
                        uv[valid] = np.nan_to_num(np.minimum(col[idx[valid]], c) / c)
                    else:
                        uv = qv
                    u += uv / nvalid
                rate = float(sn.rates[i, k]) if sn.rates is not None else 0.0
                if not np.isfinite(rate):
                    rate = 0.0
                scale = (c * rate) if np.isfinite(c) else rate
                QNt[(i, k)] = {'t': tu.copy(), 'metric': q}
                UNt[(i, k)] = {'t': tu.copy(), 'metric': u}
                TNt[(i, k)] = {'t': tu.copy(), 'metric': u * scale}

        self._tran_avg = (QNt, UNt, TNt)
        self._tran_runtime = _time.time() - t0
        return self

    def getTranAvg(self):
        """Get transient average metrics from simulation.

        Runs simulation with logging and extracts time series.

        Returns:
            Tuple of (QNt, UNt, TNt) time series dicts, or None if unavailable.
        """
        # method='replication' produced the real transient mean; return it rather
        # than the constant series the steady-state fallback below builds.
        tran = getattr(self, '_tran_avg', None)
        if tran is not None:
            return tran
        # Transient analysis runs the simulation with a finite timespan.
        if self._result is None:
            self._ensureAvgResults()

        # JMT steady-state results don't have time series natively.
        # Return steady-state values as constant time series.
        if self._result is None:
            return None, None, None

        M = self._sn.nstations
        K = self._sn.nclasses
        t = np.array([0.0, self.options.max_simulated_time if np.isfinite(self.options.max_simulated_time) else 1.0])

        QNt = {}
        UNt = {}
        TNt = {}
        for i in range(M):
            for k in range(K):
                q_val = self._result.Q[i, k] if self._result.Q is not None else 0.0
                u_val = self._result.U[i, k] if self._result.U is not None else 0.0
                t_val = self._result.T[i, k] if self._result.T is not None else 0.0
                QNt[(i, k)] = {'t': t.copy(), 'metric': np.array([q_val, q_val])}
                UNt[(i, k)] = {'t': t.copy(), 'metric': np.array([u_val, u_val])}
                TNt[(i, k)] = {'t': t.copy(), 'metric': np.array([t_val, t_val])}

        return QNt, UNt, TNt

    # =========================================================================
    # Probability Methods (Gap 3d)
    # =========================================================================

    def getProb(self, node=None, state=None):
        """Get state probability from simulation trajectory.

        Args:
            node: Node index or node object.
            state: Target state vector. If None, returns probability of current state.

        Returns:
            Float probability value, or dict of probabilities.
        """
        if self._result is None:
            self._ensureAvgResults()

        # For simulation-based solver, compute from trajectory
        result = self.sampleSysAggr(num_events=self.options.samples)
        if result is None:
            return 0.0

        if node is None:
            return self.getProbSys()

        node_idx = node if isinstance(node, int) else getattr(node, '_station_index', 0)

        if state is None:
            # Use current state from sn
            if self._sn is not None and hasattr(self._sn, 'state') and self._sn.state is not None:
                state = self._sn.state[node_idx] if node_idx < len(self._sn.state) else None

        if state is None or result is None:
            return 0.0

        # Time-weighted probability
        if 'states' in result and node_idx < len(result.get('states', [])):
            trajectory = result['states'][node_idx]
            t = result.get('t', np.array([]))
            if len(t) < 2:
                return 0.0
            target = np.asarray(state).flatten()
            total_time = 0.0
            match_time = 0.0
            for idx in range(len(t) - 1):
                dt = t[idx + 1] - t[idx]
                total_time += dt
                if np.allclose(trajectory[idx], target, atol=1e-10):
                    match_time += dt
            return match_time / total_time if total_time > 0 else 0.0

        return 0.0

    def getProbSys(self):
        """Get joint system state probability.

        Returns:
            Float probability of the current system state.
        """
        return self.getProbSysAggr()

    def getProbMarg(self, node=None, jobclass=None):
        """Get marginal state probability for a specific class at a node.

        Args:
            node: Node index or object.
            jobclass: Class index.

        Returns:
            Dict mapping state values to probabilities.
        """
        result = self.sampleSysAggr(num_events=self.options.samples)
        if result is None:
            return {}

        node_idx = node if isinstance(node, int) else getattr(node, '_station_index', 0)
        class_idx = jobclass if isinstance(jobclass, int) else getattr(jobclass, '_index', 0)

        if 'states' not in result or node_idx >= len(result.get('states', [])):
            return {}

        trajectory = result['states'][node_idx]
        t = result.get('t', np.array([]))
        if len(t) < 2:
            return {}

        # Compute time-weighted histogram for the given class
        prob_map = {}
        total_time = 0.0
        for idx in range(len(t) - 1):
            dt = t[idx + 1] - t[idx]
            total_time += dt
            state_val = int(trajectory[idx][class_idx]) if class_idx < len(trajectory[idx]) else 0
            prob_map[state_val] = prob_map.get(state_val, 0.0) + dt

        if total_time > 0:
            for k in prob_map:
                prob_map[k] /= total_time

        return prob_map

    def getProbNormConstAggr(self):
        """Log normalizing constant, from the JMVA engine only.

        A simulation computes no normalizing constant, but the analytical JMVA
        algorithms report one in the result file's <normconst logValue>, which
        MATLAB stores as result.Prob.logNormConstAggr. It is returned here for
        the jmva* methods and refused for the simulation ones rather than
        handing back the NaN placeholder.

        Raises:
            NotImplementedError: on the simulation methods.
        """
        from ....api.solvers.jmt.handler import _is_jmva_method
        if not _is_jmva_method(self.method):
            raise NotImplementedError(
                "getProbNormConstAggr() is not supported by SolverJMT with method='%s'. "
                "Use an analytical method (SolverJMT 'jmva'), SolverNC or SolverCTMC "
                "for normalizing constant computation." % self.method)
        if self._result is None:
            self._ensureAvgResults()
        return getattr(self._result, 'logNormConstAggr', float('nan'))

    # =========================================================================
    # Sampling Methods (Gap 3e)
    # =========================================================================

    def sample(self, node, numEvents: int = 1000):
        """Sample the aggregated state trajectory at a node (JMT logs only carry
        per-class counts, so this is equivalent to :meth:`sampleAggr`).

        Args:
            node: Node index or node object.
            numEvents: Number of events to sample.

        Returns:
            Dict with 't' (time vector) and 'state' (per-class counts).
        """
        return self.sampleAggr(node, numEvents)

    def _resolve_node(self, node):
        """Resolve a node argument to (node_index_0based, node_name)."""
        sn = self._sn if self._sn is not None else self.model.get_struct()
        if isinstance(node, int):
            node_idx = node
            node_name = sn.nodenames[node_idx]
        else:
            node_name = node.name if hasattr(node, 'name') else str(node)
            node_idx = self.model.get_node_index(node_name) - 1  # 0-based
        return node_idx, node_name

    def _cleanup_log_dir(self, path):
        """Remove a temporary JMT sampling log directory unless options.keep is set.

        Mirrors the JAR SolverJMT.cleanupDir and MATLAB @SolverJMT/runAnalyzer.m
        (keep=false), so log-based sampling does not accumulate scratch folders.
        """
        if getattr(self.options, 'keep', False):
            return
        if not path:
            return
        shutil.rmtree(path, ignore_errors=True)

    def _run_logged_copy(self, is_node_logged, numEvents):
        """Simulate a logged copy of the model and return (modelCopy, log_path).

        Mirrors MATLAB SolverJMT.sampleAggr: build a temporary copy, apply
        linkAndLog to the requested nodes, and run JMT so the per-node
        arrival/departure CSV logs are produced.
        """
        import tempfile
        sn = self._sn if self._sn is not None else self.model.get_struct()

        model_copy = self.model.copy()
        model_copy.reset_network()

        Plinked = None
        if hasattr(self.model, 'get_linked_routing_matrix'):
            Plinked = self.model.get_linked_routing_matrix()
        if Plinked is None and hasattr(sn, 'rtorig'):
            Plinked = sn.rtorig
        if Plinked is None:
            raise RuntimeError("JMT log-based sampling requires the routing "
                               "matrix (rtorig).")

        log_path = tempfile.mkdtemp(prefix='jmt_sample_logs_')
        model_copy.link_and_log(Plinked, is_node_logged, log_path)

        jmt = SolverJMT(model_copy, self.options)
        # JMT cannot cap events per node; scale the total event budget as MATLAB
        # does (numEvents * nnodes * nclasses) so each node yields ~numEvents.
        if numEvents and numEvents > 0:
            try:
                jmt.maxEvents = int(numEvents) * sn.nnodes * sn.nclasses
            except Exception:
                pass
        jmt.runAnalyzer()
        return model_copy, log_path

    def _node_preload(self, node_idx, K):
        """Initial per-class population at a node (INIT row of the trajectory)."""
        from ....api.sn.transforms import sn_get_state_aggr
        sn = self._sn if self._sn is not None else self.model.get_struct()
        preload = np.zeros(K)
        try:
            preload_map = sn_get_state_aggr(sn)
            isf = int(sn.nodeToStateful[node_idx]) if sn.nodeToStateful is not None else -1
            if isf in preload_map and preload_map[isf] is not None:
                pl = np.ravel(np.asarray(preload_map[isf], dtype=float))
                preload[:min(K, len(pl))] = pl[:min(K, len(pl))]
        except Exception:
            pass
        return preload

    def sampleAggr(self, node, numEvents: int = 1000):
        """Sample the aggregated (per-class count) state trajectory at a node.

        Faithful port of MATLAB SolverJMT.sampleAggr: a temporary logged copy of
        the model is simulated and its arrival/departure logs are reconstructed
        into a piecewise-constant per-class queue-length trajectory.

        Args:
            node: Node index (0-based) or node object.
            numEvents: Desired number of sampled events at the node.

        Returns:
            Dict with keys 'handle', 't' (event-boundary times), 'state'
            (len x nclasses per-class counts), 'event' (chronological event
            list) and 'isaggregate'=True.
        """
        import os
        from ....api.solvers.jmt.handler import parse_tran_state

        sn = self._sn if self._sn is not None else self.model.get_struct()
        K = sn.nclasses
        node_idx, node_name = self._resolve_node(node)

        nnodes = self.model.get_number_of_nodes()
        is_node_logged = [False] * nnodes
        is_node_logged[node_idx] = True

        model_copy, log_path = self._run_logged_copy(is_node_logged, numEvents)

        preload = self._node_preload(node_idx, K)
        class_names = [sn.classnames[r] for r in range(K)] if sn.classnames is not None else None

        arv_file = os.path.join(log_path, f"{node_name}-Arv.csv")
        dep_file = os.path.join(log_path, f"{node_name}-Dep.csv")
        state, evtype, evclass, evjob = parse_tran_state(
            arv_file, dep_file, preload, class_names)

        # see _kb/06-solver-catalog.md (Wrappers: "JMT python wrapper:
        # export/import internals") for the event-boundary time shift
        t_all = state[:, 0]
        _, uniq_idx = np.unique(t_all, return_index=True)
        uniq_idx = np.sort(uniq_idx)
        t = t_all[uniq_idx]
        qlen = state[uniq_idx, 1:1 + K]

        if len(t) > 1:
            t_shift = np.concatenate([t[1:], t[-1:]])
        else:
            t_shift = t.copy()

        m = min(len(t_shift), 1 + numEvents) if (numEvents and numEvents > 0) else len(t_shift)
        t_shift = t_shift[:m]
        qlen = qlen[:m, :]
        t_out = np.concatenate([[0.0], t_shift[:-1]]) if len(t_shift) > 0 else t_shift

        tmax = float(t_out[-1]) if len(t_out) > 0 else np.inf
        events = []
        for e in range(len(evtype)):
            te = float(state[e, 0])
            if te > tmax:
                continue
            events.append({
                'event': evtype[e],
                'node': node_idx,
                'class': (int(evclass[e]) if not np.isnan(evclass[e]) else None),
                't': te,
                'job': (int(evjob[e]) if not np.isnan(evjob[e]) else None),
            })

        self._cleanup_log_dir(log_path)
        return {
            'handle': node,
            't': t_out,
            'state': qlen,
            'event': events,
            'isaggregate': True,
        }

    def sampleSys(self, numEvents: int = 1000):
        """Sample system-wide state trajectory.

        Args:
            numEvents: Number of events to sample.

        Returns:
            Dict with 't' (time vector) and 'states' (list of per-node state matrices).
        """
        return self.sampleSysAggr(num_events=numEvents)

    # =========================================================================
    # Additional Metric Methods
    # =========================================================================

    def getAvgResidT(self) -> np.ndarray:
        """Get average residence times (M x K)."""
        if self._result is None:
            self._ensureAvgResults()
        if self._sn is not None and hasattr(self._sn, 'visits') and self._sn.visits:
            return sn_get_residt_from_respt(self._sn, self._result.R, None)
        return self._result.R.copy() if self._result.R is not None else np.array([])

    def getAvgWaitT(self) -> np.ndarray:
        """Get average waiting times (M x K). W = R - S."""
        if self._result is None:
            self._ensureAvgResults()
        R = self._result.R.copy() if self._result.R is not None else np.array([])
        if len(R) == 0:
            return R
        if hasattr(self._sn, 'rates') and self._sn.rates is not None:
            rates = np.asarray(self._sn.rates)
            S = np.zeros_like(rates)
            nonzero = rates > 0
            S[nonzero] = 1.0 / rates[nonzero]
            W = R - S
            W = np.maximum(W, 0.0)
            return W
        return R

    def getAvg(self):
        """Get all average metrics at once.

        Returns:
            Tuple of (Q, U, R, T, A, W)
        """
        if self._result is None:
            self._ensureAvgResults()
        r = self._result
        Q = r.Q if r.Q is not None else np.array([])
        U = r.U if r.U is not None else np.array([])
        R = r.R if r.R is not None else np.array([])
        T = r.T if r.T is not None else np.array([])
        A = r.A if r.A is not None else T.copy()
        W = self.getAvgResidT()
        return Q, U, R, T, A, W

    def getAvgSys(self):
        """Get system-level average metrics.

        Returns:
            Tuple of (CN, XN) - system response times and throughputs
        """
        return self.getAvgSysRespT(), self.getAvgSysTput()

    def getAvgNode(self):
        """Get average metrics per node.

        Returns:
            Tuple of (QNn, UNn, RNn, WNn, ANn, TNn) - node-level metrics
        """
        if self._result is None:
            self._ensureAvgResults()
        sn = self._sn
        I = sn.nnodes
        M = sn.nstations
        K = sn.nclasses
        QN = self._result.Q if self._result.Q is not None else np.zeros((M, K))
        UN = self._result.U if self._result.U is not None else np.zeros((M, K))
        RN = self._result.R if self._result.R is not None else np.zeros((M, K))
        TN = self._result.T if self._result.T is not None else np.zeros((M, K))
        AN = self._result.A if self._result.A is not None else TN.copy()
        WN = self.getAvgResidT() if len(RN) > 0 else np.zeros((M, K))

        QNn = np.zeros((I, K))
        UNn = np.zeros((I, K))
        RNn = np.zeros((I, K))
        WNn = np.zeros((I, K))
        TNn = np.zeros((I, K))
        ANn = np.zeros((I, K))
        for ist in range(M):
            ind = int(sn.stationToNode[ist])
            if 0 <= ind < I:
                QNn[ind, :] = QN[ist, :]
                UNn[ind, :] = UN[ist, :]
                RNn[ind, :] = RN[ist, :]
                WNn[ind, :] = WN[ist, :]
                TNn[ind, :] = TN[ist, :]
                ANn[ind, :] = AN[ist, :]
        return QNn, UNn, RNn, WNn, ANn, TNn

    def getCdfRespT(self, R=None):
        """Get response time CDF via transient simulation with logging.

        This method runs JMT twice:
        1. First run: Get steady-state queue lengths to initialize state
        2. Second run: Run with logging enabled to collect response time samples

        Ported from MATLAB's SolverJMT.getCdfRespT.

        Args:
            R: Optional response time handles (uses defaults if None)

        Returns:
            List of lists where RD[station][class] is a 2D array [cdf, time]
        """
        return self._cdfRespTPipeline(R, init_from_steady=True)

    def _cdfRespTPipeline(self, R=None, init_from_steady=True):
        """The shared logged-run pipeline behind getCdfRespT (seeded from the
        rounded steady-state queue lengths) and getTranCdfRespT (unseeded)."""
        import os
        import tempfile
        from ....api.solvers.jmt.handler import parse_tran_resp_t

        sn = self._sn
        M = sn.nstations
        K = sn.nclasses

        # Initialize result structure
        RD = [[None for _ in range(K)] for _ in range(M)]

        n = None
        if init_from_steady:
            # Step 1: Get steady-state queue lengths (first JMT run)
            QN = self.getAvgQLen()
            n = QN.copy()

            # Adjust job numbers based on network constraints
            for r in range(K):
                if np.isinf(sn.njobs[r]):
                    # Open class - use floor of queue lengths
                    for i in range(M):
                        n[i, r] = np.floor(QN[i, r])
                else:
                    # Closed class - ensure total population equals njobs
                    for i in range(M):
                        n[i, r] = np.floor(QN[i, r])
                    total_jobs = np.sum(n[:, r])
                    if total_jobs < sn.njobs[r]:
                        # Put remaining jobs on bottleneck station
                        imax = np.argmax(n[:, r])
                        n[imax, r] = n[imax, r] + sn.njobs[r] - total_jobs

        # Step 2: Copy model for CDF computation
        cdfmodel = self.model.copy()
        cdfmodel.reset_network()
        cdfmodel.reset()

        # Determine which nodes should be logged (all stations except Source/Sink)
        nnodes = cdfmodel.get_number_of_nodes()
        is_node_logged = [False] * nnodes

        for i in range(cdfmodel.get_number_of_stations()):
            station = cdfmodel.get_stations()[i]
            node_idx = cdfmodel.get_node_index(station.name) - 1  # Convert to 0-based
            # Don't log Source or Sink
            if hasattr(station, 'node_type'):
                from ....lang.base import NodeType
                if station.node_type not in (NodeType.SOURCE, NodeType.SINK):
                    is_node_logged[node_idx] = True
            else:
                is_node_logged[node_idx] = True

        # Get original routing matrix from the model
        Plinked = None
        if hasattr(self.model, 'get_linked_routing_matrix'):
            Plinked = self.model.get_linked_routing_matrix()
        if Plinked is None and hasattr(sn, 'rtorig'):
            Plinked = sn.rtorig

        if Plinked is None:
            raise RuntimeError("getCdfRespT requires routing matrix (rtorig)")

        # Step 3: Set up logging
        log_path = tempfile.mkdtemp(prefix='jmt_cdf_logs_')
        cdfmodel.link_and_log(Plinked, is_node_logged, log_path)

        # Initialize model state from marginal distribution (seeded route only;
        # the transient getter starts from the model's default initial state)
        if init_from_steady and n is not None:
            try:
                cdfmodel.init_from_marginal(n)
            except Exception:
                pass  # May not be supported for all models

        # Step 4: Run JMT on logged model (second JMT run)
        cdf_solver = SolverJMT(cdfmodel, self.options)
        cdf_solver.runAnalyzer()

        # Step 5: Parse logs to get response time samples
        node_names = self.model.get_node_names() if hasattr(self.model, 'get_node_names') else []
        station_names = sn.nodenames if sn.nodenames else []

        # Get class names in model order for deterministic mapping
        class_names = [sn.classnames[r] for r in range(K)] if sn.classnames is not None else None

        for i in range(M):
            # Get original station name
            node_idx = int(sn.stationToNode[i]) if sn.stationToNode is not None else i
            station_name = station_names[node_idx] if node_idx < len(station_names) else f'Station{i}'

            # Check if this node was logged
            if node_idx < len(is_node_logged) and is_node_logged[node_idx]:
                arv_file = os.path.join(log_path, f"{station_name}-Arv.csv")
                dep_file = os.path.join(log_path, f"{station_name}-Dep.csv")

                if os.path.exists(arv_file) and os.path.exists(dep_file):
                    class_resp_t, _, _ = parse_tran_resp_t(arv_file, dep_file, class_names=class_names)

                    for r in range(min(K, len(class_resp_t))):
                        resp_times = class_resp_t[r]
                        if len(resp_times) > 0:
                            # Create empirical CDF (ecdf equivalent, matching MATLAB ecdf)
                            sorted_times = np.sort(resp_times)
                            n = len(sorted_times)
                            # Group duplicate values (matches MATLAB ecdf and JAR createEmpiricalCDF)
                            unique_vals, counts = np.unique(sorted_times, return_counts=True)
                            cumulative = np.cumsum(counts) / n
                            # Add point at F=0 for the first data point
                            X = np.concatenate([[unique_vals[0]], unique_vals])
                            F = np.concatenate([[0.0], cumulative])
                            # Store as [F, X] format
                            RD[i][r] = np.column_stack([F, X])

        self._cleanup_log_dir(log_path)
        return RD

    def getPerctRespT(self, percentiles=None):
        """Get response time percentiles.

        Args:
            percentiles: Array of percentile values (default: [90, 95, 99])

        Returns:
            Tuple of (PercRT, PercTable) where PercRT is list of dicts
            and PercTable is a pandas DataFrame
        """
        import pandas as pd

        if percentiles is None:
            percentiles = np.array([90, 95, 99])
        else:
            percentiles = np.asarray(percentiles)

        if self._result is None:
            self._ensureAvgResults()

        R = self._result.R
        M = self._sn.nstations
        K = self._sn.nclasses

        PercRT = []
        rows = []
        perc_col_names = [f'P{int(p)}' for p in percentiles]
        percentiles_normalized = percentiles / 100.0

        station_names = self._sn.nodenames if self._sn.nodenames else [f'Station{i}' for i in range(M)]
        class_names = self._sn.classnames if self._sn.classnames else [f'Class{r}' for r in range(K)]

        for i in range(M):
            for r in range(K):
                if R is not None and i < R.shape[0] and r < R.shape[1]:
                    mean_resp_t = R[i, r]
                    if mean_resp_t > 0 and not np.isnan(mean_resp_t):
                        lambda_rate = 1.0 / mean_resp_t
                        perc_values = -np.log(1 - percentiles_normalized) / lambda_rate

                        PercRT.append({
                            'station': i + 1,
                            'class': r + 1,
                            'percentiles': percentiles.tolist(),
                            'values': perc_values.tolist(),
                        })

                        node_idx = int(self._sn.stationToNode[i]) if self._sn.stationToNode is not None else i
                        station_name = station_names[node_idx] if node_idx < len(station_names) else f'Station{i}'

                        row_data = {
                            'Station': station_name,
                            'Class': class_names[r] if r < len(class_names) else f'Class{r}',
                        }
                        for perc_col, perc_val in zip(perc_col_names, perc_values):
                            row_data[perc_col] = perc_val
                        rows.append(row_data)

        PercTable = pd.DataFrame(rows) if rows else pd.DataFrame()
        return PercRT, PercTable

    def __repr__(self) -> str:
        return f"SolverJMT(method='{self.method}', samples={self.options.samples})"

    # Snake case aliases for MATLAB compatibility
    avg_qlen = getAvgQLen
    prob_sys_aggr = getProbSysAggr
    prob_aggr = getProbAggr
    avg_util = getAvgUtil
    avg_respt = getAvgRespT
    get_avg_respt = getAvgRespT
    avg_residt = getAvgResidT
    avg_waitt = getAvgWaitT
    avg_tput = getAvgTput
    avg_arv_r = getAvgArvR
    avg_fcr = getAvgFcr
    avg_chain_table = getAvgChainTable
    avg_sys_table = getAvgSysTable
    avg_sys_resp_t = getAvgSysRespT
    avg_sys_tput = getAvgSysTput
    run_analyzer = runAnalyzer
    get_runtime = getRuntime
    get_method = getMethod
    list_valid_methods = listValidMethods
    is_available = isAvailable
    get_feature_set = getFeatureSet
    default_options = defaultOptions
    cdf_resp_t = getCdfRespT
    cdf_respt = getCdfRespT
    get_cdf_resp_t = getCdfRespT
    get_tran_cdf_respt = getTranCdfRespT
    get_tran_cdf_resp_t = getTranCdfRespT
    get_tran_cdf_pass_t = getTranCdfPassT
    perct_resp_t = getPerctRespT
    perct_respt = getPerctRespT

    def getAvgNodeTable(self) -> pd.DataFrame:
        """
        Per-node average performance metrics. Mirrors MATLAB and JAR output by
        including non-station nodes (Source, Sink, Router, ClassSwitch) alongside
        station rows. Per-node arrival rates are computed via
        sn_get_node_arvr_from_tput; passthrough nodes get Tput == ArvR; Sink Tput
        is 0; queue/delay rows reuse the station-level metrics.
        """
        if self._result is None:
            self._ensureAvgResults()

        from line_solver.api.sn import sn_get_node_arvr_from_tput
        from line_solver.lang.base import NodeType

        sn = self._sn
        I = sn.nnodes
        K = sn.nclasses
        nodenames = sn.nodenames if sn.nodenames else [f'Node{i}' for i in range(I)]
        classnames = sn.classnames if sn.classnames else [f'Class{r}' for r in range(K)]
        nodetype = list(sn.nodetype) if sn.nodetype is not None else []

        TN = self._result.T if self._result.T is not None else np.zeros((sn.nstations, K))
        # see _kb/06-solver-catalog.md (Wrappers: "JMT python wrapper: export/import internals")
        AN = self._result.A if self._result.A is not None else TN.copy()
        ANn = sn_get_node_arvr_from_tput(sn, TN, None, AN)

        # Per-node Tput: for stations, use station Tput; for passthrough nodes,
        # equate to per-node ArvR; Sinks get 0.
        TNn = np.zeros((I, K))
        stationToNode = np.asarray(sn.stationToNode).flatten() if sn.stationToNode is not None else np.array([])
        for ist in range(sn.nstations):
            ind = int(stationToNode[ist]) if ist < len(stationToNode) else ist
            if 0 <= ind < I:
                TNn[ind, :] = TN[ist, :]
        for ind in range(I):
            if ind >= len(nodetype):
                continue
            nt = nodetype[ind]
            nt_val = int(nt) if not hasattr(nt, 'value') else int(nt.value)
            sink_val = NodeType.SINK.value if hasattr(NodeType.SINK, 'value') else int(NodeType.SINK)
            source_val = NodeType.SOURCE.value if hasattr(NodeType.SOURCE, 'value') else int(NodeType.SOURCE)
            queue_val = NodeType.QUEUE.value if hasattr(NodeType.QUEUE, 'value') else int(NodeType.QUEUE)
            delay_val = NodeType.DELAY.value if hasattr(NodeType.DELAY, 'value') else int(NodeType.DELAY)
            if nt_val == sink_val:
                TNn[ind, :] = 0.0
            elif nt_val not in (source_val, queue_val, delay_val):
                # Router, ClassSwitch, Join, Cache, etc.: passthrough flow
                TNn[ind, :] = ANn[ind, :]

        # Build rows: stations first (with full Q/U/R/ResidT) then non-station nodes
        if sn.visits and self._result.R is not None:
            WN = sn_get_residt_from_respt(sn, self._result.R, None)
        else:
            WN = self._result.R.copy() if self._result.R is not None else np.zeros((sn.nstations, K))

        rows = []
        station_node_set = set()
        for ist in range(sn.nstations):
            ind = int(stationToNode[ist]) if ist < len(stationToNode) else ist
            station_node_set.add(ind)
            for r in range(K):
                qlen = self._result.Q[ist, r] if self._result.Q is not None else 0.0
                util = self._result.U[ist, r] if self._result.U is not None else 0.0
                respt = self._result.R[ist, r] if self._result.R is not None else 0.0
                residt = WN[ist, r] if ist < WN.shape[0] and r < WN.shape[1] else respt
                arvr = ANn[ind, r] if 0 <= ind < I else 0.0
                tput = TNn[ind, r] if 0 <= ind < I else 0.0
                if ind < len(nodetype):
                    nt_v = int(nodetype[ind]) if not hasattr(nodetype[ind], 'value') else int(nodetype[ind].value)
                    src_v = NodeType.SOURCE.value if hasattr(NodeType.SOURCE, 'value') else int(NodeType.SOURCE)
                    if nt_v == src_v:
                        qlen = util = respt = residt = 0.0
                        arvr = 0.0
                rows.append({
                    'Node': nodenames[ind] if ind < len(nodenames) else f'Node{ind}',
                    'JobClass': classnames[r],
                    'QLen': qlen, 'Util': util, 'RespT': respt, 'ResidT': residt,
                    'ArvR': arvr, 'Tput': tput,
                })

        # Non-station nodes: emit one row per node x class with passthrough flow
        for ind in range(I):
            if ind in station_node_set or ind >= len(nodetype):
                continue
            for r in range(K):
                arvr = ANn[ind, r]
                tput = TNn[ind, r]
                if arvr == 0 and tput == 0:
                    continue
                rows.append({
                    'Node': nodenames[ind] if ind < len(nodenames) else f'Node{ind}',
                    'JobClass': classnames[r],
                    'QLen': 0.0, 'Util': 0.0, 'RespT': 0.0, 'ResidT': 0.0,
                    'ArvR': arvr, 'Tput': tput,
                })

        # see _kb/06-solver-catalog.md (Wrappers: "JMT python wrapper: export/import internals")
        F = int(sn.nregions) if getattr(sn, 'nregions', 0) else 0
        Qfcr = getattr(self._result, 'Qfcr', None)
        if F > 0 and Qfcr is not None:
            Rfcr = getattr(self._result, 'Rfcr', None)
            Wfcr = getattr(self._result, 'Wfcr', None)
            Tfcr = getattr(self._result, 'Tfcr', None)
            region_names = []
            if self.model is not None and hasattr(self.model, 'get_regions'):
                region_names = [rg.get_name() if hasattr(rg, 'get_name') else f'FCR{f+1}'
                                for f, rg in enumerate(self.model.get_regions())]
            for f in range(F):
                rname = region_names[f] if f < len(region_names) else f'FCR{f+1}'
                for r in range(K):
                    qlen = Qfcr[f, r] if Qfcr is not None else 0.0
                    respt = Rfcr[f, r] if Rfcr is not None else 0.0
                    residt = Wfcr[f, r] if Wfcr is not None else respt
                    tput = Tfcr[f, r] if Tfcr is not None else 0.0
                    if np.nansum([qlen, respt, tput]) <= 0:
                        continue
                    rows.append({
                        'Node': rname,
                        'JobClass': classnames[r],
                        'QLen': qlen, 'Util': np.nan, 'RespT': respt, 'ResidT': residt,
                        'ArvR': np.nan, 'Tput': tput,
                    })

        df = pd.DataFrame(rows)
        if not getattr(self, '_table_silent', False):
            print(df.to_string(index=False))
        return df

    avg_node_table = getAvgNodeTable
    get_avg_node_table = getAvgNodeTable
    get_file_name = getFileName
    get_file_path = getFilePath
    get_jmt_jar_path = getJMTJarPath
    write_jmva = writeJMVA
    write_jsim = writeJSIM
