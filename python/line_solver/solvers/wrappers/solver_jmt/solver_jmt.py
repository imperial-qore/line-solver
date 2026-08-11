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
from ...base import NetworkSolver, method_label


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
        verbose = kwargs.get('verbose', False)
        keep = kwargs.get('keep', False)
        conf_int = kwargs.get('conf_int', kwargs.get('confint', 0.99))
        max_rel_err = kwargs.get('max_rel_err', 0.03)
        max_simulated_time = kwargs.get('max_simulated_time',
                                        kwargs.get('timespan', [0, float('inf')])[1]
                                        if isinstance(kwargs.get('timespan'), list) else float('inf'))

        self.options = SolverJMTOptions(
            method=self.method,
            samples=samples,
            seed=seed,
            max_simulated_time=max_simulated_time,
            conf_int=conf_int,
            max_rel_err=max_rel_err,
            verbose=verbose,
            keep=keep
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

        if getattr(self._sn, 'immfeed', None) is not None and np.any(self._sn.immfeed):
            line_warning("SolverJMT", "SolverJMT does not support immediate feedback (immfeed); no solution returned.")
            return self

        method = self.options.method
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
            keep=self.options.keep
        )

        # Call the handler (pass model for FCR region support)
        self._result = solver_jmt(self._sn, handler_options, self.model)

        # Print completion message (matches MATLAB verbose guard)
        if self.options.verbose:
            py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
            runtime = self._result.runtime if hasattr(self._result, 'runtime') else 0.0
            method = self._result.method if hasattr(self._result, 'method') else self.options.method
            print(f"JMT analysis [method: {method_label(self.options.method, method)}, lang: python, env: {py_version}] completed in {runtime:.6f}s.")

        return self

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get average performance metrics as a DataFrame.

        Returns:
            DataFrame with columns: Station, Class, QLen, Util, RespT, Tput, ArvR
        """
        if self._result is None:
            self.runAnalyzer()

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
            self.runAnalyzer()
        return self._result.Q if self._result.Q is not None else np.array([])

    def getAvgUtil(self) -> np.ndarray:
        """Get average utilizations (M x K matrix)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.U if self._result.U is not None else np.array([])

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times (M x K matrix)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.R if self._result.R is not None else np.array([])

    def getAvgTput(self) -> np.ndarray:
        """Get average throughputs (M x K matrix)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.T if self._result.T is not None else np.array([])

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates (M x K matrix)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.A if self._result.A is not None else np.array([])

    def getAvgChainTable(self) -> pd.DataFrame:
        """
        Get average performance metrics aggregated by chain.

        Returns:
            DataFrame with columns: Chain, QLen, Util, RespT, Tput
        """
        if self._result is None:
            self.runAnalyzer()

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

        return pd.DataFrame(rows)

    def getAvgSysTable(self) -> pd.DataFrame:
        """
        Get system-level average performance metrics.

        Returns:
            DataFrame with columns: Chain, SysRespT, SysTput
        """
        if self._result is None:
            self.runAnalyzer()

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
            self.runAnalyzer()
        # Sum response times across all stations for each class
        if self._result.R is not None:
            return np.nansum(self._result.R, axis=0, keepdims=True)
        return np.array([[]])

    def getAvgSysTput(self) -> np.ndarray:
        """Get system throughputs (1 x K)."""
        if self._result is None:
            self.runAnalyzer()
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

    def getProbAggr(self, station: int) -> np.ndarray:
        """Get aggregated state probabilities at station.

        Note: JMT simulation does not directly compute state probabilities.
        Returns empty array as placeholder.
        """
        if self._result is None:
            self.runAnalyzer()
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
        tokens = re.split(r'[./]', str(method).lower())
        if 'jmva' in tokens:
            return any(tok in ('ls', 'mci', 'imci', 'sampling') for tok in tokens)
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
            'SchedStrategy_SEPT', 'SchedStrategy_SRPT', 'SchedStrategy_SRPTPRIO', 'SchedStrategy_LEPT',
            'SchedStrategy_SJF', 'SchedStrategy_LJF', 'SchedStrategy_LPS',
            'SchedStrategy_POLLING', 'SchedStrategy_EXT',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND',
            'RoutingStrategy_RROBIN', 'RoutingStrategy_WRROBIN',
            'RoutingStrategy_JSQ',
            'RoutingStrategy_SQ',
            'ClosedClass', 'SelfLoopingClass', 'OpenClass',
            'Cache', 'CacheClassSwitcher',
            'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO',
            'ReplacementStrategy_SFIFO', 'ReplacementStrategy_LRU',
            'Region',
            # Exported as delayOffTime/setUpTime (_write_delayoff_strategy).
            'SetupDelayOff',
            # Exported as Impatience/Reneging and Impatience/Balking strategies.
            'Reneging', 'Balking',
        }

    def getMethodFeatureSet(self, method):
        """All JMT methods share the solver-level feature envelope.

        Defining this is what lets the base runAnalyzerChecks gate name the
        offending features (mirrors MATLAB SolverJMT.getMethodFeatureSet):
        without it the coarse supports(model) is used, which accepts every
        model. A non-Network model (e.g. a LayeredNetwork) keeps the coarse
        path and any structural checks that operate on such models."""
        from ....lang.network import Network
        model = getattr(self, 'model', None)
        if not isinstance(model, Network):
            return None
        return SolverJMT.getFeatureSet()

    get_method_feature_set = getMethodFeatureSet

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
        """Get transient CDF of response times. Alias for getCdfRespT."""
        return self.getCdfRespT(R)

    def getTranCdfPassT(self, R=None):
        """Get transient CDF of passage times. Alias for getCdfRespT."""
        return self.getCdfRespT(R)

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

    def getTranAvg(self):
        """Get transient average metrics from simulation.

        Runs simulation with logging and extracts time series.

        Returns:
            Tuple of (QNt, UNt, TNt) time series dicts, or None if unavailable.
        """
        # Transient analysis runs the simulation with a finite timespan.
        if self._result is None:
            self.runAnalyzer()

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
            self.runAnalyzer()

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
        """Get normalizing constant. Not supported by simulation solver.

        Raises:
            NotImplementedError
        """
        raise NotImplementedError(
            "getProbNormConstAggr() is not supported by SolverJMT. "
            "Use SolverNC or SolverCTMC for normalizing constant computation."
        )

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
            self.runAnalyzer()
        if self._sn is not None and hasattr(self._sn, 'visits') and self._sn.visits:
            return sn_get_residt_from_respt(self._sn, self._result.R, None)
        return self._result.R.copy() if self._result.R is not None else np.array([])

    def getAvgWaitT(self) -> np.ndarray:
        """Get average waiting times (M x K). W = R - S."""
        if self._result is None:
            self.runAnalyzer()
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
            self.runAnalyzer()
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
            self.runAnalyzer()
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
        import os
        import tempfile
        from ....api.solvers.jmt.handler import parse_tran_resp_t

        sn = self._sn
        M = sn.nstations
        K = sn.nclasses

        # Initialize result structure
        RD = [[None for _ in range(K)] for _ in range(M)]

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

        # Initialize model state from marginal distribution
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
            self.runAnalyzer()

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
            self.runAnalyzer()

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
