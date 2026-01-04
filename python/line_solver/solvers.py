"""
LINE Solver classes for queueing network analysis.

This module provides solver classes that use native Python implementations
for analyzing queueing networks.
"""

import os
import warnings
from typing import Optional, Dict, Any, List
import numpy as np
import pandas as pd

from line_solver.constants import VerboseLevel, SolverType


class OptionsDict(dict):
    """A dict that supports attribute-style access.

    Allows both `options['key']` and `options.key` access patterns
    for backward compatibility with notebooks.
    """
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


class SolverOptions:
    """Options for solvers."""

    def __init__(self, solver_type=None):
        self.solver_type = solver_type
        self.method = 'default'
        self.verbose = VerboseLevel.STD
        self.cutoff = 10
        self.seed = 23000
        self.iter_max = 200
        self.iter_tol = 1e-6
        self.tol = 1e-6
        self.samples = 10000
        self.keep = False
        self.force = False
        self.cache = True
        self.hide_immediate = False
        self.stiff = False
        self.timespan = None
        self.timestep = None
        self.init_sol = None


class _DefaultOptionsDescriptor:
    """Descriptor that allows Solver.defaultOptions to work both as property and method."""

    def __get__(self, obj, objtype=None):
        # Return an OptionsDict that can be used directly
        options = OptionsDict()
        options['keep'] = False
        options['verbose'] = VerboseLevel.STD
        options['cutoff'] = 10
        options['seed'] = 23000
        options['iter_max'] = 200
        options['samples'] = 10000
        options['method'] = 'default'
        return options


class Solver:
    """Base class for all LINE solvers.

    Provides common functionality including default options.
    """

    defaultOptions = _DefaultOptionsDescriptor()
    default_options = defaultOptions


class SampleResult:
    """Result container for sample-based analysis."""

    def __init__(self, native_result=None):
        if native_result is not None:
            self.t = getattr(native_result, 't', None)
            self.QN = getattr(native_result, 'QN', None)
            self.UN = getattr(native_result, 'UN', None)
            self.RN = getattr(native_result, 'RN', None)
            self.TN = getattr(native_result, 'TN', None)
            self.WN = getattr(native_result, 'WN', None)
            self.AN = getattr(native_result, 'AN', None)
            self.event = getattr(native_result, 'event', None)
        else:
            self.t = None
            self.QN = None
            self.UN = None
            self.RN = None
            self.TN = None
            self.WN = None
            self.AN = None
            self.event = None


class EventInfo:
    """Event information container."""

    def __init__(self, event_type, station, job_class, time):
        self.event_type = event_type
        self.station = station
        self.job_class = job_class
        self.time = time


class DistributionResult:
    """Result container for distribution analysis."""

    def __init__(self, native_result=None):
        if native_result is not None:
            self.time_points = getattr(native_result, 'time_points', None)
            self.cdf = getattr(native_result, 'cdf', None)
        else:
            self.time_points = None
            self.cdf = None


class ProbabilityResult:
    """Result container for probability analysis."""

    def __init__(self, native_result=None):
        if native_result is not None:
            self.probability = getattr(native_result, 'probability', None)
            self.log_probability = getattr(native_result, 'log_probability', None)
            self.state = getattr(native_result, 'state', None)
        else:
            self.probability = None
            self.log_probability = None
            self.state = None


class Solver:
    """Base class for all queueing network solvers."""

    @staticmethod
    def defaultOptions():
        """Get default solver options."""
        return OptionsDict({
            'keep': False,
            'verbose': VerboseLevel.STD,
            'cutoff': 10,
            'seed': 23000,
            'iter_max': 200,
            'samples': 10000,
            'method': 'default'
        })

    # Snake_case alias
    default_options = defaultOptions

    def __init__(self, options=None, *args, **kwargs):
        self.solveropt = options if options else SolverOptions()
        self._verbose_silent = False
        self._table_silent = False
        self._result = None

        # Process kwargs
        for key, value in kwargs.items():
            self._process_solver_option(key, value)

    def _process_solver_option(self, key, value):
        """Process a single solver option."""
        if hasattr(self.solveropt, key):
            setattr(self.solveropt, key, value)
        if key == 'verbose':
            self._process_verbose_option(value)

    def _process_verbose_option(self, value):
        """Process verbose option."""
        if isinstance(value, bool):
            self._verbose_silent = not value
            self._table_silent = not value
        elif value == VerboseLevel.SILENT:
            self._verbose_silent = True
            self._table_silent = True

    def getName(self):
        """Get solver name."""
        return self.__class__.__name__.replace('Solver', '')

    @property
    def name(self):
        """Get solver name (property access)."""
        return self.getName()

    def supports(self, model):
        """Check if solver supports the model."""
        return True

    get_name = getName
    default_options = defaultOptions


class EnsembleSolver(Solver):
    """Base class for ensemble solvers."""

    def __init__(self, options=None, *args, **kwargs):
        super().__init__(options, *args, **kwargs)
        self._models = []

    def getNumberOfModels(self):
        """Get number of models in ensemble."""
        return len(self._models)

    get_number_of_models = getNumberOfModels


class NetworkSolver(Solver):
    """Base class for network solvers."""

    def __init__(self, model, options=None, *args, **kwargs):
        super().__init__(options, *args, **kwargs)
        self.model = model
        self._native_solver = None

    def runAnalyzer(self):
        """Run the solver analysis."""
        if self._native_solver is not None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result
        return self

    run_analyzer = runAnalyzer

    def getAvgTable(self):
        """Get average performance metrics table."""
        # Auto-run analyzer if needed
        if self._native_solver is not None:
            if hasattr(self._native_solver, 'result') and self._native_solver.result is None:
                self._native_solver.runAnalyzer()
            elif hasattr(self._native_solver, '_result') and self._native_solver._result is None:
                self._native_solver.runAnalyzer()
            if hasattr(self._native_solver, 'getAvgTable'):
                return self._native_solver.getAvgTable()
        return self.avg_table()

    def avg_table(self):
        """Get average performance metrics table (snake_case)."""
        # Auto-run analyzer if needed
        if self._native_solver is not None:
            if hasattr(self._native_solver, 'result') and self._native_solver.result is None:
                self._native_solver.runAnalyzer()
            elif hasattr(self._native_solver, '_result') and self._native_solver._result is None:
                self._native_solver.runAnalyzer()
            # Suppress native solver output to avoid duplicate printing
            if hasattr(self._native_solver, '_table_silent'):
                self._native_solver._table_silent = True
            # Check for getAvgTable method first
            if hasattr(self._native_solver, 'getAvgTable'):
                return self._native_solver.getAvgTable()
            # Check if avg_table is a property or method
            if hasattr(self._native_solver, 'avg_table'):
                attr = getattr(self._native_solver, 'avg_table')
                if callable(attr):
                    return attr()
                else:
                    return attr  # It's a property
        if hasattr(self, '_result') and self._result is not None:
            return self._get_result_table()
        return pd.DataFrame()

    def _get_result_table(self):
        """Convert result to DataFrame."""
        if self._result is None:
            return pd.DataFrame()

        result = self._result
        if hasattr(result, 'Table') and result.Table is not None:
            return result.Table

        # Build table from QN, UN, RN, TN
        rows = []
        station_names = self._get_station_names()
        class_names = self._get_class_names()

        M = len(station_names)
        K = len(class_names)

        for i in range(M):
            for r in range(K):
                qlen = result.QN[i, r] if result.QN is not None and i < result.QN.shape[0] and r < result.QN.shape[1] else 0
                util = result.UN[i, r] if result.UN is not None and i < result.UN.shape[0] and r < result.UN.shape[1] else 0
                respt = result.RN[i, r] if result.RN is not None and i < result.RN.shape[0] and r < result.RN.shape[1] else 0
                tput = result.TN[i, r] if result.TN is not None and i < result.TN.shape[0] and r < result.TN.shape[1] else 0

                rows.append({
                    'Station': station_names[i],
                    'JobClass': class_names[r],
                    'QLen': qlen,
                    'Util': util,
                    'RespT': respt,
                    'ResidT': respt,
                    'ArvR': tput,
                    'Tput': tput,
                })

        return pd.DataFrame(rows)

    def _get_station_names(self):
        """Get station names from model."""
        if hasattr(self.model, 'get_station_names'):
            return self.model.get_station_names()
        if hasattr(self.model, '_sn') and hasattr(self.model._sn, 'nodenames'):
            return list(self.model._sn.nodenames)
        return [f'Station{i}' for i in range(self._get_num_stations())]

    def _get_class_names(self):
        """Get class names from model."""
        if hasattr(self.model, 'get_class_names'):
            return self.model.get_class_names()
        if hasattr(self.model, '_sn') and hasattr(self.model._sn, 'classnames'):
            return list(self.model._sn.classnames)
        return [f'Class{i}' for i in range(self._get_num_classes())]

    def _get_num_stations(self):
        """Get number of stations."""
        if hasattr(self.model, 'get_number_of_stations'):
            return self.model.get_number_of_stations()
        if hasattr(self.model, '_sn'):
            return self.model._sn.nstations
        return 0

    def _get_num_classes(self):
        """Get number of classes."""
        if hasattr(self.model, 'get_number_of_classes'):
            return self.model.get_number_of_classes()
        if hasattr(self.model, '_sn'):
            return self.model._sn.nclasses
        return 0

    def hasResults(self):
        """Check if solver has results."""
        return self._result is not None

    has_results = hasResults

    def getAvgQLen(self):
        """Get average queue length."""
        if self._result is not None and hasattr(self._result, 'QN'):
            return self._result.QN
        return None

    def getAvgUtil(self):
        """Get average utilization."""
        if self._result is not None and hasattr(self._result, 'UN'):
            return self._result.UN
        return None

    def getAvgRespT(self):
        """Get average response time."""
        if self._result is not None and hasattr(self._result, 'RN'):
            return self._result.RN
        return None

    def getAvgTput(self):
        """Get average throughput."""
        if self._result is not None and hasattr(self._result, 'TN'):
            return self._result.TN
        return None

    get_avg_q_len = getAvgQLen
    get_avg_util = getAvgUtil
    get_avg_resp_t = getAvgRespT
    get_avg_tput = getAvgTput
    # Additional aliases without underscores for notebook compatibility
    avg_qlen = getAvgQLen
    avg_util = getAvgUtil
    avg_respt = getAvgRespT
    avg_tput = getAvgTput
    # Node-level table aliases
    avg_node_table = avg_table
    getAvgNodeTable = getAvgTable

    def getAvgChainTable(self):
        """
        Get average performance metrics aggregated by chain.

        A chain is a group of related job classes that can switch between each other.
        This method returns metrics aggregated at the chain level.

        Returns:
            pandas.DataFrame: Table with columns Chain, QLen, Util, RespT, Tput
        """
        # Auto-run analyzer if needed
        if self._native_solver is not None:
            if hasattr(self._native_solver, '_result') and self._native_solver._result is None:
                self._native_solver.runAnalyzer()
            if hasattr(self._native_solver, 'getAvgChainTable'):
                return self._native_solver.getAvgChainTable()
            if hasattr(self._native_solver, 'avg_chain_table'):
                return self._native_solver.avg_chain_table()

        # Fallback: aggregate avg_table by chain if chains exist
        table = self.avg_table()
        if table is None or table.empty:
            return pd.DataFrame()

        # If model has chains, aggregate by chain; otherwise return class-level
        if hasattr(self.model, 'get_chains') and self.model.get_chains():
            # Group by chain and aggregate
            chains = self.model.get_chains()
            rows = []
            for chain_idx, chain in enumerate(chains):
                chain_name = chain.name if hasattr(chain, 'name') else f'Chain{chain_idx}'
                chain_classes = chain.classes if hasattr(chain, 'classes') else []
                class_names = [c.name if hasattr(c, 'name') else str(c) for c in chain_classes]

                # Filter table for this chain's classes
                if 'JobClass' in table.columns:
                    chain_data = table[table['JobClass'].isin(class_names)]
                    if not chain_data.empty:
                        rows.append({
                            'Chain': chain_name,
                            'QLen': chain_data['QLen'].sum() if 'QLen' in chain_data else 0,
                            'Util': chain_data['Util'].mean() if 'Util' in chain_data else 0,
                            'RespT': chain_data['RespT'].mean() if 'RespT' in chain_data else 0,
                            'Tput': chain_data['Tput'].sum() if 'Tput' in chain_data else 0,
                        })
            return pd.DataFrame(rows) if rows else table
        else:
            # No chains defined, return per-class data
            return table

    avg_chain_table = getAvgChainTable

    def getAvgSysTable(self):
        """
        Get system-level average performance metrics.

        Returns system-wide metrics like total response time and throughput
        aggregated across all stations.

        Returns:
            pandas.DataFrame: Table with columns Chain, SysRespT, SysTput
        """
        # Auto-run analyzer if needed
        if self._native_solver is not None:
            if hasattr(self._native_solver, '_result') and self._native_solver._result is None:
                self._native_solver.runAnalyzer()
            if hasattr(self._native_solver, 'getAvgSysTable'):
                return self._native_solver.getAvgSysTable()
            if hasattr(self._native_solver, 'avg_sys_table'):
                return self._native_solver.avg_sys_table()

        # Fallback: compute from chain table
        chain_table = self.getAvgChainTable()
        if chain_table is None or chain_table.empty:
            return pd.DataFrame()

        # System response time is sum of response times, system throughput is the throughput
        rows = []
        if 'Chain' in chain_table.columns:
            for _, row in chain_table.iterrows():
                rows.append({
                    'Chain': row.get('Chain', 'Chain0'),
                    'SysRespT': row.get('RespT', 0),
                    'SysTput': row.get('Tput', 0),
                })
        else:
            rows.append({
                'Chain': 'Chain0',
                'SysRespT': chain_table['RespT'].sum() if 'RespT' in chain_table else 0,
                'SysTput': chain_table['Tput'].mean() if 'Tput' in chain_table else 0,
            })

        return pd.DataFrame(rows)

    avg_sys_table = getAvgSysTable

    def getAvg(self):
        """
        Get all average metrics (QLen, Util, RespT, Tput, ArvR, ResidT).

        Returns:
            Tuple of (QN, UN, RN, TN, AN, WN) arrays
        """
        if self._native_solver is not None:
            if hasattr(self._native_solver, '_result') and self._native_solver._result is None:
                self._native_solver.runAnalyzer()
            if hasattr(self._native_solver, 'getAvg'):
                return self._native_solver.getAvg()
        result = self._result
        if result is not None:
            return (
                getattr(result, 'QN', None),
                getattr(result, 'UN', None),
                getattr(result, 'RN', None),
                getattr(result, 'TN', None),
                getattr(result, 'AN', None),
                getattr(result, 'WN', None),
            )
        return (None, None, None, None, None, None)

    avg = getAvg

    def getAvgChain(self):
        """Get average metrics aggregated by chain."""
        if self._native_solver is not None and hasattr(self._native_solver, 'getAvgChain'):
            return self._native_solver.getAvgChain()
        # Return from chain table
        table = self.getAvgChainTable()
        if table.empty:
            return (None, None, None, None, None, None)
        return (
            table['QLen'].values if 'QLen' in table else None,
            table['Util'].values if 'Util' in table else None,
            table['RespT'].values if 'RespT' in table else None,
            None,  # WN
            None,  # AN
            table['Tput'].values if 'Tput' in table else None,
        )

    avg_chain = getAvgChain

    def getAvgNode(self):
        """Get average metrics at node level."""
        return self.getAvg()

    avg_node = getAvgNode

    def getAvgSys(self):
        """Get system-level average metrics (response time and throughput per chain)."""
        if self._native_solver is not None and hasattr(self._native_solver, 'getAvgSys'):
            return self._native_solver.getAvgSys()
        table = self.getAvgSysTable()
        if table.empty:
            return (None, None)
        return (
            table['SysRespT'].values if 'SysRespT' in table else None,
            table['SysTput'].values if 'SysTput' in table else None,
        )

    avg_sys = getAvgSys

    def getAvgSysRespT(self):
        """Get system response time per chain."""
        sys_data = self.getAvgSys()
        return sys_data[0] if sys_data else None

    avg_sys_resp_t = getAvgSysRespT

    def getAvgSysTput(self):
        """Get system throughput per chain."""
        sys_data = self.getAvgSys()
        return sys_data[1] if sys_data else None

    avg_sys_tput = getAvgSysTput

    def getAvgQLenChain(self):
        """Get average queue length per chain."""
        table = self.getAvgChainTable()
        return table['QLen'].values if not table.empty and 'QLen' in table else None

    avg_qlen_chain = getAvgQLenChain

    def getAvgUtilChain(self):
        """Get average utilization per chain."""
        table = self.getAvgChainTable()
        return table['Util'].values if not table.empty and 'Util' in table else None

    avg_util_chain = getAvgUtilChain

    def getAvgRespTChain(self):
        """Get average response time per chain."""
        table = self.getAvgChainTable()
        return table['RespT'].values if not table.empty and 'RespT' in table else None

    avg_resp_t_chain = getAvgRespTChain

    def getAvgResidTChain(self):
        """Get average residence time per chain."""
        return self.getAvgRespTChain()  # Same as response time for chain-level

    avg_resid_t_chain = getAvgResidTChain

    def getAvgTputChain(self):
        """Get average throughput per chain."""
        table = self.getAvgChainTable()
        return table['Tput'].values if not table.empty and 'Tput' in table else None

    avg_tput_chain = getAvgTputChain

    def getAvgArvRChain(self):
        """Get average arrival rate per chain."""
        return self.getAvgTputChain()  # For closed networks, arrival = throughput

    avg_arv_r_chain = getAvgArvRChain

    def getAvgNodeChain(self):
        """Get average metrics at node level, aggregated by chain."""
        return self.getAvgChain()

    avg_node_chain = getAvgNodeChain

    def getAvgNodeChainTable(self):
        """Get node-level metrics aggregated by chain as a table."""
        return self.getAvgChainTable()

    avg_node_chain_table = getAvgNodeChainTable

    def getAvgNodeQLenChain(self):
        """Get average queue length at nodes, per chain."""
        return self.getAvgQLenChain()

    avg_node_qlen_chain = getAvgNodeQLenChain

    def getAvgNodeUtilChain(self):
        """Get average utilization at nodes, per chain."""
        return self.getAvgUtilChain()

    avg_node_util_chain = getAvgNodeUtilChain

    def getAvgNodeRespTChain(self):
        """Get average response time at nodes, per chain."""
        return self.getAvgRespTChain()

    avg_node_resp_t_chain = getAvgNodeRespTChain

    def getAvgNodeResidTChain(self):
        """Get average residence time at nodes, per chain."""
        return self.getAvgResidTChain()

    avg_node_resid_t_chain = getAvgNodeResidTChain

    def getAvgNodeTputChain(self):
        """Get average throughput at nodes, per chain."""
        return self.getAvgTputChain()

    avg_node_tput_chain = getAvgNodeTputChain

    def getAvgNodeArvRChain(self):
        """Get average arrival rate at nodes, per chain."""
        return self.getAvgArvRChain()

    avg_node_arv_r_chain = getAvgNodeArvRChain

    def getAvgResidT(self):
        """Get average residence time."""
        if self._result is not None and hasattr(self._result, 'WN'):
            return self._result.WN
        return self.getAvgRespT()  # Fallback to response time

    avg_resid_t = getAvgResidT

    def getAvgArvR(self):
        """Get average arrival rate."""
        if self._result is not None and hasattr(self._result, 'AN'):
            return self._result.AN
        return self.getAvgTput()  # Fallback to throughput

    avg_arv_r = getAvgArvR

    def getAvgT(self):
        """Alias for getAvgTable."""
        return self.getAvgTable()

    avg_t = getAvgT

    def getChainAvgT(self):
        """Alias for getAvgChainTable."""
        return self.getAvgChainTable()

    chain_avg_t = getChainAvgT

    def getNodeAvgT(self):
        """Alias for getAvgNodeTable."""
        return self.getAvgNodeTable()

    node_avg_t = getNodeAvgT

    def getNodeChainAvgT(self):
        """Alias for getAvgNodeChainTable."""
        return self.getAvgNodeChainTable()

    node_chain_avg_t = getNodeChainAvgT

    def getSysAvgT(self):
        """Alias for getAvgSysTable."""
        return self.getAvgSysTable()

    sys_avg_t = getSysAvgT

    def getPerctRespT(self, percentiles=None, jobclass=None):
        """
        Get response time percentiles.

        Args:
            percentiles: List of percentiles (e.g., [90, 95, 99])
            jobclass: Optional job class filter

        Returns:
            Tuple of (percentile_values, percentile_table)
        """
        if self._native_solver is not None and hasattr(self._native_solver, 'getPerctRespT'):
            return self._native_solver.getPerctRespT(percentiles, jobclass)
        # Fallback: return empty if not supported
        return (np.array([]), pd.DataFrame())

    perct_resp_t = getPerctRespT
    get_perct_resp_t = getPerctRespT

    def getTranAvg(self):
        """Get transient average metrics."""
        if self._native_solver is not None and hasattr(self._native_solver, 'getTranAvg'):
            return self._native_solver.getTranAvg()
        return (None, None, None)

    tran_avg = getTranAvg
    get_tran_avg = getTranAvg

    def getTranHandles(self):
        """Get transient metric handles."""
        if self._native_solver is not None and hasattr(self._native_solver, 'getTranHandles'):
            return self._native_solver.getTranHandles()
        return (None, None, None)

    tran_handles = getTranHandles
    get_tran_handles = getTranHandles

    def getAvgHandles(self):
        """Get average metric handles."""
        if self._native_solver is not None and hasattr(self._native_solver, 'getAvgHandles'):
            return self._native_solver.getAvgHandles()
        return (None, None, None, None, None, None)

    avg_handles = getAvgHandles
    get_avg_handles = getAvgHandles

    def getAvgQLenHandles(self):
        """Get queue length handles."""
        handles = self.getAvgHandles()
        return handles[0] if handles else None

    avg_qlen_handles = getAvgQLenHandles

    def getAvgUtilHandles(self):
        """Get utilization handles."""
        handles = self.getAvgHandles()
        return handles[1] if handles else None

    avg_util_handles = getAvgUtilHandles

    def getAvgRespTHandles(self):
        """Get response time handles."""
        handles = self.getAvgHandles()
        return handles[2] if handles else None

    avg_resp_t_handles = getAvgRespTHandles

    def getAvgTputHandles(self):
        """Get throughput handles."""
        handles = self.getAvgHandles()
        return handles[3] if handles else None

    avg_tput_handles = getAvgTputHandles

    def getAvgArvRHandles(self):
        """Get arrival rate handles."""
        handles = self.getAvgHandles()
        return handles[4] if handles else None

    avg_arv_r_handles = getAvgArvRHandles

    def getAvgResidTHandles(self):
        """Get residence time handles."""
        handles = self.getAvgHandles()
        return handles[5] if handles else None

    avg_resid_t_handles = getAvgResidTHandles

    def getProbNormConstAggr(self):
        """Get probability normalizing constant (aggregate)."""
        if self._native_solver is not None and hasattr(self._native_solver, 'getProbNormConstAggr'):
            return self._native_solver.getProbNormConstAggr()
        return None

    prob_norm_const_aggr = getProbNormConstAggr

    def getProb(self, node=None, state=None):
        """Get probability of a specific state at a node."""
        if self._native_solver is not None and hasattr(self._native_solver, 'getProb'):
            return self._native_solver.getProb(node, state)
        return None

    prob = getProb

    def getProbSys(self):
        """Get system state probabilities."""
        if self._native_solver is not None and hasattr(self._native_solver, 'getProbSys'):
            return self._native_solver.getProbSys()
        return None

    prob_sys = getProbSys

    def prob_sys_aggr(self):
        """
        Get system aggregate state probabilities.

        Returns:
            Tuple of (Pr, S) where:
            - Pr: Probability vector for each aggregate state
            - S: Array of aggregate state descriptors
        """
        if self._native_solver is None:
            return np.array([]), np.array([])

        if hasattr(self._native_solver, '_result') and self._native_solver._result is None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result

        # Try native solver methods
        if hasattr(self._native_solver, 'getProbSysAggr'):
            return self._native_solver.getProbSysAggr()

        if hasattr(self._native_solver, 'getSteadyState'):
            pi = self._native_solver.getSteadyState()
            space = self._native_solver.getStateSpaceAggr() if hasattr(self._native_solver, 'getStateSpaceAggr') else np.array([])
            if len(pi) > 0:
                return pi, space

        return np.array([]), np.array([])

    def prob_aggr(self, node=None):
        """
        Get aggregate state probabilities for a specific node or all nodes.

        Args:
            node: Node index (0-based) or node object. If None, returns all.

        Returns:
            Probability vector for aggregate states at the node
        """
        if self._native_solver is None:
            return np.array([])

        if hasattr(self._native_solver, '_result') and self._native_solver._result is None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result

        if hasattr(self._native_solver, 'getProb'):
            if node is not None:
                node_idx = node if isinstance(node, int) else self._get_node_index(node)
                return self._native_solver.getProb(node_idx)
            else:
                return self._native_solver.getProb()

        return np.array([])

    def _get_node_index(self, node):
        """Get node index from node object."""
        if hasattr(self.model, 'get_nodes'):
            nodes = self.model.get_nodes()
            for i, n in enumerate(nodes):
                if n == node or (hasattr(n, 'name') and hasattr(node, 'name') and n.name == node.name):
                    return i
        return 0


class SolverMVA(NetworkSolver):
    """Mean Value Analysis solver."""

    def __init__(self, model, *args, **kwargs):
        options = SolverOptions(SolverType.MVA)
        super().__init__(model, options, *args, **kwargs)

        from .solvers.solver_mva import SolverMVA
        method = kwargs.pop('method', 'exact')  # Use 'exact' for accurate results
        # Filter out unsupported kwargs for MVA solver
        supported_mva_kwargs = {'verbose', 'seed', 'max_iter', 'tol'}
        mva_kwargs = {k: v for k, v in kwargs.items() if k in supported_mva_kwargs}
        self._native_solver = SolverMVA(model, method=method, **mva_kwargs)

    @staticmethod
    def defaultOptions():
        """Get default MVA options."""
        return OptionsDict({
            'method': 'exact',  # Use exact MVA (matches CTMC results)
            'verbose': VerboseLevel.STD,
            'iter_max': 1000,
            'tol': 1e-8,
            'config': OptionsDict({
                'multiserver': 'default',
                'fork_join': 'default'
            })
        })

    default_options = defaultOptions


class SolverNC(NetworkSolver):
    """Normalizing Constant solver."""

    def __init__(self, model, *args, **kwargs):
        options = SolverOptions(SolverType.NC)
        super().__init__(model, options, *args, **kwargs)

        from .solvers.solver_nc import SolverNC
        method = kwargs.pop('method', 'exact')
        # Filter out unsupported kwargs for NC solver
        supported_nc_kwargs = {'tol', 'iter_max', 'iter_tol', 'verbose', 'seed'}
        nc_kwargs = {k: v for k, v in kwargs.items() if k in supported_nc_kwargs}
        self._native_solver = SolverNC(model, method=method, **nc_kwargs)

    def prob_sys_aggr(self):
        """
        Get system aggregate state probabilities.

        Returns:
            Tuple of (Pr, S) where:
            - Pr: Probability vector for each aggregate state
            - S: Array of aggregate state descriptors
        """
        if self._native_solver._result is None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result

        # NC solver computes marginal probabilities
        if hasattr(self._native_solver, 'getProbSysAggr'):
            return self._native_solver.getProbSysAggr()

        # Fallback: construct from steady-state if available
        if hasattr(self._native_solver, 'getSteadyState'):
            pi = self._native_solver.getSteadyState()
            space = self._native_solver.getStateSpaceAggr() if hasattr(self._native_solver, 'getStateSpaceAggr') else np.array([])
            return pi, space

        return np.array([]), np.array([])

    @staticmethod
    def defaultOptions():
        """Get default NC options."""
        return OptionsDict({
            'method': 'exact',
            'verbose': VerboseLevel.STD
        })

    default_options = defaultOptions


class SolverCTMC(NetworkSolver):
    """Continuous-Time Markov Chain solver."""

    def __init__(self, model, *args, **kwargs):
        options = SolverOptions(SolverType.CTMC)
        super().__init__(model, options, *args, **kwargs)

        from .solvers.solver_ctmc import SolverCTMC
        method = kwargs.pop('method', 'default')
        self._native_solver = SolverCTMC(model, method=method, **kwargs)

    def state_space(self):
        """
        Get the CTMC state space.

        Returns:
            Tuple of (stateSpace, nodeStateSpace) where:
            - stateSpace: 2D array of system states (nstates x nstations)
            - nodeStateSpace: Dict mapping node index to node-specific state arrays
        """
        if self._native_solver._result is None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result

        space = self._native_solver.getStateSpace()
        sn = self._native_solver._sn

        # Build nodeStateSpace dictionary
        nstations = sn.nstations if sn else space.shape[1]
        nodeStateSpace = {}

        for i in range(nstations):
            # Extract column for this station
            if space.ndim == 2 and space.shape[1] > i:
                nodeStateSpace[i] = space[:, i:i+1]
            else:
                nodeStateSpace[i] = space.reshape(-1, 1)

        return space, nodeStateSpace

    def generator(self):
        """
        Get the CTMC infinitesimal generator matrix.

        Returns:
            Tuple of (infGen, eventFilt) where:
            - infGen: Infinitesimal generator matrix Q
            - eventFilt: Event filter matrix (for now, same as infGen)
        """
        if self._native_solver._result is None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result

        infGen = self._native_solver.getInfGen()
        eventFilt = infGen.copy()  # Event filter is the same for basic CTMC

        return infGen, eventFilt

    @staticmethod
    def print_inf_gen(infGen, stateSpace):
        """
        Print the infinitesimal generator in human-readable form.

        Args:
            infGen: Infinitesimal generator matrix
            stateSpace: State space array
        """
        nstates = infGen.shape[0]
        for i in range(nstates):
            for j in range(nstates):
                if i != j and infGen[i, j] > 0:
                    from_state = stateSpace[i, :] if stateSpace.ndim > 1 else [stateSpace[i]]
                    to_state = stateSpace[j, :] if stateSpace.ndim > 1 else [stateSpace[j]]
                    from_str = "[ " + "   ".join(f"{s:.1f}" for s in from_state) + " ]"
                    to_str = "[ " + "   ".join(f"{s:.1f}" for s in to_state) + " ]"
                    print(f"{from_str} -> {to_str} : {infGen[i, j]}")

    def avg_util(self):
        """
        Get average utilization as a dictionary keyed by node.

        Returns:
            Dict mapping node to utilization array
        """
        if self._native_solver._result is None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result

        U = self._native_solver.getAvgUtil()
        result = {}
        nodes = self.model.get_nodes() if hasattr(self.model, 'get_nodes') else []

        for i, node in enumerate(nodes):
            if i < U.shape[0]:
                result[node] = U[i, :]

        return result

    def avg_wait_t(self):
        """
        Get average waiting time as a dictionary keyed by node index.

        Returns:
            Dict mapping node index to waiting time array
        """
        if self._native_solver._result is None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result

        W = self._native_solver.getAvgWaitT()
        result = {}

        for i in range(W.shape[0]):
            result[i] = W[i, :]

        return result

    def sample_sys_aggr(self, num_samples: int = 1000):
        """
        Sample system aggregate states (simulation-based).

        Args:
            num_samples: Number of events to sample

        Returns:
            SampleResult containing event trace
        """
        # For CTMC, we use Gillespie algorithm to sample from the Markov chain
        if self._native_solver._result is None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result

        from dataclasses import dataclass
        import numpy as np

        @dataclass
        class Event:
            t: float
            node: int
            event: str
            job_class: int = 0

        infGen = self._native_solver.getInfGen()
        space = self._native_solver.getStateSpace()
        nstates = infGen.shape[0]

        if nstates == 0:
            return SampleResult()

        # Initialize at steady state or state 0
        pi = self._native_solver.getSteadyState()
        if len(pi) > 0:
            current_state = np.random.choice(nstates, p=pi)
        else:
            current_state = 0

        events = []
        t = 0.0
        np.random.seed(self.solveropt.seed if hasattr(self.solveropt, 'seed') else 23000)

        for _ in range(num_samples):
            # Get transition rates from current state
            rates = infGen[current_state, :].copy()
            rates[current_state] = 0  # No self-transitions

            total_rate = np.sum(rates)
            if total_rate <= 0:
                break

            # Time to next event (exponential)
            dt = np.random.exponential(1.0 / total_rate)
            t += dt

            # Choose next state
            probs = rates / total_rate
            next_state = np.random.choice(nstates, p=probs)

            # Determine event type based on state change
            if space.ndim > 1 and space.shape[0] > max(current_state, next_state):
                old_jobs = space[current_state, :]
                new_jobs = space[next_state, :]
                diff = new_jobs - old_jobs

                # Get station-to-node mapping for proper node indices
                sn = self._native_solver._sn
                stationToNode = sn.stationToNode if hasattr(sn, 'stationToNode') and sn.stationToNode is not None else None

                # Find which station changed and convert to node index
                for station_idx in range(len(diff)):
                    if diff[station_idx] != 0:
                        # Convert station index to node index
                        if stationToNode is not None and station_idx < len(stationToNode):
                            node_idx = int(stationToNode[station_idx])
                        else:
                            node_idx = station_idx

                        if diff[station_idx] > 0:
                            events.append(Event(t=t, node=node_idx, event="ARV"))
                        elif diff[station_idx] < 0:
                            events.append(Event(t=t, node=node_idx, event="DEP"))

            current_state = next_state

        result = SampleResult()
        result.event = events
        result.t = np.array([e.t for e in events])

        return result

    def prob_sys_aggr(self):
        """
        Get system aggregate state probabilities.

        Returns:
            Tuple of (Pr, S) where:
            - Pr: Probability vector for each aggregate state
            - S: Array of aggregate state descriptors
        """
        if self._native_solver._result is None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result

        pi = self._native_solver.getSteadyState()
        space = self._native_solver.getStateSpaceAggr()

        if len(pi) == 0 or len(space) == 0:
            return np.array([]), np.array([])

        return pi, space

    def prob_aggr(self, node=None):
        """
        Get aggregate state probabilities for a specific node or all nodes.

        Args:
            node: Node index (0-based) or node object. If None, returns all.

        Returns:
            Probability vector for aggregate states at the node
        """
        if self._native_solver._result is None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result

        if node is not None:
            node_idx = node if isinstance(node, int) else self._get_node_index(node)
            return self._native_solver.getProb(node_idx)
        else:
            return self._native_solver.getProb()

    def _get_node_index(self, node):
        """Get node index from node object."""
        if hasattr(self.model, 'get_nodes'):
            nodes = self.model.get_nodes()
            for i, n in enumerate(nodes):
                if n == node or (hasattr(n, 'name') and hasattr(node, 'name') and n.name == node.name):
                    return i
        return 0

    def get_avg_reward(self):
        """
        Compute steady-state expected rewards for all reward functions defined on the model.

        The model must have reward functions defined via model.setReward(name, reward_fn).
        Each reward function takes a RewardState object and returns a scalar value.

        Returns:
            Tuple of (R, names) where:
            - R: numpy array of expected reward values
            - names: list of reward function names

        Example:
            >>> model.setReward('QueueLength', lambda state: state.at(queue, oclass))
            >>> model.setReward('Utilization', Reward.utilization(queue, oclass))
            >>> solver = CTMC(model)
            >>> R, names = solver.get_avg_reward()
        """
        if self._native_solver._result is None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result

        # Get rewards from the model
        if hasattr(self.model, 'getRewards'):
            rewards = self.model.getRewards()
        elif hasattr(self.model, 'get_rewards'):
            rewards = self.model.get_rewards()
        else:
            rewards = {}

        if not rewards:
            return np.array([]), []

        # Get state space and steady-state probabilities
        space = self._native_solver.getStateSpace()
        pi = self._native_solver.getSteadyState()

        if len(space) == 0 or len(pi) == 0:
            return np.array([0.0] * len(rewards)), list(rewards.keys())

        # Build reward state accessor
        from .lang.reward_state import RewardState
        sn = self._native_solver._sn

        # Build node-to-station and class-to-index mappings
        nodes_to_station = {}
        if hasattr(self.model, 'get_nodes'):
            for node in self.model.get_nodes():
                # Get node index - handle both wrapper and native implementations
                if hasattr(node, 'get_index'):
                    node_idx = node.get_index()  # 1-based
                elif hasattr(node, 'index'):
                    node_idx = node.index
                else:
                    continue

                if hasattr(sn, 'nodeToStation'):
                    node_idx0 = node_idx - 1  # 0-indexed
                    if node_idx0 < len(sn.nodeToStation):
                        station_idx = int(sn.nodeToStation[node_idx0])
                        if station_idx >= 0:  # Valid station
                            nodes_to_station[node_idx] = station_idx + 1  # 1-based for RewardState

        classes_to_idx = {}
        if hasattr(self.model, 'get_classes'):
            for i, jobclass in enumerate(self.model.get_classes()):
                if hasattr(jobclass, 'get_index'):
                    class_idx = jobclass.get_index()
                elif hasattr(jobclass, 'index'):
                    class_idx = jobclass.index
                else:
                    class_idx = i + 1  # 1-based
                classes_to_idx[class_idx] = i + 1  # 1-indexed

        # Compute expected reward for each reward function
        R = []
        names = list(rewards.keys())

        for name, reward_fn in rewards.items():
            expected_value = 0.0

            for state_idx, state_vec in enumerate(space):
                prob = pi[state_idx]
                if prob <= 0:
                    continue

                # Create RewardState for this state vector
                reward_state = RewardState(state_vec, sn, nodes_to_station, classes_to_idx)

                # Evaluate reward function
                try:
                    reward_value = reward_fn(reward_state)
                    expected_value += prob * reward_value
                except Exception as e:
                    # Skip if reward function fails for this state
                    pass

            R.append(expected_value)

        return np.array(R), names

    # Alias for camelCase
    getAvgReward = get_avg_reward

    @staticmethod
    def defaultOptions():
        """Get default CTMC options."""
        return OptionsDict({
            'method': 'default',
            'cutoff': 10,
            'seed': 23000,
            'verbose': VerboseLevel.STD
        })

    default_options = defaultOptions


class SolverSSA(NetworkSolver):
    """Stochastic Simulation Algorithm solver."""

    def __init__(self, model, *args, **kwargs):
        options = SolverOptions(SolverType.SSA)
        super().__init__(model, options, *args, **kwargs)

        from .solvers.solver_ssa import SolverSSA
        method = kwargs.pop('method', 'default')
        self._native_solver = SolverSSA(model, method=method, **kwargs)

    @staticmethod
    def defaultOptions():
        """Get default SSA options."""
        return OptionsDict({
            'method': 'default',
            'samples': 10000,
            'seed': 23000,
            'verbose': VerboseLevel.STD
        })

    default_options = defaultOptions


class SolverFluid(NetworkSolver):
    """Fluid approximation solver."""

    def __init__(self, model, *args, **kwargs):
        options = SolverOptions(SolverType.FLUID)
        super().__init__(model, options, *args, **kwargs)

        from .solvers.solver_fld import SolverFLD
        from .solvers.solver_fld.options import SolverFLDOptions
        method = kwargs.pop('method', 'default')
        # Filter out unsupported kwargs for FLD solver
        supported_fld_kwargs = {'tol', 'iter_max', 'iter_tol', 'stiff', 'timespan',
                                'timestep', 'pstar', 'softmin_alpha', 'hide_immediate', 'verbose'}
        fld_kwargs = {k: v for k, v in kwargs.items() if k in supported_fld_kwargs}
        # Convert kwargs to SolverFLDOptions
        fld_options = SolverFLDOptions(method=method, **fld_kwargs) if fld_kwargs else None
        self._native_solver = SolverFLD(model, method=method, options=fld_options)

    def runAnalyzer(self):
        """Run the fluid solver analysis."""
        if self._native_solver is not None:
            self._native_solver.runAnalyzer()
            # SolverFLD uses 'result' not '_result'
            self._result = self._native_solver.result
        return self

    run_analyzer = runAnalyzer

    def cdf_resp_t(self):
        """
        Get response time CDF for all station/class combinations.

        Returns:
            Matrix of CDF arrays indexed by [station][class].
            Each CDF array has columns [cdf_value, time].
        """
        if self._native_solver.result is None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver.result

        import numpy as np

        sn = self._native_solver.sn
        nstations = sn.nstations if sn else 1
        nclasses = sn.nclasses if sn else 1

        RD_matrix = [[None for _ in range(nclasses)] for _ in range(nstations)]

        # Get CDF for each station/class
        for station in range(nstations):
            for job_class in range(nclasses):
                try:
                    cdf_result = self._native_solver.getCdfRespT(station=station, job_class=job_class)
                    t = cdf_result.get('t', np.array([]))
                    cdf = cdf_result.get('cdf', np.array([]))

                    if len(t) > 0 and len(cdf) > 0:
                        # Create array with columns [cdf, time]
                        cdf_array = np.column_stack([cdf, t])
                        RD_matrix[station][job_class] = cdf_array
                except Exception as e:
                    # Skip if CDF computation fails for this station/class
                    print(f"Error in getCdfRespT: {e}")
                    continue

        return RD_matrix

    # Alias without underscore for consistency
    cdf_respt = cdf_resp_t

    @staticmethod
    def defaultOptions():
        """Get default Fluid options."""
        return OptionsDict({
            'method': 'default',
            'verbose': VerboseLevel.STD
        })

    default_options = defaultOptions


class SolverMAM(NetworkSolver):
    """Matrix Analytic Methods solver."""

    def __init__(self, model, *args, **kwargs):
        options = SolverOptions(SolverType.MAM)
        super().__init__(model, options, *args, **kwargs)

        from .solvers.solver_mam import SolverMAM
        method = kwargs.pop('method', 'default')
        self._native_solver = SolverMAM(model, method=method, **kwargs)

    @staticmethod
    def defaultOptions():
        """Get default MAM options."""
        return OptionsDict({
            'method': 'default',
            'verbose': VerboseLevel.STD
        })

    default_options = defaultOptions


class SolverJMT(NetworkSolver):
    """JMT simulation solver."""

    def __init__(self, model, *args, **kwargs):
        options = SolverOptions(SolverType.JMT)
        super().__init__(model, options, *args, **kwargs)

        # Handle options object passed as first positional arg
        if args and isinstance(args[0], (dict, OptionsDict)):
            opts = args[0]
            kwargs.update(opts)
            args = args[1:]

        from .solvers.solver_jmt import SolverJMT
        self._native_solver = SolverJMT(model, **kwargs)

    def cdf_resp_t(self):
        """
        Get response time CDF for all station/class combinations.

        Note: JMT native solver currently does not support CDF extraction.
        Returns empty matrix. Use SolverFluid (FLD) for CDF analysis.

        Returns:
            Matrix of CDF arrays indexed by [station][class].
        """
        import warnings
        warnings.warn(
            "SolverJMT native does not currently support CDF extraction. "
            "Returning empty results. Use SolverFluid (FLD) for CDF analysis.",
            UserWarning
        )

        if self._native_solver._result is None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result

        sn = self._native_solver._sn
        nstations = sn.nstations if sn else 1
        nclasses = sn.nclasses if sn else 1

        # Return empty matrix since JMT doesn't provide CDF data yet
        RD_matrix = [[None for _ in range(nclasses)] for _ in range(nstations)]
        return RD_matrix

    # Alias without underscore for consistency
    cdf_respt = cdf_resp_t

    def get_tran_cdf_respt(self):
        """Get transient response time CDF (stub - returns empty results)."""
        import warnings
        warnings.warn(
            "SolverJMT native does not currently support transient CDF extraction. "
            "Returning empty results.",
            UserWarning
        )
        return self.cdf_resp_t()

    getTranCdfRespT = get_tran_cdf_respt

    @staticmethod
    def defaultOptions():
        """Get default JMT options."""
        return OptionsDict({
            'samples': 10000,
            'seed': 23000,
            'verbose': VerboseLevel.STD
        })

    default_options = defaultOptions


class SolverDES(NetworkSolver):
    """Discrete Event Simulation solver."""

    def __init__(self, model, *args, **kwargs):
        options = SolverOptions(SolverType.DES)
        super().__init__(model, options, *args, **kwargs)

        from .solvers.solver_des import SolverDES
        self._native_solver = SolverDES(model, **kwargs)

    @staticmethod
    def defaultOptions():
        """Get default DES options."""
        return OptionsDict({
            'samples': 10000,
            'seed': 23000,
            'verbose': VerboseLevel.STD
        })

    default_options = defaultOptions


class SolverAuto(NetworkSolver):
    """Automatic solver selection."""

    def __init__(self, model, *args, **kwargs):
        options = SolverOptions(SolverType.AUTO)
        super().__init__(model, options, *args, **kwargs)

        from .solvers.solver_auto import SolverAuto
        self._native_solver = SolverAuto(model, **kwargs)

    @staticmethod
    def defaultOptions():
        """Get default AUTO options."""
        return OptionsDict({
            'verbose': VerboseLevel.STD
        })

    default_options = defaultOptions


class LINE(SolverAuto):
    """Alias for SolverAuto - automatic solver selection."""

    @staticmethod
    def load(filepath):
        """
        Load a model from a pickle file.

        Args:
            filepath: Path to the pickle file

        Returns:
            The loaded object (typically a Network model or dict)
        """
        import pickle
        with open(filepath, 'rb') as f:
            return pickle.load(f)


class SolverENV(EnsembleSolver):
    """Environment solver for random environments."""

    def __init__(self, model, solvers, *args, **kwargs):
        options = kwargs.get('options', SolverOptions(SolverType.ENV))
        if 'options' not in kwargs:
            options = SolverOptions(SolverType.ENV)
        super().__init__(options, *args, **kwargs)
        self.model = model
        self._solvers = solvers
        self._smp_method = False  # Use DTMC-based computation for Semi-Markov Processes

        # Enable SMP method if specified in options
        if hasattr(options, 'method') and options.method is not None:
            if str(options.method).lower() == 'smp':
                self._smp_method = True

    def setSMPMethod(self, flag: bool):
        """
        Enable/disable DTMC-based computation for Semi-Markov Processes.

        Args:
            flag: True to enable SMP method, False to disable
        """
        self._smp_method = flag

    def set_smp_method(self, flag: bool):
        """
        Enable/disable DTMC-based computation for Semi-Markov Processes (Pythonic name).

        Args:
            flag: True to enable SMP method, False to disable
        """
        self.setSMPMethod(flag)

    @staticmethod
    def defaultOptions():
        """Get default ENV options."""
        return OptionsDict({
            'verbose': VerboseLevel.STD
        })

    default_options = defaultOptions


class ENV(SolverENV):
    """Alias for SolverENV."""
    pass


class SolverQNS(NetworkSolver):
    """QNS solver."""

    def __init__(self, model, *args, **kwargs):
        # Extract QNS-specific kwargs before calling super()
        method = kwargs.pop('method', 'default')
        qns_kwargs = {k: kwargs.pop(k) for k in list(kwargs.keys()) if k in ['multiserver', 'samples', 'keep']}

        options = SolverOptions(SolverType.QNS)
        super().__init__(model, options, *args, **kwargs)

        from .solvers.solver_qns import SolverQNS as NativeSolverQNS
        from .solvers.solver_qns.solver_qns import QNSOptions
        qns_options = QNSOptions(method=method, **qns_kwargs)
        # Get the network structure from the model
        sn = model.getStruct(True) if hasattr(model, 'getStruct') else model._sn
        self._native_solver = NativeSolverQNS(sn, options=qns_options)

    @staticmethod
    def defaultOptions():
        """Get default QNS options."""
        return OptionsDict({
            'verbose': VerboseLevel.STD
        })

    @staticmethod
    def isAvailable():
        """Check if QNS solver is available."""
        try:
            from .solvers.solver_qns import SolverQNS
            return True
        except ImportError:
            return False

    default_options = defaultOptions
    is_available = isAvailable


class SolverLQNS(Solver):
    """LQNS solver for layered queueing networks."""

    def __init__(self, model, *args, **kwargs):
        options = SolverOptions(SolverType.LQNS)
        super().__init__(options, *args, **kwargs)
        self.model = model

        from .solvers.solver_lqns import SolverLQNS
        self._native_solver = SolverLQNS(model, **kwargs)

    def runAnalyzer(self):
        """Run the solver analysis."""
        if self._native_solver is not None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result
        return self

    run_analyzer = runAnalyzer

    def avg_table(self):
        """Get average performance metrics table."""
        if self._native_solver is not None:
            return self._native_solver.getAvgTable()
        return None

    def getAvgTable(self):
        """Get average performance metrics table."""
        return self.avg_table()

    def get_raw_avg_tables(self):
        """Get raw average tables (AvgTable, CallAvgTable)."""
        if self._native_solver is not None:
            avg_table = self._native_solver.getAvgTable()
            # CallAvgTable not yet implemented
            call_avg_table = None
            return avg_table, call_avg_table
        return None, None

    getRawAvgTables = get_raw_avg_tables

    @staticmethod
    def defaultOptions():
        """Get default LQNS options."""
        return OptionsDict({
            'verbose': VerboseLevel.STD
        })

    @staticmethod
    def isAvailable():
        """Check if LQNS solver is available."""
        try:
            from .solvers.solver_lqns import SolverLQNS
            return SolverLQNS.isAvailable()
        except ImportError:
            return False

    default_options = defaultOptions
    is_available = isAvailable


class SolverLN(EnsembleSolver):
    """Layered Network solver."""

    def __init__(self, model, solver_factory=None, options=None, **kwargs):
        solver_options = SolverOptions(SolverType.LN)
        super().__init__(solver_options, **kwargs)
        self.model = model
        self._solver_factory = solver_factory

        from .solvers.solver_ln import SolverLN
        self._native_solver = SolverLN(model, **kwargs)

    def runAnalyzer(self):
        """Run the solver analysis."""
        if self._native_solver is not None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result
        return self

    run_analyzer = runAnalyzer

    def avg_table(self):
        """Get average performance metrics table."""
        if self._native_solver is not None:
            return self._native_solver.getAvgTable()
        return None

    def getAvgTable(self):
        """Get average performance metrics table."""
        return self.avg_table()

    @staticmethod
    def defaultOptions():
        """Get default LN options."""
        return OptionsDict({
            'verbose': VerboseLevel.STD,
            'method': 'exact',
            'tol': 1e-6,
            'max_iter': 100,
        })

    default_options = defaultOptions


class SolverPosterior(EnsembleSolver):
    """Posterior solver for Bayesian inference."""

    def __init__(self, model, *args, **kwargs):
        options = SolverOptions()
        super().__init__(options, *args, **kwargs)
        self.model = model

        from .solvers.solver_posterior import SolverPosterior
        self._native_solver = SolverPosterior(model, **kwargs)

    def runAnalyzer(self):
        """Run the solver analysis."""
        if self._native_solver is not None:
            self._native_solver.runAnalyzer()
            self._result = self._native_solver._result
        return self

    run_analyzer = runAnalyzer

    @staticmethod
    def defaultOptions():
        """Get default Posterior options."""
        return OptionsDict({
            'verbose': VerboseLevel.STD
        })

    default_options = defaultOptions


# Option classes for specific solvers
class CTMCOptions(SolverOptions):
    """CTMC-specific options."""
    def __init__(self):
        super().__init__(SolverType.CTMC)


class EnvOptions(SolverOptions):
    """ENV-specific options."""
    def __init__(self):
        super().__init__(SolverType.ENV)


class FluidOptions(SolverOptions):
    """Fluid-specific options."""
    def __init__(self):
        super().__init__(SolverType.FLUID)


class JMTOptions(SolverOptions):
    """JMT-specific options."""
    def __init__(self):
        super().__init__(SolverType.JMT)


class LNOptions(SolverOptions):
    """LN-specific options."""
    def __init__(self):
        super().__init__(SolverType.LN)


class LQNSOptions(SolverOptions):
    """LQNS-specific options."""
    def __init__(self):
        super().__init__(SolverType.LQNS)


class MAMOptions(SolverOptions):
    """MAM-specific options."""
    def __init__(self):
        super().__init__(SolverType.MAM)


class MVAOptions(SolverOptions):
    """MVA-specific options."""
    def __init__(self):
        super().__init__(SolverType.MVA)


class NCOptions(SolverOptions):
    """NC-specific options."""
    def __init__(self):
        super().__init__(SolverType.NC)


class QNSOptions(SolverOptions):
    """QNS-specific options."""
    def __init__(self):
        super().__init__(SolverType.QNS)


class SSAOptions(SolverOptions):
    """SSA-specific options."""
    def __init__(self):
        super().__init__(SolverType.SSA)


class AutoOptions(SolverOptions):
    """Auto-specific options."""
    def __init__(self):
        super().__init__(SolverType.AUTO)


class SolverConfig:
    """Solver configuration container."""

    def __init__(self):
        self.solver = None
        self.options = None


# Short aliases for solver classes (MATLAB/Java API compatibility)
MVA = SolverMVA
NC = SolverNC
CTMC = SolverCTMC
SSA = SolverSSA
FLD = SolverFluid
Fluid = SolverFluid
MAM = SolverMAM
JMT = SolverJMT
DES = SolverDES
AUTO = SolverAuto
QNS = SolverQNS
LQNS = SolverLQNS
LN = SolverLN
Posterior = SolverPosterior

# Export all public classes
__all__ = [
    'Solver',
    'EnsembleSolver',
    'NetworkSolver',
    'SolverMVA',
    'SolverNC',
    'SolverCTMC',
    'SolverSSA',
    'SolverFluid',
    'SolverMAM',
    'SolverJMT',
    'SolverDES',
    'SolverAuto',
    'SolverENV',
    'SolverQNS',
    'SolverLQNS',
    'SolverLN',
    'SolverPosterior',
    'LINE',
    'ENV',
    'SolverOptions',
    'SolverConfig',
    'CTMCOptions',
    'EnvOptions',
    'FluidOptions',
    'JMTOptions',
    'LNOptions',
    'LQNSOptions',
    'MAMOptions',
    'MVAOptions',
    'NCOptions',
    'QNSOptions',
    'SSAOptions',
    'AutoOptions',
    'SampleResult',
    'EventInfo',
    'DistributionResult',
    'ProbabilityResult',
    # Short aliases
    'MVA',
    'NC',
    'CTMC',
    'SSA',
    'FLD',
    'Fluid',
    'MAM',
    'JMT',
    'DES',
    'AUTO',
    'QNS',
    'LQNS',
    'LN',
    'Posterior',
]
