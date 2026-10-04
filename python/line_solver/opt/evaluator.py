"""
LINE Model Evaluator for Optimization

This module provides the interface between the optimization solver and
LINE's SolverAUTO for model evaluation.

Key Classes:
    - LineEvaluator: Wraps LINE model evaluation for optimization

Example:
    >>> evaluator = LineEvaluator(model, variables)
    >>> result = evaluator.evaluate(x)
    >>> print(f"Response time: {result.getResponseTime('Queue')}")
"""

import time
import copy
import logging
import numpy as np
from typing import List, Dict, Any, Optional

from .variables import DecisionVariable
from .results import EvaluationResult

logger = logging.getLogger(__name__)


class LineEvaluator:
    """
    Wraps LINE model evaluation for optimization.

    Applies decision variable values to a model copy, solves with SolverAUTO,
    and extracts performance metrics.

    This class handles the mapping between continuous optimization vectors
    and discrete LINE model parameters.

    Args:
        model: LINE Network model (base model that will be copied for each evaluation)
        variables: List of DecisionVariable objects defining the optimization space
        fixed_variables: Optional list of (DecisionVariable, value) pairs applied
            to every model copy before the free variables. Used by the
            decomposition workflow to coordinate subproblems.

    Example:
        >>> from line_solver import Network, Queue, Source, Sink
        >>> from line_solver import ServerAllocation, LineEvaluator
        >>>
        >>> model = Network("Test")
        >>> # ... build model ...
        >>> vars = [ServerAllocation(queue, bounds=(1, 10))]
        >>> evaluator = LineEvaluator(model, vars)
        >>> result = evaluator.evaluate(np.array([0.5]))
    """

    def __init__(self, model: 'Network', variables: List[DecisionVariable],
                 fixed_variables: Optional[List[tuple]] = None):
        """
        Initialize the evaluator.

        Args:
            model: LINE Network model
            variables: List of decision variables
            fixed_variables: Optional list of (DecisionVariable, value) pairs
                applied before the free variables on each evaluation
        """
        self._base_model = model
        self._variables = variables
        self._fixed_variables = list(fixed_variables) if fixed_variables else []
        self._evaluation_count = 0

        # LayeredNetwork (LQN) models take a distinct evaluation path: solved
        # with SolverLN and read from the per-LQN-node average table instead of
        # SolverAUTO and the per-station table (see opt/layered.py).
        from .layered import is_layered
        self._is_layered = is_layered(model)

        # Cache dimension info
        self._total_dimension = sum(v.getDimension() for v in variables)
        self._var_offsets = []
        offset = 0
        for v in variables:
            self._var_offsets.append(offset)
            offset += v.getDimension()

        # Try to import LINE
        self._line_available = False
        try:
            if self._is_layered:
                from line_solver import SolverLN  # noqa: F401
            else:
                from line_solver import SolverAUTO  # noqa: F401
            self._line_available = True
        except ImportError:
            pass

    def getModel(self) -> 'Network':
        """Get the base model."""
        return self._base_model

    def getVariables(self) -> List[DecisionVariable]:
        """Get the decision variables."""
        return self._variables

    def getFixedVariables(self) -> List[tuple]:
        """Get the fixed (variable, value) pairs applied on each evaluation."""
        return list(self._fixed_variables)

    def getTotalDimension(self) -> int:
        """Get total dimension of the optimization vector."""
        return self._total_dimension

    def getEvaluationCount(self) -> int:
        """Get number of evaluations performed."""
        return self._evaluation_count

    def getBounds(self) -> List[tuple]:
        """
        Get bounds for the full optimization vector.

        Returns:
            List of (lower, upper) tuples for each dimension
        """
        bounds = []
        for var in self._variables:
            bounds.extend(var.getBounds())
        return bounds

    def decodeVariables(self, x: np.ndarray) -> Dict[str, Any]:
        """
        Decode optimization vector to variable values.

        Args:
            x: Continuous optimization vector

        Returns:
            Dict mapping variable names to decoded values
        """
        values = {}
        for i, var in enumerate(self._variables):
            offset = self._var_offsets[i]
            dim = var.getDimension()
            var_x = x[offset:offset + dim]
            values[var.getName()] = var.decode(var_x)
        return values

    def applyVariables(self, model: 'Network', values: Dict[str, Any]) -> None:
        """
        Apply decoded variable values to a model.

        Args:
            model: LINE Network model to modify
            values: Dict mapping variable names to values
        """
        for var in self._variables:
            value = values.get(var.getName())
            if value is not None:
                var.apply(model, value)

    def copyModel(self) -> 'Network':
        """
        Create a copy of the base model for evaluation.

        Returns:
            Copy of the base model
        """
        # see _kb/05-solvers-overview.md (LineOpt: opt/evaluator.py) for rationale
        model = self._base_model
        copier = getattr(model, 'copy', None)
        if callable(copier):
            try:
                copied = copier()
                if copied is not None:
                    return copied
            except Exception:
                pass
        try:
            return copy.deepcopy(model)
        except Exception:
            return model

    def evaluate(self, x: np.ndarray) -> EvaluationResult:
        """
        Evaluate the model at decision vector x.

        Decodes x to variable values and delegates to evaluateValues().

        Args:
            x: Continuous optimization vector

        Returns:
            EvaluationResult with performance metrics
        """
        return self.evaluateValues(self.decodeVariables(x))

    def evaluateValues(self, values: Dict[str, Any]) -> EvaluationResult:
        """
        Evaluate the model at decoded variable values.

        This method:
        1. Creates a model copy
        2. Applies fixed variables, then the given values
        3. Solves with SolverAUTO
        4. Extracts per-station and system performance metrics

        Args:
            values: Dict mapping variable names to decoded values

        Returns:
            EvaluationResult with performance metrics
        """
        self._evaluation_count += 1
        start_time = time.time()

        result = EvaluationResult()

        if not self._line_available:
            # Return mock result if LINE not available
            result.feasible = False
            result.solve_time = time.time() - start_time
            return result

        try:
            # Copy and modify model: fixed variables first, then free values
            model = self.copyModel()
            for var, fixed_value in self._fixed_variables:
                var.apply(model, fixed_value)
            self.applyVariables(model, values)

            if self._is_layered:
                # LQN path: solve with SolverLN, read the per-node LQN table.
                from .layered import solve_lqn_avg
                _solver, df = solve_lqn_avg(model)
                result.feasible = True
                result.solver_used = 'SolverLN'
                self._extractLayeredMetrics(df, result)
                self._extractLayeredSystemMetrics(model, result)
                # LQN sensitivities are expensive (they re-solve every layer),
                # so they are computed lazily by the gradient path only, not on
                # every evaluation. See LineOptSolver._lqnAnalyticGradient.
            else:
                # Flat path: SolverAUTO + per-station table.
                from line_solver import SolverAUTO
                solver = SolverAUTO(model)
                avg_table = solver.getAvgTable()

                result.feasible = True
                result.solver_used = solver.getSelectedSolverName()

                # Extract metrics from the result table
                self._extractMetrics(avg_table, result)
                self._extractSystemMetrics(solver, model, result)

                # Attach analytic performance sensitivities when the model is a
                # supported product-form network (None otherwise; the optimizer
                # then falls back to finite differences).
                from .sensitivity import compute_model_sensitivities
                result.sensitivities = compute_model_sensitivities(model)

        except Exception as e:
            # Model evaluation failed
            logger.warning("LINE evaluation failed: %s", e)
            result.feasible = False
            result.solve_time = time.time() - start_time
            return result

        result.solve_time = time.time() - start_time
        return result

    # Column-name candidates in LINE's average metric table (AvgTable).
    _STATION_COLS = ('Station', 'Node')
    _CLASS_COLS = ('JobClass', 'Class', 'Chain')
    _RESPT_COLS = ('RespT',)
    _TPUT_COLS = ('Tput',)
    _UTIL_COLS = ('Util',)
    _QLEN_COLS = ('QLen',)

    @staticmethod
    def _pickColumn(columns, candidates):
        """Return the first candidate present in columns, or None."""
        for name in candidates:
            if name in columns:
                return name
        return None

    def _extractMetrics(self, avg_table: Any,
                        result: EvaluationResult) -> None:
        """
        Extract performance metrics from LINE's average metric table.

        Parses the DataFrame returned by ``SolverAUTO.getAvgTable()``, which is
        indexed by Station/JobClass name (columns RespT, Tput, Util, QLen). This
        is robust to solver internals, unlike positional index access.

        Args:
            avg_table: IndexedTable or pandas DataFrame from getAvgTable()
            result: EvaluationResult to populate
        """
        # IndexedTable wraps the DataFrame in a .data attribute.
        df = getattr(avg_table, 'data', avg_table)
        if df is None or not hasattr(df, 'columns'):
            logger.warning("getAvgTable() returned no parseable table (%r); "
                           "no metrics extracted", type(avg_table).__name__)
            return

        columns = list(df.columns)
        station_col = self._pickColumn(columns, self._STATION_COLS)
        class_col = self._pickColumn(columns, self._CLASS_COLS)
        respt_col = self._pickColumn(columns, self._RESPT_COLS)
        tput_col = self._pickColumn(columns, self._TPUT_COLS)
        util_col = self._pickColumn(columns, self._UTIL_COLS)
        qlen_col = self._pickColumn(columns, self._QLEN_COLS)

        if station_col is None or class_col is None:
            logger.warning("AvgTable missing Station/JobClass columns "
                           "(have %s); no metrics extracted", columns)
            return

        for _, row in df.iterrows():
            station = str(row[station_col])
            jobclass = str(row[class_col])
            key = (station, jobclass)

            if respt_col is not None:
                result.response_times[key] = float(row[respt_col])
            if tput_col is not None:
                result.throughputs[key] = float(row[tput_col])
            if qlen_col is not None:
                result.queue_lengths[key] = float(row[qlen_col])
            if util_col is not None:
                # Utilization is per-station: aggregate across classes.
                result.utilizations[station] = (
                    result.utilizations.get(station, 0.0) + float(row[util_col])
                )

    def _extractSystemMetrics(self, solver: Any, model: 'Network',
                              result: EvaluationResult) -> None:
        """
        Extract end-to-end (system) metrics from getAvgSysTable().

        The system table has one row per chain (in chain-index order) with
        columns SysRespT and SysTput, labeled with generic chain names.
        Rows are mapped back to job class names via the network struct's
        chain membership matrix, so metrics are retrievable by class name.
        For open chains the response time is recomputed by Little's law
        (see inline comment); for closed chains the table value is kept.
        Absence of the table is logged, not fatal.

        Args:
            solver: SolverAUTO instance after solving
            model: The solved model (for the chain-to-class mapping)
            result: EvaluationResult to populate
        """
        try:
            sys_table = solver.getAvgSysTable()
        except Exception as e:
            logger.warning("getAvgSysTable() failed: %s; system metrics "
                           "unavailable", e)
            return

        df = getattr(sys_table, 'data', sys_table)
        if df is None or not hasattr(df, 'columns'):
            return

        columns = list(df.columns)
        respt_col = self._pickColumn(columns, ('SysRespT',))
        tput_col = self._pickColumn(columns, ('SysTput',))

        # Map chain index -> member class names via the network struct,
        # and record whether the chain is open (infinite population)
        chain_classes = []
        chain_is_open = []
        try:
            sn = model.getStruct()
            chains = np.asarray(sn.chains)
            classnames = list(sn.classnames)
            njobs = np.asarray(sn.njobs).flatten()
            for ci in range(chains.shape[0]):
                members = [k for k in range(chains.shape[1])
                           if chains[ci, k] > 0]
                chain_classes.append([classnames[k] for k in members])
                chain_is_open.append(any(np.isinf(njobs[k]) for k in members))
        except Exception as e:
            logger.warning("Chain-to-class mapping unavailable (%s); "
                           "system metrics keyed by row label", e)

        class_col = self._pickColumn(columns, ('Class', 'JobClass', 'Chain'))

        for position, (_, row) in enumerate(df.iterrows()):
            # Key by member class names when the mapping is available,
            # falling back to the table's own row label
            if position < len(chain_classes) and chain_classes[position]:
                keys = chain_classes[position]
            elif class_col is not None:
                keys = [str(row[class_col])]
            else:
                continue

            tput = float(row[tput_col]) if tput_col is not None else 0.0
            respt = float(row[respt_col]) if respt_col is not None else None

            # see _kb/05-solvers-overview.md (LineOpt: opt/evaluator.py) for rationale
            if (position < len(chain_is_open) and chain_is_open[position]
                    and tput > 0 and result.queue_lengths):
                jobs_in_system = sum(
                    qlen for (_, cls), qlen in result.queue_lengths.items()
                    if cls in keys)
                respt = jobs_in_system / tput

            for key in keys:
                if respt is not None:
                    result.system_response_times[key] = respt
                if tput_col is not None:
                    result.system_throughputs[key] = tput

    # ---- LayeredNetwork (LQN) metric extraction -------------------------

    def _extractLayeredMetrics(self, df: Any,
                               result: EvaluationResult) -> None:
        """Extract per-LQN-node metrics from SolverLN's average table.

        The LQN table has one row per node (Processor/Task/Entry/Activity) with
        columns Node, NodeType, QLen, Util, RespT, ResidT, ArvR, Tput. Metrics
        are keyed by node so objectives/constraints resolve by LQN node name:
        throughput/queue-length/response-time under (node, node) and
        utilization under node, matching the flat EvaluationResult convention
        (utilization keyed by station, the rest by (station, class)). Non-finite
        cells (e.g. a processor's NaN response time) are skipped.
        """
        if df is None or not hasattr(df, 'columns'):
            logger.warning("SolverLN.get_avg_table() returned no parseable "
                           "table; no LQN metrics extracted")
            return
        columns = list(df.columns)
        node_col = self._pickColumn(columns, ('Node', 'Station'))
        if node_col is None:
            logger.warning("LQN table missing Node column (have %s)", columns)
            return
        respt_col = self._pickColumn(columns, self._RESPT_COLS)
        tput_col = self._pickColumn(columns, self._TPUT_COLS)
        util_col = self._pickColumn(columns, self._UTIL_COLS)
        qlen_col = self._pickColumn(columns, self._QLEN_COLS)

        def finite(v):
            try:
                f = float(v)
            except (TypeError, ValueError):
                return None
            return f if np.isfinite(f) else None

        for _, row in df.iterrows():
            node = str(row[node_col])
            key = (node, node)
            if respt_col is not None:
                v = finite(row[respt_col])
                if v is not None:
                    result.response_times[key] = v
            if tput_col is not None:
                v = finite(row[tput_col])
                if v is not None:
                    result.throughputs[key] = v
            if qlen_col is not None:
                v = finite(row[qlen_col])
                if v is not None:
                    result.queue_lengths[key] = v
            if util_col is not None:
                v = finite(row[util_col])
                if v is not None:
                    result.utilizations[node] = v

    def _extractLayeredSystemMetrics(self, model: Any,
                                     result: EvaluationResult) -> None:
        """Derive end-to-end (system) metrics from the reference task.

        A closed LQN's system throughput is the reference task's throughput;
        its end-to-end response time is the sum of response times over the
        reference task's entries (one entry -> the entry's RespT). Keyed by the
        reference task's name so MinimizeSystemResponseTime /
        SystemResponseTimeConstraint resolve without a chain concept.
        """
        from .layered import elem_name
        ref_tasks = []
        for task in list(getattr(model, 'tasks', [])):
            sched = getattr(task, 'sched_strategy', None)
            sched_name = getattr(sched, 'name', str(sched))
            is_ref = getattr(task, 'is_reference', False) \
                or str(sched_name).upper() in ('REF', 'REFERENCE')
            if is_ref:
                ref_tasks.append(task)
        for task in ref_tasks:
            tname = elem_name(task)
            tput = result.throughputs.get((tname, tname), 0.0)
            if tput > 0:
                result.system_throughputs[tname] = tput
            # Sum response times over this task's entries.
            entries = list(getattr(task, 'entries', []))
            total_rt = 0.0
            have_rt = False
            for entry in entries:
                ename = elem_name(entry)
                rt = result.response_times.get((ename, ename))
                if rt is not None:
                    total_rt += rt
                    have_rt = True
            if have_rt:
                result.system_response_times[tname] = total_rt

    def evaluateLayeredSensitivities(self, values: Dict[str, Any]
                                     ) -> Optional[Dict]:
        """Compute LQN per-layer service-rate partial sensitivities on demand.

        Rebuilds the configured model copy, solves it with SolverLN, and parses
        SolverLN.getSensitivityTable into a
        ``{(station, jobclass): {metric: d/dRate}}`` dict. Called only by the
        partial-sensitivity gradient path, never on the hot evaluation loop.
        Returns None on any failure (the caller then finite-differences).
        """
        if not self._is_layered or not self._line_available:
            return None
        try:
            from .layered import make_lqn_solver, compute_lqn_sensitivities
            model = self.copyModel()
            for var, fixed_value in self._fixed_variables:
                var.apply(model, fixed_value)
            self.applyVariables(model, values)
            solver = make_lqn_solver(model)
            return compute_lqn_sensitivities(solver)
        except Exception as e:
            logger.debug("LQN sensitivity evaluation failed: %s", e)
            return None

    evaluate_layered_sensitivities = evaluateLayeredSensitivities

    @property
    def is_layered(self) -> bool:
        """True if this evaluator drives a LayeredNetwork (LQN) model."""
        return self._is_layered

    @staticmethod
    def _valuesKey(values: Dict[str, Any]) -> tuple:
        """
        Build a hashable cache key from decoded variable values.

        Keying the cache on the decoded configuration (rather than the raw
        continuous vector) collapses the rounding plateaus of integer
        variables, so distinct DE candidates that map to the same discrete
        configuration trigger a single LINE solve.
        """
        def canonical(v):
            if isinstance(v, np.ndarray):
                return tuple(np.round(v, 9).tolist())
            if isinstance(v, (list, tuple)):
                return tuple(canonical(item) for item in v)
            if isinstance(v, float):
                return round(v, 9)
            return v

        return tuple(sorted((name, canonical(v)) for name, v in values.items()))

    def evaluateWithCache(self, x: np.ndarray,
                          cache: Dict[tuple, EvaluationResult]) -> EvaluationResult:
        """
        Evaluate with result caching.

        Args:
            x: Optimization vector
            cache: Dict to cache results (keyed by decoded configuration)

        Returns:
            EvaluationResult (from cache if available)
        """
        return self.evaluateValuesWithCache(self.decodeVariables(x), cache)

    def evaluateValuesWithCache(self, values: Dict[str, Any],
                                cache: Dict[tuple, EvaluationResult]
                                ) -> EvaluationResult:
        """
        Evaluate decoded values with result caching.

        Args:
            values: Dict mapping variable names to decoded values
            cache: Dict to cache results (keyed by decoded configuration)

        Returns:
            EvaluationResult (from cache if available)
        """
        key = self._valuesKey(values)
        if key in cache:
            return cache[key]

        result = self.evaluateValues(values)
        cache[key] = result
        return result

    # Snake case aliases
    get_model = getModel
    get_variables = getVariables
    get_fixed_variables = getFixedVariables
    get_total_dimension = getTotalDimension
    get_evaluation_count = getEvaluationCount
    get_bounds = getBounds
    decode_variables = decodeVariables
    apply_variables = applyVariables
    copy_model = copyModel
    evaluate_values = evaluateValues
    evaluate_with_cache = evaluateWithCache
    evaluate_values_with_cache = evaluateValuesWithCache


# Public API
__all__ = ['LineEvaluator']
