"""
LayeredNetwork (LQN) support for line-opt.

line-opt's core (problem/evaluator/variables/solver) was written for the flat
``Network``: it copies the model, solves it with ``SolverAuto``, and reads a
per-(station, class) average table. A ``LayeredNetwork`` cannot go through that
path -- it is an ``Ensemble`` of per-layer ``Network`` submodels, solved by
``SolverLN`` into a per-LQN-node average table (Node, NodeType, QLen, Util,
RespT, ResidT, ArvR, Tput).

This module is the adapter that lets the same optimizer drive an LQN. It
provides:

* ``is_layered`` -- model-type detection used to branch problem/evaluator.
* element resolution by name in a per-evaluation model copy (processors, tasks,
  activities), mirroring ``variables._resolveNode`` for the flat case.
* ``solve_lqn_avg`` -- run ``SolverLN`` and return the parsed LQN average table.
* ``compute_lqn_sensitivities`` -- run ``SolverLN.getSensitivityTable`` and
  reshape the per-(Layer, Station, JobClass) service-rate partial derivatives
  into a ``{(station, jobclass): {metric: d/dRate}}`` dict for the gradient.

IMPORTANT (see ``SolverLN.getSensitivityTable`` docstring): the per-layer table
holds WITHIN-LAYER PARTIAL derivatives, taken with the fixed-point layer
parameters held constant. It omits the cross-layer coupling term, so it is a
biased estimate of the total derivative of the solved layered model. The
optimizer's ``lqn_gradient='fd'`` mode finite-differences the whole
``LayeredNetwork`` instead (correct total derivative); ``partial_sens`` uses
this table directly (cheap, biased); ``partial_plus_fd`` corrects it with a
periodic full-model finite difference.
"""

import logging
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


def is_layered(model: Any) -> bool:
    """True if ``model`` is a LayeredNetwork (LQN), false for a flat Network.

    Prefers an ``isinstance`` check against the native ``LayeredNetwork`` class
    and falls back to duck typing on the LQN element containers so the check
    also holds for backend-agnostic wrappers that expose the same surface.
    """
    if model is None:
        return False
    try:
        from line_solver import LayeredNetwork
        if isinstance(model, LayeredNetwork):
            return True
    except Exception:
        pass
    return (hasattr(model, 'processors') and hasattr(model, 'tasks')
            and hasattr(model, 'activities') and hasattr(model, 'getEnsemble'))


def elem_name(element: Any) -> str:
    """Name of an LQN element (Processor/Task/Entry/Activity).

    LQN elements carry a plain ``name`` attribute rather than the ``getName()``
    method that flat nodes use; support both so callers need not care.
    """
    name = getattr(element, 'name', None)
    if name is not None:
        return str(name)
    getter = getattr(element, 'getName', None)
    if callable(getter):
        return str(getter())
    return str(element)


def dist_mean(value: Any) -> Optional[float]:
    """Mean of a think-time/demand value that may be a distribution or scalar.

    Handles native distributions (``getMean``/``get_mean``), plain numbers, and
    None (returns None). Used by LQN variables to read a model's current
    parameter for layer freezing.
    """
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    for attr in ('getMean', 'get_mean'):
        getter = getattr(value, attr, None)
        if callable(getter):
            try:
                return float(getter())
            except Exception:
                return None
    return None


def _by_name(elements: List[Any], name: str) -> Optional[Any]:
    for element in elements:
        if elem_name(element) == name:
            return element
    return None


def resolve_processor(model: Any, name: str) -> Optional[Any]:
    """Resolve a Processor by name inside a (possibly copied) LQN model."""
    return _by_name(list(getattr(model, 'processors', [])), name)


def resolve_task(model: Any, name: str) -> Optional[Any]:
    """Resolve a Task by name inside a (possibly copied) LQN model."""
    return _by_name(list(getattr(model, 'tasks', [])), name)


def resolve_activity(model: Any, name: str) -> Optional[Any]:
    """Resolve an Activity by name inside a (possibly copied) LQN model."""
    return _by_name(list(getattr(model, 'activities', [])), name)


def activity_processor_name(model: Any, activity_name: str) -> Optional[str]:
    """Name of the processor an activity ultimately runs on.

    Activity -> parent Task -> deployed Processor. Returns None if the chain is
    incomplete (e.g. an activity not yet bound to a task). Used both to tag a
    HostDemand variable with the host layer it perturbs and to key its row in
    the per-layer sensitivity table, whose host-layer rows are
    (Layer=processor, Station=processor, JobClass=activity).
    """
    act = resolve_activity(model, activity_name)
    if act is None:
        return None
    task = getattr(act, 'task', None) or (
        act.getParent() if hasattr(act, 'getParent') else None)
    if task is None:
        return None
    proc = getattr(task, 'processor', None)
    if proc is None:
        return None
    return elem_name(proc)


def task_processor_name(model: Any, task_name: str) -> Optional[str]:
    """Name of the processor a task is deployed on (None if undeployed)."""
    task = resolve_task(model, task_name)
    if task is None:
        return None
    proc = getattr(task, 'processor', None)
    return elem_name(proc) if proc is not None else None


def make_lqn_solver(model: Any, options: Optional[dict] = None):
    """Construct a quiet ``SolverLN`` for a per-evaluation LQN model copy.

    Silences the layer/avg-table printing that ``SolverLN`` does by default so
    an optimization run of thousands of solves does not flood stdout.
    """
    from line_solver import SolverLN
    opts = dict(options or {})
    opts.setdefault('verbose', False)
    solver = SolverLN(model, **opts)
    # Suppress the get_avg_table print (mirrors the MATLAB _table_silent path).
    try:
        solver._table_silent = True
    except Exception:
        pass
    return solver


def solve_lqn_avg(model: Any, options: Optional[dict] = None):
    """Solve an LQN and return (solver, avg_table_DataFrame).

    The average table has one row per LQN node (processor/task/entry/activity)
    with columns Node, NodeType, QLen, Util, RespT, ResidT, ArvR, Tput.
    """
    solver = make_lqn_solver(model, options)
    table = solver.get_avg_table()
    df = getattr(table, 'data', table)
    return solver, df


# Sensitivity-table metric column -> the EvaluationResult metric kind it feeds.
_SENS_COLUMNS = {
    'Tput': 'dTput_dRate',
    'RespT': 'dRespT_dRate',
    'QLen': 'dQLen_dRate',
    'Util': 'dUtil_dRate',
}


def compute_lqn_sensitivities(solver: Any) -> Optional[Dict]:
    """Per-(Station, JobClass) within-layer service-rate partial derivatives.

    Runs ``SolverLN.getSensitivityTable`` on an already-built solver (kept at
    its converged fixed point) and reshapes the table into
    ``{(station, jobclass): {'Tput':.., 'RespT':.., 'QLen':.., 'Util':..}}``,
    each value being the derivative of that layer row's mean measure with
    respect to that station-class service RATE.

    Returns None if the table cannot be produced. See the module docstring for
    the partial-vs-total caveat that makes this a biased gradient source.
    """
    try:
        table = solver.getSensitivityTable()
    except Exception as e:
        logger.debug("LQN sensitivity table unavailable: %s", e)
        return None

    df = getattr(table, 'data', table)
    if df is None or not hasattr(df, 'columns'):
        return None

    out: Dict[tuple, Dict[str, float]] = {}
    for _, row in df.iterrows():
        key = (str(row['Station']), str(row['JobClass']))
        entry = out.setdefault(key, {})
        for kind, col in _SENS_COLUMNS.items():
            try:
                entry[kind] = float(row[col])
            except Exception:
                entry[kind] = 0.0
    return out
