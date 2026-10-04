"""
Native Python implementation of Random Environment models.

This module provides classes for defining and analyzing queueing networks
in random environments, where the network parameters change according to
an underlying Markov modulated process.

Implements full parity with MATLAB SolverENV using transient analysis
with iteration until convergence.
"""

import os as _os
import warnings

import numpy as np
from typing import Optional, List, Dict, Any, Union, Callable, Tuple
from .constants import VerboseLevel
import pandas as pd


def _snake_camel_getattr(self, name):
    """Fallback resolving documented snake_case names to the camelCase
    attributes/methods implemented on the class (e.g. prob_env -> probEnv,
    hold_time -> holdTime, get_solver -> getSolver). Only invoked when normal
    lookup fails; private/dunder names are left untouched so genuine
    ``AttributeError`` and copy/pickle behaviour are preserved.
    """
    if name.startswith('_') or '_' not in name:
        raise AttributeError(name)
    parts = name.split('_')
    camel = parts[0] + ''.join(p[:1].upper() + p[1:] for p in parts[1:])
    if camel != name:
        inst = object.__getattribute__(self, '__dict__')
        if camel in inst:
            return inst[camel]
        if hasattr(type(self), camel):
            return getattr(self, camel)
    raise AttributeError(name)


def _get_rate(dist) -> float:
    """Extract rate from a distribution object."""
    if hasattr(dist, 'getRate'):
        return dist.getRate()
    elif hasattr(dist, 'rate'):
        return dist.rate
    elif hasattr(dist, 'getMean'):
        mean = dist.getMean()
        return 1.0 / mean if mean > 0 else 0.0
    else:
        return float(dist)


def _get_map_representation(dist) -> Tuple[np.ndarray, np.ndarray]:
    """Get MAP {D0, D1} representation of a distribution.

    For phase-type distributions, returns the sub-generator matrix D0 and
    completion matrix D1. For exponential distributions, returns scalar matrices.
    """
    if hasattr(dist, 'getD0') and hasattr(dist, 'getD1'):
        D0 = np.atleast_2d(dist.getD0())
        D1 = np.atleast_2d(dist.getD1())
        return D0, D1
    elif hasattr(dist, 'T') and hasattr(dist, 'alpha'):
        # Phase-type distribution with T matrix and alpha vector
        T = np.atleast_2d(dist.T)
        alpha = np.atleast_1d(dist.alpha)
        t = -T.sum(axis=1)  # Exit rate vector
        D0 = T
        D1 = np.outer(t, alpha)
        return D0, D1
    else:
        # Fall back to exponential approximation
        rate = _get_rate(dist)
        D0 = np.array([[-rate]])
        D1 = np.array([[rate]])
        return D0, D1


def _krons(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Kronecker sum of matrices A and B.

    S = A ⊗ I_B + I_A ⊗ B
    """
    return np.kron(A, np.eye(B.shape[0])) + np.kron(np.eye(A.shape[0]), B)


def _map_prob(D0: np.ndarray, D1: np.ndarray) -> np.ndarray:
    """Compute stationary probability vector of a MAP."""
    Q = D0 + D1
    n = Q.shape[0]

    # Solve pi * Q = 0, sum(pi) = 1
    A = Q.T.copy()
    A[-1, :] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0

    try:
        pi = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        pi = np.linalg.lstsq(A, b, rcond=None)[0]

    pi = np.maximum(pi, 0)
    if pi.sum() > 0:
        pi = pi / pi.sum()
    return pi


def _map_pie(D0: np.ndarray, D1: np.ndarray) -> np.ndarray:
    """Compute equilibrium distribution of the embedded DTMC at departure epochs.

    For a MAP {D0, D1}, the embedded chain has transition matrix P = (-D0)^{-1} D1.
    The stationary vector of P is: pie = pi * D1 / (pi * D1 * e),
    where pi is the CTMC stationary vector of D0 + D1.

    This is the correct initial probability vector for CDF computation.
    Matches MATLAB's map_pie(MAP).
    """
    pi = _map_prob(D0, D1)
    A = pi @ D1
    s = A.sum()
    if s > 0:
        return A / s
    return pi  # fallback


def _map_lambda(D0: np.ndarray, D1: np.ndarray) -> float:
    """Compute arrival rate of a MAP."""
    pi = _map_prob(D0, D1)
    e = np.ones(D1.shape[1])
    return float(pi @ D1 @ e)


def _map_mean(D0: np.ndarray, D1: np.ndarray) -> float:
    """Compute mean inter-arrival time of a MAP."""
    lam = _map_lambda(D0, D1)
    return 1.0 / lam if lam > 0 else float('inf')


def _map_eval_cdf(D0: np.ndarray, D1: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Evaluate CDF of a MAP at time points t.

    For a MAP {D0, D1}, the CDF is:
    P(X <= t) = 1 - pie * exp(D0 * t) * e

    where pie is the equilibrium distribution of the embedded DTMC
    at departure epochs (matches MATLAB's map_pie).
    """
    from scipy.linalg import expm

    t = np.atleast_1d(t)
    n = D0.shape[0]
    e = np.ones(n)

    # Get initial probability vector (embedded DTMC at departure epochs)
    pie = _map_pie(D0, D1)

    cdf = np.zeros(len(t))
    for i, ti in enumerate(t):
        if ti <= 0:
            cdf[i] = 0.0
        else:
            exp_D0_t = expm(D0 * ti)
            cdf[i] = 1.0 - float(pie @ exp_D0_t @ e)

    return np.clip(cdf, 0.0, 1.0)


def _mmap_normalize(mmap: List[np.ndarray]) -> List[np.ndarray]:
    """Normalize a Marked MAP to ensure feasibility.

    MMAP format: [D0, D1, D1_class1, D1_class2, ...]
    where D1 = sum of all D1_class_k matrices.
    """
    if len(mmap) < 3:
        return mmap

    K = mmap[0].shape[0]  # Number of phases
    C = len(mmap) - 2     # Number of classes

    # Ensure non-negative off-diagonal elements in D0
    D0 = mmap[0].copy()
    for i in range(K):
        for j in range(K):
            if i != j:
                D0[i, j] = max(D0[i, j], 0)

    # Ensure non-negative elements in class-specific D1 matrices
    D1 = np.zeros_like(D0)
    for c in range(C):
        mmap[2 + c] = np.maximum(mmap[2 + c], 0)
        D1 += mmap[2 + c]
    mmap[1] = D1

    # Adjust diagonal of D0 so rows sum to zero
    for k in range(K):
        D0[k, k] = 0
        D0[k, k] = -np.sum(D0[k, :]) - np.sum(D1[k, :])
    mmap[0] = D0

    return mmap


def _mmap_count_lambda(mmap: List[np.ndarray]) -> np.ndarray:
    """Compute arrival rates for each class in a Marked MAP.

    Returns array of rates, one per class.
    """
    D0 = mmap[0]
    D1 = mmap[1]
    K = len(mmap) - 2  # Number of classes

    theta = _map_prob(D0, D1)
    e = np.ones(D0.shape[0])

    lk = np.zeros(K)
    for k in range(K):
        lk[k] = theta @ mmap[2 + k] @ e

    return lk


def _eval_cdf(dist, t) -> np.ndarray:
    """Evaluate CDF of distribution at time points t."""
    t = np.atleast_1d(t)
    if hasattr(dist, 'evalCDF'):
        return np.array([dist.evalCDF(ti) for ti in t])
    elif hasattr(dist, 'cdf'):
        return dist.cdf(t)
    else:
        # Assume exponential with rate
        rate = _get_rate(dist)
        if rate <= 0:
            return np.zeros_like(t)
        return 1.0 - np.exp(-rate * t)


def _finite_span(ts):
    """[t0, t1] as floats when both are finite and t1 > t0, else None."""
    if ts is None:
        return None
    try:
        t0 = float(ts[0])
        t1 = float(ts[1])
    except (TypeError, ValueError, IndexError, KeyError):
        return None
    if not (np.isfinite(t0) and np.isfinite(t1) and t1 > t0):
        return None
    return [t0, t1]


def _cdf_refine_grid(t_vals, mean_sojourn, n_interp):
    """Build the refinement grid the CDF-weighted exit average is summed on.

    The grid holds 90% of its points under 5*E[S], where the holding-time CDF
    carries essentially all of its mass, and spreads the rest across the tail
    up to the horizon. It is REBUILT UNCONDITIONALLY: keeping the integrator's
    own grid when it happened to be dense stops wherever its steps fell rather
    than where the sum converged, and it makes MATLAB, the JAR, C++ and this
    port sum different points. See _interpolate_for_cdf_map for the measured
    convergence sequence.
    """
    t_cdf_end = min(t_vals[-1], 5.0 * mean_sojourn)
    if t_cdf_end <= t_vals[0]:
        t_cdf_end = t_vals[-1]

    n_dense = int(0.9 * n_interp)
    n_tail = n_interp - n_dense
    t_dense = np.linspace(t_vals[0], t_cdf_end, n_dense)
    if t_cdf_end < t_vals[-1] and n_tail > 1:
        t_tail = np.linspace(t_cdf_end, t_vals[-1], n_tail + 1)[1:]  # Exclude overlap
        return np.concatenate([t_dense, t_tail])
    return t_dense


def _resample_on(t_fine, t_vals, q_metric, u_tran_ir, t_tran_ir):
    """Linearly resample the queue, utilization and throughput onto t_fine."""
    q_fine = np.interp(t_fine, t_vals, q_metric)

    _, u_metric = _get_tran_data(u_tran_ir)
    u_fine = np.interp(t_fine, t_vals, u_metric) if u_metric is not None else None

    _, t_metric = _get_tran_data(t_tran_ir)
    t_fine_data = np.interp(t_fine, t_vals, t_metric) if t_metric is not None else None

    return t_fine, q_fine, u_fine, t_fine_data


def _interpolate_for_cdf(t_vals, q_metric, u_tran_ir, t_tran_ir, dist,
                         n_interp=5000):
    """Interpolate transient data onto a finer time grid for CDF weighting.

    The ODE grid is chosen for the HORIZON, not for the holding time, so the
    CDF-weighted average reads whatever points the integrator left inside the
    sojourn support. This resamples the metrics onto a grid concentrated where
    the CDF has mass, using the distribution mean as the time scale.

    Args:
        t_vals: Original time points from ODE solver
        q_metric: Queue length metric values at t_vals
        u_tran_ir: Utilization transient result (TranResult or dict)
        t_tran_ir: Throughput transient result (TranResult or dict)
        dist: Transition distribution (for determining CDF time scale)
        n_interp: Number of interpolation points to create

    Returns:
        Tuple of (t_fine, q_fine, u_fine, t_fine_data) where each is an array
        on the finer grid, or None for u/t if not available.
    """
    rate = _get_rate(dist)
    if rate > 0:
        mean_sojourn = 1.0 / rate
    else:
        mean_sojourn = (t_vals[-1] - t_vals[0]) / 10.0

    t_fine = _cdf_refine_grid(t_vals, mean_sojourn, n_interp)
    return _resample_on(t_fine, t_vals, q_metric, u_tran_ir, t_tran_ir)


def _interpolate_for_cdf_map(t_vals, q_metric, u_tran_ir, t_tran_ir, D0, D1,
                             n_interp=5000):
    """Interpolate transient data onto a finer time grid for MAP CDF weighting.

    Same as _interpolate_for_cdf but uses MAP {D0, D1} matrices to determine
    the mean sojourn time instead of a distribution object.

    Args:
        t_vals: Original time points from ODE solver
        q_metric: Queue length metric values at t_vals
        u_tran_ir: Utilization transient result
        t_tran_ir: Throughput transient result
        D0: MAP sub-generator matrix
        D1: MAP completion matrix
        n_interp: Number of interpolation points to create

    Returns:
        Tuple of (t_fine, q_fine, u_fine, t_fine_data)
    """
    # Compute mean sojourn from MAP: mean = 1 / (alpha @ D1 @ e)
    alpha = _map_prob(D0, D1)
    e = np.ones(D0.shape[0])
    total_rate = float(alpha @ D1 @ e)
    if total_rate > 0:
        mean_sojourn = 1.0 / total_rate
    else:
        mean_sojourn = (t_vals[-1] - t_vals[0]) / 10.0

    # THE GRID IS REBUILT UNCONDITIONALLY as of 2026-08-11. The former rule kept
    # the solver's own grid when at least 50 of its points fell inside the
    # sojourn support, which stops wherever the integrator's steps happened to
    # fall rather than where the sum has converged: measured on
    # renv_node_breakdown with every stage refined, the exit average runs
    # 0.462260, 0.460580, 0.459704, 0.459272, 0.459138, 0.459122 as n_interp
    # goes 500 to 5e4. Rebuilding always is also what makes MATLAB, the JAR, C++
    # and this port sum the SAME points; each then differs only by its own
    # trajectory. The error is first order because the sum reads the RIGHT
    # endpoint and an unstable stage's queue grows linearly in t, so the mass
    # beyond 5*E[S] multiplies a large metric.

    t_fine = _cdf_refine_grid(t_vals, mean_sojourn, n_interp)
    return _resample_on(t_fine, t_vals, q_metric, u_tran_ir, t_tran_ir)


def _interp_cov(qcov, tq, n):
    """A stage's covariance trajectory read at the instants TQ.

    Clamped to its own grid rather than extrapolated: linear extrapolation of a
    covariance can leave the positive semidefinite cone, while a convex
    combination of two members stays inside it. None when the stage carries no
    covariance, which is what a stage solver with a first moment only leaves.
    """
    if not isinstance(qcov, dict):
        return None
    tv = qcov.get('t')
    C = qcov.get('C')
    if tv is None or C is None:
        return None
    tv = np.asarray(tv, dtype=float).ravel()
    C = np.asarray(C, dtype=float)
    if tv.size == 0 or C.ndim != 3 or C.shape[0] != n or C.shape[1] != n \
            or C.shape[2] != tv.size:
        return None
    tq = np.asarray(tq, dtype=float).ravel()
    if tv.size == 1:
        return np.repeat(C[:, :, :1], tq.size, axis=2)
    tc = np.clip(tq, tv[0], tv[-1])
    flat = C.reshape(n * n, tv.size)
    out = np.empty((n * n, tc.size))
    for j in range(n * n):
        out[j, :] = np.interp(tc, tv, flat[j, :])
    return out.reshape(n, n, tc.size)


def _reset_jacobian(reset_fn, Qexit, M, K):
    """The Jacobian of a reset policy at the exit mean, so a covariance can
    cross a switch as R C R'.

    A reset policy is an arbitrary map q -> q on the (station, class) mean queue
    lengths and there is no general way to push a second moment through one. The
    delta method is the first-order image, which is the order the whole
    mean-field coupling works to; the two NAMED policies are linear and R is then
    EXACT, the identity for 'keep' and zero for 'clear'. Column-major flattening
    (ir = r*M + i) matches the covariance index space.
    """
    n = M * K
    R = np.zeros((n, n))
    Qexit = np.asarray(Qexit, dtype=float)
    base = np.asarray(reset_fn(Qexit), dtype=float)
    if base.size != n:
        return R
    base = base.ravel(order='F')
    step = 1e-6 * max(1.0, float(np.max(np.abs(Qexit))) if Qexit.size else 1.0)
    for j in range(n):
        Qp = Qexit.copy()
        Qp[j % M, j // M] += step
        pert = np.asarray(reset_fn(Qp), dtype=float).ravel(order='F')
        R[:, j] = (pert - base) / step
    return R


def _set_stage_config(solver, key, value):
    """Put one entry on a stage solver's options.config, creating it if absent."""
    opts = getattr(solver, 'options', None)
    if opts is None:
        return
    if isinstance(opts, dict):
        cfg = opts.get('config')
        if not isinstance(cfg, dict):
            cfg = {}
            opts['config'] = cfg
        cfg[key] = value
        return
    cfg = getattr(opts, 'config', None)
    if not isinstance(cfg, dict):
        cfg = {}
        try:
            opts.config = cfg
        except AttributeError:
            return
    cfg[key] = value
def _finite_block(src, dst):
    """The leading block of SRC written into DST, every non-finite entry left at
    zero: a NaN is ABSENT, not a number, and must not enter the probEnv blend."""
    if src is None:
        return
    a = np.atleast_2d(np.asarray(src, dtype=float))
    m = min(dst.shape[0], a.shape[0])
    k = min(dst.shape[1], a.shape[1])
    if m == 0 or k == 0:
        return
    blk = np.array(a[:m, :k], dtype=float)
    blk[~np.isfinite(blk)] = 0.0
    dst[:m, :k] = blk


def _get_tran_data(tran_result):
    """Extract time and metric arrays from transient result.

    Handles both dict objects (with 't' and 'metric' keys) and
    TranResult objects (with .t and .metric attributes).

    Args:
        tran_result: Either a dict or TranResult object

    Returns:
        Tuple of (t_vals, metric_vals) or (None, None) if invalid
    """
    if tran_result is None:
        return None, None

    # Handle dict objects
    if isinstance(tran_result, dict):
        if 't' in tran_result and 'metric' in tran_result:
            t_vals = tran_result['t']
            metric_vals = tran_result['metric']
            if t_vals is not None and len(t_vals) > 0:
                return np.asarray(t_vals), np.asarray(metric_vals)
        return None, None

    # Handle TranResult objects (have .t and .metric attributes)
    if hasattr(tran_result, 't') and hasattr(tran_result, 'metric'):
        t_vals = tran_result.t
        metric_vals = tran_result.metric
        if t_vals is not None and len(t_vals) > 0:
            return np.asarray(t_vals), np.asarray(metric_vals)
        return None, None

    return None, None


class Environment:
    """
    A random environment model where a queueing network operates under
    different environmental conditions (stages).

    The environment switches between stages according to a Markov process,
    and each stage has its own network model with potentially different
    parameters.

    This class mirrors MATLAB's Environment class with full parity.
    """

    __getattr__ = _snake_camel_getattr


    def findSolver(self, metric: str = '', showAll: bool = False):
        """Which solvers and solver methods can analyze THIS model.

            model.findSolver()                # every (solver, method) pair that runs
            model.findSolver('cdf')           # ... that returns a passage-time law
            model.findSolver('getCdfRespT')   # the same question, asked by accessor
            model.findSolver('', True)        # also the pairs that are refused, and why

        The returned DataFrame has one row per pair, with columns Solver,
        Method, Runnable, Class ('exact', 'approx', 'bound' or 'simulation'),
        Metrics and Reason. Method is the method name to pass as a solver method, so
        a row can be acted on directly::

            T = model.findSolver('cdf')
            solver = LINE(model, T.Method[0])

        findMethod and help are aliases of this method.

        Args:
            metric: measure group ('cdf') or accessor ('getCdfRespT') to narrow
                the report to; '' or 'any' keeps every pair.
            showAll: also list the refused pairs, with the reason each was
                refused.

        Returns:
            pandas.DataFrame with the six columns above.
        """
        # The gate lives in SolverAUTO, which is the class that already knows
        # every family, how to build one and what each refuses. Asking it here
        # rather than reimplementing the walk is what keeps the model's answer
        # and AUTO's own dispatch from being two opinions.
        from .solvers.solver_auto.solver_auto import SolverAUTO
        # silenced() wraps the CONSTRUCTION too: it probes every candidate
        # with supports(model), which warns on a model one of them refuses.
        with SolverAUTO.silenced():
            auto = SolverAUTO(self, verbose=False)
        return auto.findSolver(metric, showAll)

    def findMethod(self, metric: str = '', showAll: bool = False):
        """Alias of findSolver: which solvers and solver methods can analyze
        this model.

        The two names exist because the question is asked both ways round --
        "which solver do I use" and "which method do I pass" -- and the answer
        is the same table, whose Method column carries the method name either caller
        needs.
        """
        return self.findSolver(metric, showAll)

    def help(self, metric: str = '', showAll: bool = False):
        """Alias of findSolver: what can this model be solved with?"""
        return self.findSolver(metric, showAll)

    find_solver = findSolver
    find_method = findMethod

    def __init__(self, name: str, num_stages: int = 0):
        """
        Create a random environment model.

        Args:
            name: Name of the environment model
            num_stages: Number of environmental stages (can be 0 if stages added later)
        """
        self.name = name
        self.num_stages = num_stages

        # Stage information
        self._stages: List[Dict[str, Any]] = []
        self._stage_names: List[str] = []
        self._stage_types: List[str] = []
        self._models: List[Any] = []

        # Transition information - env[e][h] is distribution from e to h
        self.env: List[List[Any]] = []
        self._transitions: Dict[tuple, Any] = {}
        self._reset_rules: Dict[tuple, Callable] = {}

        # Node breakdown/repair descriptors recorded by add_node_breakdown /
        # add_node_repair, so that the macros can be serialized declaratively.
        self._node_failures: List[Dict[str, Any]] = []

        # Environment probabilities (computed by init())
        self.probEnv: Optional[np.ndarray] = None  # steady-state probs
        self.probOrig: Optional[np.ndarray] = None  # transition origin probs

        # Hold time distributions
        self.holdTime: List[Any] = []  # hold time distribution for each stage
        self.proc: List[List[Any]] = []  # proc[e][h] = distribution from e to h

        # Reset functions - resetFun[e][h] transforms queue lengths from e to h
        self.resetFun: List[List[Callable]] = []

        # Initialize empty stages if num_stages provided
        for _ in range(num_stages):
            self._stages.append({})
            self._stage_names.append('')
            self._stage_types.append('')
            self._models.append(None)

        # Initialize env matrix
        self._init_env_matrix(num_stages)

    def _init_env_matrix(self, E: int):
        """Initialize the environment transition matrix."""
        self.env = [[None for _ in range(E)] for _ in range(E)]
        self.proc = [[None for _ in range(E)] for _ in range(E)]
        self.resetFun = [[lambda q: q for _ in range(E)] for _ in range(E)]

    def add_stage(self, index, name=None, stage_type: str = None, model: Any = None) -> None:
        """
        Add or update a stage in the environment.

        Two forms, so a MATLAB script transliterates unchanged:

        - ``add_stage(index, name, type, model)`` places the stage at a 0-based
          index, overwriting whatever was there.
        - ``add_stage(name, type, model)`` APPENDS, as MATLAB's
          ``Environment.addStage(name, type, model)`` does. Selected when the
          first argument is a string.

        Args:
            index: Stage index (0-based), or the stage NAME in the append form
            name: Name of the stage, or the type in the append form
            stage_type: Type of stage ('UP', 'DOWN', ...), or the model
            model: Network model for this stage
        """
        if isinstance(index, str):
            index, name, stage_type, model = len(self._stages), index, name, stage_type
        # Expand arrays if needed
        while len(self._stages) <= index:
            self._stages.append({})
            self._stage_names.append('')
            self._stage_types.append('')
            self._models.append(None)

        self._stages[index] = {'name': name, 'type': stage_type, 'model': model}
        self._stage_names[index] = name
        self._stage_types[index] = stage_type
        self._models[index] = model

        self.num_stages = max(self.num_stages, index + 1)

        # Reinitialize env matrix if size changed
        if len(self.env) < self.num_stages:
            self._init_env_matrix(self.num_stages)

    def _stage_index(self, stage) -> int:
        """The 0-based index of a stage given by index or by name."""
        if isinstance(stage, str):
            try:
                return self._stage_names.index(stage)
            except ValueError:
                raise ValueError("[%s] no environment stage named '%s'; the stages are %s"
                                 % (self.name, stage, self._stage_names))
        return int(stage)

    def add_transition(self, from_stage, to_stage, distribution: Any,
                       reset_rule: Optional[Callable[[np.ndarray], np.ndarray]] = None) -> None:
        """
        Add a transition between stages with an optional reset rule.

        Either endpoint may be given as a stage NAME instead of an index, as in
        MATLAB's ``Environment.addTransition(fromName, toName, dist)``; an
        unknown name is an error rather than a silently appended stage.

        Args:
            from_stage: Source stage index, or its name
            to_stage: Destination stage index, or its name
            distribution: Distribution for the transition time (e.g., Exp(rate))
            reset_rule: Optional function that transforms queue lengths when transition occurs.
        """
        from_stage = self._stage_index(from_stage)
        to_stage = self._stage_index(to_stage)
        # Ensure env matrix is large enough
        E = max(from_stage + 1, to_stage + 1, self.num_stages)
        if len(self.env) < E:
            self._init_env_matrix(E)

        self._transitions[(from_stage, to_stage)] = distribution
        self.env[from_stage][to_stage] = distribution
        self.proc[from_stage][to_stage] = distribution

        if reset_rule is not None:
            self._reset_rules[(from_stage, to_stage)] = reset_rule
            self.resetFun[from_stage][to_stage] = reset_rule
        else:
            self._reset_rules[(from_stage, to_stage)] = lambda q: q
            self.resetFun[from_stage][to_stage] = lambda q: q

    def init(self):
        """
        Initialize environment probabilities and hold time distributions.

        This method uses MMAP (Marked MAP) representation with Kronecker products
        to properly handle competing phase-type transitions, matching MATLAB's
        implementation for full parity.

        It computes:
        - probEnv: steady-state probabilities for each stage
        - probOrig: transition origin probabilities
        - holdTime: hold time MMAP representations for each stage
        """
        E = self.num_stages
        if E == 0:
            return

        Pemb = np.zeros((E, E))  # Embedded DTMC transition matrix

        # Build MMAP representations for each transition
        # emmap[e][h] is the MMAP representation for transition from e to h
        emmap = [[None for _ in range(E)] for _ in range(E)]

        for e in range(E):
            for h in range(E):
                if self.env[e][h] is not None:
                    D0, D1 = _get_map_representation(self.env[e][h])
                    # Create MMAP with E+2 elements: [D0, D1, D1_class0, ..., D1_classE-1]
                    mmap_eh = [D0, D1]
                    for j in range(E):
                        if j == h:
                            mmap_eh.append(D1.copy())
                        else:
                            mmap_eh.append(np.zeros_like(D1))
                    emmap[e][h] = mmap_eh
                else:
                    # Disabled transition - use zero matrices
                    emmap[e][h] = [np.zeros((1, 1)), np.zeros((1, 1))] + [np.zeros((1, 1)) for _ in range(E)]

        # Compute hold time MMAP for each stage by combining competing transitions
        self.holdTime = []
        for e in range(E):
            # SEEDED WITH THE SELF TRANSITION, then every other destination folded
            # in, exactly as Environment.init does. Skipping e -> e made the stage
            # sojourn the time to LEAVE rather than the time to SWITCH: on
            # renv_twostages_repairmen, whose Stage2 competes a 0.5 self arc with a
            # 0.5 arc back, the mean read 2 against MATLAB's 1, and the self origin
            # vanished from probOrig (1.0/0.0 against 0.5/0.5). A disabled arc is
            # the 1x1 zero pair and is the neutral element of the superposition, so
            # seeding with it costs nothing when there is no self transition.
            hold_mmap = [m.copy() for m in emmap[e][e]]
            for h in range(E):
                if h != e and self.env[e][h] is not None:
                    # Combine using Kronecker sums (following MATLAB algorithm)
                    n1 = hold_mmap[0].shape[0]
                    n2 = emmap[e][h][0].shape[0]

                    # D0: Kronecker sum
                    D0_new = _krons(hold_mmap[0], emmap[e][h][0])

                    # D1 and class-specific matrices: Kronecker sum, then redirect to first column
                    new_mmap = [D0_new]
                    for j in range(1, E + 2):  # indices 1 to E+1 (D1 and class-specific)
                        Dj_combined = _krons(hold_mmap[j], emmap[e][h][j])
                        completion_rates = Dj_combined @ np.ones(n1 * n2)
                        Dj_new = np.zeros((n1 * n2, n1 * n2))
                        Dj_new[:, 0] = completion_rates
                        new_mmap.append(Dj_new)

                    hold_mmap = _mmap_normalize(new_mmap)

            self.holdTime.append(hold_mmap)

            # Compute embedded transition probabilities from this stage
            count_lambda = _mmap_count_lambda(hold_mmap)
            total_lambda = np.sum(count_lambda)
            if total_lambda > 0:
                Pemb[e, :] = count_lambda / total_lambda
            else:
                Pemb[e, e] = 1.0  # Self-loop if no outgoing transitions

        # Compute holding rates lambda[e] = 1/map_mean(holdTime[e])
        lam = np.zeros(E)
        for e in range(E):
            D0 = self.holdTime[e][0]
            D1 = self.holdTime[e][1]
            mean_hold = _map_mean(D0, D1)
            lam[e] = 1.0 / mean_hold if mean_hold > 0 and np.isfinite(mean_hold) else 0.0

        # Build infinitesimal generator A[e,h] = -lambda[e] * (I[e,h] - Pemb[e,h])
        A = np.zeros((E, E))
        I = np.eye(E)
        for e in range(E):
            for h in range(E):
                A[e, h] = -lam[e] * (I[e, h] - Pemb[e, h])

        # Compute steady-state probabilities
        if np.all(lam > 0):
            self.probEnv = self._solve_ctmc(A)
        else:
            # Fall back to uniform if some stages have no outgoing transitions
            self.probEnv = np.ones(E) / E

        # Compute transition origin probabilities probOrig[h,e]
        self.probOrig = np.zeros((E, E))
        for e in range(E):
            for h in range(E):
                self.probOrig[h, e] = self.probEnv[h] * lam[h] * Pemb[h, e]
            total = np.sum(self.probOrig[:, e])
            if total > 0:
                self.probOrig[:, e] /= total

    def _solve_ctmc(self, Q: np.ndarray) -> np.ndarray:
        """Solve CTMC for steady-state probabilities."""
        E = Q.shape[0]
        if E == 0:
            return np.array([])

        # Solve pi * Q = 0, sum(pi) = 1
        A = Q.T.copy()
        A[-1, :] = 1.0
        b = np.zeros(E)
        b[-1] = 1.0

        try:
            pi = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            pi = np.linalg.lstsq(A, b, rcond=None)[0]

        # Ensure non-negative
        pi = np.maximum(pi, 0)
        pi = pi / np.sum(pi)
        return pi

    def get_stage(self, index: int) -> Dict[str, Any]:
        """Get stage information by index."""
        if 0 <= index < len(self._stages):
            return self._stages[index]
        return {}

    def get_model(self, index: int) -> Any:
        """Get the network model for a stage."""
        if 0 <= index < len(self._models):
            return self._models[index]
        return None

    def get_transition(self, from_stage: int, to_stage: int) -> Any:
        """Get the transition distribution between stages."""
        return self._transitions.get((from_stage, to_stage), None)

    def get_reset_rule(self, from_stage: int, to_stage: int) -> Optional[Callable]:
        """Get the reset rule for a transition between stages."""
        return self._reset_rules.get((from_stage, to_stage), None)

    def get_transition_rate_matrix(self) -> np.ndarray:
        """Build the transition rate matrix for the environment."""
        E = self.num_stages
        Q = np.zeros((E, E))

        for (i, j), dist in self._transitions.items():
            Q[i, j] = _get_rate(dist)

        for i in range(E):
            Q[i, i] = -np.sum(Q[i, :])

        return Q

    def get_steady_state_probs(self) -> np.ndarray:
        """Compute steady-state probabilities for the environment stages."""
        if self.probEnv is None:
            self.init()
        return self.probEnv

    def getEnsemble(self) -> List[Any]:
        """Get list of network models."""
        return self._models

    def stage_table(self) -> pd.DataFrame:
        """Get a table summarizing the environment stages."""
        data = []
        for i in range(len(self._stages)):
            model_name = 'None'
            if i < len(self._models) and self._models[i] is not None:
                if hasattr(self._models[i], 'name'):
                    model_name = self._models[i].name
            row = {
                'Stage': i,
                'Name': self._stage_names[i] if i < len(self._stage_names) else '',
                'Type': self._stage_types[i] if i < len(self._stage_types) else '',
                'Model': model_name
            }
            data.append(row)

        df = pd.DataFrame(data)
        print(df.to_string(index=False))
        return df

    def getStageTable(self):
        """Alias for stage_table (MATLAB compatibility)."""
        return self.stage_table()

    def get_reliability_table(self):
        """Compute system-wide reliability metrics for a breakdown/repair
        environment (one UP stage and one or more DOWN_* stages).

        Returns:
            dict with keys MTTF, MTTR, MTBF and Availability:
              - MTTF: mean time to failure (combined failure rate over UP->DOWN_*)
              - MTTR: mean time to repair (probability-weighted over DOWN_* states)
              - MTBF: MTTF + MTTR
              - Availability: steady-state probability of the UP state

        References:
            MATLAB: matlab/src/lang/Environment.m (getReliabilityTable)
        """
        if self.probEnv is None:
            self.init()
        E = len(self._stage_names)
        if E == 0:
            raise ValueError(
                "Environment has no stages. Add stages before computing reliability metrics.")
        up_idx = None
        down_idx = []
        for i, name in enumerate(self._stage_names):
            if name == 'UP':
                up_idx = i
            elif name.startswith('DOWN_') or name == 'DOWN':
                down_idx.append(i)
        if up_idx is None:
            raise ValueError(
                "No UP stage found. Use add_node_breakdown/add_node_repair to "
                "configure breakdown/repair transitions.")
        if not down_idx:
            raise ValueError(
                "No DOWN stages found. Use add_node_breakdown/add_node_repair to "
                "configure breakdown/repair transitions.")

        def _is_active(dist):
            return dist is not None and not (hasattr(dist, 'isDisabled') and dist.isDisabled())

        # Breakdown rates (UP -> DOWN_*): competing risks
        breakdown_rates = []
        for h in down_idx:
            dist = self.env[up_idx][h]
            if _is_active(dist):
                breakdown_rates.append(1.0 / dist.getMean())
        if not breakdown_rates:
            raise ValueError("No breakdown transitions found (UP -> DOWN_*).")
        lambda_total = sum(breakdown_rates)
        mttf = 1.0 / lambda_total

        # Repair rates (DOWN_* -> UP), weighted by steady-state DOWN probabilities
        repair_rates = []
        down_probs = []
        for e in down_idx:
            dist = self.env[e][up_idx]
            if _is_active(dist):
                repair_rates.append(1.0 / dist.getMean())
                down_probs.append(float(self.probEnv[e]))
        if not repair_rates:
            raise ValueError("No repair transitions found (DOWN_* -> UP).")
        total_down_prob = sum(down_probs)
        if total_down_prob > 0:
            mttr = sum((p / total_down_prob) / mu for p, mu in zip(down_probs, repair_rates))
        else:
            mttr = float(np.mean([1.0 / mu for mu in repair_rates]))

        mtbf = mttf + mttr
        avail_up = float(self.probEnv[up_idx])
        avail_down = float(sum(self.probEnv[e] for e in down_idx))
        availability = avail_up / (avail_up + avail_down) if (avail_up + avail_down) > 0 else 0.0

        return {'MTTF': mttf, 'MTTR': mttr, 'MTBF': mtbf, 'Availability': availability}

    getReliabilityTable = get_reliability_table

    def rel_t(self):
        """Short alias for get_reliability_table."""
        return self.get_reliability_table()

    relT = rel_t
    getRelT = get_reliability_table
    relTable = get_reliability_table
    getRelTable = get_reliability_table

    # CamelCase aliases for MATLAB API compatibility
    def addStage(self, index, name=None, stage_type: str = None, model: Any = None) -> None:
        """Alias for add_stage (MATLAB compatibility), including its append form
        ``addStage(name, type, model)``."""
        return self.add_stage(index, name, stage_type, model)

    def addTransition(self, from_stage, to_stage, distribution: Any,
                      reset_rule: Optional[Callable[[np.ndarray], np.ndarray]] = None) -> None:
        """Alias for add_transition (MATLAB compatibility); either endpoint may
        be a stage name."""
        return self.add_transition(from_stage, to_stage, distribution, reset_rule)

    def print_stage_table(self) -> None:
        """Print detailed stage table with transition information."""
        print("Stage Table:")
        print("============")
        for i in range(len(self._stages)):
            name = self._stage_names[i] if i < len(self._stage_names) else ''
            stage_type = self._stage_types[i] if i < len(self._stage_types) else ''
            model = self._models[i] if i < len(self._models) else None
            model_name = model.name if model is not None and hasattr(model, 'name') else 'None'

            n_nodes = 0
            n_classes = 0
            if model is not None:
                if hasattr(model, 'get_nodes'):
                    n_nodes = len(model.get_nodes())
                elif hasattr(model, 'nodes'):
                    n_nodes = len(model.nodes)
                if hasattr(model, 'get_classes'):
                    n_classes = len(model.get_classes())
                elif hasattr(model, 'classes'):
                    n_classes = len(model.classes)

            print(f"Stage {i + 1}: {name} (Type: {stage_type})")
            print(f"  - Network: {model_name}")
            print(f"  - Nodes: {n_nodes}")
            print(f"  - Classes: {n_classes}")

        if self._transitions:
            print("Transitions:")
            for (from_idx, to_idx), dist in sorted(self._transitions.items()):
                from_name = self._stage_names[from_idx] if from_idx < len(self._stage_names) else f'Stage{from_idx}'
                to_name = self._stage_names[to_idx] if to_idx < len(self._stage_names) else f'Stage{to_idx}'
                rate = _get_rate(dist)
                print(f"  {from_name} -> {to_name}: rate = {rate:.4f}")

    @property
    def ensemble(self) -> List[Any]:
        """Get list of network models (MATLAB compatibility)."""
        return self._models

    def _find_stage_by_name(self, name: str) -> int:
        """Find stage index by name, returns -1 if not found."""
        for i, stage_name in enumerate(self._stage_names):
            if stage_name == name:
                return i
        return -1

    def _get_node_name(self, node_or_name) -> str:
        """Extract node name from a Node object or string."""
        if hasattr(node_or_name, 'name'):
            return node_or_name.name
        return str(node_or_name)

    @staticmethod
    def resolve_reset_policy(spec):
        """Resolve a queue-length reset policy given a named policy or a callable.

        Named policies (the only serializable ones):
            'keep'  - carry the queue lengths across the transition, q -> q
            'clear' - empty the queues on the transition, q -> 0

        A callable is returned unchanged and reported as 'custom': an arbitrary
        reset function cannot be reproduced from JSON.

        Returns:
            tuple: (reset_fun, reset_name)
        """
        if spec is None:
            return (lambda q: q), 'keep'
        if isinstance(spec, str):
            name = spec.lower()
            if name == 'keep':
                return (lambda q: q), 'keep'
            if name == 'clear':
                return (lambda q: np.zeros_like(q)), 'clear'
            raise ValueError('Unknown reset policy "%s". Use "keep", "clear", or a callable q -> q.' % spec)
        if callable(spec):
            return spec, 'custom'
        raise ValueError('Reset policy must be a callable q -> q, or one of the named policies "keep" and "clear".')

    def find_node_failure(self, node_name: str) -> int:
        """Index of the node-failure descriptor for node_name, or -1 if absent."""
        for i, nf in enumerate(self._node_failures):
            if nf['node'] == node_name:
                return i
        return -1

    def register_node_failure(self, node_name, breakdown_dist, repair_dist, down_service_dist,
                              breakdown_policy='keep', repair_policy='keep'):
        """Attach a node breakdown/repair descriptor to stages that already exist.

        This is the counterpart of add_node_breakdown/add_node_repair for the case
        where the UP and DOWN_<node> stages and their transitions have already been
        built (for instance by load_model reading the expanded stages/transitions
        form). It records the descriptor and applies the queue-length reset
        policies, which the expanded form cannot carry.
        """
        down_stage_name = 'DOWN_%s' % node_name
        up_idx = self._find_stage_by_name('UP')
        down_idx = self._find_stage_by_name(down_stage_name)
        if up_idx < 0:
            raise ValueError('Cannot register a node failure on "%s": no UP stage is defined in this '
                             'environment.' % node_name)
        if down_idx < 0:
            raise ValueError('Cannot register a node failure on "%s": no "%s" stage is defined in this '
                             'environment.' % (node_name, down_stage_name))

        breakdown_fun, breakdown_name = Environment.resolve_reset_policy(breakdown_policy)
        self._reset_rules[(up_idx, down_idx)] = breakdown_fun
        self.resetFun[up_idx][down_idx] = breakdown_fun
        repair_name = ''
        if repair_dist is not None:
            repair_fun, repair_name = Environment.resolve_reset_policy(repair_policy)
            self._reset_rules[(down_idx, up_idx)] = repair_fun
            self.resetFun[down_idx][up_idx] = repair_fun

        nf = {'node': node_name, 'breakdown': breakdown_dist, 'downService': down_service_dist,
              'repair': repair_dist, 'breakdownResetPolicy': breakdown_name,
              'repairResetPolicy': repair_name}
        idx = self.find_node_failure(node_name)
        if idx >= 0:
            self._node_failures[idx] = nf
        else:
            self._node_failures.append(nf)

    def add_node_breakdown(self, base_model, node_or_name, breakdown_dist, down_service_dist,
                           reset_fun=None):
        """
        Add UP and DOWN stages for a node that can break down.

        Args:
            base_model: The base network model with normal (UP) service rates
            node_or_name: Node object or name of the node that can break down
            breakdown_dist: Distribution for time until breakdown (UP->DOWN transition)
            down_service_dist: Service distribution when the node is down
            reset_fun: Optional reset policy for queue lengths on breakdown. Either a
                callable q -> q, or one of the named policies 'keep' (identity, the
                default) and 'clear' (empty the queues). Only named policies are
                serializable.
        """
        node_name = self._get_node_name(node_or_name)

        reset_fun, reset_name = Environment.resolve_reset_policy(reset_fun)

        # Create UP stage if this is the first call
        up_idx = self._find_stage_by_name('UP')
        if up_idx < 0:
            up_model = base_model.copy() if hasattr(base_model, 'copy') else base_model
            # Use first available slot (0) if pre-allocated stages exist
            up_idx = 0
            self.add_stage(up_idx, 'UP', 'operational', up_model)

        # Create DOWN stage with modified service rate for the specified node
        down_model = base_model.copy() if hasattr(base_model, 'copy') else base_model
        nodes = down_model.get_nodes() if hasattr(down_model, 'get_nodes') else []

        node_idx = -1
        for i, n in enumerate(nodes):
            n_name = n.name if hasattr(n, 'name') else str(n)
            if n_name == node_name:
                node_idx = i
                break

        if node_idx < 0:
            raise ValueError(f'Node "{node_name}" not found in the base model.')

        # Update service distribution for the down node
        classes = down_model.get_classes() if hasattr(down_model, 'get_classes') else []
        for cls in classes:
            if hasattr(nodes[node_idx], 'set_service'):
                nodes[node_idx].set_service(cls, down_service_dist)
            elif hasattr(nodes[node_idx], 'setService'):
                nodes[node_idx].setService(cls, down_service_dist)

        # Add DOWN stage
        down_stage_name = f'DOWN_{node_name}'
        # Find next available index - either UP+1 or first empty slot after UP
        down_idx = up_idx + 1
        self.add_stage(down_idx, down_stage_name, 'failed', down_model)

        # Reinitialize env matrix if needed
        E = max(len(self._stages), self.num_stages)
        if len(self.env) < E:
            self._init_env_matrix(E)

        # Add breakdown transition (UP -> DOWN)
        self.add_transition(up_idx, down_idx, breakdown_dist, reset_fun)

        # see _kb/06-solver-catalog.md (ENV: Python environment.py additional notes) for rationale
        self._node_failures.append({
            'node': node_name,
            'breakdown': breakdown_dist,
            'downService': down_service_dist,
            'repair': None,
            'breakdownResetPolicy': reset_name,
            'repairResetPolicy': '',
        })

    def add_node_repair(self, node_or_name, repair_dist, reset_fun=None):
        """
        Add repair transition from DOWN to UP stage for a previously added breakdown.

        Args:
            node_or_name: Node object or name of the node that can be repaired
            repair_dist: Distribution for repair time (DOWN->UP transition)
            reset_fun: Optional reset policy for queue lengths on repair. Either a
                callable q -> q, or one of the named policies 'keep' (identity, the
                default) and 'clear' (empty the queues). Only named policies are
                serializable.
        """
        node_name = self._get_node_name(node_or_name)

        reset_fun, reset_name = Environment.resolve_reset_policy(reset_fun)

        down_stage_name = f'DOWN_{node_name}'
        down_idx = self._find_stage_by_name(down_stage_name)
        up_idx = self._find_stage_by_name('UP')

        if down_idx < 0:
            raise ValueError(f'DOWN stage for node "{node_name}" not found. Call add_node_breakdown first.')
        if up_idx < 0:
            raise ValueError('UP stage not found. Call add_node_breakdown first.')

        # Add repair transition (DOWN -> UP)
        self.add_transition(down_idx, up_idx, repair_dist, reset_fun)

        # Complete the descriptor recorded by add_node_breakdown.
        idx = self.find_node_failure(node_name)
        if idx >= 0:
            self._node_failures[idx]['repair'] = repair_dist
            self._node_failures[idx]['repairResetPolicy'] = reset_name

    def add_node_failure_repair(self, base_model, node_or_name, breakdown_dist, repair_dist,
                                down_service_dist, reset_breakdown=None,
                                reset_repair=None):
        """
        Convenience method to add both breakdown and repair for a node.

        Args:
            base_model: The base network model with normal (UP) service rates
            node_or_name: Node object or name of the node that can break down and repair
            breakdown_dist: Distribution for time until breakdown
            repair_dist: Distribution for repair time
            down_service_dist: Service distribution when the node is down
            reset_breakdown: Optional reset function for breakdown transition
            reset_repair: Optional reset function for repair transition
        """
        node_name = self._get_node_name(node_or_name)

        self.add_node_breakdown(base_model, node_name, breakdown_dist, down_service_dist, reset_breakdown)
        self.add_node_repair(node_name, repair_dist, reset_repair)

    def set_breakdown_reset_policy(self, node_or_name, reset_fun):
        """
        Update the reset policy for breakdown transitions (UP -> DOWN) of a node.

        Args:
            node_or_name: Node object or name of the node
            reset_fun: Reset policy for queue lengths on breakdown. Either a callable
                q -> q, or one of the named policies 'keep' and 'clear'. Only named
                policies are serializable.
        """
        node_name = self._get_node_name(node_or_name)
        reset_fun, reset_name = Environment.resolve_reset_policy(reset_fun)
        down_stage_name = f'DOWN_{node_name}'

        up_idx = self._find_stage_by_name('UP')
        down_idx = self._find_stage_by_name(down_stage_name)

        if up_idx < 0:
            raise ValueError('UP stage not found. Call add_node_breakdown first.')
        if down_idx < 0:
            raise ValueError(f'DOWN stage for node "{node_name}" not found. Call add_node_breakdown first.')

        # Update the reset function for the breakdown transition (UP -> DOWN)
        self._reset_rules[(up_idx, down_idx)] = reset_fun
        self.resetFun[up_idx][down_idx] = reset_fun
        idx = self.find_node_failure(node_name)
        if idx >= 0:
            self._node_failures[idx]['breakdownResetPolicy'] = reset_name

    def set_repair_reset_policy(self, node_or_name, reset_fun):
        """
        Update the reset policy for repair transitions (DOWN -> UP) of a node.

        Args:
            node_or_name: Node object or name of the node
            reset_fun: Reset policy for queue lengths on repair. Either a callable
                q -> q, or one of the named policies 'keep' and 'clear'. Only named
                policies are serializable.
        """
        node_name = self._get_node_name(node_or_name)
        reset_fun, reset_name = Environment.resolve_reset_policy(reset_fun)
        down_stage_name = f'DOWN_{node_name}'

        up_idx = self._find_stage_by_name('UP')
        down_idx = self._find_stage_by_name(down_stage_name)

        if up_idx < 0:
            raise ValueError('UP stage not found. Call add_node_breakdown first.')
        if down_idx < 0:
            raise ValueError(f'DOWN stage for node "{node_name}" not found. Call add_node_breakdown first.')

        # Update the reset function for the repair transition (DOWN -> UP)
        self._reset_rules[(down_idx, up_idx)] = reset_fun
        self.resetFun[down_idx][up_idx] = reset_fun
        idx = self.find_node_failure(node_name)
        if idx >= 0:
            self._node_failures[idx]['repairResetPolicy'] = reset_name

    def get_stage_table(self):
        """Get stage table (MATLAB compatibility alias)."""
        return self.stage_table()

    def get_ensemble(self) -> List[Any]:
        """Get list of network models (snake_case version)."""
        return self._models

    # CamelCase aliases (portability with the wrapper / MATLAB API)
    addNodeBreakdown = add_node_breakdown
    addNodeRepair = add_node_repair
    addNodeFailureRepair = add_node_failure_repair
    setBreakdownResetPolicy = set_breakdown_reset_policy
    setRepairResetPolicy = set_repair_reset_policy
    getSteadyStateProbs = get_steady_state_probs
    getEnsemble = get_ensemble


class _ExponentialDist:
    """Simple exponential distribution for hold times."""

    def __init__(self, rate: float):
        self.rate = rate

    def getRate(self) -> float:
        return self.rate

    def evalCDF(self, t: float) -> float:
        if self.rate <= 0:
            return 0.0
        return 1.0 - np.exp(-self.rate * t)


def _env_verbose(value) -> bool:
    """Whether this `verbose` asks for output.

    A PLAIN TRUTHINESS TEST IS WRONG HERE, because the slot holds a
    `VerboseLevel` and every Enum member is truthy -- `VerboseLevel.SILENT`
    included, value 0 notwithstanding. Asking ENV for SILENT therefore got STD:
    `renv_basic` printed its `ENV analysis [...]` banner, and with it the
    environment-averaged table stopped reading as the stage solver's (see
    `Solver._process_verbose_option`, which makes the same distinction for
    every network solver).
    """
    if isinstance(value, VerboseLevel):
        return value != VerboseLevel.SILENT
    # `verbose = 0` is the MATLAB spelling and reaches here as a plain int, so
    # anything outside the enum is read by truthiness, where it is meaningful.
    return bool(value)


def _mf_hit_by_list(xbar, lam, k):
    """Hit probability of class K resolved by cache list, from the mean-field
    occupancy XBAR (flat DDPP state, item-major over lists 0..h) and the
    isolated per-item request rates LAM. List 0 is outside the cache, so the
    hit rows are lists 1..h and sum to the total hit probability of the class."""
    if xbar is None or lam is None:
        return None
    xbar = np.asarray(xbar, dtype=float).ravel()
    lam = np.asarray(lam, dtype=float)
    if xbar.size == 0 or lam.ndim < 3 or k >= lam.shape[0]:
        return None
    nitems = lam.shape[1]
    if nitems <= 0 or xbar.size % nitems != 0:
        return None
    h = xbar.size // nitems - 1
    if h < 1:
        return None
    wpop = np.asarray(lam[k, :, 0], dtype=float).ravel()
    wpop = np.where(np.isfinite(wpop), wpop, 0.0)
    tot = float(np.sum(wpop))
    if tot <= 0:
        return None
    wpop = wpop / tot
    hbl = np.zeros(h)
    for l in range(1, h + 1):
        hbl[l - 1] = min(1.0, max(0.0, float(wpop @ xbar[l * nitems:(l + 1) * nitems])))
    return hbl


class SolverENV:
    """
    Environment solver for random environment models.

    This solver analyzes queueing networks operating in random environments
    by running transient analysis with iteration until convergence.

    Implements full parity with MATLAB SolverENV.
    """

    __getattr__ = _snake_camel_getattr

    def __init__(self, env_model: Environment, solvers: Union[List, Callable], options=None):
        """
        Create an environment solver.

        Args:
            env_model: Environment model to analyze
            solvers: Either a list of solvers (one per stage) or a factory function
            options: Solver options
        """
        self.env_model = env_model
        self.options = self._normalize_options(options)

        # Handle solvers as list or factory function
        if callable(solvers):
            self._solvers = []
            for model in env_model.ensemble:
                if model is not None:
                    self._solvers.append(solvers(model))
                else:
                    self._solvers.append(None)
        else:
            self._solvers = list(solvers)

        # Ensemble state (mirrors MATLAB)
        self.ensemble = env_model.ensemble
        self.results = {}  # results[it, e] = result for iteration it, stage e
        # stages whose unspecified horizon has been announced (see _stage_horizon)
        self._horizon_warned = set()

        # Final result
        self.result = None
        self._result = None

        # State-dependent / semi-Markov / decomposition state (MATLAB @SolverENV parity)
        self.stateDepMethod = ''
        self.SMPMethod = False
        self.newMethod = False
        self.compression = False
        self.compressionResult = None
        self.MS = None
        self.Ecompress = None
        self.E0 = None
        self.Eutil = None
        # see _kb/06-solver-catalog.md (ENV: Python environment.py additional notes) for rationale
        if str(self.options.get('method', 'default')).lower() == 'smp':
            self.SMPMethod = True

    def _normalize_options(self, options) -> Dict:
        """Normalize options to a dictionary."""
        # `lang` and `arith` are carried through because a normalization that
        # drops them would let a caller ask for lang='cpp' and be answered
        # natively, which is the one substitution the delegation forbids.
        if options is None:
            return {
                'method': 'default',
                'iter_max': 100,
                'iter_tol': 1e-4,
                # STD like every other solver; see default_options
                'verbose': VerboseLevel.STD,
                'config': None,
                'lang': _os.environ.get('LINE_SOLVER_LANG', 'python'),
                'arith': None
            }

        if isinstance(options, dict):
            opts = {
                'method': options.get('method', 'default'),
                'iter_max': options.get('iter_max', 100),
                'iter_tol': options.get('iter_tol', 1e-4),
                'verbose': options.get('verbose', VerboseLevel.STD),
                'sojourn': options.get('sojourn', None),
                'config': options.get('config', None),
                'lang': options.get('lang', _os.environ.get('LINE_SOLVER_LANG', 'python')),
                'arith': options.get('arith', None)
            }
            return opts

        # Options object
        return {
            'method': getattr(options, 'method', 'default'),
            'iter_max': getattr(options, 'iter_max', 100),
            'iter_tol': getattr(options, 'iter_tol', 1e-4),
            'verbose': getattr(options, 'verbose', VerboseLevel.STD),
            'sojourn': getattr(options, 'sojourn', None),
            'config': getattr(options, 'config', None),
            'lang': getattr(options, 'lang', _os.environ.get('LINE_SOLVER_LANG', 'python')),
            'arith': getattr(options, 'arith', None)
        }

    def getNumberOfModels(self) -> int:
        """Get number of environment stages."""
        return self.env_model.num_stages

    def getSolver(self, e: int):
        """Get solver for stage e."""
        return self._solvers[e] if e < len(self._solvers) else None

    @property
    def solvers(self):
        """Per-stage solver list (mirrors MATLAB/JAR solvers{e})."""
        return self._solvers

    def getStageResult(self, e: int):
        """Per-stage solver result for stage e from the final iteration
        (mirrors MATLAB results{e} / JAR getStageResult(e)). Returns None if
        the solver has not been run."""
        if not self.results:
            return None
        it = max(k[0] for k in self.results.keys())
        return self.results.get((it, e))

    def get_stage_result(self, e: int):
        """snake_case alias for getStageResult."""
        return self.getStageResult(e)

    # ------------------------------------------------------------------
    # State-dependent / semi-Markov / NCD-decomposition configuration
    # (mirrors MATLAB @SolverENV: setStateDepMethod, setSMPMethod,
    # setNewMethod, setCompression, ctmc_decompose, findBestPartition,
    # applyCompression). __getattr__ resolves only instance attributes, so
    # both snake_case and camelCase method spellings are provided explicitly.
    # ------------------------------------------------------------------

    def setStateDepMethod(self, method):
        """Set the state-dependent decomposition method (non-empty string)."""
        if method is None or method == '':
            raise ValueError('State-dependent method cannot be null or empty.')
        self.stateDepMethod = method
        return self

    def set_state_dep_method(self, method):
        """snake_case alias for setStateDepMethod."""
        return self.setStateDepMethod(method)

    def setSMPMethod(self, flag):
        """Enable/disable DTMC-based computation for semi-Markov environments."""
        self.SMPMethod = bool(flag)
        return self

    def set_smp_method(self, flag):
        """snake_case alias for setSMPMethod."""
        return self.setSMPMethod(flag)

    def setNewMethod(self, flag):
        """Enable/disable the DTMC-based embedded-chain environment computation."""
        self.newMethod = bool(flag)
        return self

    def set_new_method(self, flag):
        """snake_case alias for setNewMethod."""
        return self.setNewMethod(flag)

    def setCompression(self, flag):
        """Enable/disable Courtois (NCD) compression of the environment process."""
        self.compression = bool(flag)
        return self

    def set_compression(self, flag):
        """snake_case alias for setCompression."""
        return self.setCompression(flag)

    def _env_generator(self):
        """Build the environment CTMC generator Eutil from stage-transition
        rates (E0[e,h] = rate of env[e][h]); mirrors MATLAB init()."""
        from .api.mc.ctmc import ctmc_makeinfgen
        E = self.env_model.num_stages
        E0 = np.zeros((E, E))
        for e in range(E):
            row = self.env_model.env[e] if e < len(self.env_model.env) else []
            for h in range(E):
                dist = row[h] if h < len(row) else None
                if dist is not None:
                    try:
                        E0[e, h] = float(dist.getRate())
                    except Exception:
                        m = float(dist.getMean())
                        E0[e, h] = 1.0 / m if m > 0 else 0.0
        self.E0 = E0
        self.Eutil = np.asarray(ctmc_makeinfgen(E0), dtype=float)
        return self.Eutil

    def ctmc_decompose(self, Q, MS):
        """CTMC nearly-complete-decomposability (NCD) aggregation of generator
        Q under the macro-state partition MS (list of 0-based index lists).

        The algorithm is selected by options.config.da: 'courtois' (default),
        'kms', 'takahashi', 'multi'; options.config.da_iter sets the iteration
        count for kms/takahashi. Mirrors MATLAB SolverENV.ctmc_decompose.

        Returns:
            (p, eps, epsMax, q)
        """
        from .api.mc.aggregation import (ctmc_courtois, ctmc_kms,
                                         ctmc_takahashi, ctmc_multi)
        Q = np.asarray(Q, dtype=float)
        config = self.options.get('config') if isinstance(self.options, dict) else None

        def _cfg(key, default):
            if config is None:
                return default
            v = config.get(key, default) if isinstance(config, dict) else getattr(config, key, default)
            return default if v is None else v

        method = str(_cfg('da', 'courtois')).lower()
        numsteps = int(_cfg('da_iter', 10))
        if method == 'courtois':
            res = ctmc_courtois(Q, MS)
            return res.p, res.eps, res.epsMAX, res.q
        elif method == 'kms':
            res = ctmc_kms(Q, MS, numsteps)
            return res.p, res.eps, res.epsMAX, 1.05 * float(np.max(np.abs(Q)))
        elif method == 'takahashi':
            res = ctmc_takahashi(Q, MS, numsteps)
            return res.p, res.eps, res.epsMAX, 1.05 * float(np.max(np.abs(Q)))
        elif method == 'multi':
            MSS = [[i] for i in range(len(MS))]
            res = ctmc_multi(Q, MS, MSS)
            return res.p, res.eps, res.epsMAX, 1.05 * float(np.max(np.abs(Q)))
        raise ValueError('Unknown decomposition method: %s' % method)

    def findBestPartition(self, E=None):
        """Find a good NCD macro-state partition by exhaustive pairwise merges
        (for small environments). Mirrors MATLAB SolverENV.findBestPartition.

        Returns a list of 0-based index lists (the partition), and sets
        self.MS / self.Ecompress.
        """
        if self.Eutil is None:
            self._env_generator()
        if E is None:
            E = self.env_model.num_stages
        bestMS = [[i] for i in range(E)]
        _, bestEps, _, _ = self.ctmc_decompose(self.Eutil, bestMS)
        if bestEps is None or (isinstance(bestEps, float) and np.isnan(bestEps)):
            self.MS = bestMS
            self.Ecompress = len(bestMS)
            return bestMS
        for i in range(E):
            for j in range(i + 1, E):
                testMS = []
                for k in range(E):
                    if k == i:
                        testMS.append([i, j])
                    elif k != j:
                        testMS.append([k])
                _, testEps, _, _ = self.ctmc_decompose(self.Eutil, testMS)
                if testEps is not None and not np.isnan(testEps) and testEps < bestEps:
                    bestEps = testEps
                    bestMS = testMS
        self.MS = bestMS
        self.Ecompress = len(bestMS)
        return bestMS

    def applyCompression(self):
        """Compute the Courtois/NCD compression of the environment process,
        storing macro/micro stationary probabilities and the partition in
        self.compressionResult. Mirrors the environment-aggregation portion of
        MATLAB SolverENV.applyCompression.

        The ensemble is then REPLACED by the compressed one: each macro-state
        carries a copy of its first micro-state network whose service rates are
        the pmicro-weighted averages over the block, and probEnv/probOrig are
        rebuilt on the macro chain. Computing the partition and leaving the
        ensemble alone reports a compression that never happened -- every
        subsequent solve still ran the full environment.

        Returns:
            self.compressionResult (dict with p, eps, epsMax, q, MS, pMacro, pmicro).
        """
        E = self.env_model.num_stages
        self._env_generator()
        MS = self.findBestPartition(E)
        p, eps, epsMax, q = self.ctmc_decompose(self.Eutil, MS)
        p = np.asarray(p, dtype=float).ravel()
        Ecomp = len(MS)
        pMacro = np.array([float(np.sum(p[MS[i]])) for i in range(Ecomp)])
        pmicro = np.zeros(E)
        for i in range(Ecomp):
            block = p[MS[i]]
            s = float(np.sum(block))
            if s > 0:
                pmicro[np.asarray(MS[i])] = block / s
        self.Ecompress = Ecomp
        self.compressionResult = {
            'p': p, 'eps': eps, 'epsMax': epsMax, 'q': q,
            'MS': MS, 'pMacro': pMacro, 'pmicro': pmicro,
        }
        if eps is not None and epsMax is not None and eps > epsMax:
            warnings.warn('SolverENV: environment cannot be effectively '
                          'compressed (eps > epsMax).')
        if Ecomp < E:
            self._rebuild_compressed_ensemble(MS, pMacro, pmicro)
        return self.compressionResult

    def _computeMacroRate(self, MS, pmicro, from_macro, to_macro):
        """Aggregate transition rate between two macro-states, weighted by the
        micro-state conditional probabilities. Mirrors MATLAB computeMacroRate."""
        rate = 0.0
        for mi in np.asarray(MS[from_macro]).ravel():
            for mj in np.asarray(MS[to_macro]).ravel():
                rate += float(pmicro[int(mi)]) * float(self.E0[int(mi), int(mj)])
        return rate

    def _rebuild_compressed_ensemble(self, MS, pMacro, pmicro):
        """Replace the ensemble, solvers and structs by their macro-state
        aggregates. Mirrors the second half of MATLAB SolverENV.applyCompression.
        """
        from .distributions import Exp

        Ecomp = len(MS)
        sn0 = self.ensemble[0].get_struct()
        M, K = int(sn0.nstations), int(sn0.nclasses)

        # Embedding weights on the macro chain
        new_embweight = np.zeros((Ecomp, Ecomp))
        for e in range(Ecomp):
            denom = sum(pMacro[h] * self._computeMacroRate(MS, pmicro, h, e)
                        for h in range(Ecomp) if h != e)
            for k in range(Ecomp):
                if k != e and denom > 0:
                    new_embweight[k, e] = (pMacro[k]
                                           * self._computeMacroRate(MS, pmicro, k, e)
                                           / denom)
        self.env_model.probEnv = pMacro
        self.env_model.probOrig = new_embweight

        macro_ensemble, macro_solvers, macro_sn = [], [], []
        for i in range(Ecomp):
            block = np.asarray(MS[i]).ravel().astype(int)
            first = int(block[0])
            model_i = self.ensemble[first].copy()
            for m in range(M):
                for k in range(K):
                    rate_sum = 0.0
                    for micro in block:
                        rate_sum += float(pmicro[micro]) * float(
                            np.asarray(self.ensemble[micro].get_struct().rates)[m, k])
                    if not (rate_sum > 0):
                        continue
                    station = model_i.get_stations()[m]
                    jobclass = model_i.get_classes()[k]
                    cls_names = {c.__name__ for c in type(station).__mro__}
                    if cls_names & {'Queue', 'Delay'}:
                        station.setService(jobclass, Exp(rate_sum))
            model_i.refresh_struct()
            macro_ensemble.append(model_i)
            macro_sn.append(model_i.get_struct())
            macro_solvers.append(self._rebuild_stage_solver(first, model_i))

        self.ensemble = macro_ensemble
        self._solvers = macro_solvers
        self.sn = macro_sn

    def _rebuild_stage_solver(self, micro_idx, model):
        """A stage solver for a macro-state network, of the same class and with
        the same options as the micro-state solver it replaces."""
        proto = self._solvers[micro_idx] if micro_idx < len(self._solvers) else None
        if proto is None:
            from .solvers.solver_fld import SolverFLD
            return SolverFLD(model)
        try:
            return type(proto)(model, options=getattr(proto, 'options', None))
        except Exception:
            return type(proto)(model)

    def getSamplePathTable(self, sample_path):
        """Compute transient performance metrics along a sample path through
        environment stages.

        Args:
            sample_path: list of (stage, duration) pairs where ``stage`` is a
                stage name (str) or 0-based index (int) and ``duration`` is the
                positive time spent in that stage.

        Returns:
            pandas.DataFrame with columns Segment, Stage, Duration, Station,
            JobClass, InitQLen, InitUtil, InitTput, FinalQLen, FinalUtil,
            FinalTput.

        References:
            MATLAB: matlab/src/solvers/ENV/@SolverENV/SolverENV.m (getSamplePathTable)
        """
        import pandas as pd
        if not sample_path:
            raise ValueError("Sample path cannot be empty.")
        if self.env_model.probEnv is None:
            self.init()

        E = self.env_model.num_stages
        sn0 = self.ensemble[0].get_struct()
        M = sn0.nstations
        K = sn0.nclasses
        njobs = np.asarray(sn0.njobs).ravel()

        # Initial queue lengths: spread closed-class jobs uniformly over stations
        Q_current = np.zeros((M, K))
        for k in range(K):
            if k < len(njobs) and np.isfinite(njobs[k]) and njobs[k] > 0:
                Q_current[:, k] = njobs[k] / M

        rows = []
        for seg, (stage_spec, duration) in enumerate(sample_path):
            if isinstance(stage_spec, str):
                e = self.env_model._find_stage_by_name(stage_spec)
                if e is None or e < 0:
                    raise ValueError('Stage "%s" not found.' % stage_spec)
                stage_name = stage_spec
            else:
                e = int(stage_spec)
                if e < 0 or e >= E:
                    raise ValueError("Stage index %d out of range [0, %d]." % (e, E - 1))
                stage_name = self.env_model._stage_names[e]
            if duration <= 0:
                raise ValueError("Duration must be positive.")

            model_e = self.ensemble[e]
            model_e.init_from_marginal(Q_current)
            solver_e = self._solvers[e]
            solver_e.options.timespan = [0, duration]
            if hasattr(solver_e, 'reset'):
                solver_e.reset()

            Qt, Ut, Tt = model_e.getTranHandles()
            QNt, UNt, TNt = solver_e.getTranAvg(Qt, Ut, Tt)

            final_Q = np.zeros((M, K))
            for i in range(M):
                for r in range(K):
                    init_q = init_u = init_t = 0.0
                    fin_q = fin_u = fin_t = 0.0
                    qir = QNt[i][r] if QNt[i][r] is not None else None
                    uir = UNt[i][r] if UNt[i][r] is not None else None
                    tir = TNt[i][r] if TNt[i][r] is not None else None
                    if qir is not None and getattr(qir, 'metric', None) is not None and len(qir.metric) > 0:
                        init_q = float(qir.metric[0])
                        fin_q = float(qir.metric[-1])
                    if uir is not None and getattr(uir, 'metric', None) is not None and len(uir.metric) > 0:
                        init_u = float(uir.metric[0])
                        fin_u = float(uir.metric[-1])
                    if tir is not None and getattr(tir, 'metric', None) is not None and len(tir.metric) > 0:
                        init_t = float(tir.metric[0])
                        fin_t = float(tir.metric[-1])
                    final_Q[i, r] = fin_q
                    station_name = sn0.nodenames[int(sn0.stationToNode[i])] if sn0.nodenames else str(i)
                    class_name = sn0.classnames[r] if sn0.classnames else str(r)
                    rows.append({
                        'Segment': seg, 'Stage': stage_name, 'Duration': duration,
                        'Station': station_name, 'JobClass': class_name,
                        'InitQLen': init_q, 'InitUtil': init_u, 'InitTput': init_t,
                        'FinalQLen': fin_q, 'FinalUtil': fin_u, 'FinalTput': fin_t,
                    })
            # Carry final queue lengths into the next segment
            Q_current = final_Q

        return pd.DataFrame(rows)

    get_sample_path_table = getSamplePathTable

    def init(self):
        """Initialize the environment solver.

        Mirrors JAR SolverENV.init():
        - Calls envObj.init() to compute probEnv, probOrig, holdTime
        - Sets ODE max step for sub-solvers to ensure accurate integration
        - Builds E0 (infgen from rates) for embweight computation
        """
        self.env_model.init()
        self.results = {}

        E = self.getNumberOfModels()

        # see _kb/06-solver-catalog.md (ENV: Python environment.py additional notes) for rationale
        for e in range(E):
            solver = self.getSolver(e)
            if solver is None:
                continue
            opts = None
            if hasattr(solver, 'options'):
                opts = solver.options
            elif hasattr(solver, '_native_solver') and hasattr(solver._native_solver, 'options'):
                opts = solver._native_solver.options
            if opts is not None and hasattr(opts, 'timespan') and opts.timespan is not None:
                ts = opts.timespan
                if len(ts) > 1 and np.isfinite(ts[1]) and ts[1] > 0:
                    # Only set if not already set by user. Matches JAR
                    # SolverENV.init(): odemaxstep = tEnd/100 (tighter than the
                    # generic tEnd/10 default, needed for ENV convergence).
                    if hasattr(opts, 'odemaxstep') and getattr(opts, 'odemaxstep', None) is None:
                        opts.odemaxstep = ts[1] / 100.0
                        opts._ode_maxstep_set = True

        if self.SMPMethod or self.newMethod:
            self._smp_stage_probabilities()

    @staticmethod
    def _arc_cdf_grid(dist, t0, step, n):
        """F(t0 + i*step), i = 0..n-1, of one environment arc; zeros for a disabled arc.

        Read off the (D0, D1) form Environment.init uses, F(t) = 1 - pie expm(D0 t) e, which is
        evalCDF for every Markovian law, stepping the row vector by expm(D0*step)."""
        from scipy.linalg import expm
        if dist is None:
            return np.zeros(n)
        D0, D1 = _get_map_representation(dist)
        D0 = np.asarray(D0, dtype=float)
        v = _map_pie(D0, np.asarray(D1, dtype=float)) @ expm(D0 * t0)
        Pstep = expm(D0 * step)
        out = np.empty(n)
        for i in range(n):
            out[i] = 1.0 - float(v.sum())
            v = v @ Pstep
        return out

    def _smp_stage_probabilities(self):
        """Semi-Markov stage probabilities for method 'smp' (JAR SolverENV.init, MATLAB init).

        P(k,e) = int dF_ke(t) prod_{h!=k,e} (1 - F_kh(t)) is the embedded jump chain, integrated on
        N = max(1000, 100T) intervals with the survival at the midpoint, T doubling from 1 until
        F_ke(T) >= 1 - 1e-8. The mean holding time integrates the sojourn survival
        prod_{h!=k} (1 - F_kh(t)) by composite Simpson (N = 10000) over [0, U], U doubling from 10
        until the survival is <= 1e-8. Then probEnv(k) is proportional to pie_dtmc(k) * hold(k),
        and probOrig(k,e) = probEnv(k) E0(k,e) / sum_{h!=e} probEnv(h) E0(h,e)."""
        from .api.mc.dtmc import dtmc_solve
        E = self.env_model.num_stages
        env = self.env_model.env
        arc = lambda k, h: env[k][h] if k < len(env) and h < len(env[k]) else None
        self._env_generator()
        eps = 1e-8
        dtmcP = np.zeros((E, E))
        for k in range(E):
            for e in range(E):
                if k == e or arc(k, e) is None:
                    continue
                T = 1.0
                while self._arc_cdf_grid(arc(k, e), T, 0.0, 1)[0] < 1.0 - eps:
                    T *= 2.0
                    if T > 1e6:
                        break
                N = max(1000, int(round(T * 100)))
                dt = T / N
                F = self._arc_cdf_grid(arc(k, e), 0.0, dt, N + 1)
                surv = np.ones(N)
                for h in range(E):
                    if h != k and h != e and arc(k, h) is not None:
                        surv *= 1.0 - self._arc_cdf_grid(arc(k, h), dt / 2.0, dt, N)
                dtmcP[k, e] = float(np.sum(np.diff(F) * surv))
        self.dtmcP = dtmcP
        pie = np.asarray(dtmc_solve(dtmcP), dtype=float).ravel()

        hold = np.zeros(E)
        Nh = 10000
        for k in range(E):
            others = [h for h in range(E) if h != k and arc(k, h) is not None]

            def sojourn_surv(t0, step, n):
                s = np.ones(n)
                for h in others:
                    s *= 1.0 - self._arc_cdf_grid(arc(k, h), t0, step, n)
                return s

            U = 10.0
            while sojourn_surv(U, 0.0, 1)[0] > eps:
                U *= 2.0
                if U > 1e6:
                    break
            s = sojourn_surv(0.0, U / (2 * Nh), 2 * Nh + 1)
            hold[k] = float((s[0:-1:2] + 4.0 * s[1::2] + s[2::2]).sum() * (U / Nh) / 6.0)
        self.holdTimeMatrix = hold

        pi = pie * hold
        pi = pi / pi.sum()
        self.env_model.probEnv = pi
        E0 = self.E0
        emb = np.zeros((E, E))
        for e in range(E):
            denom = sum(pi[h] * E0[h, e] for h in range(E) if h != e)
            if denom > 0:
                for k in range(E):
                    if k != e:
                        emb[k, e] = pi[k] * E0[k, e] / denom
        self.env_model.probOrig = emb


    def _solve_env_limit(self, method: str):
        """Closed-form fast/slow random-environment limits (method 'avg'/'dec').

        These treat the environment (stage) process as either infinitely fast or
        infinitely slow relative to the base-model dynamics and therefore need no
        inter-stage coupling iteration, only steady-state stage solves:

        - 'avg' (fast-environment limit): the base model sees the
          stationary-probability-weighted average of the modulated rates. A
          single rate-averaged model is built and solved once. Exact as the
          stage-switching rate -> Inf.
        - 'dec' (slow-environment / quasi-stationary decomposition): each stage
          is solved independently in steady state and the per-stage metrics are
          averaged with weights probEnv(e). Exact as the stage-switching rate -> 0.

        Mirrors matlab @SolverENV/SolverENV.m solveEnvLimit().
        """
        from .api.sn.getters import sn_get_node_arvr_from_tput

        self.init()
        E = self.getNumberOfModels()
        prob_env = np.asarray(self.env_model.probEnv, dtype=float).ravel()

        def _stage_avg(solver):
            out = solver.getAvg()
            Q, U, T = out[0], out[1], out[3]
            return np.asarray(Q, dtype=float), np.asarray(U, dtype=float), np.asarray(T, dtype=float)

        # A Cache carries its result on the NODE, not in the Q/U/T tables, so an
        # aggregate that only sums those tables leaves getHitRatio() at None and
        # the caller reads a missing number as a missing class. Mirrors MATLAB
        # solveEnvLimit, which accumulates the weighted hit/miss ratios here and
        # writes them onto the stage-1 reference model.
        ref_nodes = self.ensemble[0].get_nodes()
        cache_idx = [c for c, nd in enumerate(ref_nodes)
                     if type(nd).__name__ == 'Cache' or hasattr(nd, 'set_result_hit_prob')]
        K = len(self.ensemble[0].get_classes())
        acc = {c: {'arv': None, 'hit': None, 'miss': None, 'dhit': None, 'hitl': None}
               for c in cache_idx}

        def _cache_arrival_row(model, TN, c):
            """Per-class arrival rate into cache node `c`, the blend weight.

            Read from the node table, not from the Source throughput, because
            the latter is 0 on a closed cache model."""
            arv = np.zeros(K)
            ANn = sn_get_node_arvr_from_tput(model.getStruct(), np.asarray(TN, dtype=float))
            if ANn is None or ANn.size == 0 or ANn.shape[0] <= c:
                return arv
            row = np.asarray(ANn[c, :], dtype=float).ravel()
            n = min(K, row.size)
            arv[:n] = row[:n]
            arv[~np.isfinite(arv)] = 0.0
            return arv

        def _accum_rate(a, ratio, arv):
            """One quantity of the accumulator. An ABSENT ratio contributes
            nothing and leaves the accumulator absent, so a metric no stage
            measured stays empty rather than becoming a fabricated zero."""
            if ratio is None:
                return a
            ratio = np.asarray(ratio, dtype=float).ravel()
            if ratio.size == 0:
                return a
            ratio = np.where(np.isfinite(ratio), ratio, 0.0)
            if a is None:
                a = np.zeros(ratio.size)
            n = min(ratio.size, arv.size, a.size)
            a[:n] += arv[:n] * ratio[:n]
            return a

        def _accum_cache(model, TN, w):
            """Accumulate the ARRIVAL-WEIGHTED cache RATES of `model`'s caches.

            Rates, not ratios. A stage's hit ratio is conditional on arriving in
            that stage, so a probEnv-weighted mean of ratios disagrees with the
            Sink-row hit throughput, sum_e p_e*lambda_e*h_e, in the SAME node
            table. Accumulating rates and dividing once keeps the two
            consistent, and is what the statevec and mean-field couplings do.
            Mirrors MATLAB SolverENV.accumCacheMetric."""
            nodes = model.get_nodes()
            for c in cache_idx:
                if c >= len(nodes):
                    continue
                node = nodes[c]
                arv = _cache_arrival_row(model, TN, c) * w
                if acc[c]['arv'] is None:
                    acc[c]['arv'] = np.zeros(K)
                acc[c]['arv'] += arv
                for key, getter in (('hit', 'get_hit_ratio'),
                                    ('miss', 'get_miss_ratio'),
                                    ('dhit', 'get_delayed_hit_ratio')):
                    fn = getattr(node, getter, None)
                    acc[c][key] = _accum_rate(acc[c][key], fn() if fn is not None else None, arv)
                fn = getattr(node, 'get_hit_ratio_by_list', None)
                hl = fn() if fn is not None else None
                if hl is None:
                    continue
                hl = np.atleast_2d(np.asarray(hl, dtype=float))
                if hl.size == 0:
                    continue
                hl = np.where(np.isfinite(hl), hl, 0.0)
                if acc[c]['hitl'] is None:
                    acc[c]['hitl'] = np.zeros_like(hl)
                wrow = np.zeros(hl.shape[0])
                n = min(hl.shape[0], arv.size)
                wrow[:n] = arv[:n]
                acc[c]['hitl'] += wrow[:, None] * hl

        def _populate_node_metrics(solver):
            """The cache metrics ride the NODE table, which getAvg does not build."""
            fn = getattr(solver, 'getAvgNodeTable', None) or getattr(solver, 'get_avg_node_table', None)
            if fn is not None:
                fn()

        if method == 'dec':
            Qval = Uval = Tval = None
            for e in range(E):
                solver = self.getSolver(e)
                if hasattr(solver, 'reset'):
                    solver.reset()
                Qe, Ue, Te = _stage_avg(solver)
                if Qval is None:
                    Qval = np.zeros_like(Qe)
                    Uval = np.zeros_like(Ue)
                    Tval = np.zeros_like(Te)
                Qval += prob_env[e] * Qe
                Uval += prob_env[e] * Ue
                Tval += prob_env[e] * Te
                if cache_idx:
                    _populate_node_metrics(solver)
                    _accum_cache(self.ensemble[e], Te, prob_env[e])
        else:
            avg_model = self._build_rate_averaged_model(prob_env)
            template = self.getSolver(0)
            inner = type(template)(avg_model, template.options) if hasattr(template, 'options') \
                else type(template)(avg_model)
            Qval, Uval, Tval = _stage_avg(inner)
            if cache_idx:
                _populate_node_metrics(inner)
                _accum_cache(avg_model, Tval, 1.0)

        # Divide the accumulated rates by the accumulated arrivals, ONCE. A class
        # with no arrival anywhere divides by NaN and STAYS NaN: its ratio is
        # undefined, and 0 would read as "never hits".
        cache_records = []
        for c in cache_idx:
            node = ref_nodes[c]
            arv = acc[c]['arv']
            if arv is None:
                continue
            d = np.where(arv > 0, arv, np.nan)
            rec = {'name': getattr(node, 'name', None), 'node': c,
                   'hitprob': None, 'missprob': None,
                   'delayedhitprob': None, 'hitproblist': None}
            for key, field, setter in (('hit', 'hitprob', 'set_result_hit_prob'),
                                       ('miss', 'missprob', 'set_result_miss_prob'),
                                       ('dhit', 'delayedhitprob', 'set_result_delayed_hit_prob')):
                a = acc[c][key]
                if a is None:
                    continue
                v = np.full(a.size, np.nan)
                n = min(a.size, d.size)
                v[:n] = a[:n] / d[:n]
                rec[field] = v
                fn = getattr(node, setter, None)
                if fn is not None:
                    fn(v)
            hl = acc[c]['hitl']
            if hl is not None:
                wrow = np.full(hl.shape[0], np.nan)
                n = min(hl.shape[0], d.size)
                wrow[:n] = d[:n]
                v = hl / wrow[:, None]
                rec['hitproblist'] = v
                fn = getattr(node, 'set_result_hit_prob_list', None)
                if fn is not None:
                    fn(v)
            cache_records.append(rec)

        self._cache_records = cache_records
        self.result = {'Avg': {'Q': Qval, 'U': Uval, 'T': Tval, 'Cache': cache_records}}

    def _build_rate_averaged_model(self, prob_env):
        """Fast-environment model: every stage-varying station rate replaced by
        its probEnv-weighted average, as an exponential. Non-modulated parameters
        keep their original distribution."""
        from .distributions import Exp

        E = self.getNumberOfModels()
        sn = [m.getStruct() for m in self.ensemble]
        avg_model = self.ensemble[0].copy()
        stations = avg_model.get_stations()
        classes = avg_model.get_classes()
        M = len(stations)
        K = len(classes)
        for i in range(M):
            node = stations[i]
            # a stateful or absorbing node carries no service rate to average
            if type(node).__name__ in ('Cache', 'Sink'):
                continue
            for k in range(K):
                rates = np.array([sn[e].rates[i, k] for e in range(E)], dtype=float)
                if np.any(np.isnan(rates)) or np.any(rates <= 0):
                    continue  # disabled for some stage: leave as configured
                if (rates.max() - rates.min()) <= 1e-12 * max(1.0, rates.max()):
                    continue  # not modulated: keep the original distribution
                ravg = float(np.dot(prob_env, rates))
                if callable(getattr(node, 'set_arrival', None)) and type(node).__name__ == 'Source':
                    node.set_arrival(classes[k], Exp(ravg))
                elif callable(getattr(node, 'set_service', None)):
                    node.set_service(classes[k], Exp(ravg))
        # The copy inherited a cached NetworkStruct from the stage model; force a
        # hard rebuild so the averaged rates take effect.
        avg_model.refresh_struct()
        return avg_model

    def pre(self, it: int):
        """Pre-iteration operations.

        At iteration 1, initialize each stage model from steady-state queue lengths.
        Mirrors MATLAB SolverENV.pre() which calls ensemble{e}.initFromMarginal(QN).
        """
        E = self.getNumberOfModels()

        if it == 1:
            for e in range(E):
                solver = self.getSolver(e)
                model = self.ensemble[e] if e < len(self.ensemble) else None
                if solver is not None and model is not None:
                    try:
                        # Check solver timespan to decide steady-state vs transient
                        timespan_inf = True
                        if hasattr(solver, 'options') and hasattr(solver.options, 'timespan'):
                            ts = solver.options.timespan
                            if ts is not None and len(ts) > 1 and np.isfinite(ts[1]):
                                timespan_inf = False
                        elif hasattr(solver, '_native_solver') and hasattr(solver._native_solver, 'options'):
                            opts = solver._native_solver.options
                            if hasattr(opts, 'timespan') and opts.timespan is not None:
                                ts = opts.timespan
                                if len(ts) > 1 and np.isfinite(ts[1]):
                                    timespan_inf = False

                        if timespan_inf:
                            # Steady-state: get QN from getAvg()
                            result = solver.getAvg()
                            if result is not None and result[0] is not None:
                                QN = np.atleast_2d(result[0])
                            else:
                                continue
                        else:
                            # Transient: get QN from transient analysis final values
                            QNt, UNt, TNt = solver.getTranAvg()
                            if QNt is not None:
                                M = len(QNt)
                                K = len(QNt[0]) if M > 0 else 0
                                QN = np.zeros((M, K))
                                for i in range(M):
                                    for k in range(K):
                                        t_vals, metric_vals = _get_tran_data(QNt[i][k])
                                        if metric_vals is not None and len(metric_vals) > 0:
                                            QN[i, k] = metric_vals[-1]
                            else:
                                continue

                        # see _kb/06-solver-catalog.md (ENV: Python environment.py additional notes) for rationale
                        solver_name = type(solver).__name__
                        is_lqn = type(model).__name__ == 'LayeredNetwork'
                        if 'Fluid' not in solver_name and 'FLD' not in solver_name and not is_lqn:
                            QN = self._round_marginal_for_discrete_solver(QN, e)

                        # Initialize model state from queue lengths (MATLAB: self.ensemble{e}.initFromMarginal(QN))
                        if hasattr(model, 'initFromMarginal'):
                            model.initFromMarginal(QN)
                        elif hasattr(model, 'init_from_marginal'):
                            model.init_from_marginal(QN)

                        # The model state alone does not reach a fluid stage;
                        # see the same call in post(). Without it the FIRST
                        # sweep integrates from the empty state instead of the
                        # seed just computed.
                        if hasattr(solver, 'setInitialState'):
                            solver.setInitialState(np.asarray(QN, dtype=float))
                    except Exception:
                        pass

    def _stage_steady(self, e: int, M: int, K: int):
        """The stage's own STEADY-STATE tables, the exit value of any metric for
        which its transient carries no trajectory.

        A METRIC WITH NO TRAJECTORY IS NOT A METRIC WORTH ZERO. A stage solver
        reports a transient only for the quantities its integration actually
        advances; what it leaves out it holds CONSTANT over the stage, so the
        value at the exit instant is that constant whatever the sojourn was.
        Seeding the exit tables with zeros instead and blending the zeros is
        what zeroes a whole answer whenever the coupling cannot read a cell:
        the EXT skip below, a metric the stage solver omitted, a stage that
        produced no transient at all.

        Free in the common case: the same stage solve that produced the
        transient already set the averages, so getAvg returns them cached.
        """
        Qs = np.zeros((M, K))
        Us = np.zeros((M, K))
        Ts = np.zeros((M, K))
        solver = self.getSolver(e)
        if solver is None or not hasattr(solver, 'getAvg'):
            return Qs, Us, Ts
        # A LAYERED STAGE ANSWERS getAvg IN A DIFFERENT INDEX SPACE than the
        # one this fills. (M, K) here is the coupling's BLOCK-DIAGONAL
        # station x class aggregate -- what get_tran_handles lays out and what
        # the exit blend is written in -- whereas SolverLN.getAvg is
        # getEnsembleAvg: one row over the LQN's OWN nodes (hosts, tasks,
        # entries, activities), carrying no station or class meaning at all.
        # The leading-block copy below cannot tell the two apart, so it read
        # LQN node k as class k and pasted an activity's queue length onto an
        # off-block (station, class) cell. Off-block is precisely where the
        # transient overwrites nothing, so the stray value survived the
        # probEnv blend and broke the closed-population conservation the
        # stages satisfy. getBlockAvg is the steady table in the RIGHT space.
        model = self.ensemble[e] if e < len(self.ensemble) else None
        if type(model).__name__ == 'LayeredNetwork' and hasattr(solver, 'getBlockAvg'):
            try:
                Qb, Ub, Tb = solver.getBlockAvg()
            except Exception:
                return Qs, Us, Ts
            _finite_block(Qb, Qs)
            _finite_block(Ub, Us)
            _finite_block(Tb, Ts)
            return Qs, Us, Ts
        try:
            res = solver.getAvg()
        except Exception:
            return Qs, Us, Ts
        if res is None or len(res) < 4:
            return Qs, Us, Ts
        _finite_block(res[0], Qs)
        _finite_block(res[1], Us)
        _finite_block(res[3], Ts)
        return Qs, Us, Ts

    def analyze(self, it: int, e: int) -> Tuple[Dict, float]:
        """
        Analyze stage e at iteration it using transient analysis.

        Mirrors MATLAB SolverENV.analyze():
          [Qt,Ut,Tt] = self.ensemble{e}.getTranHandles;
          self.solvers{e}.reset();
          [QNt,UNt,TNt] = self.solvers{e}.getTranAvg(Qt,Ut,Tt);

        Returns transient analysis results in MATLAB-compatible format.
        """
        import time
        t0 = time.time()

        result_e = {
            'Tran': {
                'Avg': {
                    'Q': None,
                    'U': None,
                    'T': None
                }
            }
        }

        solver = self.getSolver(e)
        if solver is None:
            return result_e, 0.0

        # The resolved horizon is set for THIS stage solve and put back after
        # it, as MATLAB's analyze_ does: it lands after pre(), which reads a
        # non-finite horizon as "start from the steady state", and before
        # getTranAvg, which would otherwise resolve it by the fluid handler's
        # own rule; and leaving it on the solver would change what a later
        # per-stage getter answers.
        ts_saved = None
        ts_stage = self._stage_horizon(e)
        if ts_stage is not None:
            opts = getattr(solver, 'options', None)
            ts_saved = opts.get('timespan') if isinstance(opts, dict) else getattr(opts, 'timespan', None)
            self._set_stage_timespan(e, ts_stage)

        # ASK FOR THE POINTS THE SOJOURN WEIGHT NEEDS, as MATLAB's analyze_
        # does. Reading them off a linear interpolation of the integrator's own
        # grid (what post() still does as a fallback) is a first-order error
        # that the ODE's continuous extension does not have. A stage solver
        # that ignores the request is resampled in post() exactly as before.
        grid = self._cdf_grid_for_stage(e)
        if grid is not None:
            opts = getattr(solver, 'options', None)
            if isinstance(opts, dict):
                opts['tranpoints'] = grid
            elif opts is not None:
                opts.tranpoints = grid

        # Reset solver so it re-reads model state (MATLAB: self.solvers{e}.reset())
        if hasattr(solver, 'reset'):
            solver.reset()

        # see _kb/06-solver-catalog.md (ENV: Python environment.py additional notes) for rationale
        try:
            if hasattr(solver, 'getTranAvg'):
                QNt, UNt, TNt = solver.getTranAvg()
                result_e['Tran']['Avg']['Q'] = QNt
                result_e['Tran']['Avg']['U'] = UNt
                result_e['Tran']['Avg']['T'] = TNt
                # 'meancov' needs the WITHIN-STAGE covariance on top of the
                # means. It costs a second integration of the same stage, which
                # is why only this method asks for it; a stage solver that
                # carries a first moment only says so through
                # supportsTransientVariance and contributes the timing variance
                # of the sojourn alone. QCov['C'] is (M*K, M*K, nt) indexed
                # ir = r*M + i, so post() needs no phase layout.
                if self._is_meancov() and getattr(
                        solver, 'supportsTransientVariance', None) is not None \
                        and solver.supportsTransientVariance():
                    tvar, _, _, QCovt = solver.getTranAvgVar()
                    if QCovt is not None and tvar is not None:
                        result_e['Tran']['Avg']['QCov'] = {
                            't': np.asarray(tvar, dtype=float).ravel(),
                            'C': np.asarray(QCovt, dtype=float)}
            else:
                # Fall back to steady-state wrapped in transient format
                result = solver.getAvg() if hasattr(solver, 'getAvg') else None
                if result is not None and result[0] is not None:
                    QN = np.atleast_2d(result[0])
                    UN = np.atleast_2d(result[1]) if result[1] is not None else np.zeros_like(QN)
                    TN = np.atleast_2d(result[3]) if len(result) > 3 and result[3] is not None else np.zeros_like(QN)
                else:
                    QN = np.array([[0.0]])
                    UN = np.array([[0.0]])
                    TN = np.array([[0.0]])
                M, K = QN.shape
                t_vals = np.array([0.0, 1000.0])
                result_e['Tran']['Avg']['Q'] = [[{'t': t_vals, 'metric': np.array([QN[i, r], QN[i, r]])} for r in range(K)] for i in range(M)]
                result_e['Tran']['Avg']['U'] = [[{'t': t_vals, 'metric': np.array([UN[i, r], UN[i, r]])} for r in range(K)] for i in range(M)]
                result_e['Tran']['Avg']['T'] = [[{'t': t_vals, 'metric': np.array([TN[i, r], TN[i, r]])} for r in range(K)] for i in range(M)]
        except Exception as ex:
            # A STAGE THAT DID NOT SOLVE IS NOT A STAGE THAT SOLVED TO ZERO.
            # Swallowing this left result_e carrying None metrics, so the
            # ensemble finished with self.result unset and avg_table() then
            # failed on an unpack of None -- an opaque TypeError in place of the
            # stage solver's own diagnosis. MATLAB's analyze_ lets the stage
            # error out, which is what lets SolverAUTO's ENV pool fall through
            # from a stage solver that cannot run the transient to one that can.
            if ts_saved is not None:
                self._set_stage_timespan(e, ts_saved)
            raise RuntimeError(
                "SolverENV: stage %d could not be analyzed by %s: %s"
                % (e, type(solver).__name__, ex)) from ex

        if ts_saved is not None:
            self._set_stage_timespan(e, ts_saved)

        runtime = time.time() - t0
        return result_e, runtime

    def post(self, it: int):
        """Post-iteration operations - compute exit metrics and update entry marginals.

        Mirrors JAR SolverENV.post():
        - Computes CDF-weighted exit metrics for each stage-to-stage transition
          using MAP CDF (map_cdf(D0, D1, t)) from proc[e][h]
        - Computes entry marginals using probOrig and reset functions
        - Normalizes QEntry to preserve closed chain populations
        - Updates state-dependent environment rates if configured
        """
        E = self.getNumberOfModels()

        # Identify EXT (Source) stations to skip in QExit computation
        # (matching JAR which skips SchedStrategy.EXT stations)
        isExtStation = [False] * 100  # will resize
        model0 = self.ensemble[0] if len(self.ensemble) > 0 else None
        if model0 is not None:
            sn0 = model0.get_struct() if hasattr(model0, 'get_struct') else None
            if sn0 is not None and hasattr(sn0, 'sched'):
                stations = model0.get_stations() if hasattr(model0, 'get_stations') else []
                isExtStation = [False] * len(stations)
                for idx, st in enumerate(stations):
                    if hasattr(sn0, 'sched') and isinstance(sn0.sched, dict):
                        sched = sn0.sched.get(st, None)
                    elif hasattr(sn0, 'sched') and hasattr(sn0.sched, '__getitem__'):
                        try:
                            sched = sn0.sched[idx] if idx < len(sn0.sched) else None
                        except (TypeError, KeyError):
                            sched = None
                    else:
                        sched = None
                    if sched is not None:
                        sched_name = sched.name if hasattr(sched, 'name') else str(sched)
                        if 'EXT' in sched_name.upper():
                            isExtStation[idx] = True

        # Compute exit metrics for each stage-to-stage transition
        Qexit = {}
        Uexit = {}
        Texit = {}

        for e in range(E):
            result_e = self.results.get((it, e))
            if result_e is None or result_e['Tran']['Avg']['Q'] is None:
                continue

            Q_tran = result_e['Tran']['Avg']['Q']
            U_tran = result_e['Tran']['Avg']['U']
            T_tran = result_e['Tran']['Avg']['T']

            M = len(Q_tran)
            K = len(Q_tran[0]) if M > 0 else 0
            # THE SEED IS THE STAGE'S OWN STEADY STATE, NEVER ZERO; see
            # _stage_steady. Every cell the weighted average below does not
            # reach -- an EXT station, a metric with no trajectory -- keeps the
            # value the stage holds constant over its sojourn, which IS its
            # exit value.
            Q_st, U_st, T_st = self._stage_steady(e, M, K)

            # The handoff averages over the SOJOURN, not over the e->h clock:
            # competing exponentials leave the exit time independent of the
            # destination, so the weight is holdTime[e] for every h and does not
            # depend on h. see _kb/06-solver-catalog.md (ENV meanfield)
            hold_mmap = self.env_model.holdTime[e]
            D0_e = np.asarray(hold_mmap[0], dtype=float)
            D1_e = np.asarray(hold_mmap[1], dtype=float)

            # ONE PASS FOR EVERY DESTINATION: the weight below reads
            # holdTime[e] and the (i, r) trajectory, so the exit averages are
            # the same matrix for all E destinations, exactly as the note above
            # says. Built inside the h loop they evaluated the sojourn CDF E
            # times per (i, r) per sweep for one answer, and
            # _interpolate_for_cdf_map hands the CDF a 5000-point grid, so a
            # 24-stage environment paid 24x for every cell of every sweep.
            Q_ex = Q_st.copy()
            U_ex = U_st.copy()
            T_ex = T_st.copy()

            for i in range(M):
                if i < len(isExtStation) and isExtStation[i]:
                    continue  # Skip Source stations
                for r in range(K):
                    Qir = Q_tran[i][r]
                    t_vals, metric_vals = _get_tran_data(Qir)
                    if t_vals is None or len(t_vals) < 2:
                        continue

                    # Interpolate transient data onto a finer grid for
                    # accurate CDF weighting
                    t_fine, q_fine, u_fine, t_fine_data = \
                        _interpolate_for_cdf_map(
                            t_vals, metric_vals,
                            U_tran[i][r], T_tran[i][r], D0_e, D1_e)

                    # Use MAP CDF for weighting (matching JAR)
                    cdf_vals = _map_eval_cdf(D0_e, D1_e, t_fine)
                    w = np.zeros(len(t_fine))
                    w[1:] = cdf_vals[1:] - cdf_vals[:-1]

                    w_sum = np.sum(w)
                    if w_sum > 0 and not np.any(np.isnan(w)):
                        Q_ex[i, r] = np.dot(q_fine, w) / w_sum
                        if u_fine is not None:
                            U_ex[i, r] = np.dot(u_fine, w) / w_sum
                        if t_fine_data is not None:
                            T_ex[i, r] = np.dot(t_fine_data, w) / w_sum
                    else:
                        # Fall back to final value
                        Q_ex[i, r] = metric_vals[-1] if len(metric_vals) > 0 else 0.0

            # A destination the environment cannot reach keeps the stage's own
            # steady state, which is the only thing the h loop ever decided.
            # COPIED per destination because the entry pass below hands
            # Qexit[(h, e)] to a reset function free to write into it.
            for h in range(E):
                if self.env_model.proc[e][h] is None:
                    Qexit[(e, h)] = Q_st.copy()
                    Uexit[(e, h)] = U_st.copy()
                    Texit[(e, h)] = T_st.copy()
                    continue
                Qexit[(e, h)] = Q_ex.copy()
                Uexit[(e, h)] = U_ex.copy()
                Texit[(e, h)] = T_ex.copy()

        # Store exit metrics for convergence check
        self._Qexit = Qexit

        # 'meancov' carries a COVARIANCE beside the mean. Independent of the
        # DESTINATION, exactly as Qexit is: the sojourn weight reads holdTime[e]
        # alone, because competing exponentials leave the exit time independent
        # of where the switch goes.
        meancov = self._is_meancov()
        Cexit = {}
        Mdim, Kdim = 0, 0
        if meancov:
            for e in range(E):
                res_e = self.results.get((it, e))
                if res_e is not None and res_e['Tran']['Avg']['Q'] is not None:
                    Mdim = len(res_e['Tran']['Avg']['Q'])
                    Kdim = len(res_e['Tran']['Avg']['Q'][0]) if Mdim > 0 else 0
                    break
            for e in range(E):
                Cexit[e] = self._stage_exit_cov(it, e, Mdim, Kdim)

        # Compute entry marginals using reset functions
        for e in range(E):
            result_e = self.results.get((it, e))
            if result_e is None or result_e['Tran']['Avg']['Q'] is None:
                continue

            Q_tran = result_e['Tran']['Avg']['Q']
            M = len(Q_tran)
            K = len(Q_tran[0]) if M > 0 else 0
            Qentry = np.zeros((M, K))
            Sentry = np.zeros((M * K, M * K))  # second moment about zero
            mentry = np.zeros(M * K)

            for h in range(E):
                if (h, e) in Qexit and self.env_model.probOrig[h, e] > 0:
                    reset_fn = self.env_model.resetFun[h][e]
                    mh = np.asarray(reset_fn(Qexit[(h, e)]), dtype=float)
                    Qentry += self.env_model.probOrig[h, e] * mh
                    if meancov and Mdim == M and Kdim == K:
                        # The reset is an arbitrary map on the means, so the
                        # covariance crosses it by the delta method. Then the
                        # MIXTURE over origins is taken on the SECOND MOMENT,
                        # not on the covariances: a convex combination of the
                        # C_h alone drops the spread of the per-origin means,
                        # which is most of the variance when the stages differ.
                        Rhe = _reset_jacobian(reset_fn, Qexit[(h, e)], M, K)
                        Che = Rhe @ Cexit.get(h, np.zeros((M * K, M * K))) @ Rhe.T
                        mhv = mh.ravel(order='F')
                        p_he = self.env_model.probOrig[h, e]
                        mentry += p_he * mhv
                        Sentry += p_he * (Che + np.outer(mhv, mhv))
            Centry = None
            if meancov and Mdim == M and Kdim == K:
                Centry = Sentry - np.outer(mentry, mentry)
                Centry = 0.5 * (Centry + Centry.T)  # drop the rounding asymmetry

            # Normalize QEntry to preserve closed chain populations
            # (matching JAR SolverENV.post() lines 1096-1119)
            sn_ref = None
            model_e = self.ensemble[e] if e < len(self.ensemble) else None
            if model_e is not None:
                sn_ref = model_e.get_struct() if hasattr(model_e, 'get_struct') else None
            if sn_ref is not None and hasattr(sn_ref, 'chains') and hasattr(sn_ref, 'njobs'):
                nchains = sn_ref.nchains if hasattr(sn_ref, 'nchains') else 0
                chains = np.atleast_2d(sn_ref.chains)
                njobs = np.atleast_1d(sn_ref.njobs).flatten()
                for c in range(nchains):
                    chain_classes = np.where(chains[c, :] > 0)[0]
                    njobs_chain = sum(njobs[k] for k in chain_classes)
                    if np.isinf(njobs_chain):
                        continue  # Open chain
                    state_chain = sum(Qentry[i, k] for i in range(M) for k in chain_classes)
                    if state_chain > 0 and abs(state_chain - njobs_chain) > 1e-10:
                        scale = njobs_chain / state_chain
                        for i in range(M):
                            for k in chain_classes:
                                Qentry[i, k] *= scale

            # Initialize model state from entry marginals and reset solver
            # MATLAB: self.solvers{e}.reset(); self.ensemble{e}.initFromMarginal(Qentry{e});
            model = self.ensemble[e] if e < len(self.ensemble) else None
            solver = self.getSolver(e)
            if solver is not None:
                if hasattr(solver, 'reset'):
                    solver.reset()
            if model is not None:
                # Round fractional marginals for discrete solvers. Skip for
                # LayeredNetwork stages (see the pre() rationale).
                solver_name = type(solver).__name__ if solver is not None else ''
                is_lqn = type(model).__name__ == 'LayeredNetwork'
                if 'Fluid' not in solver_name and 'FLD' not in solver_name and not is_lqn:
                    Qentry = self._round_marginal_for_discrete_solver(Qentry, e)

                if hasattr(model, 'initFromMarginal'):
                    model.initFromMarginal(Qentry)
                elif hasattr(model, 'init_from_marginal'):
                    model.init_from_marginal(Qentry)

                # The model state alone does not reach a fluid stage: SolverFLD
                # takes its ODE initial condition from options.init_sol, so
                # without this the transient starts at the stage's OWN steady
                # state, sits flat, and the exit average collapses to the
                # quasi-stationary blend regardless of the coupling.
                if solver is not None and hasattr(solver, 'setInitialState'):
                    solver.setInitialState(np.asarray(Qentry, dtype=float))

                # The second half of the handoff. initFromMarginal carries the
                # mean through the MODEL; a covariance has no place on a model
                # object, so it is handed to the stage SOLVER, in the same
                # (station, class) index space, together with the mean it belongs
                # to. init_qlen is not redundant with initFromMarginal: the 'kp'
                # method has its own state layout and deliberately does not read
                # options.init_sol, so without it a 'kp' stage would restart
                # empty every iteration.
                if Centry is not None and solver is not None:
                    _set_stage_config(solver, 'init_qlen',
                                      np.asarray(Qentry, dtype=float))
                    _set_stage_config(solver, 'init_qcov', Centry)

    def _is_meancov(self):
        """Whether this run carries a COVARIANCE beside the mean across a switch.

        One predicate, read by analyze(), post() and finish(), so the three
        cannot drift.
        """
        method = self.options.get('method', 'default') if isinstance(self.options, dict) \
            else getattr(self.options, 'method', 'default')
        return str(method).lower() in ('meancov', 'env.meancov')

    def _stage_exit_cov(self, it, e, M, K):
        """Covariance of the station-class queue lengths when stage E is left.

        The law of total variance over the random sojourn T:

            Cov[Q(T)] = E_T[Cov(Q(t)|t)] + Cov_T[E(Q(t)|t)]
                      = sum_n w_n (C(t_n) + m(t_n)m(t_n)')/sum(w) - m_ex m_ex'

        with the SAME weights the exit mean uses, so the two are consistent by
        construction. C(t) is the within-stage covariance the stage solver
        integrated and is absent -- hence zero -- for a stage solver carrying a
        first moment only; the timing term survives regardless.
        """
        n = M * K
        C = np.zeros((n, n))
        if n == 0:
            return C
        res = self.results.get((it, e))
        if res is None or res['Tran']['Avg']['Q'] is None:
            return C
        Q_tran = res['Tran']['Avg']['Q']
        hold_mmap = self.env_model.holdTime[e]
        D0 = np.asarray(hold_mmap[0], dtype=float)
        D1 = np.asarray(hold_mmap[1], dtype=float)

        # ONE grid for the whole stage. The exit moment multiplies station-class
        # means by one another, so they have to be read at the same instants; the
        # per-metric grids the exit MEAN uses are each self-contained and need no
        # such alignment.
        t_base = None
        for i in range(min(M, len(Q_tran))):
            for r in range(min(K, len(Q_tran[i]))):
                tv, _ = _get_tran_data(Q_tran[i][r])
                if tv is not None and len(tv) > 1 \
                        and (t_base is None or len(tv) > len(t_base)):
                    t_base = np.asarray(tv, dtype=float)
        if t_base is None or t_base.size < 2:
            return C
        alpha = _map_prob(D0, D1)
        total_rate = float(alpha @ D1 @ np.ones(D0.shape[0]))
        mean_sojourn = 1.0 / total_rate if total_rate > 0 \
            else (t_base[-1] - t_base[0]) / 10.0
        t_fine = _cdf_refine_grid(t_base, mean_sojourn, 5000)
        cdf_vals = _map_eval_cdf(D0, D1, t_fine)
        w = np.zeros(len(t_fine))
        w[1:] = cdf_vals[1:] - cdf_vals[:-1]
        w_sum = float(np.sum(w))
        if not (w_sum > 0) or np.any(np.isnan(w)):
            return C

        Mtraj = np.zeros((n, t_fine.size))
        for i in range(min(M, len(Q_tran))):
            for r in range(min(K, len(Q_tran[i]))):
                tv, mv = _get_tran_data(Q_tran[i][r])
                if tv is None or len(tv) < 2:
                    continue
                tv = np.asarray(tv, dtype=float)
                mv = np.asarray(mv, dtype=float)
                Mtraj[r * M + i, :] = np.interp(np.clip(t_fine, tv[0], tv[-1]), tv, mv)
        mexit = (Mtraj @ w) / w_sum
        S = ((Mtraj * w) @ Mtraj.T) / w_sum
        Ctraj = _interp_cov(res['Tran']['Avg'].get('QCov'), t_fine, n)
        if Ctraj is not None:
            S = S + (Ctraj @ w) / w_sum
        C = S - np.outer(mexit, mexit)
        return 0.5 * (C + C.T)

    def _round_marginal_for_discrete_solver(self, Q, stage_idx):
        """Round fractional queue lengths to integers using the largest remainder method,
        preserving closed chain populations exactly."""
        from .api.state.marginal import roundMarginalPreservingChains
        model = self.ensemble[stage_idx] if stage_idx < len(self.ensemble) else None
        if model is None:
            return Q
        sn = model.get_struct() if hasattr(model, 'get_struct') else None
        if sn is None:
            return Q
        return roundMarginalPreservingChains(Q, sn)

    def converged(self, it: int) -> bool:
        """Check if iteration has converged."""
        if it <= 1:
            return False

        E = self.getNumberOfModels()
        iter_tol = self.options.get('iter_tol', 1e-4)

        # Compare queue lengths between iterations
        for e in range(E):
            result_curr = self.results.get((it, e))
            result_prev = self.results.get((it - 1, e))

            if result_curr is None or result_prev is None:
                return False

            Q_curr = result_curr['Tran']['Avg']['Q']
            Q_prev = result_prev['Tran']['Avg']['Q']

            if Q_curr is None or Q_prev is None:
                return False

            M = len(Q_curr)
            K = len(Q_curr[0]) if M > 0 else 0

            for i in range(M):
                for k in range(K):
                    curr_val = 0.0
                    prev_val = 0.0

                    Qik_curr = Q_curr[i][k]
                    _, curr_metric = _get_tran_data(Qik_curr)
                    if curr_metric is not None and len(curr_metric) > 0:
                        curr_val = curr_metric[0]

                    Qik_prev = Q_prev[i][k]
                    _, prev_metric = _get_tran_data(Qik_prev)
                    if prev_metric is not None and len(prev_metric) > 0:
                        prev_val = prev_metric[0]

                    if prev_val > 1e-10:
                        rel_diff = abs(curr_val - prev_val) / prev_val
                        if rel_diff >= iter_tol:
                            return False

        return True

    def finish(self):
        """Compute final weighted averages using hold time CDF weighting."""
        E = self.getNumberOfModels()
        if E == 0:
            return

        # Use last iteration results
        it = max([k[0] for k in self.results.keys()]) if self.results else 0
        if it == 0:
            return

        # Compute exit metrics weighted by hold time distribution
        QExit = {}
        UExit = {}
        TExit = {}

        # Determine dimensions
        M, K = 0, 0
        for e in range(E):
            result_e = self.results.get((it, e))
            if result_e is not None and result_e['Tran']['Avg']['Q'] is not None:
                Q_tran = result_e['Tran']['Avg']['Q']
                M = len(Q_tran)
                K = len(Q_tran[0]) if M > 0 else 0
                break

        if M == 0:
            return

        for e in range(E):
            result_e = self.results.get((it, e))
            # Seeded with the stage's own steady state, not with zeros; see
            # _stage_steady for why a metric with no trajectory is not a zero.
            QExit[e], UExit[e], TExit[e] = self._stage_steady(e, M, K)

            if result_e is None or result_e['Tran']['Avg']['Q'] is None:
                continue

            Q_tran = result_e['Tran']['Avg']['Q']
            U_tran = result_e['Tran']['Avg']['U']
            T_tran = result_e['Tran']['Avg']['T']

            for i in range(M):
                for r in range(K):
                    Qir = Q_tran[i][r]
                    t_vals, metric_vals = _get_tran_data(Qir)
                    if t_vals is not None:
                        # Use hold time MAP for weighting
                        hold_mmap = self.env_model.holdTime[e]
                        D0, D1 = hold_mmap[0], hold_mmap[1]

                        # Interpolate onto finer grid for accurate CDF weighting
                        t_fine, q_fine, u_fine, t_fine_data = \
                            _interpolate_for_cdf_map(
                                t_vals, metric_vals,
                                U_tran[i][r], T_tran[i][r], D0, D1)

                        cdf_vals = _map_eval_cdf(D0, D1, t_fine)
                        w = np.zeros(len(t_fine))
                        w[1:] = cdf_vals[1:] - cdf_vals[:-1]

                        if np.sum(w) > 0:
                            QExit[e][i, r] = np.dot(q_fine, w) / np.sum(w)
                            if u_fine is not None:
                                UExit[e][i, r] = np.dot(u_fine, w) / np.sum(w)
                            if t_fine_data is not None:
                                TExit[e][i, r] = np.dot(t_fine_data, w) / np.sum(w)
                        else:
                            # Fall back to final value
                            QExit[e][i, r] = metric_vals[-1] if len(metric_vals) > 0 else 0.0

        # Compute weighted averages across stages using steady-state probabilities
        Qval = np.zeros((M, K))
        Uval = np.zeros((M, K))
        Tval = np.zeros((M, K))

        for e in range(E):
            Qval += self.env_model.probEnv[e] * QExit[e]
            Uval += self.env_model.probEnv[e] * UExit[e]
            Tval += self.env_model.probEnv[e] * TExit[e]

        self.result = {
            'Avg': {
                'Q': Qval,
                'U': Uval,
                'T': Tval
            }
        }

        # 'meancov' also REPORTS the second moment, mixed over the stages by the
        # same law of total variance the handoff uses. The term
        # sum_e p_e (m_e - m)(m_e - m)' is what the environment itself
        # contributes: two stages with identical within-stage variance but
        # different means still leave the queue length varying, and averaging the
        # per-stage covariances alone would report none of it.
        if self._is_meancov():
            n = M * K
            Sval = np.zeros((n, n))
            mval = np.zeros(n)
            for e in range(E):
                p_e = float(self.env_model.probEnv[e])
                if not (p_e > 0):
                    continue
                me = np.asarray(QExit[e], dtype=float).ravel(order='F')
                Ce = self._stage_exit_cov(it, e, M, K)
                Sval += p_e * (Ce + np.outer(me, me))
                mval += p_e * me
            QCov = Sval - np.outer(mval, mval)
            QCov = 0.5 * (QCov + QCov.T)
            self.result['Avg']['QCov'] = QCov
            self.result['Avg']['QVar'] = np.maximum(
                np.diag(QCov), 0.0).reshape((M, K), order='F')

        # Cache-hit aggregation for fluid inner solvers.
        self._aggregate_cache_meanfield()

    def _substitute_cache_blend(self, E, ref, caches, K):
        """Cache blend for stages the mean-field cache fixed point cannot
        integrate. Each stage solver has already written its own
        hit/miss/delayed ratios onto its own cache node, so blend those by the
        per-stage arrival into the cache, the same rate weighting the
        'avg'/'dec' limits use. The occupancy handoff across environment
        switches is NOT modelled here: every stage is read at its own steady
        state, which is the substitution the warning names."""
        from .api.sn.getters import sn_get_node_arvr_from_tput

        warnings.warn(
            "The mean-field cache fixed point requires fluid stage solvers; "
            "blending each stage's own cache ratios instead, without the "
            "occupancy handoff across environment switches.")

        acc = {c: {'arv': np.zeros(K), 'hit': None, 'miss': None,
                   'dhit': None, 'hitl': None} for c in caches}
        for e in range(E):
            solver = self.getSolver(e)
            fn = getattr(solver, 'getAvgNodeTable', None) or getattr(solver, 'get_avg_node_table', None)
            if fn is not None:
                fn()
            out = solver.getAvg()
            TN = np.asarray(out[3], dtype=float)
            ANn = sn_get_node_arvr_from_tput(self.ensemble[e].getStruct(), TN)
            nodes = self.ensemble[e].get_nodes()
            pe = float(self.env_model.probEnv[e])
            for c in caches:
                arv = np.zeros(K)
                if ANn is not None and ANn.size > 0 and ANn.shape[0] > c:
                    row = np.asarray(ANn[c, :], dtype=float).ravel()
                    n = min(K, row.size)
                    arv[:n] = row[:n]
                arv[~np.isfinite(arv)] = 0.0
                arv = arv * pe
                acc[c]['arv'] += arv
                node = nodes[c]
                for key, getter in (('hit', 'get_hit_ratio'),
                                    ('miss', 'get_miss_ratio'),
                                    ('dhit', 'get_delayed_hit_ratio')):
                    gf = getattr(node, getter, None)
                    v = gf() if gf is not None else None
                    if v is None:
                        continue
                    v = np.asarray(v, dtype=float).ravel()
                    if v.size == 0:
                        continue
                    v = np.where(np.isfinite(v), v, 0.0)
                    if acc[c][key] is None:
                        acc[c][key] = np.zeros(v.size)
                    n = min(v.size, arv.size, acc[c][key].size)
                    acc[c][key][:n] += arv[:n] * v[:n]
                gf = getattr(node, 'get_hit_ratio_by_list', None)
                hl = gf() if gf is not None else None
                if hl is None:
                    continue
                hl = np.atleast_2d(np.asarray(hl, dtype=float))
                if hl.size == 0:
                    continue
                hl = np.where(np.isfinite(hl), hl, 0.0)
                if acc[c]['hitl'] is None:
                    acc[c]['hitl'] = np.zeros_like(hl)
                wrow = np.zeros(hl.shape[0])
                n = min(hl.shape[0], arv.size)
                wrow[:n] = arv[:n]
                acc[c]['hitl'] += wrow[:, None] * hl

        ref_nodes = ref.get_nodes()
        records = []
        for c in caches:
            d = np.where(acc[c]['arv'] > 0, acc[c]['arv'], np.nan)
            node = ref_nodes[c]
            rec = {'name': getattr(node, 'name', None), 'node': c, 'hitprob': None,
                   'missprob': None, 'delayedhitprob': None, 'hitproblist': None}
            for key, field, setter in (('hit', 'hitprob', 'set_result_hit_prob'),
                                       ('miss', 'missprob', 'set_result_miss_prob'),
                                       ('dhit', 'delayedhitprob', 'set_result_delayed_hit_prob')):
                a = acc[c][key]
                if a is None:
                    continue
                v = np.full(a.size, np.nan)
                n = min(a.size, d.size)
                v[:n] = a[:n] / d[:n]
                rec[field] = v
                sf = getattr(node, setter, None)
                if sf is not None:
                    sf(v)
            hl = acc[c]['hitl']
            if hl is not None:
                wrow = np.full(hl.shape[0], np.nan)
                n = min(hl.shape[0], d.size)
                wrow[:n] = d[:n]
                v = hl / wrow[:, None]
                rec['hitproblist'] = v
                sf = getattr(node, 'set_result_hit_prob_list', None)
                if sf is not None:
                    sf(v)
            records.append(rec)
        self._cache_records = records

    def get_avg_cache_results(self):
        """The environment-aggregated cache surface, one record per Cache node.

        EVERY FIELD IS OPTIONAL AND ABSENT (None) MEANS NOT COMPUTED, never
        zero. Returned as a list rather than read off ensemble[0]'s nodes,
        because under mapEnvApprox those nodes belong to a throwaway copy the
        caller never sees.
        """
        recs = getattr(self, '_cache_records', None)
        if recs is None:
            self.getAvg()
            recs = getattr(self, '_cache_records', None)
        return recs if recs is not None else []

    getAvgCacheResults = get_avg_cache_results

    def _aggregate_cache_meanfield(self):
        """Cache-hit aggregation for fluid (FLD) inner solvers, mirroring the
        MATLAB aggregateCacheMeanfield_. Runs a self-contained mean-field fixed
        point that carries each cache's mean occupancy across environment
        switches (the cache analog of the queue-length handoff): stage e is
        integrated over its sojourn from an entry occupancy that mixes the exit
        occupancy of its predecessors by probOrig. At convergence the per-class
        hit throughput is the probEnv-weighted, sojourn-averaged arrival x
        hit-prob, and the reported hit ratio is hit/(hit+miss), written onto the
        reference model's cache nodes.
        """
        from .api.cache import (cache_gamma_lp, cache_miss_rmf, cache_miss_fifo_rmf,
                                cache_miss_sfifo_rmf)
        from .api.sn.network_struct import NodeType
        from .lang.base import ReplacementStrategy

        E = self.getNumberOfModels()
        if E == 0:
            return

        ref = self.ensemble[0]
        # THE BLEND IS OVER THE CACHE NODES OF FLAT STAGE NETWORKS, and a
        # LAYERED stage has none: an LQN's struct is a LayeredNetworkStruct,
        # which carries no node space at all, so there is nothing here to scan
        # and no cache surface to report. Saying so by name is what the C++
        # twin does (`has_lqn_stages()` in solver_env_meanfield.h) and what
        # MATLAB's `~isfield(sn1,'nodetype')` return amounts to. Without it
        # this reads a struct accessor a LayeredNetwork does not have, and the
        # whole LQN-in-ENV run dies at the last step of finish() -- after the
        # queue, utilization and throughput tables it was asked for are already
        # computed and correct.
        if type(ref).__name__ == 'LayeredNetwork':
            return
        sn0 = ref.get_struct()
        caches = [ind for ind in range(sn0.nnodes) if sn0.nodetype[ind] == NodeType.CACHE]
        if not caches:
            return
        K = sn0.nclasses

        # Only fluid stages expose the RMF cache transient this fixed point
        # integrates. Any other stage solver still reports its own cache ratios,
        # so blend THOSE rather than reporting nothing: a silent return leaves
        # the caller with no hit ratio and no reason why.
        # see _kb/06-solver-catalog.md (ENV: Python environment.py additional notes) for rationale
        all_fluid = all(type(self.getSolver(e)).__name__ in ('SolverFLD', 'FLD', 'SolverFluid')
                        for e in range(E))
        if not all_fluid:
            self._substitute_cache_blend(E, ref, caches, K)
            return

        def _default_routing(h_val):
            mat = np.diag(np.ones(h_val), 1)
            mat[h_val, h_val] = 1.0
            return mat

        # Per-stage isolated cache inputs (open-cache arrival = source rate) and
        # the finite integration window per stage.
        stage_info = []
        tspan = []
        for e in range(E):
            sne = self.ensemble[e].get_struct()
            rates = np.asarray(sne.rates, dtype=float)
            refstat = np.asarray(sne.refstat).ravel().astype(int)
            infos = []
            for ind in caches:
                ch = sne.nodeparam.get(ind)
                m = np.asarray(ch.itemcap, dtype=int)
                n = int(ch.nitems)
                h = len(m)
                u = K
                arate = np.zeros(K)
                for r in range(K):
                    st = refstat[r] if r < len(refstat) else 0
                    if 0 <= st < rates.shape[0] and r < rates.shape[1] and np.isfinite(rates[st, r]):
                        arate[r] = rates[st, r]
                lam = np.zeros((u, n, h + 1))
                for v in range(u):
                    pv = ch.pread[v] if v < len(ch.pread) else None
                    if pv is not None and not (np.isscalar(pv) and np.isnan(pv)):
                        for ki in range(min(n, len(pv))):
                            for l in range(h + 1):
                                lam[v, ki, l] = arate[v] * pv[ki]
                accost = getattr(ch, 'accost', None)
                Rcost = accost
                if Rcost is None:
                    Rcost = [[_default_routing(h) for _ in range(n)] for _ in range(u)]
                gamma, _, _, _, _ = cache_gamma_lp(lam, Rcost)
                infos.append(dict(node=ind, gamma=gamma, m=m, lam=lam, arate=arate,
                                  accost=accost, strat=getattr(ch, 'replacestrat', None)))
            stage_info.append(infos)

            hd = self.env_model.holdTime[e]
            mean_e = _map_mean(hd[0], hd[1])
            ts = None
            solv = self.getSolver(e)
            opt = getattr(solv, 'options', None)
            if opt is not None:
                ts = opt.get('timespan', None) if isinstance(opt, dict) else getattr(opt, 'timespan', None)
            if ts is None or not np.all(np.isfinite(np.asarray(ts, dtype=float))):
                ts = [0.0, 20.0 * mean_e]
            tspan.append([float(ts[0]), float(ts[1])])

        # RANDOM(m), FIFO(m) and strict FIFO(m) each have a drift-based
        # transient (FIFO shares RANDOM's steady state, Gast15 Thm 1, but not its
        # transient, so it runs its own position-resolved drift), exactly as in
        # MATLAB solver_fld_cacheqn_tran. LRU/HLRU/CLIMB/QLRU have none, and a
        # fluid stage cannot have solved them either.
        tran_fn = {ReplacementStrategy.RR: cache_miss_rmf,
                   ReplacementStrategy.FIFO: cache_miss_fifo_rmf,
                   ReplacementStrategy.SFIFO: cache_miss_sfifo_rmf}
        for e in range(E):
            for cc, si in enumerate(stage_info[e]):
                if si['strat'] not in tran_fn:
                    raise RuntimeError(
                        "Transient cache analysis is only available for RANDOM(m)/FIFO(m) "
                        "and strict FIFO(m) replacement via a drift-based mean field; "
                        "cache %d uses a strategy without one." % (cc + 1))

        ncaches = len(caches)
        entry = [[None] * ncaches for _ in range(E)]
        SJ = [None] * E
        Exit = [None] * E
        max_sweep = max(1, int(self.options.get('iter_max', 100)))
        tol = float(self.options.get('iter_tol', 1e-4))
        prev = None
        from .api.cache.sojourn import cache_sojourn_clock
        for _sweep in range(max_sweep):
            for e in range(E):
                hd = self.env_model.holdTime[e]
                D0, D1 = hd[0], hd[1]
                pie_e = _map_pie(D0, D1)
                Exit_e = [None] * ncaches
                SJ_e = [None] * ncaches
                for cc, si in enumerate(stage_info[e]):
                    x0 = entry[e][cc]
                    # The sojourn average int x dF(t/Lam) / int dF(t/Lam) is
                    # integrated WITH the drift (RMF per-request time mapped to
                    # real time by the cache's request rate Lam), so it does not
                    # depend on the ODE output grid.
                    # see _kb/06-solver-catalog.md (ENV: Python environment.py additional notes) for rationale
                    # tspan is REAL time, the horizon the holding time is measured
                    # in, so the drift runs over Lam times it in its own
                    # per-request time; see MATLAB solver_fld_cacheqn_tran.
                    Lam = float(np.sum(si['arate']))
                    if not (Lam > 0):
                        continue
                    tdrift = [Lam * tspan[e][0], Lam * tspan[e][1]]
                    clock = cache_sojourn_clock(D0, pie_e, Lam, tdrift[0])
                    res = tran_fn[si['strat']](si['gamma'], si['m'], si['lam'], tspan=tdrift,
                                               x0init=x0, accost=si['accost'],
                                               sojourn_clock=clock)
                    sj = res[-1]
                    if sj['xbar'] is None or not (sj['wtot'] > 0):
                        continue
                    lam = si['lam']
                    hp = np.zeros(K)
                    mp = np.zeros(K)
                    for v in range(lam.shape[0]):
                        rr = float(np.sum(lam[v, :, 0]))
                        if rr > 0 and v < K:
                            mp[v] = min(1.0, max(0.0, sj['MU'][v] / rr))
                            hp[v] = 1.0 - mp[v]
                    sj['hitprob'] = hp
                    sj['missprob'] = mp
                    SJ_e[cc] = sj
                    Exit_e[cc] = sj['xbar']
                SJ[e], Exit[e] = SJ_e, Exit_e
            # Update each stage's entry occupancy from its predecessors.
            new_entry = [[None] * ncaches for _ in range(E)]
            for e in range(E):
                for cc in range(ncaches):
                    acc = None
                    for hh in range(E):
                        po = self.env_model.probOrig[hh, e]
                        if po > 0 and Exit[hh][cc] is not None:
                            term = po * Exit[hh][cc]
                            acc = term if acc is None else acc + term
                    new_entry[e][cc] = acc
            parts = [new_entry[e][cc].ravel() for e in range(E) for cc in range(ncaches)
                     if new_entry[e][cc] is not None]
            flat = np.concatenate(parts) if parts else np.array([])
            entry = new_entry
            if prev is not None and prev.shape == flat.shape and flat.size > 0 \
                    and np.max(np.abs(flat - prev)) < tol:
                prev = flat
                break
            prev = flat

        # Aggregate converged hit/miss throughputs and write onto the ref cache.
        ref_nodes = ref.get_nodes()
        cache_records = []
        for cc, ind in enumerate(caches):
            hitT = np.zeros(K)
            missT = np.zeros(K)
            hitLT = None
            for e in range(E):
                sj = SJ[e][cc] if (SJ[e] is not None and cc < len(SJ[e])) else None
                if sj is None:
                    continue
                pe = self.env_model.probEnv[e]
                # Sojourn-averaged per-item, per-list occupancy of this stage.
                xbar = sj['xbar']
                for k in range(K):
                    a = stage_info[e][cc]['arate'][k]
                    if a <= 0:
                        continue
                    # hit/miss are affine in the occupancy, so their sojourn
                    # averages are those of xbar.
                    hbar = sj['hitprob'][k]
                    mbar = sj['missprob'][k]
                    hitT[k] += pe * a * hbar
                    missT[k] += pe * a * mbar
                    hbl = _mf_hit_by_list(xbar, stage_info[e][cc]['lam'], k)
                    if hbl is None:
                        continue
                    if hitLT is None:
                        hitLT = np.zeros((K, hbl.size))
                    hitLT[k, :] += pe * a * hbl
            hitprob = np.full(K, np.nan)
            missprob = np.full(K, np.nan)
            hitproblist = None if hitLT is None else np.full(hitLT.shape, np.nan)
            for k in range(K):
                tot = hitT[k] + missT[k]
                if tot > 0:
                    hitprob[k] = hitT[k] / tot
                    missprob[k] = missT[k] / tot
                    if hitLT is not None:
                        hitproblist[k, :] = hitLT[k, :] / tot
            node = ref_nodes[ind]
            if hasattr(node, 'set_result_hit_prob'):
                node.set_result_hit_prob(hitprob)
            if hasattr(node, 'set_result_miss_prob'):
                node.set_result_miss_prob(missprob)
            if hitproblist is not None and not np.all(np.isnan(hitproblist)) \
                    and hasattr(node, 'set_result_hit_prob_list'):
                node.set_result_hit_prob_list(hitproblist)
            cache_records.append({'name': getattr(node, 'name', None), 'node': ind,
                                  'hitprob': hitprob, 'missprob': missprob,
                                  'delayedhitprob': None, 'hitproblist': hitproblist})
        self._cache_records = cache_records

    def iterate(self):
        """Run the main iteration loop.

        Mirrors JAR SolverENV.blending():
        - init() then pre(1) for initial steady-state
        - Loop: analyze all stages, post, converged
        - If max_iter hit without convergence: average last 10% of iterations
        - finish() for CDF-weighted aggregation
        """

        method = str(self.options.get('method', 'default')).lower()
        if method in ('avg', 'dec'):
            self._solve_env_limit(method)
            return

        self.init()

        it = 0
        iter_max = self.options.get('iter_max', 100)
        verbose = _env_verbose(self.options.get('verbose', False))

        # State-vector analyzer: carries the full per-stage distribution across
        # environment switches (mirrors MATLAB/JAR options.method='statevec').
        # 'blend' is a spelling of 'default' SINCE 2026-09-13 and named the
        # state-vector coupling before that, so it must NOT be matched here:
        # falling through to the mean-field analyzer is what it now means,
        # deliberately, where it used to be the defect. Only 'statevec' selects
        # the state-vector analyzer, in all four codebases.
        if str(self.options.get('method', 'default')).lower() == 'statevec':
            # see _kb/06-solver-catalog.md (ENV: Python environment.py additional notes) for rationale
            for e in range(len(self.ensemble)):
                if type(self.ensemble[e]).__name__ == 'LayeredNetwork':
                    raise RuntimeError(
                        "The state-vector (statevec) analyzer does not support "
                        "LayeredNetwork stages (stage %d): an LQN has no single "
                        "stage generator. Use the default mean-field analyzer "
                        "(omit method='statevec')." % e)
            from .solvers.solver_env.statevec import solver_env_statevec
            QN, UN, RN, TN, cache_records = solver_env_statevec(
                self.env_model, self._solvers, self.options)
            self._cache_records = cache_records
            self.result = {'Avg': {'Q': QN, 'U': UN, 'T': TN, 'R': RN,
                                   'Cache': cache_records}}
            return

        E = self.getNumberOfModels()

        while not self.converged(it) and it < iter_max:
            it += 1
            if verbose:
                print(f"ENV solver iteration {it}")

            self.pre(it)

            # Analyze each stage
            for e in range(E):
                result_e, runtime = self.analyze(it, e)
                self.results[(it, e)] = result_e

            self.post(it)

            # If max iterations reached without convergence, average last 10%
            # (matching JAR SolverENV.blending() lines 1822-1832)
            if it == iter_max:
                it_last = int(round(iter_max * 0.9))
                for e in range(E):
                    result_curr = self.results.get((it, e))
                    if result_curr is None or result_curr['Tran']['Avg']['Q'] is None:
                        continue
                    Q_tran = result_curr['Tran']['Avg']['Q']
                    M = len(Q_tran)
                    K = len(Q_tran[0]) if M > 0 else 0
                    # Average Q values at t=0 across last 10% of iterations
                    QN_avg = np.zeros((M, K))
                    count = 0
                    for it_tmp in range(it_last, it + 1):
                        res_tmp = self.results.get((it_tmp, e))
                        if res_tmp is not None and res_tmp['Tran']['Avg']['Q'] is not None:
                            Q_tmp = res_tmp['Tran']['Avg']['Q']
                            for i in range(M):
                                for k in range(K):
                                    _, m_vals = _get_tran_data(Q_tmp[i][k])
                                    if m_vals is not None and len(m_vals) > 0:
                                        QN_avg[i, k] += m_vals[0]
                            count += 1
                    if count > 0:
                        QN_avg /= count

        self.finish()

        if verbose:
            print(f"ENV solver converged after {it} iterations")

    def runAnalyzer(self):
        """Run the environment solver (MATLAB compatibility)."""
        self.iterate()

    def generator(self) -> Tuple[np.ndarray, List[np.ndarray]]:
        """
        Get the infinitesimal generator matrices for the random environment model.

        Returns the combined infinitesimal generator for the random environment
        and the individual stage generators.

        This method requires all sub-solvers to be CTMC solvers (SolverCTMC).

        Returns:
            Tuple of (renvInfGen, stageInfGen) where:
            - renvInfGen: Combined infinitesimal generator for the random environment (scipy sparse or dense)
            - stageInfGen: List of infinitesimal generators for each stage
        """
        from scipy import sparse
        from scipy.sparse import csr_matrix, lil_matrix

        E = self.getNumberOfModels()

        # Get stage generators from CTMC solvers
        stageInfGen = []
        for e in range(E):
            solver = self.getSolver(e)
            if solver is None:
                raise ValueError(f"No solver for stage {e}")

            # Check if solver is CTMC-based
            if hasattr(solver, 'getInfGen'):
                gen = solver.getInfGen()
                stageInfGen.append(gen)
            elif hasattr(solver, 'generator'):
                gen, _ = solver.generator()
                stageInfGen.append(gen)
            else:
                raise ValueError(f"Solver for stage {e} must be a CTMC solver with getInfGen() method")

        # Get number of states for each stage
        nstates = [g.shape[0] for g in stageInfGen]

        # Get number of phases for each transition distribution
        nphases_raw = np.zeros((E, E), dtype=int)
        for e in range(E):
            for h in range(E):
                dist = self.env_model.env[e][h]
                if dist is not None and hasattr(dist, 'getNumberOfPhases'):
                    nphases_raw[e, h] = dist.getNumberOfPhases()
                else:
                    nphases_raw[e, h] = 1

        # Calculate expanded dimensions for each stage
        # Each stage's state space is expanded by the phases of its outgoing transitions
        expanded_dims = []
        for e in range(E):
            dim = nstates[e]
            for h in range(E):
                if h != e and self.env_model.env[e][h] is not None:
                    dim *= nphases_raw[e, h]
            expanded_dims.append(dim)

        total_dim = sum(expanded_dims)

        # Initialize the combined generator as a dense matrix
        renvInfGen_flat = np.zeros((total_dim, total_dim))

        # Calculate block offsets
        offsets = [0]
        for d in expanded_dims[:-1]:
            offsets.append(offsets[-1] + d)

        # Build each block
        for e in range(E):
            # Compute the expanded diagonal block for stage e
            diag_block = stageInfGen[e].copy() if hasattr(stageInfGen[e], 'copy') else np.array(stageInfGen[e])

            # Get D0 matrices for all outgoing transitions from stage e
            D0_list = []
            for h in range(E):
                if h != e:
                    dist_eh = self.env_model.env[e][h]
                    if dist_eh is not None and hasattr(dist_eh, 'getD0'):
                        D0_list.append(dist_eh.getD0())

            # Kronecker sum the diagonal with all D0 matrices
            for D0 in D0_list:
                n_A = diag_block.shape[0]
                n_B = D0.shape[0]
                I_A = np.eye(n_A)
                I_B = np.eye(n_B)
                diag_block = np.kron(diag_block, I_B) + np.kron(I_A, D0)

            # Place diagonal block
            row_start = offsets[e]
            row_end = row_start + expanded_dims[e]
            renvInfGen_flat[row_start:row_end, row_start:row_end] = diag_block

            # Build off-diagonal blocks (transitions from stage e to stage f)
            for f in range(E):
                if f != e:
                    dist_ef = self.env_model.env[e][f]
                    if dist_ef is not None:
                        # Get D1 matrix for transition completion
                        D1 = dist_ef.getD1() if hasattr(dist_ef, 'getD1') else np.array([[1.0]])
                        nph_ef = nphases_raw[e, f]

                        # Get initial phase probability for the target stage
                        dist_fe = self.env_model.env[f][e]
                        if dist_fe is not None and hasattr(dist_fe, 'getInitProb'):
                            pie_fe = dist_fe.getInitProb()
                        else:
                            # Count phases of outgoing transition from f
                            nph_fe = 1
                            for h in range(E):
                                if h != f and self.env_model.env[f][h] is not None:
                                    nph_fe = nphases_raw[f, h]
                                    break
                            pie_fe = np.zeros(nph_fe)
                            pie_fe[0] = 1.0

                        # Build the transition rate matrix
                        # onePhase vector sums over completion phases
                        onePhase = np.ones((nph_ef, 1))
                        # D1 * onePhase gives the rate of completing the transition
                        rate_vec = D1 @ onePhase  # nph_ef x 1

                        # Reset matrix: maps source state to target state
                        minStates = min(nstates[e], nstates[f])
                        reset_base = np.zeros((nstates[e], nstates[f]))
                        for i in range(minStates):
                            reset_base[i, i] = 1.0

                        # see _kb/06-solver-catalog.md (ENV: Python environment.py additional notes) for rationale
                        block = reset_base.copy()

                        # Expand source dimension by phases of OTHER outgoing transitions
                        for h in range(E):
                            if h != e and h != f:
                                dist_eh = self.env_model.env[e][h]
                                if dist_eh is not None:
                                    nph_eh = nphases_raw[e, h]
                                    # Need to select the right phase slice
                                    block = np.kron(block, np.ones((nph_eh, 1)))

                        # Expand by the completing transition phases (D1 * onePhase)
                        block = np.kron(block, rate_vec)

                        # Expand target dimension by phases of outgoing transitions from target
                        block = np.kron(block, pie_fe.reshape(1, -1))

                        # Place the block (may need reshaping)
                        col_start = offsets[f]
                        col_end = col_start + expanded_dims[f]

                        # The block dimensions should match
                        if block.shape == (expanded_dims[e], expanded_dims[f]):
                            renvInfGen_flat[row_start:row_end, col_start:col_end] = block
                        else:
                            # Reshape if necessary - this handles dimension mismatches
                            try:
                                block_reshaped = block.reshape(expanded_dims[e], expanded_dims[f])
                                renvInfGen_flat[row_start:row_end, col_start:col_end] = block_reshaped
                            except ValueError:
                                # If reshape fails, try to fit what we can
                                min_rows = min(block.shape[0], expanded_dims[e])
                                min_cols = min(block.shape[1], expanded_dims[f])
                                renvInfGen_flat[row_start:row_start+min_rows,
                                              col_start:col_start+min_cols] = block[:min_rows, :min_cols]

        # Normalize to make it a valid infinitesimal generator (rows sum to 0)
        renvInfGen_flat = self._ctmc_makeinfgen(renvInfGen_flat)

        return renvInfGen_flat, stageInfGen

    def get_generator(self) -> Tuple[np.ndarray, List[np.ndarray]]:
        """Alias of :meth:`generator`, under the name MATLAB's
        `@SolverENV/getGenerator` uses. Returns the same
        `(renvInfGen, stageInfGen)` pair; MATLAB's single-output call takes the
        first element."""
        return self.generator()

    getGenerator = get_generator

    def _ctmc_makeinfgen(self, Q: np.ndarray) -> np.ndarray:
        """
        Normalize a matrix to be a valid infinitesimal generator.

        Ensures:
        - Off-diagonal elements are non-negative
        - Diagonal elements make rows sum to 0

        Args:
            Q: Input matrix

        Returns:
            Normalized infinitesimal generator
        """
        Q = Q.copy()
        n = Q.shape[0]

        # Make off-diagonal elements non-negative
        for i in range(n):
            for j in range(n):
                if i != j and Q[i, j] < 0:
                    Q[i, j] = 0.0

        # Set diagonal to make rows sum to 0
        for i in range(n):
            Q[i, i] = 0.0
            row_sum = np.sum(Q[i, :])
            Q[i, i] = -row_sum

        return Q

    def stage_timespan(self):
        """The finite horizon each stage transient is integrated over, or None.

        THE HORIZON LIVES ON THE STAGE SOLVER, not on the ensemble: an example
        writes `ENV(env, lambda m: FLD(m, timespan=[0, 1000]))`, so a bridge that
        forwards only the ENV solver's own options states nothing and the engine
        falls back to its default 100. The stages are then integrated over a
        different interval and the exit averages are a different quadrature --
        worth 6e-3 relative on renv_twostages_repairmen, whose stages ask for
        [0, 1e3]. The ENV options are the fallback, for a caller that states the
        horizon there instead.
        """
        for e, s in enumerate(self._solvers):
            if s is None:
                continue
            opts = getattr(s, 'options', None)
            ts = None
            if isinstance(opts, dict):
                ts = opts.get('timespan')
            elif opts is not None:
                ts = getattr(opts, 'timespan', None)
            span = _finite_span(ts)
            if span is not None:
                return span
            # A NON-finite horizon is not "unstated": the native path resolves
            # it to 30/minrate and integrates that, so a bridge that dropped it
            # here would have the remote engine answer over its own default of
            # 100 instead -- a different question, and most of the parity slack
            # renv_fourstages_repairmen used to carry.
            span = _finite_span(self._stage_horizon(e))
            if span is not None:
                return span
        return _finite_span(self.options.get('timespan'))

    def _stage_horizon(self, e):
        """Stage e's transient horizon, resolved; None when it already has one.

        A stage solver left at `timespan=[0, inf]` -- what the renv examples
        write when they have no particular horizon in mind -- is otherwise
        resolved by each backend's own convention, and the conventions
        disagree: MATLAB's getTranAvg takes `30/minrate`, the native fluid
        handler takes `min(timespan[1], 10*iter_max/min_rate)`, and the C++/JAR
        `-s env` arm falls back to 100. Each codebase then integrates a
        DIFFERENT interval and the exit averages are different quadratures --
        that, and not the integrator, was most of the parity slack on
        renv_fourstages_repairmen. MATLAB is ground truth, so its rule is the
        one used here too.

        PURE ON PURPOSE, as MATLAB's stageHorizon_ is: writing it onto the
        stage solver would outlive the ENV solve and change what a later
        per-stage getter answers.
        """
        try:
            s = self._solvers[e]
        except (IndexError, TypeError):
            return None
        if s is None:
            return None
        opts = getattr(s, 'options', None)
        if opts is None:
            return None
        ts = opts.get('timespan') if isinstance(opts, dict) else getattr(opts, 'timespan', None)
        if ts is None or len(ts) < 2 or np.isfinite(ts[1]):
            return None
        sn = self._stage_struct(e)
        if sn is None:
            return None
        rates = np.asarray(sn.rates, dtype=float)
        finite = rates[np.isfinite(rates)]
        if finite.size == 0:
            return None
        minrate = float(np.min(finite))
        if not (minrate > 0):
            return None
        t0 = 0.0 if not np.isfinite(ts[0]) else float(ts[0])
        resolved = [t0, 30.0 / minrate]
        if e not in self._horizon_warned:
            self._horizon_warned.add(e)
            from .api.io.logging import line_warning
            line_warning('SolverENV',
                         "End time of transient analysis unspecified for stage %d, setting its "
                         "timespan option to [%g,%g]. Pass a stage solver with timespan=[0,T] "
                         "to customize." % (e + 1, resolved[0], resolved[1]))
        return resolved

    def _set_stage_timespan(self, e, ts):
        """Put stage e's horizon back."""
        opts = getattr(self._solvers[e], 'options', None)
        if opts is None:
            return
        if isinstance(opts, dict):
            opts['timespan'] = ts
        else:
            opts.timespan = ts

    def _cdf_grid_for_stage(self, e):
        """The instants stage e's exit average is summed over, or None.

        Mirrors cdfGrid_ in solver_env_meanfield_analyzer.m: the same 5000
        points, 90% of them below 5*E[S], built from the stage HORIZON rather
        than from a trajectory, so it can be requested before the solve. None
        when the horizon is not finite, i.e. when there is no grid to ask for.
        """
        solver = self.getSolver(e)
        if solver is None:
            return None
        opts = getattr(solver, 'options', None)
        if opts is None:
            return None
        ts = opts.get('timespan') if isinstance(opts, dict) else getattr(opts, 'timespan', None)
        if ts is None or len(ts) < 2 or not np.isfinite(ts[1]) or not (ts[1] > ts[0]):
            return None
        try:
            hold = self.env_model.holdTime[e]
        except (AttributeError, IndexError, TypeError):
            return None
        if hold is None or len(hold) < 2:
            return None
        mean_sojourn = _map_mean(np.asarray(hold[0]), np.asarray(hold[1]))
        if not (mean_sojourn > 0) or not np.isfinite(mean_sojourn):
            mean_sojourn = (float(ts[1]) - float(ts[0])) / 10.0
        return _cdf_refine_grid(np.array([float(ts[0]), float(ts[1])]), mean_sojourn, 5000)

    def _stage_struct(self, e):
        """The NetworkStruct of stage e, or None when it cannot be built."""
        try:
            model = self.ensemble[e]
        except (IndexError, TypeError):
            return None
        for name in ('getStruct', 'get_struct'):
            fn = getattr(model, name, None)
            if fn is not None:
                try:
                    return fn()
                except Exception:
                    return None
        return None

    def _assert_stages_delegable(self, lang):
        """Refuse a delegated solve whose stage solver the engine cannot run.

        THE CHOICE IS NOT IN THE MODEL, which is the root of this: model.json
        carries the stage NETWORKS and the transition process, and the ensemble's
        solver choice lives on this object. An engine that defaults to the fluid
        transient therefore answered a DIFFERENT model in silence -- on
        renv_threestages_repairmen, whose stages are SolverCTMC, `-s env`
        returned Queue1 throughput 1.5577 against the 1.3333 the CTMC stages
        give, a 17% substitution reported as that example's result.

        Both engines now take `--stage-solver`, and both run the enumerated CTMC
        as well as the fluid transient -- the JAR's `SolverENV` always handled a
        non-fluid stage (`roundMarginalForDiscreteSolver` runs for exactly that
        case) and only its CLI hardwired the factory. `_env_stage_solver` sends
        the token on either route. A MIXED ensemble is still refused, since one
        coupling runs one stage solver and sending the first stage's name would
        answer for the rest under it.
        """
        names = sorted(set(type(s).__name__ for s in self._solvers if s is not None))
        if all('FLD' in n or 'Fluid' in n for n in names):
            return
        if names and all('CTMC' in n for n in names):
            return
        raise ValueError(
            "SolverENV does not support lang='%s' with %s stages: the %s ENV engine "
            "solves every stage by the fluid transient. Use lang='python' or "
            "lang='matlab' for an ensemble built on another stage solver."
            % (lang, ', '.join(names) or 'these', 'C++' if lang == 'cpp' else 'JAR'))

    def _print_banner(self, runtime, iterations):
        """The completion line for an environment solve.

        ANNOUNCE THE SOLVER, as every network solver does: an ensemble result
        table that arrives with no banner belongs to nobody, and a reader that
        has to guess will guess -- parity-static tags it UNKNOWN and its
        agreement-gated generator then drops it as a parsing artifact. The text
        is the one MATLAB's and the JAR's SolverENV print.
        """
        opts = self.options
        if not _env_verbose(opts.get('verbose') if isinstance(opts, dict)
                            else getattr(opts, 'verbose', False)):
            return
        method = str(opts.get('method', 'default') or 'default')
        lang = str(opts.get('lang', 'python') or 'python')
        from line_solver.solvers.base import print_solver_banner
        print_solver_banner("ENV analysis [method: %s; type: approximate, deterministic; lang: %s; "
              "env: python] completed in %fs. Iterations: %d."
              % (method, lang, runtime, iterations))

    def avg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute average performance metrics across environments.

        Mirrors MATLAB SolverENV.getEnsembleAvg() which calls iterate()
        and returns environment-weighted Q, U, T from self.result.Avg.

        Returns:
            Tuple of (QN, UN, TN) - queue lengths, utilizations, throughputs
        """
        import time as _time
        # Solver console: avg() is the single point every ENV accessor reaches
        # (getAvg, getEnsembleAvg and avg_table all land here), so the narrated
        # run is opened around the WHOLE analysis. Opening it in iterate() left
        # the stage solves outside it, and each of them narrated as a run of
        # its own: 5274 lines for a two-stage model.
        from line_solver.api.io import console as _console
        _console.begin_run(self, self.options)
        _t0 = _time.time()
        try:
            return self._avg()
        finally:
            _console.close_run(self)
            # EVERY exit announces itself once: the native fixed point, the two
            # bridges and the empty-result return leave by different paths, and
            # a banner printed on only some of them labels some tables and
            # leaves others anonymous.
            iters = 0
            for key in (self.results or {}):
                if isinstance(key, tuple) and key:
                    iters = max(iters, int(key[0]))
            self._print_banner(_time.time() - _t0, iters)

    def _avg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        # lang='cpp' delegates the whole environment solve to line-cli's -s env
        # arm. It is checked before iterate() because answering natively would
        # report a python number under lang='cpp'.
        if str(self.options.get('lang', 'python')) == 'cpp':
            self._assert_stages_delegable('cpp')
            from .solvers.cpp_dispatch import env_avg_via_cpp
            Q, U, T = env_avg_via_cpp(self)
            self.result = {'Avg': {'Q': Q, 'U': U, 'T': T}}
            self._result = (Q, U, T)
            return Q, U, T

        # lang='java' delegates the WHOLE environment solve for the same reason,
        # and it has to: running this coupling with JAR stage solves sends each
        # stage's entry marginal over model.json, which carries no such state,
        # so every stage restarted from the default and the loop settled on the
        # uncoupled answer (renv_node_breakdown Server QLen 0.39634 at
        # throughput 0.71574, against a source admitting 0.8).
        if str(self.options.get('lang', 'python')) == 'java':
            self._assert_stages_delegable('java')
            from .solvers.jar_dispatch import env_avg_via_jar
            Q, U, T = env_avg_via_jar(self)
            self.result = {'Avg': {'Q': Q, 'U': U, 'T': T}}
            self._result = (Q, U, T)
            return Q, U, T

        if self.result is None:
            self.iterate()

        if self.result is None:
            return np.array([]), np.array([]), np.array([])

        Q = self.result['Avg']['Q']
        U = self.result['Avg']['U']
        T = self.result['Avg']['T']

        self._result = (Q, U, T)
        return Q, U, T

    def getAvgQLenCov(self):
        """Covariance of the station-class queue lengths across the environment.

        An (M*K, M*K) matrix indexed ir = r*M + i (column-major over (M, K), the
        same flattening MATLAB's (r-1)*M+i gives). Produced by the 'meancov'
        coupling only. It is the law of total variance over the environment
        stages,

            sum_e p_e (C_e + m_e m_e') - m m',

        so it keeps BOTH what varies inside a stage and what varies because the
        stages differ from one another. The other couplings carry a first moment
        only and this raises rather than returning a zero that reads like a
        computed answer.
        """
        self.runAnalyzer()
        result = getattr(self, 'result', None)
        if not isinstance(result, dict) or result.get('Avg', {}).get('QCov') is None:
            raise ValueError(
                "getAvgQLenCov needs options['method']='meancov'; the other SolverENV "
                "couplings carry the mean queue lengths alone and compute no second "
                "moment.")
        return result['Avg']['QCov']

    def getAvgQLenVar(self):
        """Variance of the queue length at each station and class across the
        environment, an (M, K) array.

        The diagonal of GETAVGQLENCOV reshaped, with the same requirement of
        options['method']='meancov'.
        """
        self.getAvgQLenCov()
        return self.result['Avg']['QVar']

    def getAvg(self):
        """Get average metrics (MATLAB compatibility).

        Returns (QN, UN, RN, TN, AN, WN): MATLAB's SolverENV.getAvg forwards to
        getEnsembleAvg, so throughput is the FOURTH output and RN is NaN
        throughout (ENV computes no response time). Use avg() for the bare
        (Q, U, T) triple.
        """
        return self.getEnsembleAvg()

    def getEnsembleAvg(self):
        """Get ensemble average metrics (MATLAB compatibility).

        Returns (QN, UN, RN, TN, AN, WN) matching MATLAB signature.
        """
        Q, U, T = self.avg()
        if len(Q) == 0:
            return Q, U, np.array([]), T, np.array([]), np.array([])
        W = np.where(T > 1e-10, Q / T, 0.0)
        R = np.full_like(W, np.nan)
        A = np.full_like(T, np.nan)
        return Q, U, R, T, A, W

    def avg_table(self) -> pd.DataFrame:
        """Get average metrics as a table.

        Mirrors MATLAB SolverENV.getAvgTable: iterates over stations and classes,
        uses sn{1} (first stage struct) for station/class names, filters zero rows.
        """
        if self._result is None:
            self.avg()

        QN, UN, TN = self._result

        # Ensure 2D
        QN = np.atleast_2d(QN)
        UN = np.atleast_2d(UN)
        TN = np.atleast_2d(TN)

        M = QN.shape[0]  # number of stations
        K = QN.shape[1]  # number of classes

        # Get station and class names from first ensemble model
        station_names = []
        class_names = []
        model = self.ensemble[0] if len(self.ensemble) > 0 else None
        if model is not None and type(model).__name__ == 'LayeredNetwork':
            # LQN stage: labels come from the block-diagonal layer networks
            # (prefixed by the layer name to disambiguate repeated names).
            Roff, Coff, Msz, Ksz = model._layer_blocks()
            station_names = [f'Station{i}' for i in range(M)]
            class_names = [f'Class{k}' for k in range(K)]
            for le, layer in enumerate(model._layer_ensemble()):
                lname = layer.get_name() if hasattr(layer, 'get_name') else getattr(layer, 'name', f'Layer{le}')
                lstations = layer.get_stations() if hasattr(layer, 'get_stations') else []
                lclasses = layer.get_classes() if hasattr(layer, 'get_classes') else []
                for i in range(Msz[le]):
                    nm = lstations[i].name if i < len(lstations) and hasattr(lstations[i], 'name') else f'S{i}'
                    if Roff[le] + i < M:
                        station_names[Roff[le] + i] = f'{lname}.{nm}'
                for r in range(Ksz[le]):
                    nm = lclasses[r].name if r < len(lclasses) and hasattr(lclasses[r], 'name') else f'C{r}'
                    if Coff[le] + r < K:
                        class_names[Coff[le] + r] = f'{lname}.{nm}'
        elif model is not None:
            # Get stations (not all nodes)
            stations = model.get_stations() if hasattr(model, 'get_stations') else []
            if not stations and hasattr(model, '_stations'):
                stations = model._stations
            nodes = model.get_nodes() if hasattr(model, 'get_nodes') else []
            if not nodes and hasattr(model, '_nodes'):
                nodes = model._nodes

            for ist in range(M):
                if ist < len(stations):
                    st = stations[ist]
                    # Get node name for this station (MATLAB: nodenames{stationToNode(ist)})
                    if st in nodes:
                        station_names.append(st.name if hasattr(st, 'name') else f'Station{ist}')
                    else:
                        station_names.append(st.name if hasattr(st, 'name') else f'Station{ist}')
                else:
                    station_names.append(f'Station{ist}')

            classes = model.get_classes() if hasattr(model, 'get_classes') else []
            if not classes and hasattr(model, '_classes'):
                classes = model._classes
            for k in range(K):
                if k < len(classes):
                    cls = classes[k]
                    class_names.append(cls.name if hasattr(cls, 'name') else f'Class{k}')
                else:
                    class_names.append(f'Class{k}')
        else:
            station_names = [f'Station{i}' for i in range(M)]
            class_names = [f'Class{k}' for k in range(K)]

        # Build table rows (MATLAB: filter rows where QN+UN+TN > 0)
        data = []
        for ist in range(M):
            for k in range(K):
                qlen = float(QN[ist, k])
                util = float(UN[ist, k])
                tput = float(TN[ist, k])
                if qlen + util + tput > 0:
                    respt = qlen / tput if tput > 1e-10 else 0.0
                    row = {
                        'Station': station_names[ist] if ist < len(station_names) else f'Station{ist}',
                        'JobClass': class_names[k] if k < len(class_names) else f'Class{k}',
                        'QLen': qlen,
                        'Util': util,
                        'RespT': respt,
                        'Tput': tput,
                    }
                    data.append(row)

        return pd.DataFrame(data)

    def getAvgTable(self):
        """Get average table (MATLAB compatibility)."""
        return self.avg_table()

    def avgT(self):
        """Short alias for getAvgTable (matches EnsembleSolver in MATLAB/JAR;
        native SolverENV does not extend EnsembleSolver)."""
        return self.getAvgTable()

    aT = avgT
    getAvgT = avgT

    @staticmethod
    def default_options():
        """Get default solver options."""
        return {
            'method': 'default',
            'iter_max': 100,
            'iter_tol': 1e-4,
            'verbose': VerboseLevel.STD  # STD like every other solver; pass
            # VerboseLevel.SILENT for a quiet run (see api/io/console.py)
        }

    @staticmethod
    def defaultOptions():
        """Get default solver options (MATLAB compatibility)."""
        return SolverENV.default_options()

    @staticmethod
    def getFeatureSet() -> set:
        """Return the set of features supported by SolverENV.

        Matches JAR SolverENV.getFeatureSet().
        """
        return {
            # Nodes
            'ClassSwitch', 'Delay', 'DelayStation', 'Queue',
            'Sink', 'JobSink', 'Source',
            # Distributions
            'Coxian', 'Cox2', 'Erlang', 'Exp', 'HyperExp',
            # Sections
            'StatelessClassSwitcher', 'InfiniteServer', 'SharedServer',
            'Buffer', 'Dispatcher', 'Server', 'RandomSource', 'ServiceTunnel',
            # Scheduling strategies
            'SchedStrategy_INF', 'SchedStrategy_PS', 'SchedStrategy_FCFS',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND', 'RoutingStrategy_RROBIN',
            # Customer Classes
            'ClosedClass', 'OpenClass',
        }

    def supports(self, model) -> bool:
        """Check whether the given model is supported by SolverENV.

        Args:
            model: Network model to check

        Returns:
            True if the model uses only supported features
        """
        # see _kb/06-solver-catalog.md (ENV: Python environment.py additional notes) for rationale
        from .solvers.base import supports_via_featureset
        return supports_via_featureset(type(self), model)

    @staticmethod
    def listValidMethods() -> List[str]:
        """Return list of valid solution methods.

        SolverENV.m verbatim. Each name selects a COUPLING -- what crosses an
        environment switch -- and each is dispatched here or in the analyzers:
        'default'/'meanfield'/'mean'/'blend'/'blending' carry the marginal means, 'meancov'
        carries a covariance beside them (the exit second moment over the sojourn
        law, mixed over the origins of a switch by the law of total variance and
        handed on as config['init_qlen']/['init_qcov']), 'statevec' the whole
        joint distribution, 'smp' recomputes probEnv/probOrig from the embedded
        jump chain of the arc CDFs and the mean holding times, 'statedep' makes the transition depend on the state it
        leaves, and 'avg'/'dec' are the closed-form fast/slow environment limits.
        This used to name only 'default' and 'smp' while four more were
        dispatched.

        'blend' NAMES THE MEAN-FIELD COUPLING SINCE 2026-09-13. It is kept as a
        spelling of 'default' so an existing script still runs, but it no longer
        reaches the state-vector analyzer and the numbers it returns move: ask
        for 'statevec' by name to carry the joint distribution.
        """
        return ['default', 'meanfield', 'mean', 'meancov', 'blend', 'blending', 'smp', 'statedep',
                'statevec', 'avg', 'dec']

    def getStruct(self) -> List:
        """Return network structures for all stages.

        Matches JAR SolverENV.getStruct().
        """
        E = self.getNumberOfModels()
        envsn = []
        for e in range(E):
            model = self.ensemble[e] if e < len(self.ensemble) else None
            if model is not None and hasattr(model, 'get_struct'):
                envsn.append(model.get_struct())
            else:
                envsn.append(None)
        return envsn

    def setRef(self, i: int):
        """Set reference station index for system throughput.

        Args:
            i: Station index to use as reference
        """
        self._ref = i

    def getName(self) -> str:
        """Return solver name."""
        return 'SolverENV'

    def runAnalyzerByCTMC(self) -> Dict:
        """Run analysis using direct CTMC composition.

        Builds a combined infinitesimal generator for the random environment
        by Kronecker-structuring the per-stage generators with environment
        transition rates, then solves the combined CTMC for steady-state.

        Matches JAR SolverENV.runAnalyzerByCTMC().

        Returns:
            Dict with keys 'QN', 'UN', 'TN', 'infGen'
        """
        self.init()

        E = self.getNumberOfModels()
        model0 = self.ensemble[0]
        sn0 = model0.get_struct() if hasattr(model0, 'get_struct') else None
        if sn0 is None:
            raise RuntimeError("Cannot get struct from first stage model")

        M = sn0.nstations
        K = sn0.nclasses

        # Build environment rate matrix E0
        E0 = np.zeros((E, E))
        for i in range(E):
            for j in range(E):
                dist = self.env_model.env[i][j]
                if dist is not None:
                    E0[i, j] = _get_rate(dist)
        # Make infgen
        for i in range(E):
            E0[i, i] = 0
            E0[i, i] = -np.sum(E0[i, :])

        # Get stage generators from CTMC solvers
        stage_infgen = []
        state_spaces = []
        for e in range(E):
            solver = self.getSolver(e)
            if solver is None:
                raise ValueError(f"No solver for stage {e}")
            if hasattr(solver, 'getGenerator'):
                gen_result = solver.getGenerator()
                if isinstance(gen_result, tuple):
                    stage_infgen.append(gen_result[0])
                    if len(gen_result) > 1:
                        state_spaces.append(gen_result[1])
                    else:
                        state_spaces.append(None)
                else:
                    stage_infgen.append(gen_result)
                    state_spaces.append(None)
            else:
                raise ValueError(f"Solver for stage {e} must support getGenerator()")

        nstates = [g.shape[0] for g in stage_infgen]
        states = nstates[0]  # Assume same state space size

        # Build combined generator Q
        total_states = E * states
        Q = np.zeros((total_states, total_states))
        for e in range(E):
            for h in range(E):
                if e == h:
                    block = stage_infgen[e].copy()
                    block += np.eye(states) * E0[e, e]
                    Q[e*states:(e+1)*states, h*states:(h+1)*states] = block
                else:
                    Q[e*states:(e+1)*states, h*states:(h+1)*states] = np.eye(states) * E0[e, h]

        # Make infgen
        for i in range(total_states):
            Q[i, i] = 0
            Q[i, i] = -np.sum(Q[i, :])

        # Solve for steady-state
        pi = self.env_model._solve_ctmc(Q)

        # Extract metrics
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        TN = np.zeros((M, K))

        for e in range(E):
            for s in range(states):
                p = pi[e * states + s]
                if state_spaces[e] is not None:
                    ss = state_spaces[e]
                    for m in range(M):
                        n_total = 0
                        for r in range(K):
                            n_total += ss[s, m * K + r] if ss.shape[1] > m * K + r else 0
                        sn_e = self.ensemble[e].get_struct() if hasattr(self.ensemble[e], 'get_struct') else None
                        nservers = sn_e.nservers[m] if sn_e is not None and hasattr(sn_e, 'nservers') else float('inf')
                        c = nservers

                        for k in range(K):
                            col_idx = m * K + k
                            if col_idx < ss.shape[1]:
                                prob = ss[s, col_idx]
                            else:
                                prob = 0
                            QN[m, k] += p * prob
                            if np.isinf(c):
                                UN[m, k] += p * prob
                                rate_mk = sn_e.rates[m, k] if sn_e is not None else 0
                                TN[m, k] += p * prob * rate_mk
                            else:
                                scaling = min(n_total, c) / n_total if n_total > 0 else 0
                                rate_mk = sn_e.rates[m, k] if sn_e is not None else 0
                                TN[m, k] += p * prob * rate_mk * scaling
                                u_contrib = prob * min(n_total, c) / (n_total * c) if n_total > 0 else 0
                                UN[m, k] += p * u_contrib

        return {'QN': QN, 'UN': UN, 'TN': TN, 'infGen': Q}


# Convenience aliases
ENV = SolverENV
SolverEnv = SolverENV


__all__ = [
    'Environment',
    'SolverENV',
    'SolverEnv',
    'ENV',
]
