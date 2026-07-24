"""
Native Python implementation of Random Environment models.

This module provides classes for defining and analyzing queueing networks
in random environments, where the network parameters change according to
an underlying Markov modulated process.

Implements full parity with MATLAB SolverENV using transient analysis
with iteration until convergence.
"""

import numpy as np
from typing import Optional, List, Dict, Any, Union, Callable, Tuple
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


def _interpolate_for_cdf(t_vals, q_metric, u_tran_ir, t_tran_ir, dist,
                         n_interp=500, min_points=50):
    """Interpolate transient data onto a finer time grid for CDF weighting.

    When the ODE solver produces very sparse adaptive time points (e.g., 13
    points over [0, 1000] for a linear ODE), the CDF-weighted average is
    inaccurate because most CDF mass falls in a single coarse interval.
    This function interpolates the metrics onto a denser grid concentrated
    where the CDF has significant probability mass.

    Only activates when the ODE solver produces fewer than min_points time
    points. When the solver produces enough points, adaptive placement is
    already adequate for CDF weighting.

    Args:
        t_vals: Original time points from ODE solver
        q_metric: Queue length metric values at t_vals
        u_tran_ir: Utilization transient result (TranResult or dict)
        t_tran_ir: Throughput transient result (TranResult or dict)
        dist: Transition distribution (for determining CDF time scale)
        n_interp: Number of interpolation points to create
        min_points: Minimum ODE points before interpolation is skipped

    Returns:
        Tuple of (t_fine, q_fine, u_fine, t_fine_data) where each is an array
        on the finer grid, or None for u/t if not available.
    """
    if len(t_vals) >= n_interp:
        # Already denser than interpolation target — no interpolation needed
        _, u_metric = _get_tran_data(u_tran_ir)
        _, t_metric = _get_tran_data(t_tran_ir)
        return t_vals, q_metric, u_metric, t_metric

    # Determine CDF time scale from the distribution mean
    rate = _get_rate(dist)
    if rate > 0:
        mean_sojourn = 1.0 / rate
    else:
        mean_sojourn = (t_vals[-1] - t_vals[0]) / 10.0

    # Build a fine grid concentrated where the CDF has mass (up to ~5 mean sojourns)
    # but also covering the full time range for completeness
    t_cdf_end = min(t_vals[-1], 5.0 * mean_sojourn)
    if t_cdf_end <= t_vals[0]:
        t_cdf_end = t_vals[-1]

    # Dense grid in CDF-active region, sparse beyond
    n_dense = int(0.9 * n_interp)
    n_tail = n_interp - n_dense
    t_dense = np.linspace(t_vals[0], t_cdf_end, n_dense)
    if t_cdf_end < t_vals[-1] and n_tail > 1:
        t_tail = np.linspace(t_cdf_end, t_vals[-1], n_tail + 1)[1:]  # Exclude overlap
        t_fine = np.concatenate([t_dense, t_tail])
    else:
        t_fine = t_dense

    # Interpolate metrics using linear interpolation
    q_fine = np.interp(t_fine, t_vals, q_metric)

    _, u_metric = _get_tran_data(u_tran_ir)
    u_fine = np.interp(t_fine, t_vals, u_metric) if u_metric is not None else None

    _, t_metric = _get_tran_data(t_tran_ir)
    t_fine_data = np.interp(t_fine, t_vals, t_metric) if t_metric is not None else None

    return t_fine, q_fine, u_fine, t_fine_data


def _interpolate_for_cdf_map(t_vals, q_metric, u_tran_ir, t_tran_ir, D0, D1,
                             n_interp=500, min_points=50):
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
        min_points: Minimum ODE points before interpolation is skipped

    Returns:
        Tuple of (t_fine, q_fine, u_fine, t_fine_data)
    """
    if len(t_vals) >= min_points:
        _, u_metric = _get_tran_data(u_tran_ir)
        _, t_metric = _get_tran_data(t_tran_ir)
        return t_vals, q_metric, u_metric, t_metric

    # Compute mean sojourn from MAP: mean = 1 / (alpha @ D1 @ e)
    alpha = _map_prob(D0, D1)
    e = np.ones(D0.shape[0])
    total_rate = float(alpha @ D1 @ e)
    if total_rate > 0:
        mean_sojourn = 1.0 / total_rate
    else:
        mean_sojourn = (t_vals[-1] - t_vals[0]) / 10.0

    t_cdf_end = min(t_vals[-1], 5.0 * mean_sojourn)
    if t_cdf_end <= t_vals[0]:
        t_cdf_end = t_vals[-1]

    n_dense = int(0.9 * n_interp)
    n_tail = n_interp - n_dense
    t_dense = np.linspace(t_vals[0], t_cdf_end, n_dense)
    if t_cdf_end < t_vals[-1] and n_tail > 1:
        t_tail = np.linspace(t_cdf_end, t_vals[-1], n_tail + 1)[1:]
        t_fine = np.concatenate([t_dense, t_tail])
    else:
        t_fine = t_dense

    q_fine = np.interp(t_fine, t_vals, q_metric)

    _, u_metric = _get_tran_data(u_tran_ir)
    u_fine = np.interp(t_fine, t_vals, u_metric) if u_metric is not None else None

    _, t_metric = _get_tran_data(t_tran_ir)
    t_fine_data = np.interp(t_fine, t_vals, t_metric) if t_metric is not None else None

    return t_fine, q_fine, u_fine, t_fine_data


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

    def add_stage(self, index: int, name: str, stage_type: str, model: Any) -> None:
        """
        Add or update a stage in the environment.

        Args:
            index: Stage index (0-based)
            name: Name of the stage
            stage_type: Type of stage ('UP', 'DOWN', etc.)
            model: Network model for this stage
        """
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

    def add_transition(self, from_stage: int, to_stage: int, distribution: Any,
                       reset_rule: Optional[Callable[[np.ndarray], np.ndarray]] = None) -> None:
        """
        Add a transition between stages with an optional reset rule.

        Args:
            from_stage: Source stage index
            to_stage: Destination stage index
            distribution: Distribution for the transition time (e.g., Exp(rate))
            reset_rule: Optional function that transforms queue lengths when transition occurs.
        """
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
            # Start with the first valid transition from e
            hold_mmap = None
            for h in range(E):
                if h != e and self.env[e][h] is not None:
                    if hold_mmap is None:
                        hold_mmap = [m.copy() for m in emmap[e][h]]
                    else:
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

            if hold_mmap is None:
                # No outgoing transitions from this stage
                hold_mmap = [np.zeros((1, 1)), np.zeros((1, 1))] + [np.zeros((1, 1)) for _ in range(E)]

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
    def addStage(self, index: int, name: str, stage_type: str, model: Any) -> None:
        """Alias for add_stage (MATLAB compatibility)."""
        return self.add_stage(index, name, stage_type, model)

    def addTransition(self, from_stage: int, to_stage: int, distribution: Any,
                      reset_rule: Optional[Callable[[np.ndarray], np.ndarray]] = None) -> None:
        """Alias for add_transition (MATLAB compatibility)."""
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
        if options is None:
            return {
                'method': 'default',
                'iter_max': 100,
                'iter_tol': 1e-4,
                'verbose': False,
                'config': None
            }

        if isinstance(options, dict):
            opts = {
                'method': options.get('method', 'default'),
                'iter_max': options.get('iter_max', 100),
                'iter_tol': options.get('iter_tol', 1e-4),
                'verbose': options.get('verbose', False),
                'sojourn': options.get('sojourn', None),
                'config': options.get('config', None)
            }
            return opts

        # Options object
        return {
            'method': getattr(options, 'method', 'default'),
            'iter_max': getattr(options, 'iter_max', 100),
            'iter_tol': getattr(options, 'iter_tol', 1e-4),
            'verbose': getattr(options, 'verbose', False),
            'sojourn': getattr(options, 'sojourn', None),
            'config': getattr(options, 'config', None)
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

        Note: the compressed re-solve of the ensemble (weighted-average
        macro-state networks) is not applied in the native Python solver; the
        NCD analysis is exposed via compressionResult for inspection.

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
        return self.compressionResult

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
                    except Exception:
                        pass

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
            pass

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

            for h in range(E):
                Qexit[(e, h)] = np.zeros((M, K))
                Uexit[(e, h)] = np.zeros((M, K))
                Texit[(e, h)] = np.zeros((M, K))

                # Get MAP representation of transition distribution for CDF weighting
                # JAR uses: map_cdf(proc[e][h].get(0), proc[e][h].get(1), t)
                dist_eh = self.env_model.proc[e][h]
                if dist_eh is None:
                    continue

                D0_eh, D1_eh = _get_map_representation(dist_eh)
                # Normalize MAP
                Q_map = D0_eh + D1_eh
                for row in range(D0_eh.shape[0]):
                    D0_eh[row, row] = 0
                    D0_eh[row, row] = -np.sum(D0_eh[row, :]) - np.sum(D1_eh[row, :])

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
                                U_tran[i][r], T_tran[i][r], D0_eh, D1_eh)

                        # Use MAP CDF for weighting (matching JAR)
                        cdf_vals = _map_eval_cdf(D0_eh, D1_eh, t_fine)
                        w = np.zeros(len(t_fine))
                        w[1:] = cdf_vals[1:] - cdf_vals[:-1]

                        w_sum = np.sum(w)
                        if w_sum > 0 and not np.any(np.isnan(w)):
                            Qexit[(e, h)][i, r] = np.dot(q_fine, w) / w_sum
                            if u_fine is not None:
                                Uexit[(e, h)][i, r] = np.dot(u_fine, w) / w_sum
                            if t_fine_data is not None:
                                Texit[(e, h)][i, r] = np.dot(t_fine_data, w) / w_sum
                        else:
                            # Fall back to final value
                            Qexit[(e, h)][i, r] = metric_vals[-1] if len(metric_vals) > 0 else 0.0

        # Store exit metrics for convergence check
        self._Qexit = Qexit

        # Compute entry marginals using reset functions
        for e in range(E):
            result_e = self.results.get((it, e))
            if result_e is None or result_e['Tran']['Avg']['Q'] is None:
                continue

            Q_tran = result_e['Tran']['Avg']['Q']
            M = len(Q_tran)
            K = len(Q_tran[0]) if M > 0 else 0
            Qentry = np.zeros((M, K))

            for h in range(E):
                if (h, e) in Qexit and self.env_model.probOrig[h, e] > 0:
                    reset_fn = self.env_model.resetFun[h][e]
                    Qentry += self.env_model.probOrig[h, e] * reset_fn(Qexit[(h, e)])

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
            QExit[e] = np.zeros((M, K))
            UExit[e] = np.zeros((M, K))
            TExit[e] = np.zeros((M, K))

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

        # Cache-hit aggregation for fluid inner solvers.
        self._aggregate_cache_meanfield()

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
        from .api.cache import cache_gamma_lp, cache_miss_rmf
        from .api.sn.network_struct import NodeType
        from .lang.base import ReplacementStrategy

        E = self.getNumberOfModels()
        if E == 0:
            return
        # see _kb/06-solver-catalog.md (ENV: Python environment.py additional notes) for rationale
        for e in range(E):
            if type(self.getSolver(e)).__name__ not in ('SolverFLD', 'FLD', 'SolverFluid'):
                return

        ref = self.ensemble[0]
        sn0 = ref.get_struct()
        caches = [ind for ind in range(sn0.nnodes) if sn0.nodetype[ind] == NodeType.CACHE]
        if not caches:
            return
        K = sn0.nclasses

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
                Rcost = getattr(ch, 'accost', None)
                if Rcost is None:
                    Rcost = [[_default_routing(h) for _ in range(n)] for _ in range(u)]
                gamma, _, _, _ = cache_gamma_lp(lam, Rcost)
                infos.append(dict(node=ind, gamma=gamma, m=m, lam=lam, arate=arate,
                                  strat=getattr(ch, 'replacestrat', None)))
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

        # Only RANDOM(m) has a drift-based transient.
        for e in range(E):
            for si in stage_info[e]:
                if si['strat'] != ReplacementStrategy.RR:
                    return

        ncaches = len(caches)
        entry = [[None] * ncaches for _ in range(E)]
        Hp = [None] * E
        Mp = [None] * E
        W = [None] * E
        Exit = [None] * E
        max_sweep = max(1, int(self.options.get('iter_max', 100)))
        tol = float(self.options.get('iter_tol', 1e-4))
        prev = None
        for _sweep in range(max_sweep):
            for e in range(E):
                hd = self.env_model.holdTime[e]
                D0, D1 = hd[0], hd[1]
                Hp_e = None
                Mp_e = None
                Exit_e = [None] * ncaches
                W_e = [None] * ncaches
                for cc, si in enumerate(stage_info[e]):
                    x0 = entry[e][cc]
                    _, _, _, _, tout, _, MU_t, xtraj = cache_miss_rmf(
                        si['gamma'], si['m'], si['lam'], tspan=tspan[e], x0init=x0)
                    if Hp_e is None:
                        nt = len(tout)
                        Hp_e = np.zeros((ncaches, K, nt))
                        Mp_e = np.zeros((ncaches, K, nt))
                    lam = si['lam']
                    for v in range(lam.shape[0]):
                        rr = float(np.sum(lam[v, :, 0]))
                        if rr > 0:
                            mpv = np.clip(MU_t[v, :] / rr, 0.0, 1.0)
                            Mp_e[cc, v, :] = mpv
                            Hp_e[cc, v, :] = 1.0 - mpv
                    # see _kb/06-solver-catalog.md (ENV: Python environment.py additional notes) for rationale
                    Lam = float(np.sum(si['arate']))
                    if Lam > 0:
                        treal = tout / Lam
                        cdf = _map_eval_cdf(D0, D1, treal)
                        w = np.zeros(len(tout))
                        w[1:] = cdf[1:] - cdf[:-1]
                        W_e[cc] = w
                        sw = np.sum(w)
                        if sw > 0:
                            Exit_e[cc] = xtraj @ w / sw
                Hp[e], Mp[e], W[e], Exit[e] = Hp_e, Mp_e, W_e, Exit_e
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
        for cc, ind in enumerate(caches):
            hitT = np.zeros(K)
            missT = np.zeros(K)
            for e in range(E):
                w_e = W[e][cc] if (W[e] is not None and cc < len(W[e])) else None
                sw = np.sum(w_e) if w_e is not None else 0.0
                if not (sw > 0):
                    continue
                pe = self.env_model.probEnv[e]
                for k in range(K):
                    a = stage_info[e][cc]['arate'][k]
                    if a <= 0:
                        continue
                    hbar = Hp[e][cc, k, :] @ w_e / sw
                    mbar = Mp[e][cc, k, :] @ w_e / sw
                    hitT[k] += pe * a * hbar
                    missT[k] += pe * a * mbar
            hitprob = np.full(K, np.nan)
            missprob = np.full(K, np.nan)
            for k in range(K):
                tot = hitT[k] + missT[k]
                if tot > 0:
                    hitprob[k] = hitT[k] / tot
                    missprob[k] = missT[k] / tot
            node = ref_nodes[ind]
            if hasattr(node, 'set_result_hit_prob'):
                node.set_result_hit_prob(hitprob)
            if hasattr(node, 'set_result_miss_prob'):
                node.set_result_miss_prob(missprob)

    def iterate(self):
        """Run the main iteration loop.

        Mirrors JAR SolverENV.blending():
        - init() then pre(1) for initial steady-state
        - Loop: analyze all stages, post, converged
        - If max_iter hit without convergence: average last 10% of iterations
        - finish() for CDF-weighted aggregation
        """
        self.init()

        it = 0
        iter_max = self.options.get('iter_max', 100)
        verbose = self.options.get('verbose', False)

        # State-vector analyzer: carries the full per-stage distribution across
        # environment switches (mirrors MATLAB/JAR options.method='statevec').
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
            QN, UN, RN, TN = solver_env_statevec(self.env_model, self._solvers, self.options)
            self.result = {'Avg': {'Q': QN, 'U': UN, 'T': TN, 'R': RN}}
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

    def avg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute average performance metrics across environments.

        Mirrors MATLAB SolverENV.getEnsembleAvg() which calls iterate()
        and returns environment-weighted Q, U, T from self.result.Avg.

        Returns:
            Tuple of (QN, UN, TN) - queue lengths, utilizations, throughputs
        """
        if self.result is None:
            self.iterate()

        if self.result is None:
            return np.array([]), np.array([]), np.array([])

        Q = self.result['Avg']['Q']
        U = self.result['Avg']['U']
        T = self.result['Avg']['T']

        self._result = (Q, U, T)
        return Q, U, T

    def getAvg(self):
        """Get average metrics (MATLAB compatibility).

        Returns (QN, UN, TN) matching MATLAB getAvg -> getEnsembleAvg.
        """
        return self.avg()

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
            'verbose': False
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

        Matches JAR SolverENV.listValidMethods().
        """
        return ['default', 'smp']

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
