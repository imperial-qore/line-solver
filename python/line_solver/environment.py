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
    if len(t_vals) >= min_points:
        # Already dense enough — no interpolation needed
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

    def add_node_breakdown(self, base_model, node_or_name, breakdown_dist, down_service_dist,
                           reset_fun: Optional[Callable] = None):
        """
        Add UP and DOWN stages for a node that can break down.

        Args:
            base_model: The base network model with normal (UP) service rates
            node_or_name: Node object or name of the node that can break down
            breakdown_dist: Distribution for time until breakdown (UP->DOWN transition)
            down_service_dist: Service distribution when the node is down
            reset_fun: Optional function to reset queue lengths on breakdown (default: keep jobs)
        """
        node_name = self._get_node_name(node_or_name)

        if reset_fun is None:
            reset_fun = lambda q: q

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

    def add_node_repair(self, node_or_name, repair_dist, reset_fun: Optional[Callable] = None):
        """
        Add repair transition from DOWN to UP stage for a previously added breakdown.

        Args:
            node_or_name: Node object or name of the node that can be repaired
            repair_dist: Distribution for repair time (DOWN->UP transition)
            reset_fun: Optional function to reset queue lengths on repair (default: keep jobs)
        """
        node_name = self._get_node_name(node_or_name)

        if reset_fun is None:
            reset_fun = lambda q: q

        down_stage_name = f'DOWN_{node_name}'
        down_idx = self._find_stage_by_name(down_stage_name)
        up_idx = self._find_stage_by_name('UP')

        if down_idx < 0:
            raise ValueError(f'DOWN stage for node "{node_name}" not found. Call add_node_breakdown first.')
        if up_idx < 0:
            raise ValueError('UP stage not found. Call add_node_breakdown first.')

        # Add repair transition (DOWN -> UP)
        self.add_transition(down_idx, up_idx, repair_dist, reset_fun)

    def add_node_failure_repair(self, base_model, node_or_name, breakdown_dist, repair_dist,
                                down_service_dist, reset_breakdown: Optional[Callable] = None,
                                reset_repair: Optional[Callable] = None):
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

        if reset_breakdown is None:
            reset_breakdown = lambda q: q
        if reset_repair is None:
            reset_repair = lambda q: q

        self.add_node_breakdown(base_model, node_name, breakdown_dist, down_service_dist, reset_breakdown)
        self.add_node_repair(node_name, repair_dist, reset_repair)

    def set_breakdown_reset_policy(self, node_or_name, reset_fun: Callable):
        """
        Update the reset policy for breakdown transitions (UP -> DOWN) of a node.

        Args:
            node_or_name: Node object or name of the node
            reset_fun: Function to reset queue lengths on breakdown
        """
        node_name = self._get_node_name(node_or_name)
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

    def set_repair_reset_policy(self, node_or_name, reset_fun: Callable):
        """
        Update the reset policy for repair transitions (DOWN -> UP) of a node.

        Args:
            node_or_name: Node object or name of the node
            reset_fun: Function to reset queue lengths on repair
        """
        node_name = self._get_node_name(node_or_name)
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

    def get_stage_table(self):
        """Get stage table (MATLAB compatibility alias)."""
        return self.stage_table()

    def get_ensemble(self) -> List[Any]:
        """Get list of network models (snake_case version)."""
        return self._models


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

    def _normalize_options(self, options) -> Dict:
        """Normalize options to a dictionary."""
        if options is None:
            return {
                'method': 'default',
                'iter_max': 100,
                'iter_tol': 1e-4,
                'verbose': False
            }

        if isinstance(options, dict):
            opts = {
                'method': options.get('method', 'default'),
                'iter_max': options.get('iter_max', 100),
                'iter_tol': options.get('iter_tol', 1e-4),
                'verbose': options.get('verbose', False)
            }
            return opts

        # Options object
        return {
            'method': getattr(options, 'method', 'default'),
            'iter_max': getattr(options, 'iter_max', 100),
            'iter_tol': getattr(options, 'iter_tol', 1e-4),
            'verbose': getattr(options, 'verbose', False)
        }

    def getNumberOfModels(self) -> int:
        """Get number of environment stages."""
        return self.env_model.num_stages

    def getSolver(self, e: int):
        """Get solver for stage e."""
        return self._solvers[e] if e < len(self._solvers) else None

    def init(self):
        """Initialize the environment solver."""
        self.env_model.init()
        self.results = {}

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

                        # Round fractional marginals for discrete solvers
                        solver_name = type(solver).__name__
                        if 'Fluid' not in solver_name and 'FLD' not in solver_name:
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

        # Get transient results
        # Note: timespan is NOT adjusted here — the CTMC solver's getTranAvg()
        # uses 30/minrate when timespan is unset, matching MATLAB's getTranAvg.m.
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
        """Post-iteration operations - compute exit metrics and update entry marginals."""
        E = self.getNumberOfModels()

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

                for i in range(M):
                    for r in range(K):
                        Qir = Q_tran[i][r]
                        t_vals, metric_vals = _get_tran_data(Qir)
                        if t_vals is not None:
                            # Compute weights using CDF of transition distribution
                            dist_eh = self.env_model.proc[e][h]
                            if dist_eh is not None:
                                # Interpolate transient data onto a finer grid for
                                # accurate CDF weighting. Python's LSODA may produce
                                # very few ODE time points (e.g., 13 for linear growth),
                                # but the CDF has most of its mass in a narrow range.
                                t_fine, q_fine, u_fine, t_fine_data = \
                                    _interpolate_for_cdf(
                                        t_vals, metric_vals,
                                        U_tran[i][r], T_tran[i][r], dist_eh)

                                cdf_vals = _eval_cdf(dist_eh, t_fine)
                                w = np.zeros(len(t_fine))
                                w[1:] = cdf_vals[1:] - cdf_vals[:-1]

                                if np.sum(w) > 0:
                                    Qexit[(e, h)][i, r] = np.dot(q_fine, w) / np.sum(w)
                                    if u_fine is not None:
                                        Uexit[(e, h)][i, r] = np.dot(u_fine, w) / np.sum(w)
                                    if t_fine_data is not None:
                                        Texit[(e, h)][i, r] = np.dot(t_fine_data, w) / np.sum(w)
                                else:
                                    # Use final value
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

            # Initialize model state from entry marginals and reset solver
            # MATLAB: self.solvers{e}.reset(); self.ensemble{e}.initFromMarginal(Qentry{e});
            model = self.ensemble[e] if e < len(self.ensemble) else None
            solver = self.getSolver(e)
            if solver is not None:
                if hasattr(solver, 'reset'):
                    solver.reset()
            if model is not None:
                # Round fractional marginals for discrete solvers
                solver_name = type(solver).__name__ if solver is not None else ''
                if 'Fluid' not in solver_name and 'FLD' not in solver_name:
                    Qentry = self._round_marginal_for_discrete_solver(Qentry, e)

                if hasattr(model, 'initFromMarginal'):
                    model.initFromMarginal(Qentry)
                elif hasattr(model, 'init_from_marginal'):
                    model.init_from_marginal(Qentry)

    def _round_marginal_for_discrete_solver(self, Q, stage_idx):
        """Round fractional queue lengths to integers using the largest remainder method,
        preserving closed chain populations exactly."""
        model = self.ensemble[stage_idx] if stage_idx < len(self.ensemble) else None
        if model is None:
            return Q
        sn = model.get_struct() if hasattr(model, 'get_struct') else None
        if sn is None or not hasattr(sn, 'chains') or not hasattr(sn, 'njobs'):
            return Q
        Q = Q.copy()
        nchains = sn.nchains if hasattr(sn, 'nchains') else (sn.chains.shape[0] if hasattr(sn.chains, 'shape') else 0)
        njobs = np.atleast_1d(sn.njobs).flatten()
        chains = np.atleast_2d(sn.chains)
        M = Q.shape[0]
        for c in range(nchains):
            chain_classes = np.where(chains[c, :] > 0)[0]
            njobs_chain = sum(njobs[k] for k in chain_classes)
            if np.isinf(njobs_chain):
                # Open chain: simple rounding
                for k in chain_classes:
                    for i in range(M):
                        Q[i, k] = round(Q[i, k])
            else:
                # Closed chain: largest remainder method
                vals = []
                indices = []
                for i in range(M):
                    for k in chain_classes:
                        vals.append(Q[i, k])
                        indices.append((i, k))
                vals = np.array(vals, dtype=float)
                floored = np.floor(vals)
                remainders = vals - floored
                deficit = int(round(njobs_chain - np.sum(floored)))
                if deficit > 0:
                    sort_idx = np.argsort(-remainders)
                    for d in range(min(deficit, len(sort_idx))):
                        floored[sort_idx[d]] += 1
                for j, (i, k) in enumerate(indices):
                    Q[i, k] = floored[j]
        return Q

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

    def iterate(self):
        """Run the main iteration loop."""
        self.init()

        it = 0
        iter_max = self.options.get('iter_max', 100)
        verbose = self.options.get('verbose', False)

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

                        # Build the full transition block
                        # Need to expand by phases of OTHER outgoing transitions from source stage
                        # and phases of outgoing transitions from target stage
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
        if model is not None:
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

    @staticmethod
    def default_options():
        """Get default solver options."""
        return {
            'method': 'default',
            'iter_max': 100,
            'iter_tol': 1e-4,
            'verbose': False
        }


# Convenience aliases
ENV = SolverENV


__all__ = [
    'Environment',
    'SolverENV',
    'ENV',
]
