"""
ENV Solver handler.

Native Python implementation of ENV (Ensemble Environment) solver handler
that analyzes queueing networks immersed in random environments.

The ENV solver handles CTMC-modulated queueing systems where the network
parameters change according to a continuous-time Markov chain environment.

Port from:

"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Any, Callable
import time
from scipy import linalg

from ...sn import NetworkStruct, SchedStrategy, NodeType
from ...mc import ctmc_solve


@dataclass
class SolverENVOptions:
    """Options for ENV solver."""
    method: str = 'default'
    tol: float = 1e-6
    verbose: bool = False
    iter_max: int = 100
    iter_tol: float = 1e-4
    state_dep_method: str = 'stateindep'  # 'stateindep' or 'statedep'


@dataclass
class SolverENVReturn:
    """
    Result of ENV solver handler.

    Attributes:
        Q: Mean queue lengths (M x K)
        U: Utilizations (M x K)
        R: Response times (M x K)
        T: Throughputs (M x K)
        C: Cycle times (1 x K)
        X: System throughputs (1 x K)
        pi: Environment stationary distribution
        runtime: Runtime in seconds
        method: Method used
        it: Number of iterations
    """
    Q: Optional[np.ndarray] = None
    U: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None
    T: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    X: Optional[np.ndarray] = None
    pi: Optional[np.ndarray] = None
    runtime: float = 0.0
    method: str = "default"
    it: int = 0


def _ctmc_makeinfgen(E: np.ndarray) -> np.ndarray:
    """
    Convert a rate matrix to a proper infinitesimal generator.

    Ensures row sums are zero by adjusting diagonal elements.

    Args:
        E: Rate matrix (E x E)

    Returns:
        Infinitesimal generator with zero row sums
    """
    Q = E.copy()
    n = Q.shape[0]
    for i in range(n):
        Q[i, i] = -np.sum(Q[i, :]) + Q[i, i]
    return Q


def _compute_env_stationary(env_rates: np.ndarray) -> np.ndarray:
    """
    Compute stationary distribution of environment CTMC.

    Args:
        env_rates: Environment transition rate matrix (E x E)

    Returns:
        Stationary distribution vector (1 x E)
    """
    Q = _ctmc_makeinfgen(env_rates)
    pi = ctmc_solve(Q)
    return pi.reshape(1, -1)


def solver_env_basic(
    stages: List[NetworkStruct],
    env_rates: np.ndarray,
    stage_results: List[Dict[str, np.ndarray]],
    options: Optional[SolverENVOptions] = None
) -> SolverENVReturn:
    """
    Basic ENV solver using probability blending.

    Computes performance metrics by blending stage results according to
    the stationary distribution of the environment CTMC.

    This is the simplest approach that assumes stages are independent
    and uses weighted averaging based on environment state probabilities.

    Args:
        stages: List of NetworkStruct for each environment stage
        env_rates: Environment transition rate matrix (E x E)
        stage_results: Pre-computed results for each stage
            Each dict should contain 'Q', 'U', 'R', 'T' matrices
        options: Solver options

    Returns:
        SolverENVReturn with blended performance metrics
    """
    start_time = time.time()

    if options is None:
        options = SolverENVOptions()

    E = len(stages)  # Number of environment states
    if E == 0:
        raise ValueError("At least one stage is required")

    # Get dimensions from first stage
    M = stages[0].nstations
    K = stages[0].nclasses

    # Compute environment stationary distribution
    pi = _compute_env_stationary(env_rates)

    # Initialize result matrices
    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))

    # Blend results across stages
    for e in range(E):
        if e >= len(stage_results):
            continue

        weight = pi[0, e]
        sr = stage_results[e]

        if 'Q' in sr and sr['Q'] is not None:
            Q_e = np.asarray(sr['Q'])
            if Q_e.shape == (M, K):
                QN += weight * Q_e

        if 'U' in sr and sr['U'] is not None:
            U_e = np.asarray(sr['U'])
            if U_e.shape == (M, K):
                UN += weight * U_e

        if 'R' in sr and sr['R'] is not None:
            R_e = np.asarray(sr['R'])
            if R_e.shape == (M, K):
                RN += weight * R_e

        if 'T' in sr and sr['T'] is not None:
            T_e = np.asarray(sr['T'])
            if T_e.shape == (M, K):
                TN += weight * T_e

    # Compute cycle times and system throughput
    CN = np.sum(RN, axis=0).reshape(1, -1)
    XN = np.zeros((1, K))
    for k in range(K):
        ref_stat = int(stages[0].refstat[k]) if k < len(stages[0].refstat) else 0
        if ref_stat < M:
            XN[0, k] = TN[ref_stat, k]

    result = SolverENVReturn()
    result.Q = QN
    result.U = UN
    result.R = RN
    result.T = TN
    result.C = CN
    result.X = XN
    result.pi = pi
    result.runtime = time.time() - start_time
    result.method = "basic"
    result.it = 1

    return result


def solver_env_statedep(
    stages: List[NetworkStruct],
    env_rates: np.ndarray,
    stage_solver: Callable[[NetworkStruct, Dict], Dict[str, np.ndarray]],
    options: Optional[SolverENVOptions] = None
) -> SolverENVReturn:
    """
    State-dependent ENV solver with iterative blending.

    This method accounts for state-dependence by iteratively solving
    stages with queue length estimates passed between stages.

    Args:
        stages: List of NetworkStruct for each environment stage
        env_rates: Environment transition rate matrix (E x E)
        stage_solver: Function that solves a single stage
            Should accept (sn, options_dict) and return metrics dict
        options: Solver options

    Returns:
        SolverENVReturn with blended performance metrics
    """
    start_time = time.time()

    if options is None:
        options = SolverENVOptions()

    E = len(stages)
    M = stages[0].nstations
    K = stages[0].nclasses

    # Compute environment stationary distribution
    pi = _compute_env_stationary(env_rates)

    # Compute embedding weights (probability of entering stage e from stage k)
    Q = _ctmc_makeinfgen(env_rates)
    embweight = np.zeros((E, E))
    for e in range(E):
        denom = 0.0
        for h in range(E):
            if h != e:
                denom += pi[0, h] * Q[h, e]
        if denom > 0:
            for k in range(E):
                if k != e:
                    embweight[k, e] = pi[0, k] * Q[k, e] / denom

    # Initialize stage results
    stage_results = []
    for e in range(E):
        result = stage_solver(stages[e], {})
        stage_results.append(result)

    # Iterative refinement
    prev_Q = None
    it = 0

    while it < options.iter_max:
        it += 1

        # Blend current results
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((M, K))

        for e in range(E):
            weight = pi[0, e]
            sr = stage_results[e]

            if 'Q' in sr and sr['Q'] is not None:
                QN += weight * np.asarray(sr['Q'])
            if 'U' in sr and sr['U'] is not None:
                UN += weight * np.asarray(sr['U'])
            if 'R' in sr and sr['R'] is not None:
                RN += weight * np.asarray(sr['R'])
            if 'T' in sr and sr['T'] is not None:
                TN += weight * np.asarray(sr['T'])

        # Check convergence
        if prev_Q is not None:
            diff = np.max(np.abs(QN - prev_Q))
            if diff < options.iter_tol:
                break

        prev_Q = QN.copy()

        # Re-solve stages with updated initial conditions
        for e in range(E):
            # Compute entry queue length for this stage
            Q_entry = np.zeros((M, K))
            for h in range(E):
                if h != e and 'Q' in stage_results[h]:
                    Q_entry += embweight[h, e] * np.asarray(stage_results[h]['Q'])

            # Solve with entry conditions
            stage_results[e] = stage_solver(stages[e], {'Q_entry': Q_entry})

    # Compute final metrics
    CN = np.sum(RN, axis=0).reshape(1, -1)
    XN = np.zeros((1, K))
    for k in range(K):
        ref_stat = int(stages[0].refstat[k]) if k < len(stages[0].refstat) else 0
        if ref_stat < M:
            XN[0, k] = TN[ref_stat, k]

    result = SolverENVReturn()
    result.Q = QN
    result.U = UN
    result.R = RN
    result.T = TN
    result.C = CN
    result.X = XN
    result.pi = pi
    result.runtime = time.time() - start_time
    result.method = "statedep"
    result.it = it

    return result


def solver_env(
    stages: List[NetworkStruct],
    env_rates: np.ndarray,
    stage_results: Optional[List[Dict[str, np.ndarray]]] = None,
    stage_solver: Optional[Callable] = None,
    options: Optional[SolverENVOptions] = None
) -> SolverENVReturn:
    """
    Main ENV solver handler.

    Routes to appropriate method based on options and available inputs.

    Args:
        stages: List of NetworkStruct for each environment stage
        env_rates: Environment transition rate matrix (E x E)
        stage_results: Pre-computed results for each stage (for basic method)
        stage_solver: Function to solve individual stages (for statedep method)
        options: Solver options

    Returns:
        SolverENVReturn with blended performance metrics
    """
    if options is None:
        options = SolverENVOptions()

    method = options.method.lower()

    if method == 'statedep' and stage_solver is not None:
        return solver_env_statedep(stages, env_rates, stage_solver, options)
    elif stage_results is not None:
        return solver_env_basic(stages, env_rates, stage_results, options)
    else:
        # No stage results provided - need to solve stages first
        raise ValueError("Either stage_results or stage_solver must be provided")


__all__ = [
    'solver_env',
    'solver_env_basic',
    'solver_env_statedep',
    'SolverENVReturn',
    'SolverENVOptions',
]
