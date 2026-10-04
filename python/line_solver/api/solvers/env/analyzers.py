"""
ENV Solver analyzers.

Native Python implementation of ENV solver analyzers that orchestrate
method selection and provide the main entry point for environment analysis.

Port from:

"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Callable
import time

from ...sn import NetworkStruct, SchedStrategy, NodeType
from .handler import (
    solver_env,
    solver_env_basic,
    solver_env_statedep,
    SolverENVOptions,
    SolverENVReturn,
)


@dataclass
class ENVResult:
    """
    Result of ENV solver analysis.

    Attributes:
        QN: Mean queue lengths (M x K)
        UN: Utilizations (M x K)
        RN: Response times (M x K)
        TN: Throughputs (M x K)
        CN: Cycle times (1 x K)
        XN: System throughputs (1 x K)
        pi: Environment stationary distribution
        iter: Number of iterations
        runtime: Runtime in seconds
        method: Method used
    """
    QN: Optional[np.ndarray] = None
    UN: Optional[np.ndarray] = None
    RN: Optional[np.ndarray] = None
    TN: Optional[np.ndarray] = None
    CN: Optional[np.ndarray] = None
    XN: Optional[np.ndarray] = None
    pi: Optional[np.ndarray] = None
    iter: int = 0
    runtime: float = 0.0
    method: str = ""


def solver_env_analyzer(
    stages: List[NetworkStruct],
    env_rates: np.ndarray,
    stage_results: Optional[List[Dict[str, np.ndarray]]] = None,
    stage_solver: Optional[Callable] = None,
    options: Optional[SolverENVOptions] = None
) -> ENVResult:
    """
    ENV Analyzer - main entry point for environment analysis.

    Analyzes queueing networks immersed in random environments by
    blending results across environment stages weighted by the
    stationary distribution of the environment CTMC.

    Supported methods:
        - 'default': State-independent blending
        - 'stateindep': State-independent blending (explicit)
        - 'statedep': State-dependent iterative blending

    Args:
        stages: List of NetworkStruct for each environment stage
        env_rates: Environment transition rate matrix (E x E)
        stage_results: Pre-computed results for each stage
            Each dict should contain 'Q', 'U', 'R', 'T' matrices
        stage_solver: Function to solve individual stages
            Required for state-dependent method
        options: Solver options

    Returns:
        ENVResult with all performance metrics

    Raises:
        ValueError: For unsupported configurations
    """
    start_time = time.time()

    if options is None:
        options = SolverENVOptions()

    method = options.method.lower()
    result = ENVResult()

    # Select and execute method
    if method in ['default', 'stateindep', 'basic']:
        if stage_results is None:
            raise ValueError("stage_results required for state-independent method")

        ret = solver_env_basic(stages, env_rates, stage_results, options)
        result.method = 'stateindep'

    elif method == 'statedep':
        if stage_solver is None:
            raise ValueError("stage_solver required for state-dependent method")

        ret = solver_env_statedep(stages, env_rates, stage_solver, options)
        result.method = 'statedep'

    else:
        # Unknown method - use basic if possible
        if stage_results is not None:
            if options.verbose:
                print(f"Warning: Unknown ENV method '{method}'. Using stateindep.")
            ret = solver_env_basic(stages, env_rates, stage_results, options)
            result.method = 'stateindep'
        elif stage_solver is not None:
            if options.verbose:
                print(f"Warning: Unknown ENV method '{method}'. Using statedep.")
            ret = solver_env_statedep(stages, env_rates, stage_solver, options)
            result.method = 'statedep'
        else:
            raise ValueError("Either stage_results or stage_solver must be provided")

    # Copy results
    if ret is not None:
        result.QN = ret.Q
        result.UN = ret.U
        result.RN = ret.R
        result.TN = ret.T
        result.CN = ret.C
        result.XN = ret.X
        result.pi = ret.pi
        result.iter = ret.it

    # Clean up NaN values
    if result.QN is not None:
        result.QN = np.nan_to_num(result.QN, nan=0.0)
    if result.UN is not None:
        result.UN = np.nan_to_num(result.UN, nan=0.0)
    if result.RN is not None:
        result.RN = np.nan_to_num(result.RN, nan=0.0)
    if result.TN is not None:
        result.TN = np.nan_to_num(result.TN, nan=0.0)
    if result.CN is not None:
        result.CN = np.nan_to_num(result.CN, nan=0.0)
    if result.XN is not None:
        result.XN = np.nan_to_num(result.XN, nan=0.0)

    result.runtime = time.time() - start_time

    return result


__all__ = [
    'ENVResult',
    'solver_env_analyzer',
]
