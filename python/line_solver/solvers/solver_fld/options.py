"""
FLD Solver Configuration and Results.

Defines configuration parameters (SolverFLDOptions) and result structures (FLDResult)
for SolverFLD fluid approximation solver.

This module provides dataclasses for:
- Configuring solution method, ODE tolerances, smoothing parameters
- Storing performance metrics (queue lengths, utilizations, response times)
- Accessing transient analysis data (time-dependent metrics)
- Tracking solution statistics (runtime, iterations)

Usage Examples
--------------

Basic configuration with defaults:
    >>> opts = SolverFLDOptions()
    >>> solver = SolverFLD(network, options=opts)

Custom configuration for higher accuracy:
    >>> opts = SolverFLDOptions(
    ...     method='matrix',
    ...     tol=1e-6,           # Tighter ODE tolerance
    ...     pstar=100,          # Smoother p-norm approximation
    ...     verbose=True
    ... )
    >>> solver = SolverFLD(network, options=opts).runAnalyzer()

Closing method with custom iterations:
    >>> opts = SolverFLDOptions(
    ...     method='closing',
    ...     iter_max=500,        # More iterations for convergence
    ...     iter_tol=1e-6        # Tighter iteration tolerance
    ... )
    >>> solver = SolverFLD(network, options=opts).runAnalyzer()

Diffusion method for stochastic analysis:
    >>> opts = SolverFLDOptions(
    ...     method='diffusion',
    ...     timestep=0.001,      # Small SDE time steps
    ...     tol=1e-5
    ... )
    >>> solver = SolverFLD(network, options=opts).runAnalyzer()
    >>> result = solver.result
    >>> print(f"Queue trajectory available: {len(result.QNt)} points")

Accessing results:
    >>> solver = SolverFLD(network).runAnalyzer()
    >>> result = solver.result
    >>> print(f"Queue lengths: {result.QN}")
    >>> print(f"Utilizations: {result.UN}")
    >>> print(f"Method used: {result.method}")
    >>> print(f"Runtime: {result.runtime:.3f}s")
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple
import numpy as np


@dataclass
class SolverFLDOptions:
    """Configuration options for FLD solver.

    Comprehensive configuration for SolverFLD, controlling solution method,
    ODE integration parameters, smoothing strategies, and output verbosity.

    Attributes
    ----------
    method : str, default='default'
        Solution method to use. Valid options:
        - 'default', 'matrix', 'fluid.matrix', 'pnorm': Matrix method (recommended)
        - 'softmin', 'fluid.softmin': Softmin smoothing
        - 'statedep', 'fluid.statedep': State-dependent constraints
        - 'closing', 'fluid.closing': Closing approximation with FCFS iteration
        - 'diffusion', 'fluid.diffusion': Euler-Maruyama stochastic solver
        - 'mfq', 'fluid.mfq', 'butools': Markovian fluid queue (single-queue only)

    tol : float, default=1e-4
        Tolerance for ODE integration (relative and absolute tolerance).
        Smaller values (e.g., 1e-6) increase accuracy but runtime.
        Recommended range: [1e-5, 1e-3]

    iter_max : int, default=200
        Maximum iterations for closing method outer loop (FCFS approximation).
        Increase if convergence is not achieved.

    iter_tol : float, default=1e-4
        Convergence tolerance for closing method iteration.
        Controls convergence: max(|eta - eta_prev| / eta) < iter_tol

    stiff : bool, default=True
        Enable LSODA stiffness detection for ODE solver.
        Automatically switches between Adams (non-stiff) and BDF (stiff) methods.
        Recommended to keep True for robust solution.

    timespan : tuple of float, default=(0.0, inf)
        Time interval for ODE integration (t_start, t_end).
        For steady-state analysis, t_end should be large (e.g., 1e6).
        Set t_end to finite value for transient analysis.

    timestep : float or None, default=None
        Fixed time step for diffusion method (Euler-Maruyama SDE).
        If None, automatic stepping is used.
        Smaller values increase accuracy but runtime.
        Recommended: 0.001-0.01 for diffusion method.

    pstar : float, default=20.0
        P-norm smoothing parameter for matrix method.
        Controls smoothness of capacity constraint: min(1, c/x) ≈ ghat(x, c, p)
        - p=2: Smooth approximation (fast but less accurate)
        - p=20: Default, good balance
        - p=100: Very smooth but slower
        Effect: lim_{p→∞} ghat(x, c, p) → min(1, c/x)

    softmin_alpha : float, default=20.0
        Softmin smoothing parameter for softmin/statedep methods.
        Controls smoothness via log-sum-exp: softmin(x, c, α) ≈ (log(e^(-αx) + e^(-αc)))/(-α)
        Larger values → sharper transitions, smaller values → smoother approximation.
        Recommended range: [1, 100]

    hide_immediate : bool, default=True
        Eliminate immediate transitions in network analysis.
        Recommended to keep True for most networks.

    verbose : bool, default=False
        Print detailed progress information during solution.
        Useful for debugging, includes method name and execution time.

    Notes
    -----
    Parameter Selection Guidelines:

    For **accurate solutions** (high precision):
        - tol = 1e-6 (tighter ODE tolerance)
        - pstar = 100 (smoother approximation)
        - iter_max = 500 (more closing iterations)

    For **fast solutions** (quick approximations):
        - tol = 1e-3 (looser tolerance)
        - pstar = 5 (faster but less smooth)
        - iter_max = 100 (fewer iterations)

    For **diffusion (stochastic) method**:
        - timestep = 0.001 (small SDE time step)
        - tol = 1e-5 (moderate ODE accuracy)
        - iter_max and iter_tol are ignored

    For **closing method** (FCFS networks):
        - iter_max = 200-500 (increase if not converging)
        - iter_tol = 1e-4 (controls iteration precision)

    For **MFQ method** (single queue):
        - All parameters except method are ignored
        - Analytical solution (no ODE integration)

    Examples
    --------
    Fast approximation:
        >>> opts = SolverFLDOptions(tol=1e-3, pstar=5)

    High accuracy:
        >>> opts = SolverFLDOptions(tol=1e-6, pstar=100)

    Diffusion analysis with verbose output:
        >>> opts = SolverFLDOptions(
        ...     method='diffusion',
        ...     timestep=0.01,
        ...     verbose=True
        ... )
    """
    method: str = 'default'
    tol: float = 1e-4
    iter_max: int = 200
    iter_tol: float = 1e-4
    stiff: bool = True
    timespan: Tuple[float, float] = (0.0, float('inf'))
    timestep: Optional[float] = None
    pstar: float = 20.0
    softmin_alpha: float = 20.0
    hide_immediate: bool = True
    verbose: bool = False
    # Common solver parameters (for compatibility)
    seed: Optional[int] = None  # Random seed (for compatibility)
    keep: bool = False  # Keep intermediate data (for compatibility)
    cutoff: Optional[int] = None  # State space cutoff (for compatibility)
    samples: Optional[int] = None  # Samples (for compatibility)


@dataclass
class FLDResult:
    """Result container for FLD solver analysis.

    Comprehensive result structure containing steady-state and transient metrics
    from fluid approximation analysis, including queue lengths, utilizations,
    response times, and throughputs across all stations and job classes.

    Attributes
    ----------
    QN : np.ndarray, shape (M, K)
        Average queue lengths per station and class.
        QN[i, k] = average number of customers at station i in class k
        (includes both waiting and in-service customers)

    UN : np.ndarray, shape (M, K)
        Average server utilization per station and class.
        UN[i, k] = fraction of time server at station i is busy with class k
        Range: [0, 1] for single-server, [0, c] for c-server queues
        Interpretation: ρ = λ / (c × μ) for M/M/c queue

    RN : np.ndarray, shape (M, K)
        Average response time (waiting + service) per station and class.
        RN[i, k] = mean time customer spends at station i (in time units)
        Related to queue length via Little's Law: QN[i, k] ≈ TN[i, k] × RN[i, k]

    TN : np.ndarray, shape (M, K)
        Average throughput per station and class.
        TN[i, k] = mean number of class-k customers per time unit at station i
        For open networks in equilibrium: TN[i, k] = λ_k (arrival rate)

    CN : np.ndarray, shape (1, K)
        System cycle time per class (sum of all station response times).
        CN[0, k] = Σ_i RN[i, k] = total time in system for class k
        For tandem networks: complete path response time

    XN : np.ndarray, shape (1, K)
        System throughput per class.
        XN[0, k] = overall throughput of class k through system
        For open networks: XN[0, k] = arrival rate λ_k
        For closed networks: determined by bottleneck service rate

    t : np.ndarray or None, optional
        Time points for transient analysis.
        t[n] = time value at ODE solution point n
        Used with QNt/UNt for time-dependent analysis
        None if only steady-state requested

    QNt : dict mapping (int, int) to np.ndarray, optional
        Transient queue length trajectories.
        Key: (station i, class k) → Value: array of queue lengths over time
        QNt[(i, k)][n] = queue length at station i, class k at time t[n]
        Available for methods with transient output

    UNt : dict mapping (int, int) to np.ndarray, optional
        Transient utilization trajectories.
        Key: (station i, class k) → Value: array of utilizations over time
        UNt[(i, k)][n] = utilization at station i, class k at time t[n]
        Available for methods with transient output

    xvec : np.ndarray or None, optional
        Final ODE state vector (internal representation).
        Contains raw queue length values from ODE solver.
        Length = Σ (M_i × K_i) where M_i = phases at station i, K_i = classes
        Useful for debugging and advanced analysis

    iterations : int, default=0
        Number of iterations executed (for iterative methods).
        - Matrix/MFQ/Diffusion: always 1 (not iterative)
        - Closing method: number of FCFS refinement iterations
        Useful for monitoring convergence

    runtime : float, default=0.0
        Execution time in seconds.
        Wall-clock time from start to end of solution.
        Useful for performance analysis

    method : str, default='matrix'
        Solution method used ('matrix', 'mfq', 'diffusion', 'closing', etc.)
        Identifies which algorithm produced these results

    Notes
    -----
    Array Dimensions:
    - M = number of stations
    - K = number of job classes
    - All metric arrays (QN, UN, RN, TN) are M × K

    Little's Law Validation:
    QN[i, k] ≈ TN[i, k] × RN[i, k]

    This should hold approximately (within numerical tolerance) for all
    stations and classes in stable systems.

    Transient Analysis:
    For methods supporting transient output (matrix, closing):
    - QNt and UNt provide time-dependent queue and utilization values
    - Access via: result.QNt[(i, k)] for station i, class k trajectory
    - Useful for analyzing system dynamics and convergence to steady-state

    Examples
    --------
    Accessing steady-state metrics:
        >>> solver = SolverFLD(network).runAnalyzer()
        >>> result = solver.result
        >>> qlen = result.QN  # (M, K) array
        >>> resp_time = result.RN  # (M, K) array
        >>> sys_resp = result.CN  # (1, K) array

    Analyzing bottleneck station:
        >>> util = result.UN
        >>> bottleneck_idx = np.argmax(np.mean(util, axis=1))
        >>> print(f"Bottleneck: Station {bottleneck_idx}")
        >>> print(f"Utilization: {np.mean(util[bottleneck_idx, :]):.1%}")

    Verifying Little's Law:
        >>> for i in range(result.QN.shape[0]):
        ...     for k in range(result.QN.shape[1]):
        ...         q = result.QN[i, k]
        ...         t = result.TN[i, k]
        ...         r = result.RN[i, k]
        ...         error = abs(q - t * r) / (t * r) if t > 0 else 0
        ...         print(f"Station {i}, Class {k}: error={error:.2%}")

    Accessing transient trajectory:
        >>> if result.QNt:  # Check if transient data available
        ...     qlen_traj = result.QNt[(0, 0)]  # Queue 0, class 0
        ...     import matplotlib.pyplot as plt
        ...     plt.plot(result.t, qlen_traj)
        ...     plt.xlabel('Time')
        ...     plt.ylabel('Queue Length')
        ...     plt.show()
    """
    QN: np.ndarray
    UN: np.ndarray
    RN: np.ndarray
    TN: np.ndarray
    CN: np.ndarray
    XN: np.ndarray
    t: Optional[np.ndarray] = None
    QNt: Dict[Tuple[int, int], np.ndarray] = field(default_factory=dict)
    UNt: Dict[Tuple[int, int], np.ndarray] = field(default_factory=dict)
    xvec: Optional[np.ndarray] = None
    iterations: int = 0
    runtime: float = 0.0
    method: str = 'matrix'
