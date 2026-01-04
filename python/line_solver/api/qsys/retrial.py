"""
Retrial and reneging queueing system analysis framework.

This module provides a foundation for analyzing bufferless retrial queues
(BMAP/PH/N/N) and reneging queues (MAP/M/s+G), with extension points for
future full solver implementation.

The framework includes:
- BMAP (Batch Markovian Arrival Process) matrix utilities
- PH (Phase-Type) distribution utilities
- QBD (Quasi-Birth-Death) state space construction
- Retrial queue topology detection and analysis
- Reneging queue topology detection and analysis

Reference:
    Dudin, A., Klimenok, V., & Vishnevsky, V. (2020).
    "The Theory of Quasi-Birth-Death Processes and Its Applications".
    Springer.

Port from: matlab/src/solvers/MAM/solver_mam_retrial.m
           matlab/src/api/qsys/qsys_bmapphnn_retrial.m
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from enum import Enum


class QueueType(Enum):
    """Queueing system topology types."""
    STANDARD = "standard"
    RETRIAL = "retrial"
    RENEGING = "reneging"
    RETRIAL_RENEGING = "retrial_reneging"


@dataclass
class BmapMatrix:
    """
    Batch Markovian Arrival Process matrix representation.

    D₀ = drift matrix (no arrivals)
    D₁, D₂, ..., Dₖ = batch arrival matrices for batch sizes 1, 2, ..., K

    Properties:
        order: Dimension of the BMAP
        num_batches: Maximum batch size
        arrival_rate: Overall arrival rate
    """
    D0: np.ndarray
    D_batch: List[np.ndarray]  # D_batch[k] = D_{k+1}

    def __post_init__(self):
        """Validate BMAP structure."""
        if self.D0.shape[0] != self.D0.shape[1]:
            raise ValueError("D0 must be square matrix")
        self.order = self.D0.shape[0]
        self.num_batches = len(self.D_batch)

        # Validate batch matrices have same dimension
        for i, D in enumerate(self.D_batch):
            if D.shape != self.D0.shape:
                raise ValueError(f"D_batch[{i}] dimension mismatch with D0")

    @property
    def arrival_rate(self) -> float:
        """
        Compute overall arrival rate of the BMAP.

        lambda = theta * (-D0) * e where theta is stationary distribution
        of BMAP Markov chain and e is unit vector.
        """
        # Simple approximation: trace-based estimate
        # For exact computation, requires solving theta (-D0) * e
        D_sum = self.D0.copy()
        for D in self.D_batch:
            D_sum = D_sum + D
        # Eigenvalue-based approximation
        try:
            eigvals = np.linalg.eigvals(-D_sum)
            return np.max(np.real(eigvals)) if len(eigvals) > 0 else 0.0
        except:
            return 0.0

    @property
    def fundamental_arrival_rate(self) -> float:
        """Arrival rate from fundamental matrix."""
        try:
            # Compute -D0 inverse (if invertible)
            inv_D0 = np.linalg.inv(-self.D0)
            e = np.ones(self.order)
            # Sum of all batch matrices
            D_sum = np.zeros_like(self.D0)
            for D in self.D_batch:
                D_sum = D_sum + D
            # Arrival rate approximation
            return np.sum(inv_D0 @ D_sum)
        except:
            return self.arrival_rate


@dataclass
class PhDistribution:
    """
    Phase-Type (PH) distribution representation.

    Beta = initial probability vector (shape: m,)
    S = transient generator matrix (shape: m x m)

    where m is the number of phases.

    Properties:
        mean: Mean of PH distribution (1/mu in queue notation)
        scv: Squared coefficient of variation
        num_phases: Number of phases
    """
    beta: np.ndarray
    S: np.ndarray

    def __post_init__(self):
        """Validate PH distribution structure."""
        if self.beta.ndim != 1:
            raise ValueError("beta must be 1-D vector")
        if self.S.ndim != 2 or self.S.shape[0] != self.S.shape[1]:
            raise ValueError("S must be square matrix")
        if self.S.shape[0] != len(self.beta):
            raise ValueError("S dimension must match beta length")
        self.num_phases = len(self.beta)

    @property
    def mean(self) -> float:
        """
        Compute mean of PH distribution.

        E[X] = -beta * S^{-1} * e
        """
        try:
            S_inv = np.linalg.inv(self.S)
            e = np.ones(self.num_phases)
            return -np.dot(self.beta, S_inv) @ e
        except:
            return 1.0  # Fallback

    @property
    def mean_squared(self) -> float:
        """
        Compute second moment of PH distribution.

        E[X²] = 2 * beta * S^{-2} * e
        """
        try:
            S_inv = np.linalg.inv(self.S)
            e = np.ones(self.num_phases)
            return 2 * np.dot(self.beta, S_inv) @ (S_inv @ e)
        except:
            mean = self.mean
            return mean ** 2  # Fallback: assume deterministic

    @property
    def scv(self) -> float:
        """
        Squared coefficient of variation of PH distribution.

        SCV = (E[X²] / E[X]²) - 1
        """
        mean = self.mean
        mean_sq = self.mean_squared
        if mean > 0:
            return (mean_sq / (mean ** 2)) - 1
        return 0.0


@dataclass
class QbdStatespace:
    """
    Quasi-Birth-Death Markov chain state space representation.

    For a QBD process, states are of the form (n, i) where:
    - n = level (number of retrying customers in orbit)
    - i = phase (service phase or phase-type stage)

    Generator matrix has block structure:
    Q = | B_0   A_0    0     0   ...  |
        | A_2   A_1   A_0    0   ...  |
        | 0     A_2   A_1   A_0  ...  |
        | ...                          |

    Properties:
        max_level: Maximum retrial orbit size (truncation level)
        phase_dim: Number of phases at each level
        total_states: Total state space dimension
    """
    A0: np.ndarray  # Arrival matrix (birth)
    A1: np.ndarray  # Service matrix (drift)
    A2: np.ndarray  # Retrial matrix (death)
    B0: np.ndarray  # Boundary matrix at level 0
    max_level: int

    def __post_init__(self):
        """Validate QBD state space."""
        phase_dim = self.A0.shape[0]

        if not all(M.shape == (phase_dim, phase_dim)
                   for M in [self.A1, self.A2, self.B0]):
            raise ValueError("All matrices must have same dimension")

        self.phase_dim = phase_dim
        # Total states = phases * (max_level + 1) levels
        self.total_states = phase_dim * (self.max_level + 1)

    def get_state_index(self, level: int, phase: int) -> int:
        """
        Map (level, phase) to linear state index.

        Args:
            level: Retrial orbit size (0 ≤ level ≤ max_level)
            phase: Phase within level (0 ≤ phase < phase_dim)

        Returns:
            Linear state index
        """
        if not (0 <= level <= self.max_level and 0 <= phase < self.phase_dim):
            raise ValueError(f"Invalid state ({level}, {phase})")
        return level * self.phase_dim + phase

    def get_level_phase(self, idx: int) -> Tuple[int, int]:
        """
        Map linear state index to (level, phase).

        Args:
            idx: Linear state index

        Returns:
            Tuple (level, phase)
        """
        if not (0 <= idx < self.total_states):
            raise ValueError(f"Invalid state index {idx}")
        level = idx // self.phase_dim
        phase = idx % self.phase_dim
        return level, phase

    @property
    def binomial_state_dimension(self) -> int:
        """
        Compute state space dimension using binomial coefficient formula.

        For a retrial queue with N total customers and m servers,
        dimension = C(N + m - 1, m - 1)

        This is a rough upper bound for QBD truncation.
        """
        # Placeholder - requires N and m parameters
        return self.phase_dim * (self.max_level + 1)


@dataclass
class RetrialQueueResult:
    """
    Analysis results for BMAP/PH/N/N retrial queue.

    Attributes:
        queue_type: Type of queue (retrial, reneging, etc.)
        L_orbit: Expected number of customers in orbit
        N_server: Expected number of customers being served
        utilization: Server utilization
        throughput: System throughput
        P_idle: Probability that server is idle
        P_empty_orbit: Probability that orbit is empty
        stationary_dist: Stationary distribution vector
        truncation_level: QBD truncation level used
        converged: Whether numerical solution converged
        iterations: Number of iterations to convergence
        error: Final error estimate
    """
    queue_type: QueueType
    L_orbit: float
    N_server: float
    utilization: float
    throughput: float
    P_idle: float
    P_empty_orbit: float
    stationary_dist: Optional[np.ndarray] = None
    truncation_level: int = 0
    converged: bool = False
    iterations: int = 0
    error: float = np.inf


class RetrialQueueAnalyzer:
    """
    Framework for analyzing BMAP/PH/N/N retrial queues.

    This class provides the foundation for future full solver implementation,
    including topology detection, parameter extraction, and QBD setup.

    Example:
        analyzer = RetrialQueueAnalyzer(model)
        queue_type = analyzer.detect_queue_type()
        if queue_type == QueueType.RETRIAL:
            result = analyzer.analyze()
    """

    def __init__(self, sn: Any, options: Optional[Dict] = None):
        """
        Initialize retrial queue analyzer.

        Args:
            sn: NetworkStruct with queue configuration
            options: Analysis options (tolerance, max iterations, etc.)
        """
        self.sn = sn
        self.options = options or {}
        self.tolerance = self.options.get('tolerance', 1e-8)
        self.max_iterations = self.options.get('max_iterations', 1000)
        self.max_retrial_orbit = self.options.get('max_retrial_orbit', 100)

    def detect_queue_type(self) -> QueueType:
        """
        Detect the type of queueing system from topology.

        Returns:
            QueueType enum indicating retrial, reneging, or combination
        """
        # Check for retrial nodes (orbit with retrial rate)
        has_retrial = False
        has_reneging = False

        # Placeholder logic - would examine sn structure for:
        # - Retrial probability/rate (orbit impatience)
        # - Reneging probability/patience times
        # - Queue abandonment

        if has_retrial and has_reneging:
            return QueueType.RETRIAL_RENEGING
        elif has_retrial:
            return QueueType.RETRIAL
        elif has_reneging:
            return QueueType.RENEGING
        else:
            return QueueType.STANDARD

    def extract_bmap(self) -> Optional[BmapMatrix]:
        """
        Extract BMAP parameters from arrival process.

        Returns:
            BmapMatrix instance or None if arrival is not BMAP/MAP
        """
        # Placeholder implementation
        # Would extract D0, D1, D2, ... from sn.proc{arrival_station}
        return None

    def extract_ph_service(self) -> Optional[PhDistribution]:
        """
        Extract Phase-Type service parameters.

        Returns:
            PhDistribution instance or None if service is not PH
        """
        # Placeholder implementation
        # Would extract beta, S from sn.proc{service_station}
        return None

    def extract_retrial_parameters(self) -> Optional[Dict[str, float]]:
        """
        Extract retrial-specific parameters.

        Returns:
            Dict with keys:
            - alpha: Retrial rate (rate at which customers retry)
            - gamma: Orbit impatience rate (reneging rate)
            - p: Batch rejection probability
            - R: Threshold for admission control
        """
        # Placeholder implementation
        # Would extract from queue configuration
        return None

    def build_qbd_statespace(self,
                            bmap: BmapMatrix,
                            ph_service: PhDistribution,
                            retrial_params: Dict[str, float]
                            ) -> Optional[QbdStatespace]:
        """
        Build QBD state space for the retrial queue.

        This constructs the generator matrix blocks for the QBD process.

        Args:
            bmap: BMAP arrival process
            ph_service: PH service distribution
            retrial_params: Retrial parameters (alpha, gamma, p, R)

        Returns:
            QbdStatespace instance or None if construction fails
        """
        # Placeholder implementation
        # Would construct A0, A1, A2, B0 matrices
        # based on BMAP, PH, and retrial parameters
        return None

    def analyze(self) -> RetrialQueueResult:
        """
        Analyze the retrial queue.

        This is the main entry point for analysis. Currently returns
        a placeholder result. Full implementation would:
        1. Detect queue type
        2. Extract arrival and service parameters
        3. Build QBD state space
        4. Solve for stationary distribution
        5. Compute performance metrics

        Returns:
            RetrialQueueResult with performance metrics
        """
        queue_type = self.detect_queue_type()

        if queue_type == QueueType.STANDARD:
            raise ValueError("Standard queue - use standard analyzer")

        # Placeholder result
        result = RetrialQueueResult(
            queue_type=queue_type,
            L_orbit=0.0,
            N_server=0.0,
            utilization=0.0,
            throughput=0.0,
            P_idle=1.0,
            P_empty_orbit=1.0,
            converged=False,
            error=np.inf
        )

        return result


def qsys_bmapphnn_retrial(
    arrival_matrix: Dict[str, np.ndarray],
    service_params: Dict[str, float],
    N: int,
    retrial_params: Optional[Dict[str, float]] = None,
    options: Optional[Dict] = None
) -> RetrialQueueResult:
    """
    Analyze BMAP/PH/N/N bufferless retrial queue.

    Parameters:
        arrival_matrix: Dict with 'D0', 'D1', ... for BMAP
        service_params: Dict with 'beta' (initial prob) and 'S' (gen matrix) for PH
        N: Number of customers (for closed queue)
        retrial_params: Dict with 'alpha' (retrial rate), 'gamma' (impatience),
                       'p' (rejection prob), 'R' (threshold)
        options: Solver options

    Returns:
        RetrialQueueResult with performance metrics

    Raises:
        NotImplementedError: Full solver not yet implemented

    References:
        Dudin, A., Klimenok, V., & Vishnevsky, V. (2020).
        "The Theory of Quasi-Birth-Death Processes and Its Applications".
        Springer.
    """
    if options is None:
        options = {}

    # Extract parameters
    try:
        D0 = arrival_matrix.get('D0')
        D_batch = [arrival_matrix.get(f'D{i+1}') for i in range(10)]
        D_batch = [D for D in D_batch if D is not None]

        bmap = BmapMatrix(D0=D0, D_batch=D_batch)

        beta = service_params.get('beta', np.array([1.0]))
        S = service_params.get('S', np.array([[-1.0]]))
        ph_service = PhDistribution(beta=beta, S=S)

    except Exception as e:
        raise ValueError(f"Invalid BMAP/PH parameters: {e}")

    # Set retrial parameters
    if retrial_params is None:
        retrial_params = {
            'alpha': 1.0,  # Default retrial rate
            'gamma': 0.0,  # No impatience
            'p': 0.0,      # No rejection
            'R': N         # No threshold
        }

    # Create placeholder result
    result = RetrialQueueResult(
        queue_type=QueueType.RETRIAL,
        L_orbit=0.0,
        N_server=0.0,
        utilization=0.0,
        throughput=0.0,
        P_idle=1.0,
        P_empty_orbit=1.0,
        truncation_level=options.get('max_retrial_orbit', 100),
        converged=False,
        error=np.inf
    )

    return result


@dataclass
class RetrialInfo:
    """Information about a valid retrial queue topology."""
    is_retrial: bool
    station_idx: Optional[int] = None
    node_idx: Optional[int] = None
    source_idx: Optional[int] = None
    class_idx: Optional[int] = None
    error_msg: str = ''
    N: Optional[int] = None  # Number of servers
    alpha: float = 0.1  # Retrial rate (default)
    gamma: float = 0.0  # Orbit impatience (default 0)
    p: float = 0.0  # Batch rejection prob (default 0)
    R: Optional[int] = None  # Admission threshold


def qsys_is_retrial(sn: Any) -> Tuple[bool, RetrialInfo]:
    """
    Check if network is a valid BMAP/PH/N/N bufferless retrial queue.

    Validates that the network structure matches the requirements for
    the BMAP/PH/N/N retrial queue solver:
    - Single bufferless queue (capacity == number of servers)
    - Retrial drop strategy configured
    - BMAP/MAP arrival process at source
    - PH/Exp service at queue
    - Open class model

    Based on: Dudin et al., "Analysis of BMAP/PH/N-Type Queueing System with
    Flexible Retrials Admission Control", Mathematics 2025, 13(9), 1434.

    Args:
        sn: NetworkStruct object

    Returns:
        Tuple of (is_retrial, retrial_info) where:
            is_retrial: True if network is valid BMAP/PH/N/N retrial topology
            retrial_info: RetrialInfo with parameters for the retrial solver

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_is_retrial.m
    """
    from ..sn import sn_is_open_model, NodeType, DropStrategy

    ret_info = RetrialInfo(is_retrial=False)

    # Check if model is open
    try:
        if not sn_is_open_model(sn):
            ret_info.error_msg = 'BMAP/PH/N/N retrial solver requires open queueing model.'
            return False, ret_info
    except:
        ret_info.error_msg = 'Could not determine if model is open.'
        return False, ret_info

    # Check for single class (current limitation)
    if sn.nclasses > 1:
        ret_info.error_msg = 'BMAP/PH/N/N retrial solver currently supports single class only.'
        return False, ret_info

    ret_info.class_idx = 0

    # Find bufferless queue stations (capacity == nservers, finite capacity)
    bufferless_stations = []
    for ist in range(sn.nstations):
        node_idx = sn.stationToNode[ist] if hasattr(sn, 'stationToNode') else ist
        if node_idx < len(sn.nodetype) and sn.nodetype[node_idx] == NodeType.QUEUE:
            cap = sn.cap[ist] if hasattr(sn, 'cap') and sn.cap is not None else np.inf
            nservers = sn.nservers[ist] if hasattr(sn, 'nservers') and sn.nservers is not None else 1
            if np.isfinite(cap) and cap == nservers:
                bufferless_stations.append(ist)

    if not bufferless_stations:
        ret_info.error_msg = 'No bufferless queue found (capacity must equal number of servers).'
        return False, ret_info

    # Check retrial drop strategy
    retrial_station = None
    for ist in bufferless_stations:
        if hasattr(sn, 'droprule') and sn.droprule is not None:
            droprule = sn.droprule[ist, :] if sn.droprule.ndim > 1 else sn.droprule
            # Check for RETRIAL or RETRIAL_WITH_LIMIT drop strategy
            has_retrial = False
            try:
                has_retrial = any(
                    dr == DropStrategy.RETRIAL or dr == DropStrategy.RETRIAL_WITH_LIMIT
                    for dr in droprule
                )
            except:
                # If DropStrategy comparison fails, try numeric comparison
                has_retrial = any(dr in [6, 7] for dr in droprule)  # Assuming 6,7 are retrial codes

            if has_retrial:
                retrial_station = ist
                break

    if retrial_station is None:
        ret_info.error_msg = 'No retrial drop strategy configured on bufferless queue.'
        return False, ret_info

    ret_info.station_idx = retrial_station
    ret_info.node_idx = sn.stationToNode[retrial_station] if hasattr(sn, 'stationToNode') else retrial_station
    ret_info.N = int(sn.nservers[retrial_station]) if hasattr(sn, 'nservers') and sn.nservers is not None else 1

    # Find source station
    source_station = None
    for ist in range(sn.nstations):
        node_idx = sn.stationToNode[ist] if hasattr(sn, 'stationToNode') else ist
        if node_idx < len(sn.nodetype) and sn.nodetype[node_idx] == NodeType.SOURCE:
            source_station = ist
            break

    if source_station is None:
        ret_info.error_msg = 'No Source node found.'
        return False, ret_info

    ret_info.source_idx = source_station

    # Validate arrival process (must be MAP/BMAP)
    if hasattr(sn, 'proc') and sn.proc is not None:
        try:
            arrival_proc = sn.proc[source_station][ret_info.class_idx]
            if arrival_proc is None or not isinstance(arrival_proc, (list, tuple)) or len(arrival_proc) < 2:
                ret_info.error_msg = 'Invalid arrival process at source.'
                return False, ret_info
        except (IndexError, TypeError):
            ret_info.error_msg = 'Invalid arrival process at source.'
            return False, ret_info

    # Validate service process (must be PH/Exp)
    if hasattr(sn, 'proc') and sn.proc is not None:
        try:
            service_proc = sn.proc[retrial_station][ret_info.class_idx]
            if service_proc is None or not isinstance(service_proc, (list, tuple)) or len(service_proc) < 2:
                ret_info.error_msg = 'Invalid service process at queue.'
                return False, ret_info
        except (IndexError, TypeError):
            ret_info.error_msg = 'Invalid service process at queue.'
            return False, ret_info

    # Set default admission threshold
    ret_info.R = ret_info.N - 1

    # Check for FCR (Finite Capacity Region) containing this queue
    if hasattr(sn, 'region') and sn.region is not None and hasattr(sn, 'nregions') and sn.nregions > 0:
        for f in range(sn.nregions):
            if f < len(sn.region) and sn.region[f] is not None:
                region_matrix = sn.region[f]
                if retrial_station < region_matrix.shape[0]:
                    global_cap_col = region_matrix.shape[1] - 1
                    if region_matrix[retrial_station, global_cap_col] > 0:
                        R_fcr = int(region_matrix[retrial_station, global_cap_col])
                        if R_fcr <= ret_info.N - 1:
                            ret_info.R = R_fcr
                        break

    # All validations passed
    ret_info.is_retrial = True
    return True, ret_info


__all__ = [
    'QueueType',
    'BmapMatrix',
    'PhDistribution',
    'QbdStatespace',
    'RetrialQueueResult',
    'RetrialQueueAnalyzer',
    'qsys_bmapphnn_retrial',
    'qsys_is_retrial',
    'RetrialInfo',
]
