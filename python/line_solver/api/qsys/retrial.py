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

        LINE stores arrival processes in MAP format: {D0, D1, D2, ...}
        where D0 is the "hidden" generator and D1, D2, ... are arrival matrices.

        Returns:
            BmapMatrix instance or None if arrival is not BMAP/MAP
        """
        sn = self.sn

        # Find source station
        source_idx = None
        for i in range(sn.nstations):
            node_idx = int(sn.stationToNode[i])
            # Check if this is a source node
            if hasattr(sn, 'nodetype') and sn.nodetype[node_idx] == 0:  # Source type
                source_idx = i
                break

        if source_idx is None:
            return None

        # Extract process from sn.proc
        # Assume single class for now (classIdx = 0)
        class_idx = 0
        if not hasattr(sn, 'proc') or sn.proc is None:
            return None

        try:
            proc = sn.proc[source_idx][class_idx]
        except (IndexError, KeyError, TypeError):
            return None

        if proc is None or not isinstance(proc, (list, tuple)) or len(proc) < 2:
            return None

        D0 = np.atleast_2d(proc[0])
        D_batch = []

        # Check if D0 looks like a subgenerator (negative diagonal)
        if D0.shape[0] == D0.shape[1]:
            diag_D0 = np.diag(D0)
            if np.all(diag_D0 <= 0):
                # This is in MAP format {D0, D1, ...}
                for i in range(1, len(proc)):
                    D_batch.append(np.atleast_2d(proc[i]))
                return BmapMatrix(D0=D0, D_batch=D_batch)

        # Try to interpret as PH format {alpha, T} and convert to MAP
        alpha = np.atleast_1d(proc[0]).flatten()
        T = np.atleast_2d(proc[1])

        if len(alpha) == T.shape[0] == T.shape[1]:
            n = len(alpha)
            e = np.ones(n)
            # Convert PH to MAP: D0 = T, D1 = (-T*e)*alpha
            D0_new = T
            D1_new = np.outer(-T @ e, alpha)
            return BmapMatrix(D0=D0_new, D_batch=[D1_new])

        return None

    def extract_ph_service(self) -> Optional[PhDistribution]:
        """
        Extract Phase-Type service parameters.

        LINE stores PH in MAP format: {D0, D1}
        D0 = T (subgenerator matrix)
        D1 = S0 * alpha (exit rate times initial prob)

        Returns:
            PhDistribution instance or None if service is not PH
        """
        sn = self.sn

        # Find queue station
        queue_idx = None
        for i in range(sn.nstations):
            node_idx = int(sn.stationToNode[i])
            # Check if this is a queue node
            if hasattr(sn, 'nodetype') and sn.nodetype[node_idx] == 1:  # Queue type
                queue_idx = i
                break

        if queue_idx is None:
            return None

        # Extract process from sn.proc
        class_idx = 0
        if not hasattr(sn, 'proc') or sn.proc is None:
            return None

        try:
            proc = sn.proc[queue_idx][class_idx]
        except (IndexError, KeyError, TypeError):
            return None

        if proc is None or not isinstance(proc, (list, tuple)) or len(proc) < 2:
            return None

        D0 = np.atleast_2d(proc[0])
        D1 = np.atleast_2d(proc[1])

        # Validate dimensions
        n = D0.shape[0]
        if D0.shape[1] != n or D1.shape[0] != n or D1.shape[1] != n:
            return None

        # S = D0 (subgenerator)
        S = D0

        # S0 = -S * ones (exit rates)
        e = np.ones(n)
        S0 = -S @ e

        # Extract alpha from D1 = S0 * alpha
        # Find a row with non-zero exit rate
        idx = np.where(S0 > 1e-10)[0]
        if len(idx) == 0:
            # All rows have zero exit rate - use uniform
            beta = np.ones(n) / n
        else:
            # alpha = D1(idx,:) / S0(idx)
            beta = D1[idx[0], :] / S0[idx[0]]

        # Ensure beta is properly normalized
        beta = beta.flatten()
        if abs(np.sum(beta) - 1) > 1e-6:
            if np.sum(beta) > 0:
                beta = beta / np.sum(beta)
            else:
                beta = np.ones(n) / n

        return PhDistribution(beta=beta, S=S)

    def extract_retrial_parameters(self) -> Optional[Dict[str, float]]:
        """
        Extract retrial-specific parameters.

        Returns:
            Dict with keys:
            - alpha: Retrial rate (rate at which customers retry)
            - gamma: Orbit impatience rate (reneging rate)
            - p: Batch rejection probability
            - R: Threshold for admission control
            - N: Number of servers
        """
        sn = self.sn
        params = {
            'alpha': 0.1,  # Default retrial rate
            'gamma': 0.0,  # Default orbit impatience (no abandonment)
            'p': 0.0,      # Default batch rejection probability
            'R': 0,        # Default admission threshold
            'N': 1         # Default number of servers
        }

        # Find queue station
        queue_idx = None
        for i in range(sn.nstations):
            node_idx = int(sn.stationToNode[i])
            if hasattr(sn, 'nodetype') and sn.nodetype[node_idx] == 1:
                queue_idx = i
                break

        if queue_idx is None:
            return None

        class_idx = 0

        # Extract number of servers
        if hasattr(sn, 'nservers'):
            try:
                params['N'] = int(sn.nservers[queue_idx])
            except (IndexError, TypeError):
                pass

        # Extract retrial rate from retrialDelays if available
        if hasattr(sn, 'retrialDelays') and sn.retrialDelays is not None:
            try:
                retrial_dist = sn.retrialDelays[queue_idx][class_idx]
                if retrial_dist is not None and isinstance(retrial_dist, (list, tuple)):
                    # For Exp(alpha), the rate matrix is -alpha
                    params['alpha'] = -float(retrial_dist[1][0, 0])
            except (IndexError, TypeError, KeyError):
                pass

        # Extract orbit impatience rate
        if hasattr(sn, 'orbitImpatience') and sn.orbitImpatience is not None:
            try:
                impatience_dist = sn.orbitImpatience[queue_idx][class_idx]
                if impatience_dist is not None and isinstance(impatience_dist, (list, tuple)):
                    params['gamma'] = -float(impatience_dist[1][0, 0])
            except (IndexError, TypeError, KeyError):
                pass

        # Extract batch rejection probability
        if hasattr(sn, 'batchRejectProb') and sn.batchRejectProb is not None:
            try:
                params['p'] = float(sn.batchRejectProb[queue_idx, class_idx])
            except (IndexError, TypeError):
                pass

        # Set admission threshold R = N - 1 by default
        params['R'] = params['N'] - 1

        return params

    def build_qbd_statespace(self,
                            bmap: BmapMatrix,
                            ph_service: PhDistribution,
                            retrial_params: Dict[str, float]
                            ) -> Optional[QbdStatespace]:
        """
        Build QBD state space for the retrial queue.

        This constructs the generator matrix blocks for the QBD process
        following the BMAP/PH/N/N retrial queue formulation.

        Args:
            bmap: BMAP arrival process
            ph_service: PH service distribution
            retrial_params: Retrial parameters (alpha, gamma, p, R, N)

        Returns:
            QbdStatespace instance or None if construction fails
        """
        if bmap is None or ph_service is None or retrial_params is None:
            return None

        # Extract parameters
        D0 = bmap.D0
        D1 = bmap.D_batch[0] if bmap.D_batch else np.zeros_like(D0)
        beta = ph_service.beta
        S = ph_service.S
        alpha = retrial_params.get('alpha', 0.1)
        gamma = retrial_params.get('gamma', 0.0)
        N = retrial_params.get('N', 1)

        # Dimensions
        m_a = D0.shape[0]  # BMAP phases
        m_s = S.shape[0]   # Service phases
        phase_dim = m_a * m_s * N  # Combined phase dimension

        # Exit rate vector for service
        e_s = np.ones(m_s)
        S0 = -S @ e_s

        # Identity matrices
        I_a = np.eye(m_a)
        I_s = np.eye(m_s)

        # Build QBD generator blocks
        # A0: transitions that increase level (arrivals to orbit)
        # A1: transitions at same level
        # A2: transitions that decrease level (retrials from orbit)

        # Simplified construction for single-server case
        if N == 1:
            # Phase = (arrival phase, service phase, server state)
            # Server state: 0 = idle, 1 = busy

            # Build A0 (arrivals when server busy - go to orbit)
            A0 = np.kron(D1, I_s)

            # Build A2 (retrials from orbit when server becomes free)
            A2 = alpha * np.kron(I_a, np.outer(S0, beta))

            # Build A1 (internal transitions + arrivals when idle)
            A1 = np.kron(D0, I_s) + np.kron(I_a, S)

            # Add impatience (customers leaving orbit)
            if gamma > 0:
                A1 = A1 + gamma * np.eye(phase_dim)

            # Boundary block B0
            B0 = np.kron(D0 + D1, I_s) + np.kron(I_a, S)

        else:
            # Multi-server case - more complex construction
            # Placeholder: return None for now
            return None

        return QbdStatespace(
            A0=A0,
            A1=A1,
            A2=A2,
            B0=B0,
            max_level=self.max_retrial_orbit,
            phase_dim=phase_dim
        )

    def analyze(self) -> RetrialQueueResult:
        """
        Analyze the retrial queue.

        This is the main entry point for analysis:
        1. Detect queue type
        2. Extract arrival and service parameters
        3. Build QBD state space
        4. Solve for stationary distribution using matrix-analytic methods
        5. Compute performance metrics

        Returns:
            RetrialQueueResult with performance metrics
        """
        queue_type = self.detect_queue_type()

        if queue_type == QueueType.STANDARD:
            raise ValueError("Standard queue - use standard analyzer")

        # Extract parameters
        bmap = self.extract_bmap()
        ph_service = self.extract_ph_service()
        retrial_params = self.extract_retrial_parameters()

        if bmap is None or ph_service is None or retrial_params is None:
            return RetrialQueueResult(
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

        # Build QBD state space
        qbd = self.build_qbd_statespace(bmap, ph_service, retrial_params)

        if qbd is None:
            return RetrialQueueResult(
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

        # Solve QBD using matrix-geometric method
        try:
            pi, R, converged, iterations, error = self._solve_qbd_matrix_geometric(qbd)
        except Exception:
            return RetrialQueueResult(
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

        # Compute performance metrics
        N = retrial_params.get('N', 1)
        arrival_rate = bmap.arrival_rate
        service_rate = ph_service.rate if hasattr(ph_service, 'rate') else 1.0 / ph_service.mean

        # Expected number in orbit (sum over levels weighted by level)
        L_orbit = 0.0
        for level in range(len(pi)):
            L_orbit += level * np.sum(pi[level])

        # Server utilization
        P_idle = np.sum(pi[0]) if len(pi) > 0 else 1.0
        utilization = 1.0 - P_idle

        # Throughput via departure rate
        throughput = utilization * N * service_rate

        # Expected number being served
        N_server = utilization * N

        # Probability orbit is empty
        P_empty_orbit = np.sum(pi[0]) if len(pi) > 0 else 1.0

        return RetrialQueueResult(
            queue_type=queue_type,
            L_orbit=L_orbit,
            N_server=N_server,
            utilization=utilization,
            throughput=throughput,
            P_idle=P_idle,
            P_empty_orbit=P_empty_orbit,
            stationary_dist=pi if converged else None,
            truncation_level=len(pi) - 1 if pi is not None else 0,
            converged=converged,
            iterations=iterations,
            error=error
        )

    def _solve_qbd_matrix_geometric(self, qbd: QbdStatespace) -> Tuple[np.ndarray, np.ndarray, bool, int, float]:
        """
        Solve QBD process using matrix-geometric method.

        The rate matrix R satisfies: A0 + R*A1 + R^2*A2 = 0

        Args:
            qbd: QBD state space with generator blocks

        Returns:
            Tuple of (stationary distribution, R matrix, converged, iterations, error)
        """
        A0 = qbd.A0
        A1 = qbd.A1
        A2 = qbd.A2
        n = A0.shape[0]

        # Initialize R matrix
        R = np.zeros((n, n))

        # Iteration for R: R = -A0 * (A1 + R*A2)^{-1}
        converged = False
        iterations = 0
        error = np.inf

        for it in range(self.max_iterations):
            iterations = it + 1

            try:
                # Compute (A1 + R*A2)^{-1}
                inv_term = np.linalg.inv(A1 + R @ A2)
                R_new = -A0 @ inv_term
            except np.linalg.LinAlgError:
                break

            # Check convergence
            error = np.max(np.abs(R_new - R))
            R = R_new

            if error < self.tolerance:
                converged = True
                break

        if not converged:
            return None, R, False, iterations, error

        # Compute boundary probabilities
        # pi_0 satisfies: pi_0 * (B0 + R*A2) = 0, sum(pi_0) * (I - R)^{-1} * e = 1
        try:
            B = qbd.B0 + R @ A2
            # Find null space
            eigvals, eigvecs = np.linalg.eig(B.T)
            null_idx = np.argmin(np.abs(eigvals))
            pi_0 = np.real(eigvecs[:, null_idx])
            pi_0 = np.abs(pi_0)

            # Normalize
            I_minus_R_inv = np.linalg.inv(np.eye(n) - R)
            e = np.ones(n)
            norm_factor = pi_0 @ I_minus_R_inv @ e
            if norm_factor > 0:
                pi_0 = pi_0 / norm_factor
        except np.linalg.LinAlgError:
            return None, R, False, iterations, error

        # Build stationary distribution up to truncation level
        pi = [pi_0]
        current = pi_0
        for level in range(1, qbd.max_level + 1):
            current = current @ R
            if np.sum(current) < self.tolerance:
                break
            pi.append(current.copy())

        return np.array(pi), R, converged, iterations, error


def qsys_bmapphnn_retrial(
    arrival_matrix: Dict[str, np.ndarray],
    service_params: Dict[str, np.ndarray],
    N: int,
    retrial_params: Optional[Dict[str, float]] = None,
    options: Optional[Dict] = None
) -> RetrialQueueResult:
    """
    Analyze BMAP/PH/N/N bufferless retrial queue.

    Implements the algorithm from Dudin et al., "Analysis of BMAP/PH/N-Type
    Queueing System with Flexible Retrials Admission Control",
    Mathematics 2025, 13(9), 1434.

    Parameters:
        arrival_matrix: Dict with 'D0', 'D1', ... for BMAP matrices.
            D0: hidden transition matrix (V x V).
            D1, ..., DK: arrival matrices for batch sizes 1, ..., K.
        service_params: Dict with 'beta' (initial prob vector, 1xM) and
            'S' (PH subgenerator matrix, MxM).
        N: Number of servers (also capacity, hence bufferless).
        retrial_params: Dict with 'alpha' (retrial rate per customer),
            'gamma' (impatience/abandonment rate), 'p' (batch rejection
            probability), 'R' (admission threshold, scalar or 1xV).
        options: Dict with optional keys:
            'MaxLevel': max orbit level for truncation (default: auto).
            'Tolerance': convergence tolerance (default: 1e-10).
            'Verbose': print progress (default: False).

    Returns:
        RetrialQueueResult with performance metrics.

    References:
        Dudin, A., Klimenok, V., & Vishnevsky, V. (2020).
        Port from: matlab/src/api/qsys/qsys_bmapphnn_retrial.m
    """
    if options is None:
        options = {}

    tol = options.get('Tolerance', 1e-10)
    verbose = options.get('Verbose', False)
    max_level_param = options.get('MaxLevel', None)

    # --- Extract BMAP matrices D = [D0, D1, ..., DK] ---
    D0 = np.atleast_2d(arrival_matrix.get('D0'))
    D_list = [D0]
    for i in range(1, 100):
        Dk = arrival_matrix.get(f'D{i}')
        if Dk is None:
            break
        D_list.append(np.atleast_2d(Dk))

    K = len(D_list) - 1  # max batch size
    V = D0.shape[0]      # number of BMAP states

    # Stationary distribution of fundamental process D^(1) = sum(D_k)
    D1_gen = sum(D_list)
    theta = _compute_stationary_vector(D1_gen)

    # Mean arrival rate: lambda = theta * sum(k * D_k) * e
    sum_k_Dk = np.zeros_like(D0)
    for k in range(1, K + 1):
        sum_k_Dk += k * D_list[k]
    lam = float(theta @ sum_k_Dk @ np.ones(V))

    # PH service parameters
    beta = np.atleast_1d(service_params.get('beta', np.array([1.0]))).astype(float).flatten()
    S = np.atleast_2d(service_params.get('S', np.array([[-1.0]]))).astype(float)
    M_ph = S.shape[0]
    S0 = -S @ np.ones(M_ph)
    b1 = float(beta @ np.linalg.solve(-S, np.ones(M_ph)))  # mean service time

    # Retrial parameters
    if retrial_params is None:
        retrial_params = {}
    alpha = float(retrial_params.get('alpha', 1.0))
    gamma = float(retrial_params.get('gamma', 0.0))
    p = float(retrial_params.get('p', 0.0))
    R_param = retrial_params.get('R', N)
    if np.isscalar(R_param):
        R_vec = np.full(V, float(R_param))
    else:
        R_vec = np.asarray(R_param, dtype=float).flatten()

    # T_n = C(n+M-1, M-1) = number of service states with n busy servers
    from math import comb
    T = np.array([comb(n + M_ph - 1, M_ph - 1) for n in range(N + 1)], dtype=int)
    d = int(np.sum(T))  # total dimension per BMAP state

    # Build state mapping (weak compositions)
    state_map = [_generate_compositions(n, M_ph) for n in range(N + 1)]

    # Truncation level
    rho = lam * b1 / N if N > 0 else 0.0
    if max_level_param is None:
        trunc_level = max(100, int(np.ceil(50.0 / (1.0 - min(rho, 0.99)))))
    else:
        trunc_level = int(max_level_param)

    if verbose:
        print(f"Solving BMAP/PH/N/N retrial queue...")
        print(f"  V={V}, M={M_ph}, N={N}, K={K}")
        print(f"  d={d}, block size Vd={V * d}")
        print(f"  lambda={lam:.4f}, mu={1.0/b1 if b1 > 0 else float('inf'):.4f}")
        print(f"  Offered load rho={rho:.4f}")
        print(f"  Truncation level: {trunc_level}")

    Vd = V * d
    total_dim = (trunc_level + 1) * Vd

    if verbose:
        print(f"Total matrix dimension: {total_dim} x {total_dim}")

    # Context object for helper functions
    ctx = _RetrialCtx(D_list, beta, S, S0, M_ph, N, V, K, d, T, R_vec,
                       alpha, gamma, p, state_map)

    # Build generator matrix
    use_sparse = total_dim > 5000
    if use_sparse:
        from scipy.sparse import lil_matrix
        Q = lil_matrix((total_dim, total_dim))
    else:
        Q = np.zeros((total_dim, total_dim))

    for i in range(trunc_level + 1):
        if verbose and i % 20 == 0:
            print(f"  Level {i}/{trunc_level}")
        row_s = i * Vd
        row_e = (i + 1) * Vd
        for j in range(max(0, i - 1), min(trunc_level, i + K) + 1):
            col_s = j * Vd
            col_e = (j + 1) * Vd
            Qij = _build_generator_level(ctx, i, j)
            if use_sparse:
                Q[row_s:row_e, col_s:col_e] = Qij
            else:
                Q[row_s:row_e, col_s:col_e] = Qij

    # Ensure rows sum to zero (set diagonal)
    if use_sparse:
        Q = Q.tocsr()
        Q_dense = Q.toarray()
    else:
        Q_dense = Q
    for i in range(total_dim):
        Q_dense[i, i] -= np.sum(Q_dense[i, :])

    if verbose:
        print("Solving linear system...")

    # Solve pi * Q = 0, pi * e = 1
    # Replace last column with ones for normalization
    Q_dense[:, -1] = 1.0
    b_vec = np.zeros(total_dim)
    b_vec[-1] = 1.0
    pi_flat = np.linalg.solve(Q_dense.T, b_vec)

    # Reshape to level structure: (trunc_level+1, Vd)
    pi = pi_flat.reshape(trunc_level + 1, Vd)

    # Handle numerical issues
    if np.any(pi < -1e-8):
        pi[pi < 0] = 0.0
    pi /= np.sum(pi)

    # --- Compute performance measures ---
    max_lvl = pi.shape[0] - 1

    # Mean number in orbit
    L_orbit = sum(i * np.sum(pi[i, :]) for i in range(1, max_lvl + 1))

    # Mean number of busy servers
    N_server = 0.0
    for i in range(max_lvl + 1):
        pi_level = pi[i, :]
        for nu in range(V):
            for n in range(N + 1):
                offset = nu * d + _get_block_offset(T, n)
                for t in range(T[n]):
                    idx = offset + t
                    if idx < Vd:
                        N_server += n * pi_level[idx]

    # Probability all servers idle
    P_idle = 0.0
    for i in range(max_lvl + 1):
        pi_level = pi[i, :]
        for nu in range(V):
            offset = nu * d  # n=0 block starts at offset
            P_idle += pi_level[offset]

    # Probability orbit empty
    P_empty_orbit = float(np.sum(pi[0, :]))

    # Probability system empty
    P_empty = 0.0
    pi_level0 = pi[0, :]
    for nu in range(V):
        offset = nu * d
        P_empty += pi_level0[offset]

    result = RetrialQueueResult(
        queue_type=QueueType.RETRIAL,
        L_orbit=L_orbit,
        N_server=N_server,
        utilization=N_server / N if N > 0 else 0.0,
        throughput=N_server / b1 if b1 > 0 else 0.0,
        P_idle=P_idle,
        P_empty_orbit=P_empty_orbit,
        truncation_level=trunc_level,
        converged=True,
        error=0.0
    )

    if verbose:
        print("Solution complete.")

    return result


# ========== Retrial solver helpers ==========

class _RetrialCtx:
    """Context for retrial solver helper functions (matches MATLAB ctx struct)."""
    __slots__ = ('D', 'beta', 'S', 'S0', 'M', 'N', 'V', 'K', 'd', 'T',
                 'R', 'alpha', 'gamma', 'p', 'state_map')

    def __init__(self, D, beta, S, S0, M, N, V, K, d, T, R, alpha, gamma, p, state_map):
        self.D = D
        self.beta = beta
        self.S = S
        self.S0 = S0
        self.M = M
        self.N = N
        self.V = V
        self.K = K
        self.d = d
        self.T = T
        self.R = R
        self.alpha = alpha
        self.gamma = gamma
        self.p = p
        self.state_map = state_map


def _compute_stationary_vector(Q: np.ndarray) -> np.ndarray:
    """Solve theta * Q = 0, theta * e = 1."""
    n = Q.shape[0]
    A = Q.T.copy()
    A[-1, :] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0
    return np.linalg.solve(A, b)


def _generate_compositions(n: int, M: int) -> np.ndarray:
    """Generate all weak compositions of n into M parts (reverse lexicographic)."""
    if M == 1:
        return np.array([[n]])
    from math import comb
    num_comps = comb(n + M - 1, M - 1)
    comps = np.zeros((num_comps, M), dtype=int)
    idx = 0
    for m1 in range(n, -1, -1):
        sub = _generate_compositions(n - m1, M - 1)
        num_sub = sub.shape[0]
        comps[idx:idx + num_sub, 0] = m1
        comps[idx:idx + num_sub, 1:] = sub
        idx += num_sub
    return comps


def _get_block_offset(T: np.ndarray, n: int) -> int:
    """Get starting index (0-based) for states with n busy servers."""
    return int(np.sum(T[:n])) if n > 0 else 0


def _compute_L(ctx: _RetrialCtx, n: int) -> Optional[np.ndarray]:
    """Matrix L_n: service completion transitions (T_n x T_{n-1})."""
    if n == 0:
        return None
    L = np.zeros((ctx.T[n], ctx.T[n - 1]))
    comps_n = ctx.state_map[n]
    comps_nm1 = ctx.state_map[n - 1]
    for i in range(comps_n.shape[0]):
        m = comps_n[i, :]
        for l in range(ctx.M):
            if m[l] > 0:
                m_prime = m.copy()
                m_prime[l] -= 1
                for j in range(comps_nm1.shape[0]):
                    if np.array_equal(comps_nm1[j, :], m_prime):
                        L[i, j] += m[l] * ctx.S0[l]
                        break
    return L


def _compute_A(ctx: _RetrialCtx, n: int) -> np.ndarray:
    """Matrix A_n: phase change transitions (T_n x T_n)."""
    if n == 0:
        return np.zeros((1, 1))
    A = np.zeros((ctx.T[n], ctx.T[n]))
    comps = ctx.state_map[n]
    for i in range(comps.shape[0]):
        m = comps[i, :]
        for l in range(ctx.M):
            if m[l] > 0:
                for lp in range(ctx.M):
                    if lp != l and ctx.S[l, lp] > 0:
                        m_prime = m.copy()
                        m_prime[l] -= 1
                        m_prime[lp] += 1
                        for j in range(comps.shape[0]):
                            if np.array_equal(comps[j, :], m_prime):
                                A[i, j] += m[l] * ctx.S[l, lp]
                                break
    return A


def _compute_P(ctx: _RetrialCtx, n: int) -> Optional[np.ndarray]:
    """Matrix P_n: new arrival transitions (T_n x T_{n+1})."""
    if n >= ctx.N:
        return None
    P = np.zeros((ctx.T[n], ctx.T[n + 1]))
    comps_n = ctx.state_map[n]
    comps_np1 = ctx.state_map[n + 1]
    for i in range(comps_n.shape[0]):
        m = comps_n[i, :]
        for l in range(ctx.M):
            if ctx.beta[l] > 0:
                m_prime = m.copy()
                m_prime[l] += 1
                for j in range(comps_np1.shape[0]):
                    if np.array_equal(comps_np1[j, :], m_prime):
                        P[i, j] += ctx.beta[l]
                        break
    return P


def _compute_Delta(ctx: _RetrialCtx, n: int) -> np.ndarray:
    """Diagonal matrix Delta_n: exit rates (T_n x T_n)."""
    if n == 0:
        return np.zeros((1, 1))
    comps = ctx.state_map[n]
    diag_vals = np.zeros(ctx.T[n])
    for i in range(comps.shape[0]):
        m = comps[i, :]
        diag_vals[i] = sum(m[l] * (-ctx.S[l, l]) for l in range(ctx.M))
    return np.diag(diag_vals)


def _compute_Gamma(ctx: _RetrialCtx, nu: int) -> np.ndarray:
    """Diagonal matrix Gamma^(nu): 0 for n<=R_nu, 1 for n>R_nu."""
    diag_vals = np.zeros(ctx.d)
    offset = 0
    for n in range(ctx.N + 1):
        if n > ctx.R[nu]:
            diag_vals[offset:offset + ctx.T[n]] = 1.0
        offset += ctx.T[n]
    return np.diag(diag_vals)


def _compute_G_nn(ctx: _RetrialCtx, n: int, nu: int, nu_prime: int) -> np.ndarray:
    """G_{n,n}^{(nu,nu')} matrix for batch losses."""
    if n <= ctx.N - ctx.K:
        return np.zeros((ctx.T[n], ctx.T[n]))
    total = 0.0
    for k in range(ctx.N - n + 1, ctx.K + 1):
        if 1 <= k <= ctx.K:
            total += ctx.D[k][nu, nu_prime]
    return ctx.p * total * np.eye(ctx.T[n])


def _compute_B(ctx: _RetrialCtx, nu: int) -> np.ndarray:
    """Block matrix B^(nu) of size d x d."""
    B = np.zeros((ctx.d, ctx.d))

    L_mats = [_compute_L(ctx, n) for n in range(ctx.N + 1)]
    A_mats = [_compute_A(ctx, n) for n in range(ctx.N + 1)]
    P_mats = [_compute_P(ctx, n) for n in range(ctx.N + 1)]
    Delta_mats = [_compute_Delta(ctx, n) for n in range(ctx.N + 1)]

    for n in range(ctx.N + 1):
        rs = _get_block_offset(ctx.T, n)
        re = rs + ctx.T[n]

        # Diagonal block
        G_nn = _compute_G_nn(ctx, n, nu, nu)
        if n == 0:
            B[rs, rs] = G_nn[0, 0] if G_nn.size > 0 else 0.0
        else:
            B[rs:re, rs:re] = A_mats[n] + Delta_mats[n] + G_nn

        # Subdiagonal block
        if n >= 1 and L_mats[n] is not None:
            cs = _get_block_offset(ctx.T, n - 1)
            ce = cs + ctx.T[n - 1]
            B[rs:re, cs:ce] = L_mats[n]

        # Superdiagonal blocks
        for k in range(1, ctx.K + 1):
            if n + k <= ctx.N:
                cs = _get_block_offset(ctx.T, n + k)
                ce = cs + ctx.T[n + k]
                D_k_nu_nu = ctx.D[k][nu, nu]
                Pprod = np.eye(ctx.T[n])
                for j in range(n, n + k):
                    if j < ctx.N and P_mats[j] is not None:
                        Pprod = Pprod @ P_mats[j]
                B[rs:re, cs:ce] = D_k_nu_nu * Pprod

    return B


def _compute_Bbar(ctx: _RetrialCtx, nu: int) -> np.ndarray:
    """B_bar^(nu) matrix for successful retrials."""
    Bbar = np.zeros((ctx.d, ctx.d))
    for n in range(min(int(ctx.R[nu]), ctx.N - 1) + 1):
        rs = _get_block_offset(ctx.T, n)
        re = rs + ctx.T[n]
        cs = _get_block_offset(ctx.T, n + 1)
        ce = cs + ctx.T[n + 1]
        P_n = _compute_P(ctx, n)
        if P_n is not None:
            Bbar[rs:re, cs:ce] = P_n
    return Bbar


def _compute_Btilde(ctx: _RetrialCtx, nu: int, nu_prime: int) -> np.ndarray:
    """B_tilde^(nu, nu') for BMAP state transitions."""
    Bt = np.zeros((ctx.d, ctx.d))
    P_mats = [_compute_P(ctx, n) for n in range(ctx.N + 1)]

    for n in range(ctx.N + 1):
        rs = _get_block_offset(ctx.T, n)
        re = rs + ctx.T[n]

        # Diagonal block
        G_nn = _compute_G_nn(ctx, n, nu, nu_prime)
        Bt[rs:re, rs:re] = G_nn

        # Superdiagonal blocks
        for k in range(1, ctx.K + 1):
            if n + k <= ctx.N:
                cs = _get_block_offset(ctx.T, n + k)
                ce = cs + ctx.T[n + k]
                D_k = ctx.D[k][nu, nu_prime]
                Pprod = np.eye(ctx.T[n])
                for j in range(n, n + k):
                    if j < ctx.N and P_mats[j] is not None:
                        Pprod = Pprod @ P_mats[j]
                Bt[rs:re, cs:ce] = D_k * Pprod

    return Bt


def _compute_C(ctx: _RetrialCtx, n: int, k: int, nu: int, nu_prime: int) -> Optional[np.ndarray]:
    """C_{n,k}^(nu, nu') for partial batch admission to orbit."""
    if n < ctx.N - ctx.K + k:
        return np.zeros((ctx.T[n], ctx.T[ctx.N]))
    elif n < ctx.N:
        batch_size = ctx.N - n + k
        if 1 <= batch_size <= ctx.K:
            D_batch = ctx.D[batch_size][nu, nu_prime]
            P_mats = [_compute_P(ctx, nn) for nn in range(ctx.N + 1)]
            Pprod = np.eye(ctx.T[n])
            for j in range(n, ctx.N):
                if P_mats[j] is not None:
                    Pprod = Pprod @ P_mats[j]
            return (1 - ctx.p) * D_batch * Pprod
        else:
            return np.zeros((ctx.T[n], ctx.T[ctx.N]))
    else:  # n == N
        if 1 <= k <= ctx.K:
            D_k = ctx.D[k][nu, nu_prime]
            return (1 - ctx.p) * D_k * np.eye(ctx.T[ctx.N])
        else:
            return np.zeros((ctx.T[ctx.N], ctx.T[ctx.N]))


def _build_generator_level(ctx: _RetrialCtx, i: int, j: int) -> np.ndarray:
    """Build generator block Q_{i,j}."""
    Vd = ctx.V * ctx.d
    Qij = np.zeros((Vd, Vd))

    if j < max(0, i - 1) or j > i + ctx.K:
        return Qij

    # Precompute per-BMAP-state matrices
    B_mats = [_compute_B(ctx, nu) for nu in range(ctx.V)]
    Bbar_mats = [_compute_Bbar(ctx, nu) for nu in range(ctx.V)]
    Gamma_mats = [_compute_Gamma(ctx, nu) for nu in range(ctx.V)]

    if i == j:  # Diagonal block
        for nu in range(ctx.V):
            rs = nu * ctx.d
            re = (nu + 1) * ctx.d
            for nu_prime in range(ctx.V):
                cs = nu_prime * ctx.d
                ce = (nu_prime + 1) * ctx.d
                if nu == nu_prime:
                    D0_nn = ctx.D[0][nu, nu]
                    block = (D0_nn * np.eye(ctx.d) + B_mats[nu]
                             - i * (ctx.gamma + ctx.alpha) * np.eye(ctx.d)
                             + i * ctx.alpha * Gamma_mats[nu])
                    Qij[rs:re, cs:ce] = block
                else:
                    Bt = _compute_Btilde(ctx, nu, nu_prime)
                    D0_nn_p = ctx.D[0][nu, nu_prime]
                    Qij[rs:re, cs:ce] = Bt + D0_nn_p * np.eye(ctx.d)

    elif j == i - 1 and i >= 1:  # Subdiagonal block
        for nu in range(ctx.V):
            rs = nu * ctx.d
            re = (nu + 1) * ctx.d
            cs = rs
            ce = re
            block = i * ctx.gamma * np.eye(ctx.d) + i * ctx.alpha * Bbar_mats[nu]
            Qij[rs:re, cs:ce] = block

    elif j > i and j <= i + ctx.K:  # Superdiagonal blocks
        k = j - i
        for nu in range(ctx.V):
            rs = nu * ctx.d
            re = (nu + 1) * ctx.d
            for nu_prime in range(ctx.V):
                cs = nu_prime * ctx.d
                ce = (nu_prime + 1) * ctx.d
                block = np.zeros((ctx.d, ctx.d))
                for n in range(ctx.N + 1):
                    C_nk = _compute_C(ctx, n, k, nu, nu_prime)
                    if C_nk is not None and np.any(C_nk != 0):
                        n_rs = _get_block_offset(ctx.T, n)
                        n_re = n_rs + ctx.T[n]
                        N_cs = _get_block_offset(ctx.T, ctx.N)
                        N_ce = N_cs + ctx.T[ctx.N]
                        if C_nk.shape[1] == ctx.T[ctx.N]:
                            block[n_rs:n_re, N_cs:N_ce] = C_nk
                Qij[rs:re, cs:ce] = block

    return Qij


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


@dataclass
class RenegingInfo:
    """Information about a valid reneging queue topology."""
    is_reneging: bool
    source_idx: Optional[int] = None
    queue_idx: Optional[int] = None
    class_idx: Optional[int] = None
    n_servers: Optional[int] = None
    service_rate: Optional[float] = None
    error_msg: str = ''


def detect_reneging_topology(sn: Any) -> Tuple[bool, RenegingInfo]:
    """
    Detect if model is suitable for MAP/M/s+G (MAPMsG) reneging solver.

    Requirements:
    - Open model, single class
    - Single queue station with reneging/patience configured
    - MAP/BMAP arrival at source
    - Exponential service at queue (single-phase PH)
    - FCFS scheduling

    Args:
        sn: NetworkStruct object

    Returns:
        Tuple of (is_reneging, reneging_info)

    References:
        Original MATLAB: matlab/src/solvers/MAM/solver_mam_retrial.m (detectRenegingTopology)
    """
    from ..sn import sn_is_open_model, NodeType
    from ...lang.base import ImpatienceType, SchedStrategy

    info = RenegingInfo(is_reneging=False)

    # Check open model
    try:
        if not sn_is_open_model(sn):
            info.error_msg = 'MAPMsG requires open queueing model.'
            return False, info
    except Exception:
        info.error_msg = 'Could not determine if model is open.'
        return False, info

    # Check single class (current limitation)
    if sn.nclasses > 1:
        info.error_msg = 'MAPMsG currently supports single class only.'
        return False, info
    info.class_idx = 0

    # Find source and queue stations
    source_idx = None
    queue_idx = None
    for ist in range(sn.nstations):
        node_idx = sn.stationToNode[ist] if hasattr(sn, 'stationToNode') else ist
        if node_idx < len(sn.nodetype):
            if sn.nodetype[node_idx] == NodeType.SOURCE:
                source_idx = ist
            elif sn.nodetype[node_idx] == NodeType.QUEUE:
                if queue_idx is not None:
                    info.error_msg = 'MAPMsG requires single queue station.'
                    return False, info
                queue_idx = ist

    if source_idx is None:
        info.error_msg = 'No Source node found.'
        return False, info
    if queue_idx is None:
        info.error_msg = 'No Queue node found.'
        return False, info

    info.source_idx = source_idx
    info.queue_idx = queue_idx

    # Check for reneging patience configuration
    if not hasattr(sn, 'impatienceClass') or sn.impatienceClass is None:
        info.error_msg = 'No patience/impatience configuration found.'
        return False, info

    try:
        imp_val = sn.impatienceClass[queue_idx, info.class_idx]
        if imp_val != ImpatienceType.RENEGING:
            info.error_msg = 'Queue does not have reneging configured.'
            return False, info
    except (IndexError, TypeError):
        info.error_msg = 'Queue does not have reneging configured.'
        return False, info

    # Check patience distribution exists
    if not hasattr(sn, 'patienceProc') or sn.patienceProc is None:
        info.error_msg = 'No patience distribution found.'
        return False, info
    try:
        pat_proc = sn.patienceProc[queue_idx, info.class_idx] if hasattr(sn.patienceProc, '__getitem__') else None
        if pat_proc is None or (isinstance(pat_proc, (list, tuple)) and len(pat_proc) == 0):
            info.error_msg = 'No patience distribution for this class.'
            return False, info
    except (IndexError, TypeError, KeyError):
        info.error_msg = 'No patience distribution for this class.'
        return False, info

    # Check FCFS scheduling
    sched_val = sn.sched.get(queue_idx, None) if isinstance(sn.sched, dict) else (
        sn.sched[queue_idx] if queue_idx < len(sn.sched) else None
    )
    if sched_val is not None:
        sched_id = sched_val.value if hasattr(sched_val, 'value') else int(sched_val)
        if sched_id != SchedStrategy.FCFS:
            info.error_msg = 'MAPMsG requires FCFS scheduling.'
            return False, info

    # Check exponential service (single-phase)
    try:
        service_proc = sn.proc[queue_idx][info.class_idx]
        if service_proc is None or not isinstance(service_proc, (list, tuple)) or len(service_proc) < 2:
            info.error_msg = 'Invalid service process.'
            return False, info
        D0 = np.atleast_2d(service_proc[0])
        if D0.shape[0] != 1:
            info.error_msg = 'MAPMsG requires exponential service (single-phase).'
            return False, info
        info.service_rate = -float(D0[0, 0])
    except (IndexError, TypeError):
        info.error_msg = 'Invalid service process.'
        return False, info

    info.n_servers = int(sn.nservers[queue_idx]) if hasattr(sn, 'nservers') and sn.nservers is not None else 1

    # Check MAP arrival process
    try:
        arrival_proc = sn.proc[source_idx][info.class_idx]
        if arrival_proc is None or not isinstance(arrival_proc, (list, tuple)) or len(arrival_proc) < 2:
            info.error_msg = 'Invalid arrival process.'
            return False, info
    except (IndexError, TypeError):
        info.error_msg = 'Invalid arrival process.'
        return False, info

    # All checks passed
    info.is_reneging = True
    return True, info


def has_reneging_patience(sn: Any) -> bool:
    """
    Check if model has reneging/patience configured on any queue station.

    Returns True if any queue station has ImpatienceType.RENEGING
    configured with a patience distribution.

    Args:
        sn: NetworkStruct object

    Returns:
        True if reneging patience is configured

    References:
        Original MATLAB: matlab/src/solvers/MAM/solver_mam_analyzer.m (hasRenegingPatience)
    """
    from ...lang.base import ImpatienceType

    if not hasattr(sn, 'impatienceClass') or sn.impatienceClass is None:
        return False
    if not hasattr(sn, 'patienceProc') or sn.patienceProc is None:
        return False

    for ist in range(sn.nstations):
        for r in range(sn.nclasses):
            try:
                if sn.impatienceClass[ist, r] == ImpatienceType.RENEGING:
                    pat_proc = sn.patienceProc[ist, r] if hasattr(sn.patienceProc, '__getitem__') else None
                    if pat_proc is not None:
                        return True
            except (IndexError, TypeError, KeyError):
                continue

    return False


def extract_bmap_matrices(proc: Any) -> Optional[List[np.ndarray]]:
    """
    Extract BMAP matrices {D0, D1, ...} from LINE process representation.

    LINE stores arrival processes in MAP format: {D0, D1, D2, ...}
    where D0 is the "hidden" generator and D1, D2, ... are arrival matrices.
    May also be in PH format {alpha, T} which is converted to MAP.

    Args:
        proc: Process representation from sn.proc[station][class]

    Returns:
        List of numpy arrays [D0, D1, ...] or None if extraction fails

    References:
        Original MATLAB: matlab/src/solvers/MAM/solver_mam_retrial.m (extractBMAPMatrices)
    """
    if proc is None or not isinstance(proc, (list, tuple)):
        return None

    # Check if already in BMAP format (cell of cell arrays)
    if isinstance(proc[0], (list, tuple)):
        return [np.atleast_2d(p) for p in proc]

    # For LINE, arrival processes are stored as {D0, D1, D2, ...}
    if len(proc) >= 2 and isinstance(proc[0], np.ndarray) and isinstance(proc[1], np.ndarray):
        D0 = np.atleast_2d(proc[0])
        D1 = np.atleast_2d(proc[1])
        n = D0.shape[0]

        # Check if D0 looks like a generator (negative diagonal)
        if n == 1:
            # Scalar case (exponential)
            if D0[0, 0] < 0 and D1[0, 0] > 0:
                # Already in MAP format
                return [np.atleast_2d(p) for p in proc]
        else:
            # Matrix case - check if D0 has negative diagonal
            diag_D0 = np.diag(D0)
            if np.all(diag_D0 < 0) or np.all(diag_D0 <= 0):
                # Looks like MAP format {D0, D1, ...}
                return [np.atleast_2d(p) for p in proc]

        # Try to interpret as PH format {alpha, T} and convert to MAP
        alpha = np.atleast_1d(proc[0]).flatten()
        T = np.atleast_2d(proc[1])

        if len(alpha) == T.shape[0] == T.shape[1]:
            n = len(alpha)
            e = np.ones(n)
            # Convert PH to MAP: D0 = T, D1 = (-T*e)*alpha
            D0_new = T
            D1_new = np.outer(-T @ e, alpha)
            return [D0_new, D1_new]

    # Fallback: assume it's already in MAP format
    return [np.atleast_2d(p) for p in proc] if len(proc) >= 2 else None


def extract_ph_params(proc: Any) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Extract PH parameters (beta, S) from LINE process representation.

    LINE stores PH in MAP format: {D0, D1} where D0=T (subgenerator),
    D1=S0*alpha (exit rate times initial prob).

    Args:
        proc: Process representation from sn.proc[station][class]

    Returns:
        Tuple of (beta, S) where beta is initial probability vector and S is
        subgenerator matrix. Returns (None, None) if extraction fails.

    References:
        Original MATLAB: matlab/src/solvers/MAM/solver_mam_retrial.m (extractPHParams)
    """
    if proc is None or not isinstance(proc, (list, tuple)) or len(proc) < 2:
        return None, None

    if not isinstance(proc[0], np.ndarray) or not isinstance(proc[1], np.ndarray):
        return None, None

    D0 = np.atleast_2d(proc[0])
    D1 = np.atleast_2d(proc[1])

    # Validate dimensions
    n = D0.shape[0]
    if D0.shape[1] != n or D1.shape[0] != n or D1.shape[1] != n:
        return None, None

    # S = D0 (subgenerator)
    S = D0

    # S0 = -S * ones (exit rates)
    S0 = -S @ np.ones(n)

    # Extract alpha from D1 = S0 * alpha
    idx = np.where(S0 > 1e-10)[0]
    if len(idx) == 0:
        # All rows have zero exit rate - use uniform
        beta = np.ones(n) / n
    else:
        # alpha = D1(idx,:) / S0(idx)
        beta = D1[idx[0], :] / S0[idx[0]]

    # Ensure beta is row vector and sums to 1
    beta = beta.flatten()
    if abs(np.sum(beta) - 1) > 1e-6:
        if np.sum(beta) > 0:
            beta = beta / np.sum(beta)
        else:
            beta = np.ones(n) / n

    return beta, S


def convert_patience_to_regimes(
    patience_proc: Any,
    options: Optional[Dict] = None
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Convert patience distribution to piecewise-constant abandonment regimes for MAPMsG.

    Converts a patience distribution (in MAP/PH format) to boundary levels
    and abandonment function values for the MRMFQ solver.

    Args:
        patience_proc: Patience distribution from sn.patienceProc[station, class]
        options: Dict with optional 'mapmsg_quantization' key (default 11)

    Returns:
        Tuple of (boundary_levels, ga, quantization) where:
            boundary_levels: Array of regime boundary time points
            ga: Array of abandonment probabilities at each regime
            quantization: Number of regimes

    References:
        Original MATLAB: matlab/src/solvers/MAM/solver_mam_retrial.m (convertPatienceToRegimes)
    """
    if options is None:
        options = {}

    quantization = options.get('mapmsg_quantization', 11)

    # Extract patience rate from the distribution
    # For exponential patience Exp(gamma): patienceProc = {-gamma, gamma}
    gamma = 0.1  # default
    if isinstance(patience_proc, (list, tuple)) and len(patience_proc) >= 2:
        D0 = np.atleast_2d(patience_proc[0])
        gamma = -float(D0[0, 0])
    elif isinstance(patience_proc, dict):
        gamma = patience_proc.get('rate', 0.1)

    # Generate boundary levels and abandonment probabilities
    max_time = 10.0
    boundary_levels = np.linspace(0, max_time, quantization)

    # Compute abandonment probability at each boundary
    # ga(k) = F(BoundaryLevels(k)) for piecewise constant approximation
    ga = np.zeros(quantization + 1)
    ga[0] = 0.0  # No abandonment at time 0
    for k in range(quantization):
        midpoint = (boundary_levels[k] + boundary_levels[min(k + 1, quantization - 1)]) / 2
        if k == 0:
            midpoint = boundary_levels[k] / 2
        ga[k + 1] = 1.0 - np.exp(-gamma * midpoint)

    return boundary_levels, ga, quantization


def solver_mam_retrial(sn: Any, options: Optional[Dict] = None) -> Tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray,
    np.ndarray, np.ndarray, int
]:
    """
    Solve queueing models with customer impatience (retrial or reneging).

    Dispatches to either:
    1. RETRIAL: BMAP/PH/N/N bufferless retrial solver
    2. RENEGING: MAP/M/s+G solver (MAPMsG)

    Args:
        sn: NetworkStruct object
        options: Solver options dict with optional keys:
            'iter_max': Maximum truncation level (default 150)
            'tol': Convergence tolerance (default 1e-10)
            'verbose': Print progress messages (default False)
            'config': Dict with 'mapmsg_quantization' (default 11)

    Returns:
        Tuple of (QN, UN, RN, TN, CN, XN, totiter) where:
            QN: (M, K) queue lengths
            UN: (M, K) server utilizations
            RN: (M, K) response times
            TN: (M, K) throughputs
            CN: (1, K) cycle times
            XN: (1, K) system throughputs
            totiter: iteration/truncation level count

    References:
        Original MATLAB: matlab/src/solvers/MAM/solver_mam_retrial.m
    """
    if options is None:
        options = {}

    # Check for reneging (queue abandonment) first
    is_reneging, reneging_info = detect_reneging_topology(sn)
    if is_reneging:
        return _solve_reneging(sn, options, reneging_info)

    # Check for retrial topology
    is_retrial, ret_info = qsys_is_retrial(sn)
    if not is_retrial:
        raise ValueError(
            'No valid impatience configuration detected (retrial or reneging). '
            + ret_info.error_msg
        )

    # Initialize output arrays
    M = sn.nstations
    K = sn.nclasses
    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))
    CN = np.zeros((1, K))
    XN = np.zeros((1, K))

    # Extract indices
    source_idx = ret_info.source_idx
    queue_idx = ret_info.station_idx
    class_idx = ret_info.class_idx
    N = ret_info.N

    # Extract arrival process (BMAP) from source
    arrival_proc = sn.proc[source_idx][class_idx]
    D = extract_bmap_matrices(arrival_proc)

    if D is None:
        raise ValueError('Could not extract BMAP matrices from arrival process.')

    # Extract PH service distribution from queue
    service_proc = sn.proc[queue_idx][class_idx]
    beta, S = extract_ph_params(service_proc)

    if beta is None or S is None:
        raise ValueError('Could not extract PH parameters from service process.')

    # Extract retrial rate alpha from network configuration
    alpha = 0.1  # default
    if hasattr(sn, 'retrialDelays') and sn.retrialDelays is not None:
        try:
            retrial_dist = sn.retrialDelays[queue_idx, class_idx] if hasattr(sn.retrialDelays, '__getitem__') else None
            if retrial_dist is None:
                retrial_dist = sn.retrialDelays[queue_idx][class_idx]
            if retrial_dist is not None and isinstance(retrial_dist, (list, tuple)):
                alpha = -float(np.atleast_2d(retrial_dist[1])[0, 0])
        except (IndexError, TypeError, KeyError):
            pass

    # Extract orbit impatience gamma (default 0)
    gamma = ret_info.gamma
    if hasattr(sn, 'orbitImpatience') and sn.orbitImpatience is not None:
        try:
            imp_dist = sn.orbitImpatience[queue_idx, class_idx] if hasattr(sn.orbitImpatience, '__getitem__') else None
            if imp_dist is None:
                imp_dist = sn.orbitImpatience[queue_idx][class_idx]
            if imp_dist is not None and isinstance(imp_dist, (list, tuple)):
                gamma = -float(np.atleast_2d(imp_dist[1])[0, 0])
        except (IndexError, TypeError, KeyError):
            pass

    # Extract batch rejection probability p (default 0)
    p = ret_info.p
    if hasattr(sn, 'batchRejectProb') and sn.batchRejectProb is not None:
        try:
            p = float(sn.batchRejectProb[queue_idx, class_idx])
        except (IndexError, TypeError):
            pass

    # Extract admission threshold R (from FCR or default N-1)
    R = ret_info.R

    # Solve with options
    max_level = options.get('iter_max', 150)
    tol = options.get('tol', 1e-10)
    verbose = options.get('verbose', False)

    # Build arrival_matrix dict for qsys_bmapphnn_retrial
    arrival_matrix = {}
    for i, Di in enumerate(D):
        arrival_matrix[f'D{i}'] = Di

    service_params = {'beta': beta, 'S': S}
    retrial_params = {'alpha': alpha, 'gamma': gamma, 'p': p, 'R': R}
    solver_options = {'MaxLevel': max_level, 'Tolerance': tol, 'Verbose': verbose}

    perf = qsys_bmapphnn_retrial(arrival_matrix, service_params, N,
                                  retrial_params, solver_options)

    # Map to LINE output format
    # Queue length includes both orbit and servers
    QN[queue_idx, class_idx] = perf.L_orbit + perf.N_server
    UN[queue_idx, class_idx] = perf.utilization
    TN[queue_idx, class_idx] = perf.throughput

    # Response time via Little's law
    if perf.throughput > 0:
        RN[queue_idx, class_idx] = QN[queue_idx, class_idx] / perf.throughput
    else:
        RN[queue_idx, class_idx] = np.inf

    # System-level metrics
    XN[0, class_idx] = perf.throughput
    CN[0, class_idx] = RN[queue_idx, class_idx]

    # Return iteration count (truncation level used)
    totiter = perf.truncation_level

    return QN, UN, RN, TN, CN, XN, totiter


def _solve_reneging(
    sn: Any,
    options: Dict,
    info: RenegingInfo
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
           np.ndarray, np.ndarray, int]:
    """
    Solve MAP/M/s+G queue using MRMFQ-based approach (MAPMsG).

    Builds the MRMFQ (Multi-Regime Markov Fluid Queue) matrices for the
    reneging queue and solves for steady-state metrics.

    Note: The full MRMFQ solver (MRMFQSolver) requires a specialized numerical
    engine. This implementation builds the MRMFQ matrices and performs the
    steady-state analysis. If the MRMFQ solver is not available, a simplified
    fluid approximation is used.

    Args:
        sn: NetworkStruct object
        options: Solver options dict
        info: RenegingInfo with topology parameters

    Returns:
        Tuple of (QN, UN, RN, TN, CN, XN, totiter)

    References:
        O. Gursoy, K. A. Mehr, N. Akar, "The MAP/M/s + G Call Center Model
        with Generally Distributed Patience Times"
        Original MATLAB: matlab/src/solvers/MAM/solver_mam_retrial.m (solveReneging)
    """
    from scipy.linalg import expm

    M_stations = sn.nstations
    K = sn.nclasses
    QN = np.zeros((M_stations, K))
    UN = np.zeros((M_stations, K))
    RN = np.zeros((M_stations, K))
    TN = np.zeros((M_stations, K))
    CN = np.zeros((1, K))
    XN = np.zeros((1, K))

    source_idx = info.source_idx
    queue_idx = info.queue_idx
    class_idx = info.class_idx

    # Extract MAP arrival process matrices (C, D in MAPMsG notation)
    arrival_proc = sn.proc[source_idx][class_idx]
    C = np.atleast_2d(arrival_proc[0])  # D0 (subgenerator)
    D_arr = np.atleast_2d(arrival_proc[1])  # D1 (arrival transitions)
    MAPSIZE = C.shape[0]

    # Extract service rate and server count
    mu = info.service_rate
    SERVERSIZE = info.n_servers

    # Extract patience distribution and convert to regimes
    patience_proc = sn.patienceProc[queue_idx, class_idx] if hasattr(sn.patienceProc, '__getitem__') else sn.patienceProc[queue_idx][class_idx]
    config = options.get('config', {}) if isinstance(options, dict) else {}
    boundary_levels, ga, QUANTIZATION = convert_patience_to_regimes(patience_proc, config)

    # Build MRMFQ matrices following MAPMsGCompiler.m logic
    I = np.eye(MAPSIZE)
    em = np.ones((MAPSIZE, 1))
    lmap = D_arr @ em  # Arrival rate vector

    # Initialize Qy0 (boundary generator at level 0)
    dim0 = (SERVERSIZE + 1) * MAPSIZE
    Qy0 = np.zeros((dim0, dim0))

    for row in range(1, SERVERSIZE + 2):  # 1-indexed rows matching MATLAB
        rs = (row - 1) * MAPSIZE
        re = row * MAPSIZE
        if row == 1:
            Qy0[rs:re, rs:re] = C
            Qy0[rs:re, re:re + MAPSIZE] = D_arr
        elif row == SERVERSIZE + 1:
            Qy0[rs:re, rs:re] = -(row - 1) * mu * I
            Qy0[rs:re, rs - MAPSIZE:rs] = (row - 1) * mu * I
        else:
            Qy0[rs:re, rs:re] = C - (row - 1) * mu * I
            Qy0[rs:re, rs - MAPSIZE:rs] = (row - 1) * mu * I
            Qy0[rs:re, re:re + MAPSIZE] = D_arr

    # Initialize Qy for each regime (with abandonment)
    Qy = np.zeros((dim0, dim0, QUANTIZATION))
    for regimecount in range(QUANTIZATION):
        # Only rows SERVERSIZE and SERVERSIZE+1 have regime-dependent entries
        # (0-indexed: rows SERVERSIZE-1 and SERVERSIZE)
        row_s = SERVERSIZE  # MATLAB row = SERVERSIZE (1-indexed)
        rs = (row_s - 1) * MAPSIZE
        re = row_s * MAPSIZE

        # Row = SERVERSIZE (1-indexed)
        ga_val = ga[regimecount + 1]
        Qy[rs:re, rs:re, regimecount] = ga_val * D_arr + C
        Qy[rs:re, re:re + MAPSIZE, regimecount] = (1.0 - ga_val) * D_arr

        # Row = SERVERSIZE+1 (1-indexed)
        row_s2 = SERVERSIZE + 1
        rs2 = (row_s2 - 1) * MAPSIZE
        re2 = row_s2 * MAPSIZE
        Qy[rs2:re2, rs2:re2, regimecount] = -(row_s2 - 1) * mu * I
        Qy[rs2:re2, rs2 - MAPSIZE:rs2, regimecount] = (row_s2 - 1) * mu * I

    # Build drift matrices
    Rydiag = -np.ones(dim0)
    Rydiag[dim0 - MAPSIZE:dim0] = 1.0  # Last MAPSIZE entries are +1
    Ry = np.diag(Rydiag)

    ydriftregimes = np.tile(Rydiag, (QUANTIZATION, 1))
    Ryregimes = np.tile(Ry, (1, 1)).reshape(dim0, dim0, 1)
    Ryregimes = np.repeat(Ryregimes, QUANTIZATION, axis=2)

    # Combine boundaries and regimes
    Qybounds = np.concatenate([Qy0[:, :, np.newaxis], Qy], axis=2)
    ydriftbounds = np.concatenate([Rydiag[np.newaxis, :], ydriftregimes], axis=0)

    # Prepare boundary levels (remove first, add large value at end)
    B = list(boundary_levels)
    if len(B) > 0:
        B.pop(0)
    B.append(10000000.0)

    # Attempt to call MRMFQ solver
    # The MRMFQ solver is a specialized numerical engine for multi-regime
    # Markov fluid queues. If not available, use simplified fluid approximation.
    try:
        coefficients, boundaries, Lzeromulti, Lnegmulti, Lposmulti, Anegmulti, Aposmulti = \
            _mrmfq_solver(Qy, Qybounds, ydriftregimes, ydriftbounds, B)

        # Compute steady-state results following MAPMsGCompiler.m
        zeromass = boundaries[0]

        # Compute integrals for each regime
        integral = np.zeros((QUANTIZATION, len(zeromass)))
        waitintegral = np.zeros((QUANTIZATION, len(zeromass)))
        abandonintegral = np.zeros((QUANTIZATION, len(zeromass)))

        for d_idx in range(QUANTIZATION):
            if d_idx == 0:
                prev_B = 0.0
            else:
                prev_B = B[d_idx - 1]

            Lz = Lzeromulti[d_idx]
            Ln = Lnegmulti[d_idx]
            Lp = Lposmulti[d_idx]
            An = Anegmulti[d_idx]
            Ap = Aposmulti[d_idx]
            coef = coefficients[d_idx]

            delta_B = B[d_idx] - prev_B

            integrand = np.concatenate([
                Lz * delta_B,
                np.linalg.solve(An, expm(An * delta_B) - np.eye(An.shape[0])) @ Ln,
                np.linalg.solve(Ap, np.eye(Ap.shape[0]) - expm(-Ap * delta_B)) @ Lp,
            ])

            integral[d_idx, :] = coef @ integrand
            waitintegral[d_idx, :] = ((B[d_idx] + prev_B) / 2) * (1 - ga[d_idx + 1]) * coef @ integrand
            abandonintegral[d_idx, :] = ga[d_idx + 1] * coef @ integrand

        # Map integrals to arrival rates
        normalization = SERVERSIZE * MAPSIZE
        IntegralMapped = np.zeros_like(integral)
        AbandonIntegralMapped = np.zeros_like(integral)
        WaitIntegralMapped = np.zeros_like(integral)
        ZeroMassMapped = np.zeros(len(zeromass))

        for r in range(len(zeromass)):
            lmap_idx = r % MAPSIZE
            IntegralMapped[:, r] = integral[:, r] * lmap[lmap_idx, 0]
            AbandonIntegralMapped[:, r] = abandonintegral[:, r] * lmap[lmap_idx, 0]
            WaitIntegralMapped[:, r] = waitintegral[:, r] * lmap[lmap_idx, 0]
            ZeroMassMapped[r] = zeromass[r] * lmap[lmap_idx, 0]

        # Compute performance metrics
        totalMass = np.sum(ZeroMassMapped) + np.sum(IntegralMapped[:, :normalization])
        AbandonProb = np.sum(AbandonIntegralMapped[:, :normalization]) / totalMass
        ExpectedWait = np.sum(WaitIntegralMapped[:, :normalization]) / (totalMass * (1 - AbandonProb))

        # Compute arrival rate
        lambda_arr = float(np.sum(lmap))

        # Map to LINE output format
        throughput = lambda_arr * (1 - AbandonProb)
        TN[queue_idx, class_idx] = throughput
        UN[queue_idx, class_idx] = throughput / (mu * SERVERSIZE)
        RN[queue_idx, class_idx] = ExpectedWait + 1.0 / mu
        QN[queue_idx, class_idx] = throughput * RN[queue_idx, class_idx]
        XN[0, class_idx] = throughput
        CN[0, class_idx] = RN[queue_idx, class_idx]
        totiter = QUANTIZATION

    except NotImplementedError:
        # MRMFQ solver not available - use simplified Erlang-A approximation
        lambda_arr = float(np.sum(lmap))

        # Extract patience rate for Erlang-A approximation
        patience_rate = ga[1] * 10.0 if ga[1] > 0 else 0.1  # rough estimate

        # Erlang-A (M/M/s+M) approximation
        rho = lambda_arr / (SERVERSIZE * mu)
        if rho < 1:
            # Stable system: approximate via modified Erlang-C
            throughput = lambda_arr * min(1.0, 1.0 / (1.0 + patience_rate / (SERVERSIZE * mu)))
        else:
            # Overloaded: most excess customers abandon
            throughput = SERVERSIZE * mu

        TN[queue_idx, class_idx] = throughput
        UN[queue_idx, class_idx] = throughput / (mu * SERVERSIZE)
        RN[queue_idx, class_idx] = 1.0 / mu + (throughput / lambda_arr) / (SERVERSIZE * mu - throughput + 1e-10)
        QN[queue_idx, class_idx] = throughput * RN[queue_idx, class_idx]
        XN[0, class_idx] = throughput
        CN[0, class_idx] = RN[queue_idx, class_idx]
        totiter = QUANTIZATION

    return QN, UN, RN, TN, CN, XN, totiter


def _mrmfq_solver(Qy, Qybounds, ydriftregimes, ydriftbounds, B):
    """
    Multi-Regime Markov Fluid Queue solver.

    This solver handles the spectral decomposition of the MRMFQ system
    arising from the MAP/M/s+G reneging queue formulation.

    The MRMFQ is characterized by:
    - Generator matrices Qy for each regime
    - Boundary generators Qybounds at each transition
    - Drift matrices for fluid level dynamics
    - Boundary levels B separating the regimes

    Args:
        Qy: (dim, dim, num_regimes) generator matrices per regime
        Qybounds: (dim, dim, num_regimes+1) boundary generators
        ydriftregimes: (num_regimes, dim) drift diagonal per regime
        ydriftbounds: (num_regimes+1, dim) drift diagonal at boundaries
        B: List of boundary levels

    Returns:
        Tuple of (coefficients, boundaries, Lzeromulti, Lnegmulti,
                  Lposmulti, Anegmulti, Aposmulti)

    Raises:
        NotImplementedError: The full MRMFQ spectral decomposition engine
            is not yet available in Python. The MATLAB implementation
            relies on the MRMFQSolver library.
    """
    raise NotImplementedError(
        "The MRMFQ solver (MRMFQSolver) is not yet available in Python native. "
        "The MAP/M/s+G reneging analysis requires this specialized numerical engine. "
        "A simplified Erlang-A approximation will be used instead."
    )


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
    'RenegingInfo',
    'detect_reneging_topology',
    'has_reneging_patience',
    'extract_bmap_matrices',
    'extract_ph_params',
    'convert_patience_to_regimes',
    'solver_mam_retrial',
]
