"""
Quasi-Birth-Death (QBD) Process Utilities.

Native Python implementations for QBD matrix computations including
rate matrix R computation and QBD block construction.

Key algorithms:
    qbd_R: Rate matrix via successive substitutions
    qbd_R_logred: Rate matrix via logarithmic reduction
    qbd_rg: Compute R and G matrices

References:
    Original MATLAB: matlab/src/api/mam/qbd_*.m
    Latouche & Ramaswami, "Introduction to Matrix Analytic Methods
    in Stochastic Modeling", 1999
"""

import numpy as np
from numpy.linalg import LinAlgError
from scipy import linalg
from typing import Tuple, Optional
from dataclasses import dataclass


@dataclass
class QBDResult:
    """Result of QBD analysis."""
    R: np.ndarray  # Rate matrix R
    G: Optional[np.ndarray] = None  # Rate matrix G
    U: Optional[np.ndarray] = None  # U matrix
    eta: Optional[float] = None  # Caudal characteristic


def qbd_R(B: np.ndarray, L: np.ndarray, F: np.ndarray,
          iter_max: int = 100000, tol: float = 1e-12) -> np.ndarray:
    """
    Compute QBD rate matrix R using successive substitutions.

    Solves the matrix quadratic equation:
        R^2 * A_{-1} + R * A_0 + A_1 = 0

    where A_{-1} = B, A_0 = L, A_1 = F.

    Args:
        B: Backward transition block A_{-1}
        L: Local transition block A_0
        F: Forward transition block A_1
        iter_max: Maximum iterations (default: 100000)
        tol: Convergence tolerance (default: 1e-12)

    Returns:
        Rate matrix R

    References:
        Original MATLAB: matlab/src/api/mam/qbd_R.m
    """
    B = np.asarray(B, dtype=np.float64)
    L = np.asarray(L, dtype=np.float64)
    F = np.asarray(F, dtype=np.float64)

    try:
        L_inv = linalg.inv(L)
    except LinAlgError:
        L_inv = linalg.pinv(L)

    Fil = F @ L_inv
    BiL = B @ L_inv

    R = -Fil
    Rprime = -Fil - R @ R @ BiL

    for _ in range(iter_max):
        R = Rprime
        Rprime = -Fil - R @ R @ BiL
        if linalg.norm(R - Rprime, 1) <= tol:
            break

    return Rprime


def qbd_R_logred(B: np.ndarray, L: np.ndarray, F: np.ndarray,
                 iter_max: int = 1000, tol: float = 1e-14) -> np.ndarray:
    """
    Compute QBD rate matrix R using logarithmic reduction.

    Uses the logarithmic reduction algorithm which has quadratic
    convergence compared to linear convergence of successive substitutions.

    Args:
        B: Backward transition block A_{-1}
        L: Local transition block A_0
        F: Forward transition block A_1
        iter_max: Maximum iterations (default: 1000)
        tol: Convergence tolerance (default: 1e-14)

    Returns:
        Rate matrix R

    References:
        Original MATLAB: matlab/src/api/mam/qbd_R_logred.m
        Latouche & Ramaswami, Ch. 8
    """
    B = np.asarray(B, dtype=np.float64)
    L = np.asarray(L, dtype=np.float64)
    F = np.asarray(F, dtype=np.float64)

    n = L.shape[0]

    try:
        L_inv = linalg.inv(-L)
    except LinAlgError:
        L_inv = linalg.pinv(-L)

    A = L_inv @ F  # A_1
    C = L_inv @ B  # A_{-1}

    H = A.copy()
    T = C.copy()

    for _ in range(iter_max):
        # Compute (I - AC - CA)^{-1}
        I = np.eye(n)
        M = I - A @ C - C @ A
        try:
            M_inv = linalg.inv(M)
        except LinAlgError:
            M_inv = linalg.pinv(M)

        # Update
        A_new = M_inv @ A @ A
        C_new = M_inv @ C @ C
        H_new = H + T @ M_inv @ A

        if linalg.norm(H_new - H, 1) <= tol:
            H = H_new
            break

        A = A_new
        C = C_new
        T = T @ M_inv @ (A + C)
        H = H_new

    # R = H * (-L)
    R = H @ (-L)

    return R


def qbd_rg(B: np.ndarray, L: np.ndarray, F: np.ndarray,
           method: str = 'logred', iter_max: int = 1000,
           tol: float = 1e-14) -> QBDResult:
    """
    Compute both R and G matrices for a QBD process.

    G is the minimal non-negative solution to:
        A_1 * G^2 + A_0 * G + A_{-1} = 0

    R is the minimal non-negative solution to:
        R^2 * A_{-1} + R * A_0 + A_1 = 0

    Args:
        B: Backward transition block A_{-1}
        L: Local transition block A_0
        F: Forward transition block A_1
        method: 'logred' or 'successive' (default: 'logred')
        iter_max: Maximum iterations
        tol: Convergence tolerance

    Returns:
        QBDResult with R, G, U, and eta (caudal characteristic)

    References:
        Original MATLAB: matlab/src/api/mam/qbd_rg.m
    """
    B = np.asarray(B, dtype=np.float64)
    L = np.asarray(L, dtype=np.float64)
    F = np.asarray(F, dtype=np.float64)

    n = L.shape[0]

    # Compute R
    if method == 'logred':
        R = qbd_R_logred(B, L, F, iter_max, tol)
    else:
        R = qbd_R(B, L, F, iter_max, tol)

    # Compute G using the relation: G = A_{-1} * (R*A_{-1} + A_0)^{-1}
    # Or alternatively iterate: G_{n+1} = -(A_0 + A_1*G_n^2)^{-1} * A_{-1}
    try:
        L_inv = linalg.inv(-L)
    except LinAlgError:
        L_inv = linalg.pinv(-L)

    A = L_inv @ F
    C = L_inv @ B

    # Use logarithmic reduction for G
    G = C.copy()
    T = A.copy()

    for _ in range(iter_max):
        I = np.eye(n)
        M = I - A @ C - C @ A
        try:
            M_inv = linalg.inv(M)
        except LinAlgError:
            M_inv = linalg.pinv(M)

        A_new = M_inv @ A @ A
        C_new = M_inv @ C @ C
        G_new = G + T @ M_inv @ C

        if linalg.norm(G_new - G, 1) <= tol:
            G = G_new
            break

        A = A_new
        C = C_new
        T = T @ M_inv @ (A + C)
        G = G_new

    G = G @ (-L)

    # Compute U = A_0 + A_1 * G
    U = L + F @ G

    # Caudal characteristic (spectral radius of R)
    try:
        eigvals = linalg.eigvals(R)
        eta = np.max(np.abs(eigvals))
    except:
        eta = None

    return QBDResult(R=R, G=G, U=U, eta=eta)


def qbd_blocks_mapmap1(D0_arr: np.ndarray, D1_arr: np.ndarray,
                       D0_srv: np.ndarray, D1_srv: np.ndarray
                       ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Construct QBD blocks for a MAP/MAP/1 queue.

    Builds the backward (B), local (L), and forward (F) transition
    blocks for the QBD representation of a MAP/MAP/1 queue.

    Args:
        D0_arr: Arrival MAP D0 matrix
        D1_arr: Arrival MAP D1 matrix
        D0_srv: Service MAP D0 matrix
        D1_srv: Service MAP D1 matrix

    Returns:
        Tuple of (B, L, F) QBD blocks

    References:
        Original MATLAB: matlab/src/api/mam/qbd_mapmap1.m
    """
    D0_arr = np.asarray(D0_arr, dtype=np.float64)
    D1_arr = np.asarray(D1_arr, dtype=np.float64)
    D0_srv = np.asarray(D0_srv, dtype=np.float64)
    D1_srv = np.asarray(D1_srv, dtype=np.float64)

    na = D0_arr.shape[0]
    ns = D0_srv.shape[0]

    # Forward transitions (arrivals): F = D1_arr \otimes I_ns
    F = np.kron(D1_arr, np.eye(ns))

    # Backward transitions (service completions): B = I_na \otimes D1_srv
    B = np.kron(np.eye(na), D1_srv)

    # Local transitions: L = D0_arr \otimes I_ns + I_na \otimes D0_srv
    L = np.kron(D0_arr, np.eye(ns)) + np.kron(np.eye(na), D0_srv)

    return B, L, F


def qbd_bmapbmap1(
    MAPa: Tuple[np.ndarray, np.ndarray],
    pbatcha: np.ndarray,
    MAPs: Tuple[np.ndarray, np.ndarray]
) -> Tuple[np.ndarray, np.ndarray, list, np.ndarray, list]:
    """
    Compute QBD blocks for a BMAP/BMAP/1 queue.

    Constructs the QBD (Quasi-Birth-Death) transition blocks for a
    BMAP/BMAP/1 queue with batch arrivals.

    Args:
        MAPa: Arrival process MAP as (D0, D1)
        pbatcha: Probability distribution of batch sizes (array of length maxbatch)
        MAPs: Service process MAP as (D0, D1)

    Returns:
        Tuple of (A0, A_1, A1_list, B0, B1_list) where:
            A0: Local transition block
            A_1: Downward transition block
            A1_list: List of upward transition blocks for each batch size
            B0: Initial boundary local block
            B1_list: List of boundary upward blocks for each batch size

    References:
        Original MATLAB: matlab/src/api/mam/qbd_bmapbmap1.m
    """
    D0_arr, D1_arr = MAPa
    D0_srv, D1_srv = MAPs

    D0_arr = np.asarray(D0_arr, dtype=np.float64)
    D1_arr = np.asarray(D1_arr, dtype=np.float64)
    D0_srv = np.asarray(D0_srv, dtype=np.float64)
    D1_srv = np.asarray(D1_srv, dtype=np.float64)
    pbatcha = np.asarray(pbatcha, dtype=np.float64)

    na = D0_arr.shape[0]
    ns = D0_srv.shape[0]
    maxbatch = len(pbatcha)

    # Build upward transition blocks for each batch size
    A1_list = []
    for b in range(maxbatch):
        A1_b = np.kron(D1_arr * pbatcha[b], np.eye(ns))
        A1_list.append(A1_b)

    # Local transitions: A0 = D0_arr \otimes I_ns + I_na \otimes D0_srv (Kronecker sum)
    A0 = np.kron(D0_arr, np.eye(ns)) + np.kron(np.eye(na), D0_srv)

    # Downward transitions: A_1 = I_na \otimes D1_srv
    A_1 = np.kron(np.eye(na), D1_srv)

    # Boundary blocks (for level 0)
    # MATLAB: B0 = krons(MAPa{1}, eye(ns)) = kron(D0_arr, I_ns) + kron(I_na, I_ns)
    B0 = np.kron(D0_arr, np.eye(ns)) + np.eye(na * ns)

    B1_list = []
    for b in range(maxbatch):
        B1_b = np.kron(D1_arr * pbatcha[b], np.eye(ns))
        B1_list.append(B1_b)

    return A0, A_1, A1_list, B0, B1_list


def qbd_mapmap1(
    MAPa: Tuple[np.ndarray, np.ndarray],
    MAPs: Tuple[np.ndarray, np.ndarray],
    util: Optional[float] = None
) -> Tuple[float, float, float, np.ndarray, np.ndarray, Optional[float],
           np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
           Tuple[np.ndarray, np.ndarray]]:
    """
    Analyze a MAP/MAP/1 queue using QBD methods.

    Solves a MAP/MAP/1 queue using Quasi-Birth-Death process methods,
    computing throughput, queue length, utilization, and other metrics.

    Args:
        MAPa: Arrival process MAP as (D0, D1)
        MAPs: Service process MAP as (D0, D1)
        util: Optional target utilization to scale service rate

    Returns:
        Tuple of (XN, QN, UN, pqueue, R, eta, G, A_1, A0, A1, U, MAPs_scaled)
        where:
            XN: System throughput
            QN: Mean queue length
            UN: Utilization
            pqueue: Queue length distribution
            R: Rate matrix R
            eta: Caudal characteristic (spectral radius of R)
            G: Rate matrix G
            A_1: Downward transition block
            A0: Local transition block
            A1: Upward transition block
            U: Matrix U
            MAPs_scaled: Scaled service process

    References:
        Original MATLAB: matlab/src/api/mam/qbd_mapmap1.m
    """
    from ..mam import map_scale, map_lambda

    D0_arr, D1_arr = MAPa
    D0_srv, D1_srv = MAPs

    D0_arr = np.asarray(D0_arr, dtype=np.float64)
    D1_arr = np.asarray(D1_arr, dtype=np.float64)
    D0_srv = np.asarray(D0_srv, dtype=np.float64)
    D1_srv = np.asarray(D1_srv, dtype=np.float64)

    na = D0_arr.shape[0]
    ns = D0_srv.shape[0]

    # Scale service process if target utilization provided
    if util is not None:
        lambda_a = map_lambda(D0_arr, D1_arr)
        D0_srv, D1_srv = map_scale(D0_srv, D1_srv, util / lambda_a)

    lambda_a = map_lambda(D0_arr, D1_arr)
    lambda_s = map_lambda(D0_srv, D1_srv)
    actual_util = lambda_a / lambda_s

    # Build QBD blocks
    A1 = np.kron(D1_arr, np.eye(ns))  # Forward (arrivals)
    A0 = np.kron(D0_arr, np.eye(ns)) + np.kron(np.eye(na), D0_srv)  # Local
    A_1 = np.kron(np.eye(na), D1_srv)  # Backward (services)
    A0bar = np.kron(D0_arr, np.eye(ns))  # Boundary local

    # Solve QBD using Cyclic Reduction (matching MATLAB: QBD_CR(A_1, A0, A1))
    from ..smc import qbd_cr, qbd_pi
    cr_result = qbd_cr(A_1, A0, A1)
    G = cr_result['G']
    R = cr_result['R']
    U = cr_result['U']

    # Compute caudal characteristic (eta = spectral radius of R)
    eta_val = np.max(np.abs(linalg.eigvals(R)))

    # Compute queue length distribution using QBD_pi
    # MATLAB: pqueue = QBD_pi(A_1, A0bar, R, 'MaxNumComp', 1e2)
    pi_flat = qbd_pi(A_1, A0bar, R, max_num_comp=100)
    n_phases = na * ns
    num_levels = len(pi_flat) // n_phases
    pqueue = pi_flat.reshape(num_levels, n_phases)

    # Retry with more components if needed (MATLAB line 75-77)
    if np.sum(np.sum(pqueue[1:, :], axis=1)) < actual_util * 0.99:
        pi_flat = qbd_pi(A_1, A0bar, R, max_num_comp=20000)
        num_levels = len(pi_flat) // n_phases
        pqueue = pi_flat.reshape(num_levels, n_phases)

    # Compute performance measures (matching MATLAB lines 79-88)
    if na == 1 and ns == 1:
        UN = 1.0 - pqueue[0, 0]
        QN = float(np.arange(pqueue.shape[0]) @ pqueue.flatten())
    else:
        UN = 1.0 - np.sum(pqueue[0, :])
        QN = float(np.arange(pqueue.shape[0]) @ np.sum(pqueue, axis=1))

    XN = lambda_a
    MAPs_scaled = (D0_srv, D1_srv)

    return XN, QN, UN, pqueue, R, eta_val, G, A_1, A0, A1, U, MAPs_scaled


def qbd_raprap1(
    RAPa: Tuple[np.ndarray, np.ndarray],
    RAPs: Tuple[np.ndarray, np.ndarray],
    util: Optional[float] = None
) -> Tuple[float, float, float, np.ndarray, np.ndarray, float,
           np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Analyze a RAP/RAP/1 queue using QBD methods.

    Solves a RAP/RAP/1 queue (Rational Arrival Process) using QBD methods,
    computing throughput, queue length, utilization, and other metrics.

    Note: This is a simplified implementation. For full RAP analysis,
    specialized algorithms like those in SMCSolver are needed.

    Args:
        RAPa: Arrival process RAP as (H0, H1)
        RAPs: Service process RAP as (H0, H1)
        util: Optional target utilization to scale service rate

    Returns:
        Tuple of (XN, QN, UN, pqueue, R, eta, G, B, L, F) where:
            XN: System throughput
            QN: Mean queue length
            UN: Utilization
            pqueue: Queue length distribution
            R: Rate matrix R
            eta: Caudal characteristic
            G: Rate matrix G
            B: Backward transition block
            L: Local transition block
            F: Forward transition block

    References:
        Original MATLAB: matlab/src/api/mam/qbd_raprap1.m
    """
    from ..mam import map_scale, map_lambda

    H0_arr, H1_arr = RAPa
    H0_srv, H1_srv = RAPs

    H0_arr = np.asarray(H0_arr, dtype=np.float64)
    H1_arr = np.asarray(H1_arr, dtype=np.float64)
    H0_srv = np.asarray(H0_srv, dtype=np.float64)
    H1_srv = np.asarray(H1_srv, dtype=np.float64)

    na = H0_arr.shape[0]
    ns = H0_srv.shape[0]

    # Scale service process if target utilization provided
    if util is not None:
        lambda_a = map_lambda(H0_arr, H1_arr)
        H0_srv, H1_srv = map_scale(H0_srv, H1_srv, util / lambda_a)

    lambda_a = map_lambda(H0_arr, H1_arr)
    lambda_s = map_lambda(H0_srv, H1_srv)
    actual_util = lambda_a / lambda_s

    # Build QBD blocks matching MATLAB Q_RAP_RAP_1.m
    # A0 = kron(eye(mA), D1)   -- backward (service completions)
    # A1 = kron(C0, eye(mS)) + kron(eye(mA), D0)  -- local (Kronecker sum)
    # A2 = kron(C1, eye(mS))   -- forward (arrivals)
    A0 = np.kron(np.eye(na), H1_srv)        # B (backward)
    A1 = np.kron(H0_arr, np.eye(ns)) + np.kron(np.eye(na), H0_srv)  # L (local)
    A2 = np.kron(H1_arr, np.eye(ns))        # F (forward)
    B0 = A0.copy()
    B1 = np.kron(H0_arr, np.eye(ns))        # Boundary local

    # Solve QBD using Cyclic Reduction with RAP components
    from ..smc import qbd_cr, qbd_pi
    cr_result = qbd_cr(A0, A1, A2, rap_comp=True)
    G = cr_result['G']
    R = cr_result['R']

    # Compute queue length distribution using QBD_pi with RAPComp
    pi_flat = qbd_pi(B0, B1, R, max_num_comp=50, rap_comp=True)
    n_phases = na * ns
    num_levels = len(pi_flat) // n_phases

    # Compute queue length distribution (probabilities per level)
    ql = np.zeros(num_levels)
    for i in range(num_levels):
        ql[i] = np.sum(pi_flat[i * n_phases:(i + 1) * n_phases])

    # Reshape pqueue as levels x phases
    pqueue = pi_flat.reshape(num_levels, n_phases)

    # Compute performance measures matching MATLAB qbd_raprap1.m
    eta = float(np.max(np.abs(linalg.eigvals(R))))

    if na == 1 and ns == 1:
        UN = 1.0 - pqueue[0, 0]
    else:
        UN = 1.0 - np.sum(pqueue[0, :])

    QN = float(np.arange(pqueue.shape[0]) @ np.sum(pqueue, axis=1))
    XN = lambda_a

    # Return B, L, F using original naming convention
    B = A0
    L = A1
    F = A2

    return XN, QN, UN, pqueue, R, eta, G, B, L, F


def qbd_setupdelayoff(
    lambda_val: float,
    mu: float,
    alpharate: float,
    alphascv: float,
    betarate: float,
    betascv: float
) -> float:
    """
    Analyze queue with setup delay and turn-off phases.

    Performs queue-length analysis for a queueing system with setup
    delay (warm-up) and turn-off periods using QBD methods.

    The system operates as follows:
    1. When empty and job arrives, server enters setup phase
    2. After setup, server becomes active and serves jobs
    3. When queue empties, server enters turn-off phase
    4. After turn-off, server becomes idle

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        alpharate: Rate of setup delay phase
        alphascv: Squared coefficient of variation for setup delay
        betarate: Rate of turn-off phase
        betascv: Squared coefficient of variation for turn-off period

    Returns:
        Average queue length QN

    References:
        Original MATLAB: matlab/src/api/mam/qbd_setupdelayoff.m
    """
    # Try to import AcyclicPHFromMeansAndSCVs, but it may not exist
    try:
        from ..butools.ph.baseph import AcyclicPHFromMeansAndSCVs
        HAS_ACYCLIC_PH = True
    except ImportError:
        HAS_ACYCLIC_PH = False
        AcyclicPHFromMeansAndSCVs = None

    # Fit PH distributions for setup and turn-off phases
    alpha_D0 = None
    if HAS_ACYCLIC_PH:
        try:
            alpha_ph = AcyclicPHFromMeansAndSCVs([1.0 / alpharate], [alphascv])
            alpha_D0 = alpha_ph[1]  # Subgenerator matrix
            alpha_D1 = -alpha_D0 @ np.ones((alpha_D0.shape[0], 1))
        except:
            alpha_D0 = None

    if alpha_D0 is None:
        # Fallback to Erlang approximation
        na = max(1, int(round(1.0 / alphascv)))
        rate_a = na * alpharate
        alpha_D0 = np.diag([-rate_a] * na)
        for i in range(na - 1):
            alpha_D0[i, i + 1] = rate_a
        alpha_D1 = np.zeros((na, 1))
        alpha_D1[-1, 0] = rate_a

    beta_D0 = None
    if HAS_ACYCLIC_PH:
        try:
            beta_ph = AcyclicPHFromMeansAndSCVs([1.0 / betarate], [betascv])
            beta_D0 = beta_ph[1]
            beta_D1 = -beta_D0 @ np.ones((beta_D0.shape[0], 1))
        except:
            beta_D0 = None

    if beta_D0 is None:
        # Fallback to Erlang approximation
        nb = max(1, int(round(1.0 / betascv)))
        rate_b = nb * betarate
        beta_D0 = np.diag([-rate_b] * nb)
        for i in range(nb - 1):
            beta_D0[i, i + 1] = rate_b
        beta_D1 = np.zeros((nb, 1))
        beta_D1[-1, 0] = rate_b

    na = alpha_D0.shape[0]
    nb = beta_D0.shape[0]
    n = na + nb

    # Build QBD blocks
    # States: [setup phases (1..na), active + turn-off phases (na+1..n)]
    F = np.zeros((n, n))  # Forward transitions (arrivals)
    B = np.zeros((n, n))  # Backward transitions (service)

    # Arrivals in setup phase
    for i in range(na):
        F[i, i] = lambda_val

    # Arrivals in turn-off phase (go to active state)
    for i in range(nb):
        F[na + i, na] = lambda_val
    F[na, na] = lambda_val

    # Service completions (only from active state)
    B[na, na] = mu

    # Local transitions
    L = np.zeros((n, n))

    # Setup phase transitions
    for i in range(na):
        L[i, i] = alpha_D0[i, i] - lambda_val
        if i < na - 1:
            L[i, i + 1:na] = alpha_D0[i, i + 1:na]
        else:
            # Transition from last setup phase to active
            L[na - 1, na] = -alpha_D0[na - 1, na - 1]

    # Active state
    L[na, na] = -mu - lambda_val

    # Turn-off phase (only reachable from level 0)
    for i in range(1, nb):
        L[na + i, na + i] = -lambda_val

    # Boundary block L0 for level 0
    L0 = np.zeros((n, n))

    # Setup phase at level 0
    for i in range(na):
        L0[i, i] = -lambda_val

    # Turn-off phase at level 0
    for i in range(nb):
        L0[na + i, na + i] = beta_D0[i, i] - lambda_val if i < nb else -lambda_val
        if i == nb - 1:
            L0[na + i, 0] = -beta_D0[i, i]  # Return to setup phase
        elif i < nb - 1:
            L0[na + i, na + i + 1] = -beta_D0[i, i]

    # Compute R matrix using QBD_CR equivalent
    R = qbd_R(B, L, F)

    # Compute steady-state distribution using QBD_pi algorithm
    # Follow MATLAB QBD_pi: convert to discrete time first
    I = np.eye(n)

    # Uniformization: find maximum exit rate from boundary block
    lamb = max(-np.diag(L0))
    if lamb <= 0:
        lamb = 1.0

    # Convert to discrete time stochastic matrices
    B1_dt = L0 / lamb + I  # Boundary local block
    B0_dt = B / lamb       # Backward transitions

    # Compute stochastic matrix for level 0
    stat_matrix = B1_dt + R @ B0_dt

    # Find stationary distribution using stat() approach:
    # Solve: K @ [A - I, e] = [0, ..., 0, 1]
    # This is equivalent to: K @ (A - I) = 0 and K @ e = 1
    e = np.ones((n, 1))
    aug_matrix = np.hstack([stat_matrix - I, e])
    y = np.zeros(n + 1)
    y[-1] = 1.0

    # Solve K @ aug_matrix = y using least squares (K = y @ pinv(aug_matrix))
    try:
        pi0 = linalg.lstsq(aug_matrix.T, y, cond=None)[0]
    except:
        pi0 = np.linalg.lstsq(aug_matrix.T, y, rcond=None)[0]
    pi0 = np.abs(pi0)  # Ensure non-negative

    # Normalize using QBD normalization: pi @ (I-R)^{-1} @ 1 = 1
    try:
        temp = linalg.inv(I - R)
    except:
        temp = linalg.pinv(I - R)

    norm_const = pi0 @ temp @ np.ones(n)
    if norm_const > 0:
        pi0 = pi0 / norm_const

    # Build full probability vector pn following MATLAB QBD_pi
    # Generate level probabilities until total mass approaches 1
    max_num_comp = 500
    pi_levels = [pi0]
    sum_pi = np.sum(pi0)
    numit = 1

    while sum_pi < 1 - 1e-10 and numit < max_num_comp:
        pi_next = pi_levels[-1] @ R
        pi_levels.append(pi_next)
        numit += 1
        sum_pi += np.sum(pi_next)

    # Concatenate all levels into a single vector (like MATLAB's reshape(pi', 1, []))
    pn = np.concatenate(pi_levels)

    # Compute QN using MATLAB's exact indexing
    # MATLAB code: j = n+1; QN = QN + ni*sum(pn(j:(j+n))); j = j+n
    # MATLAB's pn(j:(j+n)) includes (n+1) elements due to inclusive 1-based indexing
    # Python equivalent: pn[j:j+n+1] where j starts at n (0-based)
    QN = 0.0
    j = n  # 0-based index (equivalent to MATLAB's j=n+1 in 1-based)
    ni = 0

    while True:
        ni += 1
        # MATLAB: sum(pn(j:(j+n))) includes elements j to j+n (inclusive)
        # Python: pn[j:j+n+1] gives same range
        end_idx = min(j + n + 1, len(pn))
        QN += ni * np.sum(pn[j:end_idx])
        j += n
        if j + n >= len(pn):
            break

    return QN


__all__ = [
    'QBDResult',
    'qbd_R',
    'qbd_R_logred',
    'qbd_rg',
    'qbd_blocks_mapmap1',
    'qbd_bmapbmap1',
    'qbd_mapmap1',
    'qbd_raprap1',
    'qbd_setupdelayoff',
]
