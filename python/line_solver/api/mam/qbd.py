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

from ...constants import GlobalConstants


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

    r = L.shape[0]
    eye_r = np.eye(r)

    try:
        Linv = linalg.inv(L)
    except LinAlgError:
        Linv = linalg.pinv(L)

    iLF = -Linv @ F
    iLB = -Linv @ B
    T = iLF.copy()
    S = iLB.copy()

    for _ in range(iter_max):
        D = iLF @ iLB + iLB @ iLF
        try:
            Minv = linalg.inv(eye_r - D)
        except LinAlgError:
            Minv = linalg.pinv(eye_r - D)
        iLF = Minv @ iLF @ iLF
        iLB = Minv @ iLB @ iLB
        S = S + T @ iLB
        T = T @ iLF
        if linalg.norm(np.ones(r) - S @ np.ones(r), 1) <= tol:
            break

    # S is the G matrix; U is the taboo generator of a level, R = -F U^-1
    U = L + F @ S
    return -F @ linalg.inv(U)


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
        where::

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
        # map_scale takes the TARGET MEAN: a service mean of util/lambda_a
        # gives the requested utilization against this arrival rate.
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
    from line_solver.lib.thirdparty.smc import qbd_cr, qbd_pi
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

    References:
        N. G. Bean and B. F. Nielsen, "Quasi-Birth-and-Death Processes with
        Rational Arrival Process Components", Stochastic Models, 26(3), 2010,
        pp. 309-334. The analysis rests on the prediction-process
        interpretation of a RAP due to Asmussen and Bladt, which is what
        allows a QBD argument to be carried over to matrices that are not
        nonnegative. The same prediction process underlies the
        conditional-vector RAP sampler in RAP.sample.

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
        # As above: the argument is the target service mean.
        H0_srv, H1_srv = map_scale(H0_srv, H1_srv, util / lambda_a)

    lambda_a = map_lambda(H0_arr, H1_arr)

    # see _kb/03-api-layer.md for rationale
    B_blk = np.kron(np.eye(na), H1_srv)     # backward, service completions
    L_blk = np.kron(H0_arr, np.eye(ns)) + np.kron(np.eye(na), H0_srv)   # local
    F_blk = np.kron(H1_arr, np.eye(ns))     # forward, arrivals
    B1 = np.kron(H0_arr, np.eye(ns))        # boundary local

    # see _kb/03-api-layer.md for rationale
    core = qbd_rap(F_blk, L_blk, B_blk, F_blk, B1, 0)
    R = core.R
    G = core.G

    # see _kb/03-api-layer.md for rationale
    n_phases = na * ns
    max_num_comp = 100
    levels = [core.pi0.reshape(1, n_phases)]
    sumpi = float(np.sum(levels[0]))
    numit = 1
    while sumpi < 1.0 - 1e-10 and numit < 1 + max_num_comp:
        levels.append(levels[numit - 1] @ R)
        numit += 1
        sumpi += float(np.sum(levels[numit - 1]))
    num_levels = len(levels)
    pqueue = np.vstack(levels)

    # Compute performance measures matching MATLAB qbd_raprap1.m
    eta = core.spr

    if na == 1 and ns == 1:
        UN = 1.0 - pqueue[0, 0]
    else:
        UN = 1.0 - np.sum(pqueue[0, :])

    QN = float(np.arange(pqueue.shape[0]) @ np.sum(pqueue, axis=1))
    XN = lambda_a

    # Return B, L, F using original naming convention
    B = B_blk
    L = L_blk
    F = F_blk

    return XN, QN, UN, pqueue, R, eta, G, B, L, F


@dataclass
class QbdRapResult:
    """
    Result of the equilibrium analysis of a QBD with RAP components.

    Attributes:
        levelProb: Marginal level probabilities, levels 0..numLevels
        QN: Mean queue length, computed exactly as pi0*R*inv(I-R)^2*e
        R: Rate matrix R = A0*inv(-U)
        G: Matrix G solving A0*G^2 + A1*G + A2 = 0
        U: Matrix U = A1 + A0*G
        spr: Spectral radius Sp(R); positive recurrent iff Sp(R) < 1
        pqueue: (numLevels+1) x m array whose n-th row is the level vector pi_n
        pi0: Level-0 vector pi_0, the boundary vector of Theorem 7
    """
    levelProb: np.ndarray
    QN: float
    R: np.ndarray
    G: np.ndarray
    U: np.ndarray
    spr: float
    pqueue: np.ndarray
    pi0: np.ndarray


def qbd_rap_g(A0: np.ndarray, A1: np.ndarray, A2: np.ndarray,
              block_scale: float) -> np.ndarray:
    """
    Solve A0*G^2 + A1*G + A2 = 0 for the matrix G.

    Uses the exact rank-one closed form when A2 has rank one, and otherwise
    natural functional iteration as a warm start followed by Newton's method on
    the Sylvester-form Jacobian. Never returns an unconverged iterate.

    Args:
        A0: Level-up block
        A1: Local block
        A2: Level-down block
        block_scale: Scale of the blocks, used to set the residual tolerance

    Returns:
        The matrix G

    Raises:
        ValueError: if G cannot be computed to roundoff level
    """
    m = A1.shape[0]
    e = np.ones((m, 1))
    res_tol = 1e-10 * block_scale

    # see _kb/03-api-layer.md for rationale
    _, s2, Vh2 = linalg.svd(A2)
    if len(s2) > 1 and s2[0] > 0 and s2[1] <= 1e-10 * s2[0]:
        v = Vh2[0, :].reshape(1, m)
        ve = float((v @ e).item())
        if abs(ve) < 1e-12 * np.max(np.abs(v)):
            raise ValueError('A2 has rank one but its right factor v satisfies v*e = 0, '
                             'so the closed form G = e*v/(v*e) is undefined.')
        G = e @ (v / ve)
        res = float(linalg.norm(A0 @ G @ G + A1 @ G + A2, 'fro'))
        if res > res_tol:
            raise ValueError('The rank-one closed form for G leaves a residual '
                             '||A0*G^2 + A1*G + A2||_F = %g, which is above the roundoff '
                             'level %g.' % (res, res_tol))
        return G

    # see _kb/03-api-layer.md for rationale
    try:
        mA1inv = linalg.inv(-A1)
    except LinAlgError:
        raise ValueError('The local block A1 is singular, the iteration for G cannot be started.')
    G = np.zeros((m, m))
    for _ in range(200):
        Gnew = mA1inv @ (A2 + A0 @ G @ G)
        if not np.all(np.isfinite(Gnew)):
            break
        step = float(linalg.norm(Gnew - G, 'fro'))
        G = Gnew
        if step <= 1e-14 * max(float(linalg.norm(G, 'fro')), 1.0):
            break
    if not np.all(np.isfinite(G)):
        G = np.zeros((m, m))

    # see _kb/03-api-layer.md for rationale
    Im = np.eye(m)
    for _ in range(100):
        res_mat = A0 @ G @ G + A1 @ G + A2
        if linalg.norm(res_mat, 'fro') <= res_tol:
            break
        J = np.kron(Im, A0 @ G + A1) + np.kron(G.T, A0)
        try:
            h = linalg.solve(J, -res_mat.reshape(-1, order='F'))
        except (LinAlgError, ValueError):
            break
        H = h.reshape(m, m, order='F')
        G = G + H
        if not np.all(np.isfinite(G)):
            break

    res = np.inf
    ge_err = np.inf
    if np.all(np.isfinite(G)):
        res = float(linalg.norm(A0 @ G @ G + A1 @ G + A2, 'fro'))
        ge_err = float(np.max(np.abs(G @ e - e)))
    if not (res <= res_tol) or ge_err > 1e-8:
        raise ValueError(
            'Could not compute the matrix G for this QBD with RAP components: residual '
            '||A0*G^2 + A1*G + A2||_F = %g against a tolerance of %g, and ||G*e-e||_inf = %g. '
            'The blocks are not nonnegative, so neither the functional iteration nor '
            "Newton's method is guaranteed to converge, and the justification of algorithms "
            'for G in this setting is left as an open problem in Section 6 of N. G. Bean and '
            'B. F. Nielsen, "Quasi-Birth-and-Death Processes with Rational Arrival Process '
            'Components", Stochastic Models, 26(3), 2010, pp. 309-334. Supply a model with a '
            'rank-one A2, for which G is available in closed form.' % (res, res_tol, ge_err))
    return G


def qbd_rap(A0: np.ndarray, A1: np.ndarray, A2: np.ndarray,
            B0: Optional[np.ndarray] = None, B1: Optional[np.ndarray] = None,
            numLevels: int = 20) -> QbdRapResult:
    """
    Equilibrium analysis of a Quasi-Birth-and-Death process with Rational
    Arrival Process (RAP) components.

    The process is specified directly by its level-independent blocks
    (A0,A1,A2) and its boundary blocks (B0,B1), where A0 drives level
    increases, A2 drives level decreases and A1 the within-level evolution.
    Unlike a Markovian QBD the blocks need not be nonnegative: they are only
    required to be conservative, (A0+A1+A2)*e = 0, and to define a genuine RAP
    through the prediction-process interpretation. This makes qbd_rap strictly
    more general than qbd_raprap1, which builds a product-space QBD from two
    INDEPENDENT RAPs; here the arrival process and the sequence of service
    times may be driven from a shared phase space and therefore be
    cross-correlated.

    This is the block-level core of the RAP QBD family. qbd_raprap1 is the thin
    wrapper over it that builds the product-space blocks of two independent
    RAPs; callers with a coupled model must call qbd_rap directly because no
    product form exists to factor out.

    Algorithm (Theorem 7 of the reference)::

        1. Solve A0*G^2 + A1*G + A2 = 0 for G.
        2. U = A1 + A0*G.
        3. R = A0*inv(-U).
        4. Find the row vector pihat0 with pihat0*(B1 + R*A2) = 0, pihat0*e = 1.
        5. pi0 = K*pihat0 with K chosen so that pi0*inv(I-R)*e = 1.
        6. pi_n = pi0*R^n, and the marginal level probability is pi_n*e.

    The process is positive recurrent iff Sp(R) < 1 and step 4 has a solution.

    Computation of G: the blocks are not nonnegative, so the probabilistic
    iterations used for Markovian QBDs (logarithmic reduction, cyclic
    reduction) carry no convergence guarantee here, and the paper explicitly
    leaves the general case open ("The issue of justifying algorithms for the
    evaluation of the matrix G for such processes has not been undertaken",
    Section 6). See qbd_rap_g: the rank-one closed form is used when it
    applies, otherwise functional iteration followed by Newton's method, and an
    unconverged G is never returned.

    References:
        N. G. Bean and B. F. Nielsen, "Quasi-Birth-and-Death Processes with
        Rational Arrival Process Components", Stochastic Models, 26(3), 2010,
        pp. 309-334 (DTU technical report IMM-2007-20). The argument rests on
        the prediction-process interpretation of a RAP due to Asmussen and
        Bladt, which is what allows a QBD argument to be carried over to
        matrices that are not nonnegative; the same prediction process
        underlies the conditional-vector RAP sampler in RAP.sample.

        Original MATLAB: matlab/src/api/mam/qbd_rap.m

    Args:
        A0: Level-up block (m x m)
        A1: Local block (m x m)
        A2: Level-down block (m x m)
        B0: Boundary level-up block, defaults to A0
        B1: Boundary local block, defaults to A1
        numLevels: Highest level reported, default 20

    Returns:
        QbdRapResult with levelProb, QN, R, G, U, spr, pqueue and pi0

    Raises:
        ValueError: if the blocks are not conservative, if the process is not
            positive recurrent, or if G cannot be computed
    """
    A0 = np.asarray(A0, dtype=np.float64)
    A1 = np.asarray(A1, dtype=np.float64)
    A2 = np.asarray(A2, dtype=np.float64)
    B0 = A0.copy() if B0 is None else np.asarray(B0, dtype=np.float64)
    B1 = A1.copy() if B1 is None else np.asarray(B1, dtype=np.float64)

    m = A1.shape[0]
    for blk, name in ((A0, 'A0'), (A1, 'A1'), (A2, 'A2'), (B0, 'B0'), (B1, 'B1')):
        if blk.shape != (m, m):
            raise ValueError('All QBD blocks must be square and of the same order; '
                             '%s has shape %s.' % (name, blk.shape))
    if numLevels < 0 or int(numLevels) != numLevels:
        raise ValueError('numLevels must be a nonnegative integer.')
    numLevels = int(numLevels)

    e = np.ones((m, 1))
    I = np.eye(m)
    block_scale = max(1.0, float(linalg.norm(A0, 'fro')),
                      float(linalg.norm(A1, 'fro')), float(linalg.norm(A2, 'fro')))

    # see _kb/03-api-layer.md for rationale
    cons_a = float(np.max(np.abs((A0 + A1 + A2) @ e)))
    if cons_a > 1e-8 * block_scale:
        raise ValueError('The repeating blocks are not conservative: ||(A0+A1+A2)*e||_inf = %g. '
                         'A QBD with RAP components requires (A0+A1+A2)*e = 0.' % cons_a)
    cons_b = float(np.max(np.abs((B0 + B1) @ e)))
    if cons_b > 1e-8 * block_scale:
        raise ValueError('The boundary blocks are not conservative: ||(B0+B1)*e||_inf = %g. '
                         'A QBD with RAP components requires (B0+B1)*e = 0 at level 0.' % cons_b)

    # Step 1: matrix G.
    G = qbd_rap_g(A0, A1, A2, block_scale)

    # Steps 2 and 3: U and R.
    U = A1 + A0 @ G
    try:
        R = A0 @ linalg.inv(-U)
    except LinAlgError:
        raise ValueError('The matrix U = A1 + A0*G is singular, R = A0*inv(-U) does not exist.')

    # Corollary 8(i): positive recurrence.
    spr = float(np.max(np.abs(linalg.eigvals(R))))
    # see _kb/03-api-layer.md for rationale
    if spr >= 1.0 - 1e-12:
        raise ValueError('The process is not positive recurrent: Sp(R) = %.15g >= 1 '
                         '(Corollary 8 of Bean and Nielsen, 2010).' % spr)

    # see _kb/03-api-layer.md for rationale
    V = B1 + R @ A2
    _, sv, Wh = linalg.svd(V.T)
    null_tol = 1e-8 * max(float(sv[0]), 1.0)
    if float(sv[-1]) > null_tol:
        raise ValueError('The boundary equation x*(B1 + R*A2) = 0 has no nontrivial solution '
                         '(smallest singular value %g against tolerance %g), so the process is '
                         'not positive recurrent (Corollary 8(ii) of Bean and Nielsen, 2010).'
                         % (float(sv[-1]), null_tol))
    if m > 1 and float(sv[-2]) <= null_tol:
        raise ValueError('The boundary equation x*(B1 + R*A2) = 0 has a solution space of '
                         'dimension greater than one, the equilibrium vector is not unique.')
    pihat0 = Wh[-1, :].reshape(1, m)
    den = float((pihat0 @ e).item())
    if abs(den) < 1e-12 * np.max(np.abs(pihat0)):
        raise ValueError('The boundary vector cannot be normalised, x*e = 0.')
    pihat0 = pihat0 / den

    # Step 5: level-0 vector.
    ImRinv = linalg.inv(I - R)
    K = 1.0 / float((pihat0 @ ImRinv @ e).item())
    pi0 = K * pihat0

    # Consistency of the supplied boundary up-block: the level-0 balance
    # equation pi0*B0 + pi1*A1 + pi2*A2 = 0 must hold with pi_n = pi0*R^n.
    bal = pi0 @ B0 + pi0 @ R @ A1 + pi0 @ R @ R @ A2
    bal_norm = float(np.max(np.abs(bal)))
    if bal_norm > 1e-8 * block_scale * max(float(np.max(np.abs(pi0))), 1.0):
        raise ValueError('The boundary block B0 is inconsistent with the repeating blocks: '
                         '||pi0*B0 + pi1*A1 + pi2*A2||_inf = %g. The level-0 balance equation '
                         'of Theorem 7 requires pi0*(B0-A0) = 0.' % bal_norm)

    # Step 6: level vectors and marginal level distribution.
    pqueue = np.zeros((numLevels + 1, m))
    pin = pi0.copy()
    for n in range(numLevels + 1):
        pqueue[n, :] = pin[0, :]
        pin = pin @ R
    levelProb = pqueue @ e
    levelProb = levelProb.reshape(-1)

    # Exact mean queue length, sum_n n*pi0*R^n*e = pi0*R*inv(I-R)^2*e.
    QN = float((pi0 @ R @ ImRinv @ ImRinv @ e).item())

    return QbdRapResult(levelProb=levelProb, QN=QN, R=R, G=G, U=U, spr=spr, pqueue=pqueue,
                        pi0=pi0)


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
        from line_solver.lib.thirdparty.butools.ph.baseph import AcyclicPHFromMeansAndSCVs
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

    # see _kb/03-api-layer.md for rationale
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

    # see _kb/03-api-layer.md for rationale
    QN = 0.0
    j = n
    ni = 0
    while j + n <= len(pn):
        ni += 1
        QN += ni * np.sum(pn[j:j + n])
        j += n

    return QN


def _coxian_phase_subgen(rate: float, scv: float) -> np.ndarray:
    """
    Sub-generator of the canonical Coxian form with the given RATE and SCV,
    entered at phase 1.

    The four branches of ``Coxian.fitMeanAndSCV`` with CoarseTol 1e-3, written
    out rather than routed through the class so that the four codebases build
    the SAME phase: the entry vector has to be ``[1 0 ... 0]`` for every SCV,
    because :func:`qbd_setupdelayoff_closed` overloads the phase index by level
    and an arrival to an off server must enter the setup at phase 1.

    THIS IS A HAND COPY OF THE REFERENCE AND MUST TRACK IT. MATLAB and the JAR
    call the real fitter; native python has no ``Coxian.fit_mean_and_scv`` to
    call, so these branches are transcribed from
    ``matlab/src/lang/processes/Coxian.m`` and are only correct while that
    function is unchanged. Anything that edits the MATLAB fitter has to edit
    this too, and a divergence here is silent: it changes the phase, not the
    shape of the answer. The uncapped ``n = ceil(1/scv)`` of the low-SCV branch
    is the reference's own, deliberately reproduced rather than bounded -- a cap
    here alone would be the divergence.

    An exponential phase is built from the RATE directly rather than round-tripped
    through its mean: an Immediate setup has rate 1e8, whose mean is exactly
    FineTol, and the round trip turns it into an infinite rate.

    Args:
        rate: reciprocal of the phase mean
        scv: squared coefficient of variation of the phase

    Returns:
        The (n x n) sub-generator, n the number of Coxian phases.
    """
    if rate <= 0:
        raise ValueError("_coxian_phase_subgen: the rate must be positive")
    if scv <= 0:
        raise ValueError("_coxian_phase_subgen: the SCV must be positive")
    if scv == 1.0:
        return np.array([[-rate]], dtype=float)

    tol = 1e-3
    mean = 1.0 / rate
    if 1.0 - tol <= scv <= 1.0 + tol:
        mu = [1.0 / mean]
        phi = [1.0]
    elif 0.5 + tol < scv < 1.0 - tol:
        s = np.sqrt(1.0 + 2.0 * (scv - 1.0))
        mu = [2.0 / mean / (1.0 + s), 2.0 / mean / (1.0 - s)]
        phi = [0.0, 1.0]
    elif scv <= 0.5 + tol:
        n = int(np.ceil(1.0 / scv))
        lam = n / mean
        mu = [lam] * n
        phi = [0.0] * n
    else:
        mu1 = 2.0 / mean
        mu2 = mu1 / (2.0 * scv)
        mu = [mu1, mu2]
        phi = [1.0 - mu2 / mu1, 1.0]
    phi[-1] = 1.0

    n = len(mu)
    D0 = np.zeros((n, n))
    for i in range(n):
        D0[i, i] = -mu[i]
        if i + 1 < n:
            D0[i, i + 1] = mu[i] * (1.0 - phi[i])
    return D0


def qbd_setupdelayoff_closed(
    N: int,
    Z: float,
    mu: float,
    alpharate: float,
    alphascv: float,
    betarate: float,
    betascv: float
) -> Tuple[float, float]:
    """
    Mean queue length and throughput of a FINITE-POPULATION queue with setup
    delay and delay-off.

    The closed twin of :func:`qbd_setupdelayoff`. The population N is finite and
    the complementary delay Z is what the customers not at this station are
    passing through, so the arrival rate is state dependent, lambda(n) =
    (N - n)/Z, and the level index is bounded by N. That makes the chain a
    LEVEL-DEPENDENT QBD over finitely many levels, i.e. a finite CTMC, and it is
    solved exactly rather than by a matrix-geometric tail.

    THE SEMANTICS ARE THE SIMULATOR'S, not the mean-value shortcut's. When the
    queue empties the server begins a delay-off period; an arrival DURING it
    finds the server still warm and resumes without setup (Solver_ssj's
    cancelDelayoff), and only an arrival after the delay-off has expired pays the
    setup. That is an M/M/1 with setup time AND close-down time, which in a
    closed network is what the per-instance cold-start race
    ``p_cold * E[setup] + S`` fails to be: that formula races the delay-off
    against the per-instance idle time and carries NO queueing term, so it
    describes a serverless instance pool rather than a single-server vacation
    queue, and it left the reported response time byte-identical across a
    tenfold change in the setup mean.

    The phase index is overloaded by level exactly as in the open twin: at level
    0 phase 1 is the OFF server and the remaining phases are the delay-off; above
    level 0 the phases are the setup and the last one is the busy server.

    Args:
        N: population of the closed chain, a non-negative integer
        Z: complementary delay, the mean time a customer spends away from this station
        mu: service rate of the station
        alpharate: rate of the setup phase
        alphascv: squared coefficient of variation of the setup phase
        betarate: rate of the delay-off phase
        betascv: squared coefficient of variation of the delay-off phase

    Returns:
        ``(QN, X)``, the mean number at the station and its throughput.

    References:
        Original MATLAB: matlab/src/api/mam/qbd_setupdelayoff_closed.m
    """
    N = int(round(N))
    if N <= 0 or mu <= 0:
        return 0.0, 0.0
    Z = max(float(Z), GlobalConstants.FineTol)

    Ta = _coxian_phase_subgen(alpharate, alphascv)
    na = Ta.shape[0]
    ta = -Ta.sum(axis=1)
    Tb = _coxian_phase_subgen(betarate, betascv)
    nb = Tb.shape[0]
    tb = -Tb.sum(axis=1)

    # Only the REACHABLE states are enumerated: the open twin pads every level to
    # na+nb phases and lets the unused ones sit at zero, which a finite chain
    # cannot do -- an unreachable row is an absorbing row and makes the
    # stationary solve singular.
    off = 0                                        # level 0, server off
    base = 1 + nb                                  # level 0 delay-off: 1..nb
    def doff(j):
        return 1 + j
    def setup(n, i):
        return base + (n - 1) * (na + 1) + i
    def busy(n):
        return base + (n - 1) * (na + 1) + na
    m = base + N * (na + 1)

    Q = np.zeros((m, m))
    def lam(n):
        return (N - n) / Z if n < N else 0.0

    if lam(0) > 0:
        Q[off, setup(1, 0)] += lam(0)
    for j in range(nb):
        for j2 in range(nb):
            if j2 != j:
                Q[doff(j), doff(j2)] += Tb[j, j2]
        Q[doff(j), off] += tb[j]
        # an arrival during the delay-off cancels it and resumes WITHOUT setup
        if lam(0) > 0:
            Q[doff(j), busy(1)] += lam(0)
    for n in range(1, N + 1):
        for i in range(na):
            for i2 in range(na):
                if i2 != i:
                    Q[setup(n, i), setup(n, i2)] += Ta[i, i2]
            Q[setup(n, i), busy(n)] += ta[i]
            # an arrival during the setup joins the queue and the setup carries
            # on in the SAME phase; the level rises, the phase does not move
            if n < N and lam(n) > 0:
                Q[setup(n, i), setup(n + 1, i)] += lam(n)
        if n < N and lam(n) > 0:
            Q[busy(n), busy(n + 1)] += lam(n)
        # a completion that empties the queue starts the delay-off at its phase 1
        Q[busy(n), busy(n - 1) if n - 1 >= 1 else doff(0)] += mu
    for i in range(m):
        Q[i, i] = -Q[i].sum()

    A = np.vstack([Q.T, np.ones(m)])
    b = np.zeros(m + 1)
    b[-1] = 1.0
    pi = np.linalg.lstsq(A, b, rcond=None)[0]
    pi = np.maximum(pi, 0.0)
    total = pi.sum()
    if total <= 0:
        return 0.0, 0.0
    pi = pi / total

    QN = 0.0
    for n in range(1, N + 1):
        QN += n * (float(np.sum([pi[setup(n, i)] for i in range(na)])) + pi[busy(n)])
    X = mu * float(np.sum([pi[busy(n)] for n in range(1, N + 1)]))
    return QN, X


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
    'qbd_setupdelayoff_closed',
]
