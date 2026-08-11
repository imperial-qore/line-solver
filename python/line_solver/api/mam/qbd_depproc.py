"""
QBD departure process analysis via ETAQA truncation.

Constructs MAP approximations of departure processes from MAP/MAP/1 queues
using ETAQA-based truncation of the underlying QBD (Quasi-Birth-Death) process.

Functions:
    qbd_depproc_etaqa: Departure MAP for MAP/MAP/1-FCFS
    qbd_depproc_etaqa_ps: Departure MAP for MAP/MAP/1-PS (processor sharing)
    qbd_depproc_jointmom: Joint moments of consecutive inter-departure times
"""

import numpy as np
from numpy.linalg import inv
from scipy.special import factorial

from line_solver.lib.thirdparty.smc import qbd_cr, qbd_pi
from .map_analysis import map_lambda, _krons


def qbd_depproc_etaqa(MAPa, MAPs, n):
    """
    Construct MAP departure process for MAP/MAP/1-FCFS via ETAQA truncation.

    Builds a Markovian Arrival Process (MAP) approximation of the departure
    process from a single-server FCFS queue with MAP arrival process and MAP
    service process. The QBD is truncated at level n using ETAQA boundary
    corrections.

    Args:
        MAPa: Arrival MAP as list [D0, D1] of numpy arrays.
        MAPs: Service MAP as list [D0, D1] of numpy arrays.
        n: Truncation level (number of QBD levels beyond level 0).

    Returns:
        List [D0, D1] representing the departure MAP, where D0 and D1 are
        numpy arrays of size (n+1)*lvlsz x (n+1)*lvlsz.
    """
    na = MAPa[0].shape[0]
    ns = MAPs[0].shape[0]
    lvlsz = ns * na

    # QBD transition blocks
    F = np.kron(MAPa[1], np.eye(ns))          # forward (arrivals)
    L = _krons(MAPa[0], MAPs[0])              # local (Kronecker sum)
    B = np.kron(np.eye(na), MAPs[1])           # backward (service completions)
    L0 = np.kron(MAPa[0], np.eye(ns))          # level-0 local

    # Solve QBD for rate matrix R and G
    cr_result = qbd_cr(B, L, F)
    R = cr_result['R']
    G = inv(-L - R @ B) @ B

    # ETAQA boundary corrections
    Lhat = F + L
    Bbar = B + F @ G
    Bhat = F @ G

    # Build D0 (hidden transitions)
    vn1 = np.ones(n - 1)
    vn = np.zeros(n)
    vn[-1] = 1.0

    # Interior block: levels 1..n
    D0 = np.kron(np.diag(np.concatenate(([1.0], vn1))), L) + \
         np.kron(np.diag(vn1, 1), F)

    # Prepend rows for level 0
    D0 = np.vstack([np.zeros((L0.shape[0], D0.shape[1])), D0])

    # Prepend columns for backward transitions from level 1 to level 0
    D0 = np.hstack([np.zeros((D0.shape[0], B.shape[1])), D0])

    # Fill level-0 block: [L0, F]
    D0[:L0.shape[0], :L0.shape[1] + F.shape[1]] = np.hstack([L0, F])

    # Apply ETAQA correction at level n-1
    row_start = (n - 1) * lvlsz
    col_start = (n - 1) * lvlsz
    D0[row_start:row_start + lvlsz, col_start:col_start + lvlsz] = Lhat

    # Build D1 (observable transitions = departures)
    # ETAQA boundary: last level uses Bbar and Bhat
    D1 = np.kron(np.diag(vn, -1), Bbar) + \
         np.kron(np.diag(np.concatenate(([0.0], vn))), Bhat)

    # Interior levels use plain B
    D1_interior = np.kron(np.diag(vn1, -1), B)
    D1[:n * lvlsz, :n * lvlsz] = D1_interior

    return [D0, D1]


def qbd_depproc_etaqa_ps(MAPa, MAPs, n):
    """
    Construct MAP departure process for MAP/MAP/1-PS via ETAQA truncation.

    Same as qbd_depproc_etaqa but for Processor Sharing (PS) discipline.
    Under PS, when j customers are present, each is served at rate mu/j.
    This leads to rate-dependent service completion splitting in D0 and D1.

    Args:
        MAPa: Arrival MAP as list [D0, D1] of numpy arrays.
        MAPs: Service MAP as list [D0, D1] of numpy arrays.
        n: Truncation level (number of QBD levels beyond level 0).

    Returns:
        List [D0, D1] representing the departure MAP, where D0 and D1 are
        numpy arrays of size (n+1)*lvlsz x (n+1)*lvlsz.
    """
    na = MAPa[0].shape[0]
    ns = MAPs[0].shape[0]
    lvlsz = ns * na

    # QBD transition blocks
    F = np.kron(MAPa[1], np.eye(ns))
    L = _krons(MAPa[0], MAPs[0])
    B = np.kron(np.eye(na), MAPs[1])
    L0 = np.kron(MAPa[0], np.eye(ns))

    # Solve QBD for rate matrix R and G
    cr_result = qbd_cr(B, L, F)
    R = cr_result['R']
    G = inv(-L - R @ B) @ B

    # ETAQA boundary corrections
    Lhat = F + L
    Bbar = B + F @ G
    Bhat = F @ G

    # Build D0 (hidden transitions) - same structure as FCFS initially
    vn1 = np.ones(n - 1)
    vn = np.zeros(n)
    vn[-1] = 1.0

    D0 = np.kron(np.diag(np.concatenate(([1.0], vn1))), L) + \
         np.kron(np.diag(vn1, 1), F)

    D0 = np.vstack([np.zeros((L0.shape[0], D0.shape[1])), D0])
    D0 = np.hstack([np.zeros((D0.shape[0], B.shape[1])), D0])
    D0[:L0.shape[0], :L0.shape[1] + F.shape[1]] = np.hstack([L0, F])

    # Apply ETAQA correction at level n-1
    row_start = (n - 1) * lvlsz
    col_start = (n - 1) * lvlsz
    D0[row_start:row_start + lvlsz, col_start:col_start + lvlsz] = Lhat

    # PS-specific: add boundary terms to D0
    # D0 += kron(diag(vn,-1), Bbar) + kron(diag([0,vn]), Bhat)*(n-1)/n
    D0 += np.kron(np.diag(vn, -1), Bbar) + \
          np.kron(np.diag(np.concatenate(([0.0], vn))), Bhat) * (n - 1) / n

    # Build D1 (observable transitions = departures)
    # Boundary contribution: 1/n fraction at last level
    D1 = np.kron(np.diag(vn, -1), Bbar) + \
         np.kron(np.diag(np.concatenate(([0.0], vn))), Bhat) / n

    # PS rate-dependent splitting at interior levels
    for j in range(1, n):
        # see _kb/03-api-layer.md for rationale
        r_start = j * lvlsz
        r_end = (j + 1) * lvlsz
        c_start = (j - 1) * lvlsz
        c_end = j * lvlsz
        D0[r_start:r_end, c_start:c_end] = B * (1.0 - 1.0 / j)
        D1[r_start:r_end, c_start:c_end] = B * (1.0 / j)

    return [D0, D1]


def qbd_depproc_jointmom(MAPa, MAPs, iset):
    """
    Compute joint moments E[X_0^i * X_1^j] of consecutive inter-departure times.

    Uses matrix-analytic methods on the QBD representation of a MAP/MAP/1 queue
    to compute joint factorial moments of pairs of consecutive inter-departure
    times.

    Args:
        MAPa: Arrival MAP as list [D0, D1] of numpy arrays.
        MAPs: Service MAP as list [D0, D1] of numpy arrays.
        iset: Array of shape (K, 2) where each row [i, j] specifies the
              moment orders for E[X_0^i * X_1^j].

    Returns:
        1D numpy array of length K with the computed joint moments.
    """
    iset = np.atleast_2d(np.asarray(iset))

    na = MAPa[0].shape[0]
    ns = MAPs[0].shape[0]
    lvlsz = ns * na

    # QBD transition blocks
    F = np.kron(MAPa[1], np.eye(ns))
    L = _krons(MAPa[0], MAPs[0])
    B = np.kron(np.eye(na), MAPs[1])
    L0 = np.kron(MAPa[0], np.eye(ns))
    Z = np.zeros_like(F)

    # Solve QBD
    cr_result = qbd_cr(B, L, F)
    R = cr_result['R']

    # Stationary distribution
    pi = qbd_pi(B, L0, R)
    v0 = pi[:lvlsz].reshape(1, -1)

    # Service rate
    lambdaS = map_lambda(MAPs[0], MAPs[1])

    # Initial vectors for the 3-level representation
    I = np.eye(R.shape[0])
    v0D = (1.0 / lambdaS) * v0 @ R @ F
    v1D = (1.0 / lambdaS) * v0 @ R @ R @ F
    v2Dp = (1.0 / lambdaS) * v0 @ np.linalg.matrix_power(R, 3) @ inv(I - R) @ F

    z = np.hstack([v0D, v1D, v2Dp])
    z = z / z.sum()  # normalize to probability distribution

    # 3-level block matrices
    M0 = np.block([
        [L0, F, Z],
        [Z, L, F],
        [Z, Z, L + F]
    ])
    M1 = np.block([
        [Z, Z, Z],
        [B, Z, Z],
        [Z, B, Z]
    ])

    ones_vec = np.ones((M0.shape[0], 1))
    neg_M0_inv = inv(-M0)

    # Compute each joint moment
    JM = np.zeros(iset.shape[0])
    for k in range(iset.shape[0]):
        i_ord = int(iset[k, 0])
        j_ord = int(iset[k, 1])
        term = z @ (factorial(i_ord) *
                     np.linalg.matrix_power(neg_M0_inv, i_ord + 1) @
                     M1 @
                     (factorial(j_ord) *
                      np.linalg.matrix_power(neg_M0_inv, j_ord)) @
                     ones_vec)
        JM[k] = term.item()

    return JM
