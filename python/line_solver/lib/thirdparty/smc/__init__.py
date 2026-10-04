"""
SMC: Structured Markov Chain Solvers.

Native Python implementations of algorithms for Quasi-Birth-Death (QBD),
M/G/1, and G/I/M/1 type Markov chains.

Based on the SMCtools package by Benny Van Houdt.
"""

import warnings
import numpy as np
from numpy.typing import NDArray
from typing import Dict, Optional, Tuple, List
import scipy.linalg as la


def stat(A: np.ndarray) -> np.ndarray:
    """
    Compute the stationary distribution of a stochastic matrix.

    Args:
        A: Stochastic matrix (rows sum to 1)

    Returns:
        Stationary distribution as a row vector
    """
    A = np.asarray(A, dtype=float)
    n = A.shape[0]

    # Solve pi * (A - I) = 0 with sum(pi) = 1: matches MATLAB QMAM stat.m,
    # K = y / [A - I, e], i.e. the LEFT eigenvector for eigenvalue 1. The
    # block must use A (not A.T): lstsq below solves B^T x = y^T, which is
    # x^T [A - I | e] = [0, 1].
    B = np.hstack([A - np.eye(n), np.ones((n, 1))])
    y = np.zeros((1, n + 1))
    y[0, -1] = 1.0

    # Solve using least squares
    pi = la.lstsq(B.T, y.T)[0].T

    # No clipping: stat.m and Stat.java return the left eigenvector as solved.
    # The solution is nonnegative whenever A is a stochastic matrix, but with
    # RAP components A carries negative entries and so does its left
    # eigenvector; clipping there replaces the boundary vector of the QBD by a
    # different one and biases every level probability that follows from it.
    return pi.flatten()


def qbd_cr(A0: np.ndarray, A1: np.ndarray, A2: np.ndarray,
           max_num_it: int = 50, verbose: bool = False,
           mode: str = "Shift", rap_comp: bool = False
          ) -> Dict[str, np.ndarray]:
    """
    QBD_CR: Cyclic Reduction algorithm for QBD Markov chains.

    Computes the G, R, and U matrices for a QBD Markov chain.

    Args:
        A0: Downward transition matrix
        A1: Level transition matrix (or generator diagonal blocks)
        A2: Upward transition matrix
        max_num_it: Maximum number of iterations (default: 50)
        verbose: Print progress if True
        mode: "Shift" or "Basic" algorithm variant
        rap_comp: Set to True for RAP (Rational Arrival Process) components

    Returns:
        Dictionary with 'G', 'R', 'U' matrices
    """
    A0 = np.asarray(A0, dtype=float).copy()
    A1 = np.asarray(A1, dtype=float).copy()
    A2 = np.asarray(A2, dtype=float).copy()

    m = A1.shape[0]

    # Check if continuous time (negative diagonal)
    A1_diag = np.diag(A1)
    continues = False

    if not rap_comp:
        if np.min(A1_diag) < 0:
            continues = True
            lamb = np.max(-A1_diag)
            A0 = A0 / lamb
            A1 = A1 / lamb + np.eye(m)
            A2 = A2 / lamb
    else:
        if np.min(A1_diag) < 0:
            continues = True
            lamb = np.max(-A1_diag)
            A0 = A0 / lamb
            A1 = A1 / lamb + np.eye(m)
            A2 = A2 / lamb

    # Check whether G is known explicitly
    result = qbd_eg(A0, A1, A2, verbose)
    G = result.get('G')

    if G is not None:
        # Explicit solution found - compute R, U if not already done and return
        R = result.get('R')
        U_eg = result.get('U')
        if R is None:
            R = A2 @ la.inv(np.eye(m) - (A1 + A2 @ G))
        if U_eg is None:
            U_eg = A1 + R @ A0
        if continues:
            U_eg = lamb * (U_eg - np.eye(m))
        return {'G': G, 'R': R, 'U': U_eg}

    # Compute drift
    theta = stat(A0 + A1 + A2)
    drift = theta @ np.sum(A0, axis=1) - theta @ np.sum(A2, axis=1)

    A2_old = A2.copy()
    A0_old = A0.copy()
    uT = np.ones((1, m)) / m

    if mode == "Shift":
        if drift < 0:
            A2 = A2 - np.ones((m, 1)) @ (theta @ A2).reshape(1, -1)
            A1 = A1 + np.ones((m, 1)) @ (theta @ A0).reshape(1, -1)
        else:
            A0 = A0 - np.sum(A0, axis=1, keepdims=True) @ uT
            A1 = A1 + np.sum(A2, axis=1, keepdims=True) @ uT

    A = A1.copy()
    B = A2.copy()
    C = A0.copy()
    Ahat = A.copy()

    check = 1.0
    numit = 0

    while check > 1e-14 and numit < max_num_it:
        Atemp = la.inv(np.eye(m) - A)
        BAtemp = B @ Atemp
        CAtemp = C @ Atemp

        Ahat = Ahat + BAtemp @ C
        A = A + BAtemp @ C + CAtemp @ B
        B = BAtemp @ B
        C = CAtemp @ C

        numit += 1
        check = min(la.norm(B, np.inf), la.norm(C, np.inf))

        if verbose:
            print(f"Check after {numit} iterations: {check}")

    if numit == max_num_it and check > 1e-14:
        print("Maximum Number of Iterations reached")

    G = la.inv(np.eye(m) - Ahat) @ A0

    if mode == "Shift":
        if drift < 0:
            A1 = A1 - np.ones((m, 1)) @ (theta @ A0).reshape(1, -1)
            A2 = A2_old.copy()
        else:
            G = G + np.ones((m, 1)) @ uT
            A1 = A1 - np.sum(A2, axis=1, keepdims=True) @ uT
            A0 = A0_old.copy()

    if verbose:
        res_norm = la.norm(G - A0 - (A1 + A2 @ G) @ G, np.inf)
        print(f"Final Residual Error for G: {res_norm}")

    R = A2 @ la.inv(np.eye(m) - A1 - A2 @ G)

    if verbose:
        res_norm = la.norm(R - A2 - R @ (A1 + R @ A0), np.inf)
        print(f"Final Residual Error for R: {res_norm}")

    U = A1 + R @ A0

    if verbose:
        res_norm = la.norm(U - A1 - A2 @ la.inv(np.eye(m) - U) @ A0, np.inf)
        print(f"Final Residual Error for U: {res_norm}")

    if continues:
        U = lamb * (U - np.eye(m))

    return {'G': G, 'R': R, 'U': U}


def qbd_caudal(A0: np.ndarray, A1: np.ndarray, A2: np.ndarray,
               dual: bool = False) -> float:
    """
    QBD_CAUDAL: Compute the spectral radius of R (caudal characteristic).

    Computes the dominant eigenvalue of R, the smallest nonnegative solution
    to R = A2 + R*A1 + R^2*A0, when the QBD is recurrent.

    Args:
        A0: Downward transition matrix
        A1: Level transition matrix
        A2: Upward transition matrix
        dual: If True, return the dominant eigenvalue of the Ramaswami dual

    Returns:
        eta: The caudal characteristic (spectral radius of R)
    """
    A0 = np.asarray(A0, dtype=float)
    A1 = np.asarray(A1, dtype=float)
    A2 = np.asarray(A2, dtype=float)

    if dual:
        A0, A2 = A2.copy(), A0.copy()

    eta_min = 0.0
    eta_max = 1.0
    eta = 0.5
    while eta_max - eta_min > 1e-15:
        new_eta = np.max(np.real(la.eigvals(A2 + A1 * eta + A0 * eta**2)))
        if new_eta > eta:
            eta_min = eta
        else:
            eta_max = eta
        eta = (eta_min + eta_max) / 2.0
    return eta


def qbd_eg(A0: np.ndarray, A1: np.ndarray, A2: np.ndarray,
           verbose: bool = False) -> Dict[str, np.ndarray]:
    """
    QBD_EG: Determines G directly if rank(A0)=1 or rank(A2)=1.

    For special cases where the downward or upward transition matrix has rank 1,
    the G and R matrices can be computed in closed form without iteration.

    Args:
        A0: Downward transition matrix
        A1: Level transition matrix
        A2: Upward transition matrix
        verbose: Print residual errors if True

    Returns:
        Dictionary with 'G', 'R', 'U' matrices (empty dict values are None
        if no explicit solution exists)
    """
    A0 = np.asarray(A0, dtype=float)
    A1 = np.asarray(A1, dtype=float)
    A2 = np.asarray(A2, dtype=float)

    m = A0.shape[0]
    G = None
    R = None
    U = None

    theta = stat(A0 + A1 + A2)
    drift = theta @ np.sum(A0, axis=1) - theta @ np.sum(A2, axis=1)

    if drift > 0:  # positive recurrent case
        if np.linalg.matrix_rank(A0) == 1:
            # A0 = alpha * beta
            row_sums = np.sum(A0, axis=1)
            temp = -1
            for i in range(m):
                if row_sums[i] > 0:
                    temp = i
                    break
            if temp >= 0:
                beta = A0[temp, :] / np.sum(A0[temp, :])
                G = np.ones((m, 1)) @ beta.reshape(1, -1)
                R = A2 @ la.inv(np.eye(m) - (A1 + A2 @ G))
        elif np.linalg.matrix_rank(A2) == 1:
            eta = qbd_caudal(A0, A1, A2)
            R = A2 @ la.inv(np.eye(m) - A1 - eta * A0)
            G = la.inv(np.eye(m) - (A1 + R @ A0)) @ A0
    elif drift < 0:  # transient case
        if np.linalg.matrix_rank(A2) == 1:
            alpha = A2 @ np.ones((m, 1))
            R = (alpha @ theta.reshape(1, -1)) / (theta @ alpha)
            G = la.inv(np.eye(m) - (A1 + R @ A0)) @ A0
        elif np.linalg.matrix_rank(A0) == 1:
            # Use the dual (time-reversed) chain
            theta_inv = 1.0 / theta
            A0hat = np.diag(theta_inv) @ A2.T @ np.diag(theta)
            A1hat = np.diag(theta_inv) @ A1.T @ np.diag(theta)
            A2hat = np.diag(theta_inv) @ A0.T @ np.diag(theta)
            etahat = qbd_caudal(A0hat, A1hat, A2hat)
            Rhat = A2hat @ la.inv(np.eye(m) - A1hat - etahat * A0hat)
            G = np.diag(theta_inv) @ Rhat.T @ np.diag(theta)
            R = A2 @ la.inv(np.eye(m) - (A1 + A2 @ G))

    if R is not None:
        U = A1 + R @ A0

    if verbose:
        if G is not None:
            res_norm = la.norm(G - A0 - (A1 + A2 @ G) @ G, np.inf)
            print(f"Final Residual Error for G: {res_norm}")
        if R is not None:
            res_norm = la.norm(R - A2 - R @ (A1 + R @ A0), np.inf)
            print(f"Final Residual Error for R: {res_norm}")
        if U is not None:
            res_norm = la.norm(U - A1 - A2 @ la.inv(np.eye(m) - U) @ A0, np.inf)
            print(f"Final Residual Error for U: {res_norm}")

    return {'G': G, 'R': R, 'U': U}


def qbd_pi(B0: np.ndarray, B1: np.ndarray, R: np.ndarray,
           max_num_comp: int = 10000, verbose: int = 0,
           boundary: Optional[np.ndarray] = None,
           rap_comp: bool = False) -> np.ndarray:
    """
    QBD_pi: Stationary distribution of a QBD Markov chain.

    Computes the stationary vector for discrete or continuous time QBD
    with transition/rate matrix:

           B1  A2  0   0   0  ...
           B0  A1  A2  0   0  ...
       =   0   A0  A1  A2  0  ...
           0   0   A0  A1  A2 ...
           ...

    Args:
        B0: Boundary transition matrix (level 1 -> 0)
        B1: Boundary transition matrix (level 0)
        R: Rate matrix (minimal nonnegative solution to R = A2 + R*A1 + R^2*A0)
        max_num_comp: Maximum number of components (default: 10000)
        verbose: Print progress every verbose steps (0 = no output)
        boundary: More general boundary structure (optional)
        rap_comp: Set to True for RAP components

    Returns:
        Stationary distribution as concatenated row vectors [pi_0, pi_1, pi_2, ...]
    """
    B0 = np.asarray(B0, dtype=float).copy()
    B1 = np.asarray(B1, dtype=float).copy()
    R = np.asarray(R, dtype=float)

    m = R.shape[0]

    # Check if continuous time
    B1_diag = np.diag(B1)
    if np.min(B1_diag) < 0 or rap_comp:
        lamb = -np.min(B1_diag)
        B1 = B1 / lamb + np.eye(m)
        B0 = B0 / lamb

    # Check spectral radius of R. Positive recurrence is sp(R) < 1. For a
    # Markovian QBD, R is nonnegative and sp(R) < 1 is equivalent to
    # (I-R)^-1 >= 0, which is the test QBD_pi.m performs; with RAP components R
    # may have negative entries, (I-R)^-1 then legitimately has negative entries
    # too, and reading a single negative entry as non-recurrence rejects stable
    # RAP/RAP/1 queues. The eigenvalue test is the condition itself and agrees
    # with the entrywise one whenever R is nonnegative.
    temp = la.inv(np.eye(m) - R)
    if np.max(np.abs(np.linalg.eigvals(R))) >= 1.0:
        raise ValueError("The spectral radius of R is not below 1: QBD is not positive recurrent")

    # Compute pi_0
    pi0 = stat(B1 + R @ B0)
    normalizer = (pi0 @ temp @ np.ones((m, 1)))[0]
    pi0 = pi0 / normalizer

    # Build stationary distribution
    pi_components = [pi0.flatten()]
    sumpi = np.sum(pi0)
    numit = 1

    while sumpi < 1 - 1e-10 and numit < max_num_comp:
        pi_next = pi_components[-1] @ R
        pi_components.append(pi_next.flatten())
        numit += 1
        sumpi += np.sum(pi_next)

        if verbose > 0 and numit % verbose == 0:
            print(f"Accumulated mass after {numit} iterations: {sumpi}")

    if numit == max_num_comp:
        print(f"Maximum Number of Components {numit} reached")

    # Concatenate all components
    return np.concatenate(pi_components)


def qbd_lr(A0: np.ndarray, A1: np.ndarray, A2: np.ndarray,
           max_num_it: int = 50, verbose: bool = False,
           mode: str = "Shift") -> Dict[str, np.ndarray]:
    """
    QBD_LR: Logarithmic Reduction algorithm for QBD Markov chains
    [Latouche, Ramaswami].

    Computes the G, R, and U matrices for a QBD Markov chain using the
    Logarithmic Reduction method.

    Args:
        A0: Downward transition matrix
        A1: Level transition matrix
        A2: Upward transition matrix
        max_num_it: Maximum number of iterations (default: 50)
        verbose: Print progress if True
        mode: "Shift" or "Basic"

    Returns:
        Dictionary with 'G', 'R', 'U' matrices
    """
    A0 = np.asarray(A0, dtype=float).copy()
    A1 = np.asarray(A1, dtype=float).copy()
    A2 = np.asarray(A2, dtype=float).copy()

    m = A1.shape[0]

    # Convert to discrete time problem, if needed
    A1_diag = np.diag(A1)
    continues = False
    if np.min(A1_diag) < 0:
        continues = True
        lamb = np.max(-A1_diag)
        A0 = A0 / lamb
        A1 = A1 / lamb + np.eye(m)
        A2 = A2 / lamb

    # Check whether G is known explicitly
    result = qbd_eg(A0, A1, A2, verbose)
    G_eg = result.get('G')
    if G_eg is None or G_eg.size == 0:
        return result

    # Shift technique
    drift = 0.0
    if mode == "Shift":
        theta = stat(A0 + A1 + A2)
        drift = theta @ np.sum(A0, axis=1) - theta @ np.sum(A2, axis=1)
        if drift < 0:
            A2_old = A2.copy()
            A2 = A2 - np.ones((m, 1)) @ (theta @ A2).reshape(1, -1)
            A1 = A1 + np.ones((m, 1)) @ (theta @ A0).reshape(1, -1)
        else:
            uT = np.ones((1, m)) / m
            A0_old = A0.copy()
            A0 = A0 - np.sum(A0, axis=1, keepdims=True) @ uT
            A1 = A1 + np.sum(A2, axis=1, keepdims=True) @ uT

    # Start of Logarithmic Reduction (Basic)
    B2 = la.inv(np.eye(m) - A1)
    B0 = B2 @ A2
    B2 = B2 @ A0
    G = B2.copy()
    PI = B0.copy()

    check = 1.0
    numit = 0
    while check > 1e-14 and numit < max_num_it:
        A1star = B2 @ B0 + B0 @ B2
        A0star = B0 @ B0
        A2star = B2 @ B2
        B0 = la.inv(np.eye(m) - A1star)
        B2 = B0 @ A2star
        B0 = B0 @ A0star
        G = G + PI @ B2
        PI = PI @ B0
        check = min(la.norm(B0, np.inf), la.norm(B2, np.inf))
        numit += 1
        if verbose:
            print(f"Check after {numit} iterations: {check}")

    if numit == max_num_it and check > 1e-14:
        print("Maximum Number of Iterations reached")

    # Shift Technique restoration
    if mode == "Shift":
        if drift < 0:
            A1 = A1 - np.ones((m, 1)) @ (theta @ A0).reshape(1, -1)
            A2 = A2_old.copy()
        else:
            G = G + np.ones((m, 1)) @ uT
            A1 = A1 - np.sum(A2, axis=1, keepdims=True) @ uT
            A0 = A0_old.copy()

    if verbose:
        res_norm = la.norm(G - A0 - (A1 + A2 @ G) @ G, np.inf)
        print(f"Final Residual Error for G: {res_norm}")

    # Compute R
    R = A2 @ la.inv(np.eye(m) - (A1 + A2 @ G))
    if verbose:
        res_norm = la.norm(R - A2 - R @ (A1 + R @ A0), np.inf)
        print(f"Final Residual Error for R: {res_norm}")

    # Compute U
    U = A1 + R @ A0
    if verbose:
        res_norm = la.norm(U - A1 - A2 @ la.inv(np.eye(m) - U) @ A0, np.inf)
        print(f"Final Residual Error for U: {res_norm}")

    if continues:
        U = lamb * (U - np.eye(m))

    return {'G': G, 'R': R, 'U': U}


def mg1_g(A: np.ndarray, max_num_it: int = 100, verbose: bool = False
         ) -> Optional[np.ndarray]:
    """
    MG1_G: Compute G matrix for M/G/1-type Markov chain.

    Args:
        A: Block matrix [A0 A1 A2 ... A_max] with m rows and m*(max+1) columns
        max_num_it: Maximum number of iterations
        verbose: Print progress if True

    Returns:
        G matrix or None if computation fails
    """
    A = np.asarray(A, dtype=float)
    m = A.shape[0]
    dega = A.shape[1] // m - 1

    # Try explicit G computation first
    G = mg1_eg(A, verbose)
    if G is not None:
        return G

    # Fall back to functional iteration
    G = np.zeros((m, m))

    for numit in range(max_num_it):
        # G = A_max
        G_new = A[:, dega * m:(dega + 1) * m].copy()

        # G = A_{max-1} + G * G
        for i in range(dega - 1, -1, -1):
            G_new = A[:, i * m:(i + 1) * m] + G_new @ G

        if la.norm(G_new - G, np.inf) < 1e-14:
            if verbose:
                print(f"Converged after {numit + 1} iterations")
            return G_new

        G = G_new

    if verbose:
        print(f"Maximum iterations {max_num_it} reached")

    return G


def mg1_eg(A: np.ndarray, verbose: bool = False) -> Optional[np.ndarray]:
    """
    MG1_EG: Explicit G computation for M/G/1-type when rank(A0)=1.

    Args:
        A: Block matrix [A0 A1 A2 ... A_max]
        verbose: Print residual error if True

    Returns:
        G matrix if explicit solution exists, None otherwise
    """
    A = np.asarray(A, dtype=float)
    m = A.shape[0]
    dega = A.shape[1] // m - 1

    # Compute sumA and beta
    sumA = A[:, dega * m:(dega + 1) * m].copy()
    beta = np.sum(sumA, axis=1, keepdims=True)

    for i in range(dega - 1, 0, -1):
        sumA = sumA + A[:, i * m:(i + 1) * m]
        beta = beta + np.sum(sumA, axis=1, keepdims=True)

    sumA = sumA + A[:, 0:m]

    # Compute stationary distribution
    theta = stat(sumA)
    drift = (theta @ beta)[0]

    G = None

    if drift < 1:
        # Positive recurrent case
        A0 = A[:, 0:m]
        rank = np.linalg.matrix_rank(A0)

        if rank == 1:
            row_sums = np.sum(A0, axis=1)
            temp = -1
            for i in range(m):
                if row_sums[i] > 0:
                    temp = i
                    break

            if temp >= 0:
                row_sum = np.sum(A0[temp, :])
                beta_vec = A0[temp, :] / row_sum
                G = np.ones((m, 1)) @ beta_vec.reshape(1, -1)

    elif drift > 1:
        # Transient chain: G through the Ramaswami dual and its caudal value.
        A0 = A[:, 0:m]
        if np.linalg.matrix_rank(A0) == 1:
            At = [np.diag(1.0 / theta) @ A[:, b*m:(b+1)*m].T @ np.diag(theta)
                  for b in range(dega + 1)]
            etahat = gim1_caudal(At)
            temp_m = At[dega].copy()
            for i in range(dega - 1, 0, -1):
                temp_m = temp_m * etahat + At[i]
            M = At[0] @ la.inv(np.eye(m) - temp_m)
            G = np.diag(1.0 / theta) @ M.T @ np.diag(theta)

    if verbose and G is not None:
        # Verify solution
        Gcheck = A[:, dega * m:(dega + 1) * m].copy()
        for j in range(dega - 1, -1, -1):
            Gcheck = A[:, j * m:(j + 1) * m] + Gcheck @ G
        res_norm = la.norm(G - Gcheck, np.inf)
        print(f"Final Residual Error for G: {res_norm}")

    return G


def gim1_r(A: np.ndarray, max_num_it: int = 100, verbose: bool = False
          ) -> Optional[np.ndarray]:
    """
    GIM1_R: Compute R matrix for G/I/M/1-type Markov chain.

    Args:
        A: Block matrix [A_-max ... A_-1 A0] with m rows
        max_num_it: Maximum number of iterations
        verbose: Print progress if True

    Returns:
        R matrix or None if computation fails
    """
    A = np.asarray(A, dtype=float)
    m = A.shape[0]
    dega = A.shape[1] // m - 1

    # Functional iteration for R
    R = np.zeros((m, m))

    for numit in range(max_num_it):
        # R = A_{-max}
        R_new = A[:, 0:m].copy()

        # R = A_{-max+1} + R * R, etc.
        R_power = R.copy()
        for i in range(1, dega + 1):
            R_new = R_new + R_power @ A[:, i * m:(i + 1) * m]
            R_power = R_power @ R

        if la.norm(R_new - R, np.inf) < 1e-14:
            if verbose:
                print(f"Converged after {numit + 1} iterations")
            return R_new

        R = R_new

    if verbose:
        print(f"Maximum iterations {max_num_it} reached")

    return R


def gim1_pi(B0: np.ndarray, B1: np.ndarray, G: np.ndarray,
            max_num_comp: int = 500, verbose: int = 0) -> np.ndarray:
    """
    GIM1_pi: Stationary distribution for G/I/M/1-type Markov chain.

    Args:
        B0: Boundary block
        B1: Boundary diagonal block
        G: G matrix from gim1 solver
        max_num_comp: Maximum number of components
        verbose: Print progress every verbose steps

    Returns:
        Stationary distribution
    """
    # Similar structure to qbd_pi
    B0 = np.asarray(B0, dtype=float).copy()
    B1 = np.asarray(B1, dtype=float).copy()
    G = np.asarray(G, dtype=float)

    m = G.shape[0]

    # Compute pi_0
    pi0 = stat(B1 + B0 @ G)

    # Normalize
    temp = la.inv(np.eye(m) - G)
    normalizer = (pi0 @ temp @ np.ones((m, 1)))[0]
    pi0 = pi0 / normalizer

    # Build stationary distribution
    pi_components = [pi0.flatten()]
    sumpi = np.sum(pi0)
    numit = 1

    while sumpi < 1 - 1e-10 and numit < max_num_comp:
        pi_next = G @ pi_components[-1].reshape(-1, 1)
        pi_components.append(pi_next.flatten())
        numit += 1
        sumpi += np.sum(pi_next)

        if verbose > 0 and numit % verbose == 0:
            print(f"Accumulated mass after {numit} iterations: {sumpi}")

    if numit == max_num_comp:
        print(f"Maximum Number of Components {numit} reached")

    return np.concatenate(pi_components)


# ---------------------------------------------------------------------------
# SMCSolver engines behind ETAQA: MG1_Shifts, MG1_CR, MG1_FI, MG1_Decay,
# GIM1_Caudal and the dual-based GIM1_R, ported from
# matlab/lib/thirdparty/MG1files. These are what MG1_G_ETAQA and GIM1_R_ETAQA
# call; before they existed here the two ETAQA entry points used a plain
# functional iteration, which converges to a different accuracy and made the
# Python answers disagree with MATLAB.
# ---------------------------------------------------------------------------


def _blocks_of(A: np.ndarray, m: int) -> List[np.ndarray]:
    """Splits the wide [A0 A1 ... Amax] into its m x m blocks."""
    A = np.asarray(A, dtype=float)
    if m == 0 or A.shape[1] % m != 0:
        raise ValueError("smc: the block sequence has an incorrect number of columns")
    return [A[:, b*m:(b+1)*m].copy() for b in range(A.shape[1] // m)]


def _hcat(blocks: List[np.ndarray]) -> np.ndarray:
    """Re-assembles a block sequence into the wide [A0 A1 ... Amax]."""
    return np.hstack(blocks)


def _vblocks_of(A: np.ndarray, r: int) -> List[np.ndarray]:
    """Splits a vertical stack into its blocks of r rows."""
    A = np.asarray(A, dtype=float)
    if r == 0 or A.shape[0] % r != 0:
        raise ValueError("smc: the stacked block sequence has an incorrect number of rows")
    return [A[b*r:(b+1)*r, :].copy() for b in range(A.shape[0] // r)]


def mg1_drift(blocks: List[np.ndarray]) -> Tuple[float, np.ndarray]:
    """
    drift = theta * beta with beta = (Amax)e + (Amax+Amax-1)e + ..., the
    expected level increment per transition, and theta = stat(sum_i A_i).
    """
    dega = len(blocks) - 1
    sumA = blocks[dega].copy()
    beta = np.sum(sumA, axis=1)
    for i in range(dega - 1, 0, -1):
        sumA = sumA + blocks[i]
        beta = beta + np.sum(sumA, axis=1)
    sumA = sumA + blocks[0]
    theta = stat(sumA)
    return float(theta @ beta), theta


def _poly_at(blocks: List[np.ndarray], z: float) -> np.ndarray:
    """A(z) = A0 + A1 z + ... + Amax z^max, by Horner as the reference writes it."""
    temp = blocks[-1].copy()
    for i in range(len(blocks) - 2, -1, -1):
        temp = temp * z + blocks[i]
    return temp


def _max_eig(M: np.ndarray) -> complex:
    """
    max(eig(M)) with MATLAB's semantics on a complex spectrum: largest modulus,
    ties broken by the larger phase angle. For the nonnegative A(z) both callers
    evaluate this is the Perron-Frobenius eigenvalue and is real.
    """
    ev = np.linalg.eigvals(M)
    best = ev[0]
    for z in ev:
        if abs(z) > abs(best) or (abs(z) == abs(best) and np.angle(z) > np.angle(best)):
            best = z
    return best


def mg1_decay(blocks: List[np.ndarray]) -> float:
    """
    Decay rate of a recurrent M/G/1-type chain: the unique z > 1 with
    PF(A(z)) = z. Port of MG1_Decay.m (the eigenvector output is not returned,
    since its only caller is the 'tau' shift, which is refused below).
    """
    eta, new_eta = 1.0, 0.0
    while new_eta - eta < 0:
        eta += 1.0
        new_eta = _max_eig(_poly_at(blocks, eta)).real
    eta_min, eta_max = eta - 1.0, eta
    eta = eta_min + 0.5
    while eta_max - eta_min > 1e-15:
        new_eta = _max_eig(_poly_at(blocks, eta)).real
        if new_eta < eta:
            eta_min = eta
        else:
            eta_max = eta
        eta = (eta_min + eta_max) / 2.0
    return eta


def gim1_caudal(blocks: List[np.ndarray]) -> float:
    """
    Caudal characteristic of a GI/M/1-type chain: the spectral radius of R, the
    unique z in (0,1) with PF(A(z)) = z. Port of GIM1_Caudal.m.
    """
    eta_min, eta_max, eta = 0.0, 1.0, 0.5
    while eta_max - eta_min > 1e-15:
        new_eta = _max_eig(_poly_at(blocks, eta)).real
        if new_eta > eta:
            eta_min = eta
        else:
            eta_max = eta
        eta = (eta_min + eta_max) / 2.0
    return eta


def mg1_shifts(blocks: List[np.ndarray], shift_type: str = 'one'
               ) -> Tuple[List[np.ndarray], float]:
    """
    Shift technique for an M/G/1-type sequence. Port of MG1_Shifts.m,
    ShiftType 'one', which is the default and the only type ETAQA uses: for a
    positive recurrent chain the eigenvalue 1 of A(z) is shifted to zero, so
    cyclic reduction converges on the second largest root instead of stalling
    on the unit one; the shift is undone on G afterwards.

    'tau' and 'dbl', and 'one' at drift > 1, are REFUSED rather than
    transcribed: their reference branches write rowhatA(1,maxd*i:end) with `i`
    undefined at that point, so the line does not compute the last block it
    needs and cannot run in MATLAB either.
    """
    if shift_type != 'one':
        raise NotImplementedError(
            "MG1_Shifts: ShiftType '%s' is not ported; its reference branch writes "
            "rowhatA(1,maxd*i:end) with `i` undefined at that point. ETAQA uses 'one'"
            % shift_type)
    A = [b.copy() for b in blocks]
    m = A[0].shape[0]
    drift, _ = mg1_drift(A)
    if not drift < 1.0:
        raise NotImplementedError(
            "MG1_Shifts: the drift > 1 branch of ShiftType 'one' is not ported; it shifts "
            "one to infinity through the same defective rowhatA(1,maxd*i:end) line")

    A[1] = A[1] - np.eye(m)
    col = np.zeros(m)
    hatA = []
    for b in range(len(A)):
        col = col + np.sum(A[b], axis=1)
        hatA.append(A[b] - np.outer(col, np.ones(m)) / m)
    hatA[1] = hatA[1] + np.eye(m)
    return hatA, drift


def _cr_block_fft(blocks: List[np.ndarray], n: int, use: int) -> np.ndarray:
    """Blockwise DFT along the block index, MATLAB's fft over the sequence."""
    m = blocks[0].shape[0]
    take = min(use, len(blocks))
    arr = np.zeros((take, m, m))
    for k in range(take):
        arr[k] = blocks[k]
    return np.fft.fft(arr, n=n, axis=0)


def _cr_tail_norm(blocks: List[np.ndarray]) -> float:
    """
    The tail norm the reference measures, over blocks deg/2 .. deg-1. MATLAB's
    `for i=deg/2:deg-1` starts at a HALF-INTEGER when deg is odd and so runs
    from ceil(deg/2); reproduced, or a degree-1 sequence would be measured here
    where the reference measures nothing.
    """
    deg = len(blocks)
    start = (deg + 1) // 2
    best = 0.0
    for i in range(start, deg):
        best = max(best, np.linalg.norm(blocks[i], np.inf))
    return best


def mg1_cr(A: np.ndarray, mode: str = 'ShiftPWCR', shift_type: str = 'one',
           max_num_it: int = 50, max_num_root: int = 2048,
           epsilon: float = 1e-16) -> np.ndarray:
    """
    Cyclic reduction for M/G/1-type Markov chains [Bini, Meini]. Port of
    MG1_CR.m, default mode 'ShiftPWCR' with ShiftType 'one'.

    One step of cyclic reduction eliminates every odd level and leaves a chain
    of the same shape on the even ones, so the level distance halves per
    iteration and the iteration converges quadratically. The composition is
    done POINT-WISE at the (nj+1)-th roots of unity and interpolated back, with
    the number of roots doubling until the interpolated tail falls below
    (nj+1) eps. Everything runs on the TRANSPOSED blocks, as the reference does
    after `D=D'`, and the final G is transposed back.
    """
    A = np.asarray(A, dtype=float)
    m = A.shape[0]
    if mode not in ('ShiftPWCR', 'PWCR'):
        raise NotImplementedError("MG1_CR: Mode '%s' is not supported" % mode)

    Geg = mg1_eg(A)
    if Geg is not None:
        return Geg

    blocks = _blocks_of(A, m)
    drift = 0.0
    if mode == 'ShiftPWCR':
        blocks, drift = mg1_shifts(blocks, shift_type)

    maxd = len(blocks) - 1
    if maxd == 0:
        raise ValueError("MG1_CR: the sequence needs at least two blocks")
    target = 1
    while target < maxd:
        target <<= 1
    if target == maxd:
        target <<= 1
    target += 1

    D = [np.zeros((m, m)) for _ in range(target)]
    for b in range(maxd + 1):
        D[b] = blocks[b].T.copy()

    Aeven = D[0::2]
    Aodd = D[1::2]
    Ahatodd = Aeven[1:] + [D[-1]]
    Ahateven = list(Aodd)

    Rj = sum(D[1:])
    Rj = D[0] @ la.inv(np.eye(m) - Rj)

    G = np.zeros((m, m))
    Anew: List[np.ndarray] = []
    Ahatnew: List[np.ndarray] = []
    numit = 0
    while numit < max_num_it:
        numit += 1
        nj = len(Aodd) - 1

        if nj > 0:
            n = nj + 1
            Anew, Ahatnew = _cr_step(Aodd, Aeven, Ahatodd, Ahateven, n, n)
        else:
            temp = Aeven[0] @ la.inv(np.eye(m) - Aodd[0])
            Ahatnew = [Ahateven[0] + temp @ Ahatodd[0]]
            Anew = [temp @ Aeven[0], Aodd[0]]

        nAnew = _cr_tail_norm(Anew)
        nAhatnew = _cr_tail_norm(Ahatnew)

        while ((nAnew > (nj + 1) * epsilon or nAhatnew > (nj + 1) * epsilon)
               and nj + 1 < max_num_root):
            nj = 2 * (nj + 1) - 1
            n = nj + 1
            stopv = min(n, len(Aodd))
            Anew, Ahatnew = _cr_step(Aodd, Aeven, Ahatodd, Ahateven, n, stopv)
            nAnew = _cr_tail_norm(Anew)
            nAhatnew = _cr_tail_norm(Ahatnew)

        if nj > 1:
            keep = (nj + 1) // 2
            Anew = Anew[:min(keep, len(Anew))]
            Ahatnew = Ahatnew[:min(keep, len(Ahatnew))]

        Aeven = Anew[0::2]
        Aodd = Anew[1::2]
        Ahateven = Ahatnew[0::2]
        Ahatodd = Ahatnew[1::2]

        if mode == 'PWCR':
            Rnewj = sum(Anew[1:]) if len(Anew) > 1 else np.zeros((m, m))
            Rnewj = Anew[0] @ la.inv(np.eye(m) - Rnewj)
            U = (np.eye(m) - Anew[0] @ la.inv(np.eye(m) - Anew[1])
                 if len(Anew) > 1 else np.eye(m))
            if (np.max(np.abs(Rj - Rnewj)) < epsilon
                    or np.max(np.sum(U, axis=0)) < epsilon):
                G = Ahatnew[0].copy()
                for i in range(1, len(Ahatnew)):
                    G = G + Rnewj @ Ahatnew[i]
                G = D[0] @ la.inv(np.eye(m) - G)
                break
            Rj = Rnewj
            tail_sum = sum(float(np.sum(b)) for b in Ahatnew[1:])
            sv = la.svdvals(Anew[0])[0] if Anew[0].size else 0.0
            V = np.eye(m) - D[0] @ la.inv(np.eye(m) - Ahatnew[0])
            if (sv < epsilon or tail_sum < epsilon
                    or np.max(np.sum(V, axis=0)) < epsilon):
                G = D[0] @ la.inv(np.eye(m) - Ahatnew[0])
                break
        else:
            Gold = G
            G = D[0] @ la.inv(np.eye(m) - Ahatnew[0])
            tail = 0.0 if len(Ahatnew) < 2 else max(
                np.linalg.norm(b, np.inf) for b in Ahatnew[1:])
            if np.linalg.norm(G - Gold, np.inf) < epsilon or tail < epsilon:
                break

    if numit == max_num_it and Ahatnew:
        G = D[0] @ la.inv(np.eye(m) - Ahatnew[0])

    G = G.T
    # Undo the shift: shifting one to zero removed the rank-one term e u^T.
    if mode == 'ShiftPWCR' and drift < 1.0:
        G = G + np.ones((m, m)) / m
    return G


def _cr_step(Aodd: List[np.ndarray], Aeven: List[np.ndarray],
             Ahatodd: List[np.ndarray], Ahateven: List[np.ndarray],
             n: int, use: int) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """One point-wise evaluation of (6.20) in Meini's thesis, at n roots."""
    m = Aodd[0].shape[0]
    T1 = _cr_block_fft(Aodd, n, use)
    T2 = _cr_block_fft(Aeven, n, use)
    T3 = _cr_block_fft(Ahatodd, n, use)
    T4 = _cr_block_fft(Ahateven, n, use)
    Ah = np.zeros((n, m, m), dtype=complex)
    An = np.zeros((n, m, m), dtype=complex)
    I = np.eye(m)
    for c in range(n):
        W = la.inv(I - T1[c])
        Ah[c] = T4[c] + T2[c] @ W @ T3[c]
        An[c] = np.exp(-c * 2j * np.pi / n) * T1[c] + T2[c] @ W @ T2[c]
    Ahatnew = np.real(np.fft.ifft(Ah, n=n, axis=0))
    Anew = np.real(np.fft.ifft(An, n=n, axis=0))
    return [Anew[k] for k in range(n)], [Ahatnew[k] for k in range(n)]


def mg1_fi(A: np.ndarray, mode: str = 'U-Based', shift_type: str = 'one',
           max_num_it: int = 10000, tol: float = 1e-14) -> np.ndarray:
    """
    Functional iterations for M/G/1-type Markov chains [Neuts]. Port of
    MG1_FI.m for Natural, Traditional and U-Based plus the Shift variants; the
    reference's NonZeroBlocks option is not exposed, since it only skips
    products over vanishing A_i and converges to the same G.

    'U-Based' is the default and the one GIM1_R(...,'FI') uses.
    """
    A = np.asarray(A, dtype=float)
    m = A.shape[0]
    maxd = A.shape[1] // m - 1

    Geg = mg1_eg(A)
    if Geg is not None:
        return Geg

    blocks = _blocks_of(A, m)
    drift = 0.0
    shifted = 'Shift' in mode
    if shifted:
        blocks, drift = mg1_shifts(blocks, shift_type)

    natural = 'Natural' in mode
    traditional = 'Traditional' in mode
    ubased = 'U-Based' in mode
    if not (natural or traditional or ubased):
        raise NotImplementedError("MG1_FI: Mode '%s' is not supported" % mode)

    G = np.zeros((m, m))
    check = 1.0
    numit = 0
    while check > tol and numit < max_num_it:
        Gold = G
        if natural:
            G = blocks[maxd].copy()
            for j in range(maxd - 1, -1, -1):
                G = blocks[j] + G @ Gold
        elif traditional:
            G = blocks[maxd].copy()
            for j in range(maxd - 1, 1, -1):
                G = blocks[j] + G @ Gold
            G = blocks[0] + G @ Gold @ Gold
            G = la.inv(np.eye(m) - blocks[1]) @ G
        else:
            G = blocks[maxd].copy()
            for j in range(maxd - 1, 0, -1):
                G = blocks[j] + G @ Gold
            G = la.inv(np.eye(m) - G) @ blocks[0]
        check = np.linalg.norm(G - Gold, np.inf)
        numit += 1

    if shifted and drift < 1.0:
        G = G + np.ones((m, m)) / m
    return G


def _solve_sylv_powers(Amat: np.ndarray, B: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Solve sum_{j=1}^N B_j Y A^{j-1} = C for Y.

    Port of solveSylvPowersDirectSum.m: the equation is linear in vec(Y) with
    coefficient sum_j kron((A^{j-1})', B_j), so it is assembled and solved
    directly. The reference's Schur variants only reorganize the SAME solve for
    speed on large blocks; on the block sizes reached here the direct sum is
    both exact and cheaper than a Schur reduction.
    """
    m = Amat.shape[0]
    n = B.shape[0]
    N = B.shape[1] // n
    avec = [np.eye(m), Amat.T.copy()]
    for i in range(2, N):
        avec.append(avec[i - 1] @ avec[1])
    Z = np.zeros((m * n, m * n))
    for j in range(N):
        Aj = avec[j] if j < len(avec) else np.linalg.matrix_power(Amat.T, j)
        Z = Z + np.kron(Aj, B[:, j * n:(j + 1) * n])
    y = np.linalg.solve(Z, C.reshape(m * n, order='F'))
    return y.reshape((n, m), order='F')


def mg1_ni(A: np.ndarray, shift_type: str = 'one', max_num_it: int = 50,
           epsilon: float = 1e-14) -> np.ndarray:
    """
    G of an M/G/1-type chain by NEWTON ITERATION. Port of MG1_NI.m
    (Perez, Telek, Van Houdt), Mode 'RealSchurShift' with the direct-sum
    Sylvester solve.

    Each step linearizes G = sum_i A_i G^i about the current iterate, which
    leaves a Sylvester equation with matrix powers, sum_j B_j Y G^{j-1} = C,
    solved by _solve_sylv_powers. Newton converges quadratically, so it needs
    far fewer steps than functional iteration, at one linear solve of size
    (m^2) per step.

    Args:
        A: block matrix [A0 A1 ... A_max], m rows
        shift_type: shift applied before the iteration and undone after; only
            'one' is available, as for every other solver in this module
        max_num_it: iteration cap
        epsilon: convergence tolerance on ||G - Gold||_inf

    Returns:
        The minimal nonnegative solution G.
    """
    A = np.asarray(A, dtype=float)
    m = A.shape[0]

    explicit = mg1_eg(A)
    if explicit is not None:
        return explicit

    blocks, drift = mg1_shifts(_blocks_of(A, m), shift_type)
    A = _hcat(blocks)
    N = A.shape[1] // m - 1

    G = np.zeros((m, m))
    check = 1.0
    numit = 0
    while check > epsilon and numit < max_num_it:
        Gold = G
        # B: the matrices premultiplying Y_k; C: the right-hand side
        B = np.zeros((m, (N + 1) * m))
        B[:, N * m:(N + 1) * m] = A[:, N * m:(N + 1) * m]
        for i in range(N - 1, -1, -1):
            B[:, i * m:(i + 1) * m] = A[:, i * m:(i + 1) * m] \
                + B[:, (i + 1) * m:(i + 2) * m] @ G
        C = G - B[:, 0:m]
        B = B[:, m:]
        B[:, 0:m] = B[:, 0:m] - np.eye(m)

        Y = _solve_sylv_powers(G, B, C)
        G = G + Y
        check = float(np.linalg.norm(G - Gold, np.inf))
        numit += 1

    if numit == max_num_it and check > epsilon:
        warnings.warn('MG1_NI: maximum number of iterations %d reached' % numit)

    if drift < 1.0:
        G = G + np.ones((m, m)) / m
    return G


def mg1_is(A: np.ndarray, mode: str = 'Schur', max_num_it: int = 50,
           epsilon: float = 1e-14) -> np.ndarray:
    """
    G of an M/G/1-type chain by the INVARIANT SUBSPACE method. Port of MG1_IS.m
    (Akar, Sohraby).

    F(z) = z - A(z) is mapped by the Cayley transform z = (1+s)/(1-s) onto a
    matrix polynomial H(s), whose companion pencil is deflated so that the
    stable invariant subspace of the resulting Z carries G. The subspace comes
    either from an ordered real Schur decomposition ('Schur', the default here
    since scipy provides ordqz/ordered Schur natively) or from the matrix sign
    iteration ('MSignStandard', 'MSignBalzer').

    Args:
        A: block matrix [A0 A1 ... A_max], m rows
        mode: 'Schur', 'MSignStandard' or 'MSignBalzer'
        max_num_it: iteration cap for the matrix sign modes
        epsilon: convergence tolerance for the matrix sign modes

    Returns:
        The minimal nonnegative solution G.
    """
    import scipy.linalg as sla

    A = np.asarray(A, dtype=float)
    m = A.shape[0]
    f = A.shape[1] // m - 1

    explicit = mg1_eg(A)
    if explicit is not None:
        return explicit

    orig = _blocks_of(A, m)
    drift, _theta = mg1_drift(orig)
    # MG1_IS measures the drift as theta*beta - 1, so it is negative exactly
    # when mg1_drift's ratio is below one; only its SIGN is used below.
    drift_sign = np.sign(drift - 1.0)

    # Step 1, F(z) = z - A(z)
    F = [-b.copy() for b in orig]
    F[1] = np.eye(m) + F[1]

    # Step 2, H(s) = sum_i F_i (1-s)^(f-i) (1+s)^i = sum_i H_i s^i
    H = [np.zeros((m, m)) for _ in range(f + 1)]
    for i in range(f + 1):
        con1 = np.array([1.0])
        for _ in range(f - i):
            con1 = np.convolve(con1, np.array([1.0, -1.0]))
        con2 = np.array([1.0])
        for _ in range(i):
            con2 = np.convolve(con2, np.array([1.0, 1.0]))
        contrib = np.convolve(con1, con2)
        for j in range(f + 1):
            H[j] = H[j] + contrib[j] * F[i]

    # Step 3, hatH_i = H_f^{-1} H_i
    Hf_inv = np.linalg.inv(H[f])
    hatH = [Hf_inv @ H[i] for i in range(f)]

    # Step 4, y and xT
    y = np.concatenate([np.ones(m), np.zeros(m * (f - 1))])
    lhs = np.hstack([hatH[0], np.ones((m, 1))])
    rhs = np.concatenate([np.zeros(m), [1.0]])
    x0T = np.linalg.lstsq(lhs.T, rhs, rcond=None)[0]
    xT = np.zeros(m * f)
    for i in range(1, f):
        xT[(i - 1) * m:i * m] = x0T @ hatH[i]
    xT[(f - 1) * m:f * m] = x0T

    # Step 5, the companion matrix E_m plus the rank-one correction
    Zold = np.zeros((m * f, m * f))
    for i in range(1, f):
        Zold[(i - 1) * m:i * m, i * m:(i + 1) * m] = np.eye(m)
    for i in range(f):
        Zold[m * (f - 1):m * f, i * m:(i + 1) * m] = -hatH[i]
    denom = float(xT @ y)
    if denom != 0:
        y = y / denom
    Zold = Zold + drift_sign * np.outer(y, xT)

    if mode == 'Schur':
        # Step 6-7, the stable invariant subspace by an ordered real Schur form
        T, U, _ = sla.schur(Zold, output='real', sort='lhp')
        Tsub = U[:, :m]
    elif mode in ('MSignStandard', 'MSignBalzer'):
        # Step 6, the classic matrix sign function iteration
        Znew = (Zold + np.linalg.inv(Zold)) / 2.0
        numit = 0
        check = 1.0
        while check > epsilon and numit < max_num_it:
            numit += 1
            Zold = Znew
            if mode == 'MSignStandard':
                determ = 0.5
            else:
                determ = 1.0 / (1.0 + abs(np.linalg.det(Zold)) ** (1.0 / (m * f)))
            Znew = determ * Zold + (1.0 - determ) * np.linalg.inv(Zold)
            check = (np.linalg.norm(Znew - Zold, 1)
                     / np.linalg.norm(Zold, 1))
        if numit == max_num_it and check > epsilon:
            warnings.warn('MG1_IS: maximum number of iterations %d reached; '
                          'T may not have m columns' % numit)
        # Step 7, an orthonormal basis of the range of Znew - I
        Tsub = sla.orth(Znew - np.eye(m * f))
    else:
        raise ValueError("MG1_IS: Mode '%s' is not one of 'Schur', "
                         "'MSignStandard', 'MSignBalzer'" % mode)

    # Step 8
    top = Tsub[0:m, :]
    bot = Tsub[m:2 * m, :]
    return (top + bot) @ np.linalg.inv(top - bot)


def gim1_r_dual(A: np.ndarray, dual: str = 'A', algor: str = 'FI') -> np.ndarray:
    """
    R of a GI/M/1-type chain through the G of its DUAL. Port of GIM1_R.m for
    Dual 'A', 'R', 'B' and Algor 'FI', 'CR'.

    There is no cyclic reduction for R directly, so the chain is transposed
    into an M/G/1-type one whose G carries the same information: the Ramaswami
    dual for a transient chain, the Bright dual (which also rescales block i by
    eta^(i-1), eta the caudal characteristic) for a positive recurrent one. 'A'
    picks between them by the drift. R comes back by the inverse similarity,
    times eta in the Bright case.

    Named `gim1_r_dual` because `gim1_r` above is a plain functional iteration
    on R itself, kept for its existing callers.
    """
    A = np.asarray(A, dtype=float)
    m = A.shape[0]
    dega = A.shape[1] // m - 1
    orig = _blocks_of(A, m)
    blocks = [b.copy() for b in orig]

    drift, theta = mg1_drift(orig)
    ram = (dual == 'R') or (dual == 'A' and drift <= 1.0)
    eta = 1.0

    if ram:
        blocks = [np.diag(1.0 / theta) @ b.T @ np.diag(theta) for b in orig]
    elif dual in ('A', 'B'):
        eta = gim1_caudal(orig) if drift > 1.0 else mg1_decay(orig)
        sumAeta = orig[dega] * (eta ** dega)
        for i in range(dega - 1, -1, -1):
            sumAeta = sumAeta + orig[i] * (eta ** i)
        theta = stat(sumAeta + (1.0 - eta) * np.eye(m))
        blocks = [(eta ** (b - 1.0)) * (np.diag(1.0 / theta) @ orig[b].T @ np.diag(theta))
                  for b in range(dega + 1)]
    else:
        raise ValueError("GIM1_R: Dual '%s' is not one of 'A', 'B', 'R'" % dual)

    wide = _hcat(blocks)
    if algor == 'FI':
        G = mg1_fi(wide)
    elif algor == 'CR':
        G = mg1_cr(wide)
    elif algor == 'NI':
        G = mg1_ni(wide)
    elif algor == 'IS':
        G = mg1_is(wide)
    elif algor == 'RR':
        raise NotImplementedError(
            "GIM1_R: Algor 'RR' (Ramaswami Reduction) is not ported; use 'FI', "
            "'CR', 'NI' or 'IS', which solve the same equation")
    else:
        raise ValueError("GIM1_R: Algor '%s' is not supported" % algor)

    R = np.diag(1.0 / theta) @ G.T @ np.diag(theta)
    if not ram:
        R = R * eta
    return R


# ---------------------------------------------------------------------------
# ETAQA: matlab/lib/thirdparty/MAMSolver
# ---------------------------------------------------------------------------


def mg1_g_etaqa(A: np.ndarray) -> np.ndarray:
    """
    G of an M/G/1-type chain, uniformized first. Port of MG1_G_ETAQA.m.

    The generator is turned into the transition matrix of the uniformized chain
    by dividing through by -min(diag(A1)) and adding the identity back onto A1,
    which is what cyclic reduction expects. The reference's `isdiscrete` flag is
    write-only (it assigns `isdicrete`), so that branch is unconditional; for
    the generators LINE passes it that is the correct branch anyway.

    References:
        Riska, A., & Smirni, E. (2003). ETAQA: An Efficient Technique for the
        Analysis of QBD-Processes by Aggregation. Performance Evaluation,
        54(2):151-177.
    """
    A = np.asarray(A, dtype=float)
    r = A.shape[0]
    if A.shape[1] % r != 0:
        raise ValueError("MG1_G_ETAQA: A is not a block sequence of A's width")
    An = A.copy()
    t = np.min(np.diag(An[:, r:2*r]))
    if t > 0:
        raise ValueError(
            "MG1_G_ETAQA: this is not a stochastic matrix, neither continuous nor "
            "discrete; every row must sum to 0 or 1")
    An = An / (-t)
    for i in range(r):
        An[i, r + i] += 1.0
    return mg1_cr(An)


def mg1_pi_etaqa(B: Optional[np.ndarray], A: np.ndarray,
                 G: Optional[np.ndarray] = None,
                 C0: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Aggregated stationary vector [pi0, pi1, pi2+pi3+...] of an M/G/1-type
    chain. Port of MG1_pi_ETAQA.m.

    ETAQA replaces the infinitely many balance equations for levels 2 and above
    by their sum, so the system is finite and (mb+2m) x (mb+2m) and EXACT: no
    level is truncated and no tail is fitted. The equations are dependent by
    one, so a redundant column is dropped (found by a rank test, since it is not
    always the last) and the normalization e^T pi = 1 takes its place.

    References:
        Stathopoulos, V., et al. (2012). ETAQA Solutions for Infinite
        Markov Processes with Repetitive Structure.
    """
    A = np.asarray(A, dtype=float).copy()
    m = A.shape[0]
    dega = A.shape[1] // m - 1
    A0 = A[:, 0:m].copy()

    if B is None or np.asarray(B).size == 0:
        mb = m
        degb = dega
        Bw = A.copy()
    else:
        Bw = np.asarray(B, dtype=float).copy()
        mb = Bw.shape[0]
        if (Bw.shape[1] - mb) % m != 0:
            raise ValueError("MG1_pi_ETAQA: matrix B has an incorrect number of columns")
        degb = (Bw.shape[1] - mb) // m

    C0m = A0 if C0 is None else np.asarray(C0, dtype=float)
    if C0 is None and mb != m:
        raise ValueError(
            "MG1_pi_ETAQA: the Boundary option must be used since a dimension of B0 is "
            "not identical to A0")
    if C0 is not None and (C0m.shape[0] != m or C0m.shape[1] != mb):
        raise ValueError("MG1_pi_ETAQA: the boundary parameter value has an incorrect dimension")

    # A transition matrix is turned into a generator, as the reference tests it.
    tot = float(np.sum(np.sum(Bw, axis=1)))
    if tot > 1e-12 and (tot - mb) < 1e-12:
        Bw[:mb, :mb] -= np.eye(mb)
        A[:, m:2*m] -= np.eye(m)

    Ab = _blocks_of(A, m)
    if G is None:
        G = mg1_g_etaqa(A)
    G = np.asarray(G, dtype=float)

    drift, _ = mg1_drift(Ab)
    if drift >= 1.0:
        raise ValueError(
            "MG1_pi_ETAQA: the Markov chain characterized by A is not positive recurrent "
            "(drift = %g)" % drift)

    # Shat(j) = B(j) + B(j+1) G + B(j+2) G^2 + ..., j = 1..degb
    Shat = [Bw[:, mb + (degb-1)*m:mb + degb*m].copy()]
    for i in range(degb - 1, 0, -1):
        Shat.insert(0, Bw[:, mb + (i-1)*m:mb + i*m] + Shat[0] @ G)

    # S(j) = A(j) + A(j+1) G + A(j+2) G^2 + ..., j = 1..dega
    if dega <= 1:
        raise ValueError(
            "MG1_pi_ETAQA: the number of repetitive state blocks is less than 2, this is "
            "not an irreducible Markov chain")
    S = [Ab[dega].copy()]
    for i in range(dega - 1, 0, -1):
        S.insert(0, Ab[i] + S[0] @ G)

    firstc = np.ones((mb + 2*m, 1))
    secondc = np.vstack([Bw[:, 0:mb], C0m, np.zeros((m, mb))])

    if len(Shat) < 2:
        Shat.append(np.zeros((mb, m)))
    if len(S) < 2:
        S.append(np.zeros((m, m)))

    thirdc = np.vstack([Bw[:, mb:mb+m] + Shat[1] @ G,
                        Ab[1] + S[1] @ G,
                        np.zeros((m, m))])

    Bsum = np.zeros((mb, m))
    Shat_sum = np.zeros((mb, m))
    if degb <= 2:
        if degb == 2:
            Bsum = Bw[:, mb+m:mb+2*m].copy()
    else:
        for i in range(2, degb):
            Bsum = Bsum + Bw[:, mb + (i-1)*m:mb + i*m]
            Shat_sum = Shat_sum + Shat[i]
        Bsum = Bsum + Bw[:, mb + (degb-1)*m:]

    Asum = np.zeros((m, m))
    Ssum = np.zeros((m, m))
    if dega >= 3:
        for i in range(2, dega):
            Ssum = Ssum + S[i]
            Asum = Asum + Ab[i]
        Asum = Asum + Ab[dega]
    elif dega == 2:
        Asum = Asum + Ab[dega]
    else:
        raise ValueError(
            "MG1_pi_ETAQA: the number of repetitive state blocks is less than 3, the "
            "Markov chain is reducible")

    fourthc = np.vstack([Bsum + Shat_sum @ G,
                         Asum + Ssum @ G,
                         Asum + Ab[1] + (Ssum + S[1]) @ G])

    Xtemp = np.hstack([secondc, thirdc, fourthc])

    full = np.linalg.matrix_rank(Xtemp)
    n = mb + 2*m
    drop = n - 1
    for i in range(n):
        if np.linalg.matrix_rank(np.delete(Xtemp, i, axis=1)) == full:
            drop = i
            break
    Xtemp = np.delete(Xtemp, drop, axis=1)

    Xnew = np.hstack([firstc, Xtemp])
    rside = np.zeros(n)
    rside[0] = 1.0
    return la.solve(Xnew.T, rside)


def mg1_qlen_etaqa(B: Optional[np.ndarray], A: np.ndarray,
                   pi: np.ndarray, n: int,
                   C0: Optional[np.ndarray] = None) -> float:
    """
    n-th moment of the level (the queue length) of an M/G/1-type chain from the
    ETAQA aggregates. Port of MG1_qlen_ETAQA.m.

    The moment is NOT read off the aggregates directly: r^(k) = sum_{j>=2} j^k
    pi_j is propagated through a recurrence whose left-hand side is one fixed
    m x m system, so the n-th moment costs n solves of that size however heavy
    the tail is.
    """
    A = np.asarray(A, dtype=float)
    m = A.shape[0]
    dega = A.shape[1] // m - 1
    Ab = _blocks_of(A, m)

    if B is None or np.asarray(B).size == 0:
        mb = m
        degb = dega
        Bw = A.copy()
    else:
        Bw = np.asarray(B, dtype=float)
        mb = Bw.shape[0]
        if (Bw.shape[1] - mb) % m != 0:
            raise ValueError("MG1_qlen_ETAQA: matrix B has an incorrect number of columns")
        degb = (Bw.shape[1] - mb) // m

    if C0 is None and mb != m:
        raise ValueError(
            "MG1_qlen_ETAQA: the Boundary option must be used since the column size of B0 "
            "is not identical to A0")

    pi = np.asarray(pi, dtype=float).flatten()
    if abs(np.sum(pi) - 1.0) > 1e-10:
        raise ValueError("MG1_qlen_ETAQA: the input probability vector does not sum up to 1")
    if (len(pi) - mb) % m != 0:
        raise ValueError("MG1_qlen_ETAQA: the probability vector has an incorrect number of columns")

    pi0 = pi[:mb]
    pi1 = pi[mb:mb+m]
    pistar = pi[mb+m:mb+2*m]

    Asum = Ab[0].copy()
    for i in range(1, dega + 1):
        Asum = Asum + Ab[i]
    F11 = Ab[2].copy()
    for i in range(3, dega + 1):
        F11 = F11 + (i - 1) * Ab[i]
    lsleft = np.zeros((m, m))
    lsleft[:, :m-1] = Asum[:, :m-1]
    lsleft[:, m-1] = np.sum(F11 - Ab[0], axis=1)

    # Fhat0(j) = sum_{l>=j} B(l), j = 1..degb
    Fhat0j = [Bw[:, mb + (degb-1)*m:].copy()]
    for j in range(degb - 1, 0, -1):
        Fhat0j.insert(0, Bw[:, mb + (j-1)*m:mb + j*m] + Fhat0j[0])

    # F0(j) = sum_{l>=j+1} A(l): built from block dega down to 2 but indexed
    # from 1, which is the reference's own off-by-one and is reproduced.
    F0j = [Ab[dega].copy()]
    for j in range(dega - 1, 1, -1):
        F0j.insert(0, Ab[j] + F0j[0])

    r = [pistar]

    frestsaver = []
    fcrestsaver = []
    for l in range(1, n + 1):
        t1 = np.zeros((m, m))
        for j in range(2, dega + 1):
            t1 = t1 + (j ** l) * Ab[j]
        frestsaver.append(t1)
        t2 = np.zeros((m, m))
        for j in range(1, dega + 1):
            if j <= len(F0j):
                t2 = t2 + (j ** l) * F0j[j-1]
        fcrestsaver.append(t2)

    for k in range(1, n + 1):
        fhatkM = np.zeros((mb, m))
        for j in range(1, degb + 1):
            fhatkM = fhatkM + ((j + 1) ** k) * Bw[:, mb + (j-1)*m:mb + j*m]
        fhatk = pi0 @ fhatkM

        fkM = np.zeros((m, m))
        for j in range(2, dega + 1):
            fkM = fkM + ((j + 1) ** k) * Ab[j]
        fkM = (2 ** k) * Ab[1] + fkM
        fk = pi1 @ fkM

        frest = np.zeros(m)
        for l in range(1, k + 1):
            frest = frest + _bino(k, l) * (r[k-l] @ (Ab[1] + frestsaver[l-1]))

        bk = -fhatk - fk - frest

        fchatkM = np.zeros((mb, m))
        for j in range(2, degb + 1):
            fchatkM = fchatkM + (j ** k) * Fhat0j[j-1]
        fchatk = float(np.sum(pi0 @ fchatkM))

        fckM = np.zeros((m, m))
        for j in range(1, dega):
            if j <= len(F0j):
                fckM = fckM + ((j + 1) ** k) * F0j[j-1]
        fck = float(np.sum(pi1 @ fckM))

        fcrest = 0.0
        for l in range(1, k + 1):
            fcrest += _bino(k, l) * float(np.sum(r[k-l] @ fcrestsaver[l-1]))

        ck = -fchatk - fck - fcrest

        rside = np.zeros(m)
        rside[:m-1] = bk[:m-1]
        rside[m-1] = ck
        r.append(la.solve(lsleft.T, rside))

    return float(np.sum(r[-1]) + np.sum(pi1))


def _bino(n: int, k: int) -> float:
    """The binomial coefficient the reference computes as a ratio of factorials."""
    num = 1.0
    for i in range(1, n + 1):
        num *= i
    den = 1.0
    for i in range(1, k + 1):
        den *= i
    for i in range(1, n - k + 1):
        den *= i
    return num / den


def gim1_r_etaqa(A: np.ndarray) -> np.ndarray:
    """
    R of a GI/M/1-type chain, uniformized first. Port of GIM1_R_ETAQA.m.

    `A` is the VERTICAL stack [A0; A1; ...; Amax]; the reference transposes it
    into the horizontal sequence GIM1_R wants and asks for the automatic dual
    with functional iterations.
    """
    A = np.asarray(A, dtype=float)
    s = A.shape[1]
    if A.shape[0] % s != 0:
        raise ValueError("GIM1_R_ETAQA: A is not a vertical stack of square blocks")
    An = A.copy()
    t = np.min(np.diag(An[s:2*s, :]))
    if not t > 0:
        An = An / (-t)
        An[s:2*s, :] += np.eye(s)
    blocks = _vblocks_of(An, s)
    return gim1_r_dual(_hcat(blocks), 'A', 'FI')


def gim1_pi_etaqa(B: np.ndarray, A: np.ndarray, R: np.ndarray,
                  B0: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Aggregated stationary vector [pi0, pi1, pi2+pi3+...] of a GI/M/1-type
    chain. Port of GIM1_pi_ETAQA.m. `B` and `A` are vertical stacks; `B0` is the
    reference's 'Boundary' option (None for the default A0).
    """
    B = np.asarray(B, dtype=float).copy()
    A = np.asarray(A, dtype=float).copy()
    R = np.asarray(R, dtype=float)
    m = R.shape[0]
    mb = B.shape[1]
    if (B.shape[0] - mb) % m != 0:
        raise ValueError("GIM1_pi_ETAQA: input matrix B has an incorrect number of rows")
    degb = (B.shape[0] - mb) // m
    if A.shape[0] % m != 0:
        raise ValueError("GIM1_pi_ETAQA: input matrix A has an incorrect number of rows")
    dega = A.shape[0] // m - 1

    temp0 = la.inv(np.eye(m) - R)
    if np.all(np.any(temp0 < -100 * np.finfo(float).eps, axis=0)):
        raise ValueError(
            "GIM1_pi_ETAQA: the spectral radius of R is not below 1, GIM1 is not positive "
            "recurrent")

    B0m = A[0:m, :].copy() if B0 is None else np.asarray(B0, dtype=float)
    if B0 is not None and (B0m.shape[0] != mb or B0m.shape[1] != m):
        raise ValueError("GIM1_pi_ETAQA: Boundary has an incorrect dimension")

    Btop = B[:mb, :].copy()
    test = Btop + B0m
    if abs(float(np.sum(np.sum(test, axis=1) - 1.0))) < 1e-10:
        B[:mb, :mb] -= np.eye(mb)
        A[m:2*m, :] -= np.eye(m)
        Btop = B[:mb, :].copy()

    Ab = _vblocks_of(A, m)

    def bblk(i):
        return B[mb + (i-1)*m:mb + i*m, :]

    firstc = np.ones((mb + 2*m, 1))

    temp = np.eye(m) - R
    tempsum = np.zeros((m, mb))
    for i in range(2, degb + 1):
        tempsum = tempsum + temp @ bblk(i)
        temp = R @ temp
    secondc = np.vstack([Btop, bblk(1), tempsum])

    temp = np.eye(m) - R
    tempsum = np.zeros((m, m))
    for i in range(2, dega + 1):
        tempsum = tempsum + temp @ Ab[i]
        temp = R @ temp
    thirdc = np.vstack([B0m, Ab[1], tempsum])

    temp = R.copy()
    tempsum = np.zeros((m, m))
    for i in range(2, dega + 1):
        tempsum = tempsum + temp @ Ab[i]
        temp = R @ temp
    fourthc = np.vstack([np.zeros((mb, m)), Ab[0], Ab[0] + Ab[1] + tempsum])
    fourthc = fourthc[:, :-1]

    X = np.hstack([firstc, secondc, thirdc, fourthc])
    rside = np.zeros(mb + 2*m)
    rside[0] = 1.0
    return la.solve(X.T, rside)


def gim1_qlen_etaqa(B: np.ndarray, A: np.ndarray, R: np.ndarray,
                    pi: np.ndarray, n: int,
                    B0: Optional[np.ndarray] = None) -> float:
    """
    n-th moment of the level of a GI/M/1-type chain from the ETAQA aggregates.
    Port of GIM1_qlen_ETAQA.m.

    REPRODUCED REFERENCE DEFECT: the accumulator starts at `A(3)`, a SCALAR at
    column-major linear index 3 of the stacked A where the third BLOCK is meant,
    which MATLAB then broadcasts. It corrupts the last column of the moment
    system whenever m > 1 and can make the reported mean NEGATIVE. The port
    reproduces value and broadcast, including the case where the loop that would
    turn the scalar into a matrix never runs.
    """
    B = np.asarray(B, dtype=float).copy()
    A = np.asarray(A, dtype=float).copy()
    R = np.asarray(R, dtype=float)
    m = R.shape[0]
    mb = B.shape[1]
    if (B.shape[0] - mb) % m != 0:
        raise ValueError("GIM1_qlen_ETAQA: input matrix B has an incorrect number of rows")
    degb = (B.shape[0] - mb) // m
    if A.shape[0] % m != 0:
        raise ValueError("GIM1_qlen_ETAQA: input matrix A has an incorrect number of rows")
    dega = A.shape[0] // m - 1

    temp0 = la.inv(np.eye(m) - R)
    if np.all(np.any(temp0 < -100 * np.finfo(float).eps, axis=0)):
        raise ValueError(
            "GIM1_qlen_ETAQA: the spectral radius of R is not below 1, GIM1 is not positive "
            "recurrent")

    B0m = A[0:m, :].copy() if B0 is None else np.asarray(B0, dtype=float)
    if B0 is not None and (B0m.shape[0] != mb or B0m.shape[1] != m):
        raise ValueError("GIM1_qlen_ETAQA: Boundary has an incorrect dimension")

    pi = np.asarray(pi, dtype=float).flatten()
    pi0 = pi[:mb]
    pi1 = pi[mb:mb+m]
    pistar = pi[mb+m:mb+2*m]

    if n == 0:
        return 1.0

    # The scalar of the reproduced defect, at column-major linear index 3.
    a3 = float(A[2 % A.shape[0], 2 // A.shape[0]])

    Btop = B[:mb, :].copy()
    test = Btop + B0m
    if abs(float(np.sum(np.sum(test, axis=1) - 1.0))) < 1e-10:
        B[:mb, :mb] -= np.eye(mb)
        A[m:2*m, :] -= np.eye(m)

    Ab = _vblocks_of(A, m)

    def bblk(i):
        return B[mb + (i-1)*m:mb + i*m, :]

    Rpower = np.eye(m)
    lsum = Ab[0] + Ab[1]
    for i in range(2, dega + 1):
        lsum = lsum + Rpower @ Ab[i]
        Rpower = R @ Rpower

    leftr = np.zeros(m)
    if degb >= 2 and dega >= 2:
        loop_runs = dega >= 3
        part1 = np.full((m, m), a3)
        part2 = np.zeros((m, m))
        Rpower = R.copy()
        for i in range(1, dega - 1):
            part1 = part1 + Rpower @ Ab[i + 2]
            part2 = part2 + i * (Rpower @ Ab[i + 2])
            Rpower = R @ Rpower
        if loop_runs:
            leftr = np.sum(part1, axis=1) + np.sum(part2, axis=1) - np.sum(Ab[0], axis=1)
        else:
            leftr = a3 * np.ones(m) - np.sum(Ab[0], axis=1)
    elif degb == 1 and dega != 1:
        acc = np.zeros((m, m))
        Rpower = np.eye(m)
        for i in range(2, dega + 1):
            acc = acc + (i - 1) * (Rpower @ Ab[i])
            Rpower = R @ Rpower
        leftr = np.sum(acc, axis=1) - np.sum(Ab[0], axis=1)
    else:
        raise ValueError(
            "GIM1_qlen_ETAQA: the number of A blocks is not enough, this is a reducible "
            "Markov chain")

    lsleft = np.zeros((m, m))
    lsleft[:, :m-1] = lsum[:, :m-1]
    lsleft[:, m-1] = leftr

    r = [pistar]
    reuse = []

    for k in range(1, n + 1):
        bk = -((2 ** k) * (pi0 @ B0m)
               + pi1 @ ((2 ** k) * Ab[1] + (3 ** k) * Ab[0]))
        for l in range(1, k + 1):
            bk = bk - _bino(k, l) * (r[k-l] @ (Ab[1] + (2 ** l) * Ab[0]))

        tempsum = np.zeros((m, m))
        Rpower = R.copy()
        for i in range(1, dega - 1):
            t = sum(z ** k for z in range(1, i + 1))
            tempsum = tempsum + t * (Rpower @ Ab[i + 2])
            Rpower = Rpower @ R
        reuse.append(np.sum(Ab[0] - tempsum, axis=1))

        ck = (2 ** k) * float(np.sum(pi1 @ Ab[0]))
        tempsum2 = np.zeros((m, mb))
        Rpower = R.copy()
        for i in range(2, degb + 1):
            t = sum(z ** k for z in range(2, i + 1))
            tempsum2 = tempsum2 + t * (Rpower @ bblk(i))
            Rpower = R @ Rpower
        ck -= float(np.sum(pi1 @ tempsum2))
        for l in range(1, k + 1):
            ck += _bino(k, l) * float(r[k-l] @ reuse[l-1])

        rside = np.zeros(m)
        rside[:m-1] = bk[:m-1]
        rside[m-1] = ck
        r.append(la.solve(lsleft.T, rside))

    return float(np.sum(r[-1]) + np.sum(pi1))


__all__ = [
    'stat',
    'qbd_cr',
    'qbd_caudal',
    'qbd_eg',
    'qbd_pi',
    'qbd_lr',
    'mg1_g',
    'mg1_eg',
    'gim1_r',
    'gim1_pi',
    # ETAQA functions
    'mg1_drift',
    'mg1_decay',
    'gim1_caudal',
    'mg1_shifts',
    'mg1_cr',
    'mg1_fi',
    'gim1_r_dual',
    'gim1_r_etaqa',
    'gim1_pi_etaqa',
    'gim1_qlen_etaqa',
    'mg1_g_etaqa',
    'mg1_pi_etaqa',
    'mg1_qlen_etaqa',
]
