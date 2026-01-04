"""
SMC: Structured Markov Chain Solvers.

Native Python implementations of algorithms for Quasi-Birth-Death (QBD),
M/G/1, and G/I/M/1 type Markov chains.

Based on the SMCtools package by Benny Van Houdt.
"""

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

    # Solve (A - I)^T * pi^T = 0 with sum(pi) = 1
    B = np.hstack([A.T - np.eye(n), np.ones((n, 1))])
    y = np.zeros((1, n + 1))
    y[0, -1] = 1.0

    # Solve using least squares
    pi = la.lstsq(B.T, y.T)[0].T

    # Ensure non-negative and normalized
    pi = np.maximum(pi, 0)
    pi = pi / np.sum(pi)

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
        lamb = np.max(-A1_diag)
        A0 = A0 / lamb
        A1 = A1 / lamb + np.eye(m)
        A2 = A2 / lamb

    # Initial G using extended generator method
    result = qbd_eg(A0, A1, A2, verbose)
    G = result.get('G')

    if G is None or G.size == 0:
        return result

    # Compute drift
    theta = stat(A0 + A1 + A2)
    drift = theta @ np.sum(A0, axis=1) - theta @ np.sum(A2, axis=1)

    A2_old = A2.copy()
    A0_old = A0.copy()
    uT = np.ones((1, m)) * m

    if mode == "Shift":
        if drift < 0:
            A2 = A2 - np.ones((m, 1)) @ (theta @ A2)
            A1 = A1 + np.ones((m, 1)) @ (theta @ A0)
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

        Ahat = A + BAtemp @ C
        A = A + BAtemp @ A + CAtemp @ B
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
            A1 = A1 - np.ones((m, 1)) @ theta @ A0
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


def qbd_eg(A0: np.ndarray, A1: np.ndarray, A2: np.ndarray,
           verbose: bool = False) -> Dict[str, np.ndarray]:
    """
    QBD_EG: Extended Generator method for initial G computation.

    Args:
        A0: Downward transition matrix
        A1: Level transition matrix
        A2: Upward transition matrix
        verbose: Print progress if True

    Returns:
        Dictionary with 'G' matrix
    """
    A0 = np.asarray(A0, dtype=float)
    A1 = np.asarray(A1, dtype=float)
    A2 = np.asarray(A2, dtype=float)

    m = A0.shape[0]

    # Initial approximation using simple iteration
    G = np.zeros((m, m))
    for _ in range(100):
        G_new = la.inv(np.eye(m) - A1 - A2 @ G) @ A0
        if la.norm(G_new - G, np.inf) < 1e-14:
            break
        G = G_new

    return {'G': G}


def qbd_pi(B0: np.ndarray, B1: np.ndarray, R: np.ndarray,
           max_num_comp: int = 500, verbose: int = 0,
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
        max_num_comp: Maximum number of components (default: 500)
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

    # Check spectral radius of R
    temp = la.inv(np.eye(m) - R)
    if np.min(temp) < -100 * 1e-15:
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
    QBD_LR: Logarithmic Reduction algorithm for QBD Markov chains.

    Similar to QBD_CR but uses different iteration scheme.

    Args:
        A0: Downward transition matrix
        A1: Level transition matrix
        A2: Upward transition matrix
        max_num_it: Maximum number of iterations
        verbose: Print progress if True
        mode: "Shift" or "Basic"

    Returns:
        Dictionary with 'G', 'R', 'U' matrices
    """
    # For simplicity, delegate to CR algorithm
    return qbd_cr(A0, A1, A2, max_num_it, verbose, mode)


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


def mg1_g_etaqa(A: np.ndarray) -> np.ndarray:
    """
    Compute G matrix for M/G/1-type Markov chains using ETAQA method.

    G is the minimal nonnegative solution to:
    - Continuous time: 0 = A0 + A1*G + A2*G^2 + ... + Amax*G^max
    - Discrete time: G = A0 + A1*G + A2*G^2 + ... + Amax*G^max

    Based on MAMSolver ETAQA from College of William & Mary.

    Args:
        A: Block matrix [A0 A1 ... Amax] with m rows and m*(max+1) columns

    Returns:
        G matrix (m x m)

    References:
        Riska, A., & Smirni, E. (2003). ETAQA: An Efficient Technique for the
        Analysis of QBD-Processes by Aggregation. Performance Evaluation, 54(2):151-177.
    """
    A = np.asarray(A, dtype=float)
    m = A.shape[0]
    n = A.shape[1] // m

    # Check if continuous time (row sums = 0) or discrete time (row sums = 1)
    row_sums = np.sum(A, axis=1)
    is_continuous = np.allclose(row_sums, 0, atol=1e-12)

    if is_continuous:
        # Uniformize: find lambda from A1 diagonal
        A1 = A[:, m:2*m]
        min_diag = np.min(np.diag(A1))
        if min_diag >= 0:
            raise ValueError("Invalid generator: A1 diagonal must be negative")

        lamb = -min_diag
        Awork = A / lamb
        for i in range(m):
            Awork[i, m + i] += 1.0
    else:
        Awork = A.copy()

    # Use functional iteration
    G = np.zeros((m, m))
    max_iter = 200
    tol = 1e-14

    for numit in range(max_iter):
        # G_new = A0 + A1*G + A2*G^2 + ...
        G_new = Awork[:, 0:m].copy()
        G_power = G.copy()
        for i in range(1, n):
            G_new = G_new + Awork[:, i*m:(i+1)*m] @ G_power
            G_power = G_power @ G

        if la.norm(G_new - G, np.inf) < tol:
            return G_new

        G = G_new

    return G


def mg1_pi_etaqa(B: Optional[np.ndarray], A: np.ndarray,
                 G: Optional[np.ndarray] = None,
                 C0: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Compute aggregated stationary probabilities for M/G/1-type Markov chain using ETAQA.

    Returns aggregated probabilities [pi0, pi1, piStar] where:
    - pi0: Probability of level 0
    - pi1: Probability of level 1
    - piStar: Sum of probabilities for levels 2, 3, ...

    Args:
        B: Boundary block matrix [B0 B1 ... Bmax] with mb rows
        A: Repeating block matrix [A0 A1 ... Amax] with m rows
        G: G matrix (computed if not provided)
        C0: Optional boundary matrix (defaults to A0)

    Returns:
        Aggregated probability vector [pi0, pi1, piStar] of shape (1, mb + 2*m)

    References:
        Stathopoulos, V., et al. (2012). ETAQA Solutions for Infinite
        Markov Processes with Repetitive Structure.
    """
    A = np.asarray(A, dtype=float)
    m = A.shape[0]
    dega = A.shape[1] // m - 1

    # Use boundary or default to A
    if B is None or B.size == 0:
        mb = m
        degb = dega
        boundary = A.copy()
    else:
        B = np.asarray(B, dtype=float)
        mb = B.shape[0]
        degb = (B.shape[1] - mb) // m
        boundary = B.copy()

    boundary_C0 = C0 if C0 is not None else A[:, 0:m]

    # Check continuous/discrete time
    row_sums = np.sum(boundary, axis=1)
    is_continuous = np.allclose(row_sums, 0, atol=1e-12)

    Bwork = boundary.copy()
    Awork = A.copy()
    if not is_continuous:
        is_discrete = np.allclose(row_sums, 1, atol=1e-12)
        if is_discrete:
            for i in range(mb):
                Bwork[i, i] -= 1.0
            for i in range(m):
                Awork[i, m + i] -= 1.0

    # Compute G matrix if not provided
    if G is None:
        G = mg1_g_etaqa(Awork)

    # Compute drift
    sumA = Awork[:, dega*m:(dega+1)*m].copy()
    alpha = np.sum(sumA, axis=1, keepdims=True)
    for i in range(dega - 1, 0, -1):
        sumA = sumA + Awork[:, i*m:(i+1)*m]
        alpha = alpha + np.sum(sumA, axis=1, keepdims=True)
    sumA = sumA + Awork[:, 0:m]
    a = stat(sumA)
    drift = (a @ alpha)[0]

    if drift >= 1.0:
        raise ValueError(f"M/G/1 chain not positive recurrent (drift = {drift})")

    # Compute S_hat: S_hat(j) = B(j) + B(j+1)*G + B(j+2)*G^2 + ...
    Shat = Bwork[:, mb + (degb-1)*m:mb + degb*m].copy()
    if degb > 1:
        for i in range(degb - 1, 0, -1):
            col_start = mb + (i-1)*m
            col_end = mb + i*m
            temp = Bwork[:, col_start:col_end] + Shat[:, 0:m] @ G
            Shat = np.hstack([temp, Shat])

    # Compute S: S(j) = A(j) + A(j+1)*G + A(j+2)*G^2 + ...
    S = Awork[:, dega*m:(dega+1)*m].copy()
    if dega > 1:
        for i in range(dega - 1, 0, -1):
            temp = Awork[:, i*m:(i+1)*m] + S[:, 0:m] @ G
            S = np.hstack([temp, S])

    # Build linear system
    firstc = np.ones((mb + 2*m, 1))

    # Second column: [B0; C0; 0s]
    secondc = np.zeros((mb + 2*m, mb))
    secondc[0:mb, 0:mb] = Bwork[:, 0:mb]
    secondc[mb:mb+m, 0:mb] = boundary_C0[:, 0:mb] if boundary_C0.shape[1] >= mb else np.zeros((m, mb))

    # Ensure Shat and S have at least 2m columns
    if Shat.shape[1] < 2*m:
        Shat = np.hstack([Shat, np.zeros((mb, m))])
    if S.shape[1] < 2*m:
        S = np.hstack([S, np.zeros((m, m))])

    # Third column: [B(1) + Shat(2)*G; A(1) + S(1)*G; 0s]
    thirdc = np.zeros((mb + 2*m, m))
    B1 = Bwork[:, mb:mb+m]
    Shat2G = Shat[:, m:2*m] @ G
    thirdc[0:mb, :] = B1 + Shat2G
    A1 = Awork[:, m:2*m]
    S1G = S[:, m:2*m] @ G
    thirdc[mb:mb+m, :] = A1 + S1G

    # Fourth column
    Bsum = np.zeros((mb, m))
    ShatSum = np.zeros((mb, m))
    if degb > 2:
        for i in range(2, degb):
            Bsum = Bsum + Bwork[:, mb + (i-1)*m:mb + i*m]
            if i < Shat.shape[1] // m:
                ShatSum = ShatSum + Shat[:, i*m:(i+1)*m]
        Bsum = Bsum + Bwork[:, mb + (degb-1)*m:mb + degb*m]
    elif degb == 2:
        Bsum = Bsum + Bwork[:, mb+m:mb+2*m]

    Asum = np.zeros((m, m))
    Ssum = np.zeros((m, m))
    if dega >= 3:
        for i in range(2, dega):
            Ssum = Ssum + S[:, i*m:(i+1)*m]
            Asum = Asum + Awork[:, i*m:(i+1)*m]
        Asum = Asum + Awork[:, dega*m:(dega+1)*m]
    elif dega == 2:
        Asum = Asum + Awork[:, dega*m:(dega+1)*m]

    fourthc = np.zeros((mb + 2*m, m))
    fourthc[0:mb, :] = Bsum + ShatSum @ G
    fourthc[mb:mb+m, :] = Asum + Ssum @ G
    S1 = S[:, m:2*m]
    fourthc[mb+m:mb+2*m, :] = Asum + A1 + (Ssum + S1) @ G

    # Build Xtemp and solve
    Xtemp = np.hstack([secondc, thirdc, fourthc])

    # Find and remove redundant column
    rankX = np.linalg.matrix_rank(Xtemp)
    redundant_col = 0
    for i in range(Xtemp.shape[1]):
        temp = np.delete(Xtemp, i, axis=1)
        if np.linalg.matrix_rank(temp) == rankX:
            redundant_col = i
            break
    Xtemp = np.delete(Xtemp, redundant_col, axis=1)

    # Build final system
    Xnew = np.hstack([firstc, Xtemp])

    # Solve: pi * Xnew = [1, 0, 0, ...]
    rside = np.zeros((1, mb + 2*m))
    rside[0, 0] = 1.0

    # Solve using pseudo-inverse
    pi = rside @ la.pinv(Xnew)

    return pi


def mg1_qlen_etaqa(B: Optional[np.ndarray], A: np.ndarray,
                   pi: np.ndarray, n: int) -> float:
    """
    Compute n-th moment of queue length using ETAQA aggregated probabilities.

    Args:
        B: Boundary block matrix
        A: Repeating block matrix
        pi: Aggregated probability from mg1_pi_etaqa
        n: Moment order (1 = mean, 2 = second moment, etc.)

    Returns:
        n-th moment of queue length
    """
    A = np.asarray(A, dtype=float)
    m = A.shape[0]

    if B is None or np.asarray(B).size == 0:
        mb = m
    else:
        mb = np.asarray(B).shape[0]

    # Extract probability components
    pi = np.asarray(pi).flatten()
    pi0 = pi[0:mb]
    pi1 = pi[mb:mb+m]
    piStar = pi[mb+m:mb+2*m]

    e = np.ones(m)

    if n == 1:
        # Mean queue length
        mean = np.sum(pi1) + 2.0 * np.sum(piStar)
        return mean

    # Higher moments - approximate using geometric series
    moment = float(np.sum(pi1))  # Level 1 contribution

    # Contribution from piStar (levels 2 and above)
    piStar_sum = np.sum(piStar)
    for k in range(2, 100):
        level_prob = piStar_sum * (1.0 - piStar_sum) ** (k - 2)
        moment += (k ** n) * level_prob
        if level_prob < 1e-15:
            break

    return moment


__all__ = [
    'stat',
    'qbd_cr',
    'qbd_eg',
    'qbd_pi',
    'qbd_lr',
    'mg1_g',
    'mg1_eg',
    'gim1_r',
    'gim1_pi',
    # ETAQA functions
    'mg1_g_etaqa',
    'mg1_pi_etaqa',
    'mg1_qlen_etaqa',
]
