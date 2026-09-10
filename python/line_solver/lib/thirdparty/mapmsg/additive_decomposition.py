"""
Additive Decomposition for Multi-Regime Markov Fluid Queues.

Port of AdditiveDecomposition.m from the MAPMsG MATLAB library.

Performs ordered Schur decomposition and additive decomposition of the
generator matrix for each regime of a multi-regime Markov fluid queue (MRMFQ).

Reference:
    O. Gursoy, K. A. Mehr, N. Akar, "The MAP/M/s + G Call Center Model with
    Generally Distributed Patience Times."
"""

import numpy as np
from scipy.linalg import schur, solve_sylvester, expm
from numpy.linalg import inv


def _ordered_schur_3way(A, tol=1e-7):
    """Compute a Schur decomposition of A ordered as [zero | negative | positive].

    Uses scipy's ``schur`` with the ``sort`` parameter in two passes to
    achieve the desired 3-way ordering matching MATLAB's ordschur behaviour.

    Parameters
    ----------
    A : ndarray, shape (n, n)
        Matrix to decompose.
    tol : float
        Threshold for classifying eigenvalues as zero.

    Returns
    -------
    Z1 : ndarray, shape (n, n)
        Orthogonal Schur vectors (reordered).
    D1 : ndarray, shape (n, n)
        Upper quasi-triangular Schur form (reordered).
    zero_count : int
        Number of zero eigenvalues.
    negative_count : int
        Number of negative eigenvalues.
    positive_count : int
        Number of positive eigenvalues.
    """
    n = A.shape[0]

    # Step 1: Schur decomposition sorting zero eigenvalues first
    # "Zero" eigenvalues satisfy |Re(lambda)| <= tol
    D1, Z1, sdim_zero = schur(A, output='real',
                               sort=lambda x: abs(x.real) <= tol)
    zero_count = sdim_zero

    # Step 2: Reorder the remaining (non-zero) block so that negative
    # eigenvalues come before positive ones.
    # Extract the lower-right block corresponding to non-zero eigenvalues.
    if zero_count < n:
        A_rest = D1[zero_count:, zero_count:].copy()
        # Schur decompose the remaining block with negative first
        D_rest, Z_rest, sdim_neg = schur(A_rest, output='real',
                                          sort=lambda x: x.real < -tol)
        negative_count = sdim_neg
        positive_count = n - zero_count - negative_count

        # Update D1: replace lower-right block
        D1[zero_count:, zero_count:] = D_rest
        # Update coupling block: D1[0:zero_count, zero_count:]
        D1[:zero_count, zero_count:] = D1[:zero_count, zero_count:] @ Z_rest
        # Update Z1: Z1[:, zero_count:] = Z1[:, zero_count:] @ Z_rest
        Z1[:, zero_count:] = Z1[:, zero_count:] @ Z_rest
    else:
        negative_count = 0
        positive_count = 0

    return Z1, D1, zero_count, negative_count, positive_count


def additive_decomposition(Qmulti, driftregimes_multi, Bmulti):
    """Perform additive decomposition for each regime.

    Parameters
    ----------
    Qmulti : ndarray, shape (m, m, num_regimes)
        3-D array of infinitesimal generator matrices, one per regime.
    driftregimes_multi : ndarray, shape (num_regimes, m)
        Drift values for each regime and state.
    Bmulti : ndarray, shape (num_regimes,)
        Upper boundary levels for each regime.

    Returns
    -------
    lower_bound_solutions : list of ndarray
        f(lower_bound) matrices for each regime (coefficients = 1).
    upper_bound_solutions : list of ndarray
        f(upper_bound) matrices for each regime (coefficients = 1).
    integral_solutions : list of ndarray
        Integral F(B) matrices for each regime (coefficients = 1).
    Lzero_multi : list of ndarray
        Zero-eigenvalue projection matrices for each regime.
    Lneg_multi : list of ndarray
        Negative-eigenvalue projection matrices for each regime.
    Lpos_multi : list of ndarray
        Positive-eigenvalue projection matrices for each regime.
    Aneg_multi : list of ndarray
        Negative eigenvalue block matrices for each regime.
    Apos_multi : list of ndarray
        Positive eigenvalue block matrices for each regime.
    """
    num_regimes = Qmulti.shape[2]

    lower_bound_solutions = [None] * num_regimes
    upper_bound_solutions = [None] * num_regimes
    integral_solutions = [None] * num_regimes
    Lzero_multi = [None] * num_regimes
    Lneg_multi = [None] * num_regimes
    Lpos_multi = [None] * num_regimes
    Aneg_multi = [None] * num_regimes
    Apos_multi = [None] * num_regimes

    for multi in range(num_regimes):
        drifts = driftregimes_multi[multi, :].copy()
        Q = Qmulti[:, :, multi].copy()

        # Find zero-drift states
        zerodrift = np.where(drifts == 0)[0]

        if len(zerodrift) > 0:
            # Swap zero-drift states to the end of the matrix
            for i in range(len(zerodrift)):
                last = len(Q) - 1 - i
                nextzero = zerodrift[i]

                # Swap rows
                temp_row = Q[last, :].copy()
                Q[last, :] = Q[nextzero, :]
                Q[nextzero, :] = temp_row
                # Swap columns
                temp_col = Q[:, last].copy()
                Q[:, last] = Q[:, nextzero]
                Q[:, nextzero] = temp_col
                # Swap drift entries
                drifts[last], drifts[nextzero] = drifts[nextzero], drifts[last]

            nindex = len(drifts) - len(zerodrift)
            Qnn = Q[:nindex, :nindex]
            Qzz = Q[nindex:, nindex:]
            Qnz = Q[:nindex, nindex:]
            Qzn = Q[nindex:, :nindex]

            # Schur complement to eliminate zero-drift states
            Qin = Qnn - Qnz @ np.linalg.solve(Qzz, Qzn)
            Rn = np.diag(drifts[:nindex])
            zeroconverter = -Qnz @ inv(Qzz)
        else:
            Qin = Q
            Rn = np.diag(drifts)

        # Compute boundaries for this regime
        upper_B = Bmulti[multi]
        if multi == 0:
            lower_B = 0.0
        else:
            lower_B = Bmulti[multi - 1]

        # A = Qin * inv(Rn)  (MATLAB: A = Qin / Rn)
        A = np.linalg.solve(Rn.T, Qin.T).T

        # Ordered Schur decomposition: [zero | negative | positive]
        Z1, D1, zero_count, negative_count, positive_count = \
            _ordered_schur_3way(A)

        length_D = len(D1)
        neg_index = zero_count
        pos_index = neg_index + negative_count

        # Compute X1: Sylvester equation A0*X1 + X1*(-k2) = -k1
        A0 = D1[:zero_count, :zero_count]
        k1 = D1[:zero_count, neg_index:]
        k2 = D1[neg_index:, neg_index:]

        if zero_count > 0 and k1.size > 0 and k2.size > 0:
            X1 = solve_sylvester(A0, -k2, -k1)
        else:
            X1 = np.zeros((zero_count, length_D - zero_count))

        # Compute X2: Sylvester equation Aneg*X2 + X2*(-Apos) = -Aposneg
        pos_k2_index = negative_count
        Aneg = k2[:negative_count, :negative_count] if negative_count > 0 else np.zeros((0, 0))
        Apos = k2[pos_k2_index:, pos_k2_index:] if positive_count > 0 else np.zeros((0, 0))
        Aposneg = k2[:negative_count, pos_k2_index:] if (negative_count > 0 and positive_count > 0) else np.zeros((negative_count, positive_count))

        if negative_count > 0 and positive_count > 0 and Aposneg.size > 0:
            X2 = solve_sylvester(Aneg, -Apos, -Aposneg)
        else:
            X2 = np.zeros((negative_count, positive_count))

        # Compute Y transformation matrix
        e1 = np.eye(length_D)
        if X1.size > 0 and X1.shape[0] > 0 and X1.shape[1] > 0:
            e1[:X1.shape[0], (length_D - X1.shape[1]):] = X1

        dim_e2 = negative_count + positive_count
        if dim_e2 > 0:
            e2 = np.eye(dim_e2)
            if X2.size > 0 and X2.shape[0] > 0 and X2.shape[1] > 0:
                e2[:X2.shape[0], (dim_e2 - X2.shape[1]):] = X2
        else:
            e2 = np.zeros((0, 0))

        e3 = np.eye(length_D)
        if dim_e2 > 0:
            e3[(length_D - dim_e2):, (length_D - dim_e2):] = e2

        Y = Z1 @ e1 @ e3
        Yinv = inv(Y)

        # Extract L matrices from Yinv
        Lzero = Yinv[:zero_count, :]
        Lneg = Yinv[neg_index:zero_count + negative_count, :]
        Lpos = Yinv[pos_index:, :]

        # If there were zero-drift states, extend L matrices back to full dimension
        if len(zerodrift) > 0:
            Lzero = np.hstack([Lzero, Lzero @ zeroconverter])
            Lpos = np.hstack([Lpos, Lpos @ zeroconverter])
            Lneg = np.hstack([Lneg, Lneg @ zeroconverter])

            # Reverse the column swaps (same order as original swaps)
            for i in range(len(zerodrift)):
                last = len(Q) - 1 - i
                nextzero = zerodrift[i]
                for L in [Lzero, Lneg, Lpos]:
                    temp = L[:, last].copy()
                    L[:, last] = L[:, nextzero]
                    L[:, nextzero] = temp

        # Compute boundary values and integrals when coefficients a=1
        delta = upper_B - lower_B
        I_neg = np.eye(Aneg.shape[0]) if Aneg.shape[0] > 0 else np.zeros((0, 0))
        I_pos = np.eye(Apos.shape[0]) if Apos.shape[0] > 0 else np.zeros((0, 0))

        # Helper to get column count for empty arrays
        def _ncols(M):
            return M.shape[1] if M.ndim >= 2 else 0

        # Integral solution
        if Lzero.shape[0] > 0:
            int_zero = Lzero * delta
        else:
            int_zero = np.zeros((0, _ncols(Lzero)))

        if Aneg.shape[0] > 0:
            int_neg = np.linalg.solve(Aneg, expm(Aneg * delta) - I_neg) @ Lneg
        else:
            int_neg = np.zeros((0, _ncols(Lneg)))

        if Apos.shape[0] > 0:
            int_pos = np.linalg.solve(-Apos, expm(-Apos * delta) - I_pos) @ Lpos
        else:
            int_pos = np.zeros((0, _ncols(Lpos)))

        integral_solution = np.vstack([int_zero, int_neg, int_pos])

        # Lower bound solution: f(lower_B)
        low_zero = Lzero.copy() if Lzero.shape[0] > 0 else np.zeros((0, _ncols(Lzero)))
        low_neg = I_neg @ Lneg if Aneg.shape[0] > 0 else np.zeros((0, _ncols(Lneg)))
        if Apos.shape[0] > 0:
            low_pos = expm(-Apos * delta) @ Lpos
        else:
            low_pos = np.zeros((0, _ncols(Lpos)))
        low_solution = np.vstack([low_zero, low_neg, low_pos])

        # Upper bound solution: f(upper_B)
        up_zero = Lzero.copy() if Lzero.shape[0] > 0 else np.zeros((0, _ncols(Lzero)))
        if Aneg.shape[0] > 0:
            up_neg = expm(Aneg * delta) @ Lneg
        else:
            up_neg = np.zeros((0, _ncols(Lneg)))
        up_pos = I_pos @ Lpos if Apos.shape[0] > 0 else np.zeros((0, _ncols(Lpos)))
        up_solution = np.vstack([up_zero, up_neg, up_pos])

        lower_bound_solutions[multi] = low_solution
        upper_bound_solutions[multi] = up_solution
        integral_solutions[multi] = integral_solution
        Lzero_multi[multi] = Lzero
        Lpos_multi[multi] = Lpos
        Lneg_multi[multi] = Lneg
        Aneg_multi[multi] = Aneg
        Apos_multi[multi] = Apos

    return (lower_bound_solutions, upper_bound_solutions, integral_solutions,
            Lzero_multi, Lneg_multi, Lpos_multi, Aneg_multi, Apos_multi)
