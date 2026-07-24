"""
Multi-Regime Markov Fluid Queue (MRMFQ) Solver.

Port of MRMFQSolver.m from the MAPMsG MATLAB library.

Solves the boundary value problem for multi-regime Markov fluid queues by
assembling the global constraint matrix H, computing the null-space vector,
and normalizing to obtain probability coefficients and boundary masses.

Reference:
    O. Gursoy, K. A. Mehr, N. Akar, "The MAP/M/s + G Call Center Model with
    Generally Distributed Patience Times."
"""

import numpy as np
from numpy.linalg import inv

from .additive_decomposition import additive_decomposition


def mrmfq_solver(Qregimes, Qboundaries, driftregimes, driftboundaries, B):
    """Solve a multi-regime Markov fluid queue.

    Parameters
    ----------
    Qregimes : ndarray, shape (m, m, num_regimes)
        Generator matrices for each fluid regime.
    Qboundaries : ndarray, shape (m, m, num_boundaries)
        Generator matrices at each boundary level.
    driftregimes : ndarray, shape (num_regimes, m)
        Drift rates for each regime and state.
    driftboundaries : ndarray, shape (num_boundaries, m)
        Drift rates at each boundary and state.
    B : ndarray, shape (num_regimes,)
        Upper boundary levels for each regime.

    Returns
    -------
    coefficients : list of ndarray
        Probability density coefficients for each regime (0..num_regimes-2).
    boundaries : list of ndarray
        Boundary probability masses at each boundary level.
    Lzero_multi : list of ndarray
        Zero-eigenvalue projection matrices from additive decomposition.
    Lneg_multi : list of ndarray
        Negative-eigenvalue projection matrices from additive decomposition.
    Lpos_multi : list of ndarray
        Positive-eigenvalue projection matrices from additive decomposition.
    Aneg_multi : list of ndarray
        Negative eigenvalue block matrices from additive decomposition.
    Apos_multi : list of ndarray
        Positive eigenvalue block matrices from additive decomposition.
    """
    num_boundaries = driftboundaries.shape[0]
    num_states = driftboundaries.shape[1]
    num_regimes = driftregimes.shape[0]

    # Build R (drift) matrices for each regime
    Rregimes = np.zeros_like(Qregimes)
    for k in range(num_regimes):
        Rregimes[:, :, k] = np.diag(driftregimes[k, :])

    # Determine which boundary masses are non-zero
    nonzero_boundaries = np.ones_like(driftboundaries)
    zerolower_pdf = np.zeros_like(driftregimes)
    zeroupper_pdf = np.zeros_like(driftregimes)

    for level in range(num_boundaries):
        for state in range(num_states):
            if level == 0:
                # First boundary: zero mass if drift is positive (flow into buffer)
                if driftregimes[level, state] > 0:
                    nonzero_boundaries[level, state] = 0
            elif level == num_boundaries - 1:
                # Last boundary: zero mass if drift of previous regime is negative
                if driftregimes[level - 1, state] < 0:
                    nonzero_boundaries[level, state] = 0
            else:
                # Interior boundaries
                dr_cur = driftregimes[level, state]
                dr_prev = driftregimes[level - 1, state]
                if ((dr_cur > 0 and dr_prev > 0)
                        or (dr_cur < 0 and dr_prev < 0)
                        or (dr_cur > 0 and dr_prev < 0
                            and driftboundaries[level, state] != 0)):
                    nonzero_boundaries[level, state] = 0
                if dr_cur > 0 and driftboundaries[level, state] <= 0:
                    zerolower_pdf[level, state] = 1
                if dr_prev < 0 and driftboundaries[level, state] >= 0:
                    zeroupper_pdf[level - 1, state] = 1

    # Additive decomposition
    (lower_bound_solutions, upper_bound_solutions, integral_solutions,
     Lzero_multi, Lneg_multi, Lpos_multi, Aneg_multi, Apos_multi) = \
        additive_decomposition(Qregimes, driftregimes, B)

    # Build global constraint matrix H
    # We need to figure out total size first by doing a dry run
    # Count total rows and columns
    total_rows = 0
    total_cols = 0
    for funcindex in range(num_boundaries):
        nonzero_indices = np.where(nonzero_boundaries[funcindex, :] == 0)[0]
        Qelim = np.delete(Qboundaries[:, :, funcindex], nonzero_indices, axis=0)
        count_nonzero = int(np.sum(nonzero_boundaries[funcindex, :] == 1))

        if funcindex == 0:
            total_rows += Qelim.shape[0]
            total_rows += lower_bound_solutions[funcindex].shape[0]
            total_cols += Rregimes[:, :, funcindex].shape[1]
        elif funcindex == num_boundaries - 1:
            total_rows += upper_bound_solutions[funcindex - 1].shape[0]
            total_rows += Qelim.shape[0]
        else:
            # zeroupper columns
            for st in range(num_states):
                if zeroupper_pdf[funcindex - 1, st] == 1:
                    total_cols += 1
            total_rows += upper_bound_solutions[funcindex - 1].shape[0]
            total_rows += Qelim.shape[0]
            total_rows += lower_bound_solutions[funcindex].shape[0]
            total_cols += Rregimes[:, :, funcindex].shape[1]
            # zerolower columns
            for st in range(num_states):
                if zerolower_pdf[funcindex, st] == 1:
                    total_cols += 1

    # Actually, let's just build H incrementally (matching MATLAB logic exactly)
    # Use a list-of-lists approach, then convert to array at the end.
    # MATLAB builds H by dynamically growing it with block assignments.

    # First pass: determine dimensions
    row_sizes = []
    col_sizes = []
    for funcindex in range(num_boundaries):
        nonzero_indices = np.where(nonzero_boundaries[funcindex, :] == 0)[0]
        Qelim = np.delete(Qboundaries[:, :, funcindex], nonzero_indices, axis=0)

        if funcindex == 0:
            row_sizes.append(Qelim.shape[0])
            row_sizes.append(lower_bound_solutions[funcindex].shape[0])
            col_sizes.append(Rregimes[:, :, funcindex].shape[1])
        elif funcindex == num_boundaries - 1:
            row_sizes.append(upper_bound_solutions[funcindex - 1].shape[0])
            row_sizes.append(Qelim.shape[0])
        else:
            for st in range(num_states):
                if zeroupper_pdf[funcindex - 1, st] == 1:
                    col_sizes.append(1)
            row_sizes.append(upper_bound_solutions[funcindex - 1].shape[0])
            row_sizes.append(Qelim.shape[0])
            row_sizes.append(lower_bound_solutions[funcindex].shape[0])
            col_sizes.append(Rregimes[:, :, funcindex].shape[1])
            for st in range(num_states):
                if zerolower_pdf[funcindex, st] == 1:
                    col_sizes.append(1)

    total_r = sum(row_sizes)
    total_c = sum(col_sizes)

    # Ensure H is large enough (MATLAB grows dynamically)
    H = np.zeros((max(total_r, total_c + 10), max(total_r, total_c + 10)))

    rowindex = 0
    columnindex = 0

    for funcindex in range(num_boundaries):
        nonzero_indices = np.where(nonzero_boundaries[funcindex, :] == 0)[0]
        Qeliminated = np.delete(Qboundaries[:, :, funcindex], nonzero_indices, axis=0)

        if funcindex == 0:
            nrows_q = Qeliminated.shape[0]
            ncols_r = Rregimes[:, :, funcindex].shape[1]
            H[rowindex:rowindex + nrows_q, columnindex:columnindex + ncols_r] = -Qeliminated
            rowindex += nrows_q
            nrows_lb = lower_bound_solutions[funcindex].shape[0]
            H[rowindex:rowindex + nrows_lb, columnindex:columnindex + ncols_r] = \
                lower_bound_solutions[funcindex] @ Rregimes[:, :, funcindex]
            columnindex += ncols_r

        elif funcindex == num_boundaries - 1:
            ncols_r_prev = Rregimes[:, :, funcindex - 1].shape[1]
            nrows_ub = upper_bound_solutions[funcindex - 1].shape[0]
            H[rowindex:rowindex + nrows_ub,
              columnindex:columnindex + ncols_r_prev] = \
                upper_bound_solutions[funcindex - 1] @ Rregimes[:, :, funcindex - 1]
            rowindex += nrows_ub
            nrows_q = Qeliminated.shape[0]
            H[rowindex:rowindex + nrows_q,
              columnindex:columnindex + ncols_r_prev] = Qeliminated

        else:
            # Zero upper pdf columns
            for st in range(num_states):
                if zeroupper_pdf[funcindex - 1, st] == 1:
                    temp = upper_bound_solutions[funcindex - 1]
                    nrows_ub = temp.shape[0]
                    H[rowindex:rowindex + nrows_ub, columnindex] = temp[:, st]
                    columnindex += 1

            ncols_r_prev = Rregimes[:, :, funcindex - 1].shape[1]
            nrows_ub = upper_bound_solutions[funcindex - 1].shape[0]
            H[rowindex:rowindex + nrows_ub,
              columnindex:columnindex + ncols_r_prev] = \
                upper_bound_solutions[funcindex - 1] @ Rregimes[:, :, funcindex - 1]
            rowindex += nrows_ub

            nrows_q = Qeliminated.shape[0]
            ncols_r_cur = Rregimes[:, :, funcindex].shape[1]
            H[rowindex:rowindex + nrows_q,
              columnindex:columnindex + ncols_r_cur] = Qeliminated
            rowindex += nrows_q

            nrows_lb = lower_bound_solutions[funcindex].shape[0]
            H[rowindex:rowindex + nrows_lb,
              columnindex:columnindex + ncols_r_cur] = \
                -lower_bound_solutions[funcindex] @ Rregimes[:, :, funcindex]
            columnindex += ncols_r_cur

            # Zero lower pdf columns
            for st in range(num_states):
                if zerolower_pdf[funcindex, st] == 1:
                    temp2 = lower_bound_solutions[funcindex]
                    nrows_lb2 = temp2.shape[0]
                    H[rowindex:rowindex + nrows_lb2, columnindex] = temp2[:, st]
                    columnindex += 1

    # Trim H to actual used dimensions
    actual_rows = rowindex + (Qeliminated.shape[0] if funcindex == num_boundaries - 1
                              else lower_bound_solutions[funcindex].shape[0])
    # Actually, MATLAB doesn't trim - it uses whatever H is.
    # The key size is determined by the last rowindex + last block added.
    # Let's figure out the actual dimensions from the MATLAB logic.
    # After the loop, the matrix is accessed as H with implicit sizing.
    # We use rowindex as last row start + last block size for final row count.
    # For the last boundary (funcindex == num_boundaries - 1):
    #   rowindex was incremented by nrows_ub, then nrows_q was assigned
    #   final row = rowindex + nrows_q
    # For interior boundaries, rowindex already advanced past the last block.
    # The column dimension is columnindex (which may differ from row dimension).

    # Determine actual H dimensions
    if num_boundaries > 1 and funcindex == num_boundaries - 1:
        final_row = rowindex + Qeliminated.shape[0]
    else:
        final_row = rowindex

    final_col = columnindex
    if funcindex == num_boundaries - 1:
        final_col = columnindex + Rregimes[:, :, funcindex - 1].shape[1]

    # Use the maximum of computed dimensions
    dim = max(final_row, final_col)
    H = H[:dim, :dim]

    # Replace first column with ones for normalization, solve for z
    Hbar = H.copy()
    Hbar[:, 0] = 1.0

    # z = e1^T / Hbar  (MATLAB: z = eye(1,size(H,2))/Hbar)
    e1_vec = np.zeros(dim)
    e1_vec[0] = 1.0
    z = np.linalg.solve(Hbar.T, e1_vec)

    # Compute normalization constant
    index = 0
    normalization_coef_masses = 0.0
    normalization_coef_integrals = 0.0
    for level in range(num_boundaries):
        count_nonzero = int(np.sum(nonzero_boundaries[level, :] == 1))
        if count_nonzero > 0:
            normalization_coef_masses += np.sum(z[index:index + count_nonzero])
            index += count_nonzero
        if level < num_boundaries - 1:
            n_integral = integral_solutions[level].shape[0]
            normalization_coef_integrals += np.sum(
                z[index:index + n_integral] @ integral_solutions[level])
            index += n_integral

    normalization_coef = normalization_coef_integrals + normalization_coef_masses
    final_z = z / normalization_coef

    # Extract coefficients and boundaries
    coefficients = [None] * num_regimes
    boundaries_out = [None] * num_boundaries

    index2 = 0
    for level in range(num_boundaries):
        count_nonzero = int(np.sum(nonzero_boundaries[level, :] == 1))
        boundaries_out[level] = final_z[index2:index2 + count_nonzero]
        index2 += count_nonzero
        if level < num_boundaries - 1:
            n_integral = integral_solutions[level].shape[0]
            coefficients[level] = final_z[index2:index2 + n_integral]
            index2 += n_integral

    # Expand boundary masses back into full state vectors
    for level in range(num_boundaries):
        l_idx = 0
        full_boundary = np.zeros(num_states)
        for k in range(num_states):
            if nonzero_boundaries[level, k] == 1:
                bound = boundaries_out[level]
                full_boundary[k] = bound[l_idx]
                l_idx += 1
            else:
                full_boundary[k] = 0.0
        boundaries_out[level] = full_boundary

    return (coefficients, boundaries_out, Lzero_multi, Lneg_multi,
            Lpos_multi, Aneg_multi, Apos_multi)
