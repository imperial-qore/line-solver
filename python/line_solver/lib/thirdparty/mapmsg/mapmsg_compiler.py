"""
MAP/M/s+G Compiler.

Port of MAPMsGCompiler.m from the MAPMsG MATLAB library.

Implements steady-state and first-passage-time analysis for MAP/M/s+G queues
(call center models with MAP arrivals, exponential service, s servers, and
generally distributed patience times).

Reference:
    O. Gursoy, K. A. Mehr, N. Akar, "The MAP/M/s + G Call Center Model with
    Generally Distributed Patience Times."
"""

import numpy as np
from scipy.linalg import expm

from .cme_parameter_calculator import cme_parameter_calculator
from .mrmfq_solver import mrmfq_solver


# Solution type constants
STEADY_STATE = 1
FIRST_PASSAGE_VIRTUAL = 2
FIRST_PASSAGE_ACTUAL = 3


class MAPMsGResult:
    """Container for MAP/M/s+G computation results.

    Attributes
    ----------
    steady_state : ndarray or None
        7-element vector [Pr{W=0}, Pr{W=0|S}, Pr{A}, E{W|S}, Var{W|S},
        F_{w|s,w>0}(0.1), F_{w|s,w>0}(0.2)] for steady-state solution.
    first_passage_virtual : float or None
        First passage time probability for virtual waiting time.
    first_passage_actual : float or None
        First passage time probability for actual waiting time.
    coefficients : list
        Probability density coefficients from MRMFQ solver.
    boundaries : list
        Boundary probability masses from MRMFQ solver.
    Lzero_multi : list
        Zero-eigenvalue projection matrices.
    Lneg_multi : list
        Negative-eigenvalue projection matrices.
    Lpos_multi : list
        Positive-eigenvalue projection matrices.
    Aneg_multi : list
        Negative eigenvalue block matrices.
    Apos_multi : list
        Positive eigenvalue block matrices.
    boundary_levels : ndarray
        Boundary levels used in the computation.
    """

    def __init__(self):
        self.steady_state = None
        self.first_passage_virtual = None
        self.first_passage_actual = None
        self.coefficients = None
        self.boundaries = None
        self.Lzero_multi = None
        self.Lneg_multi = None
        self.Lpos_multi = None
        self.Aneg_multi = None
        self.Apos_multi = None
        self.boundary_levels = None


def _build_generators(server_size, map_size, mu, C, D, ga, quantization):
    """Build infinitesimal generator and drift matrices for the QBD-like structure.

    Parameters
    ----------
    server_size : int
        Number of servers (s).
    map_size : int
        Order of the MAP(C,D) arrival process.
    mu : float
        Service rate.
    C : ndarray, shape (map_size, map_size)
        MAP D0 matrix.
    D : ndarray, shape (map_size, map_size)
        MAP D1 matrix.
    ga : ndarray, shape (quantization + 2,)
        Abandonment probabilities for each regime (starts with 0).
    quantization : int
        Number of regimes.

    Returns
    -------
    Qy : ndarray, shape (state_size, state_size, quantization)
        Generator matrices for regimes.
    Qy0 : ndarray, shape (state_size, state_size)
        Generator matrix for the boundary at level 0.
    Ry : ndarray, shape (state_size, state_size)
        Base drift matrix.
    Rydiag : ndarray, shape (state_size,)
        Diagonal of the drift matrix.
    ydriftregimes : ndarray, shape (quantization, state_size)
        Drift rates for each regime.
    Ryregimes : ndarray, shape (state_size, state_size, quantization)
        Drift matrices for each regime.
    """
    I = np.eye(map_size)
    state_size = (server_size + 1) * map_size

    Qy0 = np.zeros((state_size, state_size))
    Rydiag = -np.ones(state_size)
    # Last map_size states have drift +1 (server occupancy = server_size)
    Rydiag[state_size - map_size:] = 1.0
    Ry = np.diag(Rydiag)

    # Build Qy0: boundary generator at level 0
    for row in range(1, server_size + 2):  # 1-indexed rows
        r0 = (row - 1) * map_size
        r1 = r0 + map_size
        if row == 1:
            Qy0[r0:r1, r0:r1] = C
            Qy0[r0:r1, r1:r1 + map_size] = D
        elif row == server_size + 1:
            Qy0[r0:r1, r0:r1] = -(row - 1) * mu * I
            Qy0[r0:r1, r0 - map_size:r0] = (row - 1) * mu * I
        else:
            Qy0[r0:r1, r0:r1] = C - (row - 1) * mu * I
            Qy0[r0:r1, r0 - map_size:r0] = (row - 1) * mu * I
            Qy0[r0:r1, r1:r1 + map_size] = D

    # Build Qy: regime generators (only rows for server_size and server_size+1 differ)
    Qy = np.zeros((state_size, state_size, quantization))
    ydriftregimes = np.zeros((quantization, state_size))
    Ryregimes = np.zeros((state_size, state_size, quantization))

    for regimecount in range(quantization):
        ga_val = ga[regimecount + 1]  # ga is 0-indexed, ga[0]=0, ga[1]=first regime, etc.
        for row_idx in [server_size, server_size + 1]:  # 1-indexed
            r0 = (row_idx - 1) * map_size
            r1 = r0 + map_size
            if row_idx == server_size:
                Qy[r0:r1, r0:r1, regimecount] = ga_val * D + C
                Qy[r0:r1, r1:r1 + map_size, regimecount] = (1 - ga_val) * D
            elif row_idx == server_size + 1:
                Qy[r0:r1, r0:r1, regimecount] = -(row_idx - 1) * mu * I
                Qy[r0:r1, r0 - map_size:r0, regimecount] = (row_idx - 1) * mu * I
        ydriftregimes[regimecount, :] = Rydiag
        Ryregimes[:, :, regimecount] = Ry

    return Qy, Qy0, Ry, Rydiag, ydriftregimes, Ryregimes


def solve_steady_state(server_size, map_size, mu, C, D, ga, boundary_levels,
                       quantization):
    """Compute steady-state performance metrics for a MAP/M/s+G queue.

    Parameters
    ----------
    server_size : int
        Number of servers (s).
    map_size : int
        Order of the MAP(C,D) arrival process.
    mu : float
        Service rate.
    C : ndarray, shape (map_size, map_size)
        MAP D0 matrix.
    D : ndarray, shape (map_size, map_size)
        MAP D1 matrix.
    ga : ndarray
        Abandonment probabilities for each regime.  Must start with 0 at the
        beginning.  Length should be quantization + 2.
    boundary_levels : ndarray
        Boundary levels of regimes (not including level 0).  The last element
        should be large (e.g. 10^7) to represent infinity.
    quantization : int
        Number of regimes that quantize the abandonment function.

    Returns
    -------
    result : MAPMsGResult
        Result object with ``steady_state`` attribute containing the 7-element
        vector [Pr{W=0}, Pr{W=0|S}, Pr{A}, E{W|S}, Var{W|S},
        F_{w|s,w>0}(0.1), F_{w|s,w>0}(0.2)].
    """
    (Qy, Qy0, Ry, Rydiag, ydriftregimes, Ryregimes) = \
        _build_generators(server_size, map_size, mu, C, D, ga, quantization)

    # Assemble boundary and regime arrays
    Qybounds = np.concatenate([Qy0[:, :, np.newaxis], Qy], axis=2)
    ydriftbounds = np.vstack([Rydiag, ydriftregimes])

    Qregimes = Qy
    Qbounds = Qybounds
    driftregimes = ydriftregimes
    driftbound = ydriftbounds
    B = np.asarray(boundary_levels)

    # Solve MRMFQ
    (coefficients, boundaries, Lzero_multi, Lneg_multi, Lpos_multi,
     Aneg_multi, Apos_multi) = \
        mrmfq_solver(Qregimes, Qbounds, driftregimes, driftbound, B)

    # Compute steady-state metrics
    em = np.ones((map_size, 1))
    lmap = (D @ em).flatten()

    zeromass = boundaries[0]

    # Compute integrals for each regime
    integral = np.zeros((quantization, len(zeromass)))
    waitintegral = np.zeros_like(integral)
    abandonintegral = np.zeros_like(integral)
    successfulintegral = np.zeros_like(integral)

    for d in range(quantization):
        if d == 0:
            lower = 0.0
        else:
            lower = B[d - 1]
        upper = B[d]
        delta = upper - lower

        Aneg = Aneg_multi[d]
        Apos = Apos_multi[d]
        Lzero = Lzero_multi[d]
        Lneg = Lneg_multi[d]
        Lpos = Lpos_multi[d]

        I_neg = np.eye(Aneg.shape[0]) if Aneg.shape[0] > 0 else np.zeros((0, 0))
        I_pos = np.eye(Apos.shape[0]) if Apos.shape[0] > 0 else np.zeros((0, 0))

        # Integral term
        int_zero = Lzero * delta
        if Aneg.shape[0] > 0:
            int_neg = np.linalg.solve(Aneg, expm(Aneg * delta) - I_neg) @ Lneg
        else:
            int_neg = np.zeros((0, Lneg.shape[1]))
        if Apos.shape[0] > 0:
            int_pos = np.linalg.solve(-Apos, expm(-Apos * delta) - I_pos) @ Lpos
        else:
            int_pos = np.zeros((0, Lpos.shape[1]))
        int_block = np.vstack([int_zero, int_neg, int_pos])

        coef = coefficients[d]
        integral[d, :] = coef @ int_block

        ga_val = ga[d + 1]  # regime d uses ga[d+1] (shifted because ga starts with 0)
        # For steady state, ga indexing is ga[d+1] in MATLAB (1-indexed),
        # which maps to ga[d+1] in Python (0-indexed) because ga has extra leading 0
        # Actually MATLAB uses ga(d+1+1) = ga(d+2) for d=1-indexed from 1 to QUANTIZATION
        # In the MATLAB code: ga(2) for d=1 (first regime)
        # ga array in MATLAB: [0, ga_values...]
        # so for regime d (0-indexed), MATLAB uses ga(d+2) which is ga[d+2] (0-indexed)
        # But looking at the MATLAB code more carefully:
        # ga = [0 ga] after construction, so ga has QUANTIZATION+2 elements
        # For d=1 (first regime, 1-indexed): ga(2) -> ga[1] in 0-indexed
        # Wait - let me re-read the MATLAB code:
        # ga=[0 ga]; adds a 0 at the front
        # Then in the regime loop: ga(regimecount+1) for regimecount=1:QUANTIZATION
        # So ga(2) for first regime, ga(3) for second, etc.
        # The integral uses ga(d+1) for d=1:QUANTIZATION (1-indexed)
        # = ga(2), ga(3), ... = ga[1], ga[2], ... in 0-indexed

        # For first regime (d=0 in 0-indexed): uses ga(2) in MATLAB = ga[1] in python
        # but we already set ga_val = ga[d+1] above

        # Wait, the MATLAB code actually uses:
        # For d=1 (1-indexed): ga(d+1) = ga(2)
        # For d=2: ga(d+1) = ga(3)
        # So for Python d=0: ga[1], d=1: ga[2], etc.
        # But the original MATLAB indexes into ga as ga(2) for first loop iteration
        # when abandonment_integral and successfulintegral use ga(d+1) with d=1:QUANTIZATION

        # Actually re-reading the MATLAB code carefully:
        # abandonintegral(1,:) = (ga(2)) * coefficients{1} * [...]
        # abandonintegral(d,:) = (ga(d+1)) * coefficients{d} * [...] for d=2:QUANTIZATION
        # So for d (1-indexed): ga(d+1)
        # In Python (0-indexed d): ga[d+1]
        # BUT the ga array was constructed as:
        # ga indexing: MATLAB uses ga(d+1) for d=1:QUANTIZATION (1-indexed).
        # In Python (0-indexed d): ga[d+1].
        # ga was constructed as [0, ga_raw] where ga_raw has QUANTIZATION+1 elements,
        # giving ga QUANTIZATION+2 total.  ga[0]=0, ga[1]=first regime value, etc.
        if d == 0:
            mid = B[0] / 2.0
        else:
            mid = (B[d] + B[d - 1]) / 2.0

        ga_regime = ga[d + 1]  # MATLAB ga(d+1) with d 1-indexed = ga[d+1-1+1] = ga[d+1]

        waitintegral[d, :] = mid * (1 - ga_regime) * coef @ int_block
        abandonintegral[d, :] = ga_regime * coef @ int_block
        successfulintegral[d, :] = (1 - ga_regime) * coef @ int_block

    # CDF at thresholds x = [0.1, 0.2]
    x_vals = [0.1, 0.2]
    normalization_size = server_size * map_size
    arrivals_wait_less_than_x = np.zeros((len(x_vals), len(zeromass)))

    for d_idx, x in enumerate(x_vals):
        count = 0
        while B[count] < x:
            count += 1
        if count > 0:
            lower_b = B[count - 1]
            ga_regime = ga[count + 1]
            coef_c = coefficients[count]
            Aneg_c = Aneg_multi[count]
            Apos_c = Apos_multi[count]
            Lzero_c = Lzero_multi[count]
            Lneg_c = Lneg_multi[count]
            Lpos_c = Lpos_multi[count]
            I_neg_c = np.eye(Aneg_c.shape[0]) if Aneg_c.shape[0] > 0 else np.zeros((0, 0))
            I_pos_c = np.eye(Apos_c.shape[0]) if Apos_c.shape[0] > 0 else np.zeros((0, 0))

            int_zero_c = Lzero_c * (x - lower_b)
            if Aneg_c.shape[0] > 0:
                int_neg_c = np.linalg.solve(Aneg_c, expm(Aneg_c * (x - lower_b)) - I_neg_c) @ Lneg_c
            else:
                int_neg_c = np.zeros((0, Lneg_c.shape[1]))
            if Apos_c.shape[0] > 0:
                int_pos_c = (np.linalg.solve(-Apos_c,
                    expm(-Apos_c * (B[count] - x)) - expm(-Apos_c * (B[count] - lower_b)))) @ Lpos_c
            else:
                int_pos_c = np.zeros((0, Lpos_c.shape[1]))
            int_block_c = np.vstack([int_zero_c, int_neg_c, int_pos_c])
            int_val = (1 - ga_regime) * coef_c @ int_block_c
        else:
            ga_regime = ga[count + 1]
            coef_c = coefficients[count]
            Aneg_c = Aneg_multi[count]
            Apos_c = Apos_multi[count]
            Lzero_c = Lzero_multi[count]
            Lneg_c = Lneg_multi[count]
            Lpos_c = Lpos_multi[count]
            I_neg_c = np.eye(Aneg_c.shape[0]) if Aneg_c.shape[0] > 0 else np.zeros((0, 0))
            I_pos_c = np.eye(Apos_c.shape[0]) if Apos_c.shape[0] > 0 else np.zeros((0, 0))

            int_zero_c = Lzero_c * x
            if Aneg_c.shape[0] > 0:
                int_neg_c = np.linalg.solve(Aneg_c, expm(Aneg_c * x) - I_neg_c) @ Lneg_c
            else:
                int_neg_c = np.zeros((0, Lneg_c.shape[1]))
            if Apos_c.shape[0] > 0:
                int_pos_c = (np.linalg.solve(-Apos_c,
                    expm(-Apos_c * (B[count] - x)) - expm(-Apos_c * B[count]))) @ Lpos_c
            else:
                int_pos_c = np.zeros((0, Lpos_c.shape[1]))
            int_block_c = np.vstack([int_zero_c, int_neg_c, int_pos_c])
            int_val = (1 - ga_regime) * coef_c @ int_block_c

        for l_idx in range(1, count + 1):
            arrivals_wait_less_than_x[d_idx, :] += successfulintegral[l_idx - 1, :]
        arrivals_wait_less_than_x[d_idx, :] += int_val

    # Map to arrival rates
    arrivals_mapped = np.zeros_like(arrivals_wait_less_than_x)
    integral_mapped = np.zeros_like(integral)
    abandon_mapped = np.zeros_like(abandonintegral)
    wait_mapped = np.zeros_like(waitintegral)
    zeromass_mapped = np.zeros_like(zeromass)

    for r in range(len(zeromass)):
        lmap_idx = r % map_size
        arrivals_mapped[:, r] = arrivals_wait_less_than_x[:, r] * lmap[lmap_idx]
        integral_mapped[:, r] = integral[:, r] * lmap[lmap_idx]
        abandon_mapped[:, r] = abandonintegral[:, r] * lmap[lmap_idx]
        wait_mapped[:, r] = waitintegral[:, r] * lmap[lmap_idx]
        zeromass_mapped[r] = zeromass[r] * lmap[lmap_idx]

    normalization = normalization_size
    total_integral = np.sum(integral_mapped[:, :normalization])
    total_zeromass = np.sum(zeromass_mapped)
    total = total_integral + total_zeromass

    w0 = total_zeromass / total
    abandon_prob = np.sum(abandon_mapped[:, :normalization]) / total
    w0s = total_zeromass / (total * (1 - abandon_prob))
    expected_wait = np.sum(wait_mapped[:, :normalization]) / (total * (1 - abandon_prob))

    prob_01 = np.sum(arrivals_mapped[0, :normalization]) / \
        (total * (1 - abandon_prob) * (1 - w0s))
    prob_02 = np.sum(arrivals_mapped[1, :normalization]) / \
        (total * (1 - abandon_prob) * (1 - w0s))

    # Variance computation
    variance_integral = np.zeros_like(integral)
    for d in range(quantization):
        if d == 0:
            mid = B[0] / 2.0
            lower = 0.0
        else:
            mid = (B[d] + B[d - 1]) / 2.0
            lower = B[d - 1]
        upper = B[d]
        delta = upper - lower
        ga_regime = ga[d + 1]

        Aneg = Aneg_multi[d]
        Apos = Apos_multi[d]
        Lzero = Lzero_multi[d]
        Lneg = Lneg_multi[d]
        Lpos = Lpos_multi[d]
        I_neg = np.eye(Aneg.shape[0]) if Aneg.shape[0] > 0 else np.zeros((0, 0))
        I_pos = np.eye(Apos.shape[0]) if Apos.shape[0] > 0 else np.zeros((0, 0))

        int_zero = Lzero * delta
        if Aneg.shape[0] > 0:
            int_neg = np.linalg.solve(Aneg, expm(Aneg * delta) - I_neg) @ Lneg
        else:
            int_neg = np.zeros((0, Lneg.shape[1]))
        if Apos.shape[0] > 0:
            int_pos = np.linalg.solve(-Apos, expm(-Apos * delta) - I_pos) @ Lpos
        else:
            int_pos = np.zeros((0, Lpos.shape[1]))
        int_block = np.vstack([int_zero, int_neg, int_pos])

        variance_integral[d, :] = (mid - expected_wait) ** 2 * (1 - ga_regime) * \
            coefficients[d] @ int_block

    variance_integral_mapped = np.zeros_like(variance_integral)
    for r in range(len(zeromass)):
        lmap_idx = r % map_size
        variance_integral_mapped[:, r] = variance_integral[:, r] * lmap[lmap_idx]

    variance = (np.sum(variance_integral_mapped[:, :normalization]) +
                total_zeromass * expected_wait ** 2) / (total * (1 - abandon_prob))

    result = MAPMsGResult()
    result.steady_state = np.array([
        w0, w0s, abandon_prob, expected_wait, variance, prob_01, prob_02
    ])
    result.coefficients = coefficients
    result.boundaries = boundaries
    result.Lzero_multi = Lzero_multi
    result.Lneg_multi = Lneg_multi
    result.Lpos_multi = Lpos_multi
    result.Aneg_multi = Aneg_multi
    result.Apos_multi = Apos_multi
    result.boundary_levels = B
    return result


def _build_first_passage_virtual(server_size, map_size, mu, C, D, ga,
                                 boundary_levels, quantization, b, tau,
                                 use_cme, order_of_ph, pi0, theta0):
    """Build augmented generators for virtual first-passage-time analysis.

    Returns Qregimes, Qbounds, driftregimes, driftbound, BoundaryLevelsLast.
    """
    (Qy, Qy0, Ry, Rydiag, ydriftregimes, Ryregimes) = \
        _build_generators(server_size, map_size, mu, C, D, ga, quantization)

    Qybounds = np.concatenate([Qy0[:, :, np.newaxis], Qy], axis=2)
    Rybounds = np.concatenate([Ry[:, :, np.newaxis], Ryregimes], axis=2)

    B = np.asarray(boundary_levels)

    # Determine PH representation
    if use_cme:
        me_system, _ = cme_parameter_calculator(order_of_ph, tau)
        S = me_system.A
        S0 = me_system.B
        alpha = me_system.C.flatten()
    else:
        alpha = np.zeros(order_of_ph)
        alpha[0] = 1.0
        ss_diag = -np.ones(order_of_ph)
        S = np.diag(ss_diag)
        for y in range(order_of_ph - 1):
            S[y, y + 1] = 1.0
        S = order_of_ph * S / tau
        e = np.ones((len(S), 1))
        S0 = -S @ e

    # Find boundary level containing b
    count = 0
    while B[count] < b:
        count += 1
    BoundaryLevelsLast = B[:count + 1].copy()
    BoundaryLevelsLast[count] = b

    n_ph = len(S)
    pi0_arr = np.asarray(pi0).flatten()
    theta0_arr = np.asarray(theta0).flatten()
    kron_init = np.kron(alpha, np.kron(pi0_arr, theta0_arr))
    n_aug = 1 + n_ph * Qy.shape[0]

    Qz = np.zeros((n_aug, n_aug, len(BoundaryLevelsLast)))
    Qzbounds = np.zeros((n_aug, n_aug, len(BoundaryLevelsLast) + 1))
    Rz = np.zeros((n_aug, n_aug, len(BoundaryLevelsLast)))
    Rzbounds = np.zeros((n_aug, n_aug, len(BoundaryLevelsLast) + 1))

    for reg in range(len(BoundaryLevelsLast)):
        eTilde = np.ones(Qybounds.shape[0])
        eTilde[-map_size:] = 0.0
        ITilde = np.diag(eTilde)

        # Qz regime
        Qz_block = np.kron(np.eye(n_ph), Qy[:, :, reg]) + np.kron(S, ITilde)
        Qz[0, 0, reg] = 0.0
        Qz[0, 1:, reg] = 0.0
        Qz[1:, 0, reg] = (np.kron(S0.flatten(), np.diag(ITilde)))
        Qz[1:, 1:, reg] = Qz_block

        # Qzbounds (reg+1)
        Qzb_block = np.kron(np.eye(n_ph), Qybounds[:, :, reg + 1]) + np.kron(S, ITilde)
        Qzbounds[0, 0, reg + 1] = 0.0
        Qzbounds[0, 1:, reg + 1] = 0.0
        Qzbounds[1:, 0, reg + 1] = (np.kron(S0.flatten(), np.diag(ITilde)))
        Qzbounds[1:, 1:, reg + 1] = Qzb_block

        # Drift matrices
        Ry_reg_diag = np.diag(Ryregimes[:, :, reg])
        Rz_diag = np.concatenate([[-1.0], np.kron(np.ones(n_ph), Ry_reg_diag)])
        Rz[:, :, reg] = np.diag(Rz_diag)

        Ryb_diag = np.diag(Rybounds[:, :, reg + 1])
        Rzb_diag = np.concatenate([[-1.0], np.kron(np.ones(n_ph), Ryb_diag)])
        Rzbounds[:, :, reg + 1] = np.diag(Rzb_diag)

    # Boundary at level 0 (special: has kron(alpha, kron(pi0, theta0)))
    eTilde = np.ones(Qybounds.shape[0])
    eTilde[-map_size:] = 0.0
    ITilde = np.diag(eTilde)

    Qzb0_block = np.kron(np.eye(n_ph), Qybounds[:, :, 0]) + np.kron(S, ITilde)
    Qzbounds[0, 0, 0] = -1.0
    Qzbounds[0, 1:, 0] = kron_init
    Qzbounds[1:, 0, 0] = (np.kron(S0.flatten(), np.diag(ITilde)))
    Qzbounds[1:, 1:, 0] = Qzb0_block

    Ryb0_diag = np.diag(Rybounds[:, :, 0])
    Rzb0_diag = np.concatenate([[-1.0], np.kron(np.ones(n_ph), Ryb0_diag)])
    Rzbounds[:, :, 0] = np.diag(Rzb0_diag)

    # Last boundary: absorbing
    Qzbounds[0, :, -1] = 0.0
    Qzbounds[0, 0, -1] = 0.0
    Qzbounds[1:, 0, -1] = 1.0
    Qzbounds[1:, 1:, -1] = -np.eye(n_aug - 1)
    Rzbounds[:, :, -1] = 0.0
    Rzbounds[0, 0, -1] = -1.0

    # Build drift arrays
    driftbound = np.zeros((len(BoundaryLevelsLast) + 1, n_aug))
    driftbound[0, :] = np.diag(Rzbounds[:, :, 0])
    driftregimes = np.zeros((len(BoundaryLevelsLast), n_aug))
    for reg in range(len(BoundaryLevelsLast)):
        driftregimes[reg, :] = np.diag(Rz[:, :, reg])
        driftbound[reg + 1, :] = np.diag(Rzbounds[:, :, reg + 1])

    return Qz, Qzbounds, driftregimes, driftbound, BoundaryLevelsLast


def _build_first_passage_actual(server_size, map_size, mu, C, D, ga,
                                boundary_levels, quantization, b, tau,
                                use_cme, order_of_ph, pi0, theta0):
    """Build augmented generators for actual first-passage-time analysis.

    Returns Qregimes, Qbounds, driftregimes, driftbound, BoundaryLevelsLast.
    """
    (Qy, Qy0, Ry, Rydiag, ydriftregimes, Ryregimes) = \
        _build_generators(server_size, map_size, mu, C, D, ga, quantization)

    Qybounds = np.concatenate([Qy0[:, :, np.newaxis], Qy], axis=2)
    Rybounds = np.concatenate([Ry[:, :, np.newaxis], Ryregimes], axis=2)

    B = np.asarray(boundary_levels)

    # Determine PH representation
    if use_cme:
        me_system, _ = cme_parameter_calculator(order_of_ph, tau)
        S = me_system.A
        S0 = me_system.B
        alpha = me_system.C.flatten()
    else:
        alpha = np.zeros(order_of_ph)
        alpha[0] = 1.0
        ss_diag = -np.ones(order_of_ph)
        S = np.diag(ss_diag)
        for y in range(order_of_ph - 1):
            S[y, y + 1] = 1.0
        S = order_of_ph * S / tau
        e = np.ones((len(S), 1))
        S0 = -S @ e

    # Find boundary level containing b and extend
    count = 0
    while B[count] < b:
        count += 1
    BoundaryLevelsLast = np.zeros(len(B) + 1)
    BoundaryLevelsLast[:count] = B[:count]
    BoundaryLevelsLast[count] = b
    BoundaryLevelsLast[count + 1:] = B[count:]

    n_ph = len(S)
    pi0_arr = np.asarray(pi0).flatten()
    theta0_arr = np.asarray(theta0).flatten()
    kron_init = np.kron(alpha, np.kron(pi0_arr, theta0_arr))
    n_aug = 2 + n_ph * Qy.shape[0]

    # Extend Qy and Ry arrays for the extra boundary
    Rytemp = np.concatenate([Ryregimes, Ry[:, :, np.newaxis]], axis=2)
    Ryboundstemp = np.concatenate([Rybounds, Ry[:, :, np.newaxis]], axis=2)

    Qytemp = np.zeros((Qy.shape[0], Qy.shape[1], Qy.shape[2] + 1))
    Qytemp[:, :, :count] = Qy[:, :, :count]
    Qytemp[:, :, count] = Qy[:, :, min(count, Qy.shape[2] - 1)]
    Qytemp[:, :, count + 1:] = Qy[:, :, count:]

    Qyboundstemp = np.zeros((Qybounds.shape[0], Qybounds.shape[1], Qybounds.shape[2] + 1))
    Qyboundstemp[:, :, :count + 1] = Qybounds[:, :, :count + 1]
    Qyboundstemp[:, :, count + 1] = Qybounds[:, :, count + 1] if count + 1 < Qybounds.shape[2] else Qybounds[:, :, -1]
    Qyboundstemp[:, :, count + 2:] = Qybounds[:, :, count + 1:]

    num_regimes = len(BoundaryLevelsLast)
    QzTilde = np.zeros((n_aug, n_aug, num_regimes))
    QzboundsTilde = np.zeros((n_aug, n_aug, num_regimes + 1))
    RzTilde = np.zeros((n_aug, n_aug, num_regimes))
    RzboundsTilde = np.zeros((n_aug, n_aug, num_regimes + 1))

    for reg in range(num_regimes):
        eTilde = np.ones(Qybounds.shape[0])
        eTilde[-map_size:] = 0.0
        ITilde = np.diag(eTilde)

        # QzTilde regime
        Qz_block = np.kron(np.eye(n_ph), Qytemp[:, :, reg]) + np.kron(S, ITilde)
        kron_S0_ITilde = np.kron(S0.flatten(), np.diag(ITilde))
        QzTilde[0, :, reg] = 0.0
        QzTilde[1, :, reg] = 0.0
        QzTilde[2:, 0, reg] = 0.0
        QzTilde[2:, 1, reg] = kron_S0_ITilde
        QzTilde[2:, 2:, reg] = Qz_block

        # QzboundsTilde (reg+1)
        Qzb_block = np.kron(np.eye(n_ph), Qyboundstemp[:, :, reg + 1]) + np.kron(S, ITilde)
        QzboundsTilde[0, :, reg + 1] = 0.0
        QzboundsTilde[1, :, reg + 1] = 0.0
        QzboundsTilde[2:, 0, reg + 1] = 0.0
        QzboundsTilde[2:, 1, reg + 1] = kron_S0_ITilde
        QzboundsTilde[2:, 2:, reg + 1] = Qzb_block

        # Drift matrices
        Ry_reg_diag = np.diag(Rytemp[:, :, reg])
        Rz_diag = np.concatenate([[-1.0, -1.0], np.kron(np.ones(n_ph), Ry_reg_diag)])
        RzTilde[:, :, reg] = np.diag(Rz_diag)

        Ryb_diag = np.diag(Ryboundstemp[:, :, reg + 1])
        Rzb_diag = np.concatenate([[-1.0, -1.0], np.kron(np.ones(n_ph), Ryb_diag)])
        RzboundsTilde[:, :, reg + 1] = np.diag(Rzb_diag)

    # Boundary at level 0 (special structure for actual waiting time)
    eTilde = np.ones(Qybounds.shape[0])
    eTilde[-map_size:] = 0.0
    ITilde = np.diag(eTilde)

    Qzb0_block = np.kron(np.eye(n_ph), Qybounds[:, :, 0]) + np.kron(S, ITilde)
    kron_S0_ITilde = np.kron(S0.flatten(), np.diag(ITilde))
    QzboundsTilde[0, 0, 0] = -1.0
    QzboundsTilde[0, 1, 0] = 1.0
    QzboundsTilde[0, 2:, 0] = 0.0
    QzboundsTilde[1, 0, 0] = 0.0
    QzboundsTilde[1, 1, 0] = -1.0
    QzboundsTilde[1, 2:, 0] = kron_init
    QzboundsTilde[2:, 0, 0] = 0.0
    QzboundsTilde[2:, 1, 0] = kron_S0_ITilde
    QzboundsTilde[2:, 2:, 0] = Qzb0_block

    Ryb0_diag = np.diag(Ryboundstemp[:, :, 0])
    Rzb0_diag = np.concatenate([[-1.0, -1.0], np.kron(np.ones(n_ph), Ryb0_diag)])
    RzboundsTilde[:, :, 0] = np.diag(Rzb0_diag)

    # For regimes beyond b (count+1 onward): modify to redirect arrivals to absorbing state
    qy_size = Qy.shape[1]
    for reg in range(count, num_regimes):
        for d in range(order_of_ph):
            index = 2 + d * qy_size + qy_size - 2 * map_size
            # Redirect arrivals from queue (state s) to absorbing state (column 0)
            QzboundsTilde[index:index + map_size, 0, reg + 1] = \
                QzboundsTilde[index:index + map_size,
                              index + map_size:index + 2 * map_size, reg + 1] @ \
                np.ones(map_size)
            QzboundsTilde[index:index + map_size,
                          index + map_size:index + 2 * map_size, reg + 1] = 0.0

            QzTilde[index:index + map_size, 0, reg] = \
                QzTilde[index:index + map_size,
                        index + map_size:index + 2 * map_size, reg] @ \
                np.ones(map_size)
            QzTilde[index:index + map_size,
                    index + map_size:index + 2 * map_size, reg] = 0.0

    # Also modify last boundary
    if num_regimes > 0:
        reg_last = num_regimes  # last boundary index
        for d in range(order_of_ph):
            index = 2 + d * qy_size + qy_size - 2 * map_size
            QzboundsTilde[index:index + map_size, 0, reg_last] = \
                QzboundsTilde[index:index + map_size,
                              index + map_size:index + 2 * map_size, reg_last] @ \
                np.ones(map_size)
            QzboundsTilde[index:index + map_size,
                          index + map_size:index + 2 * map_size, reg_last] = 0.0

    # Build drift arrays
    driftbound = np.zeros((num_regimes + 1, n_aug))
    driftbound[0, :] = np.diag(RzboundsTilde[:, :, 0])
    driftregimes = np.zeros((num_regimes, n_aug))
    for reg in range(num_regimes):
        driftregimes[reg, :] = np.diag(RzTilde[:, :, reg])
        driftbound[reg + 1, :] = np.diag(RzboundsTilde[:, :, reg + 1])

    return QzTilde, QzboundsTilde, driftregimes, driftbound, BoundaryLevelsLast


def solve_first_passage_virtual(server_size, map_size, mu, C, D, ga,
                                boundary_levels, quantization, b, tau,
                                use_cme=False, order_of_ph=25,
                                pi0=None, theta0=None):
    """Compute first-passage-time probability for virtual waiting time.

    Parameters
    ----------
    server_size : int
        Number of servers.
    map_size : int
        Order of the MAP(C,D) arrival process.
    mu : float
        Service rate.
    C : ndarray
        MAP D0 matrix.
    D : ndarray
        MAP D1 matrix.
    ga : ndarray
        Abandonment probabilities (starts with 0).
    boundary_levels : ndarray
        Boundary levels of regimes.
    quantization : int
        Number of regimes.
    b : float
        Threshold level for first passage.
    tau : float
        Time horizon.
    use_cme : bool, optional
        If True, use CME approximation; if False, use Erlangization.
    order_of_ph : int, optional
        Order of the PH/CME approximation (25, 51, or 101).
    pi0 : array_like, optional
        Initial server occupancy distribution.
    theta0 : array_like, optional
        Initial MAP state distribution.

    Returns
    -------
    result : MAPMsGResult
        Result object with ``first_passage_virtual`` attribute.
    """
    if pi0 is None:
        pi0 = np.zeros(server_size + 1)
        pi0[0] = 1.0
    if theta0 is None:
        theta0 = np.zeros(map_size)
        theta0[0] = 1.0

    Qregimes, Qbounds, driftregimes, driftbound, BoundaryLevelsLast = \
        _build_first_passage_virtual(
            server_size, map_size, mu, C, D, ga, boundary_levels,
            quantization, b, tau, use_cme, order_of_ph, pi0, theta0)

    (coefficients, boundaries, Lzero_multi, Lneg_multi, Lpos_multi,
     Aneg_multi, Apos_multi) = \
        mrmfq_solver(Qregimes, Qbounds, driftregimes, driftbound,
                     BoundaryLevelsLast)

    ca = boundaries[0]
    cb = boundaries[-1]
    first_passage = np.sum(cb) / ca[0]

    result = MAPMsGResult()
    result.first_passage_virtual = first_passage
    result.coefficients = coefficients
    result.boundaries = boundaries
    result.Lzero_multi = Lzero_multi
    result.Lneg_multi = Lneg_multi
    result.Lpos_multi = Lpos_multi
    result.Aneg_multi = Aneg_multi
    result.Apos_multi = Apos_multi
    result.boundary_levels = BoundaryLevelsLast
    return result


def solve_first_passage_actual(server_size, map_size, mu, C, D, ga,
                               boundary_levels, quantization, b, tau,
                               use_cme=False, order_of_ph=25,
                               pi0=None, theta0=None):
    """Compute first-passage-time probability for actual waiting time.

    Parameters
    ----------
    server_size : int
        Number of servers.
    map_size : int
        Order of the MAP(C,D) arrival process.
    mu : float
        Service rate.
    C : ndarray
        MAP D0 matrix.
    D : ndarray
        MAP D1 matrix.
    ga : ndarray
        Abandonment probabilities (starts with 0).
    boundary_levels : ndarray
        Boundary levels of regimes.
    quantization : int
        Number of regimes.
    b : float
        Threshold level for first passage.
    tau : float
        Time horizon.
    use_cme : bool, optional
        If True, use CME approximation; if False, use Erlangization.
    order_of_ph : int, optional
        Order of the PH/CME approximation (25, 51, or 101).
    pi0 : array_like, optional
        Initial server occupancy distribution.
    theta0 : array_like, optional
        Initial MAP state distribution.

    Returns
    -------
    result : MAPMsGResult
        Result object with ``first_passage_actual`` attribute.
    """
    if pi0 is None:
        pi0 = np.zeros(server_size + 1)
        pi0[0] = 1.0
    if theta0 is None:
        theta0 = np.zeros(map_size)
        theta0[0] = 1.0

    Qregimes, Qbounds, driftregimes, driftbound, BoundaryLevelsLast = \
        _build_first_passage_actual(
            server_size, map_size, mu, C, D, ga, boundary_levels,
            quantization, b, tau, use_cme, order_of_ph, pi0, theta0)

    (coefficients, boundaries, Lzero_multi, Lneg_multi, Lpos_multi,
     Aneg_multi, Apos_multi) = \
        mrmfq_solver(Qregimes, Qbounds, driftregimes, driftbound,
                     BoundaryLevelsLast)

    ca = boundaries[0]
    first_passage = ca[0] / ca[1]

    result = MAPMsGResult()
    result.first_passage_actual = first_passage
    result.coefficients = coefficients
    result.boundaries = boundaries
    result.Lzero_multi = Lzero_multi
    result.Lneg_multi = Lneg_multi
    result.Lpos_multi = Lpos_multi
    result.Aneg_multi = Aneg_multi
    result.Apos_multi = Apos_multi
    result.boundary_levels = BoundaryLevelsLast
    return result
