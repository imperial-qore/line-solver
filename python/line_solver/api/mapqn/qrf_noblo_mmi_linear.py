"""
QRF No-Blocking NLP with Linear Constraint Matrices.

Builds sparse Aeq/beq and A/b matrices for linear constraints, then
solves an NLP with nonlinear objective (MMI) subject to these linear
constraints. This is computationally more efficient than callback-based
constraints for larger problems.

The equality block must be reduced to an independent row set before it
reaches SLSQP (see qrf_noblo_common.reduce_equalities). SLSQP cannot solve a
problem with more equality constraints than variables -- that is its exit
mode 2 -- and scipy 1.16.3 does not merely refuse such a problem: its Fortran
workspace is sized by a closed form that is only valid for meq <= n, so with
meq > n it under-allocates and the solver writes past the end of the buffer.
On a 2-station N=3 instance (n=66, meq=133, mieq=4) the buffer comes to 22403
doubles against the 23051 that scipy's own comment attributes to the LSQ
stage alone, and the process dies with "double free or corruption". The
mieq == 0 branch tops the buffer up by 2*n*(n+1), which is why an otherwise
identical problem without inequality rows survives.

Port of MATLAB qrf_noblo_mmi_linear.m.

The 'linear' in the name is about HOW the constraints are built -- emitted
directly as sparse matrices rather than recovered from a residual callback --
not about which they are, and not about the objective, which is the nonlinear
MMI. Until 2026-08-29 the MATLAB reference called its own mem() here (its
mmi() survived only in a commented-out line) and all three ports mirrored
that, so the token named for mutual-information minimisation computed an
entropy extremum instead.
"""

import numpy as np
from scipy.sparse import lil_matrix

from .qrf_noblo_common import (
    sub_qrfvar,
    mmi_objective,
    extract_busy,
    extract_results,
    extract_mu_v_from_maps,
    build_q_ld,
    reduce_equalities,
    feasible_start,
    qrf_index_map,
    solve_qrf_nlp,
    mmi_gradient,
)


def _compact_dims(K, N, MR):
    """Compact (variable-phase) layout dimensions.

    The decision vector produced/consumed by sub_qrfvar uses a COMPACT layout:
    the phase indices k (for station j) and h (for station i) range only over
    1..K[j] and 1..K[i] respectively, NOT over max(K). Earlier revisions of
    _deltap2/_deltae assumed a max(K)-padded layout, which agrees with the
    compact one only when all K are equal; for mixed-phase models (e.g. an
    Erlang-2 station beside an exponential one) it indexed the wrong columns,
    silently corrupting the linear constraint matrices. These helpers reproduce
    the exact sequential order of sub_qrfvar.

    Returns (cumK, SK, P) where cumK[i]=sum_{i'<i} K[i'], SK=sum K,
    P=(N+1)^2*MR*SK^2 is the number of p2 variables.
    """
    K = np.asarray(K, dtype=int)
    cumK = np.concatenate([[0], np.cumsum(K)])
    SK = int(K.sum())
    P = (N + 1) ** 2 * MR * SK * SK
    return cumK, SK, P


def _compact_num_vars(K, N, MR):
    """Total compact decision-vector length: P p2-vars + SK e-vars."""
    _, SK, P = _compact_dims(K, N, MR)
    return P + SK


def _deltap2(j, nj, k, i, ni, h, m, M, N, K, MR):
    """Compact 0-based column index for p2[j,nj,k,i,ni,h,m].

    Matches the sequential fill order of sub_qrfvar (j,nj,k,i,ni,h,m with
    k in 0..K[j]-1, h in 0..K[i]-1). All inputs 0-based.
    """
    cumK, SK, _ = _compact_dims(K, N, MR)
    Kj = int(K[j]); Ki = int(K[i])
    return (cumK[j] * (N + 1) ** 2 * MR * SK
            + nj * Kj * (N + 1) * MR * SK
            + k * (N + 1) * MR * SK
            + cumK[i] * (N + 1) * MR
            + ni * Ki * MR
            + h * MR
            + m)


def _deltae(i, k, M, N, K, MR):
    """Compact 0-based column index for e[i,k] (after all p2 vars)."""
    cumK, SK, P = _compact_dims(K, N, MR)
    return P + cumK[i] + k


def qrf_noblo_mmi_linear(MAPs, N, rt, alpha=None):
    """QRF no-blocking NLP with linear constraint matrices.

    Builds sparse equality/inequality matrices and uses SLSQP with the MMI objective.

    Args:
        MAPs: List of [D0, D1] pairs per queue
        N: Total population
        rt: Routing matrix [M, M]
        alpha: Load-dependent scaling [M, N]. If None, defaults to ones.

    Returns:
        UN: Utilization per queue [M]
        QN: Queue length per queue [M]
    """
    M = len(MAPs)
    K = np.array([MAPs[i][0].shape[0] for i in range(M)], dtype=int)
    Kmax = int(max(K))

    if alpha is None:
        alpha = np.ones((M, N))

    mu, v = extract_mu_v_from_maps(MAPs, M, K)

    MR = 1
    BB = np.zeros((MR, M))
    F = np.full(M, N, dtype=int)

    q = build_q_ld(M, K, mu, v, rt, N, alpha)

    # Compact layout (consistent with sub_qrfvar); see _compact_dims.
    num_vars = _compact_num_vars(K, N, MR)

    # Build sparse constraint matrices
    Aeq, beq, Aub, bub = _build_linear_constraints(q, M, MR, BB, F, N, K, num_vars)

    # Bounds
    bounds = [(0.0, 1.0)] * num_vars

    # see _kb/03-api-layer.md for rationale
    Aeq_dense, beq_arr, _ = reduce_equalities(Aeq.toarray(), np.array(beq))
    Aub_dense = Aub.toarray() if Aub.shape[0] > 0 else np.zeros((0, num_vars))
    bub_arr = np.array(bub) if len(bub) > 0 else np.zeros(0)

    # see _kb/03-api-layer.md for rationale
    x0 = feasible_start(Aeq_dense, beq_arr, Aub_dense, bub_arr, num_vars)

    idx = qrf_index_map(M, N, K, MR)
    xopt = solve_qrf_nlp(
        lambda x: mmi_objective(x, M, N, K, F, MR),
        lambda x: mmi_gradient(x, M, N, K, F, MR, idx),
        x0, Aeq_dense, beq_arr, Aub_dense, bub_arr, 'qrf_noblo_mmi_linear')

    p2opt, _ = sub_qrfvar(xopt, M, N, K, MR)
    UN, QN = extract_results(p2opt, M, K, F, MR)
    BN = extract_busy(p2opt, M, K, F, MR, alpha)

    return UN, QN, BN


def _build_linear_constraints(q, M, MR, BB, F, N, K, num_vars):
    """Build sparse Aeq, beq, Aub, bub matrices.

    Port of MATLAB sub_qrfcon() linear constraint version.

    q must be the 5D load-dependent array from build_q_ld: this builder mirrors
    MATLAB qrf_noblo_mmi_linear.m, whose q always carries the population index
    of the emitting station.

    Returns:
        Aeq: Sparse equality constraint matrix
        beq: Equality RHS vector
        Aub: Sparse inequality constraint matrix (A*x <= b)
        bub: Inequality RHS vector
    """
    q = np.asarray(q)
    if q.ndim != 5:
        raise ValueError(
            "q must be 5D [M,M,Kmax,Kmax,N+1] (see build_q_ld); got ndim=%d. "
            "The 4D population-free form belongs to qrf_noblo_mmi/mem."
            % q.ndim)

    # Use lists and convert to sparse at end
    aeq_rows = []
    beq_vals = []
    aub_rows = []
    bub_vals = []

    # Temporary dense row accumulator
    def new_eq_row():
        return np.zeros(num_vars)

    # ONE: normalization
    for j in range(M):
        row = new_eq_row()
        for nj in range(N + 1):
            for k in range(K[j]):
                for m in range(MR):
                    idx = _deltap2(j, nj, k, j, nj, k, m, M, N, K, MR)
                    row[idx] += 1.0
        aeq_rows.append(row)
        beq_vals.append(1.0)

    # ZERO1: i==j and ni==nj and h!=k
    for j in range(M):
        for k in range(K[j]):
            for nj in range(N + 1):
                for i in range(M):
                    for h in range(K[i]):
                        for ni in range(N + 1):
                            for m in range(MR):
                                if i == j and nj == ni and h != k:
                                    row = new_eq_row()
                                    idx = _deltap2(j, nj, k, i, ni, h, m, M, N, K, MR)
                                    row[idx] = 1.0
                                    aeq_rows.append(row)
                                    beq_vals.append(0.0)

    # ZERO2: i==j and nj!=ni
    for j in range(M):
        for k in range(K[j]):
            for nj in range(N + 1):
                for i in range(M):
                    for h in range(K[i]):
                        for ni in range(N + 1):
                            for m in range(MR):
                                if i == j and nj != ni:
                                    row = new_eq_row()
                                    idx = _deltap2(j, nj, k, i, ni, h, m, M, N, K, MR)
                                    row[idx] = 1.0
                                    aeq_rows.append(row)
                                    beq_vals.append(0.0)

    # ZERO3: i!=j and nj+ni>N
    for j in range(M):
        for k in range(K[j]):
            for nj in range(N + 1):
                for i in range(M):
                    for h in range(K[i]):
                        for ni in range(N + 1):
                            for m in range(MR):
                                if i != j and nj + ni > N:
                                    row = new_eq_row()
                                    idx = _deltap2(j, nj, k, i, ni, h, m, M, N, K, MR)
                                    row[idx] = 1.0
                                    aeq_rows.append(row)
                                    beq_vals.append(0.0)

    # ZERO5: BB[m,j]==1 for m>=1 (0-based)
    for j in range(M):
        for k in range(K[j]):
            for i in range(M):
                for h in range(K[i]):
                    for ni in range(F[i] + 1):
                        for m in range(1, MR):
                            if BB[m, j] == 1:
                                row = new_eq_row()
                                idx = _deltap2(j, 0, k, i, ni, h, m, M, N, K, MR)
                                row[idx] = 1.0
                                aeq_rows.append(row)
                                beq_vals.append(0.0)

    # ZERO6: nj > F[j]
    for j in range(M):
        for k in range(K[j]):
            for nj in range(F[j] + 1, N + 1):
                for i in range(M):
                    for h in range(K[i]):
                        for ni in range(N + 1):
                            for m in range(MR):
                                row = new_eq_row()
                                idx = _deltap2(j, nj, k, i, ni, h, m, M, N, K, MR)
                                row[idx] = 1.0
                                aeq_rows.append(row)
                                beq_vals.append(0.0)

    # SYMMETRY
    for j in range(M):
        for nj in range(N + 1):
            for k in range(K[j]):
                for i in range(M):
                    for ni in range(N + 1):
                        for h in range(K[i]):
                            for m in range(MR):
                                row = new_eq_row()
                                idx1 = _deltap2(i, ni, h, j, nj, k, m, M, N, K, MR)
                                idx2 = _deltap2(j, nj, k, i, ni, h, m, M, N, K, MR)
                                row[idx1] += 1.0
                                row[idx2] -= 1.0
                                aeq_rows.append(row)
                                beq_vals.append(0.0)

    # MARGINALS
    for j in range(M):
        for k in range(K[j]):
            for nj in range(N + 1):
                for i in range(M):
                    for m in range(MR):
                        if i != j:
                            row = new_eq_row()
                            # LHS: p2[j,nj,k,j,nj,k,m]
                            idx = _deltap2(j, nj, k, j, nj, k, m, M, N, K, MR)
                            row[idx] += 1.0
                            # see _kb/03-api-layer.md for rationale
                            for ni in range(N + 1):
                                for h in range(K[i]):
                                    idx = _deltap2(j, nj, k, i, ni, h, m, M, N, K, MR)
                                    row[idx] -= 1.0
                            aeq_rows.append(row)
                            beq_vals.append(0.0)

    # UEFF
    for j in range(M):
        for i in range(M):
            for ki in range(K[i]):
                row = new_eq_row()
                # LHS: e[i,ki]
                idx = _deltae(i, ki, M, N, K, MR)
                row[idx] += 1.0
                # RHS: -sum p2[j,nj,kj,i,ni,ki,m]
                for nj in range(N + 1):
                    for kj in range(K[j]):
                        for m in range(MR):
                            for ni in range(1, N + 1):
                                if BB[m, i] == 0:
                                    idx = _deltap2(j, nj, kj, i, ni, ki, m, M, N, K, MR)
                                    row[idx] -= 1.0
                aeq_rows.append(row)
                beq_vals.append(0.0)

    # see _kb/03-api-layer.md for rationale
    for i in range(M):
        for k in range(K[i]):
            row = new_eq_row()
            for ni in range(1, F[i] + 1):
                for m in range(MR):
                    # LHS
                    for j in range(M):
                        for h in range(K[i]):
                            idx = _deltap2(i, ni, k, i, ni, k, m, M, N, K, MR)
                            row[idx] += q[i, j, k, h, ni]
                    # RHS
                    for j in range(M):
                        for h in range(K[i]):
                            idx = _deltap2(i, ni, h, i, ni, h, m, M, N, K, MR)
                            row[idx] -= q[i, j, h, k, ni]
            aeq_rows.append(row)
            beq_vals.append(0.0)

    # THM2
    for j in range(M):
        for k in range(K[j]):
            for nj in range(F[j] + 1):
                for m in range(MR):
                    row = new_eq_row()
                    # LHS
                    for i in range(M):
                        for ni in range(1, F[i] + 1):
                            for ki in range(K[i]):
                                idx = _deltap2(j, nj, k, i, ni, ki, m, M, N, K, MR)
                                row[idx] += ni
                    # RHS: -N*p2[j,nj,k,j,nj,k,m]
                    idx = _deltap2(j, nj, k, j, nj, k, m, M, N, K, MR)
                    row[idx] -= N
                    aeq_rows.append(row)
                    beq_vals.append(0.0)

    # COR1
    row = new_eq_row()
    for m in range(MR):
        for i in range(M):
            for j in range(M):
                for nj in range(1, F[j] + 1):
                    for ni in range(1, F[i] + 1):
                        for ki in range(K[i]):
                            for kj in range(K[j]):
                                idx = _deltap2(j, nj, kj, i, ni, ki, m, M, N, K, MR)
                                row[idx] += ni * nj
    aeq_rows.append(row)
    beq_vals.append(N ** 2)

    # see _kb/03-api-layer.md for rationale

    # see _kb/03-api-layer.md for rationale
    for i in range(M):
        for u in range(K[i]):
            row = new_eq_row()
            for j in range(M):
                if j == i:
                    continue
                for nj in range(1, F[j] + 1):
                    for k in range(K[j]):
                        coef = 0.0
                        for h in range(K[j]):
                            coef += q[j, i, k, h, nj]
                        if coef == 0.0:
                            continue
                        for m in range(MR):
                            idx = _deltap2(j, nj, k, i, 0, u, m, M, N, K, MR)
                            row[idx] += coef
            for j in range(M):
                if j == i:
                    continue
                for nj in range(F[j] + 1):
                    for k in range(K[i]):
                        # see _kb/03-api-layer.md for rationale
                        coef = q[i, j, k, u, 1]
                        if coef == 0.0:
                            continue
                        for h in range(K[j]):
                            for m in range(MR):
                                idx = _deltap2(j, nj, h, i, 1, k, m, M, N, K, MR)
                                row[idx] -= coef
            aeq_rows.append(row)
            beq_vals.append(0.0)

    # THM3 {i, ni in 0..F[i]-1}: balance across the ni -> ni+1 boundary.
    for i in range(M):
        for ni in range(F[i]):
            row = new_eq_row()
            for j in range(M):
                if j == i:
                    continue
                for nj in range(1, F[j] + 1):
                    for k in range(K[j]):
                        coef = 0.0
                        for h in range(K[j]):
                            coef += q[j, i, k, h, nj]
                        if coef == 0.0:
                            continue
                        for u in range(K[i]):
                            for m in range(MR):
                                idx = _deltap2(j, nj, k, i, ni, u, m, M, N, K, MR)
                                row[idx] += coef
            for j in range(M):
                if j == i:
                    continue
                for nj in range(F[j] + 1):
                    for k in range(K[i]):
                        coef = 0.0
                        for h in range(K[i]):
                            coef += q[i, j, k, h, ni + 1]
                        if coef == 0.0:
                            continue
                        for u in range(K[j]):
                            for m in range(MR):
                                idx = _deltap2(j, nj, u, i, ni + 1, k, m,
                                               M, N, K, MR)
                                row[idx] -= coef
            aeq_rows.append(row)
            beq_vals.append(0.0)

    # THM4 (inequality)
    for j in range(M):
        for k in range(K[j]):
            for i in range(M):
                for m in range(MR):
                    row = new_eq_row()
                    # -LHS (sign swapped: >= becomes <=)
                    for t in range(M):
                        for h in range(K[t]):
                            for nj_idx in range(N + 1):
                                for nt in range(N + 1):
                                    idx = _deltap2(j, nj_idx, k, t, nt, h, m, M, N, K, MR)
                                    row[idx] -= nt
                    # +RHS
                    for h in range(K[i]):
                        for nj_idx in range(N + 1):
                            for ni in range(1, N + 1):
                                idx = _deltap2(j, nj_idx, k, i, ni, h, m, M, N, K, MR)
                                row[idx] += N
                    aub_rows.append(row)
                    bub_vals.append(0.0)

    # Convert to sparse matrices
    if aeq_rows:
        Aeq = lil_matrix((len(aeq_rows), num_vars))
        for r, row in enumerate(aeq_rows):
            nz = np.nonzero(row)[0]
            for c in nz:
                Aeq[r, c] = row[c]
        Aeq = Aeq.tocsr()
    else:
        Aeq = lil_matrix((0, num_vars)).tocsr()

    if aub_rows:
        Aub = lil_matrix((len(aub_rows), num_vars))
        for r, row in enumerate(aub_rows):
            nz = np.nonzero(row)[0]
            for c in nz:
                Aub[r, c] = row[c]
        Aub = Aub.tocsr()
    else:
        Aub = lil_matrix((0, num_vars)).tocsr()

    return Aeq, beq_vals, Aub, bub_vals
