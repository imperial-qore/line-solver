"""
Native Python port of the FJ_codes response-time percentile algorithm
("Beyond the Mean in Fork-Join Queues", Qiu, Perez, Harrison, IFIP Performance
2015), i.e. the MATLAB matlab/lib/thirdparty/FJ_codes/returnRT2.m pipeline.

This reproduces the MATLAB reference exactly (K=2 fork-join, exponential service),
so the native SolverMAM percentile output matches MATLAB and the JAR. Uses scipy
for the ordered real Schur decomposition (computeT_NARE) and the Sylvester solve
(computePi's lyap).
"""

import numpy as np
from scipy.linalg import schur, solve_sylvester


def _build_index(m, cr):
    """Port of build_index.m: rows are the compositions of cr into m parts."""
    from math import comb
    total = comb(cr + m - 1, cr)
    idx = np.zeros((total, m), dtype=int)
    idx[0, 0] = cr
    for row in range(1, total):
        nz = np.nonzero(idx[row - 1, :] > 0)[0]
        k = nz[0] if nz.size else 0
        if k < m - 1:
            idx[row, :] = idx[row - 1, :]
            idx[row, k + 1] += 1
            idx[row, 0] = idx[row, k] - 1
            idx[row, 1:k + 1] = 0
    return idx


def _vectmatch(row, matrix):
    """Port of vectmatch.m: 0-based index of `row` in `matrix` rows."""
    for i in range(matrix.shape[0]):
        if np.array_equal(matrix[i, :], row):
            return i
    return -1


def _kronsum(A, B):
    return np.kron(A, np.eye(B.shape[1])) + np.kron(np.eye(A.shape[1]), B)


def _build_service_h(tau_st, ST):
    """Port of build_Service_h.m for the 2-node FJ job phase."""
    dim_single = len(tau_st)
    phases_single = _build_index(dim_single, 1)
    n = dim_single * dim_single
    service_phases = np.zeros((n, 2 * phases_single.shape[1]), dtype=int)
    k = 0
    for i in range(dim_single):
        for j in range(dim_single):
            service_phases[k, :] = np.concatenate([phases_single[i, :], phases_single[j, :]])
            k += 1
    beta = np.kron(tau_st, tau_st)
    Smat = _kronsum(ST, ST)
    return {'service_phases': service_phases, 'beta': beta, 'S': Smat}


def _build_sa(tau_st, ST, service_h, C):
    """Port of build_SA.m -> (S, A_jump)."""
    sp = service_h['service_phases']
    dim = len(service_h['beta'])
    m = len(tau_st)
    dim_C = C + 1
    newdim = dim_C * dim
    S = np.zeros((newdim, newdim))
    A_jump = np.zeros((newdim, newdim))
    Sh = service_h['S']
    for row in range(dim_C):
        S[row * dim:(row + 1) * dim, row * dim:(row + 1) * dim] = Sh
    STm = np.atleast_2d(ST)
    A = (-np.sum(STm, axis=1).reshape(-1, 1)) @ np.atleast_2d(tau_st)  # (m x m)

    S_Cminus1 = np.zeros((dim, dim))
    for row in range(dim):
        cv = sp[row, :]
        for i in range(m):
            if cv[i] > 0:
                for j in range(m):
                    tov = cv.copy()
                    tov[i] -= 1
                    tov[j] += 1
                    col = _vectmatch(tov, sp)
                    S_Cminus1[row, col] += cv[i] * A[i, j]
    for c in range(C, 0, -1):
        S[(C - c) * dim:(C - c + 1) * dim, (C - c + 1) * dim:(C - c + 2) * dim] = S_Cminus1

    A_Cplus1 = np.zeros((dim, dim))
    for row in range(dim):
        cv = sp[row, :]
        for i in range(m, 2 * m):
            if cv[i] > 0:
                for j in range(m, 2 * m):
                    tov = cv.copy()
                    tov[i] -= 1
                    tov[j] += 1
                    col = _vectmatch(tov, sp)
                    A_Cplus1[row, col] += A[i - m, j - m]
    for c in range(C - 1, 0, -1):
        A_jump[(C - c) * dim:(C - c + 1) * dim, (C - c - 1) * dim:(C - c) * dim] = A_Cplus1
    A_jump[0:dim, 0:dim] = A_Cplus1

    A_last = np.zeros((dim, dim))
    for row in range(dim):
        cv = sp[row, :]
        for k in range(2):
            for i in range(m):
                if cv[k * m + i] > 0:
                    tov = cv[(1 - k) * m:(2 - k) * m]
                    for j in range(m):
                        tmp = np.zeros(m, dtype=int)
                        tmp[j] = 1
                        target = np.concatenate([tov, tmp])
                        col = _vectmatch(target, sp)
                        A_last[row, col] += cv[k * m + i] * A[i, j]
    A_jump[C * dim:, (C - 1) * dim:C * dim] = A_last
    return S, A_jump


def _construct_not_all_busy(C, tau_st, ST, St, service_h):
    """Port of constructNotAllBusy.m."""
    m = len(tau_st)
    indexes_nb = _build_index(m, 1)
    dim_NB = indexes_nb.shape[0]
    dim_C = C + 1
    dim_nb = (dim_C - 1) * dim_NB + 1
    S = np.zeros((dim_nb, dim_nb))
    STm = np.atleast_2d(ST)
    for row in range(dim_C - 1):
        S[row * dim_NB:(row + 1) * dim_NB, row * dim_NB:(row + 1) * dim_NB] = STm
    A = (-np.sum(STm, axis=1).reshape(-1, 1)) @ np.atleast_2d(tau_st)
    for row in range(dim_C - 2):
        S[row * dim_NB:(row + 1) * dim_NB, (row + 1) * dim_NB:(row + 2) * dim_NB] = A
    exitcol = -np.sum(STm, axis=1).reshape(-1, 1)
    S[(dim_C - 2) * dim_NB:(dim_C - 1) * dim_NB, (dim_C - 1) * dim_NB:] = exitcol
    for row in range(dim_nb):
        S[row, row] = 0.0
        S[row, row] = -np.sum(S[row, :])
    return S


def _generate_service(tau_st, ST, St, service_h, C, S):
    """Port of generateService.m -> (T, newdim, dim_notbusy)."""
    sp = service_h['service_phases']
    dim = len(service_h['beta'])
    m = len(tau_st)
    indexes_nb = _build_index(m, 1)
    dim_NB = indexes_nb.shape[0]
    dim_C = C + 1
    newdim = dim_C * dim
    dim_notbusy = dim_C * dim_NB
    N = newdim + dim_notbusy
    T = np.zeros((N, N))
    t = np.zeros((N, 1))
    Stv = np.atleast_1d(St).reshape(-1)
    t[N - dim_NB:, 0] = Stv[:dim_NB]
    T[0:newdim, 0:newdim] = S

    S_long = np.zeros((dim, dim_NB))
    for row in range(dim):
        cv = sp[row, :]
        for i in range(m, 2 * m):
            if cv[i] > 0:
                tov = cv[0:m]
                col = _vectmatch(tov, indexes_nb)
                S_long[row, col] += cv[i] * Stv[i - m]
    for row in range(dim_C - 1):
        T[row * dim:(row + 1) * dim, newdim + row * dim_NB:newdim + (row + 1) * dim_NB] = S_long

    S_last = np.zeros((dim, dim_NB))
    for row in range(dim):
        cv = sp[row, :]
        for k in range(2):
            for i in range(k * m, (k + 1) * m):
                if cv[i] > 0:
                    tov = cv[(1 - k) * m:(2 - k) * m]
                    col = _vectmatch(tov, indexes_nb)
                    S_last[row, col] += cv[i] * Stv[i - k * m]
    T[(dim_C - 1) * dim:dim_C * dim, newdim + (dim_C - 1) * dim_NB:newdim + dim_C * dim_NB] = S_last

    STm = np.atleast_2d(ST)
    for row in range(dim_C):
        T[newdim + row * dim_NB:newdim + (row + 1) * dim_NB,
          newdim + row * dim_NB:newdim + (row + 1) * dim_NB] = STm
    A = (-np.sum(STm, axis=1).reshape(-1, 1)) @ np.atleast_2d(tau_st)
    for row in range(dim_C - 1):
        T[newdim + row * dim_NB:newdim + (row + 1) * dim_NB,
          newdim + (row + 1) * dim_NB:newdim + (row + 2) * dim_NB] = A
    for row in range(newdim, newdim + dim_notbusy):
        T[row, row] = 0.0
        T[row, row] = -np.sum(T[row, :]) - t[row, 0]
    return T, newdim, dim_notbusy


def _compute_t_nare(D0, D1, S_Arr, A_jump):
    """Port of computeT_NARE.m via scipy ordered Schur (stable eigenvalues first)."""
    m = S_Arr.shape[0]
    ma = D0.shape[0]
    ms = m // ma
    A = np.kron(np.eye(ms), D0)
    B = np.kron(np.eye(ms), D1)
    Cm = np.kron(A_jump, np.eye(ma))
    H = np.block([[A, B], [-Cm, -S_Arr]])
    # Sort eigenvalues with negative real part (stable) to the leading block.
    Tsch, Z, sdim = schur(H, output='real', sort='lhp')
    Q11 = Z[0:m, 0:m]
    Q21 = Z[m:2 * m, 0:m]
    X = np.linalg.solve(Q11.T, Q21.T).T  # X = Q21 * inv(Q11)
    return S_Arr + X @ np.kron(np.eye(ms), D1)


def _compute_t(lambda0, lambda1, tau_st, ST, St, service_h, C):
    """Port of computeT.m (NARE mode) -> (T, S, A_jump, S_Arr, sum_Ajump)."""
    S, A_jump = _build_sa(tau_st, ST, service_h, C)
    d0 = lambda0.shape[0]
    S_Arr = np.kron(S, np.eye(d0))
    A_jump_Arr = np.kron(A_jump, np.eye(d0))
    T = _compute_t_nare(lambda0, lambda1, S_Arr, A_jump)
    sum_Ajump = np.sum(A_jump_Arr, axis=1).reshape(-1, 1)
    return T, S, A_jump, S_Arr, sum_Ajump


def _compute_pi(T, lambda0, lambda1, tau_st, ST, St, service_h, C, S, A_jump):
    """Port of computePi.m (exponential-service branch, SerChoice==1)."""
    ms = S.shape[1]
    S_nab = _construct_not_all_busy(C, tau_st, ST, St, service_h)
    Q0 = _kronsum(S_nab, lambda0)
    da = lambda0.shape[0]
    n = ms * da
    # Igral = lyap(T, kron(I,lambda0), -eye(n)) solves T*X + X*B + (-eye) = 0
    # -> T*X + X*B = eye  => solve_sylvester(T, B, eye)
    B = np.kron(np.eye(ms), lambda0)
    Igral = solve_sylvester(T, B, np.eye(n))
    pi0mat = Igral @ np.kron(A_jump, np.eye(da)) @ np.linalg.inv(Q0) @ np.kron(np.eye(ms), lambda1)
    # pi0 = [0..0 -1] / [pi0mat-eye, sum(inv(T),2)]   (right division, 1 x n)
    invT_rowsum = np.sum(np.linalg.inv(T), axis=1).reshape(-1, 1)
    M = np.hstack([pi0mat - np.eye(n), invT_rowsum])  # (n x n+1)
    b = np.zeros(n + 1)
    b[-1] = -1.0
    pi0 = np.linalg.lstsq(M.T, b, rcond=None)[0].reshape(1, -1)  # x @ M = b
    En1 = (1.0 / np.sum(pi0)) * np.sum(pi0 @ pi0mat)
    return pi0, En1


def _return_wait(En1, pi0, T, phi, sum_Ajump):
    """Port of returnWait.m -> (wait_alpha, wait_Smat, prob_wait, alfa)."""
    alfa = -np.linalg.solve(T.T, pi0.reshape(-1, 1)).reshape(1, -1)  # -pi0 / T
    ds = phi.shape[0]
    phi = phi.reshape(-1, 1)
    rhos = (np.diag((phi @ alfa)).reshape(1, -1)) / float(alfa @ phi)
    En0 = float(alfa @ sum_Ajump) / float(np.sum(pi0))
    prob_wait = (En0 - 1.0) / (En0 - 1.0 + En1)
    wait_alpha = prob_wait * rhos
    wait_Smat = np.zeros((ds, ds))
    a = alfa.reshape(-1)
    for i in range(ds):
        for j in range(ds):
            wait_Smat[i, j] = a[j] * T[j, i] / a[i]
    return wait_alpha, wait_Smat, prob_wait, alfa


def _return_per(vector, Mat, pers):
    """Port of returnPer.m with bisection percentile search (fast, exact)."""
    vector = vector.reshape(1, -1)
    meanRT = -float(np.sum(np.linalg.solve(Mat.T, vector.T)))  # sum(-vector/Mat)
    c = float(np.max(-np.diag(Mat)))
    m = Mat.shape[1]
    P_res = Mat / c + np.eye(m)
    M = float(np.sum(np.linalg.solve((np.eye(m) - P_res).T, vector.T)))
    a0 = float(np.sum(vector))
    sum_a = a0
    ak = []
    vP = np.sum(P_res, axis=1).reshape(-1, 1)
    while abs(sum_a - M) >= 1e-10:
        val = float(vector @ vP)
        ak.append(val)
        sum_a += val
        vP = P_res @ vP
        if len(ak) > 200000:
            break
    ak = np.asarray(ak)
    K1 = len(ak)

    def cdf_at(t):
        pM = np.exp(-c * t)
        F = pM * a0
        for ki in range(K1):
            pM = c * t * pM / (ki + 1)
            F += pM * ak[ki]
        return 1.0 - F

    out = np.zeros((len(pers), 2))
    for p in range(len(pers)):
        out[p, 0] = pers[p]
        if pers[p] < 1.0 - float(np.sum(vector)):
            out[p, 1] = 0.0
            continue
        hi = max(3.0 * meanRT, 1e-9)
        guard = 0
        while cdf_at(hi) < pers[p] and guard < 200:
            hi *= 2.0
            guard += 1
        lo = 0.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if cdf_at(mid) < pers[p]:
                lo = mid
            else:
                hi = mid
        out[p, 1] = 0.5 * (lo + hi)
    return out


def return_rt2(lambda0, lambda1, tau_st, ST, St, pers, C):
    """
    Response-time percentiles for a K=2 fork-join queue (exponential service),
    matching MATLAB returnRT2.m. `pers` on the 0-1 probability scale.
    Returns an (len(pers) x 2) array: column 0 = pers, column 1 = RT percentiles.
    """
    lambda0 = np.atleast_2d(lambda0).astype(float)
    lambda1 = np.atleast_2d(lambda1).astype(float)
    tau_st = np.atleast_1d(tau_st).astype(float)
    ST = np.atleast_2d(ST).astype(float)
    St = np.atleast_1d(St).astype(float)

    service_h = _build_service_h(tau_st, ST)
    T, S, A_jump, S_Arr, sum_Ajump = _compute_t(lambda0, lambda1, tau_st, ST, St, service_h, C)
    mWait = A_jump.shape[0]
    phi = np.sum(T - S_Arr, axis=1).reshape(-1, 1)
    pi0, En1 = _compute_pi(T, lambda0, lambda1, tau_st, ST, St, service_h, C, S, A_jump)
    wait_alpha, wait_Smat, prob_wait, alfa = _return_wait(En1, pi0, T, phi, sum_Ajump)

    STg, dim, dim_notbusy = _generate_service(tau_st, ST, St, service_h, C, S)
    ma = lambda0.shape[0]
    STg = np.kron(STg, np.eye(ma))
    pi0 = pi0 / np.sum(pi0)
    dim_ma = ma * dim
    dim_service = dim_ma + ma * dim_notbusy

    notbusy_start = np.zeros((1, dim_service))
    notbusy_start[0, 0:dim_ma] = (1.0 - prob_wait) * pi0.reshape(-1)[0:dim_ma]

    TS = T - S_Arr
    alfaTS = alfa @ TS
    sumAlfaTS = float(np.sum(alfaTS))
    busy_start = np.zeros((1, dim_service))
    busy_start[0, 0:dim_ma] = prob_wait * alfaTS.reshape(-1)[0:dim_ma] / sumAlfaTS

    Tr = STg.shape[0]
    Sc = wait_Smat.shape[1]
    TS_rowsum = np.sum(TS, axis=1)
    TS_normalized = np.zeros_like(TS)
    for i in range(TS.shape[0]):
        if TS_rowsum[i] > 1e-12:
            TS_normalized[i, :] = TS[i, :] / TS_rowsum[i]

    negST = -STg
    stat = np.linalg.solve(negST.T, busy_start.T).T  # busy_start / (-ST)
    stat = stat.reshape(-1)
    statMax = float(np.max(stat)) if stat.size else 0.0
    statTol = 1e-10 * statMax
    busy_nz = stat > statTol
    m_tr_ST = int(np.sum(busy_nz))

    tr_start_state_full = -np.sum(STg, axis=1) * stat  # -sum(ST,2)' .* stat
    tr_ST_full = np.zeros((Tr, Tr))
    for i in range(Tr):
        if busy_nz[i]:
            for j in range(Tr):
                if busy_nz[j]:
                    tr_ST_full[i, j] = STg[j, i] * stat[j] / stat[i]
    tr_ST_exit_full = -np.sum(tr_ST_full, axis=1)
    tr_ST_exit_mat = np.tile(tr_ST_exit_full[0:Sc].reshape(-1, 1), (1, Sc))

    nColsTS = TS_normalized.shape[1]
    nRowsTS = TS_normalized.shape[0]
    TS2 = np.zeros((nColsTS, nRowsTS))
    a = alfa.reshape(-1)
    for aa in range(nColsTS):
        for bb in range(nRowsTS):
            TS2[aa, bb] = TS_normalized[bb, aa] * a[bb]
    tr_TS_norm = np.zeros((nColsTS, nRowsTS))
    for aa in range(nColsTS):
        rs = np.sum(TS2[aa, :])
        if abs(rs) < 1e-12:
            rs = 1.0
        tr_TS_norm[aa, :] = TS2[aa, :] / rs
    tildeP = np.zeros((Tr, Sc))
    for i in range(Sc):
        for j in range(Sc):
            tildeP[i, j] = tr_ST_exit_mat[i, j] * tr_TS_norm[i, j]

    nz = np.nonzero(busy_nz)[0]
    tr_ST = tr_ST_full[np.ix_(nz, nz)]
    tr_start_state = tr_start_state_full[nz]
    tildeP_red = tildeP[nz, :]

    gamma_res = np.zeros((1, dim_service + m_tr_ST + Sc))
    gamma_res[0, 0:dim_service] = notbusy_start.reshape(-1)
    gamma_res[0, dim_service:dim_service + m_tr_ST] = tr_start_state

    C_res = np.zeros((Tr + m_tr_ST + Sc, Tr + m_tr_ST + Sc))
    C_res[0:Tr, 0:Tr] = STg
    C_res[Tr:Tr + m_tr_ST, Tr:Tr + m_tr_ST] = tr_ST
    C_res[Tr:Tr + m_tr_ST, Tr + m_tr_ST:Tr + m_tr_ST + Sc] = tildeP_red
    C_res[Tr + m_tr_ST:, Tr + m_tr_ST:] = wait_Smat

    return _return_per(gamma_res, C_res, np.asarray(pers, dtype=float))
