"""
Shared utilities for QRF NLP no-blocking approximation methods.

Provides variable unflattening (sub_qrfvar), objective functions (MMI, MEM),
and constraint construction (sub_qrfcon) used by all qrf_noblo_* methods.

Port of MATLAB qrf_noblo_mmi.m / qrf_noblo_mem.m shared infrastructure.
"""

import numpy as np
from scipy.linalg import qr as scipy_qr
from scipy.optimize import linprog

LOGTOL = 1e-6


def sub_qrfvar(x, M, N, K, MR):
    """Unflatten decision vector x into 7D tensor p2 and effective utilizations e.

    Args:
        x: Flat decision vector of length M*(N+1)*Kmax*M*(N+1)*Kmax*MR + M*Kmax
        M: Number of queues
        N: Total population
        K: Array of phases per queue [M]
        MR: Number of blocking configurations

    Returns:
        p2: 7D array [M, N+1, Kmax, M, N+1, Kmax, MR] joint probabilities
        e: 2D array [M, Kmax] effective utilizations
    """
    Kmax = int(max(K))
    p2_size = M * (N + 1) * Kmax * M * (N + 1) * Kmax * MR

    # Reshape p2 from flat vector - MATLAB iteration order:
    # j, nj, k, i, ni, h, m (outermost to innermost)
    p2 = np.zeros((M, N + 1, Kmax, M, N + 1, Kmax, MR))
    ctr = 0
    for j in range(M):
        for nj in range(N + 1):
            for k in range(K[j]):
                for i in range(M):
                    for ni in range(N + 1):
                        for h in range(K[i]):
                            for m in range(MR):
                                p2[j, nj, k, i, ni, h, m] = x[ctr]
                                ctr += 1

    # Extract e variables
    e = np.zeros((M, Kmax))
    for i in range(M):
        for k in range(K[i]):
            e[i, k] = x[ctr]
            ctr += 1

    return p2, e


def mmi_objective(x, M, N, K, F, MR):
    """Mutual Information Minimization objective.

    Minimizes: sum p2[i,ni,ki,j,nj,kj,m] *
        (log(p2[i,ni,ki,j,nj,kj,m]) - log(p2[i,ni,ki,i,ni,ki,m]) - log(p2[j,nj,kj,j,nj,kj,m]))

    Port of MATLAB mmi() nested function.
    """
    p2, _ = sub_qrfvar(x, M, N, K, MR)
    fobj = 0.0
    for m in range(MR):
        for i in range(M):
            for ki in range(K[i]):
                for j in range(M):
                    if i != j:
                        for kj in range(K[j]):
                            for ni in range(1, F[i] + 1):
                                for nj in range(1, F[j] + 1):
                                    pij = p2[i, ni, ki, j, nj, kj, m]
                                    pii = p2[i, ni, ki, i, ni, ki, m]
                                    pjj = p2[j, nj, kj, j, nj, kj, m]
                                    fobj += pij * (
                                        np.log(LOGTOL + pij)
                                        - np.log(LOGTOL + pii)
                                        - np.log(LOGTOL + pjj)
                                    )
    return fobj


def mem_objective(x, M, N, K, F, MR):
    """Maximum Entropy objective (minimize negative entropy).

    Minimizes: -sum p2[i,ni,k,i,ni,k,m] * log(p2[i,ni,k,i,ni,k,m])

    Port of MATLAB mem() nested function.
    """
    p2, _ = sub_qrfvar(x, M, N, K, MR)
    fobj = 0.0
    for m in range(MR):
        for i in range(M):
            for k in range(K[i]):
                for ni in range(1, F[i] + 1):
                    pval = p2[i, ni, k, i, ni, k, m]
                    fobj -= pval * np.log(LOGTOL + pval)
    return fobj


def sub_qrfcon_noblo(x, q, M, MR, BB, F, N, K):
    """Build equality and inequality constraints for QRF no-blocking.

    Port of MATLAB sub_qrfcon() for qrf_noblo_mmi.m / qrf_noblo_mem.m and, in
    the load-dependent form, qrf_noblo_mmi_ld.m.

    The ARITY of q selects the model, exactly as the AMPL skeletons do:

    * 4D q [M,M,Kmax,Kmax] is the population-free form of
      qrboundsbas_skel.mod, used by qrf_noblo_mmi.m / qrf_noblo_mem.m. THM1 is
      stated on the aggregated e[i,k], which is exact here because there is no
      alpha to make the rate population-dependent.
    * 5D q [M,M,Kmax,Kmax,N+1] is the load-dependent form of
      qrboundsrsrd_skel.mod:11, whose fifth index is the population of the
      EMITTING station (see build_q_ld). THM1 is then stated per population
      against the station-i marginal.

    THM30 and THM3 are emitted in BOTH forms; only the q lookup differs. They
    are the marginal-balance families and are what makes the polytope depend on
    the service rates at all: without them an LP over the remaining constraints
    returns the vacuous [0,1] for every instance (verified against glpsol).

    Args:
        x: Decision vector
        q: Transition rates, 4D [M, M, Kmax, Kmax] or 5D [..., N+1]
        M, MR, N: Problem dimensions
        BB: Blocking state matrix [MR, M]
        F: Capacity per queue [M]
        K: Phases per queue [M]

    Returns:
        c_ineq: Inequality constraints (c <= 0, MATLAB convention)
        ceq: Equality constraints (ceq = 0)
    """
    q = np.asarray(q)
    load_dependent = (q.ndim == 5)
    p2, e = sub_qrfvar(x, M, N, K, MR)
    ceq = []
    c = []

    # ONE: normalization
    # sum_{nj,k,m} p2[j,nj,k,j,nj,k,m] = 1  for each j
    for j in range(M):
        val = 0.0
        for nj in range(N + 1):
            for k in range(K[j]):
                for m in range(MR):
                    val += p2[j, nj, k, j, nj, k, m]
        ceq.append(val - 1.0)

    # ZERO1: i==j and ni==nj and h!=k => p2=0
    for j in range(M):
        for k in range(K[j]):
            for nj in range(N + 1):
                for i in range(M):
                    for h in range(K[i]):
                        for ni in range(N + 1):
                            for m in range(MR):
                                if i == j and nj == ni and h != k:
                                    ceq.append(p2[j, nj, k, i, ni, h, m])

    # ZERO2: i==j and nj!=ni => p2=0
    for j in range(M):
        for k in range(K[j]):
            for nj in range(N + 1):
                for i in range(M):
                    for h in range(K[i]):
                        for ni in range(N + 1):
                            for m in range(MR):
                                if i == j and nj != ni:
                                    ceq.append(p2[j, nj, k, i, ni, h, m])

    # ZERO3: i!=j and nj+ni>N => p2=0
    for j in range(M):
        for k in range(K[j]):
            for nj in range(N + 1):
                for i in range(M):
                    for h in range(K[i]):
                        for ni in range(N + 1):
                            for m in range(MR):
                                if i != j and nj + ni > N:
                                    ceq.append(p2[j, nj, k, i, ni, h, m])

    # ZERO5: BB[m,j]==1 for m>=1 (0-based) => p2[j,0,...]=0
    for j in range(M):
        for k in range(K[j]):
            for i in range(M):
                for h in range(K[i]):
                    for ni in range(F[i] + 1):
                        for m in range(1, MR):  # m from 1 (2nd config, 0-based)
                            if BB[m, j] == 1:
                                ceq.append(p2[j, 0, k, i, ni, h, m])

    # ZERO6: nj > F[j] => p2=0
    for j in range(M):
        for k in range(K[j]):
            for nj in range(F[j] + 1, N + 1):
                for i in range(M):
                    for h in range(K[i]):
                        for ni in range(N + 1):
                            for m in range(MR):
                                ceq.append(p2[j, nj, k, i, ni, h, m])

    # ZERO7: blocking constraint (only active when MR>1 and BB has nonzeros)
    # Skipped for no-blocking (MR=1) since m range is 1:MR-1 which is empty

    # SYMMETRY: p2[i,ni,h,j,nj,k,m] = p2[j,nj,k,i,ni,h,m]
    for j in range(M):
        for nj in range(N + 1):
            for k in range(K[j]):
                for i in range(M):
                    for ni in range(N + 1):
                        for h in range(K[i]):
                            for m in range(MR):
                                ceq.append(
                                    p2[i, ni, h, j, nj, k, m]
                                    - p2[j, nj, k, i, ni, h, m]
                                )

    # MARGINALS: p2[j,nj,k,j,nj,k,m] = sum_{ni,h} p2[j,nj,k,i,ni,h,m] for i!=j
    for j in range(M):
        for k in range(K[j]):
            for nj in range(N + 1):
                for i in range(M):
                    for m in range(MR):
                        if i != j:
                            val = p2[j, nj, k, j, nj, k, m]
                            # see _kb/03-api-layer.md for rationale
                            for ni in range(N + 1):  # population 0..F[i], F[i]=N
                                for h in range(K[i]):
                                    val -= p2[j, nj, k, i, ni, h, m]
                            ceq.append(val)

    # UEFF: e[i,ki] = sum p2[j,nj,kj,i,ni,ki,m] for ni>=1 and BB[m,i]==0
    for j in range(M):
        for i in range(M):
            for ki in range(K[i]):
                val = e[i, ki]
                for nj in range(N + 1):
                    for kj in range(K[j]):
                        for m in range(MR):
                            for ni in range(1, N + 1):
                                if BB[m, i] == 0:
                                    val -= p2[j, nj, kj, i, ni, ki, m]
                ceq.append(val)

    # see _kb/03-api-layer.md for rationale
    for i in range(M):
        for k in range(K[i]):
            val = 0.0
            if load_dependent:
                for ni in range(1, F[i] + 1):
                    for m in range(MR):
                        for j in range(M):
                            for h in range(K[i]):
                                val += q[i, j, k, h, ni] * p2[i, ni, k, i, ni, k, m]
                                val -= q[i, j, h, k, ni] * p2[i, ni, h, i, ni, h, m]
            else:
                # LHS
                for j in range(M):
                    for h in range(K[i]):
                        val += q[i, j, k, h] * e[i, k]
                # RHS
                for j in range(M):
                    for h in range(K[i]):
                        val -= q[i, j, h, k] * e[i, h]
            ceq.append(val)

    # THM2: sum_{i,ni,ki} ni*p2[j,nj,k,i,ni,ki,m] = N*p2[j,nj,k,j,nj,k,m]
    for j in range(M):
        for k in range(K[j]):
            for nj in range(F[j] + 1):  # 0 to F[j]
                for m in range(MR):
                    val = 0.0
                    for i in range(M):
                        for ni in range(1, F[i] + 1):  # 1 to F[i]
                            for ki in range(K[i]):
                                val += ni * p2[j, nj, k, i, ni, ki, m]
                    val -= N * p2[j, nj, k, j, nj, k, m]
                    ceq.append(val)

    # COR1: sum ni*nj*p2 = N^2
    val = 0.0
    for m in range(MR):
        for i in range(M):
            for j in range(M):
                for nj in range(1, F[j] + 1):  # 1 to F[j]
                    for ni in range(1, F[i] + 1):  # 1 to F[i]
                        for ki in range(K[i]):
                            for kj in range(K[j]):
                                val += ni * nj * p2[j, nj, kj, i, ni, ki, m]
    val -= N ** 2
    ceq.append(val)

    # see _kb/03-api-layer.md for rationale
    def q_at(i, j, k, h, n):
        """Rate emitted by station i toward j, in phase k->h, at population n.

        The 4D population-free q of qrboundsbas_skel.mod and the 5D
        load-dependent q of qrboundsrsrd_skel.mod differ only in whether the
        emitting station's population is carried, so the two forms of these
        constraints differ only here.
        """
        return q[i, j, k, h, n] if load_dependent else q[i, j, k, h]

    # THM30 {i, u}: balance across the ni = 0 boundary of station i.
    for i in range(M):
        for u in range(K[i]):
            val = 0.0
            for j in range(M):
                if j == i:
                    continue
                for nj in range(1, F[j] + 1):
                    for k in range(K[j]):
                        coef = 0.0
                        for h in range(K[j]):
                            coef += q_at(j, i, k, h, nj)
                        if coef == 0.0:
                            continue
                        for m in range(MR):
                            val += coef * p2[j, nj, k, i, 0, u, m]
            for j in range(M):
                if j == i:
                    continue
                for nj in range(F[j] + 1):
                    for k in range(K[i]):
                        # see _kb/03-api-layer.md for rationale
                        coef = q_at(i, j, k, u, 1)
                        if coef == 0.0:
                            continue
                        for h in range(K[j]):
                            for m in range(MR):
                                val -= coef * p2[j, nj, h, i, 1, k, m]
            ceq.append(val)

    # THM3 {i, ni in 0..F[i]-1}: balance across the ni -> ni+1 boundary.
    for i in range(M):
        for ni in range(F[i]):
            val = 0.0
            for j in range(M):
                if j == i:
                    continue
                for nj in range(1, F[j] + 1):
                    for k in range(K[j]):
                        coef = 0.0
                        for h in range(K[j]):
                            coef += q_at(j, i, k, h, nj)
                        if coef == 0.0:
                            continue
                        for u in range(K[i]):
                            for m in range(MR):
                                val += coef * p2[j, nj, k, i, ni, u, m]
            for j in range(M):
                if j == i:
                    continue
                for nj in range(F[j] + 1):
                    for k in range(K[i]):
                        coef = 0.0
                        for h in range(K[i]):
                            coef += q_at(i, j, k, h, ni + 1)
                        if coef == 0.0:
                            continue
                        for u in range(K[j]):
                            for m in range(MR):
                                val -= coef * p2[j, nj, u, i, ni + 1, k, m]
            ceq.append(val)

    # THM4: inequality (>= in MATLAB, stored as <= with sign swap)
    # sum_{t,h,nj,nt} nt*p2[j,nj,k,t,nt,h,m] >= N*sum_{h,nj,ni} p2[j,nj,k,i,ni,h,m]
    for j in range(M):
        for k in range(K[j]):
            for i in range(M):
                for m in range(MR):
                    val = 0.0
                    # -LHS (sign swapped for <= form)
                    for t in range(M):
                        for h in range(K[t]):
                            for nj_idx in range(N + 1):
                                for nt in range(N + 1):
                                    val -= nt * p2[j, nj_idx, k, t, nt, h, m]
                    # +RHS (sign swapped)
                    for h in range(K[i]):
                        for nj_idx in range(N + 1):
                            for ni in range(1, N + 1):
                                val += N * p2[j, nj_idx, k, i, ni, h, m]
                    c.append(val)

    return np.array(c), np.array(ceq)


def affine_constraint_matrices(fn, n):
    """Recover (A, b) from an affine residual map fn(x) = A x - b.

    Every constraint in the QRF no-blocking inventory is linear in the decision
    vector -- only the objectives (MMI, MEM) are nonlinear -- so the residual
    callbacks are affine and their matrix form is exact, not a linearization.
    Affinity is verified at a probe point and a violation is raised rather than
    tolerated, because a silently non-affine callback would make the recovered
    matrix wrong everywhere except at the probe.

    Args:
        fn: Callable mapping a length-n vector to a residual vector.
        n: Length of the decision vector.

    Returns:
        (A, b) with fn(x) == A @ x - b for all x.
    """
    r0 = np.asarray(fn(np.zeros(n)), dtype=float)
    A = np.zeros((r0.size, n))
    basis = np.zeros(n)
    for col in range(n):
        basis[col] = 1.0
        A[:, col] = np.asarray(fn(basis), dtype=float) - r0
        basis[col] = 0.0
    b = -r0

    probe = np.linspace(0.1, 0.9, n)
    err = np.max(np.abs(np.asarray(fn(probe), dtype=float) - (A @ probe - b))) \
        if r0.size else 0.0
    if err > 1e-9:
        raise ValueError("constraint residuals are not affine in x "
                         "(probe mismatch %.3e); the matrix form recovered "
                         "here would be wrong" % err)
    return A, b


def independent_rows(A, tol=None):
    """Indices of a maximal linearly independent subset of the rows of A.

    Selection is by column-pivoted QR of A.T, so the retained rows are the
    numerically best-conditioned independent set rather than simply the first
    ones encountered.

    This exists because scipy's SLSQP cannot accept more equality constraints
    than variables: its own exit mode 2 is "More equality constraints than
    independent variables". The QRF equality block is heavily redundant --
    SYMMETRY states every station pair twice, and ZERO/MARGINALS/UEFF overlap
    -- so it carries roughly twice as many rows as its rank. Feeding that to
    SLSQP is a formulation error on our side; worse, scipy's Fortran workspace
    is under-allocated in exactly that regime and the solver corrupts the heap
    instead of returning mode 2 (see the module docstring of
    qrf_noblo_mmi_linear). Dropping dependent rows changes no feasible point:
    they are exact linear combinations of the retained ones.

    Args:
        A: 2D array.
        tol: Rank tolerance; default is max(A.shape)*eps*|R[0,0]|.

    Returns:
        Sorted index array of the retained rows.
    """
    A = np.asarray(A, dtype=float)
    if A.shape[0] == 0:
        return np.empty(0, dtype=int)
    _, R, piv = scipy_qr(A.T, mode='economic', pivoting=True)
    diag = np.abs(np.diag(R))
    if diag.size == 0:
        return np.empty(0, dtype=int)
    if tol is None:
        tol = max(A.shape) * np.finfo(float).eps * diag[0]
    rank = int(np.sum(diag > tol))
    return np.sort(piv[:rank])


def reduce_equalities(A, b):
    """Drop linearly dependent equality rows, keeping the feasible set exact.

    Raises:
        ValueError: If the dropped rows are not implied by the retained ones,
            i.e. the system is inconsistent (rank([A|b]) > rank(A)). That is a
            modelling error and must not be silently discarded.
    """
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    keep = independent_rows(A)
    if keep.size == A.shape[0]:
        return A, b, keep
    aug_rank = np.linalg.matrix_rank(np.hstack([A, b[:, None]]))
    if aug_rank > keep.size:
        raise ValueError("equality system is inconsistent: rank(A)=%d but "
                         "rank([A|b])=%d" % (keep.size, aug_rank))
    return A[keep], b[keep], keep


def feasible_start(Aeq, beq, Aub, bub, n):
    """A point of the polytope, for use as the NLP start point.

    The all-zero vector violates ONE (normalization) by a full unit and COR1 by
    N^2, and from there SLSQP terminates with mode 4, "Inequality constraints
    incompatible", returning the start point unchanged: the reported UN and QN
    are then identically zero and nothing in the result signals it. Every
    constraint of this model is linear, so a phase-1 LP over the same matrices
    lands on the polytope exactly and costs a fraction of one NLP iteration.

    Raises:
        ValueError: If the LP is infeasible. The polytope is nonempty for every
            well-posed instance, so an empty one is a modelling error and must
            surface rather than be replaced by an arbitrary point.
    """
    has_ineq = Aub is not None and len(Aub) > 0
    res = linprog(np.zeros(n),
                  A_ub=Aub if has_ineq else None,
                  b_ub=bub if has_ineq else None,
                  A_eq=Aeq, b_eq=beq,
                  bounds=[(0.0, 1.0)] * n, method='highs')
    if not res.success:
        raise ValueError("QRF no-blocking polytope is infeasible "
                         "(phase-1 LP status %d: %s)" % (res.status, res.message))
    return res.x


def compute_num_vars(M, N, K, MR):
    """Compute total number of decision variables."""
    Kmax = int(max(K))
    return M * (N + 1) * Kmax * M * (N + 1) * Kmax * MR + M * Kmax


def extract_results(p2opt, M, K, F, MR):
    """Extract UN and QN from optimal p2 tensor.

    Args:
        p2opt: Optimal 7D probability tensor
        M: Number of queues
        K: Phases per queue [M]
        F: Capacity per queue [M]
        MR: Number of blocking configurations

    Returns:
        UN: Utilization per queue [M]
        QN: Queue length per queue [M]
    """
    UN = np.zeros(M)
    QN = np.zeros(M)
    for ti in range(M):
        for m in range(MR):
            for ni in range(1, F[ti] + 1):
                for ki in range(K[ti]):
                    UN[ti] += p2opt[ti, ni, ki, ti, ni, ki, m]
                    QN[ti] += ni * p2opt[ti, ni, ki, ti, ni, ki, m]
    return UN, QN


def build_q_from_mu_v_rt(M, K, mu, v, rt):
    """Build 4D transition rate array from mu, v, rt.

    Args:
        M: Number of queues
        K: Phases per queue [M]
        mu: Completion rates [M, Kmax, Kmax]
        v: Background rates [M, Kmax, Kmax]
        rt: Routing matrix [M, M]

    Returns:
        q: Transition rates [M, M, Kmax, Kmax]
    """
    Kmax = int(max(K))
    q = np.zeros((M, M, Kmax, Kmax))
    for i in range(M):
        for j in range(M):
            for k in range(K[i]):
                for h in range(K[i]):
                    if j != i:
                        q[i, j, k, h] = rt[i, j] * mu[i, k, h]
                    else:
                        q[i, j, k, h] = v[i, k, h] + rt[i, i] * mu[i, k, h]
    return q


def build_q_ld(M, K, mu, v, rt, N, alpha=None):
    """Build the 5D load-dependent transition rate array.

    qrboundsrsrd_skel.mod:11 declares q with FIVE indices,

        q {i in 1..M, j in 1..M, k in 1..K[i], h in 1..K[i], n in 0..N}

    where n is the population of station i, the station EMITTING the
    transition, q[...,n=0] = 0, and the load-dependent scaling alpha[i,n]
    multiplies BOTH the background term v and the completion term
    r[i,i]*mu[i,k,h]. Dropping the population index (or dropping alpha from
    the v term) is what makes a load-dependent instance unrepresentable: with
    a population-free q the balance constraints cannot distinguish the rate at
    which station i empties at population n from the rate at population n', so
    the polytope stops pinning the utilization.

    The 4D form built by build_q_from_mu_v_rt is the correct one for the
    non-load-dependent variants, whose skeleton (qrboundsbas_skel.mod) declares
    q without the population index.

    Args:
        M: Number of queues
        K: Phases per queue [M]
        mu: Completion rates [M, Kmax, Kmax]
        v: Background rates [M, Kmax, Kmax]
        rt: Routing matrix [M, M]
        N: Total population
        alpha: Load-dependent scaling [M, N], 0-based in the population so that
            alpha[i, n-1] is the AMPL alpha[i,n]. None means all ones.

    Returns:
        q: Transition rates [M, M, Kmax, Kmax, N+1], indexed by the ACTUAL
           population n in 0..N (no MATLAB-style 1+n shift).
    """
    Kmax = int(max(K))
    if alpha is None:
        alpha = np.ones((M, N))
    alpha = np.asarray(alpha, dtype=float)
    q = np.zeros((M, M, Kmax, Kmax, N + 1))
    for i in range(M):
        for j in range(M):
            for k in range(K[i]):
                for h in range(K[i]):
                    for n in range(1, N + 1):
                        a = alpha[i, n - 1]
                        if j != i:
                            q[i, j, k, h, n] = rt[i, j] * mu[i, k, h] * a
                        else:
                            q[i, j, k, h, n] = (v[i, k, h] * a
                                                + rt[i, i] * mu[i, k, h] * a)
    return q


def extract_mu_v_from_maps(MAPs, M, K):
    """Extract mu and v arrays from MAPs cell array.

    Args:
        MAPs: List of [D0, D1] pairs per queue
        M: Number of queues
        K: Phases per queue [M]

    Returns:
        mu: Completion rates [M, Kmax, Kmax]
        v: Background rates [M, Kmax, Kmax]
    """
    Kmax = int(max(K))
    mu = np.zeros((M, Kmax, Kmax))
    v = np.zeros((M, Kmax, Kmax))
    for i in range(M):
        D0 = MAPs[i][0]
        D1 = MAPs[i][1]
        for h in range(K[i]):
            for k in range(K[i]):
                mu[i, h, k] = D1[h, k]
                if h == k:
                    v[i, k, h] = 0.0
                else:
                    v[i, k, h] = D0[h, k]
    return mu, v
