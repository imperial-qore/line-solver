"""
Shared utilities for QRF NLP no-blocking approximation methods.

Provides variable unflattening (sub_qrfvar), objective functions (MMI, MEM,
BETHE),
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

    D1 FIX (2026-09-02): the population sums run from n = 0. They used to
    start at 1, dropping the idle/idle cells -- and since U_i = 1 - p_i(0),
    that excluded the strongest correlation in a closed chain from the very
    functional meant to measure coupling. The AMPL source reads
    `sum {ni, nj in 0..F}`. Note the structurally zero entries now inside
    the sum contribute 0*log(LOGTOL) = 0 to the VALUE but log(LOGTOL) to the
    GRADIENT, so the descent direction is LOGTOL-sensitive where the value is
    not; the optimum moves < 1e-7 over LOGTOL in [1e-8, 1e-4].
    """
    p2, _ = sub_qrfvar(x, M, N, K, MR)
    fobj = 0.0
    for m in range(MR):
        for i in range(M):
            for ki in range(K[i]):
                for j in range(M):
                    if i != j:
                        for kj in range(K[j]):
                            for ni in range(0, F[i] + 1):
                                for nj in range(0, F[j] + 1):
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
    """Maximum-entropy objective, as the NEGATIVE entropy a minimiser wants.

    Returns: +sum p2[i,ni,k,i,ni,k,m] * log(p2[i,ni,k,i,ni,k,m])

    The AMPL model this comes from states the objective as `maximize H` with
    H = -sum p log p, and every port hands it to a MINIMISER, so the entropy
    must be returned negated. Returning +H (as all four ports did until
    2026-08-29) selects the MINIMUM-entropy point of the polytope -- the
    opposite face -- under a method advertised as maximum-entropy.

    Port of MATLAB mem() nested function.
    """
    p2, _ = sub_qrfvar(x, M, N, K, MR)
    fobj = 0.0
    for m in range(MR):
        for i in range(M):
            for k in range(K[i]):
                for ni in range(1, F[i] + 1):
                    pval = p2[i, ni, k, i, ni, k, m]
                    fobj += pval * np.log(LOGTOL + pval)
    return fobj


def bethe_objective(x, M, N, K, F, MR):
    """Tree-reweighted (Bethe) free entropy at the uniform spanning-tree weight.

    With lambda = 1/M this is

        lambda * sum_m sum_{i!=j} sum_{ki,kj} sum_{ni,nj>=0}
                     p_ij*(log p_ij - log p_ii - log p_jj)
              + sum_m sum_i sum_k sum_{n>=0} p_ii*log p_ii

    i.e. lambda*sum_{i!=j} I(n_i;n_j) - sum_i H(n_i), the NEGATIVE of a
    tree-reweighted entropy with uniform edge weight rho_ij = 2*lambda on the
    complete station graph.

    WHY lambda = 1/M AND NOT THE BETHE 1/2. H_rho is a convex combination of
    tree entropies, hence concave on the local marginal polytope, exactly when
    rho lies in the spanning tree polytope of K_M. Its uniform point is
    rho_ij = 2/M, so lambda = 1/M is the LARGEST uniform weight for which
    minimising this is a CONVEX program -- every local optimum global, the
    answer a property of the model rather than of the start point. The Bethe
    weight lambda = 1/2 (total edge mass C(M,2) against the M-1 a spanning tree
    can carry) is outside it for every M > 2 and coincides with 1/M at M = 2.

    TWO DIFFERENCES FROM mmi_objective, BOTH DELIBERATE. The population loops
    start at n = 0, the range the AMPL source states (`ni, nj in 0..F`) and the
    one mmi_objective does not use, so the idle/idle cell -- the strongest
    correlation in a closed chain -- is inside the sum; and the entropy term is
    mem_objective's body over the same restored range, which already carries
    the sign a minimiser needs. Neither repair touches qrf.mmi or qrf.mem,
    whose values are pinned by tests.

    NUMERICAL NOTE. Restoring the n = 0 cells brings the structurally zero
    entries inside the sum: they contribute 0*log(LOGTOL) = 0 to the VALUE but
    log(LOGTOL) ~ -13.8 to the GRADIENT, so the value is insensitive to LOGTOL
    while the descent direction is not. Measured over LOGTOL in
    {1e-4, 1e-6, 1e-8} the optimum moves by under 1e-7 in U on the fixtures of
    test_qrf_noblo_bethe.py, but do not freeze test digits below that.
    """
    p2, _ = sub_qrfvar(x, M, N, K, MR)
    lam = 1.0 / M
    fobj = 0.0
    for m in range(MR):
        for i in range(M):
            for ki in range(K[i]):
                for j in range(M):
                    if i != j:
                        for kj in range(K[j]):
                            for ni in range(0, F[i] + 1):
                                for nj in range(0, F[j] + 1):
                                    pij = p2[i, ni, ki, j, nj, kj, m]
                                    pii = p2[i, ni, ki, i, ni, ki, m]
                                    pjj = p2[j, nj, kj, j, nj, kj, m]
                                    fobj += lam * pij * (
                                        np.log(LOGTOL + pij)
                                        - np.log(LOGTOL + pii)
                                        - np.log(LOGTOL + pjj)
                                    )
    for m in range(MR):
        for i in range(M):
            for k in range(K[i]):
                for ni in range(0, F[i] + 1):
                    pval = p2[i, ni, k, i, ni, k, m]
                    fobj += pval * np.log(LOGTOL + pval)
    return fobj


def qrf_index_map(M, N, K, MR):
    """Flat position of every p2 entry, in the fill order of sub_qrfvar.

    The layout is not a plain strided tensor: the phase loops run over K[j] and
    K[i], not over Kmax, so a station with fewer phases leaves gaps. Entries
    that carry no variable keep -1, and the gradient scatter must skip them.
    """
    Kmax = int(max(K))
    idx = -np.ones((M, N + 1, Kmax, M, N + 1, Kmax, MR), dtype=np.int64)
    ctr = 0
    for j in range(M):
        for nj in range(N + 1):
            for k in range(K[j]):
                for i in range(M):
                    for ni in range(N + 1):
                        for h in range(K[i]):
                            for m in range(MR):
                                idx[j, nj, k, i, ni, h, m] = ctr
                                ctr += 1
    return idx


def mmi_gradient(x, M, N, K, F, MR, idx=None):
    """Gradient of mmi_objective.

    For a term t = pij (log pij - log pii - log pjj) the three partials are
    dt/dpij = log pij - log pii - log pjj + pij/pij', dt/dpii = -pij/pii' and
    dt/dpjj = -pij/pjj', a primed denominator standing for the LOGTOL-shifted
    value. i != j throughout, so no term aliases its own partials.
    """
    if idx is None:
        idx = qrf_index_map(M, N, K, MR)
    p2, _ = sub_qrfvar(x, M, N, K, MR)
    g = np.zeros_like(np.asarray(x, dtype=float))
    for m in range(MR):
        for i in range(M):
            ri = np.arange(0, F[i] + 1)          # D1: from n = 0
            for ki in range(K[i]):
                pii = p2[i, ri, ki, i, ri, ki, m]
                log_pii = np.log(LOGTOL + pii)
                idx_ii = idx[i, ri, ki, i, ri, ki, m]
                for j in range(M):
                    if i == j:
                        continue
                    rj = np.arange(0, F[j] + 1)          # D1: from n = 0
                    for kj in range(K[j]):
                        pjj = p2[j, rj, kj, j, rj, kj, m]
                        log_pjj = np.log(LOGTOL + pjj)
                        idx_jj = idx[j, rj, kj, j, rj, kj, m]
                        pij = p2[i, ri[:, None], ki, j, rj[None, :], kj, m]
                        idx_ij = idx[i, ri[:, None], ki, j, rj[None, :], kj, m]
                        gij = (np.log(LOGTOL + pij) - log_pii[:, None]
                               - log_pjj[None, :] + pij / (LOGTOL + pij))
                        np.add.at(g, idx_ij.ravel(), gij.ravel())
                        np.add.at(g, idx_ii, -pij.sum(axis=1) / (LOGTOL + pii))
                        np.add.at(g, idx_jj, -pij.sum(axis=0) / (LOGTOL + pjj))
    return g


def mem_gradient(x, M, N, K, F, MR, idx=None):
    """Gradient of mem_objective: d/dp of p log(p') is log p' + p/p'."""
    if idx is None:
        idx = qrf_index_map(M, N, K, MR)
    p2, _ = sub_qrfvar(x, M, N, K, MR)
    g = np.zeros_like(np.asarray(x, dtype=float))
    for m in range(MR):
        for i in range(M):
            ri = np.arange(1, F[i] + 1)
            for k in range(K[i]):
                p = p2[i, ri, k, i, ri, k, m]
                np.add.at(g, idx[i, ri, k, i, ri, k, m],
                          np.log(LOGTOL + p) + p / (LOGTOL + p))
    return g


def bethe_gradient(x, M, N, K, F, MR, idx=None):
    """Gradient of bethe_objective.

    df/dp_ij = lam*(log p_ij' + p_ij/p_ij' - log p_ii' - log p_jj') for i != j,
    df/dp_ii = log p_ii' + p_ii/p_ii'
               - (lam/p_ii') * sum_{j!=i} sum_{kj,nj} (p_ij + p_ji),

    a primed denominator standing for the LOGTOL-shifted value. The second sum
    is accumulated by the scatter below, which visits both orderings of every
    pair: (i,j) contributes p_ij to idx_ii and (j,i) contributes p_ji.
    """
    if idx is None:
        idx = qrf_index_map(M, N, K, MR)
    p2, _ = sub_qrfvar(x, M, N, K, MR)
    lam = 1.0 / M
    g = np.zeros_like(np.asarray(x, dtype=float))
    for m in range(MR):
        for i in range(M):
            ri = np.arange(0, F[i] + 1)
            for ki in range(K[i]):
                pii = p2[i, ri, ki, i, ri, ki, m]
                log_pii = np.log(LOGTOL + pii)
                idx_ii = idx[i, ri, ki, i, ri, ki, m]
                for j in range(M):
                    if i == j:
                        continue
                    rj = np.arange(0, F[j] + 1)
                    for kj in range(K[j]):
                        pjj = p2[j, rj, kj, j, rj, kj, m]
                        log_pjj = np.log(LOGTOL + pjj)
                        idx_jj = idx[j, rj, kj, j, rj, kj, m]
                        pij = p2[i, ri[:, None], ki, j, rj[None, :], kj, m]
                        idx_ij = idx[i, ri[:, None], ki, j, rj[None, :], kj, m]
                        gij = lam * (np.log(LOGTOL + pij) - log_pii[:, None]
                                     - log_pjj[None, :] + pij / (LOGTOL + pij))
                        np.add.at(g, idx_ij.ravel(), gij.ravel())
                        np.add.at(g, idx_ii, -lam * pij.sum(axis=1) / (LOGTOL + pii))
                        np.add.at(g, idx_jj, -lam * pij.sum(axis=0) / (LOGTOL + pjj))
    for m in range(MR):
        for i in range(M):
            ri = np.arange(0, F[i] + 1)
            for k in range(K[i]):
                p = p2[i, ri, k, i, ri, k, m]
                np.add.at(g, idx[i, ri, k, i, ri, k, m],
                          np.log(LOGTOL + p) + p / (LOGTOL + p))
    return g


def solve_qrf_nlp(objective, gradient, x0, Aeq, beq, Aub, bub, name):
    """Minimise a QRF no-blocking objective over its polytope.

    The equalities are ELIMINATED, not passed to the solver. They are nearly
    square -- 219 independent rows for 228 variables on a 3-station, N=4
    instance -- so the feasible set is a handful of dimensions inside a large
    ambient space. In the FULL space SLSQP stops at iteration 1 with mode 4,
    "Inequality constraints incompatible", and returns the start point: every
    token backed by it reported the phase-1 LP vertex, identically for all four
    objectives, with nothing in the result to say so. Parameterising
    x = x0 + Z t with Z an orthonormal basis of null(Aeq) leaves an
    unconstrained-in-equalities problem of about nine dimensions, on which the
    same SLSQP converges in 14 iterations and under a second. (trust-constr
    also solves the reduced problem, to the same objective, but spends ~2 s per
    iteration on the mapped box rows: 95 s against 0.6 s.)

    Raises:
        RuntimeError: If the solver makes no progress at all. Returning the
            start point as though it were a solution is what this function
            exists to prevent.
    """
    from scipy.linalg import null_space
    from scipy.optimize import minimize
    x0 = np.asarray(x0, dtype=float)
    n = len(x0)
    Aeq = np.asarray(Aeq, dtype=float)
    Z = null_space(Aeq) if Aeq.size else np.eye(n)
    d = Z.shape[1]
    if d == 0:
        # The equalities pin a single point; there is nothing left to optimise
        # and x0 IS the answer, not a failure to optimise.
        return x0

    # The objective is p log p, so a trial point with p < -LOGTOL would take the
    # log of a negative number. It is evaluated on the CLIPPED point, and the
    # clipped coordinates carry zero sensitivity.
    def obj_t(t):
        return objective(np.clip(x0 + Z @ t, 0.0, 1.0))

    def grad_t(t):
        x = x0 + Z @ t
        g = gradient(np.clip(x, 0.0, 1.0))
        g[(x < 0.0) | (x > 1.0)] = 0.0
        return Z.T @ g

    # Box rows whose Z row vanishes are coordinates the equalities PIN (the
    # model pins most p2 entries to zero); in t they read 0 <= 0 <= 0, which
    # constrains nothing and only degrades the active-set factorisation.
    row_norm = np.linalg.norm(Z, axis=1)
    varies = row_norm > 1e-12 * max(1.0, float(np.max(row_norm)))
    pinned = ~varies
    if np.any(pinned) and (np.min(x0[pinned]) < -1e-9 or np.max(x0[pinned]) > 1 + 1e-9):
        raise RuntimeError("%s: the equality system pins a variable outside "
                           "[0,1]; the polytope is empty." % name)
    A = np.vstack([Z[varies], -Z[varies]])
    b = np.concatenate([1.0 - x0[varies], x0[varies]])
    if Aub is not None and len(Aub) > 0:
        Aub = np.asarray(Aub, dtype=float)
        AubZ = Aub @ Z
        rhs = np.asarray(bub, dtype=float) - Aub @ x0
        keep = np.linalg.norm(AubZ, axis=1) > 1e-12
        if np.any(~keep) and np.max(rhs[~keep]) < -1e-9:
            raise RuntimeError("%s: an inequality independent of the free "
                               "directions is violated at the feasible start."
                               % name)
        if np.any(keep):
            A = np.vstack([A, AubZ[keep]])
            b = np.concatenate([b, rhs[keep]])

    constraints = [{'type': 'ineq',
                    'fun': lambda t: b - A @ t,
                    'jac': lambda t: -A}]
    def _run(t_start):
        return minimize(obj_t, t_start, jac=grad_t, method='SLSQP',
                        constraints=constraints,
                        options={'maxiter': 300, 'ftol': 1e-12, 'disp': False})

    result = _run(np.zeros(d))

    # ESCAPE A STATIONARY POINT OF A CONCAVE OBJECTIVE.
    #
    # SLSQP stopping without moving is not by itself a defect: on these
    # polytopes it commonly lands on a point where NO feasible direction is
    # descending TO FIRST ORDER, and it reports status 0 there. The LP probe
    # below asks exactly that first-order question, and on the ragged
    # two-phase / one-phase instance at N = 2 it answers 0.0 -- the start is
    # genuinely first-order stationary, so the guard correctly declines to
    # fire and the old code returned the phase-1 vertex.
    #
    # But the MI objective is CONCAVE on the phase-type polytopes, and for a
    # concave f a point with zero directional derivative is a local MAXIMUM or
    # a saddle, not a minimum: f(t + s*u) <= f(t) + s*<g,u> = f(t) for every
    # feasible u. A first-order test is therefore structurally blind exactly
    # where it matters, and no tightening of it can help. The minimisers of a
    # concave program are VERTICES, so the escape has to take a FINITE step to
    # one and compare values, which is what a conditional-gradient method does
    # implicitly and what the C++ port has always done.
    #
    # Measured on that instance: the start has f = 1.900535 and this escape
    # reaches 0.824607, the value the C++ Frank-Wolfe port reports; both are
    # feasible, and U_0 = 0.8 and 1.0 are respectively the LP minimum and the
    # LP maximum of the same polytope.
    #
    # Deterministic: the probe costs are drawn from a fixed seed, so the answer
    # is a property of the model and not of the run.
    rng = np.random.default_rng(20260902)
    bounds_t = [(-1.0, 1.0)] * d
    for _ in range(8):
        t_cur = result.x
        f_cur = float(result.fun)
        gred = grad_t(t_cur)
        probe = linprog(gred, A_ub=A, b_ub=b, bounds=bounds_t, method='highs')
        moved = np.max(np.abs(t_cur)) > 0.0
        # A feasible descent direction exists and SLSQP did not move, so the
        # start would be reported as the optimum. That is not yet a failure:
        # the phase-1 LP returns a VERTEX, and SLSQP declines a degenerate one
        # with status 4 ("Inequality constraints incompatible") while a finite
        # step along a segment into the polytope descends perfectly well. The
        # finite step is below, and it is what the C++ port does at every
        # iterate rather than once at the start. Only an escape that cannot
        # improve either makes this a failure, which is where it is raised.
        stuck = (probe.success and not moved
                 and probe.fun - float(gred @ t_cur) < -1e-8)

        # Candidate vertices: the steepest-descent one, plus random-cost ones
        # for the case the reduced gradient vanishes and the first carries no
        # information.
        cands = []
        if probe.success:
            cands.append(np.asarray(probe.x, dtype=float))
        for _k in range(2 * d + 4):
            v = linprog(rng.standard_normal(d), A_ub=A, b_ub=b,
                        bounds=bounds_t, method='highs')
            if v.success:
                cands.append(np.asarray(v.x, dtype=float))

        best_t, best_f = None, f_cur
        for v in cands:
            for gam in (1.0, 0.5, 0.25, 0.1, 0.01, 1e-3):
                t_try = t_cur + gam * (v - t_cur)
                f_try = obj_t(t_try)
                if np.isfinite(f_try) and f_try < best_f - 1e-10 * (1.0 + abs(best_f)):
                    best_t, best_f = t_try, float(f_try)
        if best_t is None:
            if stuck:
                raise RuntimeError(
                    "%s: the NLP did not move from the feasible start although "
                    "a feasible descent direction exists there (directional "
                    "derivative %.3e, status %d: %s), and no finite step to a "
                    "vertex of the polytope improved the objective either. The "
                    "reported metrics would be those of the phase-1 LP vertex, "
                    "not of the QRF optimum."
                    % (name, probe.fun, result.status, result.message))
            break
        result = _run(best_t)
        if float(result.fun) > best_f:      # keep the better of the two
            class _R:
                pass
            r = _R()
            r.x, r.fun, r.status, r.message = best_t, best_f, 0, 'escape step'
            result = r

    return np.clip(x0 + Z @ result.x, 0.0, 1.0)


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


def extract_busy(p2opt, M, K, F, MR, alpha):
    """Alpha-weighted diagonal marginal mean: the mean number of jobs in service.

    ``E[min(n,c)]`` at a c-server station, ``E[n]`` at a delay and ``P(n >= 1)``
    where alpha is 1. It is what the departure rate is proportional to, since
    alpha(i,n) scales the completion rate, so a station's throughput is
    ``BN / stime`` exactly at the relaxed point. See ``sn_to_qrf_alpha``.

    Args:
        p2opt: Optimal 7D probability tensor
        M: Number of queues
        K: Phases per queue [M]
        F: Capacity per queue [M]
        MR: Number of blocking configurations
        alpha: Load-dependent scaling [M, N], indexed by population n = 1..N

    Returns:
        BN: Mean number in service per queue [M]
    """
    alpha = np.asarray(alpha, dtype=float)
    BN = np.zeros(M)
    for ti in range(M):
        for m in range(MR):
            for ni in range(1, F[ti] + 1):
                a = alpha[ti, ni - 1]
                for ki in range(K[ti]):
                    BN[ti] += a * p2opt[ti, ni, ki, ti, ni, ki, m]
    return BN


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
                # (from, to), the order mu and the q assembly already use: writing
                # v[i, k, h] transposes D0 alone and reverses an Erlang.
                if h == k:
                    v[i, h, k] = 0.0
                else:
                    v[i, h, k] = D0[h, k]
    return mu, v
