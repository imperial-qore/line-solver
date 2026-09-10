"""
Per-type waiting times of the MMAP[K]/G[K]/1 FCFS queue.

THE METHOD, which is He's, theorem for theorem. FCFS makes the actual waiting
time of a customer the WORKLOAD it finds on arrival, so everything follows from
the joint transform of workload and arrival phase,
f(s)_j = E[exp(-s V) 1{phase = j}], which by He's Theorem 4.1 (eq. 4.6)
satisfies

    f(s) [ s I + D0 + sum_k Dk gk(s) ] = s v0,                          (*)

with v0 the idle-phase vector, his y0. The unknown v0 needs NO search for the
roots of the determinant: the matrix U solving

    U = D0 + sum_k Dk Fk(U),      Fk(U) = int_0^inf exp(U t) dFk(t),

is his eq. (4.4), the generator of the underlying Markov process obtained by
EXCISING the busy periods, and eq. (4.5) with Theorem 4.2 give y0 Q = 0 and
y0 e = 1 - rho, i.e. v0 = (1 - rho) pi_U.

The same vector is what the analyticity of (*) forces: for every left eigenpair
(w, u) of U one has w [D0 + sum_k Dk gk(-u) - u I] = 0, so the roots of (*) in
the closed right half plane are exactly s = -u over the spectrum of U, and
imposing v0 r_i = 0 at each right null vector reproduces the stationary vector
to 2.5e-13. The stationary route is the one taken, as it needs no complex
eigenvector and no rule for telling the structural root at the origin from a
genuine one.

The per-type actual waiting time is the workload seen by a type-k arrival,
biased by that type's own arrival block, which is his Theorem 5.1 eq. (5.1)
summed over the post-arrival phase:

    E[exp(-s Wk)] = f(s) Dk e / lambda_k.

SCOPE. He allows an arrival to be a BATCH carrying a sequence of types, and his
Theorem 5.3 then multiplies the transform by prod_{i<n} f*_{h_i}(s), the service
of the customers ahead of the tagged one WITHIN its own batch. This module
covers the single-customer-per-arrival case, his Special case 3.3, where that
product is empty -- which is exactly the MMAP convention LINE carries,
[D0, D1, D^(1), ..., D^(K)] with one customer per epoch.

Reference:
    Qi-Ming He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1 queue",
    Queueing Systems 38(4):397-418, 2001.
    Original MATLAB: matlab/src/api/qsys/qsys_mmapgk1.m
"""

from dataclasses import dataclass
from math import comb
from typing import Any, Callable, List, Optional, Sequence

import numpy as np
from scipy.linalg import expm, null_space


@dataclass
class MmapGk1Result:
    """Result structure for the MMAP[K]/G[K]/1 analysis."""
    lambdas: np.ndarray
    arrivalRate: float
    utilization: float
    idleVector: np.ndarray
    waitLST: Callable[[complex], np.ndarray]
    waitMoments: np.ndarray
    meanWaitingTime: np.ndarray
    meanSojournTime: np.ndarray
    meanQueueLength: float
    waitCDF: Optional[np.ndarray]
    waitPoints: Optional[np.ndarray]
    analyzer: str


def _stat_left_null(G: np.ndarray) -> np.ndarray:
    n = G.shape[0]
    B = np.hstack([G, np.ones((n, 1))])
    y = np.zeros(n + 1)
    y[n] = 1.0
    x, _, _, _ = np.linalg.lstsq(B.T, y, rcond=None)
    return x


def _is_ph_pair(law: Any) -> bool:
    return isinstance(law, (list, tuple)) and len(law) >= 2 and np.ndim(law[0]) == 2


def _is_transform_law(law: Any) -> bool:
    """
    A law given as its TRANSFORM plus raw moments, {'lst': f, 'moments': [...]}.

    This is what a solver hands in: sn.lst already carries the transform of the
    ORIGINAL law, while sn.proc carries only its phase-type fit, and it is the
    original this analysis needs. The handle must admit a COMPLEX argument, since
    the matrix transform is read off the spectrum of U.
    """
    return isinstance(law, dict) and 'lst' in law


def _svc_mean(law: Any) -> float:
    if _is_transform_law(law):
        return float(law['moments'][0])
    if _is_ph_pair(law):
        D0s = np.atleast_2d(np.asarray(law[0], dtype=float))
        n = D0s.shape[0]
        alpha = _ph_alpha(law)
        return float(-alpha @ np.linalg.solve(D0s, np.ones(n)))
    return float(law.getMean())


def _ph_alpha(law: Any) -> np.ndarray:
    D0s = np.atleast_2d(np.asarray(law[0], dtype=float))
    D1s = np.atleast_2d(np.asarray(law[1], dtype=float))
    n = D0s.shape[0]
    # embedded stationary vector at arrival epochs
    P = np.linalg.solve(-D0s, D1s)
    w, V = np.linalg.eig(P.T)
    i = int(np.argmin(np.abs(w - 1.0)))
    v = np.real(V[:, i])
    return v / np.sum(v)


def _raw_moment(law: Any, j: int) -> float:
    if _is_transform_law(law):
        moms = law['moments']
        if j > len(moms):
            raise ValueError(
                "The service law was given as a transform with %d moments, but order %d "
                "is needed; supply at least num_w_moms+1 of them" % (len(moms), j))
        return float(moms[j - 1])
    if _is_ph_pair(law):
        D0s = np.atleast_2d(np.asarray(law[0], dtype=float))
        n = D0s.shape[0]
        alpha = _ph_alpha(law)
        M = np.linalg.inv(-D0s)
        acc = np.eye(n)
        for _ in range(j):
            acc = acc @ M
        from math import factorial
        return float(factorial(j) * alpha @ acc @ np.ones(n))
    name = type(law).__name__
    if name == 'Det':
        return float(law.getMean()) ** j
    # A Markovian law has the exact moment j! alpha (-D0)^-j e; the quadrature
    # fallback below is only for a law with no phase-type representation.
    if hasattr(law, 'getD0') and hasattr(law, 'getInitProb'):
        try:
            D0s = np.atleast_2d(np.asarray(law.getD0(), dtype=float))
            alpha = np.asarray(law.getInitProb(), dtype=float).ravel()
            n = D0s.shape[0]
            if alpha.size == n and np.all(np.isfinite(D0s)):
                return _raw_moment([D0s, np.outer(-(D0s @ np.ones(n)), alpha)], j)
        except Exception:
            pass
    if name == 'Uniform':
        a, b = (float(v) for v in law.getSupport())
        return (b ** (j + 1) - a ** (j + 1)) / ((b - a) * (j + 1))
    if j == 1:
        return float(law.getMean())
    if j == 2:
        return float(law.getVar() + law.getMean() ** 2)
    if j == 3:
        try:
            m1 = float(law.getMean())
            v = float(law.getVar())
            m2 = v + m1 ** 2
            return float(law.getSkewness()) * v ** 1.5 + 3 * m1 * m2 - 2 * m1 ** 3
        except Exception:
            pass
    # A HEAVY TAIL HAS NO MOMENT, and the truncated sum below would hand back a
    # finite number for one that diverges: a Pareto of shape <= j has
    # E[S^j] = Inf, and E[Wq] is then genuinely infinite rather than merely
    # large. The quadrature integrates over a cut support and cannot see that,
    # so the divergence is decided from the tail index first.
    if name == 'Pareto':
        shape = getattr(law, '_alpha', None)
        if shape is not None and float(shape) <= j:
            return float('inf')
    x, w = _stieltjes_nodes(law)
    return float(np.sum(w * x ** j))


def _stieltjes_nodes(law: Any, n_grid: int = 2400):
    """Midpoint nodes with true CDF increments, a proper measure for any law."""
    lo, hi = 0.0, float('inf')
    try:
        lo, hi = (float(v) for v in law.getSupport())
    except Exception:
        pass
    if not np.isfinite(hi):
        hi = float(law.getMean()) * 60.0
        try:
            hi = max(hi, float(law.getMean()) + 12.0 * float(law.getVar()) ** 0.5)
        except Exception:
            pass
    edges = np.linspace(lo, hi, n_grid + 1)
    x = 0.5 * (edges[:-1] + edges[1:])
    cdf = np.array([float(law.evalCDF(t)) for t in edges])
    w = np.diff(cdf)
    mass = np.sum(w)
    if mass > 0:
        w = w / mass
    return x, w


def _scalar_lst(law: Any, s: complex) -> complex:
    if _is_transform_law(law):
        return complex(law['lst'](s))
    if _is_ph_pair(law):
        D0s = np.atleast_2d(np.asarray(law[0], dtype=float))
        n = D0s.shape[0]
        alpha = _ph_alpha(law)
        s0 = -D0s @ np.ones(n)
        return complex(alpha @ np.linalg.solve(s * np.eye(n) - D0s, s0))
    try:
        return complex(law.evalLST(s))
    except Exception:
        x, w = _stieltjes_nodes(law)
        return complex(np.sum(w * np.exp(-s * x)))


def _matrix_lst(law: Any, U: np.ndarray) -> np.ndarray:
    """int_0^inf exp(U t) dF(t)."""
    n = U.shape[0]
    if _is_transform_law(law):
        # The transform handle alone suffices: diagonalizing U turns the MATRIX
        # transform into the SCALAR one at the eigenvalues, which is why sn.lst
        # has to admit a complex argument.
        uv, Vd = np.linalg.eig(U)
        gv = np.array([complex(law['lst'](-u)) for u in uv])
        return np.real(Vd @ np.diag(gv) @ np.linalg.inv(Vd))
    if _is_ph_pair(law):
        # The density is the SCALAR beta exp(St) s0, so the integral is exact on
        # the Kronecker sum: int exp(Ut) x exp(St) dt = -(U (+) S)^-1.
        D0s = np.atleast_2d(np.asarray(law[0], dtype=float))
        ms = D0s.shape[0]
        beta = _ph_alpha(law)
        s0 = -D0s @ np.ones(ms)
        KS = np.kron(U, np.eye(ms)) + np.kron(np.eye(n), D0s)
        return np.kron(np.eye(n), beta.reshape(1, -1)) @ (
            -np.linalg.solve(KS, np.kron(np.eye(n), s0.reshape(-1, 1))))
    if type(law).__name__ == 'Det':
        return expm(U * float(law.getMean()))
    # A GENERAL law needs no quadrature: diagonalizing U turns the matrix
    # transform into the SCALAR transform at the eigenvalues.
    try:
        uv, Vd = np.linalg.eig(U)
        if np.linalg.cond(Vd) < 1e12:
            gv = np.array([_scalar_lst(law, -u) for u in uv])
            if np.all(np.isfinite(gv)):
                return np.real(Vd @ np.diag(gv) @ np.linalg.inv(Vd))
    except Exception:
        pass
    x, w = _stieltjes_nodes(law)
    F = np.zeros((n, n))
    for xi, wi in zip(x, w):
        F += wi * expm(U * xi)
    return F


def _euler_invert(lst_fun: Callable[[complex], complex], t: float) -> float:
    """Abate-Whitt Euler inversion of a CDF from its Laplace-Stieltjes transform."""
    if t <= 0:
        return float(np.real(lst_fun(1e12)))
    A = 18.4
    n_euler, m_euler = 15, 11
    u = np.exp(A / 2) / t
    x = A / (2 * t)
    terms = np.zeros(n_euler + m_euler + 1)
    terms[0] = float(np.real(lst_fun(x))) / x / 2.0
    for k in range(1, n_euler + m_euler + 1):
        sk = x + 1j * np.pi * k / t
        terms[k] = ((-1) ** k) * float(np.real(lst_fun(sk) / sk))
    partial = np.cumsum(terms)
    wts = np.array([comb(m_euler, j) / 2 ** m_euler for j in range(m_euler + 1)])
    F = u * float(np.sum(wts * partial[n_euler:n_euler + m_euler + 1]))
    return min(max(F, 0.0), 1.0)


def qsys_mmapgk1(MMAP: Sequence, svc: Sequence, w_points: Optional[Sequence[float]] = None,
                 num_w_moms: int = 3, tol: float = 1e-12, iter_max: int = 10000) -> MmapGk1Result:
    """
    Analyze an MMAP[K]/G[K]/1 FCFS queue.

    Args:
        MMAP: LINE convention [D0, D1, D^(1), ..., D^(K)] with D1 = sum_k D^(k)
        svc: K service laws, each a LINE Distribution, a [D0, D1] phase-type
            pair, or a dict {'lst': handle, 'moments': [E[S], E[S^2], ...]} with
            at least num_w_moms+1 moments
        w_points: times at which to evaluate the per-type waiting time CDF
        num_w_moms: how many per-type waiting time moments to return

    Returns:
        MmapGk1Result
    """
    K = len(MMAP) - 2
    if K < 1:
        raise ValueError("The MMAP must carry at least one marked arrival block")
    if len(svc) != K:
        raise ValueError("One service law per marked type is required "
                         "(%d given, %d types)" % (len(svc), K))
    D0 = np.atleast_2d(np.asarray(MMAP[0], dtype=float))
    ma = D0.shape[0]
    Dk = [np.atleast_2d(np.asarray(MMAP[k + 2], dtype=float)) for k in range(K)]
    Dsum = D0.copy()
    for k in range(K):
        Dsum = Dsum + Dk[k]

    theta = _stat_left_null(Dsum)
    lambdas = np.array([float(theta @ Dk[k] @ np.ones(ma)) for k in range(K)])
    mean_s = np.array([_svc_mean(svc[k]) for k in range(K)])
    rho = float(np.sum(lambdas * mean_s))
    if rho >= 1.0:
        raise ValueError("The load %g of the system is not below one" % rho)

    U = D0.copy()
    for _ in range(iter_max):
        Unew = D0.copy()
        for k in range(K):
            Unew = Unew + Dk[k] @ _matrix_lst(svc[k], U)
        done = np.max(np.abs(Unew - U)) <= tol
        U = Unew
        if done:
            break
    # U e = 0 EXACTLY. It is a property of the fixed point, not of the iterate:
    # the iteration converges linearly, so the row sums still carry O(1e-9) at
    # the tolerance above, and that residue moves the eigenvalue that belongs at
    # the origin off it, where it would be mistaken for a real condition.
    U = U - np.diag(U.sum(axis=1))

    # The idle vector. U is a proper generator: its off-diagonals are nonnegative
    # (D0's are, Dk >= 0 and Fk(U) >= 0) and its rows sum to zero, and it governs
    # the arrival phase at the epochs the level drops by one. Its stationary
    # vector carries the idle mass,
    #
    #     v0 = (1 - rho) pi_U,     pi_U U = 0,  pi_U e = 1,
    #
    # which is what the analyticity conditions of (*) deliver: imposing
    # v0 r_i = 0 at the right null vector of the bracket for each of the ma-1
    # roots s = -u off the origin, plus v0 e = 1 - rho, reproduces this vector to
    # 2.5e-13 on an Erlang case. The stationary route needs no complex
    # eigenvector and no rule for telling the structural root at the origin from
    # a genuine one.
    v0 = (1 - rho) * _stat_left_null(U)

    def _bracket(s: complex) -> np.ndarray:
        M = s * np.eye(ma) + D0.astype(complex)
        for k in range(K):
            M = M + Dk[k] * _scalar_lst(svc[k], s)
        return M

    def wait_lst(s: complex) -> np.ndarray:
        if abs(s) < 1e-14:
            return np.ones(K)
        f = np.linalg.solve(_bracket(s).T, (s * v0).astype(complex))
        return np.array([complex(f @ Dk[k] @ np.ones(ma)) / lambdas[k] for k in range(K)])

    wait_moments = _moments(D0, Dk, svc, theta, v0, lambdas, num_w_moms)
    mean_wt = wait_moments[:, 0]
    mean_st = mean_wt + mean_s
    mean_ql = float(np.sum(lambdas * mean_st))

    wait_cdf = None
    pts = None
    if w_points is not None and len(np.atleast_1d(w_points)) > 0:
        pts = np.atleast_1d(np.asarray(w_points, dtype=float))
        wait_cdf = np.zeros((K, pts.size))
        for k in range(K):
            for it, t in enumerate(pts):
                wait_cdf[k, it] = _euler_invert(lambda s, kk=k: wait_lst(s)[kk], float(t))

    return MmapGk1Result(
        lambdas=lambdas,
        arrivalRate=float(np.sum(lambdas)),
        utilization=rho,
        idleVector=v0,
        waitLST=wait_lst,
        waitMoments=wait_moments,
        meanWaitingTime=mean_wt,
        meanSojournTime=mean_st,
        meanQueueLength=mean_ql,
        waitCDF=wait_cdf,
        waitPoints=pts,
        analyzer="LINE:MMAP[%d]/G[%d]/1" % (K, K),
    )


def _moments(D0, Dk, svc, theta, v0, lambdas, num_w_moms):
    """
    Derivatives of f(s) M(s) = s v0 at s = 0. M_0 is SINGULAR with right null
    vector e, so each order fixes f_j only up to a multiple of theta; that
    multiple is what the NEXT order's solvability condition supplies. At j = 0
    the same condition reads theta M_1 e = v0 e, i.e. 1 - rho = 1 - rho, which
    is the identity that validates the setup.
    """
    ma = D0.shape[0]
    K = len(Dk)
    Mder = []
    for j in range(num_w_moms + 2):
        if j == 0:
            Mj = D0.copy()
            for k in range(K):
                Mj = Mj + Dk[k]
        else:
            Mj = np.zeros((ma, ma))
            for k in range(K):
                Mj = Mj + Dk[k] * (((-1) ** j) * _raw_moment(svc[k], j))
            if j == 1:
                Mj = Mj + np.eye(ma)
        Mder.append(Mj)
    e = np.ones(ma)
    denom = float(theta @ Mder[1] @ e)
    fder = [theta]
    Abase = np.hstack([Mder[0], e.reshape(-1, 1)])
    for j in range(1, num_w_moms + 1):
        rhs = np.zeros(ma)
        if j == 1:
            rhs = rhs + v0
        for i in range(j):
            rhs = rhs - comb(j, i) * (fder[i] @ Mder[j - i])
        fp = np.linalg.lstsq(Abase.T, np.concatenate([rhs, [0.0]]), rcond=None)[0]
        acc = 0.0
        for i in range(j):
            acc += comb(j + 1, i) * float(fder[i] @ Mder[j + 1 - i] @ e)
        cfree = (-acc / (j + 1) - float(fp @ Mder[1] @ e)) / denom
        fder.append(fp + cfree * theta)
    moms = np.zeros((K, num_w_moms))
    for k in range(K):
        for j in range(1, num_w_moms + 1):
            moms[k, j - 1] = ((-1) ** j) * float(fder[j] @ Dk[k] @ e) / lambdas[k]
    return moms
