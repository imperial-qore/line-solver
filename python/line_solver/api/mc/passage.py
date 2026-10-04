"""
First passage times in Markov and semi-Markov chains.

Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions in
Large Markov Chains", 2002. Eqs. 1-3 for the Markov case, Eqs. 4-8 for the
semi-Markov one.

Twin of the MATLAB matlab/src/api/mc/ctmc_passage_*.m and smp_passage_*.m.
TARGET indices are 0-based here and 1-based there, as elsewhere between the two
codebases.
"""

from math import comb, factorial
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.linalg import expm

from ..lti import laplace_invert_cdf, laplace_invert_pdf


def _as_target(target, n: int) -> np.ndarray:
    t = np.unique(np.atleast_1d(np.asarray(target, dtype=int)).ravel())
    if t.size == 0:
        raise ValueError("The target state set is empty: a first passage time "
                         "into no state is undefined.")
    if t.min() < 0 or t.max() >= n:
        raise ValueError("A target state index is outside the state space.")
    return t


def ctmc_passage_ph(Q: np.ndarray, pi0: Optional[np.ndarray],
                    target) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                     np.ndarray, float]:
    """
    Phase-type representation (alpha, S, s0, keep, atom) of the first passage
    time from pi0 into the target state set.

    With A the complement of the target, S = Q(A,A), s0 = -S*1 and alpha =
    pi0(A), so L(s) = alpha (sI-S)^{-1} s0 + atom and F(t) = 1 - alpha exp(St)1.
    That is Eqs. 1-2 written as a phase-type law rather than as n scalar
    equations.

    ALPHA IS DELIBERATELY NOT NORMALIZED. Its mass is 1 - atom; the missing mass
    is the ATOM AT ZERO carried by initial states already inside the target. A
    caller that normalizes alpha and forgets the atom reports F(0) = 0 for a
    passage that has already completed with probability atom.
    """
    Q = np.asarray(Q, dtype=float)
    n = Q.shape[0]
    if Q.shape[1] != n:
        raise ValueError("The generator must be square.")
    tgt = _as_target(target, n)
    scale = max(1.0, float(np.max(np.abs(Q))) if Q.size else 1.0)
    if np.max(np.abs(Q.sum(axis=1))) > 1e-8 * scale:
        raise ValueError("Q is not an infinitesimal generator: its rows do not "
                         "sum to zero. Pass it through ctmc_makeinfgen first.")

    is_target = np.zeros(n, dtype=bool)
    is_target[tgt] = True
    keep = np.flatnonzero(~is_target)

    S = Q[np.ix_(keep, keep)]
    s0 = -S.dot(np.ones(len(keep)))

    if pi0 is None:
        from .ctmc import ctmc_solve
        p = np.asarray(ctmc_solve(Q), dtype=float).ravel()
        mass = p[keep].sum()
        if mass <= 0:
            raise ValueError("The stationary law puts no mass outside the target "
                             "set, so there is no passage to time.")
        alpha = p[keep] / mass
        atom = 0.0
    else:
        pi0 = np.asarray(pi0, dtype=float).ravel()
        if pi0.size != n:
            raise ValueError("pi0 must be a distribution over the state space, "
                             "one entry per state.")
        alpha = pi0[keep]
        atom = float(pi0[tgt].sum())
    return alpha, S, s0, keep, atom


def ctmc_passage_lst(Q: np.ndarray, pi0: Optional[np.ndarray], target,
                     s) -> np.ndarray:
    """
    Laplace-Stieltjes transform of the first passage time, L(s) = alpha
    (sI-S)^{-1} s0 + atom, at the (possibly complex) points s. Eqs. 1-2: one
    linear system per value of s.

    ONE SOLVE PER s, NOT PER (s,t) PAIR. The saving over a dense matrix
    exponential is that the solves are sparse, so this route reaches chains a
    dense expm cannot hold. It is NOT a saving in the number of time points:
    every Abate-Whitt inverter places its nodes at s = beta/t, so a grid of T
    points costs T*|beta| solves. On a small chain ctmc_passage_time's default
    'expm' route is faster.
    """
    alpha, S, s0, _, atom = ctmc_passage_ph(Q, pi0, target)
    sv = np.atleast_1d(np.asarray(s))
    I = np.eye(S.shape[0])
    out = np.zeros(sv.shape, dtype=complex)
    flat = sv.ravel()
    res = np.zeros(flat.size, dtype=complex)
    for i, si in enumerate(flat):
        res[i] = alpha.dot(np.linalg.solve(si * I - S, s0)) + atom
    out = res.reshape(sv.shape)
    if np.all(np.abs(out.imag) < 1e-14):
        out = out.real
    return out


def _reaches_target(Q: np.ndarray, keep: np.ndarray, n: int,
                    tgt: np.ndarray) -> np.ndarray:
    """Backward reachability closure over the transition graph."""
    adj = Q != 0
    np.fill_diagonal(adj, False)
    seen = np.zeros(n, dtype=bool)
    seen[tgt] = True
    frontier = list(tgt)
    while frontier:
        pred = np.flatnonzero(adj[:, frontier].any(axis=1))
        pred = pred[~seen[pred]]
        seen[pred] = True
        frontier = list(pred)
    return seen[keep]


def ctmc_passage_moments(Q: np.ndarray, pi0: Optional[np.ndarray], target,
                         nmax: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    """
    Moments of order 1..nmax of the first passage time into the target set.

    Returns (mall, m): mall is (nstates, nmax), row i for a passage started in
    state i, zero on target states and inf where the target cannot be reached;
    m is the pi0-weighted moment vector.

    This is Eq. 3, -q_ii M_i(n) = sum_{k not in B} q_ik M_k(n) + n M_i(n-1),
    i.e. (-S) M(n) = n M(n-1) with M(0) = 1: nmax linear solves and no
    transform inversion at all. The equivalent closed form n! alpha (-S)^{-n} 1
    is NOT how it is evaluated -- forming the inverse of the sub-generator
    destroys the sparsity the recursion preserves.
    """
    alpha, S, _, keep, atom = ctmc_passage_ph(Q, pi0, target)
    Q = np.asarray(Q, dtype=float)
    n = Q.shape[0]
    tgt = _as_target(target, n)
    nA = len(keep)
    mall = np.zeros((n, nmax))
    if nA == 0:
        return mall, np.zeros(nmax)

    # A state that cannot reach the target has an infinite passage time; the
    # sub-generator is singular on that block and a least-squares solve would
    # return a finite number instead of saying so.
    reach = _reaches_target(Q, keep, n, tgt)

    A = -S
    x = np.ones(nA)
    X = np.zeros((nA, nmax))
    for k in range(1, nmax + 1):
        if reach.all():
            x = np.linalg.solve(A, k * x)
        else:
            x = np.linalg.lstsq(A, k * x, rcond=None)[0]
        X[:, k - 1] = x
    X[~reach, :] = np.inf

    mall[keep, :] = X
    if np.any((~reach) & (alpha > 0)):
        m = np.full(nmax, np.inf)
    else:
        sel = reach & (alpha != 0)
        m = np.array([alpha[sel].dot(X[sel, k]) for k in range(nmax)])
    if atom >= 1:
        m = np.zeros(nmax)
    return mall, m


def ctmc_hitting_time(Q: np.ndarray, target_states) -> np.ndarray:
    """
    Mean time to reach any state in target_states from each state of a CTMC.

    Continuous-time twin of dtmc_hitting_time and the first-moment special case
    of ctmc_passage_moments: (-S) h = 1 on the non-target block, where
    dtmc_hitting_time solves (I - P_NT) h = 1. Unreachable states give inf.
    """
    Q = np.asarray(Q, dtype=float)
    n = Q.shape[0]
    # mall does not depend on the initial law, so a uniform one is passed
    # rather than paying for the stationary solve that pi0 = None would trigger.
    mall, _ = ctmc_passage_moments(Q, np.ones(n) / n, target_states, 1)
    return mall[:, 0]


def ctmc_passage_time(Q: np.ndarray, pi0: Optional[np.ndarray], target, tset,
                      method: str = 'expm', lti_method: str = 'euler'
                      ) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    CDF and density of the first passage time into the target state set:
    F(t) = 1 - alpha exp(St) 1 and f(t) = alpha exp(St) s0.

    method='expm' (default) is exact and reuses one matrix exponential along a
    uniform grid. method='lt' inverts the transform of Eqs. 1-2 through api/lti;
    it exists for chains whose non-target block is too large for a dense
    exp(St), not because it needs fewer time points (see ctmc_passage_lst). On a
    small chain 'expm' is both faster and more accurate, hence the default.

    The returned dict carries 'atom', the mass of pi0 already inside the target,
    which is F(0).
    """
    alpha, S, s0, keep, atom = ctmc_passage_ph(Q, pi0, target)
    tv = np.atleast_1d(np.asarray(tset, dtype=float)).ravel()
    F = np.zeros(tv.size)
    f = np.zeros(tv.size)
    out = {'atom': atom, 'alpha': alpha, 'S': S, 's0': s0, 'keep': keep,
           'method': method}

    if method == 'expm':
        e = np.ones(S.shape[0])
        dt = np.diff(tv)
        uniform = (tv.size > 2 and dt[0] > 0
                   and np.all(np.abs(dt - dt[0]) < 1e-12 * max(1.0, abs(dt[0]))))
        if uniform:
            # One exponential, then propagate: recomputing expm(S*t) at every
            # grid point is the same answer at a cost linear in the grid.
            E = expm(S * dt[0])
            v = alpha.dot(expm(S * tv[0]))
            for i in range(tv.size):
                if i > 0:
                    v = v.dot(E)
                F[i] = 1.0 - v.dot(e)
                f[i] = v.dot(s0)
        else:
            for i, t in enumerate(tv):
                if t < 0:
                    continue
                v = alpha.dot(expm(S * t))
                F[i] = 1.0 - v.dot(e)
                f[i] = v.dot(s0)
    elif method == 'lt':
        I = np.eye(S.shape[0])

        def Lfun(sv):
            return alpha.dot(np.linalg.solve(sv * I - S, s0)) + atom

        F = laplace_invert_cdf(Lfun, tv, lti_method)
        f = laplace_invert_pdf(lambda sv: Lfun(sv) - atom, tv, lti_method)
    else:
        raise ValueError("Unknown passage-time method: %s. Supported: expm, lt."
                         % method)

    return np.clip(F, 0.0, 1.0), np.maximum(f, 0.0), out


# --------------------------------------------------------------------------
# Semi-Markov chains (Sec. 3)
# --------------------------------------------------------------------------

def smp_passage_moments(P: np.ndarray, hmom, pi0: Optional[np.ndarray], target,
                        nmax: int = 1) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """
    Moments of order 1..nmax of the first passage time into the target state set
    for a semi-Markov chain with embedded transition matrix P.

    hmom selects which of the paper's two recursions runs, and they are NOT the
    same computation:

      (nstates, nmax) array   m_i(r), the holding time in i depends only on i.
                              Eq. 7 with the u_i(r) recurrence of Eq. 8,

                                  u_i(r) = -sum_{j=1..r} C(r,j) m_i(j) u_i(r-j),
                                  u_i(0) = 1,

                              which are the derivatives at the origin of
                              1/h*_i(s). Cheaper: no per-pair moments.

      (nstates, nstates) list-of-lists   hmom[i][k] = [m_ik(1) ... m_ik(nmax)],
                              the r-th moment of the holding time in i WHEN THE
                              NEXT STATE IS k. Eq. 6, the full Markov-renewal
                              kernel. m_ik(0) = P[i,k] is implied and must not
                              be supplied.

    Unlike the Markov case the n-th moment needs every moment from 1 to n, so
    nmax cannot be raised for free.
    """
    P = np.asarray(P, dtype=float)
    n = P.shape[0]
    if P.shape[1] != n:
        raise ValueError("The embedded transition matrix must be square.")
    if np.max(np.abs(P.sum(axis=1) - 1.0)) > 1e-8:
        raise ValueError("The embedded transition matrix rows must sum to one.")
    tgt = _as_target(target, n)

    is_target = np.zeros(n, dtype=bool)
    is_target[tgt] = True
    A = np.flatnonzero(~is_target)
    nA = len(A)
    mall = np.zeros((n, nmax))
    if nA == 0:
        return mall, (np.zeros(nmax) if pi0 is not None else None)

    PAA = P[np.ix_(A, A)]
    I = np.eye(nA)
    Msub = np.zeros((nA, nmax))

    kernel = not isinstance(hmom, np.ndarray) or hmom.dtype == object
    if kernel:
        # --- Eq. 6, the full kernel ---------------------------------------
        def mik(i, k):
            v = hmom[i][k]
            return np.zeros(nmax) if v is None or len(v) == 0 else np.asarray(v, float)

        for q in range(1, nmax + 1):
            b = np.zeros(nA)
            for a, i in enumerate(A):
                acc = 0.0
                for c, k in enumerate(A):
                    mv = mik(i, k)
                    for r in range(1, q + 1):
                        term = comb(q, r) * mv[r - 1]
                        acc += term * (Msub[c, q - r - 1] if r < q else 1.0)
                for k in tgt:
                    acc += mik(i, k)[q - 1]
                b[a] = acc
            Msub[:, q - 1] = np.linalg.solve(I - PAA, b)
    else:
        # --- Eq. 7 with the Eq. 8 recurrence -------------------------------
        hm = np.asarray(hmom, dtype=float)
        if hm.shape[0] != n:
            raise ValueError("An array hmom must carry one row per state.")
        if hm.shape[1] < nmax:
            raise ValueError("hmom must carry at least nmax holding-time "
                             "moments per state.")
        u = _smp_u(hm[A, :], nmax)
        for q in range(1, nmax + 1):
            b = np.zeros(nA)
            for r in range(1, q + 1):
                base = Msub[:, q - r - 1] if r < q else np.ones(nA)
                b -= comb(q, r) * u[:, r - 1] * base
            Msub[:, q - 1] = np.linalg.solve(I - PAA, b)

    mall[A, :] = Msub
    if pi0 is None:
        return mall, None
    pi0 = np.asarray(pi0, dtype=float).ravel()
    if pi0.size != n:
        raise ValueError("pi0 must be a distribution over the state space, one "
                         "entry per state.")
    return mall, pi0.dot(mall)


def _smp_u(mrows: np.ndarray, nmax: int) -> np.ndarray:
    """
    Eq. 8: u_i(r) = -sum_{j=1..r} C(r,j) m_i(j) u_i(r-j), u_i(0) = 1. These are
    the derivatives at the origin of 1/h*_i(s), obtained from the moments of
    h*_i alone by differentiating h*(s) y(s) = 1 r times.
    """
    nA = mrows.shape[0]
    u = np.zeros((nA, nmax))
    u0 = np.ones(nA)
    for r in range(1, nmax + 1):
        acc = np.zeros(nA)
        for j in range(1, r + 1):
            base = u0 if r - j == 0 else u[:, r - j - 1]
            acc += comb(r, j) * mrows[:, j - 1] * base
        u[:, r - 1] = -acc
    return u


def smp_passage_lst(P: np.ndarray, hlst, pi0: Optional[np.ndarray], target,
                    s) -> np.ndarray:
    """
    Laplace-Stieltjes transform of the semi-Markov first passage time, Eqs. 4-5:

        L_i(s) = sum_{k not in B} r*_ik(s) L_k(s) + sum_{k in B} r*_ik(s)

    so (I - R*_AA(s)) L_A(s) = R*_AB(s) 1, one linear system per value of s.

    hlst is either a length-nstates sequence of handles h*_i(s), in which case
    r*_ik(s) = P[i,k] h*_i(s) and the complex numbers stay on the DIAGONAL of
    the system (Eq. 5), or an (nstates, nstates) nested sequence of handles
    r*_ik(s) for the full Markov-renewal kernel, the harder case the paper
    flags.

    Distribution objects supply their own transform: Markovian.evalLST gives the
    closed form pie (sI-D0)^{-1} (-D0) e for the phase-type family, so
    hlst[i] = dist.evalLST is the intended way to build these.
    """
    P = np.asarray(P, dtype=float)
    n = P.shape[0]
    tgt = _as_target(target, n)
    is_target = np.zeros(n, dtype=bool)
    is_target[tgt] = True
    A = np.flatnonzero(~is_target)
    nA = len(A)

    if pi0 is None:
        pi0 = np.ones(n) / n
    pi0 = np.asarray(pi0, dtype=float).ravel()
    if pi0.size != n:
        raise ValueError("pi0 must be a distribution over the state space, one "
                         "entry per state.")
    atom = float(pi0[tgt].sum())

    per_state = not hasattr(hlst[0], '__len__')
    sv = np.atleast_1d(np.asarray(s))
    flat = sv.ravel()
    res = np.zeros(flat.size, dtype=complex)
    I = np.eye(nA)
    for idx, si in enumerate(flat):
        if per_state:
            h = np.array([hlst[i](si) for i in A], dtype=complex)
            RAA = h[:, None] * P[np.ix_(A, A)]
            RAB = h * P[np.ix_(A, tgt)].sum(axis=1)
        else:
            Rs = np.zeros((n, n), dtype=complex)
            for i in range(n):
                for k in range(n):
                    fn = hlst[i][k]
                    if fn is not None:
                        Rs[i, k] = fn(si)
            RAA = Rs[np.ix_(A, A)]
            RAB = Rs[np.ix_(A, tgt)].sum(axis=1)
        res[idx] = pi0[A].dot(np.linalg.solve(I - RAA, RAB)) + atom
    out = res.reshape(sv.shape)
    if np.all(np.abs(out.imag) < 1e-14):
        out = out.real
    return out


def smp_passage_time(P: np.ndarray, hlst, pi0: Optional[np.ndarray], target,
                     tset, lti_method: str = 'euler'
                     ) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    CDF and density of the semi-Markov first passage time, by inverting
    smp_passage_lst through api/lti.

    There is no matrix-exponential route here: a semi-Markov chain has no
    generator to exponentiate, which is exactly the case uniformization does not
    reach and the transform does.

    lti_method defaults to 'euler' RATHER THAN 'weeks'. Semi-Markov passage
    densities are the case Sec. 4.2 singles out as slow-converging for a
    Laguerre series: a kernel with a deterministic or discontinuous holding time
    gives a density whose derivatives jump, and laplace_weeks_scaling then
    refuses by name rather than returning noise.
    """
    P = np.asarray(P, dtype=float)
    n = P.shape[0]
    tgt = _as_target(target, n)
    if pi0 is None:
        pi0 = np.ones(n) / n
    pi0 = np.asarray(pi0, dtype=float).ravel()
    atom = float(pi0[tgt].sum())

    def Lfun(sv):
        return complex(np.atleast_1d(smp_passage_lst(P, hlst, pi0, tgt, sv))[0])

    tv = np.atleast_1d(np.asarray(tset, dtype=float)).ravel()
    F = laplace_invert_cdf(Lfun, tv, lti_method)
    f = laplace_invert_pdf(lambda sv: Lfun(sv) - atom, tv, lti_method)
    return F, f, {'atom': atom, 'lti_method': lti_method}
