"""
Continuous-Time Markov Chain (CTMC) analysis algorithms.

Native Python implementations for CTMC steady-state, transient analysis,
and related methods. Uses scipy for linear algebra and numerical integration.

Key algorithms:
    ctmc_solve: Steady-state distribution
    ctmc_makeinfgen: Construct valid infinitesimal generator
    ctmc_transient: Transient probabilities via matrix exponential
    ctmc_uniformization: Uniformization for transient analysis
    ctmc_stochcomp: Stochastic complementation
"""

import numpy as np
from numpy.linalg import LinAlgError
from scipy import linalg
import scipy.sparse as sp
from scipy.sparse import csc_matrix, issparse, diags
from scipy.sparse.linalg import spsolve, gmres, splu
from scipy.sparse.csgraph import connected_components
from scipy.integrate import solve_ivp
from typing import Dict, Any, Optional, List, Tuple, Union
from dataclasses import dataclass
import json
import os
import tempfile
import time
import platform


# see _kb/03-api-layer.md for rationale
GMRES_MIN_STATES = 6000


def _raise_no_recurrent_state():
    """Raised when the active-set elimination in ctmc_solve leaves no state.

    The elimination drops every state whose row is all-zero, which is precisely
    an ABSORBING state; ctmc_makeinfgen then re-zeroes the diagonal of the
    survivors that only fed it, so the elimination cascades until nothing is
    left. Returning a uniform vector here does NOT satisfy p*Q=0 (it is not a
    stationary distribution, just a shape of the right size), and a caller
    cannot tell it apart from a real answer: a generator missing all its
    arrivals reads back as a plausible mean of cutoff/2. Fail instead. A
    genuinely absorbing chain has no unique stationary distribution without an
    initial vector, so it belongs in ctmc_solve_reducible.
    """
    raise ValueError(
        "The infinitesimal generator has no recurrent state: every state was "
        "eliminated as absorbing. This generator admits no unique stationary "
        "distribution. It usually means the generator is malformed -- e.g. a "
        "state with no outgoing transitions that absorbs the whole chain, as "
        "happens when a class of transitions was dropped while building it. "
        "Use ctmc_solve_reducible for a genuinely absorbing chain.")


def issym(Q):
    """Check if matrix contains sympy symbolic expressions."""
    try:
        import sympy
        if isinstance(Q, sympy.MatrixBase):
            return True
        if isinstance(Q, np.ndarray) and Q.dtype == object:
            return any(isinstance(x, sympy.Basic) for x in Q.flat)
    except ImportError:
        pass
    return False


def ctmc_makeinfgen(Q):
    """
    Convert a matrix into a valid infinitesimal generator for a CTMC.

    An infinitesimal generator has:
    - Row sums equal to zero
    - Non-positive diagonal elements
    - Non-negative off-diagonal elements

    Args:
        Q: Candidate infinitesimal generator matrix

    Returns:
        Valid infinitesimal generator matrix with corrected diagonal
    """
    if issym(Q):
        import sympy
        if isinstance(Q, sympy.MatrixBase):
            M = Q.copy()
        else:
            M = sympy.Matrix(Q.tolist())
        n = M.rows
        for i in range(n):
            M[i, i] = 0
        for i in range(n):
            M[i, i] = -sum(M[i, j] for j in range(n))
        return M

    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]

    # Extract off-diagonal elements
    result = Q.copy()
    np.fill_diagonal(result, 0.0)

    # Set diagonal to negative row sum (ensures row sums = 0)
    row_sums = result.sum(axis=1)
    np.fill_diagonal(result, -row_sums)

    return result


def _find_weakly_connected_components(Q: np.ndarray) -> List[List[int]]:
    """
    Find weakly connected components in the CTMC graph.

    Args:
        Q: Generator matrix

    Returns:
        List of lists, each containing state indices for a component
    """
    n = Q.shape[0]

    # Build adjacency matrix: B = |Q + Q'| > 0
    B = np.abs(Q + Q.T) > 0

    # Find connected components
    n_components, labels = connected_components(
        csc_matrix(B), directed=False, return_labels=True
    )

    if n_components == 1:
        return [list(range(n))]

    # Group states by component
    components = [[] for _ in range(n_components)]
    for i, label in enumerate(labels):
        components[label].append(i)

    return components


def ctmc_solve(Q: np.ndarray, method: Optional[str] = None) -> np.ndarray:
    """
    Solve for steady-state probabilities of a CTMC.

    Computes the stationary distribution π by solving πQ = 0
    with normalization constraint Σπ = 1.

    Handles reducible CTMCs by decomposing into strongly connected
    components and solving each separately.

    Args:
        Q: Infinitesimal generator matrix (row sums should be zero)
        method: 'gmres' or 'direct' to force a solution method; None or
            'default' selects by size, GMRES above GMRES_MIN_STATES states

    Returns:
        Steady-state probability distribution (1D array)
    """
    if issym(Q):
        import sympy
        if isinstance(Q, sympy.MatrixBase):
            Qs = Q.copy()
        else:
            Qs = sympy.Matrix(Q.tolist())
        n = Qs.rows
        if n == 1:
            return sympy.Matrix([[1]])
        Qs = ctmc_makeinfgen(Qs)

        # see _kb/03-api-layer.md for rationale
        from ..sym import prefers_sage, resolve as _resolve_sym
        if prefers_sage():
            engine = _resolve_sym()
            if engine is not None:
                local = dict((str(s), s) for s in Qs.free_symbols)
                symbols = sorted(local)
                r = engine.solve_ctmc(
                    [[str(Qs[i, j]) for j in range(n)] for i in range(n)], symbols)
                # see _kb/03-api-layer.md for rationale
                return sympy.Matrix([[sympy.sympify(p.replace("^", "**"), locals=local)
                                      for p in r["pi"]]])

        def _sym_to_numeric(M):
            # Structure of the generator with all symbols replaced by 1.0
            # (MATLAB: double(subs(Q, symvar(Q), ones(...))))
            return np.array(
                M.subs({s: 1.0 for s in M.free_symbols}).evalf().tolist(),
                dtype=np.float64,
            )

        # Reducible generator: solve each weakly connected component
        # recursively and renormalize (MATLAB ctmc_solve lines 66-82)
        components = _find_weakly_connected_components(_sym_to_numeric(Qs))
        if len(components) > 1:
            pi = sympy.zeros(1, n)
            for comp in components:
                Qc = ctmc_makeinfgen(Qs[comp, comp])
                pc = ctmc_solve(Qc)
                for k, s in enumerate(comp):
                    pi[s] = pc[k]
            return pi / sum(pi)

        # Iteratively trim states with zero row or column sums
        nnzel = np.arange(n)
        Qnnz = Qs
        prev_size = -1
        while Qnnz.rows != prev_size:
            prev_size = Qnnz.rows
            Qabs = np.abs(_sym_to_numeric(Qnnz))
            active = np.where((Qabs.sum(axis=1) > 0) & (Qabs.sum(axis=0) > 0))[0]
            if len(active) == 0:
                _raise_no_recurrent_state()
            if len(active) < Qnnz.rows:
                nnzel = nnzel[active]
                Qnnz = ctmc_makeinfgen(Qnnz[list(active), list(active)])

        # Solve pi * Qmod = e_m where Qmod is Qnnz with its last column
        # replaced by ones (normalization), as in the numeric branch
        m = Qnnz.rows
        A = Qnnz.copy()
        for i in range(m):
            A[i, m - 1] = 1
        b = sympy.zeros(m, 1)
        b[m - 1] = 1
        try:
            x = A.T.solve(b)
        except Exception:
            # Singular system: re-check reducibility of the trimmed generator
            components = _find_weakly_connected_components(_sym_to_numeric(Qnnz))
            if len(components) > 1:
                pi = sympy.zeros(1, n)
                for comp in components:
                    Qc = ctmc_makeinfgen(Qnnz[comp, comp])
                    pc = ctmc_solve(Qc)
                    for k, s in enumerate(comp):
                        pi[nnzel[s]] = pc[k]
                return pi / sum(pi)
            # NaN signals failure to the caller, as in the numeric branch
            return sympy.nan * sympy.ones(1, n)
        pi = sympy.zeros(1, n)
        for k, s in enumerate(nnzel):
            pi[s] = x[k]
        return pi

    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]

    # Trivial case
    if n == 1:
        return np.array([1.0])

    # Ensure valid generator
    Q = ctmc_makeinfgen(Q)

    # No transitions at all: every distribution satisfies p*Q=0, so the
    # stationary distribution is not unique and uniform is as good as any.
    if np.all(Q == 0):
        return np.ones(n) / n

    # see _kb/03-api-layer.md for rationale
    nnzel = np.arange(n)
    Qnnz = Q.copy()
    Qnnz_prev_size = -1

    while len(Qnnz) != Qnnz_prev_size:
        Qnnz_prev_size = len(Qnnz)
        row_sums = np.abs(Qnnz).sum(axis=1)
        col_sums = np.abs(Qnnz).sum(axis=0)
        active = (row_sums > 0) & (col_sums > 0)

        if not np.any(active):
            _raise_no_recurrent_state()

        active_local = np.where(active)[0]
        nnzel = nnzel[active_local]
        Qnnz = Qnnz[np.ix_(active_local, active_local)]
        Qnnz = ctmc_makeinfgen(Qnnz)

    if Qnnz.size == 0:
        _raise_no_recurrent_state()

    active_idx = nnzel

    # see _kb/03-api-layer.md for rationale
    n_active = len(active_idx)

    # Check for multiple weakly connected components first (matches JAR lines 46-73)
    sub_components = _find_weakly_connected_components(Qnnz)
    if len(sub_components) > 1:
        pi_active = np.zeros(n_active)
        for comp in sub_components:
            idx = np.array(comp)
            Qc = Qnnz[np.ix_(idx, idx)]
            Qc = ctmc_makeinfgen(Qc)
            pc = ctmc_solve(Qc, method)
            for i, state in enumerate(comp):
                pi_active[state] = pc[i]
        if pi_active.sum() > 0:
            pi_active /= pi_active.sum()
        else:
            pi_active = np.ones(n_active) / n_active

        p = np.zeros(n)
        p[active_idx] = pi_active
        if p.sum() > 0:
            p /= p.sum()
        else:
            p = np.ones(n) / n
        return p

    # Single component: solve Q'x = b with last column replaced by 1s
    # This matches JAR ctmc_solve lines 147-157
    Qsol = Qnnz.copy()
    Qsol[:, -1] = 1.0  # Replace last column with 1s
    b = np.zeros(n_active)
    b[-1] = 1.0  # Normalization constraint

    pi_active = np.full(n_active, np.nan)
    QsolT = Qsol.T

    # see _kb/03-api-layer.md for rationale
    method = (method or 'default').lower()
    if method == 'gmres' or (method != 'direct' and n_active > GMRES_MIN_STATES):
        from .gmres import ctmc_gmres
        x_gmres, gflag, _, _ = ctmc_gmres(sp.csc_matrix(QsolT), b)
        if gflag == 0:
            p = np.zeros(n)
            p[active_idx] = x_gmres
            total = p.sum()
            return p / total if total > 0 else np.ones(n) / n

    # see _kb/03-api-layer.md for rationale
    rcond_threshold = 1e-10
    fast_pi = None
    try:
        lu, piv = linalg.lu_factor(QsolT)
        anorm = linalg.norm(QsolT, 1)
        rcond, _info = linalg.lapack.dgecon(lu, anorm)
        if rcond > rcond_threshold:
            candidate = linalg.lu_solve((lu, piv), b)
            if np.all(np.isfinite(candidate)):
                fast_pi = candidate
    except LinAlgError:
        fast_pi = None

    if fast_pi is not None:
        pi_active = fast_pi
    else:
        qsol_rank = np.linalg.matrix_rank(QsolT)
        if qsol_rank < n_active:
            # Matrix is singular - return NaN to trigger fallback in caller
            # (matches MATLAB ctmc_solve behavior for reducible generators)
            p = np.zeros(n)
            p[active_idx] = pi_active  # NaN values
            return p

        try:
            pi_active = np.linalg.solve(QsolT, b)
        except np.linalg.LinAlgError:
            pass  # Leave as NaN - will trigger fallback in caller

    # If solution has NaN, try weakly connected component decomposition
    # (matches JAR lines 159-201)
    if np.any(np.isnan(pi_active)):
        # Already checked components above and found only 1, so just return NaN
        p = np.zeros(n)
        p[active_idx] = pi_active
        return p

    # Map back to full state space
    p = np.zeros(n)
    p[active_idx] = pi_active

    # see _kb/03-api-layer.md for rationale
    if p.sum() > 0:
        p /= p.sum()
    else:
        p = np.ones(n) / n

    return p


def ctmc_sens(Q: np.ndarray, dQ: np.ndarray,
              pi: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Sensitivity of the steady-state distribution of a CTMC to a scalar
    parameter theta, given the generator Q, its derivative dQ = dQ/dtheta,
    and the steady-state vector pi.

    Differentiating the balance equations pi*Q = 0 and pi*e = 1 with respect
    to theta gives the linear system

      (dpi/dtheta) * Q = -pi * (dQ/dtheta),   sum_i dpi_i/dtheta = 0,

    i.e. Trivedi and Bobbio (2017), Eq. (9.81). The system has the same
    coefficient matrix as the steady-state solve itself, so obtaining a
    sensitivity costs one extra solve against a matrix that is already
    assembled. The normalization replaces one row of the singular Q',
    exactly as in the steady-state solve.

    Args:
        Q: Generator matrix (n x n)
        dQ: Derivative of the generator with respect to theta (n x n)
        pi: Steady-state distribution (1 x n); computed if omitted

    Returns:
        Derivative of the steady-state distribution (1D array of length n)
    """
    Q = np.asarray(Q.todense() if issparse(Q) else Q, dtype=np.float64)
    dQ = np.asarray(dQ.todense() if issparse(dQ) else dQ, dtype=np.float64)
    n = Q.shape[0]
    if pi is None:
        pi = ctmc_solve(Q)
    pi = np.asarray(pi, dtype=np.float64).flatten()

    if dQ.shape != Q.shape:
        raise ValueError("dQ must have the same size as Q")

    # Right-hand side of Eq. (9.81)
    b = -pi @ dQ

    # Solve dpi * Q = b subject to sum(dpi) = 0. Transpose to column form and
    # replace the last equation by the normalization, mirroring ctmc_solve.
    A = Q.T.copy()
    A[n - 1, :] = 1.0
    b = b.flatten().copy()
    b[n - 1] = 0.0

    return np.linalg.solve(A, b)


def ctmc_solve_reducible(Q: np.ndarray, pin: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Solve reducible CTMCs by converting to DTMC via uniformization.

    Port of MATLAB ctmc_solve_reducible.m, which is a thin delegate to
    dtmc_solve_reducible(ctmc_randomization(Q), pin, tol=1e-12). The whole
    convention for a reducible chain therefore lives in dtmc_solve_reducible:
    the limiting vector of a chain with several closed communicating classes is
    NOT unique, and it is resolved by starting uniformly over the SCCs and
    propagating through the lumped limiting matrix (so, when every SCC is
    closed, each class carries equal weight).

    This formerly delegated to ctmc_solve, which just returns whichever null
    vector the linear solver happens to land on -- for two closed classes that
    is all the mass on the first, disagreeing with MATLAB and the JAR.

    Args:
        Q: Infinitesimal generator matrix (possibly reducible)
        pin: Initial probability vector, or None when not available

    Returns:
        Steady-state probability vector
    """
    # Deferred import: dtmc imports ctmc at module level, so a top-level import
    # here would be circular.
    from .dtmc import dtmc_solve_reducible, dtmc_makestochastic

    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]
    if n == 1:
        return np.array([1.0])

    # see _kb/03-api-layer.md for rationale
    qmax = np.max(np.abs(Q))
    q = 1.05 * qmax if qmax > 0 else 1.0
    P = np.eye(n) + Q / q
    P = dtmc_makestochastic(P)

    return dtmc_solve_reducible(P, pin)


def ctmc_solve_reducible_blkdecomp(
    Q: np.ndarray,
    pin: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Solve reducible CTMCs via direct block decomposition on the generator.

    Algorithm:
      1. Decompose states into transient and recurrent classes via SCC
      2. For transient states: solve n * Q_tt = -p0_t for expected sojourn
      3. Compute hitting probabilities: h = n * Q_ta + p0_r
      4. For each recurrent class: solve pi_c * Q_cc = 0, scale by hitting prob

    This avoids the randomization to DTMC used in ctmc_solve_reducible.

    Args:
        Q: Infinitesimal generator matrix (possibly reducible)
        pin: Initial probability distribution (optional)

    Returns:
        Steady-state probability vector
    """
    Q = np.asarray(Q, dtype=np.float64)
    N = Q.shape[0]

    if N == 1:
        return np.array([1.0])

    # Ensure valid generator
    Q = ctmc_makeinfgen(Q)

    # Find strongly connected components
    Adj = Q.copy()
    np.fill_diagonal(Adj, 0.0)
    n_components, scc_labels = connected_components(
        csc_matrix(Adj > 0), directed=True, connection='strong',
        return_labels=True
    )

    # Irreducible case
    if n_components == 1:
        return ctmc_solve(Q)

    num_scc = n_components
    scc_idx = [np.where(scc_labels == i)[0] for i in range(num_scc)]

    # Classify SCCs as recurrent (no outgoing edges) or transient
    is_rec = np.zeros(num_scc, dtype=bool)
    for i in range(num_scc):
        states_i = scc_idx[i]
        outgoing = 0.0
        for s in states_i:
            for j in range(num_scc):
                if j != i:
                    outgoing += np.sum(Adj[s, scc_idx[j]])
        is_rec[i] = outgoing < 1e-10

    trans_scc_ids = np.where(~is_rec)[0]
    rec_scc_ids = np.where(is_rec)[0]

    # Gather ordered state indices
    trans_states = np.sort(np.concatenate([scc_idx[i] for i in trans_scc_ids])
                           ) if len(trans_scc_ids) > 0 else np.array([], dtype=int)
    rec_states = np.sort(np.concatenate([scc_idx[i] for i in rec_scc_ids])
                         ) if len(rec_scc_ids) > 0 else np.array([], dtype=int)
    nt = len(trans_states)
    nr = len(rec_states)

    # Extract Q sub-blocks
    Q_tt = None
    Q_ta = None
    if nt > 0 and nr > 0:
        Q_tt = Q[np.ix_(trans_states, trans_states)]
        Q_ta = Q[np.ix_(trans_states, rec_states)]

    # Compute per-SCC limiting distributions
    pis = np.zeros((num_scc, N))

    for s in range(num_scc):
        class_states = scc_idx[s]
        class_size = len(class_states)

        # Build initial distribution: uniform within SCC s
        p0 = np.zeros(N)
        p0[class_states] = 1.0 / class_size

        # Compute absorption probabilities into recurrent states
        hit = np.zeros(nr)

        if nt > 0 and Q_tt is not None and Q_ta is not None:
            p0_t = p0[trans_states]
            if np.any(np.abs(p0_t) > 0):
                # Solve n * Q_tt = -p0_t  =>  Q_tt' * n' = -p0_t'
                # Q_tt is non-singular (Hurwitz) for transient states
                try:
                    # Above the dispatch threshold the transient block is what
                    # the direct factorization cannot hold; it stays the fallback.
                    sojourn = None
                    if Q_tt.shape[0] > GMRES_MIN_STATES:
                        from .gmres import ctmc_gmres
                        x_gmres, gflag, _, _ = ctmc_gmres(Q_tt.T, -p0_t)
                        if gflag == 0:
                            sojourn = x_gmres
                    if sojourn is None:
                        sojourn = np.linalg.solve(Q_tt.T, -p0_t)
                except np.linalg.LinAlgError:
                    sojourn = np.linalg.lstsq(Q_tt.T, -p0_t, rcond=None)[0]
                hit = sojourn @ Q_ta

        # Add initial mass already in recurrent states
        hit = hit + p0[rec_states]

        # Solve steady-state per recurrent class, scaled by hitting probability
        for c in rec_scc_ids:
            idx_c = scc_idx[c]
            # Map class states to positions in rec_states
            loc = np.searchsorted(rec_states, idx_c)
            reachprob = np.sum(hit[loc])

            if reachprob < 1e-15:
                continue

            if len(idx_c) == 1:
                # Absorbing state: hitting probability IS the final probability
                pis[s, idx_c[0]] = reachprob
            else:
                # Solve pi_c * Q_cc = 0 within this recurrent class
                Q_cc = Q[np.ix_(idx_c, idx_c)]
                pi_c = ctmc_solve(Q_cc)
                pis[s, idx_c] = pi_c * reachprob

    # Compute initial SCC probabilities for weighted average
    if pin is None:
        pinl = np.ones(num_scc)
        # Zero out SCCs containing states with zero column sums (no incoming)
        col_sums = np.sum(np.abs(Q), axis=0)
        for j in np.where(col_sums < 1e-12)[0]:
            pinl[scc_labels[j]] = 0.0
        total_pinl = np.sum(pinl)
        if total_pinl > 0:
            pinl /= total_pinl
        else:
            pinl = np.ones(num_scc) / num_scc
    else:
        pinl = np.zeros(num_scc)
        for i in range(num_scc):
            pinl[i] = np.sum(pin[scc_idx[i]])

    # Weighted average over starting SCCs
    pi = np.zeros(N)
    for i in range(num_scc):
        if pinl[i] > 0:
            pi += pis[i, :] * pinl[i]

    # Special case: single transient SCC without explicit initial distribution
    if len(trans_scc_ids) == 1 and pin is None:
        pi = pis[trans_scc_ids[0], :]

    # Normalize
    total = np.sum(pi)
    if total > 0:
        pi /= total

    return pi


def ctmc_transient(
    Q: np.ndarray,
    initial_dist: np.ndarray,
    time_points: Union[float, np.ndarray],
    method: str = 'expm'
) -> np.ndarray:
    """
    Compute transient probabilities of a CTMC.

    Calculates time-dependent state probabilities π(t) for
    specified time points using matrix exponential methods.

    Args:
        Q: Infinitesimal generator matrix
        initial_dist: Initial probability distribution π(0)
        time_points: Array of time points to evaluate, or single time value
        method: 'expm' for matrix exponential, 'ode' for ODE solver

    Returns:
        Transient probabilities at each time point.
        Shape: (len(time_points), n) if multiple times, (n,) if single time
    """
    Q = np.asarray(Q, dtype=np.float64)
    initial_dist = np.asarray(initial_dist, dtype=np.float64).flatten()

    if np.isscalar(time_points):
        time_points = np.array([time_points])
        single_time = True
    else:
        time_points = np.asarray(time_points)
        single_time = False

    n = Q.shape[0]
    results = np.zeros((len(time_points), n))

    if method == 'expm':
        # Use matrix exponential: π(t) = π(0) * exp(Qt)
        for i, t in enumerate(time_points):
            if t == 0:
                results[i] = initial_dist
            else:
                expQt = linalg.expm(Q * t)
                results[i] = initial_dist @ expQt
    else:
        # Use ODE solver: dπ/dt = π * Q
        def ode_func(t, pi):
            return pi @ Q

        for i, t in enumerate(time_points):
            if t == 0:
                results[i] = initial_dist
            else:
                sol = solve_ivp(
                    ode_func, [0, t], initial_dist,
                    method='LSODA', dense_output=False
                )
                results[i] = sol.y[:, -1]

    return results[0] if single_time else results


def ctmc_timeaverage(pi0: np.ndarray, Q: np.ndarray, t: float,
                     tol: float = 1e-12, maxiter: int = 100):
    """
    Time-averaged transient distribution of a CTMC over [0, t] via uniformization.

    Companion of the endpoint pi0*exp(Q*t); additionally returns the time average

        piTimeAvg = pi0 * (1/t) * \\int_0^t exp(Q*tau) d(tau)

    as well as the endpoint piExit = pi0*exp(Q*t), both from the same Jensen
    uniformization series. Used by the SolverENV state-vector analyzer
    (deterministic-sojourn option). Mirrors matlab ctmc_timeaverage.m.

    Returns:
        (piTimeAvg, piExit) as 1D arrays.
    """
    Q = np.asarray(Q, dtype=np.float64)
    pi0 = np.asarray(pi0, dtype=np.float64).ravel()
    n = Q.shape[0]
    q = 1.1 * np.max(np.abs(np.diag(Q)))
    Qs = np.eye(n) + Q / q
    qt = q * t

    # Number of Poisson terms needed (right-tail below tol).
    k = 0
    s = 1.0
    r = 1.0
    it = 0
    kmax = 1
    while it < maxiter:
        it += 1
        k += 1
        r = r * qt / k
        s += r
        if 1 - np.exp(-qt) * s <= tol:
            kmax = k
            break

    w = np.exp(-qt)   # Poisson PMF w_0
    W = w             # Poisson CDF W_0
    P = pi0.copy()    # pi0*P^0
    piExit = w * P
    piIntSum = max(1 - W, 0.0) * P
    for j in range(1, kmax + 1):
        P = P @ Qs
        w = w * qt / j
        W += w
        piExit = piExit + w * P
        piIntSum = piIntSum + max(1 - W, 0.0) * P
    piTimeAvg = piIntSum / qt
    return piTimeAvg, piExit


def ctmc_uniformization(
    Q: np.ndarray,
    lambda_rate: Optional[float] = None
) -> Dict[str, Any]:
    """
    Uniformize CTMC generator matrix.

    Converts CTMC to an equivalent uniformized discrete-time chain
    for numerical analysis and simulation purposes.

    The uniformized DTMC has transition matrix P = I + Q/λ where
    λ is the uniformization rate (max exit rate).

    Args:
        Q: Infinitesimal generator matrix
        lambda_rate: Uniformization rate (optional, auto-computed if None)

    Returns:
        dict containing:
            - 'P': Uniformized transition matrix
            - 'lambda': Uniformization rate
    """
    Q = np.asarray(Q, dtype=np.float64)

    if lambda_rate is None:
        # Use max exit rate (max |diagonal|)
        lambda_rate = -np.min(np.diag(Q))
        if lambda_rate <= 0:
            lambda_rate = 1.0

    n = Q.shape[0]
    I = np.eye(n)
    P = I + Q / lambda_rate

    return {
        'P': P,
        'lambda': lambda_rate
    }


def ctmc_randomization(
    Q: np.ndarray,
    initial_dist: np.ndarray,
    time_points: np.ndarray,
    precision: float = 1e-10
) -> np.ndarray:
    """
    Compute CTMC transient probabilities using randomization.

    Uses Jensen's randomization method (uniformization) to compute
    transient probabilities by converting the CTMC to a uniformized DTMC.

    This method is numerically stable and avoids matrix exponentials.

    Args:
        Q: Infinitesimal generator matrix
        initial_dist: Initial probability distribution
        time_points: Array of time points to evaluate
        precision: Numerical precision for truncation (Poisson tail)

    Returns:
        Transient probabilities at each time point
    """
    Q = np.asarray(Q, dtype=np.float64)
    initial_dist = np.asarray(initial_dist, dtype=np.float64).flatten()
    time_points = np.asarray(time_points)

    n = Q.shape[0]

    # Uniformization rate
    lambda_rate = -np.min(np.diag(Q))
    if lambda_rate <= 0:
        lambda_rate = 1.0

    # Uniformized transition matrix
    P = np.eye(n) + Q / lambda_rate

    results = np.zeros((len(time_points), n))

    for idx, t in enumerate(time_points):
        if t == 0:
            results[idx] = initial_dist
            continue

        # Compute Poisson probabilities and truncation point
        q = lambda_rate * t

        # Find truncation point k_max such that sum of Poisson tail < precision
        from scipy.stats import poisson
        k_max = int(poisson.ppf(1 - precision, q)) + 10

        # Compute Poisson probabilities
        poisson_probs = poisson.pmf(np.arange(k_max + 1), q)

        # Compute π(t) = Σ_k P(N(t)=k) * π(0) * P^k
        pi_t = np.zeros(n)
        pi_k = initial_dist.copy()  # π(0) * P^0 = π(0)

        for k in range(k_max + 1):
            pi_t += poisson_probs[k] * pi_k
            pi_k = pi_k @ P  # π(0) * P^(k+1)

        results[idx] = pi_t

    return results


def ctmc_stochcomp(
    Q: np.ndarray,
    I: Optional[np.ndarray] = None
) -> Dict[str, np.ndarray]:
    """
    Compute stochastic complement of CTMC.

    Reduces the CTMC by eliminating states not in I while preserving
    the steady-state distribution restricted to the kept states.

    Args:
        Q: Infinitesimal generator matrix
        I: States to retain (array of indices). If None, defaults to
           0..ceil(n/2)-1 (matching JAR/MATLAB).

    Returns:
        dict containing:
            - 'S': Stochastic complement (reduced generator)
            - 'Q11': Submatrix for kept states
            - 'Q12': Transitions from kept to eliminated
            - 'Q21': Transitions from eliminated to kept
            - 'Q22': Submatrix for eliminated states
            - 'T': Transient contribution matrix
    """
    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]

    if I is None:
        I = np.arange(int(np.ceil(n / 2)))
    else:
        I = np.asarray(I, dtype=int).flatten()

    keep_set = set(I.tolist())
    Ic = np.array([j for j in range(n) if j not in keep_set], dtype=int)

    if len(Ic) == 0:
        return {
            'S': Q[np.ix_(I, I)],
            'Q11': Q[np.ix_(I, I)],
            'Q12': np.zeros((len(I), 0)),
            'Q21': np.zeros((0, len(I))),
            'Q22': np.zeros((0, 0)),
            'T': np.zeros((len(I), len(I)))
        }

    Q11 = Q[np.ix_(I, I)]
    Q12 = Q[np.ix_(I, Ic)]
    Q21 = Q[np.ix_(Ic, I)]
    Q22 = Q[np.ix_(Ic, Ic)]

    # S = Q11 + Q12 * (-Q22)^{-1} * Q21
    Q22_neg = -Q22

    # see _kb/03-api-layer.md for rationale
    zero_diag = np.where(np.abs(np.diag(Q22_neg)) < 1e-10)[0]
    if zero_diag.size > 0:
        Q22_neg = Q22_neg.copy()
        Q22_neg[zero_diag, zero_diag] = 1.0

    # see _kb/03-api-layer.md for rationale
    T = None
    if Q22_neg.shape[0] > GMRES_MIN_STATES:
        from .gmres import ctmc_gmres_multi
        T_it, gflag = ctmc_gmres_multi(Q22_neg, Q21)
        if gflag == 0:
            T = T_it
    if T is None:
        try:
            T = linalg.solve(Q22_neg, Q21)
        except LinAlgError:
            # see _kb/03-api-layer.md for rationale
            T = linalg.lstsq(Q22_neg, Q21)[0]

    T = Q12 @ T
    S = Q11 + T

    return {
        'S': S,
        'Q11': Q11,
        'Q12': Q12,
        'Q21': Q21,
        'Q22': Q22,
        'T': T
    }


def ctmc_timereverse(
    Q: np.ndarray,
    pi: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Compute time-reversed CTMC generator.

    The time-reversed generator Q* has elements:
    Q*_{ij} = π_j * Q_{ji} / π_i

    Args:
        Q: Original infinitesimal generator matrix
        pi: Steady-state distribution (optional, computed if None)

    Returns:
        Time-reversed generator matrix
    """
    Q = np.asarray(Q, dtype=np.float64)

    if pi is None:
        pi = ctmc_solve(Q)
    else:
        pi = np.asarray(pi, dtype=np.float64).flatten()

    n = Q.shape[0]

    Q_rev = np.zeros_like(Q)

    for i in range(n):
        for j in range(n):
            if pi[i] > 0:
                Q_rev[i, j] = pi[j] * Q[j, i] / pi[i]

    # Ensure valid generator
    Q_rev = ctmc_makeinfgen(Q_rev)

    return Q_rev


def ctmc_rand(
    n: int,
    density: float = 0.3,
    max_rate: float = 10.0
) -> np.ndarray:
    """
    Generate random CTMC generator matrix.

    Args:
        n: Number of states
        density: Sparsity density (0 to 1, default 0.3)
        max_rate: Maximum transition rate (default 10.0)

    Returns:
        Random infinitesimal generator matrix
    """
    # Generate random off-diagonal entries
    Q = np.random.rand(n, n) * max_rate

    # Apply density mask
    mask = np.random.rand(n, n) < density
    Q = Q * mask

    # Zero out diagonal
    np.fill_diagonal(Q, 0.0)

    # Make valid generator
    Q = ctmc_makeinfgen(Q)

    return Q


def ctmc_simulate(
    Q: np.ndarray,
    initial_state: int,
    max_time: float,
    max_events: int = 10000,
    seed: Optional[int] = None
) -> Dict[str, np.ndarray]:
    """
    Simulate CTMC sample path using Gillespie algorithm.

    Generates a realization of the continuous-time Markov chain
    using the next-reaction method.

    Args:
        Q: Infinitesimal generator matrix
        initial_state: Starting state (integer index)
        max_time: Maximum simulation time
        max_events: Maximum number of transitions (default: 10000)
        seed: Random seed for reproducibility (optional)

    Returns:
        dict with:
            - 'states': Array of visited states
            - 'times': Array of transition times
            - 'sojourn_times': Time spent in each state
    """
    if seed is not None:
        np.random.seed(seed)

    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]

    states = [initial_state]
    times = [0.0]
    sojourn_times = []

    current_state = initial_state
    current_time = 0.0

    for _ in range(max_events):
        # Exit rate from current state
        exit_rate = -Q[current_state, current_state]

        if exit_rate <= 0:
            # Absorbing state
            sojourn_times.append(max_time - current_time)
            break

        # Time to next transition (exponential)
        sojourn = np.random.exponential(1.0 / exit_rate)

        if current_time + sojourn > max_time:
            sojourn_times.append(max_time - current_time)
            break

        sojourn_times.append(sojourn)
        current_time += sojourn

        # Determine next state
        rates = Q[current_state, :].copy()
        rates[current_state] = 0.0
        probs = rates / rates.sum()
        next_state = np.random.choice(n, p=probs)

        states.append(next_state)
        times.append(current_time)
        current_state = next_state

    return {
        'states': np.array(states),
        'times': np.array(times),
        'sojourn_times': np.array(sojourn_times)
    }


def ctmc_isfeasible(Q: np.ndarray, tolerance: float = 1e-10) -> bool:
    """
    Check if matrix is a valid CTMC infinitesimal generator.

    Validates:
    - Off-diagonal elements are non-negative
    - Row sums are zero
    - Diagonal elements are non-positive

    Args:
        Q: Candidate generator matrix
        tolerance: Numerical tolerance (default: 1e-10)

    Returns:
        True if matrix is valid CTMC generator
    """
    Q = np.asarray(Q)
    n = Q.shape[0]

    if Q.shape[0] != Q.shape[1]:
        return False

    # Check off-diagonal non-negative
    for i in range(n):
        for j in range(n):
            if i != j and Q[i, j] < -tolerance:
                return False

    # Check diagonal non-positive
    if np.any(np.diag(Q) > tolerance):
        return False

    # Check row sums are zero
    row_sums = Q.sum(axis=1)
    if np.any(np.abs(row_sums) > tolerance):
        return False

    return True


# ============================================================================
# State Space Generation
# ============================================================================

@dataclass
class CtmcSsgResult:
    """Result from CTMC state space generation."""
    state_space: np.ndarray  # Complete state space matrix
    state_space_aggr: np.ndarray  # Aggregated state space (per station-class)
    state_space_hashed: np.ndarray  # Hashed state indices
    node_state_space: Dict[int, np.ndarray]  # State space per node
    sn: Any  # Updated network structure


def ctmc_ssg(sn: Any, options: Optional[Dict] = None) -> CtmcSsgResult:
    """
    Generate complete CTMC state space for a queueing network.

    Creates all possible network states including those not reachable from
    the initial state. For open classes, a cutoff parameter limits the
    maximum population to keep state space finite.

    The state space is aggregated to show per-station-class job counts.

    Args:
        sn: NetworkStruct object (from getStruct())
        options: Solver options dict with fields:
            - cutoff: Population cutoff for open classes (required if open)
            - config.hide_immediate: Hide immediate transitions (default True)

    Returns:
        CtmcSsgResult containing:
            - state_space: Complete state space matrix (rows=states, cols=state components)
            - state_space_aggr: Aggregated state space (rows=states, cols=stations*classes)
            - state_space_hashed: Hashed state indices for lookup
            - node_state_space: Dictionary of per-node state spaces
            - sn: Updated network structure with space field populated

    References:
        MATLAB: matlab/src/api/mc/ctmc_ssg.m
    """
    from ..state import spaceGenerator, toMarginal

    if options is None:
        options = {}

    # Get cutoff from options
    cutoff = options.get('cutoff', None)

    # Generate state space
    state_space, state_space_hashed, sn, adj, st = spaceGenerator(sn, cutoff, options)

    # Get node state space from sn.space
    node_state_space = {}
    if hasattr(sn, 'space') and sn.space is not None:
        if isinstance(sn.space, dict):
            node_state_space = sn.space
        elif isinstance(sn.space, (list, np.ndarray)):
            for i, space in enumerate(sn.space):
                if space is not None:
                    node_state_space[i] = space

    # Set default hide_immediate
    if 'hide_immediate' not in options:
        options['hide_immediate'] = True

    # Compute aggregated state space
    nstateful = sn.nstateful if hasattr(sn, 'nstateful') else sn.nstations
    nclasses = sn.nclasses if hasattr(sn, 'nclasses') else 1

    # Initialize aggregated state space
    if state_space_hashed.size > 0:
        state_space_aggr = np.zeros((state_space_hashed.shape[0], sn.nstations * nclasses))
    else:
        state_space_aggr = np.zeros((0, sn.nstations * nclasses))

    # Process synchronizations to compute aggregated states
    sync = sn.sync if hasattr(sn, 'sync') else []
    A = len(sync) if sync else 0

    for s in range(state_space_hashed.shape[0] if state_space_hashed.size > 0 else 0):
        state = state_space_hashed[s, :] if state_space_hashed.ndim > 1 else state_space_hashed

        # Update state cell array
        state_cell = {}
        for ind in range(sn.nnodes if hasattr(sn, 'nnodes') else sn.nstations):
            isstateful = sn.isstateful(ind) if hasattr(sn, 'isstateful') else True
            if isstateful:
                isf = sn.nodeToStateful[ind] if hasattr(sn, 'nodeToStateful') else ind
                if isf < len(node_state_space) and node_state_space.get(isf) is not None:
                    space_isf = node_state_space[isf]
                    state_idx = int(state[isf]) if isf < len(state) else 0
                    if state_idx < len(space_isf):
                        state_cell[isf] = space_isf[state_idx]

                isstation = sn.isstation(ind) if hasattr(sn, 'isstation') else True
                if isstation:
                    ist = sn.nodeToStation[ind] if hasattr(sn, 'nodeToStation') else ind
                    if isf in state_cell:
                        # Compute marginal job counts
                        try:
                            ni, nir, sir, kir = toMarginal(sn, ind, state_cell[isf])
                            # Store in aggregated state space
                            start_col = ist * nclasses
                            end_col = (ist + 1) * nclasses
                            if nir.ndim == 1:
                                state_space_aggr[s, start_col:end_col] = nir
                            else:
                                state_space_aggr[s, start_col:end_col] = nir[0, :]
                        except Exception:
                            pass

    return CtmcSsgResult(
        state_space=state_space,
        state_space_aggr=state_space_aggr,
        state_space_hashed=state_space_hashed,
        node_state_space=node_state_space,
        sn=sn
    )


def ctmc_ssg_reachability(sn: Any, options: Optional[Dict] = None) -> CtmcSsgResult:
    """
    Generate reachable CTMC state space for a queueing network.

    Creates only the states reachable from the initial state through valid
    transitions. This is more efficient than ctmc_ssg for networks with
    constrained reachability.

    Args:
        sn: NetworkStruct object (from getStruct())
        options: Solver options dict with fields:
            - config.hide_immediate: Hide immediate transitions (default True)

    Returns:
        CtmcSsgResult containing:
            - state_space: Reachable state space matrix
            - state_space_aggr: Aggregated state space (per station-class)
            - state_space_hashed: Hashed state indices
            - node_state_space: Dictionary of per-node state spaces
            - sn: Updated network structure

    References:
        MATLAB: matlab/src/api/mc/ctmc_ssg_reachability.m
    """
    from ..state import spaceGenerator, toMarginal

    if options is None:
        options = {}

    # Set default hide_immediate in config
    config = options.get('config', {})
    if 'hide_immediate' not in config:
        config['hide_immediate'] = True
        options['config'] = config

    # see _kb/03-api-layer.md for rationale
    state_space, state_space_hashed, sn, adj, st = spaceGenerator(sn, options.get('cutoff'), options)

    # Get node state space
    node_state_space = {}
    if hasattr(sn, 'space') and sn.space is not None:
        if isinstance(sn.space, dict):
            node_state_space = sn.space
        elif isinstance(sn.space, (list, np.ndarray)):
            for i, space in enumerate(sn.space):
                if space is not None:
                    node_state_space[i] = space

    # Compute aggregated state space
    nstateful = sn.nstateful if hasattr(sn, 'nstateful') else sn.nstations
    nclasses = sn.nclasses if hasattr(sn, 'nclasses') else 1

    if state_space_hashed.size > 0:
        state_space_aggr = np.zeros((state_space_hashed.shape[0], sn.nstations * nclasses))
    else:
        state_space_aggr = np.zeros((0, sn.nstations * nclasses))

    # Process states to compute aggregated representation
    sync = sn.sync if hasattr(sn, 'sync') else []
    A = len(sync) if sync else 0

    for s in range(state_space_hashed.shape[0] if state_space_hashed.size > 0 else 0):
        state = state_space_hashed[s, :] if state_space_hashed.ndim > 1 else state_space_hashed

        state_cell = {}
        for ind in range(sn.nnodes if hasattr(sn, 'nnodes') else sn.nstations):
            isstateful = sn.isstateful(ind) if hasattr(sn, 'isstateful') else True
            if isstateful:
                isf = sn.nodeToStateful[ind] if hasattr(sn, 'nodeToStateful') else ind
                if isf < len(node_state_space) and node_state_space.get(isf) is not None:
                    space_isf = node_state_space[isf]
                    state_idx = int(state[isf]) if isf < len(state) else 0
                    if state_idx < len(space_isf):
                        state_cell[isf] = space_isf[state_idx]

                isstation = sn.isstation(ind) if hasattr(sn, 'isstation') else True
                if isstation:
                    ist = sn.nodeToStation[ind] if hasattr(sn, 'nodeToStation') else ind
                    if isf in state_cell:
                        try:
                            ni, nir, sir, kir = toMarginal(sn, ind, state_cell[isf])
                            start_col = ist * nclasses
                            end_col = (ist + 1) * nclasses
                            if nir.ndim == 1:
                                state_space_aggr[s, start_col:end_col] = nir
                            else:
                                state_space_aggr[s, start_col:end_col] = nir[0, :]
                        except Exception:
                            pass

    return CtmcSsgResult(
        state_space=state_space,
        state_space_aggr=state_space_aggr,
        state_space_hashed=state_space_hashed,
        node_state_space=node_state_space,
        sn=sn
    )


def ctmc_memory_gate(log_nstates: float, force: bool = False,
                      verbose: bool = False, safety_fraction: float = 0.6) -> Tuple[bool, str]:
    """
    Hardware-aware, profiling-calibrated CTMC memory pre-gate.

    Decides whether a CTMC steady-state solve of a state space of worst-case
    size exp(log_nstates) is safe on the current host. The budget is a fraction
    of available memory; the per-state cost is calibrated by profiling sparse LU
    factorization and cached per machine.

    Args:
        log_nstates: log(number of states) in the CTMC
        force: If True, override the memory limit and proceed anyway
        verbose: If True, print calibration and memory predictions
        safety_fraction: Fraction of available memory to use as safe budget (default 0.6)

    Returns:
        Tuple (ok, msg) where:
            - ok (bool): True if solve is safe, False if memory exceeded (and force=False)
            - msg (str): Status or warning message
    """
    try:
        import psutil
        avail = psutil.virtual_memory().available
    except ImportError:
        import resource
        avail = resource.getrlimit(resource.RLIMIT_AS)[0]
        if avail <= 0:
            avail = 16 * 1024**3

    budget = safety_fraction * avail
    calib = _ctmc_get_calibration(verbose)

    log_pred = np.log(calib['alpha_mem']) + calib['beta_mem'] * log_nstates
    log_budget = np.log(max(budget, 1.0))
    pred_gb = np.exp(min(log_pred, 700)) / (1024**3)
    budget_gb = budget / (1024**3)

    if log_pred > log_budget:
        msg = (f"CTMC predicted peak memory ~{pred_gb:.2f} GB exceeds the safe "
               f"budget ~{budget_gb:.2f} GB ({100*safety_fraction:.0f}% of {avail/(1024**3):.2f} GB available). "
               f"Reduce the state space (e.g. lower cutoff), use another solver (MVA/NC/FLD), or "
               f"set force=True to override.")
        if not force:
            return (False, msg)
        if verbose:
            print(f"Warning (forced): {msg}")
    elif verbose and log_pred > np.log(max(0.5 * budget, 1.0)):
        print(f"CTMC predicted peak memory ~{pred_gb:.2f} GB (budget ~{budget_gb:.2f} GB).")

    return (True, "")


def _ctmc_get_calibration(verbose: bool = False) -> Dict[str, float]:
    """Get or compute calibrated power-law coefficients for memory/time prediction."""
    BYTES_PER_NZ = 16
    G1, G2 = 40, 80
    FALLBACK_ALPHA = BYTES_PER_NZ * 8
    FALLBACK_BETA = 1.3

    sig = _ctmc_machine_signature()
    cachefile = os.path.join(tempfile.gettempdir(), 'line_ctmc_calib_python.json')

    if os.path.exists(cachefile):
        try:
            with open(cachefile, 'r') as f:
                data = json.load(f)
            if data.get('sig') == sig:
                return data
        except Exception:
            pass

    try:
        n1, b1, t1 = _ctmc_profile_point(G1, BYTES_PER_NZ)
        n2, b2, t2 = _ctmc_profile_point(G2, BYTES_PER_NZ)
        am, bm = _ctmc_fit_power_law(n1, b1, n2, b2)
        at, bt = _ctmc_fit_power_law(n1, max(t1, 1e-9), n2, max(t2, 1e-9))
        calib = {
            'sig': sig,
            'alpha_mem': float(am),
            'beta_mem': float(bm),
            'alpha_t': float(at),
            'beta_t': float(bt),
            'timestamp': time.time()
        }
        try:
            with open(cachefile, 'w') as f:
                json.dump(calib, f)
        except Exception:
            pass
        if verbose:
            print(f"CTMC calibration: bytes ~ {am:.3g}*N^{bm:.3f}")
        return calib
    except Exception as e:
        if verbose:
            print(f"CTMC calibration failed ({str(e)}); using fallback model")
        return {
            'sig': sig,
            'alpha_mem': FALLBACK_ALPHA,
            'beta_mem': FALLBACK_BETA,
            'alpha_t': 0.0,
            'beta_t': 1.0,
            'timestamp': time.time()
        }


def _ctmc_machine_signature() -> str:
    """Generate a machine signature for caching calibration data."""
    try:
        ncores = os.cpu_count() or 1
    except Exception:
        ncores = 1
    arch = platform.machine()
    return f"{arch}|{ncores}|{platform.platform()}"


def _ctmc_profile_point(g: int, bytes_per_nz: int = 16) -> Tuple[int, float, float]:
    """
    Profile sparse LU factorization on a g*g lattice generator.

    Returns: (n_states, bytes_used, time_seconds)
    """
    A = _ctmc_lattice_generator(g)
    n = A.shape[0]

    A_csc = A.tocsc()
    t0 = time.time()
    lu_result = splu(A_csc)
    elapsed = time.time() - t0

    L_nnz = lu_result.L.nnz
    U_nnz = lu_result.U.nnz
    bytes_used = bytes_per_nz * (L_nnz + U_nnz)

    return n, bytes_used, elapsed


def _ctmc_lattice_generator(g: int):
    """Create a g*g nearest-neighbor lattice generator (QBD-like structure)."""
    n = g * g
    idx = np.arange(n).reshape(g, g)

    src_list, dst_list = [], []

    for row in range(g):
        for col in range(g - 1):
            src_list.append(idx[row, col])
            dst_list.append(idx[row, col + 1])
            src_list.append(idx[row, col + 1])
            dst_list.append(idx[row, col])

    for row in range(g - 1):
        for col in range(g):
            src_list.append(idx[row, col])
            dst_list.append(idx[row + 1, col])
            src_list.append(idx[row + 1, col])
            dst_list.append(idx[row, col])

    Q = csc_matrix((np.ones(len(src_list)), (src_list, dst_list)), shape=(n, n))
    rowsum = np.array(Q.sum(axis=1)).flatten()
    Q = Q - diags(rowsum)

    return Q[:-1, :-1]


def _ctmc_fit_power_law(x1: float, y1: float, x2: float, y2: float) -> Tuple[float, float]:
    """Fit power-law coefficients: y = alpha * x^beta."""
    if x1 <= 0 or x2 <= 0 or y1 <= 0 or y2 <= 0 or x1 == x2:
        raise ValueError("Degenerate power-law fit")
    beta = np.log(y2 / y1) / np.log(x2 / x1)
    alpha = y1 / (x1 ** beta)
    return alpha, beta


__all__ = [
    'ctmc_solve',
    'ctmc_solve_reducible',
    'ctmc_solve_reducible_blkdecomp',
    'ctmc_makeinfgen',
    'ctmc_transient',
    'ctmc_timeaverage',
    'ctmc_uniformization',
    'ctmc_randomization',
    'ctmc_stochcomp',
    'ctmc_timereverse',
    'ctmc_rand',
    'ctmc_simulate',
    'ctmc_isfeasible',
    'ctmc_ssg',
    'ctmc_ssg_reachability',
    'ctmc_memory_gate',
    'CtmcSsgResult',
]
