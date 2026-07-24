"""
DAC (Distribution Analysis by Chain) Method for Joint Queue-Length Distributions.

Native Python implementation of the DAC algorithm, which computes the whole set
of joint queue-length probabilities of a closed product-form queueing network
with single-server fixed-rate, infinite-server and queue-dependent centers.

Key functions:
    pfqn_dac: DAC method for joint queue-length distributions

References:
    E. de Souza e Silva, "Distribution Analysis of Product Form Queueing
    Networks", UCLA Computer Science Department, CSD-870023, April 1987.

    Original MATLAB: matlab/src/api/pfqn/pfqn_dac.m
"""

from typing import List, Optional, Tuple

import numpy as np


def _dac_compositions(J: int, k: int) -> np.ndarray:
    """All J-part compositions of k, in lexicographic order."""
    if J == 1:
        return np.array([[k]], dtype=np.int64)
    blocks = []
    for v in range(k + 1):
        B = _dac_compositions(J - 1, k - v)
        blocks.append(np.hstack((np.full((B.shape[0], 1), v, dtype=np.int64), B)))
    return np.vstack(blocks)


def _dac_lattice(J: int, Nt: int) -> Tuple[List[np.ndarray], List[Optional[np.ndarray]]]:
    """Enumerate the aggregate state space level by level and precompute, for
    each state at level k, the index of that state with one more job at center j."""
    lvstates = [_dac_compositions(J, k) for k in range(Nt + 1)]

    # Binomial table, C[a, b] = binom(a, b), for ranking compositions.
    C = np.zeros((Nt + J + 1, J + 2), dtype=np.int64)
    for a in range(Nt + J):
        for b in range(min(a, J) + 1):
            if b == 0:
                C[a, b] = 1
            else:
                C[a, b] = C[a - 1, b - 1] + C[a - 1, b]

    def rank(n, k):
        idx = 0
        rem = k
        for j in range(J - 1):
            for v in range(n[j]):
                idx += C[(rem - v) + (J - 1 - j) - 1, (J - 1 - j) - 1]
            rem -= n[j]
        return idx

    succ: List[Optional[np.ndarray]] = [None] * (Nt + 1)
    for k in range(Nt):
        S = lvstates[k]
        ix = np.zeros((S.shape[0], J), dtype=np.int64)
        for i in range(S.shape[0]):
            for j in range(J):
                t = S[i].copy()
                t[j] += 1
                ix[i, j] = rank(t, k + 1)
        succ[k] = ix
    return lvstates, succ


def _dac_step(p, r, mu, k, lvstates, succ):
    """One step of the recursion: add a single customer with demands r to a
    network holding k customers. Returns the distribution at level k+1, the
    throughput of the added customer, and its per-center presence probabilities."""
    J = len(r)
    S = lvstates[k]
    # Marginal queue lengths of the network with k customers.
    marg = np.zeros((J, k + 1))
    for j in range(J):
        marg[j] = np.bincount(S[:, j], weights=p, minlength=k + 1)
    c = np.zeros(J)
    for j in range(J):
        for n in range(1, k + 2):
            c[j] += (n / mu[j, n - 1]) * marg[j, n - 1]
    lam = 1.0 / float(np.dot(r, c))
    Lq = lam * r * c

    ix = succ[k]
    pn = np.zeros(lvstates[k + 1].shape[0])
    for j in range(J):
        if r[j] <= 0:
            continue
        nj = S[:, j] + 1
        w = lam * r[j] * (nj / mu[j, nj - 1]) * p
        pn += np.bincount(ix[:, j], weights=w, minlength=pn.shape[0])
    return pn, lam, Lq


def pfqn_dac(L: np.ndarray, N: np.ndarray,
             Z: Optional[np.ndarray] = None,
             mu: Optional[np.ndarray] = None):
    """
    DAC (Distribution Analysis by Chain) method for joint queue-length distributions.

    Computes the joint queue-length distribution of a closed product-form
    network by a chain-by-chain recursion over a related network in which every
    chain holds a single customer, a transformation that leaves the aggregate
    queue-length distribution unchanged. Given the distribution of a network
    with k-1 such chains, adding one customer of a chain with demands r gives

        c_j      = sum_{n=1..k} (n/mu_j(n)) * P_j^{k-1}(n-1)
        lambda_k = 1 / sum_j r_j c_j
        P^k(n)   = lambda_k * sum_j r_j (n_j/mu_j(n_j)) * P^{k-1}(n-e_j)

    where lambda_k is the throughput of the customer being added and 1/c_j is
    the throughput of a chain visiting center j only. The recursion conserves
    probability mass by construction, hence it is numerically stable.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,), optional (default: zeros). If sum(Z)>0 an
           extra infinite-server station is appended, so that states has M+1
           columns and its last column holds the think-station population.
        mu: Load-dependent rate matrix (M x Nt), Nt=sum(N), optional
            (default: ones, i.e. single-server fixed rate). Use mu[j,:]=1..Nt
            for infinite server and mu[j,n]=min(n+1,c) for a c-server station.

    Returns:
        Tuple of (Pjoint, states, XN, QN, UN, CN, pi):
            Pjoint: Probability of the aggregate state in the corresponding row
                    of states, of length nchoosek(Nt+J-1, J-1)
            states: Aggregate states (S x J), states[s,j] = jobs at center j
            XN: Throughput of chain r (R,)
            QN: Mean number of chain-r customers at station j (M x R)
            UN: Utilization of station j (M,), i.e. 1-P_j(0)
            CN: Cycle time of chain r, exclusive of think time (R,)
            pi: Marginal probabilities, pi[j,n] = P(n jobs at station j),
                shape (M, Nt+1)

    References:
        E. de Souza e Silva, "Distribution Analysis of Product Form Queueing
        Networks", UCLA Computer Science Department, CSD-870023, April 1987.
    """
    L = np.asarray(L, dtype=float)
    if L.ndim == 1:
        L = L.reshape(-1, 1)
    M, R = L.shape
    N = np.asarray(N, dtype=float).flatten()
    if Z is None:
        Z = np.zeros(R)
    Z = np.asarray(Z, dtype=float).flatten()
    Nt = int(round(N.sum()))

    if np.any(N < 0):
        raise ValueError("pfqn_dac: population vector must be non-negative.")
    if mu is None:
        mu = np.ones((M, max(Nt, 1)))
    mu = np.asarray(mu, dtype=float)
    if mu.ndim == 1:
        mu = mu.reshape(1, -1)
    if mu.shape[0] != M:
        raise ValueError("pfqn_dac: the mu matrix must have one row per station.")
    if Nt > 0 and mu.shape[1] < Nt:
        raise ValueError("pfqn_dac: the mu matrix must have at least sum(N) columns.")

    # A non-zero think time is modelled as an appended infinite-server center.
    hasZ = Z.sum() > 0
    W = max(Nt, 1)
    if hasZ:
        Lx = np.vstack((L, Z.reshape(1, -1)))
        mux = np.vstack((mu[:, :W], np.arange(1, W + 1, dtype=float).reshape(1, -1)))
    else:
        Lx = L
        mux = mu[:, :W]
    J = Lx.shape[0]

    if Nt == 0:
        pi = np.zeros((M, 1))
        pi[:, 0] = 1.0
        return (np.array([1.0]), np.zeros((1, J), dtype=np.int64), np.zeros(R),
                np.zeros((M, R)), np.zeros(M), np.zeros(R), pi)

    active = [r for r in range(R) if N[r] > 0]
    D = len(active)
    for r in active:
        if not np.any(Lx[:, r] > 0):
            raise ValueError("pfqn_dac: chain %d has null demand at every center." % (r + 1))

    lvstates, succ = _dac_lattice(J, Nt)

    # Chains are ordered so that the last D added are distinct, which lets all
    # per-chain measures reuse the common prefix of the recursion.
    prefix = []
    for r in active:
        prefix += [r] * (int(N[r]) - 1)
    tail = list(active)

    # Prefix of the recursion, shared by every per-chain run.
    p = np.array([1.0])
    k = 0
    for c in prefix:
        p, _, _ = _dac_step(p, Lx[:, c], mux, k, lvstates, succ)
        k += 1

    # Base run over the tail, saving the intermediate distributions S_0..S_{D-1}.
    Sp = [None] * D
    Sp[0] = p
    pb = p
    kb = k
    lam = 0.0
    Lq = np.zeros(J)
    for idx in range(D):
        pb, lam, Lq = _dac_step(pb, Lx[:, tail[idx]], mux, kb, lvstates, succ)
        kb += 1
        if idx < D - 1:
            Sp[idx + 1] = pb
    Pjoint = pb
    states = lvstates[Nt]

    XN = np.zeros(R)
    QN = np.zeros((M, R))
    # The base run already places chain tail[D-1] last.
    r = tail[D - 1]
    XN[r] = N[r] * lam
    QN[:, r] = N[r] * Lq[:M]

    # Re-run the tail with chain tail[idx] moved last, restarting from S_idx.
    for idx in range(D - 1):
        pc = Sp[idx]
        kc = k + idx
        order = tail[idx + 1:] + [tail[idx]]
        for c in order:
            pc, lam, Lq = _dac_step(pc, Lx[:, c], mux, kc, lvstates, succ)
            kc += 1
        r = tail[idx]
        XN[r] = N[r] * lam
        QN[:, r] = N[r] * Lq[:M]

    # Marginal queue-length probabilities at the full population.
    pi = np.zeros((M, Nt + 1))
    for j in range(M):
        pi[j] = np.bincount(states[:, j], weights=Pjoint, minlength=Nt + 1)
    UN = 1.0 - pi[:, 0]
    CN = np.zeros(R)
    for r in active:
        CN[r] = N[r] / XN[r] - Z[r]

    return Pjoint, states, XN, QN, UN, CN, pi
