"""
MVAC (Mean Value Analysis by Chain) for closed product-form queueing networks.

Native Python implementation of the exact mean value analysis by chain of
Conway, de Souza e Silva and Lavenberg, which recurs on the chains rather than
on the population vector and is therefore attractive for networks with few
service centers and many chains.

Key functions:
    pfqn_mvac: MVAC exact solution of a closed multichain product-form network

References:
    A. E. Conway, E. de Souza e Silva and S. S. Lavenberg, "Mean Value Analysis
    by Chain of Product Form Queueing Networks", IEEE Trans. Computers,
    38(3):432-442, 1989.
    Original MATLAB: matlab/src/api/pfqn/pfqn_mvac.m
"""

import numpy as np
from typing import Optional, Tuple

from .utils import multichoose

__all__ = ['pfqn_mvac']


def _part1(k0, K, J, J1, a, ak, Vlist, vsum, succ, Lall, Ljkall, lamall):
    """
    Part 1 of the MVAC basic step: evaluate (10), (9a)-(9b) and (7) for
    k = k0,...,K over v in V_k = I_k u ... u I_K.
    """
    nv = Vlist.shape[0]
    for k in range(k0, K + 1):
        idxV = np.flatnonzero(vsum <= K - k)
        Lp = Lall[k - 1]  # L^{k-1}
        Vv = Vlist[idxV, :]
        Lpv = Lp[:, idxV].T  # L^{k-1}_j(v)
        ajk = a[:, k - 1]  # a_jk
        # (10): the chain-k customer is at some center, so sum_j L^k_{jk}(v) = 1
        den = ak[k - 1] + np.sum((Lpv[:, :J1] + Vv[:, :J1]) * ajk[:J1], axis=1)
        lam = 1.0 / den
        Ljkv = np.zeros((len(idxV), J))
        Ljkv[:, :J1] = (lam[:, None] * (1.0 + Lpv[:, :J1] + Vv[:, :J1])) * ajk[:J1]  # (9a)
        Ljkv[:, J1:J] = lam[:, None] * ajk[J1:J]  # (9b)
        # (7)
        Lkv = Ljkv.copy()
        for j in range(J):
            Lkv = Lkv + Ljkv[:, j][:, None] * Lp[:, succ[idxV, j]].T
        Lk = np.zeros((J, nv))
        Lk[:, idxV] = Lkv.T
        Ljk = np.zeros((J, nv))
        Ljk[:, idxV] = Ljkv.T
        lamv = np.zeros(nv)
        lamv[idxV] = lam
        Lall[k] = Lk
        Ljkall[k] = Ljk
        lamall[k] = lamv
    return Lall, Ljkall, lamall


def pfqn_mvac(L: np.ndarray, N: np.ndarray,
              Z: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray,
                                                       np.ndarray, np.ndarray]:
    """
    MVAC (Mean Value Analysis by Chain) for closed product-form networks.

    Exact mean value analysis of a closed multichain product-form queueing
    network composed of single-server fixed-rate (SSFR) queues and
    infinite-server (IS) centers. Unlike the classic MVA recursion of pfqn_mva,
    which recurs on the population vector and costs O(prod(N+1)), MVAC recurs on
    the chains: each class is reduced to single-customer chains and the removed
    chains are replaced by self-looping single-customer (SCSL) chains pinned at
    a service center, so the multiplicity vector v = (v_1,...,v_J), with v_j the
    number of SCSL chains at center j, indexes the recursion in place of the
    population vector. MVAC is thus attractive for networks with few centers and
    many chains, and is the mean-value counterpart of the RECAL
    normalizing-constant recursion of pfqn_recal. Since no normalizing constant
    is formed, MVAC does not suffer the floating-point underflow/overflow that
    complicates RECAL and convolution.

    With j, i indexing centers (j = 1,...,J1 SSFR and j = J1+1,...,J IS), k, l
    indexing single-customer chains, a_jk the relative utilization (service
    demand) of chain k at center j, a_k = sum_j a_jk, K the total number of
    single-customer chains, and I_k = {v : sum_j v_j = K - k}, the recursion of
    Section II of the paper reads

        lambda^k_k(v) = 1 / (a_k + sum_{j=1}^{J1} a_jk (L^{k-1}_j(v) + v_j))  (10)
        L^k_{jk}(v)   = lambda^k_k(v) a_jk (1 + L^{k-1}_j(v) + v_j), j SSFR   (9a)
        L^k_{jk}(v)   = lambda^k_k(v) a_jk,                          j IS     (9b)
        L^k_i(v)      = sum_j L^k_{jk}(v) L^{k-1}_i(v + 1_j) + L^k_{ik}(v)    (7)
        L^k_{il}(v)   = sum_j L^k_{jk}(v) L^{k-1}_{il}(v + 1_j), l = 1..k-1   (6)

    with L^0_j(v) = 0, where L^k_j(v) is the mean number of customers at center
    j (SCSL customers excluded), L^k_{jl}(v) the mean number of chain-l
    customers at center j, and lambda^k_k(v) the throughput of chain k, all for
    the network with normalizing constant G_k(v). Equation (10) is the
    arrival-theorem closure obtained by summing (9a)-(9b) over all centers,
    since chain k holds a single customer. The measures of the original network
    are read off at k = K and v = 0.

    Part 1 of the basic step evaluates (10), (9) and (7) and yields the measures
    of chain K; part 2 evaluates (6) and yields those of the chains that visit
    at least one IS center, whose throughput then follows from Little's law at
    that center. Chains visiting only SSFR centers require a re-execution of
    part 1 with their label interchanged with K, which is cheap because the
    levels below the interchanged label are unaffected and are reused. Classes
    with N_r > 1, and classes with identical demand columns, collapse into a
    single subset of identical single-customer chains: only one representative
    per subset is analyzed and its per-chain measures are scaled by the class
    population, so the cost depends on the number D of distinct chains, not K.

    Args:
        L: Service demand matrix of the SSFR queues (M x R)
        N: Population vector (R,), finite and nonnegative
        Z: Service demand matrix of the IS centers, (R,) or (Mz x R), one row
           per IS center, optional (default: zeros)

    Returns:
        Tuple of (XN, QN, UN, CN):
            XN: Per-class throughput at the reference station (R,)
            QN: Per-class mean queue-length at the SSFR queues (M x R)
            UN: Per-class utilization, XN[r] * L[i,r] (M x R)
            CN: Per-class residence time, QN[i,r] / XN[r] (M x R)

    References:
        A. E. Conway, E. de Souza e Silva and S. S. Lavenberg, "Mean Value
        Analysis by Chain of Product Form Queueing Networks", IEEE Trans.
        Computers, 38(3):432-442, 1989.
        Original MATLAB: matlab/src/api/pfqn/pfqn_mvac.m
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    M, R = L.shape

    if Z is None:
        Z = np.zeros((1, R))
    else:
        Z = np.atleast_2d(np.asarray(Z, dtype=float))
    if Z.shape[1] != R:
        raise ValueError('the think time matrix and the demand matrix have a '
                         'different number of classes')

    N = np.rint(np.asarray(N, dtype=float).flatten()).astype(int)
    if len(N) != R:
        raise ValueError('the demand matrix and the population vector have a '
                         'different number of classes')
    if np.any(N < 0):
        raise ValueError('the population vector must be finite and nonnegative')

    XN = np.zeros(R)
    QN = np.zeros((M, R))
    UN = np.zeros((M, R))
    CN = np.zeros((M, R))

    K = int(np.sum(N))  # number of single-customer chains in the original network
    if K == 0:
        return XN, QN, UN, CN

    # Discard the centers that no chain visits: they carry no customers and
    # would only inflate the multiplicity vector v.
    ssfr_idx = np.flatnonzero(np.any(L > 0, axis=1))
    is_idx = np.flatnonzero(np.any(Z > 0, axis=1))
    J1 = len(ssfr_idx)
    J = J1 + len(is_idx)
    if J == 0:
        raise ValueError('all service demands are zero, the throughput is unbounded')
    # A[j,r] = a_jr, centers ordered SSFR first as required by (9a)-(9b), (10)
    A = np.vstack([L[ssfr_idx, :], Z[is_idx, :]])

    # see _kb/03-api-layer.md for rationale
    posr = np.flatnonzero(N > 0)
    Acols = A[:, posr].T
    Adist, grp_of_class = _unique_rows_stable(Acols)
    Dall = Adist.shape[0]
    visits_is = np.array([np.any(Adist[g, J1:J] > 0) for g in range(Dall)])
    # see _kb/03-api-layer.md for rationale
    gorder = np.concatenate([np.flatnonzero(visits_is), np.flatnonzero(~visits_is)])
    D = Dall
    S = int(np.sum(~visits_is))
    Agrp = Adist[gorder, :].T  # (J x D) demands of each subset representative
    mult = np.array([int(np.sum(N[posr[grp_of_class == gorder[g]]])) for g in range(D)])

    # see _kb/03-api-layer.md for rationale
    chain_group = np.zeros(K, dtype=int)
    chain_group[K - D:K] = np.arange(D)
    p = 0
    for g in range(D):
        for _ in range(mult[g] - 1):
            chain_group[p] = g
            p += 1
    a = Agrp[:, chain_group]  # (J x K) relative utilizations a_jk
    ak = np.sum(a, axis=0)  # (K,) a_k = sum_j a_jk, over ALL centers

    # see _kb/03-api-layer.md for rationale
    Vlist = np.vstack([multichoose(J, t) for t in range(K + 1)])
    nv = Vlist.shape[0]
    vsum = np.sum(Vlist, axis=1)
    # succ[vi,j] is the index of v + 1_j, or 0 when out of range (sum(v) == K),
    # in which case it is never read
    key_of = {tuple(Vlist[i, :]): i for i in range(nv)}
    succ = np.zeros((nv, J), dtype=int)
    for vi in range(nv):
        v = list(Vlist[vi, :])
        for j in range(J):
            v[j] += 1
            succ[vi, j] = key_of.get(tuple(v), 0)
            v[j] -= 1
    z0 = int(np.flatnonzero(vsum == 0)[0])  # index of v = 0

    # MVAC recursion
    Lall = [None] * (K + 1)  # Lall[k][i,vi] = L^k_i(v)
    Ljkall = [None] * (K + 1)  # Ljkall[k][j,vi] = L^k_{jk}(v)
    lamall = [None] * (K + 1)  # lamall[k][vi] = lambda^k_k(v)
    Lall[0] = np.zeros((J, nv))  # L^0_i(v) = 0

    lam_chain = np.zeros(K)  # per-chain throughput, by ORIGINAL chain label
    L_chain = np.zeros((J, K))  # per-chain mean queue-length, same indexing

    # First execution of the basic step, k0 = 1
    Lall, Ljkall, lamall = _part1(1, K, J, J1, a, ak, Vlist, vsum, succ,
                                  Lall, Ljkall, lamall)
    lam_chain[K - 1] = lamall[K][z0]
    L_chain[:, K - 1] = Ljkall[K][:, z0]

    # Part 2: measures of the chains that visit at least one IS center
    lmaxK = min(K - 1, K - S)
    if D >= 2 and lmaxK >= K - D + 1:
        L2prev = [None] * (K + 1)
        for k in range(K - D + 2, K + 1):
            idxI = np.flatnonzero(vsum == K - k)  # v in I_k
            L2cur = [None] * (K + 1)
            for l in range(K - D + 1, min(k - 1, K - S) + 1):
                acc = np.zeros((len(idxI), J))
                for j in range(J):
                    if l == k - 1:
                        # base case of (6): L^{k-1}_{i,k-1} comes from (9)
                        prev = Ljkall[k - 1][:, succ[idxI, j]].T
                    else:
                        prev = L2prev[l][:, succ[idxI, j]].T
                    acc = acc + Ljkall[k][j, idxI][:, None] * prev
                L2cur[l] = np.zeros((J, nv))
                L2cur[l][:, idxI] = acc.T
            if k == K:
                for l in range(K - D + 1, lmaxK + 1):
                    L_chain[:, l - 1] = L2cur[l][:, z0]
                    # Little's law at an IS center visited by chain l
                    jIS = int(np.flatnonzero(a[J1:J, l - 1] > 0)[0]) + J1
                    lam_chain[l - 1] = L_chain[jIS, l - 1] / a[jIS, l - 1]
            L2prev = L2cur

    # see _kb/03-api-layer.md for rationale
    perm = np.arange(K)  # perm[k] is the original label of the chain now at k
    for l in range(1, S):
        perm[[K - l - 1, K - 1]] = perm[[K - 1, K - l - 1]]
        a[:, [K - l - 1, K - 1]] = a[:, [K - 1, K - l - 1]]
        ak[[K - l - 1, K - 1]] = ak[[K - 1, K - l - 1]]
        Lall, Ljkall, lamall = _part1(K - l, K, J, J1, a, ak, Vlist, vsum, succ,
                                      Lall, Ljkall, lamall)
        lam_chain[perm[K - 1]] = lamall[K][z0]
        L_chain[:, perm[K - 1]] = Ljkall[K][:, z0]

    # see _kb/03-api-layer.md for rationale
    for g in range(D):
        kg = K - D + g  # original label of the representative of subset g
        for r in posr[grp_of_class == gorder[g]]:
            XN[r] = N[r] * lam_chain[kg]
            QN[ssfr_idx, r] = N[r] * L_chain[:J1, kg]
            UN[ssfr_idx, r] = XN[r] * L[ssfr_idx, r]
            CN[ssfr_idx, r] = QN[ssfr_idx, r] / XN[r]
    # An empty class has zero throughput and queue-length, so its residence time
    # is reported as the bare service demand, as in pfqn_mva.
    CN[:, N == 0] = L[:, N == 0]

    return XN, QN, UN, CN


def _unique_rows_stable(A: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Distinct rows of A in order of first appearance, with the group index of
    each row (equivalent to MATLAB unique(...,'rows','stable')).
    """
    seen = {}
    order = []
    grp = np.zeros(A.shape[0], dtype=int)
    for i in range(A.shape[0]):
        key = tuple(A[i, :])
        if key not in seen:
            seen[key] = len(order)
            order.append(i)
        grp[i] = seen[key]
    return A[order, :], grp
