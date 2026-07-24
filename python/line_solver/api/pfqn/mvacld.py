"""
MVAC (Mean Value Analysis by Chain) for load-dependent closed product-form networks.

Native Python implementation of the Section V extension of MVAC to queue-length
dependent (QLD) service centers, which propagates the marginal queue-length
distributions instead of the mean queue-lengths of pfqn_mvac.

Key functions:
    pfqn_mvacld: MVAC for closed multichain networks with QLD centers

References:
    A. E. Conway, E. de Souza e Silva and S. S. Lavenberg, "Mean Value Analysis
    by Chain of Product Form Queueing Networks", IEEE Trans. Computers,
    38(3):432-442, 1989, Section V.
    Original MATLAB: matlab/src/api/pfqn/pfqn_mvacld.m
"""

import numpy as np
from typing import Optional, Tuple

from .utils import multichoose
from .mvac import _unique_rows_stable

__all__ = ['pfqn_mvacld']


def _part1(k0, K, J, J1, a, MU, Vlist, off, cnt, succ, Pall, Ljkall, lamall):
    """
    Part 1 of the basic step for networks with QLD centers: evaluate (21)-(25)
    for k = k0,...,K over v in I_k, the block of multiplicity vectors of sum K-k.
    """
    for k in range(k0, K + 1):
        t = K - k
        nvk = cnt[t]
        Pp = Pall[k - 1]  # P^{k-1}, over I_{k-1}, last axis n = 0,...,k-1
        Lk = np.zeros((J, nvk))
        lam = np.zeros(nvk)
        Pk = np.zeros((J1, nvk, k + 1))
        for vloc in range(nvk):
            vi = off[t] + vloc
            v = Vlist[vi, :]
            sloc = succ[vi, :] - off[t + 1]  # local index of v + 1_i within I_{k-1}
            # see _kb/03-api-layer.md for rationale
            c = np.ones(J)
            for i in range(J1):
                nn = np.arange(k)
                c[i] = np.sum(Pp[i, sloc[i], :k] * MU[i, nn + v[i]] / (nn + v[i] + 1))
            # (23)-(24) in reference-station-free form:
            # theta_jk/tau_k(v,j) = a_jk/c_j
            w = np.zeros(J)
            nz = a[:, k - 1] > 0
            w[nz] = a[nz, k - 1] / c[nz]
            sw = w.sum()
            lam[vloc] = 1.0 / sw
            Lk[:, vloc] = w / sw
            # (25): condition on the center holding the single chain-k customer
            for j in range(J1):
                for n in range(k + 1):
                    s = 0.0
                    if n >= 1:
                        # chain k is at j, so n-1 of the k-1 others are there too
                        s = Lk[j, vloc] * Pp[j, sloc[j], n - 1]
                    if n <= k - 1:
                        # chain k is elsewhere, so all n are from the k-1 others
                        for m in range(J):
                            if m != j:
                                s += Lk[m, vloc] * Pp[j, sloc[m], n]
                    Pk[j, vloc, n] = s
        Pall[k] = Pk
        Ljkall[k] = Lk
        lamall[k] = lam
    return Pall, Ljkall, lamall


def pfqn_mvacld(L: np.ndarray, N: np.ndarray,
                Z: Optional[np.ndarray] = None,
                mu: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray,
                                                          np.ndarray, np.ndarray,
                                                          np.ndarray]:
    """
    MVAC for closed product-form networks with queue-length dependent centers.

    Exact mean value analysis by chain of a closed multichain product-form
    queueing network that may contain queue-length dependent (QLD) service
    centers. This is the Section V extension of Conway, de Souza e Silva and
    Lavenberg (1989); pfqn_mvac implements Sections II-IV, which cover
    single-server fixed-rate (SSFR) and infinite-server (IS) centers only.

    Where pfqn_mvac propagates the MEAN queue-lengths L^k_j(v) through eq. (7)
    and closes the recursion with the arrival-theorem identity (10), the QLD
    extension propagates the MARGINAL queue-length DISTRIBUTIONS P^k_j(n,v)
    instead. That is forced by load dependence -- the rate seen by a job depends
    on the whole occupancy, so a mean no longer suffices -- but it also
    SIMPLIFIES the recursion: eq. (21)-(25) read level k-1 only at the shifted
    vectors v + 1_i, so the basic step sweeps v in I_k alone, where pfqn_mvac
    must sweep the larger I_k u ... u I_K. The marginals come almost for free and
    are returned as a first-class output.

    Notation follows pfqn_mvac: j, i index centers (j = 1,...,J1 the QLD centers
    of L, j = J1+1,...,J the IS centers of Z), k, l index the single-customer
    chains, a_jk = theta_jk T_jk is the demand of chain k at center j, K = sum(N)
    and I_k = {v : sum_j v_j = K - k} with v_j the number of self-looping
    single-customer (SCSL) chains pinned at center j. P^k_j(n,v) is the
    probability of n customers at center j -- EXCLUDING the v_j SCSL customers
    there -- in the network with normalizing constant G_k(v). Writing tau_k(v,i)
    for the throughput of an SCSL chain that replaces chain k at center i:

        tau_k(v,i)  = T_ik^-1 sum_{n=0}^{k-1} P^{k-1}_i(n,v+1_i)
                                 * mu_i(n+v_i+1)/(n+v_i+1)                   (21)
        tau_k(v,i)  = T_ik^-1,                                   i IS        (22)
        L^k_{jk}(v) = theta_jk tau_k(v,j)^-1 / sum_m theta_mk tau_k(v,m)^-1  (23)
        lambda^k_k(v) = tau_k(v,j(k)) L^k_{j(k)k}(v)                         (24)
        P^k_j(n,v)  = L^k_{jk}(v) P^{k-1}_j(n-1,v+1_j)
                      + sum_{m != j} L^k_{mk}(v) P^{k-1}_j(n,v+1_m)          (25)

    with P^0_j(0,v) = 1 and P^{k-1}_j(n,.) = 0 for n < 0 or n > k-1. Eq. (21) is
    just "the mean rate at which an SCSL chain is served": given n other
    customers the processor-sharing rate share is mu_i(n+v_i+1)/(n+v_i+1),
    averaged over the distribution of those others. The queueing discipline may
    be assumed PS with no loss of generality, since product-form measures do not
    depend on it.

    This implementation writes (23)-(24) in the reference-station-free form

        c_i(k,v)      = sum_{n=0}^{k-1} P^{k-1}_i(n,v+1_i)
                            * mu_i(n+v_i+1)/(n+v_i+1)
        L^k_{jk}(v)   = (a_jk/c_j) / sum_m (a_mk/c_m)
        lambda^k_k(v) = 1 / sum_m (a_mk/c_m)

    which follows from theta_jk tau_k(v,j)^-1 = a_jk/c_j and theta_{j(k)k} = 1,
    so only the demands a_jk are needed and the visit ratios never appear
    separately. For an IS center c_i = 1 identically, which is exactly (22).
    Eq. (25) is self-normalizing, sum_n P^k_j(n,v) = sum_m L^k_{mk}(v) = 1, so no
    normalizing constant is formed and the recursion involves only positive
    quantities: unlike the classic load-dependent MVA of pfqn_mvald it cannot
    produce negative probabilities and needs no stabilization.

    Parts 2 and 3 are unchanged from pfqn_mvac, since eq. (6) holds verbatim in
    the presence of QLD centers: part 2 resolves the chains that visit at least
    one IS center, and the chains that visit no IS center are resolved by
    re-executing part 1 with their label interchanged with K.

    Args:
        L: Service demand matrix of the QLD centers (M x R)
        N: Population vector (R,), finite and nonnegative
        Z: Demand matrix of the IS centers, (R,) or (Mz x R), one row per
           center, optional (default: zeros)
        mu: Load-dependent rates (M x Nt) with Nt >= sum(N); mu[j,n-1] is the
            total service rate of center j with n jobs present. mu[j,:] = 1 is a
            single-server fixed-rate queue, mu[j,n-1] = min(n,c) a c-server
            queue, mu[j,n-1] = n an infinite server. Optional (default: ones,
            i.e. all centers SSFR, in which case results agree with pfqn_mvac)

    Returns:
        Tuple of (XN, QN, UN, CN, pij):
            XN: Per-class throughput at the reference station (R,)
            QN: Per-class mean queue-length at the QLD centers (M x R)
            UN: Utilization of each center (M,), 1 - P_j(0). PER-STATION, not
                per-class, as in pfqn_mvald and pfqn_dac: for a load-dependent
                center the per-class product XN[r]*L[j,r] of pfqn_mvac is NOT
                the utilization
            CN: Per-class cycle time exclusive of think time (R,),
                N[r]/XN[r]-Z[r], as in pfqn_mvald and pfqn_dac. NOT the (M x R)
                per-station residence time of pfqn_mvac: the whole LD family
                reports a cycle time here
            pij: Marginal queue-length probabilities (M x (sum(N)+1)),
                 pij[j,n] = P(n jobs at center j)

    References:
        A. E. Conway, E. de Souza e Silva and S. S. Lavenberg, "Mean Value
        Analysis by Chain of Product Form Queueing Networks", IEEE Trans.
        Computers, 38(3):432-442, 1989, Section V.
        Original MATLAB: matlab/src/api/pfqn/pfqn_mvacld.m
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

    K = int(np.sum(N))  # number of single-customer chains in the original network
    if mu is None:
        mu = np.ones((M, max(K, 1)))
    else:
        mu = np.atleast_2d(np.asarray(mu, dtype=float))
    if mu.shape[0] != M:
        raise ValueError('the rate matrix and the demand matrix have a '
                         'different number of centers')
    if K > 0 and mu.shape[1] < K:
        raise ValueError('the rate matrix must supply a rate for every '
                         'population up to sum(N)')

    XN = np.zeros(R)
    QN = np.zeros((M, R))
    UN = np.zeros(M)
    CN = np.zeros(R)
    pij = np.zeros((M, K + 1))
    pij[:, 0] = 1.0  # an unvisited center holds no jobs
    if K == 0:
        return XN, QN, UN, CN, pij

    # Discard the centers that no chain visits: they carry no customers and
    # would only inflate the multiplicity vector v.
    ld_idx = np.flatnonzero(np.any(L > 0, axis=1))
    is_idx = np.flatnonzero(np.any(Z > 0, axis=1))
    J1 = len(ld_idx)
    J = J1 + len(is_idx)
    if J == 0:
        raise ValueError('all service demands are zero, the throughput is unbounded')
    # A[j,r] = a_jr, centers ordered QLD first, then IS, as required by (21)-(22)
    A = np.vstack([L[ld_idx, :], Z[is_idx, :]])
    MU = mu[ld_idx, :K]  # (J1 x K) rates of the QLD centers
    if np.any(MU <= 0):
        raise ValueError('the service rates must be strictly positive for every '
                         'population up to sum(N)')

    # Partition the chains into subsets of identical single-customer chains.
    posr = np.flatnonzero(N > 0)
    Adist, grp_of_class = _unique_rows_stable(A[:, posr].T)
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

    # see _kb/03-api-layer.md for rationale
    blocks = []
    cnt = np.zeros(K + 1, dtype=int)
    off = np.zeros(K + 1, dtype=int)
    total = 0
    for t in range(K + 1):
        Vt = multichoose(J, t)
        off[t] = total
        cnt[t] = Vt.shape[0]
        total += Vt.shape[0]
        blocks.append(Vt)
    Vlist = np.vstack(blocks)
    nv = Vlist.shape[0]
    # succ[vi,i] is the global index of v + 1_i, or 0 when out of range
    # (sum(v) == K), in which case it is never read
    key_of = {tuple(Vlist[i, :]): i for i in range(nv)}
    succ = np.zeros((nv, J), dtype=int)
    for vi in range(nv):
        v = list(Vlist[vi, :])
        for j in range(J):
            v[j] += 1
            succ[vi, j] = key_of.get(tuple(v), 0)
            v[j] -= 1

    # MVAC recursion
    Pall = [None] * (K + 1)  # Pall[k][j,vloc,n] = P^k_j(n,v), v the vloc-th of I_k
    Ljkall = [None] * (K + 1)  # Ljkall[k][j,vloc] = L^k_{jk}(v)
    lamall = [None] * (K + 1)  # lamall[k][vloc] = lambda^k_k(v)
    Pall[0] = np.ones((J1, cnt[K], 1))  # P^0_j(0,v) = 1, I_0 is the block of sum K

    lam_chain = np.zeros(K)  # per-chain throughput, by ORIGINAL chain label
    L_chain = np.zeros((J, K))  # per-chain mean queue-length, same indexing

    # First execution of the basic step, k0 = 1
    Pall, Ljkall, lamall = _part1(1, K, J, J1, a, MU, Vlist, off, cnt, succ,
                                  Pall, Ljkall, lamall)
    lam_chain[K - 1] = lamall[K][0]
    L_chain[:, K - 1] = Ljkall[K][:, 0]
    # The marginals of the ORIGINAL network are read at k = K and v = 0; the
    # label interchanges below overwrite this level, so capture them now.
    pij[ld_idx, :] = Pall[K][:, 0, :]

    # Part 2: measures of the chains that visit at least one IS center. Eq. (6)
    # is unchanged by load dependence.
    lmaxK = min(K - 1, K - S)
    if D >= 2 and lmaxK >= K - D + 1:
        L2prev = [None] * (K + 1)
        for k in range(K - D + 2, K + 1):
            t = K - k
            L2cur = [None] * (K + 1)
            for l in range(K - D + 1, min(k - 1, K - S) + 1):
                acc = np.zeros((J, cnt[t]))
                for vloc in range(cnt[t]):
                    vi = off[t] + vloc
                    for j in range(J):
                        sloc = succ[vi, j] - off[t + 1]
                        if l == k - 1:
                            # base case of (6): L^{k-1}_{i,k-1} comes from (23)
                            prev = Ljkall[k - 1][:, sloc]
                        else:
                            prev = L2prev[l][:, sloc]
                        acc[:, vloc] += Ljkall[k][j, vloc] * prev
                L2cur[l] = acc
            if k == K:
                for l in range(K - D + 1, lmaxK + 1):
                    L_chain[:, l - 1] = L2cur[l][:, 0]
                    # Little's law at an IS center visited by chain l
                    jIS = int(np.flatnonzero(a[J1:J, l - 1] > 0)[0]) + J1
                    lam_chain[l - 1] = L_chain[jIS, l - 1] / a[jIS, l - 1]
            L2prev = L2cur

    # see _kb/03-api-layer.md for rationale
    perm = np.arange(K)  # perm[k] is the original label of the chain now at k
    for l in range(1, S):
        perm[[K - l - 1, K - 1]] = perm[[K - 1, K - l - 1]]
        a[:, [K - l - 1, K - 1]] = a[:, [K - 1, K - l - 1]]
        Pall, Ljkall, lamall = _part1(K - l, K, J, J1, a, MU, Vlist, off, cnt,
                                      succ, Pall, Ljkall, lamall)
        lam_chain[perm[K - 1]] = lamall[K][0]
        L_chain[:, perm[K - 1]] = Ljkall[K][:, 0]

    # see _kb/03-api-layer.md for rationale
    for g in range(D):
        kg = K - D + g  # original label of the representative of subset g
        for r in posr[grp_of_class == gorder[g]]:
            XN[r] = N[r] * lam_chain[kg]
            QN[ld_idx, r] = N[r] * L_chain[:J1, kg]
    # Utilization of a load-dependent center is 1-P_j(0), NOT XN[r]*L[j,r].
    UN = 1.0 - pij[:, 0]
    # Cycle time exclusive of think time, as in pfqn_mvald and pfqn_dac. An empty
    # class has zero throughput, hence no cycle time.
    CN[posr] = N[posr] / XN[posr] - Z[:, posr].sum(axis=0)

    return XN, QN, UN, CN, pij
