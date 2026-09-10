"""
Horizontal-cut mean value analysis for a MAP server (method 'amva.mapqn').

Closed multiclass network of an exponential infinite-server station (think rate
mu_r for class r) and one FCFS single-server station whose class-r service is
the MAP (D0_r, D1_r); the MAP of class r moves only while a class-r job is in
service and is frozen otherwise (the SolverCTMC convention). The recursion
walks the population lattice n <= N in lexicographic order and solves ONE
linear R x R system per point. Its unknowns are the per-phase means
Q_r^k = E[n_r 1{k}] over the joint phase k = (k_1..k_R), the busy laws
U_r^k = P[serving r, k], the phase law pi_k and the throughputs X_r.

Exact relations: the joint phase balance, the class marginals
U_r = X_r E[S_r] theta_r, and the per-class horizontal cut (the generator
balance of n_r 1{k}) of Casale-Smirni, DSN 2009. Closures: the product busy
law theta_r(k_r) prod_{s != r} phi_s(k_s), phi the post-completion law of a
frozen MAP, which solves the phase balance identically; the service-age
closure of the cross term E[n_r 1{serving s} 1{k}] (class r accumulates at its
throughput over the elapsed class-s service, whose mean given the phase is
theta_s (-D0_s)^{-1} / theta_s, the wait before service fixed so the phase
average is the product-form value); Little's law resolved by arrival phase
with the exact FCFS response of the queue composition seen at n - e_r (the
multiclass arrival theorem). K_r = 1 for every class reproduces multiclass
FCFS MVA on class means. Mirrors matlab/src/api/mapqn/mapqn_amva.m.
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class MapqnAmvaResult:
    """X: class throughputs; Qq: mean queue lengths at the MAP station (job in
    service included); U: busy probability of the server per class, X E[S];
    pi: joint phase law at N (class R fastest); ES: mean service times."""
    X: np.ndarray
    Qq: np.ndarray
    U: np.ndarray
    pi: np.ndarray
    ES: np.ndarray


def _stationary(G):
    K = G.shape[0]
    A = np.vstack([G.T, np.ones(K)])
    b = np.zeros(K + 1)
    b[-1] = 1.0
    th, *_ = np.linalg.lstsq(A, b, rcond=None)
    th = np.maximum(th, 0.0)
    return th / th.sum()


def mapqn_amva(mu, D0s, D1s, N) -> MapqnAmvaResult:
    """Horizontal-cut MVA of an exponential delay plus one FCFS MAP queue.

    Args:
        mu: think rates of the R classes at the infinite server
        D0s: list of the R hidden-transition matrices (K_r x K_r)
        D1s: list of the R completion matrices (K_r x K_r)
        N: populations (a class with N_r = 0 is absent)

    Returns:
        MapqnAmvaResult

    Reference:
        G. Casale, E. Smirni, "MAP-AMVA: Approximate Mean Value Analysis of
        Bursty Systems", IEEE/IFIP DSN 2009, pp. 409-418 (the horizontal cut).
    """
    mu = np.asarray(mu, dtype=float).ravel()
    N = np.asarray(np.round(np.asarray(N, dtype=float)), dtype=int).ravel()
    R = N.size
    D0s = [np.atleast_2d(np.asarray(D, dtype=float)) for D in D0s]
    D1s = [np.atleast_2d(np.asarray(D, dtype=float)) for D in D1s]
    Ks = [D0s[r].shape[0] for r in range(R)]
    K = int(np.prod(Ks))
    stride = [int(np.prod(Ks[r + 1:])) for r in range(R)]          # class R-1 fastest
    kr_of = np.array([[(k // stride[r]) % Ks[r] for r in range(R)] for k in range(K)], dtype=int)
    Z = 1.0 / mu
    Ntot = int(N.sum())

    G = [D0s[r] + D1s[r] for r in range(R)]
    th = [_stationary(G[r]) for r in range(R)]
    ES = np.array([1.0 / float(th[r] @ D1s[r].sum(1)) for r in range(R)])
    phi = [th[r] @ D1s[r] * ES[r] for r in range(R)]              # post-completion phase law
    T, age, abar, Ainv = [], [], [], []
    for r in range(R):
        Kr = Ks[r]
        s = np.linalg.solve(-D0s[r], np.ones(Kr))                   # mean time to the next completion
        P = np.linalg.solve(-D0s[r], D1s[r])                        # embedded phase transition
        Tr = np.zeros((Ntot + 2, Kr))                               # T[j,k] = e_k'(I+P+..+P^(j-1)) s
        acc = np.zeros(Kr)
        v = s.copy()
        for j in range(1, Ntot + 2):
            acc = acc + v
            Tr[j] = acc
            v = P @ v
        T.append(Tr)
        w = th[r] @ np.linalg.inv(-D0s[r])
        age.append(w / th[r])                                       # mean elapsed service given the phase
        abar.append(float(w.sum()))
        Ainv.append(np.linalg.inv(G[r] - mu[r] * np.eye(Kr)))

    F = np.ones(K)                                                  # idle law shape
    u = np.ones((R, K))                                             # class-r busy law shape
    for k in range(K):
        for r in range(R):
            F[k] *= phi[r][kr_of[k, r]]
            for s in range(R):
                u[r, k] *= th[s][kr_of[k, s]] if s == r else phi[s][kr_of[k, s]]

    def apply_axis(V, M, r):
        """out[k] = sum_h V[k with k_r -> h] M[h, k_r]"""
        out = np.zeros(K)
        Kr = Ks[r]
        for k in range(K):
            kr = kr_of[k, r]
            base = k - kr * stride[r]
            acc = 0.0
            for h in range(Kr):
                acc += V[base + h * stride[r]] * M[h, kr]
            out[k] = acc
        return out

    def t_at(r, b, k):
        Tr = T[r]
        kr = kr_of[k, r]
        j0 = int(np.floor(b))
        j0 = max(0, min(j0, Tr.shape[0] - 2))
        f = min(max(b - j0, 0.0), 1.0)
        return (1.0 - f) * Tr[j0, kr] + f * Tr[j0 + 1, kr]

    lstride = [int(np.prod([N[s] + 1 for s in range(r + 1, R)])) for r in range(R)]
    L = int(np.prod(N + 1))
    Qs = np.zeros((L, R, K))
    pis = np.zeros((L, K))
    Xs = np.zeros((L, R))
    pis[0] = F
    for l in range(1, L):
        n = np.array([(l // lstride[r]) % (N[r] + 1) for r in range(R)], dtype=int)
        present = [r for r in range(R) if n[r] >= 1]
        # arrival-theorem conditionals at n - e_r, class responses, age closure
        b = {}
        bN = {}
        Rk = np.zeros((R, K))
        for r in present:
            lp = l - lstride[r]
            pip = pis[lp]
            b[r] = np.where(pip > 0, Qs[lp] / np.maximum(pip, 1e-300), 0.0)
            for k in range(K):
                Rk[r, k] = sum(t_at(t, b[r][t, k] + (1.0 if t == r else 0.0), k) for t in range(R))
            bN[r] = np.zeros((R, K))
            for t in range(R):
                if t == r or n[t] == 0:
                    continue
                Xt = Xs[lp, t]
                Qt = Qs[lp, t].sum()
                if Xt > 0:
                    W = max(Qt / Xt - abar[r], 0.0)
                    for k in range(K):
                        bN[r][t, k] = Xt * min(W + age[r][kr_of[k, r]], n[t] / Xt)
        # the cut, linear in X: Q_r = c0[r] + sum_s X_s c1[r][s]
        c0 = {}
        c1 = {}
        for r in present:
            c0[r] = apply_axis(-mu[r] * n[r] * F, Ainv[r], r)
            c1[r] = {}
            for s in range(R):
                term = -mu[r] * n[r] * ES[s] * (u[s] - F)
                if s == r:
                    term = term + ES[r] * apply_axis(u[r], D1s[r], r)
                elif s in present:
                    W = ES[s] * u[s] * bN[s][r]
                    term = term + apply_axis(W, G[r], r) - apply_axis(W, G[s], s)
                c1[r][s] = apply_axis(term, Ainv[r], r)
        # Little's law by arrival phase: one R x R solve
        M = np.eye(R)
        v = np.zeros(R)
        for r in present:
            v[r] = n[r] - mu[r] * float(((n[r] * F - c0[r]) * Rk[r]).sum())
            M[r, r] = Z[r]
            for s in present:
                M[r, s] += mu[r] * float(((n[r] * ES[s] * (u[s] - F) - c1[r][s]) * Rk[r]).sum())
        X = np.linalg.solve(M, v)
        pi = F + sum(X[s] * ES[s] * (u[s] - F) for s in range(R))
        pi = np.maximum(pi, 0.0)
        pi /= pi.sum()
        for r in present:
            Ur = X[r] * ES[r] * u[r]
            Qr = c0[r] + sum(X[s] * c1[r][s] for s in range(R))
            # project onto Q >= U keeping the flow-balance total n_r - X_r/mu_r
            Wr = np.maximum(Qr - Ur, 0.0)
            tot = max(n[r] - X[r] / mu[r] - Ur.sum(), 0.0)
            Qs[l, r] = Ur + (Wr * tot / Wr.sum() if Wr.sum() > 0 else 0.0)
        pis[l] = pi
        Xs[l] = X
    X = Xs[L - 1]
    return MapqnAmvaResult(X=X, Qq=Qs[L - 1].sum(1), U=X * ES, pi=pis[L - 1], ES=ES)
