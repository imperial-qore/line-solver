"""FPI-based delayed-hit count and expected latency for a delayed-hit cache.

Native-Python port of retrieval_fpi_latency.m / Retrieval_fpi_latency.java
(paper thm:di, eq:moments, eq:latency tot). Supported station types: IS, PS,
SIRO, FCFS, LCFSPR. SIRO/FCFS require exponential service with identical rates.

Representation: ``alpha[s][i]`` is the (fsz_s,) PH entry vector of item i at
station s; ``T[s][i]`` is its (fsz_s, fsz_s) PH subgenerator; ``R[i]`` is the
(S+1, S+1) routing matrix of item i (index 0 = outside/cache, 1..S = stations).
"""
import numpy as np

from .retrieval_fpi import retrieval_fpi


def _ph_mean(al, Tm):
    al = np.asarray(al, dtype=float)
    Tm = np.asarray(Tm, dtype=float)
    z = np.linalg.solve(-Tm, np.ones(Tm.shape[0]))
    return float(al @ z)


def _visits(R, S):
    R = np.asarray(R, dtype=float)
    a = R[0, 1:S + 1]
    P = R[1:S + 1, 1:S + 1]
    return np.linalg.solve(np.eye(S) - P.T, a)


def retrieval_fpi_latency(m, lambda_, gamma, alpha, T, R, station_type):
    """Return (Z, d, phi, pi0)."""
    m = np.asarray(m, dtype=float).ravel()
    lambda_ = np.asarray(lambda_, dtype=float).ravel()
    gamma = np.asarray(gamma, dtype=float)
    n = len(lambda_)
    S = len(station_type)
    fsz = [np.asarray(T[s][0]).shape[0] for s in range(S)]

    for s in range(S):
        st = station_type[s]
        if st not in ("IS", "PS", "SIRO", "FCFS", "LCFSPR"):
            raise RuntimeError("retrieval_fpi_latency supports only IS, PS, SIRO, FCFS and LCFSPR stations; got %s" % st)
        if st in ("SIRO", "FCFS") and fsz[s] > 1:
            raise RuntimeError("retrieval_fpi_latency supports SIRO/FCFS only with exponential (single-phase) service; station %d" % s)
        if st in ("SIRO", "FCFS"):
            tau0 = _ph_mean(alpha[s][0], T[s][0])
            for i in range(1, n):
                if abs(_ph_mean(alpha[s][i], T[s][i]) - tau0) > 1e-9 * max(tau0, 1e-300):
                    raise RuntimeError("retrieval_fpi_latency requires identical per-class rates at SIRO/FCFS station %d" % s)

    is_is = [station_type[s] == "IS" for s in range(S)]
    ps_idx = [s for s in range(S) if not is_is[s]]
    r = len(ps_idx)

    eta = np.zeros((n, r + 1))
    for i in range(n):
        visits = _visits(R[i], S)
        tau = np.array([_ph_mean(alpha[s][i], T[s][i]) for s in range(S)])
        for s in range(S):
            if is_is[s]:
                eta[i, 0] += visits[s] * tau[s]
        for p in range(r):
            eta[i, 1 + p] = visits[ps_idx[p]] * tau[ps_idx[p]]

    pi0, _, pdh_full = retrieval_fpi(m, lambda_, eta, gamma)
    phi = pdh_full.sum(axis=0)

    d = np.zeros(n)
    for i in range(n):
        keep = [k for k in range(n) if k != i]
        _, _, pdh_i = retrieval_fpi(m, lambda_[keep], eta[keep, :], gamma[keep, :])
        phitilde = np.zeros(S)
        for p in range(r):
            phitilde[ps_idx[p]] = pdh_i[1 + p, :].sum()

        Phi = int(sum(fsz))
        off = np.zeros(S, dtype=int)
        for s in range(1, S):
            off[s] = off[s - 1] + fsz[s - 1]
        D0 = np.zeros((Phi, Phi))
        pe = np.zeros(Phi)
        Ri = np.asarray(R[i], dtype=float)
        for s in range(S):
            scale = 1.0 if is_is[s] else 1.0 / (1.0 + phitilde[s])
            blk = scale * np.asarray(T[s][i], dtype=float)
            D0[off[s]:off[s] + fsz[s], off[s]:off[s] + fsz[s]] += blk
            compl = -blk.sum(axis=1)
            al_s = [np.asarray(alpha[sp][i], dtype=float) for sp in range(S)]
            for sp in range(S):
                rprob = Ri[s + 1, sp + 1]
                if rprob == 0:
                    continue
                D0[off[s]:off[s] + fsz[s], off[sp]:off[sp] + fsz[sp]] += np.outer(compl, rprob * al_s[sp])
            pe[off[s]:off[s] + fsz[s]] = Ri[0, s + 1] * al_s[s]

        A = -D0
        e = np.ones(Phi)
        x1 = np.linalg.solve(A, e)
        x2 = np.linalg.solve(A, x1)
        M1 = float(pe @ x1)
        M2 = 2.0 * float(pe @ x2)
        d[i] = phi[i] * lambda_[i] * M2 / (2 * M1)

    num = float(np.sum(phi + d))
    den = float(np.sum(lambda_ * (phi + pi0)))
    return num / den, d, phi, pi0
