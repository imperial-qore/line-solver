"""
Normalizing constant and mean measures for mixed open/closed networks with
limited load dependence.

The closed-conditional normalizing constant of a mixed limited load-dependent (LLD)
network equals a purely closed load-dependent normalizing constant in which every
queueing station i carries the Bruell-Balbo-Afshari effective capacity rate
mu_i^eff(n) = 1/EC_i(n), where EC is returned by pfqn_ldmx_ec and folds the open
classes into the closed subnetwork. The open classes contribute the separable
prefactor lGopen = sum_i log E_i(0), which reduces to -sum_i log(1-rho_i) in the
load-independent limit.

Mean measures follow from that identification without ever enumerating the closed
population lattice, which is what makes this the normalizing-constant counterpart
of pfqn_mvaldmx rather than a rename of it:

  - closed throughputs are the ratios X_r = G(N-e_r)/G(N);
  - closed queue lengths are the conditional normalizing-constant recursion of the
    load-dependent closed network (pfqn_mushift / pfqn_fnc), applied to the
    effective-capacity rates;
  - open queue lengths are the Bruell-Balbo-Afshari sum
    Q_ir = lambda_r D_ir sum_n (n+1) EC_i(n+1) P_i(n) with its SATURATED TAIL
    FOLDED ONTO THE CLOSED MEAN. EC_i(n) is constant for n >= b_i, the level where
    the rate row stops growing, so writing EC_i(n) = EC_i^inf + delta_i(n) with
    delta_i(n) = 0 for n >= b_i leaves

      Q_ir = lambda_r D_ir [ EC_i^inf (Q_i^closed + 1)
                             + sum_{n=0}^{b_i-2} (n+1) delta_i(n+1) P_i(n) ]

    using sum_n P_i(n) = 1 and sum_n n P_i(n) = Q_i^closed. Only the first b_i-1
    marginal probabilities survive, and b_i is the number of servers, not the
    population: a single-server station needs none at all, and the formula
    collapses to the classical lambda_r D_ir (1+Q_i^closed)/(1-rho_i).

The marginals that remain are themselves normalizing-constant ratios,
P_i(n) = sum_{|k|=n} F_i(k) G_{-i}(N-k) / G(N).

Reproduces pfqn_mvaldmx to machine precision on multiserver, nonlinear-mu,
think-time and multi-chain models.

References:
    MATLAB: matlab/src/api/pfqn/pfqn_ncldmx.m
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from scipy.special import gammaln

from .mvaldmx import pfqn_ldmx_ec
from .ncld import pfqn_ncld, pfqn_mushift, pfqn_fnc


@dataclass
class PfqnNcldmxResult:
    """Result of a mixed limited load-dependent normalizing constant computation."""
    G: float
    lG: float
    lGopen: float
    method: str = "default"
    XN: Optional[np.ndarray] = None
    QN: Optional[np.ndarray] = None


def _ncld_lg(L: np.ndarray, N: np.ndarray, Z: np.ndarray, mu: np.ndarray,
             options: Optional[Dict[str, Any]]) -> float:
    """pfqn_ncld, with the empty-station residual network handled explicitly.

    A network reduced to its think times alone has G(N) = prod_r Z_r^N_r / N_r!,
    and no G at all when a class has jobs but neither demand nor think time.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    if L.shape[0] == 0:
        if np.all(N <= 0):
            return 0.0
        Zt = np.sum(np.atleast_2d(np.asarray(Z, dtype=float)), axis=0)
        acc = 0.0
        for r in range(len(N)):
            if N[r] <= 0:
                continue
            if r >= len(Zt) or Zt[r] <= 0:
                return -np.inf
            acc += N[r] * np.log(Zt[r]) - gammaln(N[r] + 1.0)
        return acc
    return float(pfqn_ncld(L, N, Z, mu, options).lG)


def _compositions(n: int, Nc: Sequence[int]) -> List[Tuple[int, ...]]:
    """Non-negative integer vectors k with sum(k) == n and k <= Nc."""
    C = len(Nc)
    if C == 0:
        return [()] if n == 0 else []
    if C == 1:
        return [(n,)] if n <= Nc[0] else []
    out: List[Tuple[int, ...]] = []
    for v in range(0, min(n, int(Nc[0])) + 1):
        for sub in _compositions(n - v, Nc[1:]):
            out.append((v,) + sub)
    return out


def _lld_level(murow: np.ndarray) -> int:
    """First column of the trailing constant run of a limited load-dependence row.

    This is the level b with mu(n) = mu(b) for every n >= b, i.e. the level
    pfqn_ldmx_ec infers and hence the one past which EC is constant.
    """
    b = len(murow)
    if b == 0:
        return 1
    while b > 1 and murow[b - 2] == murow[b - 1]:
        b -= 1
    return b


def _marginal(n: int, ist: int, Dc: np.ndarray, Nc: np.ndarray, Zc: np.ndarray,
              muEff: np.ndarray, Dminus: np.ndarray, muminus: np.ndarray,
              lG: float, options: Optional[Dict[str, Any]]) -> float:
    """P_i(n) = sum_{|k|=n, k<=Nc} F_i(k) G_{-i}(Nc-k) / G(Nc)."""
    Cc = len(Nc)
    P = 0.0
    for k in _compositions(n, [int(x) for x in Nc]):
        kk = np.asarray(k, dtype=float)
        lF = 0.0 if n == 0 else _ncld_lg(Dc[ist:ist + 1, :], kk, np.zeros(Cc),
                                         muEff[ist:ist + 1, :], options)
        lGbar = _ncld_lg(Dminus, Nc - kk, Zc, muminus, options)
        P += np.exp(lF + lGbar - lG)
    return float(P)


def pfqn_ncldmx(lam: np.ndarray, D: np.ndarray, N: np.ndarray,
                Z: Optional[np.ndarray] = None,
                mu: Optional[np.ndarray] = None,
                S: Optional[np.ndarray] = None,
                options: Optional[Dict[str, Any]] = None) -> PfqnNcldmxResult:
    """
    Normalizing constant and mean measures for mixed open/closed networks with
    limited load dependence.

    Args:
        lam: Arrival rate vector (R,) - 0 on closed classes
        D: Service demand matrix (M x R)
        N: Population vector (R,) - inf for open classes
        Z: Think time vector (R,), optional
        mu: Load-dependent rate matrix (M x >= sum(N_closed)), optional
        S: Number of servers per station (M,), kept for signature parity
        options: Solver options forwarded to pfqn_ncld

    Returns:
        PfqnNcldmxResult with the closed-conditional constant (G, lG), the open-class
        normalizing prefactor lGopen, the method used for the closed-conditional
        solve, and the mean throughputs XN (1 x R) and queue lengths QN (M x R).
    """
    D = np.atleast_2d(np.asarray(D, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    lam = np.asarray(lam, dtype=float).flatten()
    M, R = D.shape
    if Z is None:
        Z = np.zeros(R)
    Z = np.asarray(Z, dtype=float).flatten()

    openClasses = np.where(np.isinf(N))[0]
    closedClasses = np.array([r for r in range(R) if r not in openClasses], dtype=int)
    for r in closedClasses:
        if lam[r] != 0 and N[r] > 0:
            raise ValueError("pfqn_ncldmx: Arrival rate cannot be specified on closed classes.")

    Kc = int(np.sum(N[closedClasses])) if closedClasses.size > 0 else 0

    if mu is None:
        mu = np.ones((M, max(1, Kc)))
    mu = np.atleast_2d(np.asarray(mu, dtype=float))

    # pad mu to at least max(1,Kc) columns, then append one extra column as in
    # pfqn_mvaldmx so pfqn_ldmx_ec returns EC with Nt = max(1,Kc)+1 columns.
    min_cols = max(1, Kc)
    pad_cols = max(mu.shape[1], min_cols) + 1
    mup = np.empty((M, pad_cols))
    ncol_mu = mu.shape[1]
    for j in range(pad_cols):
        mup[:, j] = mu[:, min(j, ncol_mu - 1)]

    lamo = np.zeros(R)
    lamo[openClasses] = lam[openClasses]
    EC, E, _Eprime, _Lo = pfqn_ldmx_ec(lamo, D, mup)
    lGopen = float(np.sum(np.log(E[:, 0])))

    Dc = D[:, closedClasses] if closedClasses.size > 0 else np.zeros((M, 0))
    Nc = N[closedClasses] if closedClasses.size > 0 else np.zeros(0)
    Zc = Z[closedClasses] if closedClasses.size > 0 else np.zeros(0)
    muEff = 1.0 / EC[:, :max(1, Kc)]

    if Kc == 0:
        lG, G, method = 0.0, 1.0, "exact"
    else:
        cc = pfqn_ncld(Dc, Nc, Zc, muEff, options)
        lG, G, method = float(cc.lG), cc.G, cc.method

    # ---- mean measures ----
    XN = np.zeros(R)
    QN = np.zeros((M, R))
    for r in openClasses:
        XN[r] = lam[r]

    Cc = int(closedClasses.size)
    lGr = np.zeros(max(Cc, 1))
    if Kc > 0:
        for rc in range(Cc):
            if Nc[rc] <= 0:
                continue
            Ncr = Nc.copy()
            Ncr[rc] -= 1
            lGr[rc] = _ncld_lg(Dc, Ncr, Zc, muEff, options)
            XN[closedClasses[rc]] = np.exp(lGr[rc] - lG)
        # closed queue lengths: conditional normalizing-constant recursion of the
        # load-dependent closed network, on the effective-capacity rates
        for ist in range(M):
            if not np.any(Dc[ist, :] > 0):
                continue
            muhat = pfqn_mushift(muEff, ist)
            fres = pfqn_fnc(muhat[ist, :])
            muhat_f = np.atleast_2d(fres.mu)
            cshift = float(np.asarray(fres.c).flatten()[0])
            Dminus = np.delete(Dc, ist, axis=0)
            muminus = np.delete(muEff, ist, axis=0)
            for rc in range(Cc):
                if Nc[rc] <= 0 or Dc[ist, rc] <= 0:
                    continue
                Ncr = Nc.copy()
                Ncr[rc] -= 1
                lGhat = _ncld_lg(Dc, Ncr, Zc, muhat, options)
                lGhatf = _ncld_lg(np.vstack([Dc, Dc[ist:ist + 1, :]]), Ncr, Zc,
                                  np.vstack([muhat, muhat_f]), options)
                lGminus = _ncld_lg(Dminus, Ncr, Zc, muminus, options)
                CQ = (np.exp(lGhatf - lGhat) - 1.0) + cshift * (np.exp(lGminus - lGhat) - 1.0)
                ld = np.log(Dc[ist, rc]) + lGhat - np.log(muEff[ist, 0]) - lGr[rc]
                QN[ist, closedClasses[rc]] = np.exp(ld) * XN[closedClasses[rc]] * (1.0 + CQ)

    # open queue lengths, with the saturated tail of EC folded onto the closed mean
    if openClasses.size > 0:
        Qtot = np.sum(QN[:, closedClasses], axis=1) if Cc > 0 else np.zeros(M)
        for ist in range(M):
            b = _lld_level(mup[ist, :])
            ECinf = EC[ist, min(b, EC.shape[1]) - 1]
            acc = ECinf * (Qtot[ist] + 1.0)
            if b >= 2:
                if Kc > 0:
                    Dminus = np.delete(Dc, ist, axis=0)
                    muminus = np.delete(muEff, ist, axis=0)
                for n in range(0, b - 1):
                    delta = EC[ist, n] - ECinf   # EC_i(n+1), 1-based
                    if delta == 0.0:
                        continue
                    # WITH NO CLOSED POPULATION THE MARGINAL IS DEGENERATE, not
                    # absent: P_i(0)=1 and P_i(n)=0 above it, so only the n=0
                    # term survives and acc collapses to EC_i(1), which is the
                    # exact open load-dependent mean. Skipping the loop instead
                    # left acc at EC_i^inf, i.e. read a c-server station as if
                    # every arrival found it saturated -- an M/M/3 at lambda=1.5
                    # came back with mean 1 against the exact 1.7368.
                    if Kc == 0:
                        Pn = 1.0 if n == 0 else 0.0
                    else:
                        Pn = _marginal(n, ist, Dc, Nc, Zc, muEff,
                                       Dminus, muminus, lG, options)
                    acc += (n + 1) * delta * Pn
            for r in openClasses:
                QN[ist, r] = lam[r] * D[ist, r] * acc

    return PfqnNcldmxResult(G=G, lG=lG, lGopen=lGopen, method=method, XN=XN, QN=QN)
