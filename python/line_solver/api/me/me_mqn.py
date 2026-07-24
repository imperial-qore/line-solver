"""
Maximum Entropy Method for Mixed Queueing Networks.

Extension of the Kouvatsos (1994) Maximum Entropy Method to mixed
open/closed multiclass networks. The 1994 survey notes that the closed
two-stage treatment carries over to mixed networks (Section 2.3) but
gives no algorithm; this implementation composes the open (Section 3.2)
and closed (Section 3.3) ME algorithms by product-form-style conditioning.
"""

from typing import Any, Dict, Optional, Tuple

import numpy as np

from .me_cqn import me_cqn
from .me_oqn import me_oqn


def me_mqn(
    M: int,
    R: int,
    openClasses,
    lambda0,
    Ca0,
    N,
    mu,
    Cs,
    P,
    c: Optional[np.ndarray] = None,
    refstat: Optional[np.ndarray] = None,
    insens: Optional[np.ndarray] = None,
    options: Optional[Dict[str, Any]] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
           np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Maximum Entropy algorithm for Mixed Queueing Networks.

    Composes the open (Section 3.2) and closed (Section 3.3) ME algorithms:
      1. the open classes are solved by the open GE-type fixed point on the
         station set, ignoring the closed classes;
      2. the closed classes are solved by the two-stage pseudo-open plus
         convolution algorithm on servers whose capacity is reduced by the
         open-class utilization, mu_c[i,r] = mu[i,r]*(1-rho_o[i]);
      3. the open mean queue lengths are inflated by the closed occupancy,
         L_o[i,r] *= (1 + Lc[i]) at single-server stations.
    Steps 2-3 are exact in the BCMP product-form limit (where they reduce
    to the classical mixed MVA treatment) and GE-type approximations
    otherwise. Only single-server and infinite-server stations are
    supported.

    Args:
        M: Number of queues (stations).
        R: Number of job classes.
        openClasses: Boolean vector [R], True for open classes.
        lambda0: External arrival rates [M x R], zero columns for closed.
        Ca0: External arrival scv [M x R].
        N: Class populations [R vector], np.inf for open classes.
        mu: Service rates [M x R matrix].
        Cs: Service scv [M x R matrix].
        P: Routing probability matrix [M x M x R], P[j,i,r] = p_ji,r.
        c: Optional servers per queue [M vector]; np.inf marks an IS queue;
           finite values must be 1 (default: all ones).
        refstat: Optional reference station per class [R vector, 0-based].
        options: Optional dict with tol / maxiter / verbose fields.

    Returns:
        Tuple of (L, W, Ca, Cd, lambda_arr, rho, X, iters).

    Reference:
        Kouvatsos (1994), Sections 3.2-3.3.
    """
    if options is None:
        options = {}

    openClasses = np.asarray(openClasses, dtype=bool).ravel()
    lambda0 = np.atleast_2d(np.asarray(lambda0, dtype=float))
    Ca0 = np.atleast_2d(np.asarray(Ca0, dtype=float))
    N = np.asarray(N, dtype=float).ravel()
    mu = np.atleast_2d(np.asarray(mu, dtype=float))
    Cs = np.atleast_2d(np.asarray(Cs, dtype=float))
    P = np.asarray(P, dtype=float)
    if c is None:
        c = np.ones(M)
    else:
        c = np.asarray(c, dtype=float).ravel()
    if refstat is None:
        refstat = -np.ones(R, dtype=int)
    else:
        refstat = np.asarray(refstat, dtype=int).ravel()
    if insens is None:
        insens = np.zeros(M, dtype=bool)
    else:
        insens = np.asarray(insens, dtype=bool).ravel()

    oc = np.where(openClasses)[0]
    cc = np.where(~openClasses)[0]
    Ro = len(oc)
    Rc = len(cc)

    L = np.zeros((M, R))
    W = np.zeros((M, R))
    Ca = np.ones((M, R))
    Cd = np.ones((M, R))
    lam = np.zeros((M, R))
    rho = np.zeros((M, R))
    X = np.zeros(R)
    iters = 0

    # Step 1: open classes by the Section 3.2 GE-type fixed point
    rho_o = np.zeros(M)  # per-station aggregate open utilization
    if Ro > 0:
        Lo, _, Cao, Cdo, lamo, rhoo, itero = me_oqn(
            M, Ro, lambda0[:, oc], Ca0[:, oc], mu[:, oc], Cs[:, oc],
            P[:, :, oc], c, insens, options)
        iters += itero
        L[:, oc] = Lo
        Ca[:, oc] = Cao
        Cd[:, oc] = Cdo
        lam[:, oc] = lamo
        rho[:, oc] = rhoo
        for i in range(M):
            if not np.isinf(c[i]):
                rho_o[i] = np.sum(rhoo[i, :])
        for k in range(Ro):
            X[oc[k]] = np.sum(lambda0[:, oc[k]])

    # Step 2: closed classes by the Section 3.3 algorithm on servers with
    # capacity reduced by the open-class utilization
    if Rc > 0:
        mu_c = mu[:, cc].copy()
        for i in range(M):
            if not np.isinf(c[i]):
                mu_c[i, :] = mu_c[i, :] * max(1.0 - rho_o[i], 0.0)
        Lc, Wc, Cac, Cdc, lamc, rhoc, Xc, iterc = me_cqn(
            M, Rc, N[cc], mu_c, Cs[:, cc], P[:, :, cc], c, refstat[cc], insens, options)
        iters += iterc
        L[:, cc] = Lc
        W[:, cc] = Wc
        Ca[:, cc] = Cac
        Cd[:, cc] = Cdc
        lam[:, cc] = lamc
        X[cc] = Xc
        # Closed utilizations are relative to the reduced capacity; rescale
        # to the busy fraction of the physical server
        for i in range(M):
            for k in range(Rc):
                if np.isinf(c[i]):
                    rho[i, cc[k]] = rhoc[i, k]
                else:
                    rho[i, cc[k]] = rhoc[i, k] * max(1.0 - rho_o[i], 0.0)

    # see _kb/03-api-layer.md for rationale
    if Ro > 0 and Rc > 0:
        for i in range(M):
            if not np.isinf(c[i]):
                Lc_i = np.sum(L[i, cc])
                L[i, oc] = L[i, oc] * (1.0 + Lc_i)
    for i in range(M):
        for k in range(Ro):
            if lam[i, oc[k]] > 0:
                W[i, oc[k]] = L[i, oc[k]] / lam[i, oc[k]]

    return L, W, Ca, Cd, lam, rho, X, iters


__all__ = ['me_mqn']
