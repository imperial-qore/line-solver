"""Whittle balance check for a globally state-dependent rate scaling.

Twin of MATLAB matlab/src/api/sn/sn_gd_balance.m and of the JAR
jline.api.sn.Sn_gd_balance.
"""
import numpy as np

__all__ = ['sn_gd_balance']


def sn_gd_balance(phi, cutoffs):
    """Worst relative violation of the Whittle balance property by phi.

    For every state n of the lattice 0..cutoffs and every pair of stations (s,t)
    populated in n, the property requires

        phi_s(n) phi_t(n - e_s) = phi_t(n) phi_s(n - e_t).

    When it holds, the chain is reversible with pi(n) ~ Phi(n) prod rho**n for
    the balance function Phi implied by phi, and the stationary law is
    insensitive to the service-time distribution beyond its mean. When it fails,
    the model is still solvable by SolverCTMC but has no product form and is
    sensitive.

    phi is evaluated on an (nstations,) population vector, i.e. the single-class
    reading of the (nstations, nclasses) contract of set_global_dependence, and
    must return a scalar or an (nstations,) vector.

    Args:
        phi: the scaling callable
        cutoffs: scalar (same bound at every station) or (nstations,) vector

    Returns:
        (worst relative violation, the state attaining it)

    Reference: P. Whittle, "Partial balance and insensitivity", J. Appl. Prob.
    22(1), 1985; T. Bonald, A. Proutiere, "Insensitivity in processor-sharing
    networks", Perf. Eval. 49, 2002.
    """
    if not callable(phi):
        raise ValueError("phi must be callable.")
    cut = np.atleast_1d(np.asarray(cutoffs, dtype=int)).ravel()
    S = cut.size
    if S < 2:
        raise ValueError(
            "cutoffs must have one entry per station (at least two stations are "
            "needed for a balance pair).")

    def _eval(n):
        v = np.atleast_1d(np.asarray(phi(n), dtype=float)).ravel()
        if v.size == 1:
            v = np.full(S, v[0])
        if v.size != S:
            raise ValueError("phi must return a scalar or a vector of length %d." % S)
        return v

    base = cut + 1
    viol = 0.0
    nworst = np.zeros(S, dtype=int)
    for idx in range(int(np.prod(base))):
        n = np.zeros(S, dtype=int)
        rem = idx
        for s in range(S):
            n[s] = rem % base[s]
            rem //= base[s]
        for s in range(S):
            if n[s] == 0:
                continue
            for t in range(s + 1, S):
                if n[t] == 0:
                    continue
                es = np.zeros(S, dtype=int)
                es[s] = 1
                et = np.zeros(S, dtype=int)
                et[t] = 1
                xn = _eval(n.astype(float))
                xs = _eval((n - es).astype(float))
                xt = _eval((n - et).astype(float))
                lhs = xn[s] * xs[t]
                rhs = xn[t] * xt[s]
                scale = max(abs(lhs), abs(rhs))
                if scale > 0:
                    v = abs(lhs - rhs) / scale
                    if v > viol:
                        viol = v
                        nworst = n.copy()
    return viol, nworst
