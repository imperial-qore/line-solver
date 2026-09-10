"""
Lattice-Poisson inversion of the multichain generating function for stations
whose rate reads the whole per-class occupancy vector.

Native Python port of the MATLAB routines:
    pfqn_clwoi - order-independent (OI) stations, i.e. rates that depend on the
                 occupancy only through its support. Transform counterpart of
                 the convolution routine pfqn_ncoi.
    pfqn_clwjd - limited joint-dependent (LJD) stations, i.e. rates that
                 saturate coordinatewise past a cutoff vector. Generalization
                 of pfqn_clwoi, which is the case lcut = 1.

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_clwoi.m, pfqn_clwjd.m
    G. L. Choudhury, K. K. Leung and W. Whitt, "Calculating normalization
    constants of closed queueing networks by numerically inverting their
    generating functions", J. ACM 42(5):935-970, 1995.
    T. Bonald and A. Proutiere, "Insensitive bandwidth sharing in data
    networks", Queueing Systems 44, 2003.
"""

import numpy as np
from typing import Callable, List, Optional, Sequence, Tuple, Union


def _clw_default_lattice(options, R: int) -> Tuple[np.ndarray, np.ndarray]:
    """Inner lattice and aliasing parameters (CLW Section 2.2, page 962)."""
    lpar = None
    gam = None
    if options is not None:
        if isinstance(options, dict):
            lpar = options.get('l', None)
            gam = options.get('gamma', None)
        else:
            lpar = getattr(options, 'l', None)
            gam = getattr(options, 'gamma', None)
    if lpar is None or len(np.atleast_1d(lpar)) == 0:
        lpar = np.full(R, 3)
        lpar[0] = 1
        if R >= 2:
            lpar[1] = 2
        if R >= 3:
            lpar[2] = 2
    lpar = np.round(np.asarray(lpar, dtype=float).ravel()).astype(int)
    if gam is None or len(np.atleast_1d(gam)) == 0:
        gam = np.full(R, 15.0)
        gam[0] = 11.0
        if R >= 2:
            gam[1] = 13.0
        if R >= 3:
            gam[2] = 13.0
    gam = np.asarray(gam, dtype=float).ravel()
    return lpar, gam


def _clw_prepare(who: str, Z, N, mu, visits) -> tuple:
    """Shared argument normalization of the two inversion routines."""
    if mu is None:
        mu = []
    elif callable(mu):
        mu = [mu]
    else:
        mu = list(mu)
    M = len(mu)

    N = np.asarray(N, dtype=float).ravel()
    R = N.size
    if np.any(~np.isfinite(N)):
        raise ValueError('%s requires finite (closed) populations.' % who)
    N = np.round(N).astype(int)

    if Z is None or np.asarray(Z).size == 0:
        Z = np.zeros(R)
    Z = np.asarray(Z, dtype=float).ravel()
    if Z.size != R:
        raise ValueError('%s: Z and N must have the same number of classes.' % who)

    if visits is None:
        vis = [np.ones(R) for _ in range(M)]
    elif isinstance(visits, np.ndarray) and visits.ndim == 2:
        vis = [np.asarray(visits[m, :], dtype=float).ravel() for m in range(M)]
    else:
        vis = [np.asarray(visits[m], dtype=float).ravel() for m in range(M)]
    for vm in vis:
        if vm.size != R:
            raise ValueError('%s: each visit vector must have one entry per class.' % who)
    return mu, M, N, R, Z, vis


def _clw_scaling(Lt: np.ndarray, Nk: np.ndarray, lpar: np.ndarray,
                 rad: np.ndarray, Zk: np.ndarray) -> np.ndarray:
    """Restrictive static scaling of CLW eqs. 5.41-5.46 with unit multiplicities.

    Lt lists one row per singular hyperplane, holding its unit-pole intensities.
    """
    nrow, Rk = Lt.shape
    alpha = np.ones(Rk)
    used = np.zeros(nrow)
    etaMat = (Lt != 0).astype(float)
    for j in range(Rk):
        Kj = int(Nk[j])
        lj = int(lpar[j])
        denom = 1.0 - used
        denom[denom <= 0] = np.finfo(float).eps
        e = Lt[:, j] / denom
        posq = np.where(Lt[:, j] > 0)[0]
        aj = np.inf
        if posq.size > 0:
            order = np.argsort(-e[posq], kind='stable')
            qs = posq[order]
            es = e[qs]
            cumrho = np.cumsum(es) / np.arange(1, es.size + 1)
            for n in range(es.size):
                qi = qs[n]
                # N_{ij} = n - 1 + sum_{k>j} K_k eta_{k,qi} (eq. 5.43, m_i = 1)
                Nn = int(round(n + np.sum(Nk[j + 1:Rk] * etaMat[qi, j + 1:Rk])))
                if Nn <= 0:
                    an = 1.0
                else:
                    # in the log domain: the product runs over N_{ij} factors
                    # below one and underflows to zero at a few hundred of them,
                    # which would silently set alpha_j = 0 and lG = NaN
                    ll = np.arange(1, Nn + 1)
                    an = np.exp(np.sum(np.log((Kj + ll) / (Kj + 2 * lj * Kj + ll)))
                                / (2 * lj * Kj))
                aj = min(aj, an / cumrho[n])
        if Zk[j] > 0:
            aj = min(aj, Kj / Zk[j])          # IS/Poisson term K_j/rho_{j0}
        if not np.isfinite(aj):
            aj = 1.0                          # chain with no demand anywhere
        alpha[j] = aj
        used = used + aj * Lt[:, j] * rad[j]
    return alpha


def _clw_invert(j: int, wfixed: np.ndarray, ctx: dict, gbar) -> complex:
    """One-dimensional lattice-Poisson inversion (CLW eq. 2.3), scaled.

    Extracts the coefficient of w_j^{N_j} from g^(j), recursing on inner chains.
    """
    Kj = int(ctx['N'][j])
    lj = int(ctx['l'][j])
    rj = float(ctx['r'][j])
    p = ctx['p']
    kk = np.arange(-Kj, Kj)
    signs = (-1.0) ** kk
    acc = 0.0 + 0.0j
    for k1 in range(lj):
        ph = np.exp(-1j * np.pi * k1 / lj)
        theta = np.pi * (k1 + lj * kk) / (lj * Kj)
        wj = rj * np.exp(1j * theta)          # 2Kj contour points
        inner = 0.0 + 0.0j
        if j == p - 1:
            nk = wj.size
            chunk = ctx['chunk']
            for a in range(0, nk, chunk):
                b = min(a + chunk, nk)
                W = np.empty((b - a, p), dtype=complex)
                if j > 0:
                    W[:, :j] = wfixed
                W[:, j] = wj[a:b]
                inner += np.sum(signs[a:b] * gbar(W, ctx))
        else:
            for t in range(wj.size):
                inner += signs[t] * _clw_invert(j + 1, np.concatenate([wfixed, [wj[t]]]), ctx, gbar)
        acc += ph * inner
    val = acc / (2 * lj * Kj * rj ** Kj)
    if j == 0:
        val = val.real
    return val


def _clw_recover(gbarN, arho0: np.ndarray, Nk: np.ndarray,
                 alpha: np.ndarray) -> Tuple[float, float]:
    """G(N) = exp(sum alpha_r Z_r) prod alpha_r^{-N_r} gbar(N) (CLW eq. 7.1)."""
    lG = float(np.log(gbarN) + np.sum(arho0) - np.sum(Nk * np.log(alpha)))
    G = np.inf if lG > 709 else float(np.exp(lG))
    return G, lG


def _clwoi_supportrate(murate: Callable, chi: np.ndarray, m: int) -> float:
    """Support rate of an OI station, read at the support indicator.

    chi is the 0/1 indicator of the support, itself a lattice point of that
    support since every retained chain has N_r >= 1.
    """
    rate = float(murate(chi))
    if not (rate > 0):
        raise ValueError('pfqn_clwoi: station %d has a non-positive rate on a reachable support.'
                         % (m + 1))
    return rate


def _clwoi_checksupport(mu: List[Callable], muS: np.ndarray, keep: np.ndarray,
                        Nk: np.ndarray, R: int) -> None:
    """Exhaustive support-only check over the count lattice 0 < n <= N.

    The scan costs prod_r (N_r+1) rate evaluations per station, below the
    prod_r 2 l_r N_r contour points the inversion itself spends.
    """
    M = len(mu)
    if M == 0:
        return
    Rk = keep.size
    L = int(np.prod(Nk + 1))
    for m in range(M):
        n = np.zeros(R)
        for idx in range(1, L):           # idx 0 is the empty support, unused
            rem = idx
            mask = 0
            n[:] = 0.0
            for j in range(Rk):
                nj = rem % (int(Nk[j]) + 1)
                rem //= (int(Nk[j]) + 1)
                n[keep[j]] = nj
                if nj > 0:
                    mask += 1 << j
            rate = float(mu[m](n))
            ref = muS[m, mask]
            if abs(rate - ref) > 1e-9 * max(1.0, abs(ref)):
                chi = np.zeros(R)
                chi[keep[[b for b in range(Rk) if (mask >> b) & 1]]] = 1.0
                raise ValueError(
                    'pfqn_clwoi: station %d has a rate that varies within a support: '
                    'mu=%g at n=%s but mu=%g at the indicator %s of the same support. '
                    'pfqn_clwoi requires order-independent (support-only) rates, '
                    'mu(n)=mu(supp(n)); a rate that varies inside a support is a general '
                    'balanced-fairness station and must be solved with pfqn_ncoi.'
                    % (m + 1, rate, np.array2string(n), ref, np.array2string(chi)))


def _clwoi_gbar(W: np.ndarray, ctx: dict) -> np.ndarray:
    """Scaled generating function Gbar evaluated at the rows of W.

    Gbar(w) = exp(sum_r alpha_r Z_r (w_r - 1)) prod_i F_i(alpha_r v_{i,r} w_r),
    with F_i given by the support recursion.
    """
    expo = (W - 1.0) @ ctx['arho0']
    logF = np.zeros(W.shape[0], dtype=complex)
    for i in range(ctx['M']):
        X = W * ctx['vs'][i, :]
        FS = np.zeros((W.shape[0], ctx['nmask']), dtype=complex)
        FS[:, 0] = 1.0                    # empty support: Phi(0) = 1
        for mask in range(1, ctx['nmask']):
            b = ctx['bits'][mask]
            sc = ctx['subcol'][mask]
            num = np.zeros(W.shape[0], dtype=complex)
            den = np.full(W.shape[0], ctx['muS'][i, mask], dtype=complex)
            for t in range(len(b)):
                xt = X[:, b[t]]
                num = num + xt * FS[:, sc[t]]
                den = den - xt
            FS[:, mask] = num / den
        logF = logF + np.log(np.sum(FS, axis=1))
    return np.exp(expo + logF)


def pfqn_clwoi(Z: Sequence[float],
               N: Sequence[int],
               mu: Optional[Union[Callable, List[Callable]]] = None,
               visits=None,
               options=None) -> Tuple[float, float]:
    """Normalizing constant of a closed delay + order-independent network.

    Inverts the multichain generating function with the lattice-Poisson
    algorithm of Choudhury, Leung and Whitt (J. ACM 42(5):935-970, 1995). This
    is the transform counterpart of the convolution routine :func:`pfqn_ncoi`;
    both return the same G(N) and differ in cost.

    An OI station factor is rational and available in closed form: splitting the
    count lattice by support S, on which mu_i(n) = mu_{i,S} is constant,

        (mu_{i,S} - sum_{r in S} v_{i,r} z_r) F_{i,S}(z)
            = sum_{r in S} v_{i,r} z_r F_{i,S-r}(z),   F_{i,{}} = 1,
        F_i(z) = sum_S F_{i,S}(z).

    The singularities are the |S| hyperplanes sum_{r in S} v_{i,r} z_r =
    mu_{i,S}, one per support, and the restrictive static scaling of CLW
    eqs. 5.41-5.46 runs on the expanded constraint matrix that lists one row per
    (station, nonempty support) pair.

    Cost: prod_r 2 l_r N_r contour points, each O(M R 2^R), against
    O(M prod_r (N_r+1)(N_r+2)/2) for the convolution of :func:`pfqn_ncoi`. The
    inversion is linear rather than quadratic in each population and returns G
    at the single population N.

    Parameters
    ----------
    Z : (R,) think-time demand vector of the aggregated delay node.
    N : (R,) closed population vector, finite.
    mu : list of callables, one per OI station. Each ``mu[m](n)`` returns the
        total service rate of station m at the per-class occupancy vector n and
        must depend on n only through its support. May be None for a pure delay
        network.
    visits : (M x R) array or list of (R,) vectors of class visit ratios
        weighting the balance recursion. Default: unit visits.
    options : dict with optional keys ``l`` (inner lattice parameters) and
        ``gamma`` (aliasing parameters). Defaults follow CLW.

    Returns
    -------
    (G, lG) : normalizing constant G(N), inf on overflow, and its natural log.

    Raises
    ------
    ValueError
        If a rate varies inside a support. Every handle is verified
        exhaustively on the count lattice before the inversion, because a
        violation would otherwise return a plausible but wrong G(N).
    """
    mu, M, N, R, Z, vis = _clw_prepare('pfqn_clwoi', Z, N, mu, visits)

    if np.any(N < 0):
        return 0.0, -np.inf
    if np.all(N == 0):
        return 1.0, 0.0

    lpar, gam = _clw_default_lattice(options, R)

    # Drop zero-population chains: the coefficient of z_r^0 is the generating
    # function restricted to z_r = 0, which kills every F_{i,S} with r in S.
    keep = np.where(N > 0)[0]
    Rk = keep.size
    Nk = N[keep]
    Zk = Z[keep]
    lpar = lpar[keep]
    gam = gam[keep]
    nmask = 2 ** Rk

    # Support rate table mu_{i,S}, S encoded as a bitmask over retained chains.
    muS = np.zeros((max(M, 1), nmask))
    for m in range(M):
        for mask in range(1, nmask):
            chi = np.zeros(R)
            chi[keep[[b for b in range(Rk) if (mask >> b) & 1]]] = 1.0
            muS[m, mask] = _clwoi_supportrate(mu[m], chi, m)

    _clwoi_checksupport(mu, muS, keep, Nk, R)

    # Per-mask chain lists and the S-minus-r column indices used by the recursion.
    bits = [None] * nmask
    subcol = [None] * nmask
    for mask in range(1, nmask):
        b = [j for j in range(Rk) if (mask >> j) & 1]
        bits[mask] = b
        subcol[mask] = [mask - (1 << j) for j in b]

    V = np.ones((max(M, 1), Rk))
    for m in range(M):
        V[m, :] = vis[m][keep]

    rad = 10.0 ** (-gam / (2.0 * lpar * Nk))

    # One constraint row per (station, nonempty support), holding the unit-pole
    # intensities v_{i,r}/mu_{i,S} of that singular hyperplane. Dominated
    # hyperplanes are dropped: support S of station i is implied by a superset
    # S' with mu_{i,S'} <= mu_{i,S}.
    rows = []
    for m in range(M):
        for mask in range(1, nmask):
            dominated = False
            for mask2 in range(1, nmask):
                if mask2 != mask and (mask & mask2) == mask \
                        and muS[m, mask2] <= muS[m, mask] * (1 + 1e-12):
                    dominated = True
                    break
            if dominated:
                continue
            row = np.zeros(Rk)
            for b in bits[mask]:
                row[b] = V[m, b] / muS[m, mask]
            rows.append(row)
    Lt = np.array(rows) if rows else np.zeros((1, Rk))

    alpha = _clw_scaling(Lt, Nk, lpar, rad, Zk)

    ctx = {'N': Nk, 'l': lpar, 'r': rad, 'p': Rk, 'M': M,
           'arho0': alpha * Zk, 'vs': V * alpha, 'muS': muS,
           'bits': bits, 'subcol': subcol, 'nmask': nmask,
           'chunk': max(1, int(2e6 // nmask))}

    gbarN = _clw_invert(0, np.zeros(0, dtype=complex), ctx, _clwoi_gbar)
    return _clw_recover(gbarN, ctx['arho0'], Nk, alpha)


def _clwjd_regionrate(murate: Callable, t: np.ndarray, keep: np.ndarray, R: int,
                      N: np.ndarray, lrow: np.ndarray, m: int) -> float:
    """Rate of one clipped region, with the constancy check.

    The region {n : t(n) = t} pins every unsaturated coordinate and leaves the
    saturated ones free above the cutoff, so the rate is probed at the region
    representative and at two larger occupancies of the same region.
    """
    nrep = np.zeros(R)
    nrep[keep] = t
    rate = float(murate(nrep))
    if not (rate > 0):
        raise ValueError('pfqn_clwjd: station %d has a non-positive rate on a reachable region.'
                         % (m + 1))
    sat = np.where(t == lrow)[0]
    if sat.size == 0:
        return rate
    for pass_ in range(2):
        nprobe = nrep.copy()
        if pass_ == 0:
            nprobe[keep[sat]] = N[keep[sat]]
        else:
            nprobe[keep[sat]] = np.maximum(t[sat], np.floor((t[sat] + N[keep[sat]]) / 2.0))
        if np.any(nprobe != nrep):
            rt = float(murate(nprobe))
            if abs(rt - rate) > 1e-9 * max(1.0, abs(rate)):
                raise ValueError(
                    'pfqn_clwjd: station %d has a rate that varies within a clipped region '
                    '(mu=%g at the region representative, %g above the cutoff). pfqn_clwjd '
                    'requires the rate to saturate at lcut; raise lcut or use pfqn_ncjd.'
                    % (m + 1, rate, rt))
    return rate


def _clwjd_gbar(W: np.ndarray, ctx: dict) -> np.ndarray:
    """Scaled generating function Gbar with the clipped-region recursion."""
    expo = (W - 1.0) @ ctx['arho0']
    logF = np.zeros(W.shape[0], dtype=complex)
    for i in range(ctx['M']):
        X = W * ctx['vs'][i, :]
        nt = ctx['ntreg'][i]
        FT = np.zeros((W.shape[0], nt), dtype=complex)
        FT[:, 0] = 1.0                    # empty region: Phi(0) = 1
        for tl in range(1, nt):
            num = np.zeros(W.shape[0], dtype=complex)
            for (chain, col) in ctx['regDec'][i][tl]:
                num = num + X[:, chain] * FT[:, col]
            den = np.full(W.shape[0], ctx['muT'][i][tl], dtype=complex)
            for chain in ctx['regSat'][i][tl]:
                den = den - X[:, chain]
            FT[:, tl] = num / den
        logF = logF + np.log(np.sum(FT, axis=1))
    return np.exp(expo + logF)


def pfqn_clwjd(Z: Sequence[float],
               N: Sequence[int],
               mu: Optional[Union[Callable, List[Callable]]] = None,
               visits=None,
               lcut=None,
               options=None) -> Tuple[float, float]:
    """Normalizing constant of a closed delay + limited joint-dependent network.

    Joint-dependent generalization of :func:`pfqn_clwoi`, which is the case
    lcut = 1. Station i has a rate that reads the whole per-class occupancy but
    saturates coordinatewise: with a cutoff vector l_i,

        mu_i(n) = c_{i,t},  t = (min(n_1, l_{i,1}), ..., min(n_R, l_{i,R})),

    so the clipped vector t ranges over a finite box and the station factor is
    rational, with

        (mu_{i,t} - sum_{r: t_r = l_{i,r}} v_{i,r} z_r) F_{i,t}(z)
            = sum_{r: t_r >= 1} v_{i,r} z_r F_{i,t-e_r}(z),   F_{i,0} = 1.

    The singular hyperplanes are indexed by the SATURATED sets, at most 2^R per
    station however large the cutoffs are. G(N) is then recovered by R nested
    lattice-Poisson inversions exactly as in :func:`pfqn_clwoi`.

    Cost: prod_r 2 l_r N_r contour points, each O(M R prod_r (lcut_{i,r}+1)).
    With lcut = N the region box is the whole lattice and the convolution of
    :func:`pfqn_ncjd` wins outright; the inversion pays off when the joint
    dependence saturates early.

    Parameters
    ----------
    Z : (R,) think-time demand vector of the aggregated delay node.
    N : (R,) closed population vector, finite.
    mu : list of callables, one per LJD station, each taking the per-class
        occupancy vector. May be None for a pure delay network.
    visits : (M x R) array or list of (R,) vectors of class visit ratios.
    lcut : (M x R) matrix of per-station per-class saturation cutoffs >= 1, or a
        scalar/row broadcast to every station. Entries are clipped to N, which
        is exact. Default: N (no truncation).
    options : dict with optional keys ``l`` and ``gamma``.

    Returns
    -------
    (G, lG) : normalizing constant G(N), inf on overflow, and its natural log.
    """
    mu, M, N, R, Z, vis = _clw_prepare('pfqn_clwjd', Z, N, mu, visits)

    if np.any(N < 0):
        return 0.0, -np.inf
    if np.all(N == 0):
        return 1.0, 0.0

    # saturation cutoffs, broadcast and clipped to the reachable lattice
    if lcut is None or np.asarray(lcut).size == 0:
        L = np.tile(N.astype(float), (max(M, 1), 1))
    else:
        lc = np.asarray(lcut, dtype=float)
        if lc.ndim == 0:
            L = float(lc) * np.ones((max(M, 1), R))
        elif lc.ndim == 1:
            L = np.tile(lc.ravel(), (max(M, 1), 1))
        else:
            L = lc
    if M > 0 and (L.shape[0] != M or L.shape[1] != R):
        raise ValueError('pfqn_clwjd: lcut must be (M x R), a (R,) row, or a scalar.')
    L = np.minimum(np.maximum(np.round(L), 1), np.tile(np.maximum(N, 1), (L.shape[0], 1)))
    L = L.astype(int)

    lpar, gam = _clw_default_lattice(options, R)

    keep = np.where(N > 0)[0]
    Rk = keep.size
    Nk = N[keep]
    Zk = Z[keep]
    lpar = lpar[keep]
    gam = gam[keep]
    Lk = L[:, keep]

    # Per-station region tables over the clipped box, in mixed radix so that
    # t - e_r always precedes t.
    ntreg = np.ones(max(M, 1), dtype=int)
    muT = [None] * max(M, 1)
    regDec = [None] * max(M, 1)
    regSat = [None] * max(M, 1)
    satMask = [None] * max(M, 1)
    for m in range(M):
        rad_m = Lk[m, :] + 1
        st = np.concatenate(([1], np.cumprod(rad_m[:-1]))).astype(int)
        nt = int(np.prod(rad_m))
        ntreg[m] = nt
        muT[m] = np.zeros(nt)
        regDec[m] = [None] * nt
        regSat[m] = [None] * nt
        satMask[m] = np.zeros(nt, dtype=int)
        for tl in range(nt):
            t = (tl // st) % rad_m
            dec = []
            sat = []
            smask = 0
            for b in range(Rk):
                if t[b] >= 1:
                    dec.append((b, tl - int(st[b])))
                if t[b] == Lk[m, b]:
                    sat.append(b)
                    smask += 1 << b
            regDec[m][tl] = dec
            regSat[m][tl] = sat
            satMask[m][tl] = smask
            if tl == 0:
                muT[m][0] = 1.0           # F_0 = Phi(0) = 1, rate unused
            else:
                muT[m][tl] = _clwjd_regionrate(mu[m], t, keep, R, N, Lk[m, :], m)

    V = np.ones((max(M, 1), Rk))
    for m in range(M):
        V[m, :] = vis[m][keep]

    rad = 10.0 ** (-gam / (2.0 * lpar * Nk))

    # Binding rate of each saturated set: regions sharing a saturated set share
    # the hyperplane, so the smallest rate constrains.
    nmask = 2 ** Rk
    muS = np.full((max(M, 1), nmask), np.inf)
    for m in range(M):
        for tl in range(ntreg[m]):
            sm = satMask[m][tl]
            if sm > 0:
                muS[m, sm] = min(muS[m, sm], muT[m][tl])

    rows = []
    for m in range(M):
        for mask in range(1, nmask):
            if not np.isfinite(muS[m, mask]):
                continue
            dominated = False
            for mask2 in range(1, nmask):
                if mask2 != mask and (mask & mask2) == mask \
                        and np.isfinite(muS[m, mask2]) \
                        and muS[m, mask2] <= muS[m, mask] * (1 + 1e-12):
                    dominated = True
                    break
            if dominated:
                continue
            row = np.zeros(Rk)
            for b in range(Rk):
                if (mask >> b) & 1:
                    row[b] = V[m, b] / muS[m, mask]
            rows.append(row)
    Lt = np.array(rows) if rows else np.zeros((1, Rk))

    alpha = _clw_scaling(Lt, Nk, lpar, rad, Zk)

    ctx = {'N': Nk, 'l': lpar, 'r': rad, 'p': Rk, 'M': M,
           'arho0': alpha * Zk, 'vs': V * alpha, 'muT': muT,
           'regDec': regDec, 'regSat': regSat, 'ntreg': ntreg,
           'chunk': max(1, int(2e6 // max(int(np.max(ntreg)), 1)))}

    gbarN = _clw_invert(0, np.zeros(0, dtype=complex), ctx, _clwjd_gbar)
    return _clw_recover(gbarN, ctx['arho0'], Nk, alpha)
