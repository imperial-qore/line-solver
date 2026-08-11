"""
Order-independent (OI) normalizing-constant tools.

Native Python port of the MATLAB OI functional-server routines:
    pfqn_ncoi  - balanced-fairness normalizing constant for a closed network
                 of order-independent stations and one aggregated delay node,
                 tabulated on the count lattice (macrostate). The microstate
                 counterpart, needed only when a swap graph imposes a placement
                 order, is pfqn_pas_nc in pas.py.
    pfqn_oi_fnc - OI generalization of the load-dependent functional server
                 (FNC) of Casale, "On Single-Class Load-Dependent Normalizing
                 Constant Equations", QEST 2006 (Theorem 3, Corollary 1).

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_ncoi.m, pfqn_oi_fnc.m
    Bonald and Proutiere, "Insensitive bandwidth sharing in data networks"
    (2003); Casale, QEST 2006.
"""

import numpy as np
from typing import Callable, List, Optional, Sequence, Tuple, Union


def pfqn_ncoi(Z: Sequence[float],
              N: Sequence[int],
              mu: Optional[Union[Callable, List[Callable]]] = None,
              visits=None,
              options=None) -> Tuple[float, float, np.ndarray]:
    """Normalizing constant of a closed OI + single-delay product-form network.

    The OI stations are analyzed by the balanced-fairness recursion of Bonald
    and Proutiere (2003) combined with the multichain convolution over
    stations. For a single OI station with rank rate mu(supp(n)) the balance
    function is Phi(0) = 1, Phi(n) = (1/mu(n)) sum_{r: n_r>0} Phi(n - e_r), and
    G(N) is the convolution of the per-station balance functions with the
    multinomial delay factor F_Z(n) = prod_r Z_r^{n_r}/n_r!,
    g_0 = F_Z, g_m(n) = sum_{0<=x<=n} Phi_m(x) g_{m-1}(n-x), G(N) = g_M(N).

    This is a MACROSTATE routine: everything is tabulated over the count
    lattice 0 <= n <= N, never over orderings, which is legitimate because an
    OI rate is permutation-invariant so Phi closes on the count vector. With a
    non-empty swap graph that closure fails; use :func:`pfqn_pas_nc` instead.

    Cost: O(M R L) for the balance functions and O(M prod_r (N_r+1)(N_r+2)/2)
    for the convolutions, with L = prod_r (N_r+1).

    Parameters
    ----------
    Z : (R,) think-time demand vector of the aggregated delay node.
    N : (R,) closed population vector, finite.
    mu : list of callables, one per OI station. Each ``mu[m](n)`` returns the
        total service rate of station m given the per-class occupancy (count)
        vector n. A station state with non-positive rate is unreachable and is
        assigned a zero balance value. May be empty/None for a pure delay
        network.
    options : accepted for signature parity; unused.

    Returns
    -------
    (G, lG, Gtab) : normalizing constant, its natural log, and the WHOLE
        lattice of normalizing constants, ``Gtab[dot(n, strides)] = G(n)`` for
        every ``0 <= n <= N`` with ``strides = cumprod([1, N[:-1]+1])``. The
        convolution produces this table anyway, so a caller needing G at more
        than one population must index this output, NOT re-call the routine per
        population -- the latter costs a needless factor ``prod_r (N_r+1)``.
    """
    N = np.round(np.asarray(N, dtype=float)).astype(int).ravel()
    R = N.size

    if mu is None:
        mu = []
    elif callable(mu):
        mu = [mu]
    else:
        mu = list(mu)

    if Z is None or len(np.asarray(Z).ravel()) == 0:
        Z = np.zeros(R)
    Z = np.asarray(Z, dtype=float).ravel()
    if Z.size != R:
        raise ValueError('pfqn_ncoi: Z and N must have the same number of classes.')
    if np.any(N < 0):
        raise ValueError('pfqn_ncoi requires finite (closed) populations.')

    # Count lattice 0 <= n <= N, flattened column-major; lin(v) = sum(v*strides)
    # is linear, so lin(x+y) = lin(x) + lin(y), which the convolution exploits.
    dims = N + 1
    ngrid = int(np.prod(dims))
    strides = np.concatenate(([1], np.cumprod(dims[:-1]))).astype(int)
    counts = np.stack(np.unravel_index(np.arange(ngrid), tuple(dims), order='F'), axis=1).astype(int)
    order = np.argsort(counts.sum(axis=1), kind='stable')

    # Delay balance function: the multinomial factor F_Z(n). A class with
    # population but no delay demand makes the state infeasible (weight zero).
    from scipy.special import gammaln
    g = np.zeros(ngrid)
    for k in range(ngrid):
        v = counts[k]
        logf = 0.0
        feas = True
        for r in range(R):
            if v[r] > 0:
                if Z[r] <= 0:
                    feas = False
                    break
                logf += v[r] * np.log(Z[r]) - gammaln(v[r] + 1)
        if feas:
            g[k] = float(np.exp(logf))

    # Per-OI-station class visit ratios: visits[m] is the 1xR vector entering
    # the v-weighted balanced-fairness balance Phi^v(n)=(1/mu(n)) sum_r v_r
    # Phi^v(n-e_r); None -> unit visits.
    # Convolve in one OI station at a time.
    for m in range(len(mu)):
        vism = None if visits is None else visits[m]
        Phim = _oi_nc_balance(counts, order, strides, mu[m], R, ngrid, vism)
        gnext = np.zeros(ngrid)
        for kx in range(ngrid):
            px = Phim[kx]
            if px == 0.0:
                continue
            x = counts[kx]
            rem = N - x
            # Linear indices of the sub-box 0 <= y <= rem, built by mixed radix.
            ylin = np.zeros(1, dtype=int)
            for d in range(R):
                ylin = (ylin[:, None] + strides[d] * np.arange(rem[d] + 1)[None, :]).ravel()
            base = int(np.dot(x, strides))
            np.add.at(gnext, base + ylin, px * g[ylin])
        g = gnext

    # g now holds G(n) for EVERY n on the lattice, not just n = N.
    G = float(g[ngrid - 1])
    lG = float(np.log(G)) if G > 0 else -np.inf
    return G, lG, g


def _oi_nc_balance(counts: np.ndarray, order: np.ndarray, strides: np.ndarray,
                   murate: Callable, R: int, ngrid: int, vis=None) -> np.ndarray:
    """v-weighted balanced-fairness recursion over the count lattice.

    Phi(n) = (1/mu(n)) sum_{r: n_r>0} v_r Phi(n-e_r); vis is the 1xR class visit
    vector (None -> unit visits).
    """
    Phi = np.zeros(ngrid)
    for k in order:
        n = counts[k]
        if not n.any():
            Phi[k] = 1.0
            continue
        mun = murate(n)
        if not (mun > 0):
            continue  # unreachable station state: zero balance value
        acc = 0.0
        for r in range(R):
            if n[r] > 0:
                vr = 1.0 if vis is None else vis[r]
                acc += vr * Phi[k - strides[r]]
        Phi[k] = acc / mun
    return Phi


def pfqn_mvaoi(Z: Sequence[float],
               N: Sequence[int],
               mu: Union[Callable, List[Callable]],
               Dli: Optional[Sequence[Sequence[float]]] = None,
               visits=None,
               options=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Mean-value analysis of a closed product-form OI network.

    Mean-value counterpart of :func:`pfqn_ncoi` and the marginal form
    :func:`pfqn_mvaoi_marg`: for a closed product-form network of an aggregated
    infinite-server (delay) node, any number of load-independent (LI) single-server
    product-form queues, and any number of order-independent (OI) stations, it
    returns the same exact per-class throughput and queue-lengths WITHOUT computing
    any normalizing constant or joint marginal, using only mean quantities. It is
    the composition-dependent generalization of the Conditional MVA (CMVA) of
    Casale, "A Note on Stable Flow-Equivalent Aggregation in Closed Networks"
    (QUESTA 2009), extended to MULTIPLE OI stations by carrying one rate-shift
    vector ``s_i`` per OI station i (row ``i`` of the shift matrix ``S``).

    Throughout, ``r`` and ``s`` index job classes; ``i`` indexes OI stations; ``j``
    indexes LI queues. State ``(S, Nn)`` is processed by increasing ``sum(Nn)``;
    each OI station keeps its own ``D^i``, ``rho^i`` and ``Q^i`` recursions driven
    by the common throughput ``X^{(S)}(Nn)``, and the population conservation
    aggregates every station's contribution::

        Nn_r = X_r Z_r + sum_j Q^{(j)}_r + sum_i Q^{(i)}_r,

    with the LI queue term ``Q^{(j)}_r = X_r D_{j,r}(1 + sum_s Q^{(j)}_s(Nn - e_r))``.

    Parameters
    ----------
    Z : (R,) think-time demand vector of the aggregated delay node.
    N : (R,) closed population vector, finite.
    mu : callable or list of callables ``mu_i(n)`` returning the OI total service
        rate of station i for the per-class occupancy (count) vector n. A bare
        callable is accepted as the single-station shorthand.
    Dli : (J, R) per-class demand matrix of the LI single-server queues; None or
        empty when J = 0.
    options : accepted for signature parity; unused.

    Returns
    -------
    (X, Qoi, Qli, Qdelay, Soi) : per-class throughput (R,), OI queue-lengths
        (K, R), LI queue-lengths (J, R), delay queue-length (R,) = ``X * Z``, and
        the per-class mean number of IN-SERVICE jobs at each OI station (K, R).
        ``Soi[i, r] = E[sir_r]`` counts the class-r jobs receiving a strictly
        positive rank rate (see :func:`pfqn_oi_insvc`); the utilization of OI
        station i is ``Soi[i, r] / c_i``. Unlike X/Qoi/Qli, which are pure
        mean-value quantities, Soi is a distributional statistic and is obtained
        from the OI count marginal assembled from the zero-shift throughputs
        ``X^{(0)}(k)`` already cached by the mean-value recursion above (no
        normalizing constant is formed).
    """
    if callable(mu):
        mu = [mu]
    mu = list(mu)
    if len(mu) == 0 or not all(callable(m) for m in mu):
        raise ValueError('mu must be a (nonempty) list of OI rate callables.')
    K = len(mu)

    Z = np.asarray(Z, dtype=float).ravel()
    N = np.round(np.asarray(N, dtype=float)).astype(int).ravel()
    R = N.size
    if np.any(~np.isfinite(N)):
        raise ValueError('pfqn_mvaoi requires finite (closed) populations.')
    if Dli is None:
        Dli = np.zeros((0, R))
    Dli = np.asarray(Dli, dtype=float).reshape(-1, R) if np.size(Dli) else np.zeros((0, R))
    J = Dli.shape[0]
    ei = np.eye(R, dtype=int)

    Xc: dict = {}          # X^{(S)}(Nn)
    Qlc: dict = {}         # Qli^{(S)}(Nn) (J, R)
    Dc = [dict() for _ in range(K)]   # D_i^{(S)}(Nn)
    Qc = [dict() for _ in range(K)]   # Q_i^{(S)}(Nn)
    rho_cache = [dict() for _ in range(K)]

    def skey(S, Nn):
        return (tuple(np.asarray(S, dtype=int).ravel()), tuple(np.asarray(Nn, dtype=int)))

    def shift_plus(S, i, r):
        Sp = np.array(S, dtype=int, copy=True)
        Sp[i] = Sp[i] + ei[r]
        return Sp

    def rho(i, r, S, M):
        M = np.asarray(M, dtype=int)
        key = (skey(S, M), r)
        v = rho_cache[i].get(key)
        if v is not None:
            return v
        if M.sum() == 0:
            rho_cache[i][key] = 1.0
            return 1.0
        s = next(ss for ss in range(R) if ss != r and M[ss] > 0)
        xu = Xc[skey(S, M)][s]
        xu2 = Xc[skey(shift_plus(S, i, r), M)][s]
        val = rho(i, r, S, M - ei[s]) * (xu / xu2 if xu2 > 0 else 0.0)
        rho_cache[i][key] = val
        return val

    zeroS = np.zeros((K, R), dtype=int)
    k0 = skey(zeroS, np.zeros(R, dtype=int))
    Xc[k0] = np.zeros(R)
    Qlc[k0] = np.zeros((J, R))
    for i in range(K):
        Dc[i][k0] = np.zeros(R)
        Qc[i][k0] = np.zeros(R)

    # see _kb/03-api-layer.md for rationale
    from itertools import product
    perclass = []
    for r in range(R):
        rows = [comp[:K + 1] for comp in _compositions(int(N[r]), K + 2)]
        perclass.append(rows)
    states = []
    for combo in product(*perclass):
        S = np.zeros((K, R), dtype=int)
        Nn = np.zeros(R, dtype=int)
        for r in range(R):
            row = combo[r]
            for i in range(K):
                S[i, r] = row[i]
            Nn[r] = row[K]
        states.append((S, Nn))
    states.sort(key=lambda st: int(st[1].sum()))

    for S, Nn in states:
        kk = skey(S, Nn)
        if Nn.sum() == 0:
            Xc[kk] = np.zeros(R)
            Qlc[kk] = np.zeros((J, R))
            for i in range(K):
                Dc[i][kk] = np.zeros(R)
                Qc[i][kk] = np.zeros(R)
            continue

        Dt = np.zeros((K, R))
        Qsub = [np.zeros((R, R)) for _ in range(K)]
        for i in range(K):
            for r in range(R):
                if Nn[r] == 0:
                    continue
                if Nn[r] == 1:
                    vfac = 1.0 if visits is None else visits[i][r]
                    mur = mu[i](S[i] + ei[r])
                    if mur > 0:
                        Dt[i, r] = (vfac / mur) * rho(i, r, S, Nn - ei[r])
                else:
                    Nmr = Nn - ei[r]
                    xs = Xc[skey(S, Nmr)][r]
                    xs2 = Xc[skey(shift_plus(S, i, r), Nmr)][r]
                    if xs2 > 0:
                        Dt[i, r] = (xs / xs2) * Dc[i][skey(S, Nmr)][r]
            for s in range(R):
                if Nn[s] > 0:
                    Qsub[i][s] = Qc[i][skey(shift_plus(S, i, s), Nn - ei[s])]

        betaLI = np.zeros((J, R))
        for r in range(R):
            if Nn[r] == 0:
                continue
            QsubLI = Qlc[skey(S, Nn - ei[r])]
            for j in range(J):
                betaLI[j, r] = Dli[j, r] * (1.0 + float(np.sum(QsubLI[j, :])))

        idx = [r for r in range(R) if Nn[r] > 0]
        A = np.zeros((len(idx), len(idx)))
        for a, r in enumerate(idx):
            for b, s in enumerate(idx):
                if s == r:
                    val = Z[r] + float(np.sum(betaLI[:, r]))
                    for i in range(K):
                        val += Dt[i, r] * (1.0 + Qsub[i][r, r])
                    A[a, b] = val
                else:
                    val = 0.0
                    for i in range(K):
                        val += Dt[i, s] * Qsub[i][s, r]
                    A[a, b] = val
        # see _kb/03-api-layer.md for rationale
        Xk = np.zeros(R)
        for a, r in enumerate(idx):
            denom = A[a, a]
            for b, s in enumerate(idx):
                if b != a:
                    Xner = Xc[skey(S, Nn - ei[r])]   # X(S, Nn-e_r)
                    Xnes = Xc[skey(S, Nn - ei[s])]   # X(S, Nn-e_s)
                    if Xnes[r] > 0:
                        denom += A[a, b] * (Xner[s] / Xnes[r])
            if denom > 0:
                Xk[r] = Nn[r] / denom

        QkLI = np.zeros((J, R))
        for r in range(R):
            if Nn[r] == 0:
                continue
            QkLI[:, r] = Xk[r] * betaLI[:, r]
        for i in range(K):
            U = Dt[i] * Xk
            Qi = np.array([U[r] + float(np.dot(U, Qsub[i][:, r])) for r in range(R)])
            Qc[i][kk] = Qi
            Dc[i][kk] = Dt[i].copy()
        Xc[kk] = Xk
        Qlc[kk] = QkLI

    keyN = skey(zeroS, N)
    X = Xc[keyN]
    Qoi = np.array([Qc[i][keyN] for i in range(K)]) if K > 0 else np.zeros((0, R))
    Qli = Qlc[keyN]
    # see _kb/03-api-layer.md for rationale
    Soi = _LazyArray(lambda: _oi_insvc_means(N, mu, Xc, skey, zeroS, K, R))
    return X, Qoi, Qli, X * Z, Soi


class _LazyArray(object):
    """Deferred ndarray: ``factory()`` runs on first use, then is cached.

    Python has no ``nargout``, so a routine whose last output is expensive and
    frequently unused cannot skip it by inspecting the call site. This proxy
    supplies the equivalent: it is returned in place of the array, forwards
    every attribute, item access, iteration and operator to the materialized
    value, and converts through ``__array__`` so ``np.asarray`` and numpy
    ufuncs see a plain ndarray.
    """

    __slots__ = ('_factory', '_value')

    def __init__(self, factory):
        self._factory = factory
        self._value = None

    def value(self) -> np.ndarray:
        """Materialize (once) and return the underlying ndarray."""
        if self._value is None:
            self._value = self._factory()
        return self._value

    def __array__(self, dtype=None):
        arr = self.value()
        return arr if dtype is None else arr.astype(dtype)

    def __getattr__(self, name):
        # The two slots are resolved by the descriptor protocol; reaching here
        # for them means they are unset, and delegating would recurse.
        if name in ('_factory', '_value'):
            raise AttributeError(name)
        return getattr(self.value(), name)

    def __getitem__(self, key):
        return self.value()[key]

    def __setitem__(self, key, val):
        self.value()[key] = val

    def __len__(self):
        return len(self.value())

    def __iter__(self):
        return iter(self.value())

    def __bool__(self):
        return bool(self.value())

    def __repr__(self):
        return repr(self.value())


def _install_lazyarray_operators():
    """Forward the arithmetic/comparison protocol to the materialized array.

    ``__getattr__`` is not consulted for operators invoked through the type, so
    each dunder has to exist on the class itself.
    """
    binary = ['add', 'sub', 'mul', 'matmul', 'truediv', 'floordiv', 'mod',
              'pow', 'and', 'or', 'xor', 'lshift', 'rshift']
    for op in binary:
        name = '__%s__' % op
        rname = '__r%s__' % op

        def _fwd(self, other, _name=name):
            return getattr(self.value(), _name)(other)

        def _rfwd(self, other, _rname=rname):
            return getattr(self.value(), _rname)(other)

        setattr(_LazyArray, name, _fwd)
        setattr(_LazyArray, rname, _rfwd)
    for op in ['neg', 'pos', 'abs', 'invert']:
        name = '__%s__' % op

        def _ufwd(self, _name=name):
            return getattr(self.value(), _name)()

        setattr(_LazyArray, name, _ufwd)
    for op in ['eq', 'ne', 'lt', 'le', 'gt', 'ge']:
        name = '__%s__' % op

        def _cfwd(self, other, _name=name):
            return getattr(self.value(), _name)(other)

        setattr(_LazyArray, name, _cfwd)


_install_lazyarray_operators()


def _oi_insvc_means(N, mu, Xc, skey, zeroS, K, R):
    """Mean number of in-service jobs per class at each OI station, E[sir_r].

    Built from the OI count marginal pM_i(n|k) on the cached zero-shift
    throughputs X^{(0)}(k). Exact because in product form
    ``pM_i(n|k) = Phi_i(n) G_{-i}(k-n)/G(k)`` and ``X_r(k) = G(k-e_r)/G(k)``, so
    the recursion below reproduces the balanced-fairness identity for Phi_i::

        pM_i(n|k) = (1/mu_i(n)) sum_r X_r(k) pM_i(n-e_r|k-e_r)
        pM_i(0|k) = 1 - sum_{n != 0} pM_i(n|k)
    """
    N = np.asarray(N, dtype=int).ravel()
    shp = N + 1
    stride = np.ones(R, dtype=int)
    for d in range(1, R):
        stride[d] = stride[d - 1] * shp[d - 1]
    total = int(np.prod(shp))
    subs = np.zeros((total, R), dtype=int)
    for i in range(total):
        li = i
        for d in range(R):
            subs[i, d] = li % shp[d]
            li //= shp[d]
    order = sorted(range(total), key=lambda i: subs[i].sum())

    Xk = np.zeros((total, R))
    for i in range(total):
        Xk[i, :] = Xc[skey(zeroS, subs[i])]

    Soi = np.zeros((K, R))
    for m in range(K):
        gm, _, _ = pfqn_oi_insvc(mu[m], N)
        muv = np.zeros(total)
        for i in range(total):
            if subs[i].sum() > 0:
                muv[i] = float(mu[m](subs[i].copy()))
        pMv = np.zeros((total, total))
        pMv[0, 0] = 1.0                       # pM(0|0) = 1
        for b in order:
            k = subs[b]
            if k.sum() == 0:
                continue
            acc0 = 0.0
            for a in order:
                n = subs[a]
                if n.sum() == 0 or np.any(n > k) or muv[a] <= 0:
                    continue
                acc = 0.0
                for r in range(R):
                    if n[r] > 0:
                        acc += Xk[b, r] * pMv[a - stride[r], b - stride[r]]
                pMv[a, b] = acc / muv[a]
                acc0 += pMv[a, b]
            pMv[0, b] = 1.0 - acc0            # empty state by complement
        idxN = int(np.sum(N * stride))
        for r in range(R):
            Soi[m, r] = float(pMv[:, idxN] @ gm[:, r])
    return Soi


def _compositions(m: int, p: int) -> List[Tuple[int, ...]]:
    """All nonnegative integer tuples of length p summing to m."""
    if p == 1:
        return [(m,)]
    out = []
    for first in range(m + 1):
        for rest in _compositions(m - first, p - 1):
            out.append((first,) + rest)
    return out


def pfqn_mvaoi_marg(D: Sequence[Sequence[float]],
                    N: Sequence[int],
                    isDelay: Sequence[bool],
                    mu: List[Optional[Callable]]
                    ) -> Tuple[np.ndarray, np.ndarray]:
    """Exact marginal load-dependent MVA for OI networks.

    Marginal-distribution counterpart of :func:`pfqn_mvaoi`. Carries, for each OI
    station, its joint count-vector marginal ``pM_i(n | k)`` and closes the per-class
    throughput by population conservation. Handles delay + LI product-form queues +
    any number of OI stations.

    Parameters
    ----------
    D : (M, R) per-class demand at every station (OI rows ignored).
    N : (R,) closed population vector, finite.
    isDelay : (M,) True for infinite-server (delay) stations.
    mu : length-M list; ``mu[i]`` is the OI rate callable of the count vector n, or
        None for non-OI stations.

    Returns
    -------
    (XN, QN) : per-class throughput (R,) and per-station queue-lengths (M, R).
    """
    D = np.asarray(D, dtype=float)
    M = D.shape[0]
    N = np.round(np.asarray(N, dtype=float)).astype(int).ravel()
    R = N.size
    isDelay = np.asarray(isDelay, dtype=bool).ravel()
    oi_list = [i for i in range(M) if mu[i] is not None]
    nOI = len(oi_list)
    if nOI == 0:
        raise ValueError('pfqn_mvaoi_marg requires at least one OI station.')
    ei = np.eye(R, dtype=int)
    from itertools import product

    def vecs_upto(k):
        return [np.array(t, dtype=int)
                for t in product(*[range(int(k[r]) + 1) for r in range(R)])]

    zero = tuple(np.zeros(R, dtype=int))
    Xc = {zero: np.zeros(R)}
    Qc = {zero: np.zeros((M, R))}
    pM = [{zero: {zero: 1.0}} for _ in range(nOI)]

    def _muM(fun, n):
        cls = np.repeat(np.arange(R), np.asarray(n, dtype=int))
        if cls.size == 0:
            return 0.0
        return float(fun(cls))

    def pM_given(o, k, Xk):
        fun = mu[oi_list[o]]
        pmk = {}
        for n in vecs_upto(k):
            if n.sum() == 0:
                continue
            rate = _muM(fun, n)
            if rate <= 0:
                continue
            acc = 0.0
            for r in range(R):
                if n[r] >= 1 and k[r] >= 1:
                    acc += Xk[r] * pM[o][tuple(k - ei[r])].get(tuple(n - ei[r]), 0.0)
            pmk[tuple(n)] = acc / rate
        pmk[zero] = 1.0 - sum(pmk.values())
        return pmk

    pops = sorted(vecs_upto(N), key=lambda v: (v.sum(), tuple(v)))
    for k in pops:
        kt = tuple(k)
        if k.sum() == 0:
            continue
        Rfix = np.zeros((M, R))
        A = np.zeros(R)
        for r in range(R):
            if k[r] == 0:
                continue
            Qkr = Qc[tuple(k - ei[r])]
            for i in range(M):
                if mu[i] is not None:
                    continue
                if isDelay[i]:
                    Rfix[i, r] = D[i, r]
                else:
                    Rfix[i, r] = D[i, r] * (1.0 + float(np.sum(Qkr[i, :])))
                A[r] += Rfix[i, r]

        Xk = np.array([k[r] / (A[r] + 1.0) if k[r] > 0 else 0.0 for r in range(R)])
        pmk = [None] * nOI
        for _ in range(2000):
            QMtot = np.zeros(R)
            for o in range(nOI):
                pmk[o] = pM_given(o, k, Xk)
                for n_t, p in pmk[o].items():
                    QMtot += np.array(n_t) * p
            Xnew = np.array([max((k[r] - QMtot[r]) / A[r], 0.0) if (k[r] > 0 and A[r] > 0)
                             else 0.0 for r in range(R)])
            if np.max(np.abs(Xnew - Xk)) < 1e-13:
                Xk = Xnew
                break
            Xk = 0.5 * Xk + 0.5 * Xnew

        Qk = np.zeros((M, R))
        for o in range(nOI):
            pmk[o] = pM_given(o, k, Xk)
            QM = np.zeros(R)
            for n_t, p in pmk[o].items():
                QM += np.array(n_t) * p
            Qk[oi_list[o], :] = QM
        for r in range(R):
            for i in range(M):
                if mu[i] is not None:
                    continue
                Qk[i, r] = Xk[r] * Rfix[i, r]
        Xc[kt] = Xk
        Qc[kt] = Qk
        for o in range(nOI):
            pM[o][kt] = pmk[o]

    Nt = tuple(N)
    return Xc[Nt].copy(), Qc[Nt].copy()


def pfqn_oi_insvc(oirate: Callable,
                  N: Sequence[int],
                  options=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Conditional mean number of in-service jobs per class at an OI station.

    This is the quantity underlying the LINE utilization convention at
    order-independent stations, ``U_r = E[sir_r] / c`` with ``c`` the number of
    servers and ``sir_r`` the number of class-r jobs receiving a strictly
    positive service rate.

    In an OI station the state is the ordered list ``c = (c_1,...,c_n)`` of job
    classes (position 1 = head) and the job in position p is served at the rank
    rate increment ``Delta_p(c) = mu(c_1..c_p) - mu(c_1..c_{p-1})``, so the total
    rate telescopes to ``mu(c)``. Position p is in service when
    ``Delta_p(c) > 0``, and ``sir_r(c) = #{p : c_p = r, Delta_p(c) > 0}``. Note
    that ``sir_r`` counts JOBS, not servers: a single job served concurrently by
    several compatible servers counts once. This matches the definition used by
    the exact CTMC solver (``State.to_marginal``, PAS branch) and by LDES.

    Because ``mu`` is permutation-invariant, the unnormalized weight of an
    ordering c of the multiset n factorizes over its prefixes as
    ``w(c) = prod_p 1/mu(n(c_1..c_p))``, and ``Phi(n) = sum_c w(c)`` obeys the
    balanced-fairness recursion (condition on the tail element)::

        Phi(0) = 1,   Phi(n) = (1/mu(n)) sum_{r: n_r>0} Phi(n - e_r)

    Conditioning the same way and using
    ``sir_r(c) = sir_r(c_1..c_{|n|-1}) + [c_{|n|} = r] * 1{mu(n) > mu(n - e_r)}``
    gives the companion recursion for ``Xi_r(n) = sum_c w(c) sir_r(c)``::

        Xi_r(0) = 0
        Xi_r(n) = (1/mu(n)) [ sum_{s: n_s>0} Xi_r(n - e_s)
                              + 1{n_r > 0} 1{mu(n) > mu(n - e_r)} Phi(n - e_r) ]

    Given n every ordering carries the same class-weight factor, so the
    conditional law of the ordering is ``w(c)/Phi(n)`` and
    ``E[sir_r | n] = Xi_r(n)/Phi(n) =: g_r(n)``, a function of the count vector
    alone. The station mean then follows from the count marginal pM as
    ``E[sir_r] = sum_n pM(n) g_r(n)``, or in normalizing-constant form from the
    functional-server identity of :func:`pfqn_oi_fnc` applied to ``f(n) = g_r(n)``
    (note ``g_r(0) = 0``, as required).

    Parameters
    ----------
    oirate : function ``mu(n)`` returning the OI total service rate for the
        per-class count vector n (length R). ``mu(0)`` is taken as 0.
    N : (R,) closed population vector, finite.
    options : accepted for signature parity, currently unused.

    Returns
    -------
    g : (prod(N+1), R) table, column-major over the lattice 0 <= n <= N, with
        ``g[1 + sum(n * stride), r] = E[sir_r | n]``.
    Xi : (prod(N+1), R) table with the sir-weighted balance ``Xi_r(n)``.
    Phi : (prod(N+1),) table with the OI balance function ``Phi(n)``.

    See Also
    --------
    pfqn_oi_fnc, pfqn_ncoi, pfqn_mvaoi, pfqn_mvaoi_marg
    """
    if not callable(oirate):
        raise ValueError('oirate must be a callable mu(n).')
    N = np.asarray(N, dtype=int).ravel()
    if np.any(N < 0):
        raise ValueError('pfqn_oi_insvc requires finite nonnegative populations.')
    R = N.size
    shp = N + 1
    stride = np.ones(R, dtype=int)
    for d in range(1, R):
        stride[d] = stride[d - 1] * shp[d - 1]
    total = int(np.prod(shp))

    subs = np.zeros((total, R), dtype=int)
    for i in range(total):
        li = i
        for d in range(R):
            subs[i, d] = li % shp[d]
            li //= shp[d]
    muv = np.zeros(total)
    for i in range(total):
        if subs[i].sum() > 0:
            muv[i] = float(oirate(subs[i].copy()))

    Phi = np.zeros(total)
    Xi = np.zeros((total, R))
    for i in range(total):
        n = subs[i]
        if n.sum() == 0:
            Phi[i] = 1.0
            continue
        mun = muv[i]
        if mun <= 0:
            # Unreachable state (no server can serve this composition).
            continue
        sPhi = 0.0
        sXi = np.zeros(R)
        for s in range(R):
            if n[s] > 0:
                j = i - stride[s]
                sPhi += Phi[j]
                sXi += Xi[j, :]
        Phi[i] = sPhi / mun
        for r in range(R):
            acc = sXi[r]
            if n[r] > 0:
                j = i - stride[r]
                if mun > muv[j]:
                    acc += Phi[j]        # the tail class-r job is in service
            Xi[i, r] = acc / mun

    g = np.zeros((total, R))
    for i in range(total):
        if Phi[i] > 0:
            g[i, :] = Xi[i, :] / Phi[i]
    return g, Xi, Phi


def pfqn_oi_fnc(Phi: Sequence[float],
               N: Optional[Sequence[int]] = None,
               f: Optional[Callable] = None,
               options=None) -> Tuple[Callable, np.ndarray, np.ndarray]:
    """OI generalization of the load-dependent functional server.

    Builds an auxiliary OI station whose balance function Psi satisfies the
    convolution identity ``(Psi * Phi)(n) = (1 + f(n)) Phi(n)``, then inverts
    Psi to the FNC rate ``mu_f(n) = (sum_{r: n_r>0} Psi(n-e_r)) / Psi(n)``.

    Parameters
    ----------
    Phi : balance function of the existing OI station over the lattice, an
        R-dimensional array of shape (N_1+1, ..., N_R+1), or a flat
        column-major vector.
    N : (R,) closed population vector. Optional when Phi is a full
        R-dimensional array (then N = shape(Phi) - 1).
    f : target queue-dependent function f(n), f(0)=0 (default f = sum(n)).
    options : accepted for signature parity; unused.

    Returns
    -------
    (muf, Psi, mu) : callable rate handle (Inf outside the lattice), FNC
        balance array (lattice shape), and tabulated FNC rate array.
    """
    Phi_arr = np.asarray(Phi, dtype=float)

    if f is None:
        f = lambda n: float(np.sum(n))
    if N is None:
        N = np.asarray(Phi_arr.shape) - 1
    N = np.round(np.asarray(N, dtype=float)).astype(int).ravel()
    R = N.size
    shp = (N + 1).astype(int)

    # Flatten Phi in column-major (Fortran) order to match the MATLAB Phi(:).
    Phiv = Phi_arr.ravel(order='F')
    total = int(np.prod(shp))
    if Phiv.size != total:
        raise ValueError('numel(Phi) must equal prod(N+1).')

    # Column-major strides and decoded subscripts.
    stride = np.ones(R, dtype=int)
    for d in range(1, R):
        stride[d] = stride[d - 1] * shp[d - 1]
    subs = np.zeros((total, R), dtype=int)
    for i in range(total):
        li = i
        for d in range(R):
            subs[i, d] = li % shp[d]
            li //= shp[d]

    # Step 1: deconvolve (Psi * Phi)(n) = (1 + f(n)) Phi(n) for Psi.
    Psiv = np.zeros(total)
    for i in range(total):
        n = subs[i, :]
        acc = (1.0 + f(n)) * Phiv[i]
        for j in range(i):
            k = subs[j, :]
            if np.all(k <= n):
                idx = int(np.sum((n - k) * stride))  # linear index of n-k
                acc -= Psiv[j] * Phiv[idx]
        Psiv[i] = acc

    # Step 2: balanced-fairness inversion of Psi to the FNC rate.
    muv = np.full(total, np.inf)
    for i in range(total):
        n = subs[i, :]
        if np.all(n == 0):
            muv[i] = 0.0            # empty state
            continue
        denom = Psiv[i]
        if denom == 0:
            muv[i] = np.inf         # non-physical / undefined rate
            continue
        num = 0.0
        for r in range(R):
            if n[r] > 0:
                num += Psiv[i - stride[r]]
        muv[i] = num / denom

    if R == 1:
        Psi = Psiv.copy()
        mu = muv.copy()
    else:
        Psi = Psiv.reshape(tuple(shp), order='F')
        mu = muv.reshape(tuple(shp), order='F')

    muf = _make_oi_fnc_eval(muv, shp, stride)
    return muf, Psi, mu


def _make_oi_fnc_eval(muv: np.ndarray, shp: np.ndarray, stride: np.ndarray) -> Callable:
    """Return the tabulated-lattice rate handle muf(n) (Inf out of lattice)."""
    def oi_fnc_eval(n):
        n = np.round(np.asarray(n, dtype=float)).astype(int).ravel()
        if n.size != shp.size or np.any(n < 0) or np.any(n > shp - 1):
            return np.inf
        return float(muv[int(np.sum(n * stride))])
    return oi_fnc_eval
