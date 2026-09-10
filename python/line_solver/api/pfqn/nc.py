"""
Normalizing Constant methods for Product-Form Queueing Networks.

Native Python implementations of methods for computing normalizing constants:
- Convolution Algorithm (pfqn_ca)
- Related utility functions

References:
    Buzen, J.P. "Computational algorithms for closed queueing networks with
    exponential servers." Communications of the ACM 16.9 (1973): 527-531.
"""

import numpy as np
from math import log, exp, factorial, lgamma, ldexp
# largest log-scale value whose exponential is still a finite double
_LOG_DBL_MAX = log(np.finfo(float).max)
from typing import Tuple, Dict, Any
from functools import lru_cache
import itertools
from math import comb


def _factln(n: float) -> float:
    """Compute log(n!) using log-gamma function."""
    if n <= 0:
        return 0.0
    return lgamma(n + 1)


def _population_lattice_pprod(n: np.ndarray, N: np.ndarray = None) -> np.ndarray:
    """
    Generate next population vector in lexicographic order.

    Iterates through all population vectors from (0,0,...,0) to N.

    Args:
        n: Current population vector
        N: Maximum population per class (for bounds)

    Returns:
        Next population vector, or (-1,...,-1) when exhausted
    """
    R = len(n)
    n_next = n.copy()

    if N is None:
        # Just increment
        n_next[-1] += 1
        return n_next

    # Find rightmost position that can be incremented
    for i in range(R - 1, -1, -1):
        if n_next[i] < N[i]:
            n_next[i] += 1
            # Reset positions to the right
            for j in range(i + 1, R):
                n_next[j] = 0
            return n_next

    # All exhausted
    return -np.ones(R, dtype=int)


def _hashpop(n: np.ndarray, N: np.ndarray) -> int:
    """
    Compute linear index for population vector.

    Maps population vector n to unique integer index in
    the flattened population lattice [0, prod(N+1)).

    Args:
        n: Population vector
        N: Maximum population per class

    Returns:
        Linear index
    """
    R = len(n)
    idx = 0
    mult = 1
    for i in range(R - 1, -1, -1):
        idx += int(n[i]) * mult
        mult *= int(N[i]) + 1
    return idx


def _pfqn_pff_delay(Z: np.ndarray, n: np.ndarray) -> float:
    """
    Product-form factor for delay stations (think times).

    Computes contribution to normalizing constant from delay stations.

    Args:
        Z: Think times per class
        n: Population vector

    Returns:
        Product-form factor
    """
    R = len(n)
    if np.sum(n) == 0:
        return 1.0
    # Log-space accumulation (MATLAB pfqn_ca Fz): naive Z^n/n! overflows the
    # int->float conversion for n beyond ~170.
    f = 0.0
    for r in range(R):
        if Z[r] > 0:
            f += log(Z[r]) * n[r]
            f -= lgamma(1.0 + n[r])
        elif n[r] > 0:
            return 0.0
    return exp(f)


def pfqn_ca(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None
            ) -> Tuple[float, float]:
    """
    Convolution Algorithm for normalizing constant computation.

    Computes the normalizing constant G(N) for a closed product-form
    queueing network using Buzen's convolution algorithm.

    Args:
        L: Service demand matrix (M x R) where M is stations, R is classes
        N: Population vector (1 x R or R,) - number of jobs per class
        Z: Think time vector (1 x R or R,) - think time per class (default 0)

    Returns:
        Tuple (G, lG) where:
            - G: Normalizing constant
            - lG: log(G)
    """
    L = np.asarray(L, dtype=np.float64)
    N = np.asarray(N, dtype=np.float64).flatten()
    N = np.ceil(N).astype(int)

    R = len(N)

    if L.ndim == 1:
        L = L.reshape(-1, 1) if R == 1 else L.reshape(1, -1)

    M = L.shape[0]  # Number of stations

    # Handle Z
    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()

    # Special case: no stations (only delay)
    if M == 0:
        # G = prod_r (Z[r]^N[r] / N[r]!)
        lGn = 0.0
        for r in range(R):
            lGn += -_factln(N[r])  # -log(N[r]!)
            if Z[r] > 0 and N[r] > 0:
                lGn += N[r] * log(Z[r])
        Gn = exp(lGn)
        return Gn, lGn

    # Check for negative populations
    if np.any(N < 0):
        return 0.0, float('-inf')

    # Check for zero population
    if N.sum() == 0:
        return 1.0, 0.0

    # see _kb/03-api-layer.md for rationale
    Nt = int(N.sum())
    # Each class independently takes whichever station -- or the delay -- gives it
    # its largest factor. The mixed state so named has term at least the product of
    # those factors, because a station holding several classes carries a multinomial
    # coefficient of at least one, so this is still a LOWER bound on log G. It
    # dominates the per-configuration maximum it replaces, which asked ONE station
    # (or the delay) to hold every class at once and so dropped the delay entirely
    # as soon as a single class had no think time. That collapse is what made the
    # scaling scale UP: on L=[1e-9,1], N=[99,1], Z=[1,0] the old estimate was the
    # all-at-the-queue -2051.6 against a true log G of -359.1, giving kscale=-30,
    # and Z/2^-30 = 1.07e9 overflowed the delay column Z^n/n! at n=[40,0].
    lGest = 0.0
    for r in range(R):
        if N[r] <= 0:
            continue
        best = float('-inf')
        for i in range(M):
            if L[i, r] > 0:
                best = max(best, N[r] * log(L[i, r]))
        Zr = float(Z[r])
        if Zr > 0:
            best = max(best, N[r] * log(Zr) - _factln(N[r]))
        if not np.isfinite(best):
            # no station and no delay can hold class r, so G(N) is exactly zero
            lGest = float('-inf')
            break
        lGest += best
    if not np.isfinite(lGest):
        kscale = 0
    else:
        kscale = int(round(lGest / (Nt * log(2))))
    cscale = ldexp(1.0, kscale)  # exact power of two
    L = L / cscale
    Z = Z / cscale

    # Compute total number of population vectors
    product_N_plus_one = int(np.prod(N + 1))

    # G[m, idx] = G_m(n) where idx = hashpop(n, N)
    G = np.ones((M + 1, product_N_plus_one))

    # Iterate through all population vectors
    n = np.zeros(R, dtype=int)  # Start at (0, 0, ..., 0)

    while True:
        # Check if done (n becomes all -1)
        if np.all(n < 0):
            break

        idxn = _hashpop(n, N)

        # Base case: delay station contribution
        G[0, idxn] = _pfqn_pff_delay(Z, n)

        # Convolution recursion: G_m(n) = G_{m-1}(n) + sum_r L[m-1,r] * G_m(n - e_r)
        for m in range(1, M + 1):
            G[m, idxn] = G[m - 1, idxn]
            for r in range(R):
                if n[r] >= 1:
                    n[r] -= 1
                    idxn_1r = _hashpop(n, N)
                    n[r] += 1
                    G[m, idxn] += L[m - 1, r] * G[m, idxn_1r]

        # Next population vector
        n = _population_lattice_pprod(n, N)

    # see _kb/03-api-layer.md for rationale
    G_final = G[M, product_N_plus_one - 1]
    if G_final > 0:
        lGn = log(G_final) + Nt * kscale * log(2)
    else:
        lGn = float('-inf')
    # see _kb/03-api-layer.md for rationale
    # Undoing the scaling can leave the double range; lGn is then the only usable
    # form, so report the linear-scale G as infinite rather than overflowing ldexp.
    if lGn > _LOG_DBL_MAX:
        Gn = float('inf')
    else:
        Gn = float(np.ldexp(G_final, int(Nt * kscale)))

    return Gn, lGn


def pfqn_is(L, N, Z=None, options=None):
    """
    Importance-sampling (IS) estimate of the normalizing constant of a closed
    LOAD-INDEPENDENT product-form queueing network with M single-server queues of
    per-class demand L and an aggregated delay of think time Z.

    This is the load-independent case of :func:`pfqn_ld_is` (capacities
    mu_i(k)=1), and the ordinary-network counterpart of the order-independent
    pfqn_oi_is and the pass-and-swap pfqn_pas_is: all four are the same
    sample-an-ordering estimator, differing only in the per-position factor of
    each station's balance function. For a single-server queue that factor is the
    demand of the class at that position, L(i,q_p); for the delay it is Z(q_p)/p;
    for an OI/P&S station it is the reciprocal rank rate 1/mu_i(supp(q_1..q_p)).

    With ell = sum(N), an ordering c of all ell jobs is drawn by placing a
    uniformly random present class at each step (probability p(c) = product of the
    reciprocal branching factors), and the sum over ALL ways of cutting c into
    contiguous per-station segments is computed exactly by dynamic programming::

        G(N) = E_{C~p}[ S(C)/p(C) ],
        S(c) = sum_{cuts} prod_m prod_p L(m, seg_m(p))

    which is unbiased for the exact constant of :func:`pfqn_nc`.

    Parameters
    ----------
    L : (M, R) array
        Per-class service demands at the M single-server queues.
    N : (R,) array
        Closed population vector, finite.
    Z : (R,) array, optional
        Aggregated think time (delay) demand; None or zeros if none.
    options : dict or options object, optional
        Fields ``samples`` (default 1e4) and ``seed`` (optional).

    Returns
    -------
    (G, lG) : tuple of float
        IS estimate of the normalizing constant and its logarithm.

    Examples
    --------
    >>> L = np.array([[0.5, 0.3], [0.2, 0.4]]); N = np.array([3, 2]); Z = np.array([1.0, 1.0])
    >>> G, lG = pfqn_is(L, N, Z, {'samples': 100000})

    See Also
    --------
    pfqn_ld_is, pfqn_nc
    """
    from .ncld import pfqn_ld_is
    res = pfqn_ld_is(L, N, Z, None, options)
    return res.G, res.lG


def _nc_G(lG: float) -> float:
    """exp(lG) for the un-logged constant, saturating instead of raising.

    math.exp raises OverflowError above ~709, so a perfectly good lG (862 nats
    on a 12-station model with N = 1100) turned a correct result into an
    exception. G simply does not fit a double there; inf is the honest value and
    lG carries the answer.
    """
    if not np.isfinite(lG):
        return 0.0
    try:
        return exp(lG)
    except OverflowError:
        return float('inf')


def pfqn_nc_resolved_method(method: str) -> str:
    """The algorithm pfqn_nc actually runs for METHOD.

    pfqn_nc substitutes a different algorithm for some method names, so the
    requested name is not always the one that ran. Callers that report the
    method to the user must resolve it through here, otherwise the banner
    names an algorithm that never executed.

    Currently only 'comom' substitutes: the native pfqn_comomrm port is not
    numerically robust for R>1, so convolution is used instead, which is exact
    for the single-station product-form models CoMoM-RM targets. The
    substitution is unconditional, hence resolvable without solving. Note this
    is BROADER than MATLAB, whose pfqn_nc reports 'ca' only for the R==1 case
    and genuinely runs CoMoM-RM for R>1 with a single queue.
    """
    if method == 'comom':
        return 'ca'
    return method


def pfqn_nc(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
            method: str = 'ca', options=None) -> Tuple[float, float]:
    """
    Normalizing constant computation dispatcher.

    Selects appropriate algorithm based on method parameter.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,) (default: zeros)
        method: Algorithm to use:
            - 'ca', 'exact': Convolution algorithm
            - 'default': Auto-select based on problem size
            - 'le': Logistic expansion (Cas17 eq. 34 as published)
            - 'ble': Logistic expansion plus the empirical eps->0 correction (BLE)
            - 'aghq': adaptive Gauss-Hermite over the simplex; q=1 is 'le'
            - 'cub': Controllable upper bound
            - 'imci': Importance sampling Monte Carlo integration
            - 'pana': PANACEA asymptotic expansion (load-independent)
            - 'propfair': Proportionally fair allocation
            - 'rgf': Recursion by generating functions (grouped stations; multiclass
              by iterated residues, with think times)
            - 'mmint2': Gauss-Legendre quadrature
            - 'gleint': Gauss-Legendre integration
            - 'sampling': Monte Carlo sampling
            - 'kt': Knessl-Tier expansion
            - 'bkt': Knessl-Tier expansion minus the Stirling remainder of each Laplaced class (BKT)
            - 'lekt': the estimator 'ble' and 'bkt' both compute, on the cheaper side
            - 'comom': Conditional moments
            - 'rd': Reduction heuristic
            - 'ls': Logistic sampling
            - 'mcmc': Chen-O'Cinneide regularization (Markov chain Monte Carlo on the
              regularized network); supplies X and Q, not a constant of its own

    Returns:
        Tuple (G, lG) - normalizing constant and its log
    """
    method = method.lower() if method else 'ca'

    # Preprocessing matching MATLAB pfqn_nc:
    # 1. Convert inputs
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()

    # Early returns
    if np.any(N < 0) or len(N) == 0:
        return 0.0, float('-inf')
    if np.sum(N) == 0:
        return 1.0, 0.0

    if Z is None:
        Z = np.zeros(len(N))
    else:
        Z = np.asarray(Z, dtype=float).ravel()

    # 2. Erase open classes (inf population -> 0)
    N = np.where(np.isinf(N), 0.0, N)

    # 3. Remove zero-population classes
    nnz_classes = np.where(N > 0)[0]
    if len(nnz_classes) == 0:
        # All classes have zero population
        Z_sum = np.sum(Z)
        if Z_sum == 0:
            return 1.0, 0.0
        else:
            lG = float(-np.sum([_factln(int(N[r])) for r in range(len(N))]) +
                        np.sum([N[r] * log(max(Z[r], 1e-300)) for r in range(len(N)) if N[r] > 0]))
            return exp(lG) if np.isfinite(lG) else 0.0, lG

    L = L[:, nnz_classes]
    N = N[nnz_classes]
    Z = Z[nnz_classes]

    # 4. Scale demands to improve numerical stability
    M, R = L.shape
    scalevec = np.ones(R)
    for r in range(R):
        scalevec[r] = max(np.max(L[:, r]), Z[r]) if (np.max(L[:, r]) > 0 or Z[r] > 0) else 1.0
    L = L / scalevec
    Z = Z / scalevec

    # 5. Remove zero-demand stations
    Lsum = np.sum(L, axis=1)
    dem_stations = np.where(Lsum > 1e-12)[0]
    L = L[dem_stations, :]

    # Check degenerate: no demand at any station
    if L.size == 0 or np.sum(L) < 1e-12:
        Z_sum = np.sum(Z)
        if Z_sum < 1e-12:
            lG = 0.0
        else:
            lG = float(-np.sum([_factln(int(N[r])) for r in range(R)]) +
                        np.sum([N[r] * log(max(Z[r], 1e-300)) for r in range(R) if N[r] > 0]) +
                        np.dot(N, np.log(scalevec)))
        return exp(lG) if np.isfinite(lG) else 0.0, lG

    M, R = L.shape
    Ntot = int(np.sum(N))

    # The algorithm the 'default' branch settles on is not derivable from the
    # arguments alone, so it is recorded here and handed back through OPTIONS
    # rather than through the return tuple, whose arity callers depend on.
    resolved = {'method': method}

    # Dispatch to method
    def _compute_nc(L, N, Z, method):
        M, R = L.shape
        Ntot = int(np.sum(N))
        Z_row = Z  # Z is already 1D

        if method in ['ca', 'exact']:
            return pfqn_ca(L, N, Z_row)
        elif method == 'clw':
            # Choudhury-Leung-Whitt generating function inversion: each
            # single-server station is a multiplicity-1 queue, delay is the IS term
            return pfqn_clw(L, N, Z_row)
        elif method == 'ger':
            # Residue closed form of the same generating function 'clw' inverts
            # numerically. A class eliminated by residues enters only as a pole ORDER,
            # so its population is free: this is the cheap route when one population
            # dwarfs the others, and the expensive one when the classes are many, since
            # the term count grows as C(S+M-1,M-1) per further elimination. maxterms
            # REFUSES rather than truncating, so an oversized model errors here instead
            # of returning a wrong lG. The solver options are deliberately not
            # forwarded: pfqn_gerasimov's tol is a pole-merging threshold, not the
            # iterative tolerance options.tol carries.
            from .gerasimov import pfqn_gerasimov
            return pfqn_gerasimov(L, N, Z_row)
        elif method == 'pana':
            return pfqn_panacea(L, N, Z_row)
        elif method == 'propfair':
            G, lG, _ = pfqn_propfair(L, N, Z_row)
            return G, lG
        elif method == 'divdiff':
            # Divided-difference closed form, Casale (SIGMETRICS 2017), Eqs. (15)
            # and (16). Load-independent single-server queues only: a think time
            # needs the integral form of Corollary 3.4, which is not implemented.
            # Unlike the default route below this one keeps pfqn_explicit's
            # warnings, since a caller that named the method has no fallback.
            from .explicit import pfqn_explicit
            if float(np.sum(Z_row)) > 0:
                raise ValueError("pfqn_nc: the 'divdiff' method requires a model without think "
                                 "time, which needs the integral form of Corollary 3.4. Use "
                                 "'ca' or 'default'.")
            lG, G, expr, _ = pfqn_explicit(L, N)
            resolved['method'] = 'divdiff/' + expr
            return G, lG
        elif method == 'rgf':
            # Recursion by generating functions. Single class: one sequence per
            # group of identically loaded stations (Coury-Harrison 1997, Property 1).
            # Multiclass: iterated residues (Harrison-Coury 2002, Thm 1) with think
            # times carried by the Bertozzi-McKenna truncation, which neither RGF
            # paper has. That sum is ALTERNATING, so pfqn_rgfmc refuses when the
            # cancellation leaves no significant digits rather than returning a
            # confidently wrong lG; the exact convolution answers those, and says so.
            from .rgf import pfqn_rgf, pfqn_rgfmc
            if L.shape[1] > 1:
                try:
                    return pfqn_rgfmc(L, N, Z_row)
                except ValueError:
                    resolved['method'] = 'rgf/ca'
                    return pfqn_ca(L, N, Z_row)
            G, lG, _ = pfqn_rgf(L[:, 0], N[0], float(np.sum(Z_row)))
            return G, lG
        elif method == 'le':
            # pfqn_le returns (Gn, lGn): taking [0] fed Gn in as lGn and exp() overflowed
            from .asymptotic import pfqn_le
            result = pfqn_le(L, N, Z_row)
            lG = result[1] if isinstance(result, tuple) else result
            G = _nc_G(lG)
            return G, lG
        elif method == 'ble':
            # LE plus the empirical eps->0 correction; see _kb/03-api-layer.md
            from .asymptotic import pfqn_ble
            result = pfqn_ble(L, N, Z_row)
            lG = result[1] if isinstance(result, tuple) else result
            G = _nc_G(lG)
            return G, lG
        elif method == 'aghq':
            # adaptive Gauss-Hermite over the simplex; q=1 would be 'le'.
            # options.config.aghq_nodes overrides the node count; the rule costs
            # q^(M-1) evaluations, so the default stays small.
            from .asymptotic import pfqn_aghq
            aghq_nodes = 3
            cfg = getattr(options, 'config', None) if options is not None else None
            if cfg is None and isinstance(options, dict):
                cfg = options.get('config')
            if cfg is not None:
                q = cfg.get('aghq_nodes') if isinstance(cfg, dict) else getattr(cfg, 'aghq_nodes', None)
                if q is not None:
                    aghq_nodes = max(1, int(round(float(q))))
            result = pfqn_aghq(L, N, Z_row, aghq_nodes)
            lG = result[1] if isinstance(result, tuple) else result
            G = _nc_G(lG)
            return G, lG
        elif method in ['cub', 'gm']:
            from .asymptotic import pfqn_cub
            order = int(np.ceil((Ntot - 1) / 2))
            result = pfqn_cub(L, N, Z_row, order=order, atol=1e-8)
            lG = result[1] if isinstance(result, tuple) else result
            G = _nc_G(lG)
            return G, lG
        elif method == 'is':
            # see _kb/03-api-layer.md for rationale
            return pfqn_is(L, N, Z_row, options)
        elif method == 'mci':
            from .asymptotic import pfqn_mci
            # plain Monte Carlo integration: pfqn_mci defaults to the 'imci' tilt
            result = pfqn_mci(L, N, Z_row, variant='mci')
            lG = result[1] if isinstance(result, tuple) else result
            G = _nc_G(lG)
            return G, lG
        elif method == 'imci':
            from .asymptotic import pfqn_mci
            # pfqn_mci returns (G, lG, lZ): taking [0] fed G in as lG
            result = pfqn_mci(L, N, Z_row)
            lG = result[1] if isinstance(result, tuple) else result
            G = _nc_G(lG)
            return G, lG
        elif method in ['mmint2', 'gleint']:
            from .quadrature import pfqn_mmint2, pfqn_mmint2_gausslegendre
            if L.shape[0] > 1:
                # THE EMPTY CONSTANT IS THE REFERENCE'S ANSWER, NOT A WORKAROUND.
                # pfqn_nc.m:287 has this exact arm: on more than one queueing
                # station it warns, sets lG = [] and returns, and getAvg then
                # renders a table of ZEROS while reporting a completed analysis.
                # The warning is verbose-gated there and the return is NOT, which
                # is why line_warning (verbosity-aware) sits above an
                # unconditional return here.
                #
                # A zero table is a poor answer and the reference knows it. It is
                # still the right one to reproduce, for three reasons. It is what
                # MATLAB, the ground truth, does. The support gate no longer
                # OFFERS this pair -- nc_method_refusal refuses it for the report
                # -- so nobody reaches it through model.help() or SolverAUTO; only
                # a caller naming the method does, and that caller gets what the
                # reference gives. And the alternative is what this port did
                # until now: quadrature on a multi-station L returned a
                # PLAUSIBLE WRONG NUMBER (QLen 2.9037/0.6578/0.4385 against the
                # exact 1.5548/1.6109/0.8343), which is strictly worse because
                # nothing distinguishes it from an answer.
                #
                # So do not "fix" this back into a silent wrong number, and do not
                # turn it forward into a refusal either: the run path is ruled
                # (2026-07-25), and cpp/tests/test_nc.cpp pins the zero table.
                from ..io.logging import line_warning
                line_warning('pfqn_nc', 'The %s method requires a model with a delay '
                                        'and a single queueing station.' % method)
                return float('nan'), float('nan')
            if method == 'gleint':
                lG, _ = pfqn_mmint2_gausslegendre(L, N, Z_row)
            else:
                lG, _ = pfqn_mmint2(L, N, Z_row)
            G = _nc_G(lG)
            return G, lG
        elif method == 'sampling':
            from .quadrature import pfqn_mmsample2
            lG, _ = pfqn_mmsample2(L, N, Z_row)
            G = _nc_G(lG)
            return G, lG
        elif method == 'kt':
            from .kt import pfqn_kt
            # pfqn_kt returns (G, lG, X, Q): unpacking two names raised
            # ValueError, so the advertised method could never run
            G, lG = pfqn_kt(L, N, Z_row)[:2]
            return G, lG
        elif method == 'bkt':
            # KT minus the exact Stirling remainder of each Laplaced class; see _kb/03-api-layer.md
            from .kt import pfqn_bkt
            G, lG = pfqn_bkt(L, N, Z_row)
            return G, lG
        elif method == 'lekt':
            # the estimator ble and bkt both compute, on the cheaper side; see _kb/03-api-layer.md
            from .asymptotic import pfqn_lekt
            G, lG = pfqn_lekt(L, N, Z_row)
            return G, lG
        elif method == 'bk':
            from .bk import pfqn_bk
            G, lG = pfqn_bk(L, N, Z_row)[:2]
            return G, lG
        elif method == 'bkue':
            # The uniform expansion is single chain by construction; the
            # multichain fallback is the saddle point of the same paper, which
            # is also how the analyzer reaches this branch when it conditions
            # on a station population with an auxiliary class.
            from .bk import pfqn_bk, pfqn_bkue
            if L.shape[1] > 1:
                G, lG = pfqn_bk(L, N, Z_row)[:2]
            else:
                G, lG = pfqn_bkue(L, float(N[0]), float(np.sum(Z_row)))
            return G, lG
        elif method in ('lc', 'lc.ue'):
            # Algorithm 2 returns mean values, not a multichain constant; the
            # saddle point that seeds it supplies lG on the same asymptotics
            from .bk import pfqn_bk
            G, lG = pfqn_bk(L, N, Z_row)[:2]
            return G, lG
        elif method == 'comom':
            # see _kb/03-api-layer.md for rationale
            return pfqn_ca(L, N, Z_row)
        elif method == 'rd':
            from .rd import pfqn_rd
            result = pfqn_rd(L, N, Z_row)
            # pfqn_rd returns a tuple (lGN, Cgamma)
            lG = result[0] if isinstance(result, tuple) else result.lGN
            G = _nc_G(lG)
            return G, lG
        elif method == 'mcmc':
            # Chen-O'Cinneide REGULARIZATION (TOMACS 8(3), 1998) estimates the RATIOS
            # G(N-e_r)/G(N) and the queue lengths, never G, so it supplies no constant of
            # its own. lG here is the BLE expansion and is NOT part of the paper: it
            # cancels out of every mean value the analyzer reports, and only
            # getProbNormConstAggr reads it. The mean values themselves are taken from
            # pfqn_mcmc by the NC handler, which is where the simulation is run.
            if M > 1:
                from .asymptotic import pfqn_ble
                _, lG = pfqn_ble(L, N, Z_row)[:2]
            else:
                # MATLAB reaches CoMoM-RM here. This port substitutes convolution for the
                # same reason the 'comom' branch above does (see
                # pfqn_nc_resolved_method): the native pfqn_comomrm is not numerically
                # robust for R>1, and pfqn_ca is EXACT on the single-station models
                # CoMoM-RM targets, so the substitution can only improve this lG.
                return pfqn_ca(L, N, Z_row)
            return _nc_G(lG), lG
        elif method == 'ls':
            from .ls import pfqn_ls
            return pfqn_ls(L, N, Z_row)
        elif method == 'nrl':
            from .laplace import pfqn_nrl
            lG = pfqn_nrl(L, N, Z_row)
            G = _nc_G(lG)
            return G, lG
        elif method == 'nrp':
            from .laplace import pfqn_nrp
            lG = pfqn_nrp(L, N, Z_row)
            G = _nc_G(lG)
            return G, lG
        elif method == 'nre':
            from .nre import pfqn_nre
            lG = pfqn_nre(L, N, Z_row)
            G = _nc_G(lG)
            return G, lG
        elif method == 'default':
            from math import comb
            from .asymptotic import pfqn_cub, pfqn_le, pfqn_cub_evals, CUB_MAX_EVALS
            # Single class and no think time: the divided-difference closed form
            # of Casale (SIGMETRICS 2017), Eqs. (15) and (16), is exact, costs
            # O(M^2) at any population and is evaluated in the log domain, so it
            # neither truncates like the cubature past sum(N)~33 nor leaves the
            # range of a double like the convolution. Demands that are exactly
            # tied fix a multiplicity and cost nothing; demands that are close
            # but distinct spend digits on cancellation instead, so the closed
            # form is told to REFUSE rather than warn once it would burn half of
            # double precision, and the general route below picks those up.
            if M > 1 and R == 1 and float(np.sum(Z_row)) == 0:
                from .explicit import pfqn_explicit
                maxLossDigits = 8  # half the ~16 digits a double carries
                lGe, Ge, expr, _ = pfqn_explicit(L, N, None, 'auto', maxLossDigits)
                if np.isfinite(lGe):
                    resolved['method'] = 'divdiff/' + expr
                    return Ge, lGe
            if M > 1:
                order = -1
                if Ntot < 1000:
                    # CUB with order selection matching MATLAB cost budget
                    Cmax = M * R * (50 ** 3)
                    maxorder = min(int(np.ceil((Ntot - 1) / 2)), 16)
                    order = 0
                    totCost = 0
                    while order < maxorder:
                        nextCost = R * comb(M + 2 * (order + 1), M - 1)
                        if totCost + nextCost <= Cmax:
                            order += 1
                            totCost += nextCost
                        else:
                            break
                # Cmax above prices neither the Grundmann-Moeller node count nor
                # the think-time v-integration that repeats the whole rule, so the
                # order is re-priced here and lowered until it fits the budget.
                # Lowering the order keeps the cubature; switching to le instead
                # would hand these models to a Laplace expansion whose mode sits
                # on the simplex boundary when L has near-zero rows, which is
                # exactly the flat-layer case that trips the budget.
                while order > 0 and pfqn_cub_evals(M, order, Z_row) > CUB_MAX_EVALS:
                    order -= 1
                if order >= 0:
                    result = pfqn_cub(L, N, Z_row, order=order, atol=1e-8)
                    lG = result[1] if isinstance(result, tuple) else result
                    G = _nc_G(lG)
                    resolved['method'] = 'cub'
                    return G, lG
                else:
                    # BLE on the default path: strictly better on lG and it
                    # cancels in G(N-e_r)/G(N). 'le' stays the published form.
                    from .asymptotic import pfqn_ble
                    result = pfqn_ble(L, N, Z_row)
                    lG = result[1] if isinstance(result, tuple) else result
                    G = _nc_G(lG)
                    resolved['method'] = 'ble'
                    # Birman-Kogan Algorithm 2 supplies the MEAN VALUES here.
                    # The caller's fallback differences lG at R+M*R reduced
                    # populations, which on many stations is both dearer and
                    # ~300x less accurate than the load concealment fixed point. Gated
                    # on the station count, since load concealment is mean
                    # field in M: see _kb/06-solver-catalog.md. The fixed point
                    # itself runs in the NC handler, which owns X and Q.
                    if M >= 10 and R > 1 and np.all(N >= 0) and np.sum(N) > 0:
                        resolved['method'] = 'ble/lc'
                    return G, lG
            elif M == 1:
                Z_sum = np.sum(Z_row)
                if Z_sum < 1e-12:
                    # Single queue, no delay: exact formula
                    lG = float(-np.dot(N, np.log(L[0, :])))
                    G = _nc_G(lG)
                    resolved['method'] = 'exact'
                    return G, lG
                else:
                    if Ntot < 10000:
                        resolved['method'] = 'ca'
                        return pfqn_ca(L, N, Z_row)
                    else:
                        from .asymptotic import pfqn_ble
                        result = pfqn_ble(L, N, Z_row)
                        lG = result[1] if isinstance(result, tuple) else result
                        G = _nc_G(lG)
                        resolved['method'] = 'ble'
                        return G, lG
            else:
                resolved['method'] = 'ca'
                return pfqn_ca(L, N, Z_row)
        else:
            # An unrecognized method name must not alias onto exact convolution: the
            # caller reports the analysis under the requested method name, so
            # the substitution is invisible in the result.
            raise ValueError("pfqn_nc: unrecognized method: %s" % method)

    G, lG = _compute_nc(L, N, Z, method)
    if isinstance(options, dict):
        options['resolved_method'] = resolved['method']

    # Scale back: lG += N * log(scalevec)
    lG = lG + float(np.dot(N, np.log(scalevec)))
    G = _nc_G(lG)
    return G, lG


def pfqn_panacea(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
                 terms: int = 3) -> Tuple[float, float]:
    """
    PANACEA asymptotic expansion for load-independent closed networks.

    McKenna-Mitra normal-usage expansion whose coefficients are linear
    combinations of pseudonetwork partition functions evaluated by convolution.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,)
        terms: Number of terms in the normal-usage asymptotic series
            (1, 2, or 3; default 3), as selectable in the original PANACEA
            package (Ramakrishnan-Mitra, BSTJ 61(10):2849-2872, 1982)

    Returns:
        Tuple (G, lG) - normalizing constant and its log, both NaN when the
        model is not in normal usage
    """
    if terms not in (1, 2, 3):
        raise ValueError("The terms parameter must be 1, 2, or 3 "
                         "(higher-order coefficients are not implemented).")
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    q, p = L.shape
    if Z is None or np.size(Z) == 0:
        Z = np.full(p, 1e-8)
    else:
        Z = np.asarray(Z, dtype=float)
        Z = np.sum(Z, axis=0) if Z.ndim > 1 else Z.flatten()

    if L.size == 0 or np.all(np.sum(L, axis=0) == 0):
        lGn = -float(np.sum([_factln(n) for n in N])) + float(np.sum(N * np.log(Z)))
        return exp(lGn), lGn

    r = L / np.tile(Z, (q, 1))
    Nt = np.max(1.0 / r[r > 0])   # ignore structural zeros
    beta = N / Nt
    gamma = r * Nt
    alpha = 1.0 - N.dot(r.T)
    if np.min(alpha) < 0:
        return float('nan'), float('nan')
    gammatilde = gamma / np.tile(alpha.reshape(-1, 1), (1, p))

    A0 = 1.0
    A1 = 0.0
    if terms >= 2:
        for j in range(p):
            m = np.zeros(p)
            m[j] = 2
            A1 -= beta[j] * pfqn_ca(gammatilde, m)[0]

    A2 = 0.0
    if terms >= 3:
        for j in range(p):
            m = np.zeros(p)
            m[j] = 3
            A2 += 2 * beta[j] * pfqn_ca(gammatilde, m)[0]
            m = np.zeros(p)
            m[j] = 4
            A2 += 3 * beta[j] ** 2 * pfqn_ca(gammatilde, m)[0]
            for k in range(p):
                if k != j:
                    m = np.zeros(p)
                    m[j] = 2
                    m[k] = 2
                    A2 += 0.5 * beta[j] * beta[k] * pfqn_ca(gammatilde, m)[0]

    I = [A0, A1 / Nt, A2 / Nt ** 2][:terms]
    sumI = sum(I)
    if sumI <= 0:
        return float('nan'), float('nan')

    lGn = (-float(np.sum([_factln(n) for n in N]))
           + float(np.sum(N * np.log(Z))) + log(sumI)
           - float(np.sum(np.log(alpha))))
    if not np.isfinite(lGn):
        return float('nan'), float('nan')
    return exp(lGn), lGn


def pfqn_propfair(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None
                  ) -> Tuple[float, float, np.ndarray]:
    """
    Proportionally Fair allocation approximation for normalizing constant.

    Estimates the normalizing constant using a convex optimization program
    that is asymptotically exact in models with single-server PS queues only.

    This method is based on Schweitzer's approach and Walton's proportional
    fairness theory for multi-class networks.

    Args:
        L: Service demand matrix (M x R) where M is stations, R is classes
        N: Population vector (1 x R or R,) - number of jobs per class
        Z: Think time vector (1 x R or R,) - think time per class (default 0)

    Returns:
        Tuple (G, lG, X) where:
            - G: Estimated normalizing constant
            - lG: log(G)
            - X: Asymptotic throughputs per class (1 x R)

    References:
        Schweitzer, P. J. (1979). Approximate analysis of multiclass closed networks
        of queues. In Proceedings of the International Conference on Stochastic
        Control and Optimization.

        Walton, N. (2009). Proportional fairness and its relationship with
        multi-class queueing networks.
    """
    from scipy.optimize import minimize

    L = np.asarray(L, dtype=np.float64)
    N = np.asarray(N, dtype=np.float64).flatten()

    R = len(N)

    if L.ndim == 1:
        L = L.reshape(-1, 1) if R == 1 else L.reshape(1, -1)

    M = L.shape[0]  # Number of stations

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()

    FineTol = 1e-12

    # Objective function: maximize sum_r (N[r] - x[r]*Z[r]) * log(x[r])
    # We minimize the negative
    def objective(x):
        obj = 0.0
        for r in range(R):
            obj += (N[r] - x[r] * Z[r]) * log(abs(x[r]) + FineTol)
        return -obj  # Minimize negative

    # Constraints: sum_r L[i,r] * x[r] <= 1 for all stations i
    # And x[r] >= 0 for all r
    constraints = []

    # Capacity constraints
    for i in range(M):
        constraints.append({
            'type': 'ineq',
            'fun': lambda x, i=i: 1.0 - sum(L[i, r] * x[r] for r in range(R))
        })

    # Non-negativity constraints
    for r in range(R):
        constraints.append({
            'type': 'ineq',
            'fun': lambda x, r=r: x[r]
        })

    # Initial guess - balanced throughput
    x0 = np.zeros(R)
    for r in range(R):
        D_max = L[:, r].max() if M > 0 else 0
        if D_max > 0:
            x0[r] = min(N[r] / (Z[r] + 1), 1.0 / D_max)
        elif Z[r] > 0:
            x0[r] = N[r] / Z[r]
        else:
            x0[r] = N[r]

    # Ensure positive initial guess
    x0 = np.maximum(x0, FineTol)

    # Run optimization with COBYLA
    result = minimize(
        objective,
        x0,
        method='COBYLA',
        constraints=constraints,
        options={'maxiter': 10000, 'rhobeg': 1.0}
    )

    Xasy = result.x

    # Compute lG
    lG = 0.0
    for r in range(R):
        x = Xasy[r]
        if x > FineTol:
            lG += (N[r] - x * Z[r]) * log(1.0 / (x + FineTol))

    # Factorial correction for think times
    for r in range(R):
        thinking = Xasy[r] * Z[r]
        if thinking > 0:
            lG -= _factln(thinking)

    G = exp(lG) if lG > -700 else 0.0  # Avoid underflow

    Xa = Xasy.reshape(1, -1)

    return G, lG, Xa


def pfqn_ls(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
            I: int = 100000) -> Tuple[float, float]:
    """
    Logistic sampling approximation for normalizing constant.

    Approximates the normalizing constant using importance sampling from
    a multivariate normal distribution fitted at the leading eigenvalue mode.

    This method is particularly effective for large networks where
    convolution becomes computationally expensive.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,) (default: zeros)
        I: Number of samples for Monte Carlo integration (default: 100000)

    Returns:
        Tuple (G, lG) where:
            G: Estimated normalizing constant
            lG: log(G)

    Reference:
        G. Casale. "Accelerating performance inference over closed systems by
        asymptotic methods." ACM SIGMETRICS 2017.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()

    # Filter out zero-demand stations
    Lsum = np.sum(L, axis=1)
    L = L[Lsum > 1e-4, :]

    M, R = L.shape

    # Handle empty network
    if L.size == 0 or np.sum(L) < 1e-4 or N.size == 0 or np.sum(N) == 0:
        # Per-class Z, as MATLAB's sum(Z,1) is, and an empty class contributes 0
        # rather than 0*log(0). Z may also be absent altogether on this branch.
        Zt = np.zeros(len(N)) if Z is None else np.asarray(Z, dtype=float).ravel()
        lGn = -np.sum([_factln(n) for n in N])
        for r in range(len(N)):
            if N[r] > 0:
                lGn += N[r] * np.log(Zt[r] if r < Zt.size else 0.0)
        return np.exp(lGn), lGn

    if Z is None or len(Z) == 0:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float).flatten()

    # Find the mode using fixed-point iteration
    u, converged = _pfqn_le_fpi(L, N, Z)
    if not converged:
        # Fall back to convolution if mode-finding fails
        return pfqn_ca(L, N, Z)

    Ntot = np.sum(N)

    if np.sum(Z) <= 0:
        # Case without think times
        # Compute Hessian at the mode
        A = _pfqn_le_hessian(L, N, u)
        A = (A + A.T) / 2  # Ensure symmetry

        try:
            iA = np.linalg.inv(A)
        except np.linalg.LinAlgError:
            return pfqn_ca(L, N, Z)

        x0 = np.log(u[:M-1] / u[M-1])

        # Sample from multivariate normal
        samples = np.random.multivariate_normal(x0, iA, I)

        # Evaluate function at samples
        T = np.zeros(I)
        for i in range(I):
            T[i] = _simplex_fun(samples[i, :], L, N)

        # Evaluate PDF at samples
        dpdf = np.zeros(I)
        for i in range(I):
            diff = samples[i, :] - x0
            dpdf[i] = np.exp(-0.5 * diff @ np.linalg.inv(iA) @ diff) / np.sqrt((2 * np.pi) ** (M-1) * np.linalg.det(iA))

        # Compute normalizing constant
        valid = dpdf > 0
        if np.sum(valid) == 0:
            return pfqn_ca(L, N, Z)

        lGn = _multinomialln(np.append(N, M-1)) + _factln(M-1) + np.log(np.mean(T[valid] / dpdf[valid]))
        Gn = np.exp(lGn)

    else:
        # Case with think times Z > 0
        u, v, converged = _pfqn_le_fpiZ(L, N, Z)
        if not converged:
            return pfqn_ca(L, N, Z)

        # Compute Hessian
        A = _pfqn_le_hessianZ(L, N, Z, u, v)
        A = (A + A.T) / 2

        try:
            iA = np.linalg.inv(A)
        except np.linalg.LinAlgError:
            return pfqn_ca(L, N, Z)

        x0 = np.append(np.log(u[:M-1] / u[M-1]), np.log(v))

        # Sample from multivariate normal
        samples = np.random.multivariate_normal(x0, iA, I)

        # Evaluate function at samples
        epsilon = 1e-10
        eN = epsilon * np.sum(N)
        eta = np.sum(N) + M * (1 + eN)
        K = M

        T = np.zeros(I)
        for i in range(I):
            x = samples[i, :]
            term1 = -np.exp(x[K-1]) + K * (1 + eN) * x[M-1]
            term2 = 0.0
            for r in range(R):
                inner = L[K-1, r] * np.exp(x[K-1]) + Z[r]
                for k in range(K-1):
                    inner += np.exp(x[k]) * (L[k, r] * np.exp(x[K-1]) + Z[r])
                term2 += N[r] * np.log(np.maximum(inner, 1e-300))
            term3 = np.sum(x[:K-1])
            term4 = -eta * np.log(1 + np.sum(np.exp(x[:K-1])))
            T[i] = np.exp(term1 + term2 + term3 + term4)

        # Evaluate PDF at samples
        dpdf = np.zeros(I)
        for i in range(I):
            diff = samples[i, :] - x0
            try:
                dpdf[i] = np.exp(-0.5 * diff @ np.linalg.inv(iA) @ diff) / np.sqrt((2 * np.pi) ** len(x0) * np.linalg.det(iA))
            except:
                dpdf[i] = 0

        valid = dpdf > 0
        if np.sum(valid) == 0:
            return pfqn_ca(L, N, Z)

        Gn = np.exp(-np.sum([lgamma(1 + n) for n in N])) * np.mean(T[valid] / dpdf[valid])
        lGn = np.log(Gn) if Gn > 0 else float('-inf')

    return Gn, lGn


def _pfqn_le_fpi(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None):
    """
    Fixed-point iteration to find mode of Gaussian approximation.

    Returns:
        (u, converged): Mode vector and convergence flag
    """
    M, R = L.shape
    Ntot = np.sum(N)
    u = np.ones(M) / M
    u_prev = np.ones(M) * np.inf

    for iteration in range(1000):
        u_prev = u.copy()
        for i in range(M):
            u[i] = 1 / (Ntot + M)
            for r in range(R):
                denom = np.dot(u_prev, L[:, r])
                if denom > 0:
                    u[i] += N[r] / (Ntot + M) * L[i, r] * u_prev[i] / denom

        if np.linalg.norm(u - u_prev, 1) < 1e-10:
            return u, True

    return u, False


def _pfqn_le_fpiZ(L: np.ndarray, N: np.ndarray, Z: np.ndarray):
    """
    Fixed-point iteration with think times.

    Returns:
        (u, v, converged): Mode vector, scale factor, and convergence flag
    """
    M, R = L.shape
    eta = np.sum(N) + M
    u = np.ones(M) / M
    # eq. (35) in the SIGMETRICS 2017 paper has a spurious +1 in the v equation;
    # the correct stationary point is v = eta - sum_r xi_r*Z_r.
    v = eta

    for iteration in range(1000):
        u_prev = u.copy()
        v_prev = v

        for ist in range(M):
            u[ist] = 1 / eta
            for r in range(R):
                denom = Z[r] + v * np.dot(u_prev, L[:, r])
                if denom > 0:
                    u[ist] += (N[r] / eta) * (Z[r] + v * L[ist, r]) * u_prev[ist] / denom

        xi = np.zeros(R)
        for r in range(R):
            denom = Z[r] + v * np.dot(u_prev, L[:, r])
            if denom > 0:
                xi[r] = N[r] / denom

        v = eta - np.dot(xi, Z)

        if np.linalg.norm(u - u_prev, 1) + abs(v - v_prev) < 1e-10:
            return u, v, True

    return u, v, False


def _pfqn_le_hessian(L: np.ndarray, N: np.ndarray, u: np.ndarray):
    """
    Compute Hessian of Gaussian approximation (without think times).
    """
    M, R = L.shape
    Ntot = np.sum(N)
    hu = np.zeros((M-1, M-1))

    for i in range(M-1):
        for j in range(M-1):
            if i != j:
                hu[i, j] = -(Ntot + M) * u[i] * u[j]
                for r in range(R):
                    denom = np.dot(u, L[:, r]) ** 2
                    if denom > 0:
                        hu[i, j] += N[r] * L[i, r] * L[j, r] * u[i] * u[j] / denom
            else:
                sum_other = np.sum(u) - u[i]
                hu[i, j] = (Ntot + M) * u[i] * sum_other
                for r in range(R):
                    denom = np.dot(u, L[:, r]) ** 2
                    L_other = np.sum(L[:, r]) - L[i, r]
                    if denom > 0:
                        hu[i, j] -= N[r] * L[i, r] * u[i] * (sum_other * L_other) / denom

    return hu


def _pfqn_le_hessianZ(L: np.ndarray, N: np.ndarray, Z: np.ndarray, u: np.ndarray, v: float):
    """
    Compute Hessian of Gaussian approximation (with think times).
    """
    K, R = L.shape
    Ntot = np.sum(N)
    A = np.zeros((K, K))

    csi = np.zeros(R)
    for r in range(R):
        denom = Z[r] + v * np.dot(u, L[:, r])
        if denom > 0:
            csi[r] = N[r] / denom

    Lhat = np.zeros((K, R))
    for k in range(K):
        for r in range(R):
            Lhat[k, r] = Z[r] + v * L[k, r]

    eta = Ntot + K

    for i in range(K):
        for j in range(K):
            if i != j:
                A[i, j] = -eta * u[i] * u[j]
                for r in range(R):
                    if N[r] > 0:
                        A[i, j] += csi[r]**2 * Lhat[i, r] * Lhat[j, r] * u[i] * u[j] / N[r]

    for i in range(K):
        A[i, i] = -np.sum(A[i, :]) + A[i, i]

    # Reduce to (K-1) x (K-1) and add v column
    A_reduced = A[:K-1, :K-1]
    A_result = np.zeros((K, K))
    A_result[:K-1, :K-1] = A_reduced

    A_result[K-1, K-1] = 1
    for r in range(R):
        if N[r] > 0:
            A_result[K-1, K-1] -= (csi[r]**2 / N[r]) * Z[r] * np.dot(u, L[:, r])
    A_result[K-1, K-1] *= v

    for i in range(K-1):
        A_result[i, K-1] = 0
        for r in range(R):
            if N[r] > 0:
                A_result[i, K-1] += v * u[i] * ((csi[r]**2 / N[r]) * Lhat[i, r] * np.dot(u, L[:, r]) - csi[r] * L[i, r])
        A_result[K-1, i] = A_result[i, K-1]

    return A_result


def _simplex_fun(x: np.ndarray, L: np.ndarray, N: np.ndarray) -> float:
    """
    Evaluate simplex function for LS algorithm.
    """
    M = len(x) + 1
    v = np.zeros(M)
    for i in range(M-1):
        v[i] = np.exp(x[i])
    v[M-1] = 1

    term1 = np.sum(N * np.log(np.dot(v, L)))
    term2 = np.sum(x)
    term3 = -(np.sum(N) + M) * np.log(np.sum(v))

    return np.exp(term1 + term2 + term3)


def _multinomialln(n: np.ndarray) -> float:
    """
    Compute log of multinomial coefficient.
    """
    return _factln(np.sum(n)) - np.sum([_factln(ni) for ni in n])


def _clw_euler_weights(n: int, mm: int) -> np.ndarray:
    """
    Euler weights of eq. (2.22): E(m,n) = sum_i w_i (-1)^i a_i.

    E(m,n) = 2^-m sum_{k=0}^{m} C(m,k) S_{n+k} with S_t = sum_{i<=t} (-1)^i a_i,
    so a_i carries the mass of every partial sum that contains it.
    """
    b = np.empty(mm + 1)
    b[0] = 1.0
    for k in range(1, mm + 1):
        b[k] = b[k - 1] * (mm - k + 1) / k         # C(mm,k)
    b = b * (2.0 ** (-mm))
    tail = np.cumsum(b[::-1])[::-1]                # tail[k] = 2^-mm sum_{j>=k} C(mm,j)
    w = np.ones(n + mm + 1)
    w[n + 1:] = tail[1:mm + 1]
    return w


def _clw_components(adj: np.ndarray, mask: np.ndarray):
    """Connected components of adj with the nodes in mask removed."""
    p = adj.shape[0]
    lab = -np.ones(p, dtype=int)
    nc = 0
    for s in range(p):
        if mask[s] or lab[s] >= 0:
            continue
        lab[s] = nc
        stack = [s]
        while stack:
            v = stack.pop()
            for u in range(p):
                if adj[v, u] and not mask[u] and lab[u] < 0:
                    lab[u] = nc
                    stack.append(u)
        nc += 1
    return [np.where(lab == c)[0] for c in range(nc)]


def _clw_dimred(L: np.ndarray, enabled: bool, maxd: int):
    """
    Interdependence graph of the factors of (4.5) and the subset D minimizing
    the inversion dimension |D| + max_i |S_i(D)| (eqs. 3.1-3.3), by enumeration
    in increasing cardinality.
    """
    p = L.shape[1]
    if not enabled or p <= 2:
        return np.zeros(0, dtype=int), [np.arange(p)]
    adj = np.zeros((p, p), dtype=bool)
    for i in range(L.shape[0]):
        v = np.where(L[i, :] != 0)[0]
        if v.size:
            adj[np.ix_(v, v)] = True               # each factor is a clique (eq. 3.1)
    np.fill_diagonal(adj, False)
    bestC = _clw_components(adj, np.zeros(p, dtype=bool))
    best = max(len(c) for c in bestC)              # |D| = 0
    bestD = np.zeros(0, dtype=int)
    for dd in range(1, min(maxd, p - 1) + 1):
        if dd >= best:
            break                                  # dimension is at least |D|
        if comb(p, dd) > 2e5:
            break                                  # (3.3) is solved by enumeration only
        for sub in itertools.combinations(range(p), dd):
            mask = np.zeros(p, dtype=bool)
            mask[list(sub)] = True
            cc = _clw_components(adj, mask)
            mx = max((len(c) for c in cc), default=0)
            if dd + mx < best:
                best = dd + mx
                bestD = np.array(sub, dtype=int)
                bestC = cc
    if best < p:
        return bestD, bestC
    return np.zeros(0, dtype=int), [np.arange(p)]


def pfqn_clw(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
             m: np.ndarray = None, l: np.ndarray = None,
             gamma: np.ndarray = None, euler: bool = True, euler_n: int = 11,
             euler_m: int = 20, euler_tol: float = 1e-10, euler_maxm: int = 160,
             beta: np.ndarray = None, dimred: bool = True,
             dimred_maxd: int = 4) -> Tuple[float, float]:
    """
    Choudhury-Leung-Whitt normalization constant by numerical inversion of the
    generating function (JACM 42(5):935-970, 1995).

    Computes g(K) of a multichain closed product-form network with single-server
    and (optionally) infinite-server queues by numerically inverting its
    p-dimensional generating function (eq. 4.5)

        G(z) = exp(sum_j rho_{j0} z_j) / prod_i (1 - sum_j rho_{ji} z_j)^{m_i}

    where j=1..p indexes chains, i=1..q' the distinct single-server queues with
    multiplicity m_i. g(K) is recovered by nested one-dimensional lattice-Poisson
    inversions (eq. 2.3) with restrictive static scaling (eqs. 5.41-5.46) and
    log-domain recovery (eq. 7.1).

    Both of the paper's accelerations are applied. Dimension reduction by
    decomposition (Sec. 3, Sec. 5.4) removes from the interdependence graph of
    the factors the subset D minimizing |D| + max_i |S_i(D)| (eq. 3.3): with the
    D variables fixed on their contours the remaining factors share no variable,
    so each connected component is inverted separately and the results
    multiplied. Euler summation (Sec. 2.4, eq. 2.22) replaces the nearly
    alternating inner sum of (2.3) by the Euler sum of its first n+m+1 terms,
    applied once for k >= 0 and once for k < 0, so prod_j K_j becomes
    prod_j min(n+m+1, K_j) in the cost (eq. 2.26); the order m is doubled until
    the paper's own estimate |E(m,n) - E(m,n+1)| falls under euler_tol.

    Args:
        L: (q' x p) single-server relative traffic intensities, L[i,j]=rho_{ji}.
        N: (p,) closed-chain population vector K.
        Z: (p,) aggregate infinite-server relative intensities rho_{j0}. Default 0.
        m: (q',) queue multiplicities m_i. Default ones.
        l: (p,) inner lattice parameters l_j indexed by chain. Default by
           inversion depth: 1 at depth 1, 2 at depths 2-3, 3 deeper.
        gamma: (p,) aliasing parameters gamma_j. Default by depth: 11, 13, 13, 15.
        euler: apply Euler summation where K_j > euler_n + euler_m.
        euler_n: terms summed exactly before averaging (n in eq. 2.22).
        euler_m: starting order of the Euler averaging (m in eq. 2.22).
        euler_tol: relative tolerance on |E(m,n) - E(m,n+1)|.
        euler_maxm: largest Euler order reached by doubling.
        beta: (p,) multipliers on the scale parameters alpha_j, the manual tuning
           of page 956 (the paper uses 0.8 <= beta <= 1.2 on its largest
           examples). Default ones.
        dimred: apply dimension reduction by decomposition.
        dimred_maxd: largest |D| examined when minimizing (3.3).

    Returns:
        Tuple (G, lG): normalization constant (inf if it overflows double) and
        its natural logarithm (always finite).
    """
    L = np.asarray(L, dtype=np.float64)
    if L.ndim == 1:
        L = L.reshape(-1, 1)
    qd, p = L.shape
    N = np.round(np.asarray(N, dtype=np.float64).flatten()).astype(int)
    if Z is None:
        Z = np.zeros(p)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()
    if m is None:
        m = np.ones(qd)
    else:
        m = np.asarray(m, dtype=np.float64).flatten()

    if beta is None:
        beta = np.ones(p)
    else:
        beta = np.asarray(beta, dtype=np.float64).flatten()

    if np.any(N < 0):
        return 0.0, -np.inf
    if np.all(N == 0):
        return 1.0, 0.0

    # dimension reduction (Section 3): D is inverted first, then each connected
    # component of the interdependence graph minus D, independently
    ordD, comps = _clw_dimred(L, dimred, dimred_maxd)
    d = ordD.size
    order = np.concatenate([ordD] + comps).astype(int)

    # depth of each chain: D occupies depths 1..d and every component restarts at
    # depth d+1, since components are inverted in parallel
    depth = np.zeros(p, dtype=int)
    depth[ordD] = np.arange(1, d + 1)
    for c in comps:
        depth[c] = d + np.arange(1, c.size + 1)

    if l is None:
        l = np.where(depth == 1, 1, np.where(depth <= 3, 2, 3))
    else:
        l = np.round(np.asarray(l, dtype=np.float64).flatten()).astype(int)
    l = np.asarray(l, dtype=int)
    if gamma is None:
        gamma = np.where(depth == 1, 11.0, np.where(depth <= 3, 13.0, 15.0))
    else:
        gamma = np.asarray(gamma, dtype=np.float64).flatten()

    # contour radii r_j = 10^{-gamma_j/(2 l_j K_j)} (eq. 2.7)
    r = np.ones(p)
    for j in range(p):
        if N[j] > 0:
            r[j] = 10.0 ** (-gamma[j] / (2 * l[j] * N[j]))

    # restrictive static scaling (eqs. 5.41-5.46), outer vars at |z_k| = r_k.
    # The loop runs in inversion order, which is what dimension reduction changes
    # (Section 5.4); chains of distinct components never deflate one another,
    # because they share no queue.
    alpha = np.ones(p)
    used = np.zeros(qd)
    eta = (L != 0).astype(float)
    for t in range(p):
        j = int(order[t])
        Kj = int(N[j])
        lj = int(l[j])
        if Kj == 0:
            continue        # empty chain: no lattice, and 2*lj*Kj = 0 below
        denom = 1.0 - used
        denom[denom <= 0] = np.finfo(float).eps
        e = L[:, j] / denom
        posq = np.where(L[:, j] > 0)[0]
        aj = np.inf
        if posq.size > 0:
            srt = np.argsort(-e[posq])
            qs = posq[srt]
            es = e[qs]
            ms = m[qs]
            cumrho = np.cumsum(es) / np.arange(1, es.size + 1)
            cummb = np.cumsum(ms)
            inner = order[t + 1:]
            for n in range(es.size):
                qi = qs[n]
                # N_{ij} = mbar_n - 1 + sum_{k inside j} K_k eta_{k,qi} (eq. 5.43)
                Nn = int(round(cummb[n] - 1 + np.sum(N[inner] * eta[qi, inner])))
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
        if Z[j] > 0:
            aj = min(aj, Kj / Z[j])
        if not np.isfinite(aj):
            aj = 1.0
        alpha[j] = beta[j] * aj          # page 956: manual tuning of alpha_j
        used = used + alpha[j] * L[:, j] * r[j]

    arho0 = alpha * Z          # (p,)
    rhoS = L * alpha           # (q' x p)
    chunk = 2000000

    # split the factors of (4.5) over D and the components: a queue whose chains
    # all lie in D is constant during the component inversions, and every other
    # queue has all of its non-D chains inside a single component
    inD = np.zeros(p, dtype=bool)
    inD[ordD] = True
    compOf = np.zeros(p, dtype=int)
    for c in range(len(comps)):
        compOf[comps[c]] = c
    qBucket = -np.ones(qd, dtype=int)
    for i in range(qd):
        v = np.where(L[i, :] != 0)[0]
        v = v[~inD[v]]
        if v.size:
            qBucket[i] = compOf[v[0]]
    qD = np.where(qBucket == -1)[0]
    qC = [np.where(qBucket == c)[0] for c in range(len(comps))]

    # per-group normalization (Section 2.2, page 944): the scaling normalizes the
    # whole generating function, not each group, and a decomposition multiplies
    # the groups together. Every factor has nonnegative coefficients, so a group's
    # modulus is maximized at w = r; that constant cancels in the recovery and is
    # left at zero on the undecomposed path, which is thus unchanged.
    def group_bound(queues, chains):
        off = 0.0
        if chains.size:
            off += float(np.sum(arho0[chains] * (r[chains] - 1.0)))
        if queues.size:
            pole = 1.0 - rhoS[queues, :] @ r
            pole = np.maximum(pole, np.finfo(float).tiny)
            off -= float(np.sum(m[queues] * np.log(pole)))
        return off

    offD = 0.0
    offC = np.zeros(len(comps))
    if d > 0 or len(comps) > 1:
        offD = group_bound(qD, ordD)
        offC = np.array([group_bound(qC[c], comps[c]) for c in range(len(comps))])

    def gbar_sub(W, queues, chains, off=0.0):
        # Gbar_S(w) = exp(sum_{j in chains} arho0_j (w_j-1))
        #             / prod_{i in queues} (1 - sum_j rhoS_ij w_j)^{m_i}
        if chains.size:
            expo = (W[:, chains] - 1.0) @ arho0[chains]
        else:
            expo = np.zeros(W.shape[0])
        if queues.size:
            A = W @ rhoS[queues, :].T
            logden = np.log(1.0 - A) @ m[queues]
        else:
            logden = np.zeros(W.shape[0])
        return np.exp(expo - logden - off)

    def evalpts(j, kk, Kj, lj, k1, rj, W, fn, vec):
        theta = np.pi * (k1 + lj * kk) / (lj * Kj)
        wj = rj * np.exp(1j * theta)
        nk = wj.size
        out = np.empty(nk, dtype=complex)
        if vec:
            for a in range(0, nk, chunk):
                b = min(a + chunk, nk)
                Wn = np.repeat(W, b - a, axis=0)
                Wn[:, j] = wj[a:b]
                out[a:b] = fn(Wn)
        else:
            for t in range(nk):
                Wn = W.copy()
                Wn[:, j] = wj[t]
                out[t] = fn(Wn)
        return out

    def inner_sum(j, Kj, lj, k1, rj, W, fn, vec):
        # the sum splits at k = 0 into two nearly alternating series (Section 2.4);
        # each is replaced by its Euler sum, refined until it stops moving
        mCur = euler_m
        while True:
            T = euler_n + mCur
            if (not euler) or Kj <= T + 1:
                kk = np.arange(-Kj, Kj)
                v = evalpts(j, kk, Kj, lj, k1, rj, W, fn, vec)
                return np.sum(((-1.0) ** kk) * v)
            kk = np.arange(-(T + 2), T + 2)
            v = evalpts(j, kk, Kj, lj, k1, rj, W, fn, vec)
            vpos = v[T + 2:]                   # k =  s     , s = 0..T+1
            vneg = v[T + 1::-1]                # k = -(s+1) , s = 0..T+1
            w1 = _clw_euler_weights(euler_n, mCur)
            w2 = _clw_euler_weights(euler_n + 1, mCur)
            sg = (-1.0) ** np.arange(T + 2)
            dv = vpos - vneg
            E1 = np.sum(sg[:T + 1] * w1 * dv[:T + 1])
            E2 = np.sum(sg * w2 * dv)
            if abs(E1 - E2) <= euler_tol * abs(E2) or mCur >= euler_maxm:
                return E2
            mCur *= 2

    def lattice(j, W, fn, vec):
        Kj = int(N[j])
        lj = int(l[j])
        rj = r[j]
        if Kj == 0:
            # [w_j^0] Gbar = Gbar(w_j=0): the K=0 lattice is the single point 0,
            # and exp(-arho0_j) there cancels the +arho0_j added back into lG
            Wn = W.copy()
            Wn[:, j] = 0.0
            v = fn(Wn)
            return np.sum(v) if vec else v
        acc = 0.0 + 0.0j
        for k1 in range(lj):
            ph = np.exp(-1j * np.pi * k1 / lj)
            acc += ph * inner_sum(j, Kj, lj, k1, rj, W, fn, vec)
        return acc / (2 * lj * Kj * rj ** Kj)

    def invert_c(c, s, W):
        vars_ = comps[c]
        if s >= vars_.size:
            return gbar_sub(W, qC[c], vars_, offC[c])[0]
        j = int(vars_[s])
        if s == vars_.size - 1:
            val = lattice(j, W, lambda Wn: gbar_sub(Wn, qC[c], vars_, offC[c]), True)
        else:
            val = lattice(j, W, lambda Wn: invert_c(c, s + 1, Wn), False)
        if d == 0 and s == 0:
            # with D empty every component is an independent subnetwork, so its
            # coefficient is real
            val = val.real
        return val

    def invert_d(t, W):
        if t >= d:
            # D fixed: the remaining factors have no variable in common, so the
            # coefficient of the inner monomial is the product of the components'
            val = gbar_sub(W, qD, ordD, offD)[0]
            for c in range(len(comps)):
                val = val * invert_c(c, 0, W)
            return val
        j = int(ordD[t])
        val = lattice(j, W, lambda Wn: invert_d(t + 1, Wn), False)
        if t == 0:
            val = val.real
        return val

    gbar = invert_d(0, np.zeros((1, p), dtype=complex))
    lG = np.log(gbar) + offD + float(np.sum(offC)) + np.sum(arho0) - np.sum(N * np.log(alpha))
    G = np.inf if lG > 709 else np.exp(lG)
    return float(G), float(lG)


def pfqn_clw_lld(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
                 mu: np.ndarray = None, l: np.ndarray = None,
                 gamma: np.ndarray = None) -> Tuple[float, float]:
    """
    Choudhury-Leung-Whitt normalization constant by numerical inversion of the
    generating function (JACM 42(5):935-970, 1995), extended to limited
    load-dependent (LLD) stations via the per-center transforms of Bertozzi
    and McKenna (SIAM Review 35(2):239-268, 1993).

    The generating function is (Bertozzi-McKenna eqs. 2.17/2.23)

        G(z) = exp(sum_j rho_{j0} z_j) prod_i F_i(sum_j rho_{ji} z_j)

    where F_i is the transform of the station factor of queue i (eq. 2.16)
    with load-dependent rate scalings S_i(k) = mu[i,k]. For an LLD queue,
    S_i(k) = c_i constant for k >= l_i, and F_i is the rational function
    (eq. 2.19)

        F_i(x) = [c_i + sum_{n=1}^{l_i-1} (c_i - S_i(n))
                        / prod_{k=1}^n S_i(k) * x^n] / (c_i - x),

    analytic except for a simple pole at x = c_i. Multiserver and
    load-independent queues are special cases. Since g(K) depends on S_i(k)
    only for k <= sum(K), general load-dependent input is truncated to LLD
    at sum(K) without loss of exactness.

    g(K) is recovered by p nested one-dimensional lattice-Poisson inversions
    (CLW eq. 2.3) with restrictive static scaling adapted from CLW eqs.
    5.41-5.46 (each queue normalized by its pole c_i, simple pole) and
    log-domain recovery (eq. 7.1).

    Args:
        L: (q' x p) single-server relative traffic intensities, L[i,j]=rho_{ji}.
        N: (p,) closed-chain population vector K.
        Z: (p,) aggregate infinite-server relative intensities rho_{j0}. Default 0.
        mu: (q' x n) load-dependent rate scalings mu[i,k] = S_i(k+1); if fewer
            than sum(N) columns are given the last column is extended (LLD
            assumption). Default ones (all queues load-independent).
        l: (p,) inner lattice parameters l_j. Default 1,2,2,3,3,...
        gamma: (p,) aliasing parameters gamma_j. Default 11,13,13,15,15,...

    Returns:
        Tuple (G, lG): normalization constant (inf if it overflows double) and
        its natural logarithm (always finite).

    Note: cost is prod_j 2 l_j K_j contour points, each of cost O(sum_i l_i);
    practical for moderate populations and few chains.
    """
    L = np.asarray(L, dtype=np.float64)
    if L.ndim == 1:
        L = L.reshape(-1, 1)
    qd, p = L.shape
    N = np.round(np.asarray(N, dtype=np.float64).flatten()).astype(int)
    if Z is None:
        Z = np.zeros(p)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()
    ntot = int(np.sum(N))
    if mu is None:
        mu = np.ones((qd, max(ntot, 1)))
    else:
        mu = np.asarray(mu, dtype=np.float64)
        if mu.ndim == 1:
            mu = mu.reshape(qd, -1)
    if l is None:
        l = np.full(p, 3)
        l[0] = 1
        if p >= 2:
            l[1] = 2
        if p >= 3:
            l[2] = 2
    else:
        l = np.round(np.asarray(l, dtype=np.float64).flatten()).astype(int)
    l = np.asarray(l, dtype=int)
    if gamma is None:
        gamma = np.full(p, 15.0)
        gamma[0] = 11
        if p >= 2:
            gamma[1] = 13
        if p >= 3:
            gamma[2] = 13
    else:
        gamma = np.asarray(gamma, dtype=np.float64).flatten()

    if np.any(N < 0):
        return 0.0, -np.inf
    if np.all(N == 0):
        return 1.0, 0.0

    # extend/truncate mu to sum(N) columns (LLD extension of last column)
    if mu.shape[1] < ntot:
        mu = np.hstack([mu, np.tile(mu[:, -1:], (1, ntot - mu.shape[1]))])
    else:
        mu = mu[:, :ntot]
    if np.any(mu <= 0):
        raise ValueError('pfqn_clw_lld: load-dependent rates mu[i,k] must be positive.')

    # drop zero-population chains: the coefficient of z_j^0 equals the pgf
    # restricted to z_j = 0, so chain j is removed exactly
    keep = N > 0
    L = L[:, keep]
    N = N[keep]
    Z = Z[keep]
    l = l[keep]
    gamma = gamma[keep]
    p = N.size

    # pole c_i and LLD cutoff l_i of each queue: S_i(k) = c_i for k >= l_i
    cpole = mu[:, -1].copy()
    numc = []
    for i in range(qd):
        mism = np.where(mu[i, :] != cpole[i])[0]
        li = int(mism[-1]) + 2 if mism.size > 0 else 1
        a = np.zeros(li)
        a[0] = cpole[i]
        if li > 1:
            cp = np.cumprod(mu[i, :li - 1])   # prod_{k=1}^n S_i(k)
            a[1:] = (cpole[i] - mu[i, :li - 1]) / cp
        numc.append(a)

    # contour radii r_j = 10^{-gamma_j/(2 l_j K_j)} (CLW eq. 2.7)
    r = 10.0 ** (-gamma / (2 * l * N))

    # see _kb/03-api-layer.md for rationale
    Lt = L / cpole[:, None]
    alpha = np.ones(p)
    used = np.zeros(qd)
    eta = (L != 0).astype(float)
    for j in range(p):
        Kj = int(N[j])
        lj = int(l[j])
        if Kj == 0:
            continue        # empty chain: no lattice, and 2*lj*Kj = 0 below
        denom = 1.0 - used
        denom[denom <= 0] = np.finfo(float).eps
        e = Lt[:, j] / denom
        posq = np.where(Lt[:, j] > 0)[0]
        aj = np.inf
        if posq.size > 0:
            order = np.argsort(-e[posq])
            qs = posq[order]
            es = e[qs]
            cumrho = np.cumsum(es) / np.arange(1, es.size + 1)
            for n in range(es.size):
                qi = qs[n]
                # N_{ij} = n - 1 + sum_{k>j} K_k eta_{k,qi} (eq. 5.43, m_i = 1)
                Nn = int(round(n + np.sum(N[j + 1:p] * eta[qi, j + 1:p])))
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
        if Z[j] > 0:
            aj = min(aj, Kj / Z[j])
        if not np.isfinite(aj):
            aj = 1.0
        alpha[j] = aj
        used = used + aj * Lt[:, j] * r[j]

    arho0 = alpha * Z          # (p,)
    rhoS = L * alpha           # (q' x p)
    chunk = 2000000

    def gbar_eval(W):
        # see _kb/03-api-layer.md for rationale
        expo = (W - 1.0) @ arho0
        X = W @ rhoS.T
        logf = np.zeros(W.shape[0], dtype=complex)
        for i in range(qd):
            xi = X[:, i]
            a = numc[i]
            num = np.full(xi.shape, a[-1], dtype=complex)   # Horner on N_i(x)
            for k in range(a.size - 2, -1, -1):
                num = num * xi + a[k]
            logf += np.log(num) - np.log(cpole[i] - xi)
        return np.exp(expo + logf)

    def invert(j, wfixed):
        Kj = int(N[j])
        lj = int(l[j])
        rj = r[j]
        if Kj == 0:
            # [w_j^0] Gbar = Gbar(w_j=0): the K=0 lattice is the single point 0,
            # and exp(-arho0_j) there cancels the +arho0_j added back into lG
            if j == p - 1:
                W = np.zeros((1, p), dtype=complex)
                if j > 0:
                    W[:, :j] = wfixed
                val = complex(np.sum(gbar_eval(W)))
            else:
                val = invert(j + 1, np.concatenate([wfixed, [0.0 + 0.0j]]))
            return val.real if j == 0 else val
        kk = np.arange(-Kj, Kj)
        signs = (-1.0) ** kk
        acc = 0.0 + 0.0j
        for k1 in range(lj):
            ph = np.exp(-1j * np.pi * k1 / lj)
            theta = np.pi * (k1 + lj * kk) / (lj * Kj)
            wj = rj * np.exp(1j * theta)
            if j == p - 1:
                inner = 0.0 + 0.0j
                nk = wj.size
                for a in range(0, nk, chunk):
                    b = min(a + chunk, nk)
                    W = np.empty((b - a, p), dtype=complex)
                    if j > 0:
                        W[:, :j] = wfixed
                    W[:, j] = wj[a:b]
                    inner += np.sum(signs[a:b] * gbar_eval(W))
            else:
                inner = 0.0 + 0.0j
                for t in range(wj.size):
                    inner += signs[t] * invert(j + 1, np.concatenate([wfixed, [wj[t]]]))
            acc += ph * inner
        val = acc / (2 * lj * Kj * rj ** Kj)
        if j == 0:
            val = val.real
        return val

    gbar = invert(0, np.array([], dtype=complex))
    # No group offsets here: the lld form keeps ONE contour over all stations,
    # so there is no decomposition bound to add back (matlab pfqn_clw_lld.m:211).
    lG = np.log(gbar) + np.sum(arho0) - np.sum(N * np.log(alpha))
    G = np.inf if lG > 709 else np.exp(lG)
    return float(G), float(lG)


def pfqn_perm(A: np.ndarray, m: np.ndarray = None) -> float:
    """
    Permanent of a demand matrix, with optional column multiplicities.

    The pfqn_ entry point of the permanent library. It exists because the
    product-form joint queue-length probability of the per-station TOTAL
    populations is a permanent of the demand matrix replicated once per job,
    which is a normalizing-constant quantity rather than a general-purpose
    linear algebra one; see pfqn_jointmarg.

    Orientation is chosen before repeated lines are grouped. That is a
    correctness concern, not an optimisation: perm(A) is transpose-invariant
    but the Ryser sum is not, and exploiting repeated rows silently expands the
    transpose.

    Args:
        A: Square matrix, or the matrix of distinct columns when m is given
        m: Multiplicity of each column of A; sum(m) must be the order of the
           expanded matrix

    Returns:
        The permanent value; 1.0 for the empty matrix

    References:
        H. J. Ryser, "Combinatorial Mathematics", Carus Mathematical
        Monographs 14, Mathematical Association of America, 1963.
    """
    from ..perm import compute_permanent

    A = np.atleast_2d(np.asarray(A, dtype=float))
    if A.shape[0] == 0 and (A.shape[1] == 0 or m is None):
        return 1.0

    if m is not None:
        m = np.asarray(m, dtype=int).ravel()
        if m.size != A.shape[1]:
            raise ValueError(
                "pfqn_perm: the multiplicity vector has %d entries but A has %d columns."
                % (m.size, A.shape[1]))
        cols = [A[:, j] for j in range(A.shape[1]) for _ in range(int(m[j]))]
        A = np.column_stack(cols) if cols else np.zeros((A.shape[0], 0))
        if A.shape[0] == 0 and A.shape[1] == 0:
            return 1.0

    if A.shape[0] != A.shape[1]:
        raise ValueError("pfqn_perm: the matrix must be square, got %d by %d."
                         % (A.shape[0], A.shape[1]))
    return float(compute_permanent(A))


__all__ = [
    'pfqn_ca',
    'pfqn_nc',
    'pfqn_panacea',
    'pfqn_propfair',
    'pfqn_ls',
    'pfqn_clw',
    'pfqn_clw_lld',
    'pfqn_perm',
]
