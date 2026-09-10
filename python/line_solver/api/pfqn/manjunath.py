"""
Exact normalizing constant of a closed multiclass product-form network whose
state space carries arbitrary linear integer constraints.

Native Python port of matlab/src/api/pfqn/pfqn_manjunath.m. This is the
queueing-network half of the transform technique of Manjunath and Sikdar, of
which `lossn_manjunath` is the loss-network half. The two solve the same problem
-- sum a product form over an irregular integer state space -- from opposite
ends of the paper: `lossn_manjunath` implements Section 2.2, a set of '<=' rows
over the Poisson terms nu^n/n!, while this routine implements Section 3 together
with Section 5.3, a MIXED set of '=', '<=' and '>' rows over the BCMP terms,
where the population constraint of a closed network is itself one of the
equalities.

THE MODEL. M queueing stations (FCFS, PS or LCFS; rows of L) and Mz delay
stations (rows of Z) serve R closed classes with populations N. Writing n_ir for
the class r jobs at station i and n_i = sum_r n_ir, the BCMP product form is

    p(n) = (1/G) prod_{i queueing} n_i! prod_r L_ir^{n_ir}/n_ir!
                 prod_{i delay}         prod_r Z_ir^{n_ir}/n_ir!

Every state obeys the R population equalities sum_i n_ir = N_r; on top of those
the caller may impose any number of further rows

    sum_{i,r} A[j, i + S*r] n_ir  {=, <=, >}  b[j],   S = M + Mz,

i.e. A acts on n.ravel(order='F'), the (M+Mz, R) occupancy read column by column
with the queueing stations first. With no extra rows the routine returns exactly
the normalizing constant of `pfqn_ca`, which is the parity oracle used by the
tests; with extra rows it answers a question no other routine in the pfqn family
can, the convolution and MVA recursions having nowhere to carry a second
constraint.

WHY THE GENERATING FUNCTION IS A PRODUCT, AND WHERE THE n_i! GOES. Marking class
r by z_r and row j by y_j, and abbreviating the monomial one class r job at
station i contributes as u_ir = z_r prod_j y_j^{A[j, i + S*r]}, the sum over the
occupancies of a single QUEUEING station is, by the multinomial theorem,

    sum_{n_i.} n_i! prod_r (L_ir u_ir)^{n_ir}/n_ir!
      = sum_k (sum_r L_ir u_ir)^k = 1 / (1 - sum_r L_ir u_ir),

so the n_i! that couples the classes is exactly what turns the station's factor
from an exponential into a geometric one. The paper reaches the same place
through the Euler integral n! = int_0^inf e^-t t^n dt (Eqns 16-18), which is that
geometric series evaluated; the closed form is used here because there is then no
quadrature to discretize. A DELAY station has no n_i! and keeps its exponential.
Hence

    F(z,y) = prod_i 1/(1 - sum_r L_ir u_ir) * prod_k prod_r exp(Z_kr u_kr)

and G is read off F as a coefficient: degree exactly N_r in z_r, and for row j
the degree its sense dictates -- exactly b_j for '=', the sum of degrees 0..b_j
for '<=' (the multiplier (y^{b+1}-1)/(y-1) of Eqn 5, whose residue is that
partial sum), and the complement of the latter for '>' (Eqn 6).

WHY IT IS A COEFFICIENT COMPUTATION AND NOT A QUADRATURE. The contour integrals
of Eqn 9 all have their only pole at the origin, of order one more than the
right-hand side, so each is a residue and hence a Taylor coefficient. The routine
therefore never integrates: it carries F as a multivariate power series truncated
at degree N_r in z_r and b_j in y_j. Truncation is exact because A is nonnegative
-- no monomial above a cut can be brought back down by a later factor.

Each queueing station is applied by SOLVING (1 - sum_r L_ir u_ir) x = ser rather
than by expanding the geometric series, which keeps the cost at one pass. Every
monomial of the operator raises the total class degree, so sweeping the lattice
in increasing total class degree lets each coefficient read only coefficients
already final: a Gauss-Seidel sweep whose result is the exact solve, not an
iterate. A delay station has no such recurrence and is convolved with exp term by
term, which is where its extra factor of the population in the cost comes from.

THE ELIMINATION ORDER IS THE MEMORY BOUND. Variable y_j is created when the first
station its row touches is multiplied in and discharged immediately after the
last, so peak memory is prod_r (N_r+1) times the product of (b_j+1) over the
SIMULTANEOUSLY LIVE rows, not over all rows. A row constraining one station
therefore costs essentially nothing. The class axes are live throughout, so
prod_r (N_r+1) is a floor -- the same lattice `pfqn_ca` walks.

SCOPE. Load-dependent and multiserver stations are NOT covered: their
per-station term is not geometric, and while the paper admits an arbitrary
f_i(n_i) in the single-class case (Section 2), the multiclass n_i! coupling used
above then breaks. Use `pfqn_gld` or `pfqn_conwayms` for those. A and b must be
integer valued and A nonnegative, since the residue argument counts whole units;
a fractional entry is refused rather than rounded.

Reference: D. Manjunath and B. Sikdar, Integral Expressions for the Numerical
Evaluation of Product Form Expressions Over Irregular Multidimensional Integer
Spaces. Sections 3 and 5.3.
"""

import math
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

__all__ = ['pfqn_manjunath', 'PfqnManjunathStats']


@dataclass
class PfqnManjunathStats:
    """
    Per-class solution of the constrained closed network, all (R,) arrays.

    `Q + think + blocked == N` exactly: a refused admission is a DELETED
    transition, so a blocked job never leaves the delay, and because the think
    time is exponential a held job is indistinguishable from one still thinking.
    Little's law is what separates the two.
    """
    Q: np.ndarray        # mean class r jobs at the queueing station
    X: np.ndarray        # class r cycle throughput
    U: np.ndarray        # class r utilization of the queueing station
    think: np.ndarray    # class r jobs genuinely thinking, X_r * Z_r
    blocked: np.ndarray  # class r jobs held at the delay by the constraint
    delay: np.ndarray    # class r jobs at the delay, think + blocked


def _empty_stats(R):
    z = np.zeros(R)
    return PfqnManjunathStats(Q=z.copy(), X=z.copy(), U=z.copy(), think=z.copy(),
                              blocked=z.copy(), delay=z.copy())


def _ret(G, lG, peak, stats, R):
    """Early exit, carrying an all-zero decomposition when one was asked for."""
    if not stats:
        return G, lG, peak
    return G, lG, peak, _empty_stats(R)


def _stats(L, N, Z, A, b, sense, lG, M, Mz, S, R):
    """
    Per-class decomposition, for the one configuration in which the truncated
    product form is the exact stationary law: a single queueing station inside
    the region and a single delay station outside it. See `pfqn_manjunath`.
    """
    if Mz != 1:
        raise ValueError("pfqn_manjunath: the per-class decomposition needs exactly one "
                         "delay station, got %d. Pass Z as a 1xR row of think times" % Mz)
    if M != 1:
        raise ValueError(
            "pfqn_manjunath: the per-class decomposition needs exactly one queueing "
            "station, got %d. With two or more the delay->q1->q2->delay cycle makes the "
            "chain irreversible, Kelly truncation no longer holds, and the truncated "
            "product form is not the stationary law (measured at 131%% error). G and lG "
            "are still returned and still correct as a sum over the admissible set" % M)
    # The delay must sit OUTSIDE the region: its columns are S*r + (S-1).
    for r in range(R):
        if np.any(A[:, S * r + (S - 1)] != 0):
            raise ValueError(
                "pfqn_manjunath: constraint row(s) reference the delay station in class %d "
                "(column %d). The delay must lie OUTSIDE the finite capacity region, "
                "because the decomposition charges every held job to it"
                % (r + 1, S * r + (S - 1)))

    # Ratios of normalizing constants are taken in the LOG domain, so the
    # internal power-of-two rescaling cancels without ever being reconstructed.
    Q = np.zeros(R)
    X = np.zeros(R)
    for r in range(R):
        qcol = S * r                       # column of (queueing station, class r)

        # Throughput. One class r job removed from the queue leaves a state of
        # population N - e_r whose admission rule is shifted by that job's own
        # requirement column, exactly as the loss network's g(C - A e_r).
        Nr = np.array(N, dtype=np.int64)
        Nr[r] -= 1
        if Nr[r] >= 0:
            _, lGr, _ = pfqn_manjunath(L, Nr, Z, A, b - A[:, qcol], sense)
            X[r] = math.exp(lGr - lG) if np.isfinite(lGr) else 0.0

        # Mean queue length from the marginal law. An '=' row is discharged by
        # picking a single coefficient, so each call returns the mass of exactly
        # that occupancy.
        for k in range(1, int(N[r]) + 1):
            row = np.zeros((1, S * R), dtype=np.int64)
            row[0, qcol] = 1
            Ak = np.vstack([A, row]) if A.size else row
            bk = np.concatenate([b, [k]])
            _, lGk, _ = pfqn_manjunath(L, N, Z, Ak, bk, sense + 'E')
            if np.isfinite(lGk):
                Q[r] += k * math.exp(lGk - lG)

    U = X * np.asarray(L[0, :], dtype=np.float64)      # one server, one visit
    think = X * np.asarray(Z[0, :], dtype=np.float64)  # Little's law at the delay
    delay = np.asarray(N, dtype=np.float64) - Q        # all that is not at the queue
    blocked = delay - think                            # the remainder is held there
    return PfqnManjunathStats(Q=Q, X=X, U=U, think=think, blocked=blocked, delay=delay)


def _factln(n):
    return math.lgamma(n + 1.0)


def _matlab_round(x):
    """MATLAB's round: half away from zero. Python's round() is half to EVEN, so
    it disagrees on every half-integer, and the reference is MATLAB."""
    return math.floor(x + 0.5) if x >= 0 else math.ceil(x - 0.5)


def _lattice(dims, R):
    """
    Column-major strides and total class degree of every lattice point, plus the
    subscript table used to test whether a shift stays in range.
    """
    dims = np.asarray(dims, dtype=np.int64)
    q = dims.size
    P = int(np.prod(dims))
    sub = np.zeros((P, q), dtype=np.int64)
    rep = 1
    for k in range(q):
        col = np.repeat(np.arange(dims[k], dtype=np.int64), rep)
        sub[:, k] = np.tile(col, P // (rep * int(dims[k])))
        rep *= int(dims[k])
    stride = np.ones(q, dtype=np.int64)
    for k in range(1, q):
        stride[k] = stride[k - 1] * dims[k - 1]
    lev = sub[:, :R].sum(axis=1)
    return sub, stride, lev


def _shifts(A, sub, stride, i, S, R, J, dims):
    """
    Flat offset and in-range mask of the monomial one class r job at station i
    contributes: z_r gains one degree and y_j gains A[j, i + S*r]. A row not live
    at this station has a zero entry here by construction of first/last, so a
    dead axis is never shifted.
    """
    offs = []
    oks = []
    for r in range(R):
        delta = np.zeros(R + J, dtype=np.int64)
        delta[r] = 1
        for j in range(J):
            delta[R + j] = A[j, i + S * r]
        if np.any(delta > np.asarray(dims, dtype=np.int64) - 1):
            offs.append(None)          # a single job already breaks the cut
            oks.append(None)
            continue
        offs.append(int(delta @ stride))
        oks.append(np.all(sub >= delta, axis=1))
    return offs, oks


def _expand(ser, dims, k, newdim):
    """Create marker k, keeping the existing content at degree zero: nothing
    multiplied in so far carries any power of it."""
    pre = int(np.prod(dims[:k]))
    post = int(np.prod(dims[k + 1:]))
    grown = np.zeros((pre, newdim, post), dtype=np.float64)
    grown[:, 0, :] = ser.reshape((pre, 1, post), order='F')[:, 0, :]
    return grown.reshape(-1, order='F')


def _reduce(ser, dims, k, sense_j, rhs):
    """Discharge marker k. The multiplier (y^{b+1}-1)/(y-1) of a '<=' row turns
    its residue into the partial sum of the coefficients of degrees 0..b, and the
    multiplier 1/y^{b+1} of an '=' row picks the single coefficient of degree b."""
    pre = int(np.prod(dims[:k]))
    dk = int(dims[k])
    post = int(np.prod(dims[k + 1:]))
    T = ser.reshape((pre, dk, post), order='F')
    if sense_j == 'E':
        return np.ascontiguousarray(T[:, rhs, :]).reshape(-1, order='F')
    return T.sum(axis=1).reshape(-1, order='F')


def _series(L, Z, N, A, b, sense):
    """
    Coefficient-domain evaluation of the multiple contour integral of Eqn 9 for a
    set of '=' and '<=' rows. The series is a truncated multivariate polynomial
    whose first R axes are the class markers z_r, of extent N_r+1 throughout, and
    whose remaining J axes are the row markers y_j, of extent 1 while row j is not
    live and b_j+1 while it is.
    """
    M = L.shape[0]
    Mz = Z.shape[0]
    S = M + Mz
    R = len(N)
    J = len(b)

    # Row j is created at the first station it touches and discharged after the
    # last, so only an induced width of rows is ever live. A row that reached
    # here touches at least one station, the trivial ones having been decided.
    first = np.zeros(J, dtype=np.int64)
    last = np.zeros(J, dtype=np.int64)
    for j in range(J):
        touched = [i for i in range(S) if any(A[j, i + S * r] != 0 for r in range(R))]
        first[j] = touched[0]
        last[j] = touched[-1]

    dims = list(np.asarray(N, dtype=np.int64) + 1) + [1] * J
    ser = np.zeros(int(np.prod(dims)), dtype=np.float64)
    ser[0] = 1.0
    peak = ser.size
    sub, stride, lev = _lattice(dims, R)

    for i in range(S):
        for j in np.nonzero(first == i)[0]:
            newdim = int(b[j]) + 1
            ser = _expand(ser, dims, R + int(j), newdim)
            dims[R + int(j)] = newdim
            peak = max(peak, ser.size)
            sub, stride, lev = _lattice(dims, R)

        offs, oks = _shifts(A, sub, stride, i, S, R, J, dims)

        if i < M:
            # Queueing station: solve (1 - sum_r L_ir u_ir) x = ser in place.
            # Every monomial of the operator raises the total class degree by
            # one, so a sweep in increasing total class degree reads only final
            # coefficients and the sweep IS the solve.
            coef = L[i, :]
            maxlev = int(lev.max()) if lev.size else 0
            for ell in range(1, maxlev + 1):
                at = np.nonzero(lev == ell)[0]
                if at.size == 0:
                    continue
                for r in range(R):
                    if coef[r] == 0.0 or offs[r] is None:
                        continue
                    sel = at[oks[r][at]]
                    if sel.size == 0:
                        continue
                    ser[sel] += coef[r] * ser[sel - offs[r]]
        else:
            # Delay station: no n_i! coupling, so the factor is a product of
            # exponentials, one per class, each convolved in term by term. There
            # is no first-order recurrence to exploit here, which is why the
            # delay costs a factor of the population that the queueing station
            # does not.
            coef = Z[i - M, :]
            for r in range(R):
                if coef[r] == 0.0 or offs[r] is None or N[r] == 0:
                    continue
                nxt = ser.copy()
                term = ser
                okr = np.nonzero(oks[r])[0]
                for n in range(1, int(N[r]) + 1):
                    shifted = np.zeros_like(ser)
                    shifted[okr] = term[okr - offs[r]]
                    term = (coef[r] / n) * shifted
                    if not np.any(term):
                        break
                    nxt = nxt + term
                ser = nxt

        for j in np.nonzero(last == i)[0]:
            ser = _reduce(ser, dims, R + int(j), sense[int(j)], int(b[j]))
            dims[R + int(j)] = 1
            sub, stride, lev = _lattice(dims, R)

    if ser.size != int(np.prod(np.asarray(N) + 1)):
        raise RuntimeError(
            "pfqn_manjunath: a constraint row was never discharged; the "
            "elimination order is inconsistent")
    # The closed network's own equalities: degree exactly N_r in every class.
    flat = 0
    st = 1
    for r in range(R):
        flat += int(N[r]) * st
        st *= int(N[r]) + 1
    return float(ser[flat]), peak


def pfqn_manjunath(L, N, Z=None, A=None, b=None, sense: Optional[str] = None,
                   stats: bool = False):
    """
    Exact normalizing constant of a constrained closed product-form network.

    Args:
        L: Service demand of class r at queueing station i, (M, R).
        N: Population of class r, (R,) nonnegative integers.
        Z: Think time of class r at delay station k, (Mz, R) or (R,); default none.
        A: Extra constraint coefficients on n.ravel(order='F'),
            (J, (M+Mz)*R) nonnegative integers; default none.
        b: Extra constraint right-hand sides, (J,) integers; default none.
        sense: One character per row, 'E' (=), 'L' (<=) or 'G' (>); default all 'L'.
        stats: Also return the per-class decomposition. Requires the ONE
            configuration in which the truncated product form is the EXACT
            stationary law: a single queueing station inside the region and a
            SINGLE DELAY STATION OUTSIDE IT. Anything else is refused by name.

    Returns:
        Tuple (G, lG, peak), or (G, lG, peak, PfqnManjunathStats) when `stats`.

    WHY THE CONFIGURATION IS NOT A CONVENIENCE. With one queueing station the
    state is the queue occupancy alone (the delay holds the complement) and
    every transition moves one job of one class by one unit, so the chain is a
    multidimensional birth-death process. That process is reversible, and
    Kelly's truncation theorem then applies verbatim: restricting it to the
    coordinate-convex set A n <= b and renormalizing gives exactly the truncated
    product form. Add a second queueing station and the delay -> q1 -> q2 ->
    delay cycle destroys reversibility; truncation no longer preserves the
    product form, measured at 131% relative error on the stationary law of a
    2-class, N = [2 2] instance. G and lG stay correct as a sum over the
    admissible set in every configuration; only the metrics are withheld.

    Everything follows from two ratios of normalizing constants, both taken in
    the log domain so the internal rescaling cancels without being
    reconstructed:

        X_r         = G(N - e_r ; b - A[:, qcol_r]) / G(N ; b)
        P(n_qr = k) = G(N ; b, with the added row n_qr = k) / G(N ; b)

    The first is the loss network's g(C - A e_r) in another guise: removing one
    class r job from the queue leaves a state whose admission rule is shifted by
    that job's own requirement column. The second is what an '=' row is for.
    """
    N = np.atleast_1d(np.asarray(N, dtype=np.float64)).ravel()
    R = N.size
    L = np.zeros((0, R)) if L is None else np.atleast_2d(np.asarray(L, dtype=np.float64))
    if L.size == 0:
        L = np.zeros((0, R))
    if L.shape[1] != R:
        raise ValueError("pfqn_manjunath: L must have %d columns, one per class" % R)
    if Z is None:
        Z = np.zeros((0, R))
    else:
        Z = np.atleast_2d(np.asarray(Z, dtype=np.float64))
        if Z.size == 0:
            Z = np.zeros((0, R))
    if Z.shape[1] != R:
        raise ValueError("pfqn_manjunath: Z must have %d columns, one per class" % R)
    M, Mz = L.shape[0], Z.shape[0]
    S = M + Mz

    b = np.zeros(0) if b is None else np.atleast_1d(np.asarray(b, dtype=np.float64)).ravel()
    J = b.size
    if A is None or J == 0:
        A = np.zeros((J, S * R))
    else:
        A = np.atleast_2d(np.asarray(A, dtype=np.float64))
    if A.shape != (J, S * R):
        raise ValueError("pfqn_manjunath: A must be %dx%d (J x (M+Mz)*R), acting on "
                         "n.ravel(order='F'), got %s" % (J, S * R, A.shape))
    if sense is None:
        sense = 'L' * J
    sense = str(sense).upper()
    if len(sense) != J:
        raise ValueError("pfqn_manjunath: sense must have one character per row of b")
    if any(c not in 'ELG' for c in sense):
        raise ValueError("pfqn_manjunath: sense must contain only 'E' (=), 'L' (<=) "
                         "or 'G' (>)")

    if np.any(N < 0):
        return _ret(0.0, -np.inf, 0, stats, R)
    if np.any(np.abs(N - np.round(N)) > 1e-9):
        raise ValueError("pfqn_manjunath: N must contain nonnegative integers")
    N = np.round(N).astype(np.int64)
    if np.any(A < 0) or np.any(np.abs(A - np.round(A)) > 1e-9):
        raise ValueError("pfqn_manjunath: A must contain nonnegative integers; the "
                         "residue argument counts whole units")
    if np.any(np.abs(b - np.round(b)) > 1e-9):
        raise ValueError("pfqn_manjunath: b must contain integers; the residue "
                         "argument counts whole units")
    A = np.round(A).astype(np.int64)
    b = np.round(b).astype(np.int64)

    # A row of zeros constrains nothing, so it is decided here rather than
    # carried as a one-coefficient dimension: 0 = b, 0 <= b and 0 > b are each
    # settled by the sign of b alone. Same for a negative right-hand side, which
    # no nonnegative combination can meet ('E', 'L') or fail to beat ('G').
    keep = np.ones(J, dtype=bool)
    for j in range(J):
        trivial = not np.any(A[j, :] != 0)
        s = sense[j]
        if s == 'E':
            if b[j] < 0 or (trivial and b[j] != 0):
                return _ret(0.0, -np.inf, 0, stats, R)
            if trivial:
                keep[j] = False
        elif s == 'L':
            if b[j] < 0:
                return _ret(0.0, -np.inf, 0, stats, R)
            if trivial:
                keep[j] = False
        else:                       # 'G'
            if b[j] < 0:
                keep[j] = False     # 0 > negative always holds
            elif trivial:
                return _ret(0.0, -np.inf, 0, stats, R)  # 0 > nonneg never holds
    A = A[keep, :]
    b = b[keep]
    sense = ''.join(c for c, k in zip(sense, keep) if k)
    J = b.size

    if S == 0:
        # No station: the only state is empty, admissible when every class is too.
        if np.all(N == 0):
            return _ret(1.0, 0.0, 1, stats, R)
        return _ret(0.0, -np.inf, 0, stats, R)

    # Every monomial that survives the extraction has total degree sum(N) in the
    # demands, so a common rescaling of L and Z moves lG by a known amount and
    # nothing else. The exponent is chosen as in pfqn_ca, from the largest single
    # term the network can produce, so that the series is centred near 1.
    Nt = int(N.sum())
    cscale = 1.0
    if Nt > 0:
        lGest = -np.inf
        for i in range(M):
            t, ok = 0.0, True
            for r in range(R):
                if N[r] > 0:
                    if L[i, r] > 0:
                        t += N[r] * math.log(L[i, r])
                    else:
                        ok = False
                        break
            if ok:
                lGest = max(lGest, t)
        if Mz > 0:
            t, ok = 0.0, True
            for r in range(R):
                if N[r] > 0:
                    zs = float(Z[:, r].sum())
                    if zs > 0:
                        t += N[r] * math.log(zs) - _factln(int(N[r]))
                    else:
                        ok = False
                        break
            if ok:
                lGest = max(lGest, t)
        if np.isfinite(lGest):
            cscale = float(math.ldexp(1.0, _matlab_round(lGest / (Nt * math.log(2)))))
    Ls = L / cscale
    Zs = Z / cscale

    # A '>' row is the complement of a '<=' row at the same right-hand side,
    # which is how the paper discharges it (Eqn 6). With several such rows the
    # product of the complements expands by inclusion-exclusion, so the series is
    # evaluated once per subset of them, with the subset's rows re-entered as
    # '<=' and the rest dropped. Exact, and the only place the cost is
    # exponential -- in the number of '>' rows, which is normally zero.
    gt = [j for j in range(J) if sense[j] == 'G']
    base = [j for j in range(J) if sense[j] != 'G']
    K = len(gt)
    Gs = 0.0
    peak = 0
    for mask in range(1 << K):
        sel = sorted(base + [gt[t] for t in range(K) if mask & (1 << t)])
        subsense = ''.join('L' if sense[j] == 'G' else sense[j] for j in sel)
        g, pk = _series(Ls, Zs, N, A[sel, :] if sel else np.zeros((0, S * R), dtype=np.int64),
                        b[sel] if sel else np.zeros(0, dtype=np.int64), subsense)
        Gs += (1 - 2 * (bin(mask).count('1') % 2)) * g
        peak = max(peak, pk)

    if Gs <= 0.0:
        # Either the admissible set is empty or the '>' complements cancelled it.
        return _ret(0.0, -np.inf, peak, stats, R)
    lG = math.log(Gs) + Nt * math.log(cscale)
    G = math.exp(lG)
    if not stats:
        return G, lG, peak
    return G, lG, peak, _stats(L, N, Z, A, b, sense, lG, M, Mz, S, R)
