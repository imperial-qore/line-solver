"""
Product form of a stochastic Petri net: decide whether one exists and derive
the per-level factors g_l that mdd_rec and spn_metrics take as input.

THIS IS THE PART THE MDD-REC PAPER DECLARES OUT OF SCOPE (FGCS Sec. 3.2).
Every other function in api/spn/ receives the g_l already formed; this one
derives them from the net, which is what lets a solver reach them.

The theory, in one paragraph. Write I(t), O(t) for the input and output vectors
of mode t and lambda_t for its rate constant. Henderson-Taylor and
Coleman-Henderson-Taylor show that a net whose firing rate has the form

    r_t(m) = lambda_t psi(m - I(t)) / psi(m),     m >= I(t)

has invariant measure pi(m) = psi(m) prod_l y_l^{m_l} whenever the positive
vector y satisfies COMPLEX BALANCE: reading the distinct vectors appearing as
some I(t) or O(t) as the COMPLEXES of the net, the flow into every complex must
equal the flow out of it,

    sum_{t : O(t)=v} lambda_t y^{I(t)} = ( sum_{t : I(t)=v} lambda_t ) y^v.

Two choices of psi are realisable in LINE's own rate law, and they are the two
this module tests for:

    psi = 1             r_t = lambda_t, the rate of a SINGLE-SERVER mode.
                        pi(m) = prod_l y_l^{m_l}, so g_l(k) = y_l^k.
    psi = prod 1/m_l!   r_t = lambda_t prod_l m_l!/(m_l-I_l)!, MASS ACTION,
                        reached through Transition.set_firing_rate_dependence
                        or, for a mode drawing one token from one place, by
                        infinite-server semantics.
                        pi(m) = prod_l y_l^{m_l}/m_l!, so g_l(k) = y_l^k/k!.

Which one holds is not guessed from the model API: the effective rate LINE
would use, lambda_t min(enabling degree, servers) g(m), is EVALUATED at every
reachable marking and compared against both laws. A net that matches neither
under one common psi is refused by name, never approximated.

Solving for y. Complex balance reads A_lambda Psi(y) = 0 with A_lambda the
Laplacian of the weighted digraph on complexes and Psi(y)_v = y^v. That Laplacian
is the TRANSPOSED GENERATOR of a Markov chain that hops from complex to complex
at the rate of the mode joining them, so its kernel on one linkage class is that
chain's stationary distribution and ctmc_solve returns it -- strictly positive
exactly when the class is strongly connected, which is weak reversibility. With
that positive vector kappa in hand y follows from the LINEAR system in x = log y

    (v - v0) x = log kappa_v - log kappa_v0,   v, v0 in the same linkage class,

solved in minimum norm. Feinberg's Deficiency Zero Theorem says this system is
consistent for every choice of rate constants when the net is weakly reversible
and its deficiency c - l - s is zero, which is why those two numbers are
reported; but consistency is CHECKED rather than assumed, so a net of positive
deficiency whose particular rates still admit a complex-balanced point is
accepted on the evidence.

The gauge, and why the minimum-norm solution is the canonical one. Complex
balance fixes y only up to y -> y * exp(u) for any u orthogonal to the
stoichiometric subspace S. Such a shift multiplies pi(m) by exp(u.m), which is
CONSTANT on one compatibility class, so every reported measure is invariant under
it -- but the normalising constant G itself is not, it scales by that constant. A
gauge must therefore be FIXED, or the four codebases would report four different G
on the same net. The one fixed here is x in the row space of the constraint
matrix, i.e. the minimum-norm solution, reached in a form that is unique whichever
least-squares primitive a codebase carries: solve (rows rows^T) w = rhs and set
x = rows^T w. Any two solutions w of that system give the SAME rows^T w, so the
answer does not depend on how the rank-deficient solve breaks its tie.

References
----------
J. L. Coleman, W. Henderson, P. G. Taylor, "Product form equilibrium
distributions and a convolution algorithm for stochastic Petri nets",
Performance Evaluation 26(3), 1996.
M. Feinberg, "Complex balancing in general kinetic systems", Arch. Rational
Mech. Anal. 49, 1972.
D. F. Anderson, G. Craciun, T. G. Kurtz, "Product-form stationary distributions
for deficiency zero chemical reaction networks", Bull. Math. Biol. 72, 2010.

See also: mdd_rec, spn_metrics, spn_mdd, spn_conv.
"""

from math import factorial, log, exp, isfinite
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..io.logging import line_error, line_printf
from ...constants import SchedStrategy
from ..mc.ctmc import ctmc_solve
from .mdd import spn_mdd

__all__ = ['spn_pf']


def spn_pf(model, options: Optional[Dict] = None) -> Dict:
    """Derive the product form of a stochastic Petri net.

    Parameters
    ----------
    model : a Network holding Places and Transitions
    options : dict with optional keys
        bound   - per-place-level token bound, passed to spn_mdd
        tol     - relative tolerance of the rate-law and complex-balance checks
                  (default 1e-9)
        verbose - print the certificate (default False)

    Returns
    -------
    dict with keys g, y, kind, complexes, deficiency, linkage, srank,
    weaklyreversible, residual, mdds and info.
    """
    if options is None:
        options = {}
    tol = float(options.get('tol', 1e-9))
    verbose = bool(options.get('verbose', False))

    mddopt = {'descriptor': False}
    if options.get('bound') is not None:
        mddopt['bound'] = options['bound']
    mdds, _, info = spn_mdd(model, mddopt)

    L = int(info['nplacelevels'])
    md = info['modes']
    E = len(md)
    if E == 0:
        line_error('spn_pf', 'the net has no timed mode')

    # ---- a queueing place holds an embedded server, not a token container
    sn = model.getStruct()
    node_to_station = np.ravel(np.asarray(sn.nodeToStation)).astype(int)
    for pp, node in enumerate(info['places']):
        ist = int(node_to_station[node])
        if ist < 0:
            continue
        sched = sn.sched[ist] if isinstance(sn.sched, dict) else np.ravel(sn.sched)[ist]
        if SchedStrategy(sched) != SchedStrategy.INF:
            line_error('spn_pf',
                       'place %s is a QUEUEING place (scheduling %s): its embedded service is '
                       'state that the marking does not carry, so the net is not the '
                       'token-container Petri net this product form is written for'
                       % (sn.nodenames[node], SchedStrategy(sched).name))

    # ---- the rate constants and the structural vectors
    lam = np.zeros(E)
    Iv = np.zeros((E, L))
    Ov = np.zeros((E, L))
    for e in range(E):
        lam[e] = float(np.asarray(md[e]['D1']).ravel()[0])
        Iv[e, :] = md[e]['enab']
        Ov[e, :] = md[e]['fire']
        if np.any(np.isfinite(md[e]['inhib'])):
            line_error('spn_pf',
                       'mode %d of node %d has an inhibitor arc. An inhibitor zeroes the firing '
                       'rate on markings that still satisfy m >= I(t), so the rate is not '
                       'lambda*psi(m-I)/psi(m) on any psi and the net has no product form of '
                       'this kind' % (md[e]['mode'] + 1, md[e]['trans'] + 1))
        if md[e]['srv'] != 1 and not np.any(md[e]['enab'] > 0):
            line_error('spn_pf',
                       'mode %d of node %d has %g servers but consumes from no place, so its '
                       'enabling degree is unbounded and its firing rate undefined'
                       % (md[e]['mode'] + 1, md[e]['trans'] + 1, md[e]['srv']))
        if not lam[e] > 0:
            line_error('spn_pf', 'mode %d of node %d has a non-positive firing rate'
                       % (md[e]['mode'] + 1, md[e]['trans'] + 1))

    # ---- which psi does LINE's own rate law follow on this net?
    states = info['mdd'].enumerate()
    kind = _ratelaw(states, md, lam, info, tol)

    # ---- complexes and the weighted digraph on them
    C, src, dst = _complexes(Iv, Ov)
    c = C.shape[0]

    # ---- Laplacian of the complex graph: A[j, i] is the rate of the arc i -> j
    A = np.zeros((c, c))
    for e in range(E):
        if src[e] == dst[e]:
            continue                                   # a mode that moves nothing
        A[dst[e], src[e]] += lam[e]
        A[src[e], src[e]] -= lam[e]

    # ---- linkage classes, weak reversibility, deficiency
    lclass, nlink = _linkage(c, src, dst)
    wr = _weakly_reversible(c, src, dst, lclass, nlink)
    srank = int(np.linalg.matrix_rank(Ov - Iv, tol=1e-9))
    deficiency = c - nlink - srank

    # ---- kappa: the positive kernel of the Laplacian on each linkage class
    kappa = np.zeros(c)
    for b in range(nlink):
        idx = np.nonzero(lclass == b)[0]
        if idx.size == 1:
            kappa[idx[0]] = 1.0
            continue
        # A[idx, idx] is the transposed generator of the complex-hopping chain,
        # so its kernel is that chain's stationary law.
        v = np.ravel(ctmc_solve(A[np.ix_(idx, idx)].T))
        if np.any(v <= 0):
            line_error('spn_pf',
                       'linkage class %d of the complex graph carries no flow through complex '
                       '%d, so the net admits no positive complex-balanced point. A weakly '
                       'reversible net has a strictly positive balance flow on every linkage '
                       'class; this one is %s'
                       % (b + 1, idx[int(np.argmin(v))] + 1, _wrtext(wr)))
        kappa[idx] = v / v.max()

    # ---- x = log y from the linear system on each linkage class
    rows: List[np.ndarray] = []
    rhs: List[float] = []
    for b in range(nlink):
        idx = np.nonzero(lclass == b)[0]
        v0 = idx[0]
        for v in idx[1:]:
            rows.append(C[v, :] - C[v0, :])
            rhs.append(log(kappa[v]) - log(kappa[v0]))
    if not rows:
        x = np.zeros(L)
    else:
        Ar = np.asarray(rows, dtype=float)
        br = np.asarray(rhs, dtype=float)
        # minimum norm through the row space; see the gauge note in the docstring
        x = Ar.T @ (np.linalg.pinv(Ar @ Ar.T) @ br)
        res = float(np.max(np.abs(Ar @ x - br))) if br.size else 0.0
        scale = max(1.0, float(np.max(np.abs(br))) if br.size else 1.0)
        if res > tol * scale:
            line_error('spn_pf',
                       'the complex-balance equations are inconsistent (residual %.3e): this net '
                       'has no product form of the tested kind at these rates. Its deficiency is '
                       '%d and it is %s; the Deficiency Zero Theorem guarantees a solution only '
                       'at deficiency 0 with weak reversibility'
                       % (res, deficiency, _wrtext(wr)))
    y = np.exp(x)

    # ---- verify complex balance itself, which is what makes pi stationary
    Psi = np.array([np.prod(y ** C[v, :]) for v in range(c)])
    resb = float(np.max(np.abs(A @ Psi))) if c else 0.0
    scaleb = max(float(np.max(np.abs(A) @ Psi)) if c else 0.0, np.finfo(float).tiny)
    if resb > tol * scaleb:
        line_error('spn_pf',
                   'complex balance fails at the computed point (relative residual %.3e), so the '
                   'product form would not be stationary' % (resb / scaleb))

    # ---- the per-level factors, tabulated over the reachable domain
    domain = np.asarray(mdds.domain if hasattr(mdds, 'domain') else mdds['domain'], dtype=int)
    g = []
    for l in range(L):
        k = np.arange(int(domain[l]))
        if kind == 'geometric':
            g.append(y[l] ** k)
        else:
            g.append((y[l] ** k) / np.array([float(factorial(int(j))) for j in k]))

    pf = {
        'g': g,
        'y': y,
        'kind': kind,
        'complexes': C,
        'deficiency': int(deficiency),
        'linkage': int(nlink),
        'srank': srank,
        'weaklyreversible': bool(wr),
        'residual': resb / scaleb,
        'mdds': mdds,
        'info': info,
    }
    if verbose:
        line_printf('\nSPN product form: %s, %d complexes, %d linkage classes, rank %d, '
                    'deficiency %d, %s\n' % (kind, c, nlink, srank, deficiency, _wrtext(wr)))
        line_printf('  y = %s (complex-balance residual %.2e)\n'
                    % (np.array2string(y, precision=6), resb / scaleb))
    return pf


def _ratelaw(states, md, lam, info, tol) -> str:
    """Which psi reproduces the rate LINE would actually use, everywhere.

    Both candidates are tried on every mode at every reachable marking; a net
    matching neither, or matching different ones on different modes, has no
    product form of this family and is refused by name.
    """
    E = len(md)
    L = int(info['nplacelevels'])
    okgeo = True
    okma = True
    for s in range(states.shape[0]):
        m = np.asarray(states[s, :L], dtype=float)
        marc = None
        for e in range(E):
            if np.any(m < md[e]['enab']):
                continue
            actual = lam[e] * _servers(m, md[e])
            if md[e]['dep'] is not None:
                if marc is None:
                    marc = _arcmatrix(m, info)
                actual = actual * float(md[e]['dep'](marc))
            okgeo = okgeo and _close(actual, lam[e], tol)
            okma = okma and _close(actual, lam[e] * _massaction(m, md[e]['enab']), tol)
            if not okgeo and not okma:
                line_error('spn_pf',
                           "mode %d of node %d fires at rate %g in a reachable marking, which is "
                           "neither its rate constant (single-server, psi = 1) nor its mass-action "
                           "rate %g (psi = prod 1/m!). LINE's rate law on this mode is "
                           "lambda*min(enabling degree, servers)*g(m), and no psi puts that in "
                           "the form lambda*psi(m-I)/psi(m)"
                           % (md[e]['mode'] + 1, md[e]['trans'] + 1, actual,
                              lam[e] * _massaction(m, md[e]['enab'])))
    return 'geometric' if okgeo else 'massaction'


def _servers(m, mde) -> float:
    """min(enabling degree, servers): the sets of tokens firing at once."""
    deg = np.inf
    for l in np.nonzero(np.asarray(mde['enab']) > 0)[0]:
        deg = min(deg, np.floor(m[l] / mde['enab'][l]))
    if not isfinite(deg):
        deg = 1.0                                   # consumes nothing: always one set
    return float(min(deg, mde['srv']))


def _massaction(m, enab) -> float:
    """prod_l m_l!/(m_l - I_l)!, the ordered ways to pick the input tokens."""
    r = 1.0
    enab = np.asarray(enab)
    for l in np.nonzero(enab > 0)[0]:
        for j in range(int(enab[l])):
            r *= (m[l] - j)
    return r


def _close(a, b, tol) -> bool:
    return abs(a - b) <= tol * max(1.0, max(abs(a), abs(b)))


def _arcmatrix(m, info) -> np.ndarray:
    """Place-major level vector -> the (nnodes x nclasses) marking matrix that a
    Transition.set_firing_rate_dependence handle is written against."""
    R = int(info['nclasses'])
    M = np.zeros((int(info['nnodes']), R))
    for pp, node in enumerate(info['places']):
        for k in range(R):
            M[node, k] = m[pp * R + k]
    return M


def _complexes(Iv, Ov) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The distinct input and output vectors, in first-seen order, and the arc
    each mode draws. First-seen rather than sorted so that the complex indices
    agree with the MATLAB, Java and C++ twins."""
    E = Iv.shape[0]
    stack = np.vstack([Iv, Ov])
    seen: Dict[Tuple, int] = {}
    order: List[np.ndarray] = []
    idx = np.zeros(stack.shape[0], dtype=int)
    for i in range(stack.shape[0]):
        key = tuple(stack[i, :])
        if key not in seen:
            seen[key] = len(order)
            order.append(stack[i, :])
        idx[i] = seen[key]
    C = np.asarray(order, dtype=float) if order else np.zeros((0, Iv.shape[1]))
    return C, idx[:E], idx[E:]


def _linkage(c, src, dst) -> Tuple[np.ndarray, int]:
    """Connected components of the UNDIRECTED complex graph."""
    lclass = -np.ones(c, dtype=int)
    nlink = 0
    for v in range(c):
        if lclass[v] >= 0:
            continue
        stack = [v]
        lclass[v] = nlink
        while stack:
            u = stack.pop()
            nb = list(dst[src == u]) + list(src[dst == u])
            for w in nb:
                if lclass[w] < 0:
                    lclass[w] = nlink
                    stack.append(w)
        nlink += 1
    return lclass, nlink


def _weakly_reversible(c, src, dst, lclass, nlink) -> bool:
    """Every linkage class strongly connected in the DIRECTED complex graph."""
    for b in range(nlink):
        idx = np.nonzero(lclass == b)[0]
        fwd = _reach(idx[0], src, dst, c)
        bwd = _reach(idx[0], dst, src, c)
        if not fwd[idx].all() or not bwd[idx].all():
            return False
    return True


def _reach(v0, frm, to, c) -> np.ndarray:
    seen = np.zeros(c, dtype=bool)
    seen[v0] = True
    stack = [v0]
    while stack:
        u = stack.pop()
        for w in to[frm == u]:
            if not seen[w]:
                seen[w] = True
                stack.append(w)
    return seen


def _wrtext(wr) -> str:
    return 'weakly reversible' if wr else 'not weakly reversible'
