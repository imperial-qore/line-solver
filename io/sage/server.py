#!/usr/bin/env sage-python
"""
LINE symbolic REST service (SageMath backend).

Exposes the computer algebra operations that LINE needs to all three
codebases over a small JSON/HTTP protocol, so that the JAR (which has no CAS
at all), MATLAB (Symbolic Math Toolbox licence) and Python (sympy) can share
one exact engine and one canonical normal form.

Endpoints. The four meta routes are GET, as in the other imperialqore
line-*-rest services; the solve routes are POST:

  /api/v1/health             -                                -> {status, timestamp}
  /api/v1/ready              -                                -> {ready, checks}
  /api/v1/info               -                                -> {name, version, ...}
  /api/v1/citation           -                                -> {notice, references}
  /api/v1/ctmc/solve         {Q, symbols, ...}                -> {pi, num, den, ...}
  /api/v1/ctmc/sensitivity   {Q, symbols, theta, reward}      -> {dpi, S, SS, pi}
  /api/v1/simplify           {exprs, form}                    -> {results}
  /api/v1/diff               {exprs, var, order}              -> {results}
  /api/v1/eval               {exprs, values}                  -> {exact, values}
  /api/v1/fluid/odes         {rhs, vars, want}                -> {jacobian, equilibria, latex}

Every response is a JSON object with a "status" field, "ok" or "error". On
error it carries "message" and a machine-readable "code".

Expressions travel as plain ASCII infix strings, e.g. "2*x1 - 3*x2". They are
parsed by an explicit AST walker (see ExprParser), never by eval/sage_eval:
that keeps the service safe to expose on a port and, more importantly, keeps
the arithmetic exact. Decimal literals are read off the source text and
converted to exact rationals, so "0.1" becomes 1/10 rather than the binary
double nearest to it.

Concurrency: the server forks per request. That gives each computation its own
process, so signal.alarm can enforce the per-request timeout and a runaway
symbolic solve can be killed without taking the service down. It also means a
long solve never blocks the health probe used by container start-up.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import ast
import json
import os
import signal
import socketserver
import sys
import traceback
from datetime import datetime
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from http.server import BaseHTTPRequestHandler, HTTPServer

from sage.all import (QQ, SR, ZZ, PolynomialRing, copy, diff, exp, gcd, latex,
                      log, matrix, max_symbolic, min_symbolic, solve, sqrt,
                      var, version)

API_PREFIX = "/api/v1"
API_VERSION = "v1"
SERVICE_VERSION = "1.0.0"
DEFAULT_TIMEOUT_S = 300
MAX_BODY_BYTES = 64 * 1024 * 1024

# This service only exposes SageMath. All scientific credit for the algebra it
# performs belongs to the Sage developers; results obtained through it must
# cite SageMath, not this service or the image that carries it.
NOTICE = ("This service is a REST packaging of SageMath and adds no computer "
          "algebra of its own. Publications reporting results obtained through "
          "it must cite SageMath, and must not cite the imperialqore repository "
          "or image.")

CITATION = {
    "tool": "SageMath",
    "full_name": "SageMath, the Sage Mathematics Software System",
    "url": "https://www.sagemath.org",
    "source": "https://github.com/sagemath/sage",
    "authors": ["The Sage Developers"],
    "license": "GPL-2.0-or-later",
    "references": [
        {"text": "The Sage Developers, SageMath, the Sage Mathematics Software "
                 "System (Version %s), https://www.sagemath.org",
         "doi": "10.5281/zenodo.593563",
         "note": "Substitute the version reported by /api/v1/info."},
    ],
}


class LineSymError(Exception):
    """Error carrying a machine-readable code for the JSON response."""

    def __init__(self, code, message):
        Exception.__init__(self, message)
        self.code = code
        self.message = message


# ---------------------------------------------------------------------------
# Expression parsing
# ---------------------------------------------------------------------------

_SR_FUNCTIONS = {
    "sqrt": sqrt,
    "exp": exp,
    "log": log,
    "abs": abs,
    "min": min_symbolic,
    "max": max_symbolic,
}


class ExprParser:
    """Parses an infix expression string into a Sage ring element.

    The walker accepts +, -, *, /, ** and parentheses, plus, when calls are
    allowed, the whitelisted functions in _SR_FUNCTIONS. Anything else is a
    parse error rather than an evaluation: no name lookup escapes the
    environment passed in.

    Numeric literals are converted through their *source text*, not through
    the float that ast produces, so a decimal literal becomes the exact
    rational with the same decimal expansion.
    """

    def __init__(self, resolve_name, one, allow_calls=False):
        self.resolve_name = resolve_name
        self.one = one
        self.allow_calls = allow_calls

    def parse(self, text):
        if not isinstance(text, str):
            # Numbers may arrive unquoted in the JSON payload.
            text = repr(text)
        text = text.strip()
        if not text:
            raise LineSymError("parse", "empty expression")
        # Sage prints powers as x^2 and MATLAB writes them the same way, but
        # Python reads ^ as bitwise xor, which also binds looser than *: read
        # as Python, "2*x^2" would silently become (2*x)^2. The substitution
        # is textual and unconditional because xor has no meaning in this
        # protocol.
        text = text.replace("^", "**")
        try:
            tree = ast.parse(text, mode="eval")
        except SyntaxError as e:
            raise LineSymError("parse", "cannot parse expression %r: %s" % (text, e))
        return self._eval(tree.body, text)

    def _number(self, node, source):
        segment = ast.get_source_segment(source, node)
        if segment is None:
            segment = repr(node.value)
        segment = segment.strip()
        try:
            return self.one * QQ(Fraction(Decimal(segment)))
        except (InvalidOperation, ValueError, ArithmeticError):
            raise LineSymError("parse", "cannot read numeric literal %r" % segment)

    def _eval(self, node, source):
        if isinstance(node, ast.Constant):
            if isinstance(node.value, bool) or node.value is None:
                raise LineSymError("parse", "unsupported literal %r" % (node.value,))
            if isinstance(node.value, (int, float)):
                return self._number(node, source)
            if isinstance(node.value, str):
                # A quoted rational such as "3/2" nested in the expression.
                return self.parse(node.value)
            raise LineSymError("parse", "unsupported literal %r" % (node.value,))
        if isinstance(node, ast.Name):
            return self.resolve_name(node.id)
        if isinstance(node, ast.UnaryOp):
            operand = self._eval(node.operand, source)
            if isinstance(node.op, ast.UAdd):
                return operand
            if isinstance(node.op, ast.USub):
                return -operand
            raise LineSymError("parse", "unsupported unary operator")
        if isinstance(node, ast.BinOp):
            left = self._eval(node.left, source)
            right = self._eval(node.right, source)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                if right == 0:
                    raise LineSymError("parse", "division by zero in expression")
                return left / right
            if isinstance(node.op, ast.Pow):
                return left ** self._exponent(node.right, source)
            raise LineSymError("parse", "unsupported binary operator")
        if isinstance(node, ast.Call) and self.allow_calls:
            if not isinstance(node.func, ast.Name):
                raise LineSymError("parse", "unsupported call target")
            fn = _SR_FUNCTIONS.get(node.func.id)
            if fn is None or node.keywords:
                raise LineSymError("parse", "unsupported function %r" % node.func.id)
            args = [self._eval(a, source) for a in node.args]
            return fn(*args)
        raise LineSymError("parse", "unsupported expression element %s"
                           % type(node).__name__)

    def _exponent(self, node, source):
        """Exponent of a power.

        In the symbolic ring anything constant is allowed, including the 1/p of
        a p-norm smoothed fluid drift. In a rational function field only an
        integer is: x^(1/8) is not a rational function. int() must not be used
        to decide that, because it TRUNCATES: int(1/8) is 0, which turned
        (1 + t^8)^(1/8) into 1 and returned a drift that was wrong by a
        thousandth with nothing reported.
        """
        value = self._eval(node, source)
        if self.allow_calls:
            return value
        try:
            if QQ(value).denominator() == 1:
                return int(QQ(value))
        except (TypeError, ValueError, ArithmeticError):
            pass
        raise LineSymError("parse",
                           "non-integer exponent in expression: a rational function field "
                           "has no such element")


def _ctmc_field(symbols):
    """Rational function field QQ(x1, ..., xE) the generator lives in.

    With no symbols the generator is numeric and QQ itself is the right field;
    PolynomialRing rejects an empty generator list, so that case is split out.
    """
    if not symbols:
        return QQ, {}
    R = PolynomialRing(QQ, names=[str(s) for s in symbols])
    F = R.fraction_field()
    env = dict(zip([str(s) for s in symbols], F.gens()))
    return F, env


def _ctmc_parser(symbols):
    F, env = _ctmc_field(symbols)

    def resolve(name):
        if name not in env:
            raise LineSymError("parse", "unknown symbol %r; declare it in \"symbols\"" % name)
        return env[name]

    return F, ExprParser(resolve, F.one(), allow_calls=False)


def _sr_parser():
    def resolve(name):
        return var(name)

    return ExprParser(resolve, SR.one(), allow_calls=True)


def _parse_matrix(rows, parser, F):
    if not isinstance(rows, list) or not rows:
        raise LineSymError("input", "\"Q\" must be a non-empty list of rows")
    n = len(rows)
    entries = []
    for i, row in enumerate(rows):
        if not isinstance(row, list) or len(row) != n:
            raise LineSymError("input", "row %d of \"Q\" is not of length %d" % (i, n))
        entries.append([parser.parse(e) for e in row])
    return matrix(F, n, n, entries)


def _to_str(e):
    return str(e)


def _coeff_lcm_denominator(p):
    """Least common multiple of the denominators of a polynomial's coefficients."""
    try:
        return ZZ(p.denominator())
    except (AttributeError, TypeError):
        return ZZ(1)


def _integer_coefficients(p):
    try:
        return [ZZ(c) for c in p.coefficients()]
    except (AttributeError, TypeError):
        return [ZZ(p)]


def _leading_coefficient(p):
    try:
        return p.lc()
    except AttributeError:
        return p


def _primitive(num, den):
    """Puts num/den in the primitive normal form: integer coefficients, no
    common integer factor, positive leading coefficient in the denominator.

    The fraction field only cancels common polynomial factors, so the same
    probability comes back scaled differently from one state to the next
    (6*x2^3/(36*x1^3 + ...) next to 6*x1^2*x2/(6*x1^3 + ...)): a rational
    constant is a unit, so it is never part of a gcd. Two expressions that
    print differently but are equal defeat the whole point of returning a
    normal form, hence this step.
    """
    d = _coeff_lcm_denominator(num).lcm(_coeff_lcm_denominator(den))
    if d != 1:
        num, den = num * d, den * d
    coeffs = _integer_coefficients(num) + _integer_coefficients(den)
    g = gcd(coeffs)
    if g and g != 1:
        num, den = num / g, den / g
    if _leading_coefficient(den) < 0:
        num, den = -num, -den
    return num, den


def _canonical(e, F):
    """Reduced rational function in the primitive normal form of _primitive."""
    if F is QQ:
        # QQ already prints its own normal form, sign and all.
        return F(e)
    e = F(e)
    num, den = _primitive(e.numerator(), e.denominator())
    return F(num) / F(den)


# ---------------------------------------------------------------------------
# CTMC steady state
# ---------------------------------------------------------------------------

def _makeinfgen(Q):
    """Zeroes the diagonal and sets it to minus the off-diagonal row sum.

    Mirrors ctmc_makeinfgen in all three codebases, so a caller may send
    either the raw rate matrix or an already-formed generator.
    """
    n = Q.nrows()
    A = copy(Q)
    for i in range(n):
        A[i, i] = 0
        A[i, i] = -sum(A[i, j] for j in range(n))
    return A


def _weak_components(Q):
    """Weakly connected components of the generator's support.

    ctmc_solve.m tests connectivity after substituting 1.0 for every symbol,
    which can cancel two rates that are structurally distinct. Testing the
    entry against the zero element of the rational function field instead is
    exact and never manufactures a spurious disconnection.
    """
    n = Q.nrows()
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    for i in range(n):
        for j in range(i + 1, n):
            if Q[i, j] != 0 or Q[j, i] != 0:
                union(i, j)

    labels = {}
    comp = [0] * n
    for i in range(n):
        r = find(i)
        if r not in labels:
            labels[r] = len(labels) + 1
        comp[i] = labels[r]
    return len(labels), comp


def _nonzero_support(Q):
    """Iteratively drops states with an all-zero row or an all-zero column.

    This is the elimination loop of ctmc_solve.m: an absorbing state leaves a
    zero row, removing it can zero the diagonal of the states that only fed
    it, and the elimination therefore cascades to a fixed point.
    """
    n = Q.nrows()
    keep = list(range(n))
    A = Q
    while True:
        rows = [i for i in range(A.nrows())
                if any(A[i, j] != 0 for j in range(A.ncols()))]
        cols = [j for j in range(A.ncols())
                if any(A[i, j] != 0 for i in range(A.nrows()))]
        idx = [i for i in range(A.nrows()) if i in set(rows) and i in set(cols)]
        if len(idx) == A.nrows():
            return A, keep
        if not idx:
            return A.matrix_from_rows_and_columns([], []), []
        keep = [keep[i] for i in idx]
        A = _makeinfgen(A.matrix_from_rows_and_columns(idx, idx))


def _solve_block(Q, F):
    """Solves pi * Q = 0 with sum(pi) = 1 on an irreducible block.

    The normalization replaces the last column of Q with ones and the right
    hand side with e_n, exactly as ctmc_solve.m does.
    """
    n = Q.nrows()
    if n == 1:
        return [F.one()]
    A = copy(Q)
    for i in range(n):
        A[i, n - 1] = F.one()
    b = [F.zero()] * n
    b[n - 1] = F.one()
    try:
        pi = A.solve_left(matrix(F, 1, n, b))
    except (ValueError, ZeroDivisionError, ArithmeticError):
        raise LineSymError("singular",
                           "the normalized generator is singular; the chain has no "
                           "unique stationary distribution")
    return [pi[0, j] for j in range(n)]


def _ctmc_solve(Q, F):
    """Symbolic stationary distribution, following ctmc_solve.m.

    Returns (pi, nConnComp, connComp) with pi a list over the full state
    space. A reducible generator is decomposed into weakly connected
    components, each solved on its own and the result renormalized, which is
    the documented (initial-vector free) MATLAB behaviour.
    """
    n = Q.nrows()
    if n == 1:
        return [F.one()], 1, [1]

    Q = _makeinfgen(Q)
    ncomp, comp = _weak_components(Q)
    if ncomp > 1:
        pi = [F.zero()] * n
        for c in range(1, ncomp + 1):
            idx = [i for i in range(n) if comp[i] == c]
            Qc = _makeinfgen(Q.matrix_from_rows_and_columns(idx, idx))
            pic, _, _ = _ctmc_solve(Qc, F)
            for k, i in enumerate(idx):
                pi[i] = pic[k]
        total = sum(pi)
        if total == 0:
            raise LineSymError("singular", "the decomposed solution sums to zero")
        return [_canonical(p / total, F) for p in pi], ncomp, comp

    if all(Q[i, j] == 0 for i in range(n) for j in range(n)):
        # No transitions at all: every distribution is stationary, and uniform
        # is as defensible as any other choice.
        return [F.one() / n] * n, 1, comp

    A, keep = _nonzero_support(Q)
    if not keep:
        raise LineSymError(
            "absorbing",
            "the infinitesimal generator has no recurrent state: every state was "
            "eliminated as absorbing. This generator admits no unique stationary "
            "distribution and usually means it is malformed, e.g. a class of "
            "transitions was dropped while building it.")

    block = _solve_block(A, F)
    pi = [F.zero()] * n
    for k, i in enumerate(keep):
        pi[i] = _canonical(block[k], F)
    return pi, 1, comp


def _common_fraction(pi, F):
    """Puts the whole vector over one common denominator.

    Clients compare symbolic results across codebases by substitution, but the
    (num, den) split is what makes a printed result readable and is what the
    MATLAB and Python callers expose as [pi_i, num, den].
    """
    if F is QQ:
        den = ZZ(1)
        for p in pi:
            den = den.lcm(QQ(p).denominator())
        return [str(QQ(p) * den) for p in pi], str(den)

    den = F.one().numerator()
    for p in pi:
        den = den.lcm(F(p).denominator())
    num = []
    for p in pi:
        q = F(p) * F(den)
        # Exact by construction: den is a common multiple of every denominator,
        # so each product is a polynomial and the division below is exact.
        num.append(q.numerator() // q.denominator())
    # One primitive normal form for the pair, so num and den scale together.
    scale = ZZ(1)
    for d in [_coeff_lcm_denominator(x) for x in num] + [_coeff_lcm_denominator(den)]:
        scale = scale.lcm(d)
    if scale != 1:
        num, den = [x * scale for x in num], den * scale
    coeffs = []
    for x in num:
        coeffs += _integer_coefficients(x)
    coeffs += _integer_coefficients(den)
    g = gcd(coeffs)
    if g and g != 1:
        num, den = [x / g for x in num], den / g
    if _leading_coefficient(den) < 0:
        num, den = [-x for x in num], -den
    return [str(x) for x in num], str(den)


# ---------------------------------------------------------------------------
# Handlers
# ---------------------------------------------------------------------------

def handle_health(_payload):
    return {"status": "ok",
            "timestamp": datetime.utcnow().isoformat() + "Z"}


def handle_ready(_payload):
    """Readiness: the CAS is reachable and gives the right answer.

    A liveness probe only proves the socket is open. This actually solves a
    two-state chain, so a Sage installation that imports but cannot compute is
    reported as not ready rather than discovered on the first real request.
    """
    checks = {}
    try:
        F, parser = _ctmc_parser(["x1", "x2"])
        Q = _parse_matrix([["-x1", "x1"], ["x2", "-x2"]], parser, F)
        pi, _, _ = _ctmc_solve(Q, F)
        checks["ctmc_solve"] = (str(pi[0]) == "x2/(x1 + x2)"
                                or pi[0] - F.gen(1) / (F.gen(0) + F.gen(1)) == 0)
    except Exception:  # noqa: BLE001 - readiness must report, not raise
        checks["ctmc_solve"] = False
    try:
        checks["symbolic_ring"] = bool(diff(var("t") ** 2, var("t")) == 2 * var("t"))
    except Exception:  # noqa: BLE001
        checks["symbolic_ring"] = False
    return {"status": "ok", "ready": all(checks.values()), "checks": checks}


def handle_info(_payload):
    return {"status": "ok",
            "name": "LINE symbolic REST API (SageMath)",
            "version": SERVICE_VERSION,
            "api_version": API_VERSION,
            "sage_version": version(),
            "python_version": sys.version.split()[0],
            "endpoints": sorted(list(META_ROUTES) + list(ROUTES)),
            "citation": CITATION,
            "notice": NOTICE}


def handle_citation(_payload):
    result = {"status": "ok", "notice": NOTICE}
    result.update(CITATION)
    return result


def handle_ctmc_solve(payload):
    symbols = payload.get("symbols", [])
    F, parser = _ctmc_parser(symbols)
    Q = _parse_matrix(payload.get("Q"), parser, F)
    pi, ncomp, comp = _ctmc_solve(Q, F)
    if payload.get("normalize", True):
        total = sum(pi)
        if total != 0:
            pi = [_canonical(p / total, F) for p in pi]
    num, den = _common_fraction(pi, F)
    return {"status": "ok",
            "pi": [_to_str(p) for p in pi],
            "num": num,
            "den": den,
            "nConnComp": int(ncomp),
            "connComp": [int(c) for c in comp]}


def handle_ctmc_sensitivity(payload):
    """Exact parametric sensitivity of a steady-state reward.

    The steady-state vector is solved symbolically and then differentiated
    with respect to theta, which is exact where @SolverCTMC/getSensitivity.m
    uses a central difference accurate to O(h^2). As there, dr/dtheta is taken
    to be zero: a reward whose rates themselves depend on theta needs the
    second term of Trivedi and Bobbio (2017), Eq. (9.83).
    """
    symbols = payload.get("symbols", [])
    theta = payload.get("theta")
    if not theta:
        raise LineSymError("input", "\"theta\" must name the parameter to differentiate")
    if str(theta) not in [str(s) for s in symbols]:
        raise LineSymError("input", "\"theta\" must be one of \"symbols\"")
    F, parser = _ctmc_parser(symbols)
    Q = _parse_matrix(payload.get("Q"), parser, F)
    pi, _, _ = _ctmc_solve(Q, F)
    total = sum(pi)
    if total != 0:
        pi = [_canonical(p / total, F) for p in pi]

    # Differentiation of a rational function is done in the symbolic ring: the
    # fraction field of a multivariate polynomial ring has no derivation that
    # covers the quotient rule in every Sage version.
    t = var(str(theta))
    pi_sr = [SR(str(p)) for p in pi]
    dpi = [diff(p, t) for p in pi_sr]

    result = {"status": "ok",
              "pi": [_to_str(p) for p in pi],
              "dpi": [_to_str(d.simplify_rational()) for d in dpi]}

    reward = payload.get("reward")
    if reward is not None:
        sr_parser = _sr_parser()
        r = [sr_parser.parse(e) for e in reward]
        if len(r) != len(pi_sr):
            raise LineSymError("input", "\"reward\" must have one entry per state")
        Er = sum(p * ri for p, ri in zip(pi_sr, r))
        S = sum(d * ri for d, ri in zip(dpi, r))
        result["Er"] = _to_str(Er.simplify_rational())
        result["S"] = _to_str(S.simplify_rational())
        if Er != 0:
            result["SS"] = _to_str((t / Er * S).simplify_rational())
    return result


_SIMPLIFY_FORMS = {
    "simplify": lambda e: e.simplify_full(),
    "factor": lambda e: e.factor(),
    "together": lambda e: e.combine(),
    "cancel": lambda e: e.simplify_rational(),
    "expand": lambda e: e.expand(),
    "latex": lambda e: latex(e),
}


def handle_simplify(payload):
    form = payload.get("form", "cancel")
    fn = _SIMPLIFY_FORMS.get(form)
    if fn is None:
        raise LineSymError("input", "unknown form %r; expected one of %s"
                           % (form, ", ".join(sorted(_SIMPLIFY_FORMS))))
    parser = _sr_parser()
    exprs = payload.get("exprs")
    if not isinstance(exprs, list):
        raise LineSymError("input", "\"exprs\" must be a list")
    out = []
    for e in exprs:
        parsed = parser.parse(e)
        try:
            out.append(_to_str(fn(parsed)))
        except (ArithmeticError, ValueError, TypeError) as err:
            raise LineSymError("compute", "cannot apply %s to %r: %s" % (form, e, err))
    return {"status": "ok", "results": out, "form": form}


def handle_diff(payload):
    parser = _sr_parser()
    exprs = payload.get("exprs")
    if not isinstance(exprs, list):
        raise LineSymError("input", "\"exprs\" must be a list")
    name = payload.get("var")
    if not name:
        raise LineSymError("input", "\"var\" must name the differentiation variable")
    order = int(payload.get("order", 1))
    if order < 1:
        raise LineSymError("input", "\"order\" must be a positive integer")
    t = var(str(name))
    out = [_to_str(diff(parser.parse(e), t, order).simplify_rational()) for e in exprs]
    return {"status": "ok", "results": out}


def handle_eval(payload):
    """Substitutes exact values for the symbols and evaluates.

    This is how the three codebases compare symbolic results with each other.
    Symbol numbering x1..xE follows event enumeration order, which is matched
    in practice but not structurally guaranteed, so comparing expression
    strings across rows is unsound; substituting each event's own rate and
    comparing numbers is not. The rationals are returned alongside the floats
    because the comparison is only exact in the former.
    """
    parser = _sr_parser()
    exprs = payload.get("exprs")
    if not isinstance(exprs, list):
        raise LineSymError("input", "\"exprs\" must be a list")
    values = payload.get("values", {})
    if not isinstance(values, dict):
        raise LineSymError("input", "\"values\" must be an object mapping symbol to value")
    subs = {}
    for k, v in values.items():
        subs[var(str(k))] = parser.parse(v)
    exact, numeric = [], []
    for e in exprs:
        val = parser.parse(e).subs(subs)
        try:
            val = val.simplify_rational()
        except (ArithmeticError, ValueError, TypeError):
            pass
        exact.append(_to_str(val))
        try:
            numeric.append(float(val))
        except (TypeError, ValueError):
            # A residual free symbol leaves the value non-numeric; report it as
            # such rather than coercing it to a number that is not one.
            numeric.append(None)
    return {"status": "ok", "exact": exact, "values": numeric}


def handle_fluid_odes(payload):
    """Jacobian, equilibria and LaTeX of a fluid ODE right hand side.

    The vector field is read in the symbolic ring because the fluid drift
    contains min() terms (see solver_fluid_symodes.m), which no polynomial
    ring can hold.
    """
    parser = _sr_parser()
    rhs_in = payload.get("rhs")
    vars_in = payload.get("vars")
    if not isinstance(rhs_in, list) or not rhs_in:
        raise LineSymError("input", "\"rhs\" must be a non-empty list")
    if not isinstance(vars_in, list) or not vars_in:
        raise LineSymError("input", "\"vars\" must be a non-empty list")
    if len(rhs_in) != len(vars_in):
        raise LineSymError("input", "\"rhs\" and \"vars\" must have the same length")
    xs = [var(str(v)) for v in vars_in]
    f = [parser.parse(e) for e in rhs_in]
    want = payload.get("want", ["jacobian", "latex"])
    result = {"status": "ok"}
    if "jacobian" in want:
        result["jacobian"] = [[_to_str(diff(fi, xj)) for xj in xs] for fi in f]
    if "latex" in want:
        result["latex"] = [_to_str(latex(fi)) for fi in f]
    if "equilibria" in want:
        try:
            sols = solve([fi == 0 for fi in f], xs, solution_dict=True)
        except (ValueError, TypeError, NotImplementedError) as err:
            raise LineSymError("compute", "cannot solve for equilibria: %s" % err)
        result["equilibria"] = [{str(k): _to_str(v) for k, v in s.items()} for s in sols]
    return result


# GET routes: metadata only, no request body.
META_ROUTES = {
    API_PREFIX + "/health": handle_health,
    API_PREFIX + "/ready": handle_ready,
    API_PREFIX + "/info": handle_info,
    API_PREFIX + "/citation": handle_citation,
}

ROUTES = {
    API_PREFIX + "/ctmc/solve": handle_ctmc_solve,
    API_PREFIX + "/ctmc/sensitivity": handle_ctmc_sensitivity,
    API_PREFIX + "/simplify": handle_simplify,
    API_PREFIX + "/diff": handle_diff,
    API_PREFIX + "/eval": handle_eval,
    API_PREFIX + "/fluid/odes": handle_fluid_odes,
}


# ---------------------------------------------------------------------------
# HTTP plumbing
# ---------------------------------------------------------------------------

class SymHandler(BaseHTTPRequestHandler):
    server_version = "line-sage-rest/1.0"

    def log_message(self, fmt, *args):
        if os.environ.get("LINE_SAGE_QUIET", "") not in ("", "0", "false"):
            return
        sys.stderr.write("[line-sage] %s - %s\n" % (self.address_string(), fmt % args))

    def _send(self, code, obj):
        body = json.dumps(obj).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        path = self.path.split("?")[0].rstrip("/") or "/"
        # The bare /health form is what a container HEALTHCHECK usually probes.
        handler = META_ROUTES.get(path) or META_ROUTES.get(API_PREFIX + path)
        if handler is not None:
            self._send(200, handler({}))
            return
        self._send(404, {"status": "error", "code": "route",
                         "message": "unknown endpoint %s" % path})

    def do_POST(self):
        path = self.path.split("?")[0].rstrip("/")
        handler = ROUTES.get(path) or META_ROUTES.get(path)
        if handler is None:
            self._send(404, {"status": "error", "code": "route",
                             "message": "unknown endpoint %s" % path})
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            length = 0
        if length > MAX_BODY_BYTES:
            self._send(413, {"status": "error", "code": "size",
                             "message": "request body exceeds %d bytes" % MAX_BODY_BYTES})
            return
        raw = self.rfile.read(length) if length else b"{}"
        try:
            payload = json.loads(raw.decode("utf-8")) if raw.strip() else {}
        except (ValueError, UnicodeDecodeError) as err:
            self._send(400, {"status": "error", "code": "json",
                             "message": "malformed JSON body: %s" % err})
            return
        if not isinstance(payload, dict):
            self._send(400, {"status": "error", "code": "json",
                             "message": "request body must be a JSON object"})
            return

        timeout = payload.get("timeout_s", DEFAULT_TIMEOUT_S)
        try:
            timeout = int(timeout)
        except (TypeError, ValueError):
            timeout = DEFAULT_TIMEOUT_S
        # The server forks per request, so this alarm is armed in the main
        # thread of a process that handles exactly this computation.
        if timeout > 0:
            signal.signal(signal.SIGALRM, _on_alarm)
            signal.alarm(timeout)
        try:
            self._send(200, handler(payload))
        except LineSymError as err:
            self._send(200, {"status": "error", "code": err.code, "message": err.message})
        except _Timeout:
            self._send(200, {"status": "error", "code": "timeout",
                             "message": "computation exceeded %d s; symbolic solves grow "
                                        "superpolynomially in the number of states" % timeout})
        except Exception as err:  # noqa: BLE001 - the service must not die on one request
            self._send(200, {"status": "error", "code": "internal",
                             "message": "%s: %s" % (type(err).__name__, err),
                             "traceback": traceback.format_exc()})
        finally:
            if timeout > 0:
                signal.alarm(0)


class _Timeout(Exception):
    pass


def _on_alarm(_signum, _frame):
    raise _Timeout()


class ForkingHTTPServer(socketserver.ForkingMixIn, HTTPServer):
    """One process per request, so a runaway solve can be killed by its alarm."""
    daemon_threads = False
    allow_reuse_address = True


def main(argv):
    port = int(os.environ.get("PORT", "8080"))
    host = os.environ.get("HOST", "0.0.0.0")
    args = list(argv[1:])
    while args:
        a = args.pop(0)
        if a in ("-p", "--port"):
            port = int(args.pop(0))
        elif a in ("-h", "--host"):
            host = args.pop(0)
        elif a in ("--help",):
            sys.stdout.write(__doc__)
            return 0
        else:
            sys.stderr.write("unknown argument: %s\n" % a)
            return 2
    server = ForkingHTTPServer((host, port), SymHandler)
    sys.stderr.write("[line-sage] listening on %s:%d (Sage %s)\n"
                     % (host, port, version()))
    sys.stderr.flush()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
