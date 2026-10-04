"""
Native Python client for the LINE symbolic REST service (SageMath backend).

Same JSON protocol as the JAR (``jline.api.sym``) and MATLAB (``SAGE.m``)
clients; see ``io/sage/server.py`` for the endpoints. Uses only the standard
library, so it adds no dependency to the native Python implementation and in
particular no JVM.

Two reasons to route symbolic work through it rather than sympy: it is
exact over the rational function field and considerably faster there on a
state space of any size, and it gives all three codebases one normal form, so
a symbolic result computed in one can be compared with another.

Backend resolution, in order:

1. an explicit URL passed in, or set with :func:`set_backend`;
2. the ``LINE_SAGE_URL`` environment variable;
3. a line-sage-rest service already listening on a conventional port;
4. a container started here from a locally present image;
5. nothing, in which case the caller falls back to sympy.

Step 3 checks identity through ``/api/v1/info`` rather than trusting the port:
every imperialqore line-*-rest service listens on 8080 by convention, so a
health probe alone would accept the LQNS service.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import atexit
import json
import os
import socket
import subprocess
import sys
import urllib.error
import urllib.request

DOCKER_IMAGES = ("imperialqore/line-sage-rest:latest", "imperialqore/line-sage-rest")
URL_ENV = "LINE_SAGE_URL"
MODE_ENV = "LINE_SYMBOLIC_BACKEND"
PROBE_PORTS = (8085, 8080)
STARTUP_TIMEOUT_S = 120
DEFAULT_TIMEOUT_S = 300

# The POST routes this client calls; a service that does not serve them all is
# older than the client and is refused by name. See missing_routes.
REQUIRED_ROUTES = ("/api/v1/ctmc/solve", "/api/v1/ctmc/sensitivity",
                   "/api/v1/ctmc/measures", "/api/v1/ctmc/passage",
                   "/api/v1/simplify", "/api/v1/diff", "/api/v1/eval",
                   "/api/v1/fluid/odes")

# A service that answers /health can still be unable to COMPUTE. The
# line-sage-rest image ships a FLINT built for CPUs that have BMI2 and ADX; on
# an older host the first multi-limb exact operation raises SIGILL, the worker
# dies mid-request and the call returns no bytes at all. Health is pure Python
# and keeps answering, so it cannot see this. The canary is the
# weighted-average softmin form -- what the fluid export actually sends -- and
# is the smallest expression observed to trigger it; see
# _kb/11-conventions-and-gotchas.md.
CANARY_EXPR = "(x*exp(-x) + exp(-1))/(exp(-x) + exp(-1))"
CANARY_ARG = "0.68999999999999995"
CANARY_VALUE = 0.82116556904906557

_configured_url = None
_resolved = None
_started_container = None
_usable = {}


class SageError(RuntimeError):
    """Error reported by the service, carrying its machine-readable code."""

    def __init__(self, code, message):
        RuntimeError.__init__(self, "line-sage-rest [%s]: %s" % (code, message))
        self.code = code
        self.message = message


class SageRestEngine:
    """Client for one line-sage-rest service."""

    def __init__(self, base_url, timeout_s=DEFAULT_TIMEOUT_S):
        self.base_url = base_url.rstrip("/")
        self.timeout_s = timeout_s

    # -- plumbing ----------------------------------------------------------

    def _request(self, path, payload=None, timeout=None):
        url = self.base_url + path
        timeout = self.timeout_s if timeout is None else timeout
        if payload is None:
            req = urllib.request.Request(url, method="GET")
        else:
            body = dict(payload)
            if self.timeout_s > 0:
                body.setdefault("timeout_s", self.timeout_s)
            data = json.dumps(body).encode("utf-8")
            req = urllib.request.Request(
                url, data=data, method="POST",
                headers={"Content-Type": "application/json; charset=utf-8"})
        with urllib.request.urlopen(req, timeout=max(timeout + 30, 60)) as resp:
            response = json.loads(resp.read().decode("utf-8"))
        if response.get("status") != "ok":
            raise SageError(response.get("code", "error"),
                            response.get("message", "unspecified error"))
        return response

    def is_available(self):
        try:
            self._request("/api/v1/health", timeout=5)
            return True
        except (urllib.error.URLError, OSError, ValueError, SageError):
            return False

    def is_usable(self):
        """True if the service answers AND can evaluate.

        Verdicts are cached per URL: this costs one small request the first
        time a service is considered, and nothing after.
        """
        cached = _usable.get(self.base_url)
        if cached is not None:
            return cached
        ok = False
        try:
            r = self._request("/api/v1/eval",
                              {"exprs": [CANARY_EXPR], "values": {"x": CANARY_ARG},
                               "timeout_s": 30},
                              timeout=30)
            values = r.get("values") or []
            ok = (len(values) == 1 and values[0] is not None
                  and abs(float(values[0]) - CANARY_VALUE) < 1e-9)
        except (urllib.error.URLError, OSError, ValueError, TypeError, SageError):
            # A dead worker closes the connection without a reply, which
            # surfaces as a transport error rather than a service one. Either
            # way the backend cannot serve us.
            ok = False
        if not ok:
            import sys
            sys.stderr.write(
                "[LINE] Ignoring symbolic backend at %s: it did not return the "
                "usability canary. On a CPU without BMI2/ADX the image's FLINT "
                "raises SIGILL mid-request.\n" % self.base_url)
        _usable[self.base_url] = ok
        return ok

    def info(self):
        return self._request("/api/v1/info", timeout=5)

    def missing_routes(self):
        """Name the routes this client calls that the service does not serve.

        A SERVICE CAN BE HEALTHY AND STILL TOO OLD, and that case used to reach
        the caller as a request failure rather than as a resolution failure.
        The roster ran a line-sage-rest image predating the ctmc/measures and
        ctmc/passage routes on 2026-09-11: it answered /api/v1/health, passed
        the is_usable arithmetic canary, and then failed every call to those
        routes with "unknown endpoint /api/v1/ctmc/measures" -- an error that
        reads as a defect in the solver rather than as an out-of-date
        container. /api/v1/info lists what the running server actually routes,
        so asking it turns that into a named refusal at resolution time.

        Returns:
            the missing routes, empty when the service serves them all or is
            too old to report its route list at all
        """
        served = self.info().get("endpoints")
        if not isinstance(served, list):
            # A service that does not list its routes cannot be interrogated;
            # it is left to the request itself to fail, as before.
            return []
        return [r for r in REQUIRED_ROUTES if r not in served]

    # -- operations --------------------------------------------------------

    def solve_ctmc(self, Q, symbols, normalize=True):
        """Symbolic stationary distribution, pi Q = 0 with sum(pi) = 1.

        Args:
            Q: square nested sequence of expression strings, row major
            symbols: the symbols occurring in Q
            normalize: renormalize the solution, as ctmc_solve does

        Returns:
            dict with pi, num, den, nConnComp and connComp
        """
        payload = {"Q": _as_expr_matrix(Q), "symbols": [str(s) for s in symbols],
                   "normalize": bool(normalize)}
        return self._request("/api/v1/ctmc/solve", payload)

    def ctmc_sensitivity(self, Q, symbols, theta, reward=None):
        """Exact d(pi)/d(theta) and, given a reward, d(E[r])/d(theta)."""
        payload = {"Q": _as_expr_matrix(Q), "symbols": [str(s) for s in symbols],
                   "theta": str(theta)}
        if reward is not None:
            payload["reward"] = [_as_expr(r) for r in reward]
        return self._request("/api/v1/ctmc/sensitivity", payload)

    def ctmc_measures(self, Q, symbols, weights, names=None, ratios=None,
                      normalize=True):
        """Stationary distribution and a set of weighted sums of it.

        Every mean SolverCTMC reports is a NUMERIC linear functional of pi, so
        the caller builds one weight vector per measure with no algebra at all
        and this one round trip returns each ``w . pi`` as a rational function.
        A throughput is ``pi . depRates``, a queue length is
        ``pi . stateSpaceAggr``, a marginal probability is ``pi . indicator``
        and a mean reward is ``pi . r``.

        Args:
            Q: square nested sequence of expression strings, row major
            symbols: the symbols occurring in Q
            weights: one row per measure, each of length len(Q). Entries are
                parsed in the same field as Q, so an expression is admissible
                and not only a number.
            names: a name per weight row, defaulting to w0, w1, ...
            ratios: dicts {name, num, den} whose num/den index into weights.
                A ratio whose denominator measure is identically zero comes
                back with expr "undefined" and a reason, never 0.
            normalize: renormalize pi, as ctmc_solve does

        Returns:
            dict with pi, num, den, nConnComp, connComp, measures and ratios
        """
        payload = {"Q": _as_expr_matrix(Q), "symbols": [str(s) for s in symbols],
                   "weights": _as_expr_rows(weights), "normalize": bool(normalize)}
        if names is not None:
            payload["names"] = [str(n) for n in names]
        if ratios:
            payload["ratios"] = [dict(r) for r in ratios]
        return self._request("/api/v1/ctmc/measures", payload)

    def ctmc_passage(self, S, s0, alpha, atom="0", symbols=(), svar="s",
                     want=("moments",), nmax=1):
        """First passage time into a target set: the transform and the moments.

        With A the complement of the target, S = Q(A,A), s0 = -S*1 and alpha
        the initial law restricted to A,

            L(s) = alpha (sI - S)^-1 s0 + atom
            (-S) M(n) = n M(n-1),  M(0) = 1

        Moments come from the recursion, not from differentiating L: it
        carries no transform symbol, stays in the same fraction field, and is
        nmax right solves against one matrix.

        Args:
            S: the sub-generator on the non-target states, row major
            s0: the exit vector -S*1
            alpha: the initial law restricted to the non-target states
            atom: the initial mass already inside the target set
            symbols: the rate symbols occurring in S
            svar: the transform symbol, adjoined only when a transform is asked
                for and required to differ from every rate symbol
            want: any of lst, lstall, moments, momall
            nmax: highest moment order

        Returns:
            dict with the requested items plus unreachable
        """
        payload = {"S": _as_expr_matrix(S), "s0": _as_expr_list(s0),
                   "alpha": _as_expr_list(alpha), "atom": _as_expr(atom),
                   "symbols": [str(x) for x in symbols], "svar": str(svar),
                   "want": [str(w) for w in want], "nmax": int(nmax)}
        return self._request("/api/v1/ctmc/passage", payload)

    def simplify(self, exprs, form="cancel"):
        """Rewrite expressions: simplify, factor, together, cancel, expand, latex."""
        payload = {"exprs": [_as_expr(e) for e in exprs], "form": form}
        return self._request("/api/v1/simplify", payload)["results"]

    def diff(self, exprs, variable, order=1):
        """Differentiate expressions with respect to ``variable``."""
        payload = {"exprs": [_as_expr(e) for e in exprs],
                   "var": str(variable), "order": int(order)}
        return self._request("/api/v1/diff", payload)["results"]

    def eval(self, exprs, assignment):
        """Substitute values for symbols and evaluate.

        This is how symbolic results are compared across codebases: symbol
        numbering follows event enumeration order, so comparing expression
        text is unsound while comparing substituted values is not.

        Returns:
            (values, exact) with values a list of floats, None where a free
            symbol remains, and exact the same values as strings
        """
        payload = {"exprs": [_as_expr(e) for e in exprs],
                   # Sent as text so the server reads the decimal exactly.
                   "values": dict((str(k), _as_expr(v))
                                  for k, v in assignment.items())}
        r = self._request("/api/v1/eval", payload)
        return r["values"], r["exact"]

    def fluid_odes(self, rhs, variables, want=("jacobian", "latex")):
        """Jacobian, LaTeX form and equilibria of a fluid drift."""
        payload = {"rhs": [_as_expr(e) for e in rhs],
                   "vars": [str(v) for v in variables],
                   "want": list(want)}
        return self._request("/api/v1/fluid/odes", payload)


def _as_expr(e):
    """Expression as text the service reads exactly.

    A float is written with full precision and read server side as the exact
    rational with that decimal expansion, never as the binary double.
    """
    if isinstance(e, str):
        return e
    if isinstance(e, float):
        if e == int(e) and abs(e) < 1e15:
            return "%d" % int(e)
        return repr(e)
    if isinstance(e, int):
        return "%d" % e
    return str(e)


def _as_expr_matrix(Q):
    try:
        rows = Q.tolist()
    except AttributeError:
        rows = [list(row) for row in Q]
    return [[_as_expr(e) for e in row] for row in rows]


def _as_expr_rows(rows):
    """Encodes a possibly non-square list of rows, e.g. a weight block.

    _as_expr_matrix is for the generator and takes a square matrix; a weight
    block is (number of measures) by (number of states) and is rarely square.
    """
    try:
        rows = rows.tolist()
    except AttributeError:
        rows = [list(r) for r in rows]
    return [[_as_expr(e) for e in r] for r in rows]


def _as_expr_list(values):
    try:
        values = values.tolist()
    except AttributeError:
        values = list(values)
    return [_as_expr(v) for v in values]


# ---------------------------------------------------------------------------
# Backend resolution
# ---------------------------------------------------------------------------

def set_backend(url_or_mode):
    """Pin the backend for this process.

    Args:
        url_or_mode: a service URL, "auto" to search, or None/"none" to
            disable the service and keep sympy
    """
    global _configured_url, _resolved
    _configured_url = url_or_mode
    _resolved = None


def backend_mode():
    """Requested backend: "auto" (default), "sympy", "sage", or a URL.

    Read from :func:`set_backend` if called, else from ``LINE_SYMBOLIC_BACKEND``.
    """
    if _configured_url is not None:
        return _configured_url
    return os.environ.get(MODE_ENV, "auto")


def prefers_sage():
    """True if the caller asked for the service rather than sympy.

    ``auto`` keeps sympy: switching engine changes the printed normal form of
    every symbolic result, so it is opt in.
    """
    mode = str(backend_mode()).strip().lower()
    return mode == "sage" or mode.startswith("http://") or mode.startswith("https://")


def resolve(requested=None):
    """Resolve an engine, starting a container if necessary.

    Args:
        requested: "" or "auto" to search, a URL, "none" to disable, or an
            image name to run

    Returns:
        a :class:`SageRestEngine`, or None if no backend could be resolved
    """
    global _resolved, _started_container
    req = (requested if requested is not None else backend_mode())
    req = "auto" if req is None else str(req).strip()
    if req.lower() in ("none", "off", "sympy"):
        return None
    if req.startswith("http://") or req.startswith("https://"):
        engine = SageRestEngine(req)
        return (engine if engine.is_available() and _serves_required_routes(engine)
                and engine.is_usable() else None)

    env = os.environ.get(URL_ENV, "").strip()
    if env:
        engine = SageRestEngine(env)
        if engine.is_available() and _serves_required_routes(engine) and engine.is_usable():
            return engine

    if _resolved is not None and _resolved.is_available() and _resolved.is_usable():
        return _resolved

    for port in PROBE_PORTS:
        engine = SageRestEngine("http://localhost:%d" % port)
        if _is_sage_service(engine) and _serves_required_routes(engine) and engine.is_usable():
            _resolved = engine
            return engine

    search = req.lower() in ("", "auto", "true", "sage")
    image = _find_image() if search else req
    # Auto-pull only on an explicit opt-in: the "sage" keyword or a named image.
    # Bare "auto"/"true"/"" keep sympy unless the image is already local, so
    # leaving the backend on auto never triggers a pull.
    if image is None and req.lower() == "sage":
        image = _pull_image(DOCKER_IMAGES[0])
    elif not search and image is not None and not _has_local_image(image):
        image = _pull_image(image)
    if image is None:
        return None
    engine = _start_container(image)
    _resolved = engine
    return engine


def _is_sage_service(engine):
    try:
        return "sage_version" in engine.info()
    except (urllib.error.URLError, OSError, ValueError, SageError):
        return False


def _serves_required_routes(engine):
    """True if the service is new enough to use, saying what is missing if not.

    A service that is line-sage-rest, healthy and arithmetically sound can
    still be OLDER THAN THE CLIENT; accepting it defers the failure to the
    first request needing a route it does not have. Refusing here sends the
    caller down the same path as "no service at all" -- sympy, or a skip --
    and says why. See SageRestEngine.missing_routes for the run that motivated
    it.
    """
    try:
        missing = engine.missing_routes()
    except (urllib.error.URLError, OSError, ValueError, SageError):
        return False
    if not missing:
        return True
    sys.stderr.write(
        "[LINE] Ignoring the line-sage-rest service at %s: it is older than this "
        "client and does not serve %s. Refresh it with 'docker pull %s'.\n"
        % (engine.base_url, ", ".join(missing), DOCKER_IMAGES[0]))
    return False


def _has_local_image(image):
    from ..io import docker_util
    return docker_util.has_local_image(image)


def _pull_image(target):
    """Pull the SAGE image if the Docker storage location has room.

    Returns the tag on success, else None. Storage-guarded via docker_util (the
    same guard as the LQNS/JMT wrappers), so an opt-in symbolic request never
    silently fills the Docker disk; on refusal the caller keeps sympy.
    """
    import sys
    from ..io import docker_util
    if not docker_util.has_storage_for(target):
        print("[LINE] Skipping docker pull of %s: insufficient free space at the Docker "
              "storage location; keeping the native symbolic backend." % target, file=sys.stderr)
        return None
    print("[LINE] Pulling Docker image %s (this may take a while)..." % target)
    if docker_util.pull(target) and docker_util.has_local_image(target):
        return target
    return None


def _find_image():
    for image in DOCKER_IMAGES:
        try:
            out = subprocess.run(["docker", "images", "-q", image],
                                 capture_output=True, text=True, timeout=30)
        except (OSError, subprocess.SubprocessError):
            return None
        if out.returncode == 0 and out.stdout.strip():
            return image
    return None


def _free_port():
    s = socket.socket()
    try:
        s.bind(("", 0))
        return s.getsockname()[1]
    finally:
        s.close()


def _start_container(image):
    """Run the service and wait for it to answer, or return None."""
    global _started_container
    port = _free_port()
    name = "line-sage-rest-%d" % port
    try:
        out = subprocess.run(
            ["docker", "run", "-d", "--rm", "--name", name,
             "-p", "%d:8080" % port, image],
            capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0 or not out.stdout.strip():
        return None
    _started_container = name
    atexit.register(stop_container)

    engine = SageRestEngine("http://localhost:%d" % port)
    import time
    deadline = time.time() + STARTUP_TIMEOUT_S
    while time.time() < deadline:
        if engine.is_available():
            if _serves_required_routes(engine) and engine.is_usable():
                return engine
            # Booted, but its arithmetic dies on this CPU. Keeping it running
            # would only cost memory, and returning it would hand the caller a
            # backend that kills every request.
            stop_container()
            return None
        time.sleep(0.5)
    stop_container()
    return None


def stop_container():
    """Stop the container started by this process, if any."""
    global _started_container, _resolved
    if _started_container:
        try:
            subprocess.run(["docker", "stop", "-t", "1", _started_container],
                           capture_output=True, timeout=60)
        except (OSError, subprocess.SubprocessError):
            pass
        _started_container = None
        _resolved = None


def require(requested=None):
    """Resolve a backend or raise with how to get one."""
    engine = resolve(requested)
    if engine is None:
        raise RuntimeError(
            "No symbolic backend is available. Start one with\n"
            "  docker run -d -p 8080:8080 %s\n"
            "point %s at a running service, or pass its URL to set_backend()."
            % (DOCKER_IMAGES[0], URL_ENV))
    return engine
