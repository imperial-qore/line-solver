"""
Native Python client for the LINE symbolic REST service (SageMath backend).

Same JSON protocol as the JAR (``jline.api.sym``) and MATLAB (``SAGE.m``)
clients; see ``sage/server.py`` for the endpoints. Uses only the standard
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
import urllib.error
import urllib.request

DOCKER_IMAGES = ("imperialqore/line-sage-rest:latest", "imperialqore/line-sage-rest")
URL_ENV = "LINE_SAGE_URL"
MODE_ENV = "LINE_SYMBOLIC_BACKEND"
PROBE_PORTS = (8085, 8080)
STARTUP_TIMEOUT_S = 120
DEFAULT_TIMEOUT_S = 300

_configured_url = None
_resolved = None
_started_container = None


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

    def info(self):
        return self._request("/api/v1/info", timeout=5)

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
        return engine if engine.is_available() else None

    env = os.environ.get(URL_ENV, "").strip()
    if env:
        engine = SageRestEngine(env)
        if engine.is_available():
            return engine

    if _resolved is not None and _resolved.is_available():
        return _resolved

    for port in PROBE_PORTS:
        engine = SageRestEngine("http://localhost:%d" % port)
        if _is_sage_service(engine):
            _resolved = engine
            return engine

    image = _find_image() if req.lower() in ("", "auto", "true", "sage") else req
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
            return engine
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
