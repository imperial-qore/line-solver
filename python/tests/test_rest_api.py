"""Integration tests for the rest-api module (jline.rest.LineRestServer).

The server is a Java process; these tests drive it over HTTP only, with the
standard library, so nothing here depends on a JVM binding. They are skipped
unless the module has been packaged:

    cd io/rest-api && mvn clean package

Routes covered: health/ready/info, the solver catalogue and its detail and
methods routes, synchronous and asynchronous solve, validate, convert,
describe, the analysis routes, the metrics routes, and the API-key and
rate-limit filters.
"""

import base64
import json
import os
import shutil
import socket
import subprocess
import time
import urllib.error
import urllib.request
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
REST_JAR = REPO_ROOT / "io" / "rest-api" / "target" / "line-rest.jar"
JLINE_JAR = REPO_ROOT / "common" / "jline.jar"
MODEL_FILE = Path(__file__).resolve().parent / "fixtures" / "example.jsimg"

STARTUP_TIMEOUT = 90.0
SOLVE_TIMEOUT = 300.0

pytestmark = pytest.mark.skipif(
    shutil.which("java") is None
    or not REST_JAR.exists()
    or not JLINE_JAR.exists()
    or not MODEL_FILE.exists(),
    reason=(
        "needs java, common/jline.jar, tests/fixtures/example.jsimg and a packaged "
        "io/rest-api/target/line-rest.jar (cd io/rest-api && mvn clean package)"
    ),
)


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _request(base, path, payload=None, method=None, headers=None, timeout=SOLVE_TIMEOUT):
    """Issue a request and return (status, parsed_body). Never raises on 4xx/5xx."""
    url = base + path
    data = None
    hdrs = dict(headers or {})
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        hdrs.setdefault("Content-Type", "application/json")
    req = urllib.request.Request(url, data=data, headers=hdrs, method=method)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read().decode("utf-8")
            status = resp.status
    except urllib.error.HTTPError as e:
        raw = e.read().decode("utf-8")
        status = e.code
    try:
        return status, json.loads(raw)
    except ValueError:
        return status, raw


def _start_server(extra_args=()):
    port = _free_port()
    cmd = [
        "java", "-cp", "%s%s%s" % (REST_JAR, os.pathsep, JLINE_JAR),
        "jline.rest.LineRestServer", "--port", str(port),
    ]
    cmd.extend(extra_args)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    base = "http://127.0.0.1:%d/api/v1" % port

    deadline = time.time() + STARTUP_TIMEOUT
    while time.time() < deadline:
        if proc.poll() is not None:
            out = proc.stdout.read().decode("utf-8", "replace")
            raise RuntimeError("server exited with %s:\n%s" % (proc.returncode, out))
        try:
            status, _ = _request(base, "/health", timeout=2.0)
            if status == 200:
                return proc, base
        except (urllib.error.URLError, socket.timeout, ConnectionError):
            pass
        time.sleep(0.25)

    proc.kill()
    raise RuntimeError("server did not become healthy within %.0fs" % STARTUP_TIMEOUT)


def _stop_server(proc):
    proc.terminate()
    try:
        proc.wait(timeout=15)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait(timeout=15)


@pytest.fixture(scope="module")
def server():
    proc, base = _start_server()
    try:
        yield base
    finally:
        _stop_server(proc)


@pytest.fixture(scope="module")
def secured_server():
    """A second instance with authentication and a tight rate limit."""
    proc, base = _start_server(["--api-keys", "k1,k2", "--rate-limit", "3"])
    try:
        yield base
    finally:
        _stop_server(proc)


@pytest.fixture(scope="module")
def model_input():
    content = base64.b64encode(MODEL_FILE.read_bytes()).decode("ascii")
    return {"format": "jsimg", "content": content, "base64": True}


# --------------------------------------------------------------------------
# Health, readiness and server information
# --------------------------------------------------------------------------

def test_health_reports_healthy(server):
    status, body = _request(server, "/health")
    assert status == 200
    assert body["status"] == "healthy"
    assert isinstance(body["timestamp"], int)


def test_ready_reports_heap_usage(server):
    status, body = _request(server, "/ready")
    assert status in (200, 503)
    assert body["status"] in ("ready", "not_ready")
    assert 0 <= body["memoryUsagePercent"] <= 100


def test_info_pins_the_api_version(server):
    status, body = _request(server, "/info")
    assert status == 200
    assert body["apiVersion"] == "v1"
    assert body["name"] == "LINE Solver"
    assert "jsimg" in body["capabilities"]["inputFormats"]
    assert body["memory"]["heapUsed"] > 0


# --------------------------------------------------------------------------
# Solver catalogue
# --------------------------------------------------------------------------

def test_solvers_catalogue_advertises_ldes_not_des(server):
    status, body = _request(server, "/solvers")
    assert status == 200
    ids = [s["id"] for s in body["solvers"]]
    assert "ldes" in ids
    assert "des" not in ids
    for expected in ("mva", "ctmc", "fluid", "jmt", "nc", "ssa", "ln", "lqns"):
        assert expected in ids


def test_solver_detail_and_unknown_id(server):
    status, body = _request(server, "/solvers/mva")
    assert status == 200
    assert body["id"] == "mva"
    assert body["type"] == "analytical"
    assert "jsimg" in body["formats"]

    status, body = _request(server, "/solvers/nosuchsolver")
    assert status == 404
    assert body["solverId"] == "nosuchsolver"


def test_solver_methods_for_a_posted_model(server, model_input):
    status, body = _request(server, "/solvers/mva/methods", payload=model_input)
    assert status == 200
    assert body["solverId"] == "mva"
    assert isinstance(body["methods"], list) and body["methods"]


# --------------------------------------------------------------------------
# Synchronous solve
# --------------------------------------------------------------------------

def _assert_solved(body, solver):
    assert body["status"] == "completed", body.get("error")
    assert body["solver"] == solver
    assert body["runtime"] >= 0.0
    avg = body["results"]["avgTable"]
    assert avg["stations"] and avg["classes"]
    for metric in ("QLen", "Util", "RespT", "ResidT", "ArvR", "Tput"):
        assert metric in avg["metrics"]
        rows = avg["metrics"][metric]
        assert len(rows) == len(avg["stations"])
        assert all(len(r) == len(avg["classes"]) for r in rows)


def test_solve_with_mva(server, model_input):
    status, body = _request(
        server, "/models/solve",
        payload={"model": model_input, "solver": "mva", "analysis": "all"},
    )
    assert status == 200
    _assert_solved(body, "mva")
    sys_table = body["results"]["sysTable"]
    assert sys_table["chains"]
    assert "SysRespT" in sys_table["metrics"]
    assert "SysTput" in sys_table["metrics"]


def test_solve_with_ldes(server, model_input):
    """Regression: ldes was rejected by the validator while /solvers advertised it."""
    status, body = _request(
        server, "/models/solve",
        payload={"model": model_input, "solver": "ldes", "analysis": "avg"},
    )
    assert status == 200
    _assert_solved(body, "ldes")


def test_ldes_and_mva_agree_on_utilization(server, model_input):
    """A simulated and an analytical run of the same model must broadly agree."""
    _, mva = _request(
        server, "/models/solve",
        payload={"model": model_input, "solver": "mva", "analysis": "avg"},
    )
    _, ldes = _request(
        server, "/models/solve",
        payload={"model": model_input, "solver": "ldes", "analysis": "avg"},
    )
    mva_util = mva["results"]["avgTable"]["metrics"]["Util"]
    ldes_util = ldes["results"]["avgTable"]["metrics"]["Util"]
    assert len(mva_util) == len(ldes_util)
    for mrow, lrow in zip(mva_util, ldes_util):
        for m, l in zip(mrow, lrow):
            assert abs(m - l) < 0.15, "Util %.4f (mva) vs %.4f (ldes)" % (m, l)


def test_solve_rejects_des_and_unknown_solvers(server, model_input):
    for bad in ("des", "nosuchsolver"):
        status, body = _request(
            server, "/models/solve",
            payload={"model": model_input, "solver": bad},
        )
        assert status == 400, bad
        assert "Invalid solver" in json.dumps(body), bad


def test_solve_rejects_a_bad_analysis_type(server, model_input):
    status, body = _request(
        server, "/models/solve",
        payload={"model": model_input, "solver": "mva", "analysis": "nonsense"},
    )
    assert status == 400
    assert "Invalid analysis type" in json.dumps(body)


def test_solve_rejects_an_unparseable_model(server):
    status, body = _request(
        server, "/models/solve",
        payload={"model": {"format": "jsimg", "content": "<not-a-model/>"},
                 "solver": "mva"},
    )
    assert status in (400, 500)
    assert "error" in json.dumps(body).lower()


# --------------------------------------------------------------------------
# Validate, convert, describe
# --------------------------------------------------------------------------

def test_validate_accepts_the_model(server, model_input):
    # The route takes a bare ModelInput, not a {"model": ...} envelope.
    status, body = _request(server, "/models/validate", payload=model_input)
    assert status == 200
    assert body["valid"] is True


def test_convert_to_json_interchange(server, model_input):
    status, body = _request(
        server, "/models/convert",
        payload={"model": model_input, "targetFormat": "json"},
    )
    assert status == 200
    assert body["format"] == "json"
    assert json.loads(body["content"])


def test_describe_returns_the_interchange_document(server, model_input):
    status, body = _request(server, "/models/describe", payload=model_input)
    assert status == 200
    assert isinstance(body, dict) and body


# --------------------------------------------------------------------------
# Asynchronous solve and job management
# --------------------------------------------------------------------------

def test_async_solve_runs_to_completion(server, model_input):
    status, body = _request(
        server, "/models/solve/async",
        payload={"model": model_input, "solver": "mva", "analysis": "avg"},
    )
    assert status == 202
    job_id = body["jobId"]
    assert job_id

    deadline = time.time() + SOLVE_TIMEOUT
    job = None
    while time.time() < deadline:
        _, job = _request(server, "/jobs/%s" % job_id)
        if job["status"] in ("completed", "failed"):
            break
        time.sleep(0.5)

    assert job is not None and job["status"] == "completed", job
    _, listing = _request(server, "/jobs")
    assert any(j["jobId"] == job_id for j in listing["jobs"])


def test_unknown_job_is_not_found(server):
    status, _ = _request(server, "/jobs/does-not-exist")
    assert status == 404


# --------------------------------------------------------------------------
# Analysis routes
# --------------------------------------------------------------------------

def test_bottleneck_analysis(server, model_input):
    status, body = _request(
        server, "/analysis/bottleneck",
        payload={"model": model_input, "solver": "mva", "includeAll": True},
    )
    assert status == 200
    assert body["status"] == "completed", body.get("error")
    # allStations is only returned when includeAll is requested.
    assert body["allStations"]
    assert isinstance(body["bottlenecks"], list)
    for st in body["allStations"]:
        assert 0.0 <= st["utilization"] <= 1.0 + 1e-9


def test_whatif_sweep(server, model_input):
    status, body = _request(
        server, "/analysis/whatif",
        payload={
            "model": model_input,
            "solver": "mva",
            # The class field is bound to the JSON key "class", not "className".
            "parameters": [{"type": "arrival_rate", "station": "Source 1",
                            "class": "Class A", "values": [0.1, 0.2]}],
            "metrics": ["Util"],
        },
    )
    assert status == 200
    assert body["status"] == "completed", body.get("error")
    assert body["totalPoints"] == 2
    assert body["completedPoints"] == 2
    assert body["failedPoints"] == 0
    # Each point carries the full result tables: WhatIfRequest.metrics has a
    # getter no caller reads, so the requested filter is silently ignored.
    for point in body["results"]:
        assert "arrivalRate" in json.dumps(point["parameters"])
        assert "Util" in point["metrics"]["avgTable"]["metrics"]
    # ParameterSweep.getDescription() keys each point as station.class.property.
    key = "Source 1.Class A.arrivalRate"
    assert sorted(p["parameters"][key] for p in body["results"]) == [0.1, 0.2]


# --------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------

def test_prometheus_metrics_are_exposed(server, model_input):
    _request(server, "/models/solve",
             payload={"model": model_input, "solver": "mva", "analysis": "avg"})
    status, body = _request(server, "/metrics")
    assert status == 200
    assert "line_requests_total" in body
    assert "line_solve_requests_total" in body


def test_json_metrics_count_requests(server):
    status, body = _request(server, "/metrics/json")
    assert status == 200
    assert body["requests"]["total"] > 0
    assert body["requests"]["solve"] > 0


# --------------------------------------------------------------------------
# Security filters
# --------------------------------------------------------------------------

def test_api_key_is_required_and_health_is_exempt(secured_server):
    status, _ = _request(secured_server, "/solvers")
    assert status == 401

    status, _ = _request(secured_server, "/solvers", headers={"X-API-Key": "k1"})
    assert status == 200

    # A missing key is 401 (UNAUTHORIZED), a wrong one 403 (FORBIDDEN).
    status, _ = _request(secured_server, "/solvers", headers={"X-API-Key": "wrong"})
    assert status == 403

    status, _ = _request(secured_server, "/health")
    assert status == 200


def test_rate_limit_returns_429(secured_server):
    seen = []
    for _ in range(8):
        status, _ = _request(secured_server, "/info", headers={"X-API-Key": "k2"})
        seen.append(status)
    assert 429 in seen, seen
