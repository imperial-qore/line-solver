"""
The MCP server's tool surface and its two output formats.

Covers what doc/latex/LINE-mcp.tex promises a client: the nine tools, the text and json
shapes of a successful analysis, and -- the part that had no coverage and no consumer
until sweep_parameter grew one -- that a run which never reaches a result answers in the
format the caller asked for. sweep_parameter reads its own points back with json.loads,
so a plain-text error there surfaced as "Expecting value: line 1 column 1" instead of the
reason the point failed.
"""

import json

import pytest

pytest.importorskip("mcp", reason="the MCP SDK is what the server is written against")

from line_solver.mcp_server import (analyze_closed_network, analyze_queue, clear_session,
                                    compare_solvers, get_example_code, list_examples,
                                    solve_line_model, sweep_parameter, visualize_model)

MM1 = """model = Network('M/M/1')
source = Source(model, 'Source')
queue = Queue(model, 'Queue', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')
oclass = OpenClass(model, 'Class1')
source.setArrival(oclass, Exp(1.0))
queue.setService(oclass, Exp(2.0))
model.link(Network.serial_routing([source, queue, sink]))
"""

# The routing form the guide shows: one hop per set() call. Indexing a routing matrix by
# node (P[jobclass, src, dst] = p) is not a native-Python idiom and raises ValueError.
WEBAPP = """model = Network('WebApp')
source = Source(model, 'Arrivals')
lb = Queue(model, 'LoadBalancer', SchedStrategy.PS)
app = Queue(model, 'AppServer', SchedStrategy.FCFS)
sink = Sink(model, 'Departures')
oclass = OpenClass(model, 'Requests')
source.setArrival(oclass, Exp(5.0))
lb.setService(oclass, Exp(100.0))
app.setService(oclass, Exp(10.0))
P = model.initRoutingMatrix()
P.set(oclass, oclass, source, lb, 1.0)
P.set(oclass, oclass, lb, app, 1.0)
P.set(oclass, oclass, app, sink, 1.0)
model.link(P)
"""


def _queue_row(payload):
    """The Queue station's record out of a json reply."""
    return next(r for r in payload["results"] if r["Station"] == "Queue")


def test_nine_tools_are_registered():
    from line_solver.mcp_server import mcp

    names = {t.name for t in mcp._tool_manager.list_tools()}
    assert names == {"solve_line_model", "list_examples", "get_example_code", "analyze_queue",
                     "analyze_closed_network", "sweep_parameter", "visualize_model",
                     "compare_solvers", "clear_session"}


def test_analyze_queue_text_carries_the_model_header():
    out = analyze_queue(arrival_rate=2.0, service_rate=5.0)
    assert "Model: M/M/1" in out
    assert "utilization (rho)     = 0.400000" in out
    assert "solver                = MVA" in out


def test_analyze_queue_json_shape_and_values():
    payload = json.loads(analyze_queue(arrival_rate=2.0, service_rate=5.0, format="json"))
    assert payload["metadata"] == {"model": "M/M/1", "arrival_rate": 2.0, "service_rate": 5.0,
                                   "servers": 1, "utilization": 0.4, "solver": "MVA"}
    row = _queue_row(payload)
    assert row["QLen"] == pytest.approx(2.0 / 3.0, rel=1e-9)
    assert row["Util"] == pytest.approx(0.4, rel=1e-9)
    assert row["RespT"] == pytest.approx(1.0 / 3.0, rel=1e-9)


def test_analyze_queue_finite_capacity_multiserver_notation():
    out = analyze_queue(arrival_rate=2.0, service_rate=5.0, servers=2, queue_capacity=10,
                        solver="CTMC")
    assert "Model: M/M/2/10" in out


def test_analyze_closed_network_json():
    payload = json.loads(analyze_closed_network(population=10, think_time=5.0, service_rate=2.0,
                                                format="json"))
    assert payload["metadata"]["model"] == "Closed / N=10 / k=1"
    # Utilisation at the Delay is the mean number thinking, so the two stations hold N jobs.
    total = sum(r["QLen"] for r in payload["results"])
    assert total == pytest.approx(10.0, rel=1e-9)


@pytest.mark.parametrize("kwargs, fragment", [
    (dict(arrival_rate=6.0, service_rate=5.0), "system is unstable"),
    (dict(arrival_rate=0.0, service_rate=5.0), "arrival_rate must be > 0"),
    (dict(arrival_rate=1.0, service_rate=0.0), "service_rate must be > 0"),
    (dict(arrival_rate=1.0, service_rate=5.0, servers=0), "servers must be >= 1"),
    (dict(arrival_rate=1.0, service_rate=5.0, solver="BOGUS"), "Unknown solver"),
    (dict(arrival_rate=1.0, service_rate=5.0, scheduling="BOGUS"), "Unknown scheduling"),
])
def test_analyze_queue_errors_honour_the_requested_format(kwargs, fragment):
    text = analyze_queue(**kwargs)
    assert fragment in text
    with pytest.raises(json.JSONDecodeError):
        json.loads(text)

    payload = json.loads(analyze_queue(format="json", **kwargs))
    assert "results" not in payload
    assert fragment in payload["error"]


@pytest.mark.parametrize("kwargs, fragment", [
    (dict(population=0, think_time=1.0, service_rate=2.0), "population must be >= 1"),
    (dict(population=2, think_time=-1.0, service_rate=2.0), "think_time must be >= 0"),
    (dict(population=2, think_time=1.0, service_rate=0.0), "service_rate must be > 0"),
    (dict(population=2, think_time=1.0, service_rate=2.0, solver="BOGUS"), "Unknown solver"),
])
def test_analyze_closed_network_errors_honour_the_requested_format(kwargs, fragment):
    assert fragment in analyze_closed_network(**kwargs)
    payload = json.loads(analyze_closed_network(format="json", **kwargs))
    assert fragment in payload["error"]


def test_solve_line_model_error_is_json_when_json_was_asked_for():
    payload = json.loads(solve_line_model(code="raise RuntimeError('boom')", format="json"))
    assert "RuntimeError: boom" in payload["error"]
    assert "RuntimeError: boom" in solve_line_model(code="raise RuntimeError('boom')")


def test_sweep_reports_the_reason_a_point_failed():
    # servers=1 saturates (rho = 4/3); 2 and 3 are stable. The failing point must name
    # instability rather than the json parse error it used to surface as.
    out = sweep_parameter(param="servers", values="1:3:1", arrival_rate=2.0, service_rate=1.5)
    assert "servers=1.0: Error: system is unstable" in out
    assert "Expecting value" not in out
    assert out.count("Queue") == 2  # only the two stable points produced rows


def test_sweep_with_no_successful_run_lists_every_reason():
    out = sweep_parameter(param="arrival_rate", values="1,2", service_rate=5.0, solver="BOGUS")
    assert out.startswith("No successful runs.")
    assert out.count("Unknown solver 'BOGUS'") == 2


def test_sweep_open_and_closed_produce_a_row_per_point():
    out = sweep_parameter(param="arrival_rate", values="0.5:1.5:0.5", service_rate=5.0)
    assert [v in out for v in ("0.5", "1.0", "1.5")] == [True, True, True]
    closed = sweep_parameter(param="population", values="1,2,5", model_type="closed",
                             think_time=1.0, service_rate=2.0)
    assert "population" in closed.splitlines()[0]


def test_sweep_rejects_a_bad_range_and_an_unknown_parameter():
    assert "Error parsing values" in sweep_parameter(param="arrival_rate", values="1:2")
    assert "Unknown parameter" in sweep_parameter(param="nosuch", values="1,2")
    assert "maximum 100 sweep points" in sweep_parameter(param="arrival_rate", values="0:200:1")


def test_solve_line_model_session_persists_the_namespace():
    first = solve_line_model(code=MM1 + "avg_table = MVA(model).avg_table()\n",
                             session_id="pytest_mcp")
    assert "Error" not in first
    second = solve_line_model(code="avg_table = NC(model).avg_table()\n", format="json",
                              session_id="pytest_mcp")
    assert _queue_row(json.loads(second))["Util"] == pytest.approx(0.5, rel=1e-9)
    assert "cleared" in clear_session(session_id="pytest_mcp")
    assert "not found" in clear_session(session_id="pytest_mcp")


def test_compare_solvers_agrees_on_a_product_form_model():
    out = compare_solvers(code=MM1, solvers="MVA,NC,MAM")
    assert [s in out for s in ("MVA", "NC", "MAM")] == [True, True, True]
    # One header plus one row per station per solver.
    assert out.count("Queue") == 3


def test_visualize_model_emits_mermaid_for_the_documented_routing_idiom():
    out = visualize_model(code=WEBAPP)
    assert out.startswith("graph LR")
    for hop in ("N0 --> N1", "N1 --> N2", "N2 --> N3"):
        assert hop in out
    assert "Requests (open)" in out


def test_example_gallery_is_browsable():
    listing = list_examples()
    assert "gallery_mm1" in listing
    code = get_example_code(name="gallery_mm1")
    assert "Network('M/M/1')" in code
    assert "not found" in get_example_code(name="no_such_example")
