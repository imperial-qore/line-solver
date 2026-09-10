#!/usr/bin/env python3
"""LINE Solver MCP Server

Exposes LINE's Python native queueing network solver as an MCP server,
allowing Claude (or any MCP client) to create, analyze, and solve
queueing models interactively during conversations.
"""

import sys
import os
import io
import json
import signal
import threading
import time
import traceback
import math

# Add python/ to path so line_solver is importable
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from mcp.server.fastmcp import FastMCP

# ---------------------------------------------------------------------------
# Richer server instructions (enhancement 10)
# ---------------------------------------------------------------------------
_INSTRUCTIONS = """\
LINE is an open-source queueing network solver for performance and capacity modeling.

## Available Tools

| Tool | Purpose |
|------|---------|
| `solve_line_model` | Execute arbitrary LINE Python code (sandbox-restricted) |
| `analyze_queue` | Quick M/M/k or M/M/k/N open-queue analysis |
| `analyze_closed_network` | Quick closed-network (think time + queue) analysis |
| `sweep_parameter` | Sweep one parameter over a range and collect metrics |
| `compare_solvers` | Run the same model with multiple solvers side-by-side |
| `visualize_model` | Generate a Mermaid diagram of a queueing network |
| `list_examples` | Browse available example models by category |
| `get_example_code` | Retrieve source code of a named example |
| `clear_session` | Remove a persistent code session |

## Example Taxonomy

- **Open queues**: `gallery_mm1`, `oqn_oneline`, `oqn_multiclass`
- **Closed networks**: `gallery_cqn`, `cqn_multiserver`, `cqn_repairmen`
- **Layered QN**: `lqn_basic`, `lqn_workflows`
- **Fork-join**: `fj_basic_open`, `fj_basic_nesting`
- **Caches**: `cache_replc_lru`, `cache_replc_routing`
- **Random environments**: `renv_threestages_repairmen`, `renv_node_breakdown`

## Solver Guide

| Solver | Best For |
|--------|----------|
| **MVA** | Product-form networks with exponential service (fast, exact) |
| **NC** | Closed networks needing exact normalizing-constant analysis |
| **CTMC** | Small state-space models; supports non-product-form features |
| **MAM** | Phase-type / MAP service distributions (matrix analytic) |
| **SSA** | Stochastic simulation; works for any supported model |
| **FLD** | Fluid / ODE approximation for transient and large models |

## Output Formats

All analysis tools accept `format="text"` (default, human-readable) or \
`format="json"` (structured, machine-readable with metadata + results keys).

## Code Execution

`solve_line_model` runs user code in a sandboxed namespace with \
`from line_solver import *`, numpy, and pandas pre-loaded.  Dangerous \
builtins (`__import__`, `exec`, `eval`, `open`, `compile`) are blocked.  \
Use `session_id` to persist variables across calls.

All code runs with 'from line_solver import *' pre-loaded.
"""

mcp = FastMCP("LINE Solver", instructions=_INSTRUCTIONS)

EXAMPLES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "examples")

# ---------------------------------------------------------------------------
# 1. Sandbox helpers
# ---------------------------------------------------------------------------
_SAFE_BUILTINS = {
    "abs": abs, "all": all, "any": any, "bin": bin, "bool": bool,
    "bytearray": bytearray, "bytes": bytes, "callable": callable,
    "chr": chr, "complex": complex, "dict": dict, "dir": dir,
    "divmod": divmod, "enumerate": enumerate, "filter": filter,
    "float": float, "format": format, "frozenset": frozenset,
    "getattr": getattr, "hasattr": hasattr, "hash": hash, "hex": hex,
    "id": id, "int": int, "isinstance": isinstance, "issubclass": issubclass,
    "iter": iter, "len": len, "list": list, "map": map, "max": max,
    "min": min, "next": next, "object": object, "oct": oct, "ord": ord,
    "pow": pow, "print": print, "property": property, "range": range,
    "repr": repr, "reversed": reversed, "round": round, "set": set,
    "setattr": setattr, "slice": slice, "sorted": sorted,
    "staticmethod": staticmethod, "str": str, "sum": sum, "super": super,
    "tuple": tuple, "type": type, "vars": vars, "zip": zip,
    "True": True, "False": False, "None": None,
    "ArithmeticError": ArithmeticError, "AssertionError": AssertionError,
    "AttributeError": AttributeError, "EOFError": EOFError,
    "Exception": Exception, "IndexError": IndexError,
    "KeyError": KeyError, "NameError": NameError,
    "NotImplementedError": NotImplementedError, "OSError": OSError,
    "RuntimeError": RuntimeError, "StopIteration": StopIteration,
    "TypeError": TypeError, "ValueError": ValueError,
    "ZeroDivisionError": ZeroDivisionError,
}

# Explicitly blocked names (for clarity)
_BLOCKED_BUILTINS = {
    "__import__", "exec", "eval", "compile", "open",
    "globals", "locals", "breakpoint", "exit", "quit",
    "input", "memoryview", "help",
}


def _make_exec_namespace(extra: dict | None = None) -> dict:
    """Create sandboxed execution namespace with line_solver pre-imported.

    Phase 1: import with full builtins so import machinery works.
    Phase 2: replace __builtins__ with restricted dict.
    """
    ns = {"__builtins__": __builtins__}
    if extra:
        ns.update(extra)

    # Phase 1 — imports need full builtins
    exec("from line_solver import *", ns)
    exec("import numpy as np", ns)
    exec("import pandas as pd", ns)
    exec("import time as _time", ns)

    # Phase 2 — lock down builtins
    ns["__builtins__"] = dict(_SAFE_BUILTINS)
    return ns


# ---------------------------------------------------------------------------
# 2. Execution timeout
# ---------------------------------------------------------------------------
class _TimeoutError(Exception):
    pass


def _exec_with_timeout(code: str, namespace: dict, timeout: int = 120):
    """Execute *code* in *namespace* with a wall-clock timeout (seconds)."""
    # Prefer SIGALRM on Linux (works only in main thread)
    if threading.current_thread() is threading.main_thread():
        def _handler(signum, frame):
            raise _TimeoutError(
                f"Execution timed out after {timeout}s. "
                "Consider simplifying the model or reducing state-space size."
            )
        old = signal.signal(signal.SIGALRM, _handler)
        signal.alarm(timeout)
        try:
            exec(code, namespace)
        finally:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old)
    else:
        # Fallback for non-main threads: run in a worker thread
        result = {"exc": None}

        def _worker():
            try:
                exec(code, namespace)
            except Exception as e:
                result["exc"] = e

        t = threading.Thread(target=_worker, daemon=True)
        t.start()
        t.join(timeout)
        if t.is_alive():
            raise _TimeoutError(
                f"Execution timed out after {timeout}s. "
                "Consider simplifying the model or reducing state-space size."
            )
        if result["exc"] is not None:
            raise result["exc"]


# ---------------------------------------------------------------------------
# 4. Structured JSON output helper
# ---------------------------------------------------------------------------
def _avg_table_to_json(avg_table, metadata: dict | None = None) -> str:
    """Convert an avg_table (IndexedTable or DataFrame) to JSON string."""
    import pandas as pd

    df = avg_table.data if hasattr(avg_table, "data") else avg_table
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame({"result": [str(avg_table)]})

    records = df.to_dict(orient="records")
    # Sanitise NaN / Inf for JSON
    for rec in records:
        for k, v in rec.items():
            if isinstance(v, float):
                if math.isnan(v):
                    rec[k] = None
                elif math.isinf(v):
                    rec[k] = "Inf" if v > 0 else "-Inf"

    payload = {"results": records}
    if metadata:
        payload["metadata"] = metadata
    return json.dumps(payload, indent=2, default=str)


# ---------------------------------------------------------------------------
# 9. Session state
# ---------------------------------------------------------------------------
_sessions: dict[str, dict] = {}
_SESSION_TTL = 3600  # 1 hour


def _get_or_create_session(session_id: str) -> dict:
    """Return (and lazily GC) a persistent namespace for *session_id*."""
    now = time.time()
    # Lazy cleanup of expired sessions
    expired = [k for k, v in _sessions.items() if now - v["last_access"] > _SESSION_TTL]
    for k in expired:
        del _sessions[k]

    if session_id not in _sessions:
        _sessions[session_id] = {
            "namespace": _make_exec_namespace(),
            "last_access": now,
        }
    else:
        _sessions[session_id]["last_access"] = now
    return _sessions[session_id]["namespace"]


# ---------------------------------------------------------------------------
# Tools
# ---------------------------------------------------------------------------

@mcp.tool()
def solve_line_model(code: str, format: str = "text", session_id: str = "") -> str:
    """Execute LINE Python code and return results.

    The code runs in a namespace with `from line_solver import *` already loaded,
    plus numpy as np and pandas as pd. Assign results to `avg_table` or print them
    to see output.

    Example:
        model = Network('M/M/1')
        source = Source(model, 'Source')
        queue = Queue(model, 'Queue', SchedStrategy.FCFS)
        sink = Sink(model, 'Sink')
        oclass = OpenClass(model, 'Class1')
        source.setArrival(oclass, Exp(1.0))
        queue.setService(oclass, Exp(2.0))
        model.link(Network.serial_routing([source, queue, sink]))
        avg_table = MVA(model).avg_table()
        print(avg_table)
    """
    # Build or retrieve namespace
    try:
        if session_id:
            exec_globals = _get_or_create_session(session_id)
        else:
            exec_globals = _make_exec_namespace()
    except Exception as e:
        return f"Error importing line_solver: {e}"

    # Capture stdout
    old_stdout = sys.stdout
    sys.stdout = captured = io.StringIO()

    try:
        _exec_with_timeout(code, exec_globals)
    except _TimeoutError as te:
        sys.stdout = old_stdout
        return str(te)
    except Exception:
        sys.stdout = old_stdout
        output = captured.getvalue()
        tb = traceback.format_exc()
        parts = []
        if output.strip():
            parts.append(output.strip())
        parts.append(f"Error:\n{tb}")
        return "\n".join(parts)
    finally:
        sys.stdout = old_stdout

    output = captured.getvalue()

    # Check for avg_table in the namespace
    parts = []
    if output.strip():
        parts.append(output.strip())

    avg = exec_globals.get("avg_table")
    if avg is not None:
        if format == "json":
            return _avg_table_to_json(avg)
        table_str = str(avg)
        if table_str not in output:
            parts.append(f"avg_table:\n{table_str}")

    if format == "json" and not avg:
        return json.dumps({"output": output.strip()}, indent=2, default=str)

    return "\n".join(parts) if parts else "(no output)"


@mcp.tool()
def list_examples() -> str:
    """List available LINE solver example models by category.

    Returns a categorized list of example files that can be retrieved
    with get_example_code.
    """
    if not os.path.isdir(EXAMPLES_DIR):
        return f"Examples directory not found: {EXAMPLES_DIR}"

    lines = ["Available LINE Solver Examples", "=" * 40, ""]

    for category in sorted(os.listdir(EXAMPLES_DIR)):
        cat_path = os.path.join(EXAMPLES_DIR, category)
        if not os.path.isdir(cat_path):
            continue

        examples = []
        for root, _dirs, files in os.walk(cat_path):
            for f in sorted(files):
                if f.endswith(".py") and not f.startswith("__"):
                    name = f[:-3]  # strip .py
                    rel = os.path.relpath(os.path.join(root, f), EXAMPLES_DIR)
                    examples.append((name, rel))

        if examples:
            lines.append(f"## {category}")
            for name, rel in examples:
                lines.append(f"  - {name}  ({rel})")
            lines.append("")

    return "\n".join(lines)


@mcp.tool()
def get_example_code(name: str) -> str:
    """Return the source code of a named example.

    Args:
        name: Example name without .py extension (e.g., "gallery_mm1", "oqn_oneline").
              Use list_examples to see available names.
    """
    target = name if name.endswith(".py") else name + ".py"

    for root, _dirs, files in os.walk(EXAMPLES_DIR):
        for f in files:
            if f == target:
                path = os.path.join(root, f)
                with open(path, "r") as fh:
                    return fh.read()

    return f"Example '{name}' not found. Use list_examples() to see available examples."


# ---------------------------------------------------------------------------
# 3. Input validation + analyze_queue
# ---------------------------------------------------------------------------

@mcp.tool()
def analyze_queue(
    arrival_rate: float,
    service_rate: float,
    servers: int = 1,
    queue_capacity: int = -1,
    scheduling: str = "FCFS",
    solver: str = "MVA",
    format: str = "text",
) -> str:
    """Quick single-queue analysis without writing code.

    Builds and solves an M/M/k model (or M/M/k/N with finite capacity).

    Args:
        arrival_rate: Job arrival rate (lambda).
        service_rate: Service rate per server (mu).
        servers: Number of servers (k). Default 1.
        queue_capacity: Total queue capacity including servers. -1 for infinite (default).
        scheduling: Scheduling strategy: FCFS, PS, LCFS, INF. Default FCFS.
        solver: Solver to use: MVA, CTMC, NC, SSA, FLD, MAM. Default MVA.
        format: Output format: "text" (default) or "json" (structured).

    Returns:
        Formatted performance metrics table.
    """
    # --- Input validation ---
    warnings: list[str] = []
    if arrival_rate <= 0:
        return "Error: arrival_rate must be > 0."
    if service_rate <= 0:
        return "Error: service_rate must be > 0."
    if servers < 1:
        return "Error: servers must be >= 1."

    rho = arrival_rate / (service_rate * servers)
    if queue_capacity < 0 and rho >= 1.0:
        return (
            f"Error: system is unstable (rho = {rho:.4f} >= 1). "
            "Reduce arrival_rate, increase service_rate/servers, "
            "or set a finite queue_capacity."
        )

    solver_upper = solver.upper()
    if solver_upper == "FLD":
        warnings.append("Note: FLD (fluid) is an approximation; results may differ from exact solvers.")
    if solver_upper == "CTMC" and queue_capacity < 0:
        warnings.append("Note: CTMC on infinite-capacity queue uses a truncated state space (cutoff heuristic).")

    try:
        from line_solver import (
            Network, Source, Queue, Sink, OpenClass,
            Exp, SchedStrategy,
        )
        from line_solver.solvers import (
            MVA as SolverMVA_cls, CTMC as SolverCTMC_cls,
            NC as SolverNC_cls, SSA as SolverSSA_cls,
            FLD as SolverFLD_cls, MAM as SolverMAM_cls,
        )
    except Exception as e:
        return f"Error importing line_solver: {e}"

    sched_map = {
        "FCFS": SchedStrategy.FCFS,
        "PS": SchedStrategy.PS,
        "LCFS": SchedStrategy.LCFS,
        "INF": SchedStrategy.INF,
    }
    sched = sched_map.get(scheduling.upper())
    if sched is None:
        return f"Unknown scheduling strategy '{scheduling}'. Use: {', '.join(sched_map)}"

    solver_map = {
        "MVA": SolverMVA_cls,
        "CTMC": SolverCTMC_cls,
        "NC": SolverNC_cls,
        "SSA": SolverSSA_cls,
        "FLD": SolverFLD_cls,
        "MAM": SolverMAM_cls,
    }
    solver_cls = solver_map.get(solver_upper)
    if solver_cls is None:
        return f"Unknown solver '{solver}'. Use: {', '.join(solver_map)}"

    try:
        notation = f"M/M/{servers}" if servers > 1 else "M/M/1"
        if queue_capacity > 0:
            notation += f"/{queue_capacity}"

        model = Network(notation)
        source = Source(model, "Source")
        queue = Queue(model, "Queue", sched)
        sink = Sink(model, "Sink")

        oclass = OpenClass(model, "Class1")
        source.setArrival(oclass, Exp(arrival_rate))
        queue.setService(oclass, Exp(service_rate))

        if servers > 1:
            queue.setNumberOfServers(servers)
        if queue_capacity > 0:
            queue.setCapacity(queue_capacity)

        model.link(Network.serial_routing([source, queue, sink]))

        if solver_upper == "CTMC":
            cutoff = queue_capacity if queue_capacity > 0 else max(50, int(10 * rho))
            s = solver_cls(model, cutoff=cutoff)
        else:
            s = solver_cls(model)

        avg_table = s.avg_table()

        metadata = {
            "model": notation,
            "arrival_rate": arrival_rate,
            "service_rate": service_rate,
            "servers": servers,
            "utilization": round(rho, 6),
            "solver": solver_upper,
        }

        if format == "json":
            if warnings:
                metadata["warnings"] = warnings
            return _avg_table_to_json(avg_table, metadata)

        header = []
        for w in warnings:
            header.append(f"WARNING: {w}")
        if warnings:
            header.append("")
        header += [
            f"Model: {notation}",
            f"  arrival_rate (lambda) = {arrival_rate}",
            f"  service_rate (mu)     = {service_rate}",
            f"  servers (k)           = {servers}",
            f"  utilization (rho)     = {rho:.6f}",
            f"  solver                = {solver_upper}",
            "",
        ]
        return "\n".join(header) + str(avg_table)

    except Exception:
        return f"Error solving model:\n{traceback.format_exc()}"


# ---------------------------------------------------------------------------
# 5. Closed-network quick analysis
# ---------------------------------------------------------------------------

@mcp.tool()
def analyze_closed_network(
    population: int,
    think_time: float,
    service_rate: float,
    servers: int = 1,
    scheduling: str = "FCFS",
    solver: str = "MVA",
    format: str = "text",
) -> str:
    """Quick closed-network analysis without writing code.

    Builds a Delay (think time) + Queue model with a fixed job population.

    Args:
        population: Number of jobs circulating in the network (>= 1).
        think_time: Mean think time at the Delay node (>= 0). 0 means no think time.
        service_rate: Service rate per server (mu) at the Queue.
        servers: Number of servers (k). Default 1.
        scheduling: Scheduling strategy: FCFS, PS, LCFS, INF. Default FCFS.
        solver: Solver to use: MVA, CTMC, NC, SSA, FLD, MAM. Default MVA.
        format: Output format: "text" (default) or "json" (structured).

    Returns:
        Formatted performance metrics table.
    """
    if population < 1:
        return "Error: population must be >= 1."
    if think_time < 0:
        return "Error: think_time must be >= 0."
    if service_rate <= 0:
        return "Error: service_rate must be > 0."
    if servers < 1:
        return "Error: servers must be >= 1."

    try:
        from line_solver import (
            Network, Delay, Queue, ClosedClass,
            Exp, SchedStrategy,
        )
        from line_solver.solvers import (
            MVA as SolverMVA_cls, CTMC as SolverCTMC_cls,
            NC as SolverNC_cls, SSA as SolverSSA_cls,
            FLD as SolverFLD_cls, MAM as SolverMAM_cls,
        )
    except Exception as e:
        return f"Error importing line_solver: {e}"

    sched_map = {
        "FCFS": SchedStrategy.FCFS,
        "PS": SchedStrategy.PS,
        "LCFS": SchedStrategy.LCFS,
        "INF": SchedStrategy.INF,
    }
    sched = sched_map.get(scheduling.upper())
    if sched is None:
        return f"Unknown scheduling strategy '{scheduling}'. Use: {', '.join(sched_map)}"

    solver_map = {
        "MVA": SolverMVA_cls,
        "CTMC": SolverCTMC_cls,
        "NC": SolverNC_cls,
        "SSA": SolverSSA_cls,
        "FLD": SolverFLD_cls,
        "MAM": SolverMAM_cls,
    }
    solver_upper = solver.upper()
    solver_cls = solver_map.get(solver_upper)
    if solver_cls is None:
        return f"Unknown solver '{solver}'. Use: {', '.join(solver_map)}"

    try:
        notation = f"Closed / N={population} / k={servers}"
        model = Network("ClosedModel")
        delay = Delay(model, "Think")
        queue = Queue(model, "Queue", sched)

        # think_time = 0 edge case: use very small value to avoid division by zero
        effective_think = think_time if think_time > 0 else 1e-9
        cclass = ClosedClass(model, "Class1", population, delay)
        delay.setService(cclass, Exp(1.0 / effective_think))
        queue.setService(cclass, Exp(service_rate))

        if servers > 1:
            queue.setNumberOfServers(servers)

        model.link(Network.serial_routing([delay, queue]))

        if solver_upper == "CTMC":
            s = solver_cls(model, cutoff=population)
        else:
            s = solver_cls(model)

        avg_table = s.avg_table()

        metadata = {
            "model": notation,
            "population": population,
            "think_time": think_time,
            "service_rate": service_rate,
            "servers": servers,
            "solver": solver_upper,
        }

        if format == "json":
            return _avg_table_to_json(avg_table, metadata)

        header = [
            f"Model: {notation}",
            f"  population (N)        = {population}",
            f"  think_time (Z)        = {think_time}",
            f"  service_rate (mu)     = {service_rate}",
            f"  servers (k)           = {servers}",
            f"  solver                = {solver_upper}",
            "",
        ]
        return "\n".join(header) + str(avg_table)

    except Exception:
        return f"Error solving model:\n{traceback.format_exc()}"


# ---------------------------------------------------------------------------
# 6. Parameter sweep
# ---------------------------------------------------------------------------

def _parse_values(values_str: str) -> list[float]:
    """Parse comma-separated or start:stop:step range notation."""
    values_str = values_str.strip()
    if ":" in values_str:
        parts = values_str.split(":")
        if len(parts) != 3:
            raise ValueError("Range notation must be start:stop:step (e.g. '0.1:2.0:0.1')")
        start, stop, step = float(parts[0]), float(parts[1]), float(parts[2])
        if step <= 0:
            raise ValueError("Step must be > 0")
        vals = []
        v = start
        while v <= stop + step * 1e-9:
            vals.append(round(v, 10))
            v += step
        return vals
    else:
        return [float(x.strip()) for x in values_str.split(",")]


@mcp.tool()
def sweep_parameter(
    param: str,
    values: str,
    arrival_rate: float = 1.0,
    service_rate: float = 2.0,
    servers: int = 1,
    queue_capacity: int = -1,
    scheduling: str = "FCFS",
    solver: str = "MVA",
    model_type: str = "open",
    population: int = 10,
    think_time: float = 5.0,
) -> str:
    """Sweep one parameter over a range and collect performance metrics.

    Args:
        param: Parameter to sweep. One of: arrival_rate, service_rate, servers,
               population, think_time.
        values: Values to sweep — comma-separated ("0.5,1.0,1.5") or
                range notation ("0.1:2.0:0.1" meaning start:stop:step).
        arrival_rate: Base arrival rate (for open models). Default 1.0.
        service_rate: Base service rate. Default 2.0.
        servers: Base number of servers. Default 1.
        queue_capacity: Queue capacity (-1 = infinite). Default -1.
        scheduling: Scheduling strategy. Default FCFS.
        solver: Solver to use. Default MVA.
        model_type: "open" or "closed". Default "open".
        population: Base population (for closed models). Default 10.
        think_time: Base think time (for closed models). Default 5.0.

    Returns:
        Table of metrics for each swept value.
    """
    import pandas as pd

    try:
        vals = _parse_values(values)
    except ValueError as e:
        return f"Error parsing values: {e}"

    if len(vals) > 100:
        return "Error: maximum 100 sweep points allowed."
    if len(vals) == 0:
        return "Error: no values to sweep."

    valid_params = {"arrival_rate", "service_rate", "servers", "population", "think_time"}
    if param not in valid_params:
        return f"Unknown parameter '{param}'. Use one of: {', '.join(sorted(valid_params))}"

    if model_type not in ("open", "closed"):
        return f"Unknown model_type '{model_type}'. Use 'open' or 'closed'."

    rows = []
    errors = []
    for v in vals:
        kwargs = {
            "service_rate": service_rate,
            "servers": servers,
            "scheduling": scheduling,
            "solver": solver,
            "format": "json",
        }
        if model_type == "open":
            kwargs["arrival_rate"] = arrival_rate
            kwargs["queue_capacity"] = queue_capacity
        else:
            kwargs["population"] = population
            kwargs["think_time"] = think_time

        # Override the swept parameter
        if param == "servers":
            kwargs["servers"] = int(v)
        elif param == "population":
            kwargs["population"] = int(v)
        else:
            kwargs[param] = v

        try:
            if model_type == "open":
                raw = analyze_queue(**kwargs)
            else:
                raw = analyze_closed_network(**kwargs)

            data = json.loads(raw)
            if "results" in data:
                for rec in data["results"]:
                    rec[param] = v
                    rows.append(rec)
            else:
                errors.append(f"{param}={v}: {data.get('output', raw)}")
        except Exception as exc:
            errors.append(f"{param}={v}: {exc}")

    if not rows:
        return "No successful runs.\n" + "\n".join(errors)

    df = pd.DataFrame(rows)
    # Move swept param to front
    if param in df.columns:
        cols = [param] + [c for c in df.columns if c != param]
        df = df[cols]

    result = df.to_string(index=False)
    if errors:
        result += "\n\nErrors:\n" + "\n".join(errors)
    return result


# ---------------------------------------------------------------------------
# 7. Model visualization (Mermaid)
# ---------------------------------------------------------------------------

# NodeType int → Mermaid shape template
_NODE_SHAPES = {
    0: ("([{label}])", "source"),      # SOURCE — stadium
    1: ("([{label}])", "sink"),        # SINK — stadium
    2: ("[{label}]", "queue"),         # QUEUE — rectangle
    3: ("({label})", "delay"),         # DELAY — rounded
    4: ("{{{{{label}}}}}", "fork"),    # FORK — diamond (hexagon)
    5: ("{{{{{label}}}}}", "join"),    # JOIN — diamond (hexagon)
    6: ("[({label})]", "cache"),       # CACHE — cylinder
    7: ("({label})", "router"),        # ROUTER — rounded
    8: ("[{label}]", "classswitch"),   # CLASSSWITCH — rectangle
}

_NODE_COLORS = {
    "source": "#4CAF50",
    "sink": "#F44336",
    "queue": "#2196F3",
    "delay": "#FF9800",
    "fork": "#9C27B0",
    "join": "#9C27B0",
    "cache": "#00BCD4",
    "router": "#795548",
    "classswitch": "#607D8B",
}


def _generate_mermaid(sn) -> str:
    """Generate a Mermaid graph definition from a NetworkStruct."""
    import numpy as np

    lines = ["graph LR"]
    styles_used: set[str] = set()
    nnodes = sn.nnodes

    # Nodes
    for i in range(nnodes):
        ntype = int(sn.nodetype[i]) if i < len(sn.nodetype) else 2
        shape_tpl, style_class = _NODE_SHAPES.get(ntype, ("[{label}]", "queue"))

        name = sn.nodenames[i] if i < len(sn.nodenames) else f"Node{i}"
        label = name

        # Add server count for station nodes
        n2s = sn.nodeToStation
        if n2s is not None and i < len(n2s):
            ist = int(n2s[i])
            if ist >= 0 and sn.nservers is not None and ist < len(sn.nservers):
                ns_raw = float(sn.nservers[ist])
                if np.isinf(ns_raw):
                    label += " (IS)"
                elif ns_raw > 1:
                    label += f" (k={int(ns_raw)})"

        node_def = shape_tpl.format(label=label)
        lines.append(f"    N{i}{node_def}")
        styles_used.add(style_class)

    # Edges from connmatrix
    if sn.connmatrix is not None:
        cm = np.array(sn.connmatrix)
        for i in range(nnodes):
            for j in range(nnodes):
                if cm[i, j]:
                    lines.append(f"    N{i} --> N{j}")

    # Styles
    for sc in sorted(styles_used):
        color = _NODE_COLORS.get(sc, "#9E9E9E")
        # Collect node IDs for this style
        ids = []
        for i in range(nnodes):
            ntype = int(sn.nodetype[i]) if i < len(sn.nodetype) else 2
            _, sc2 = _NODE_SHAPES.get(ntype, ("[{label}]", "queue"))
            if sc2 == sc:
                ids.append(f"N{i}")
        if ids:
            lines.append(f"    style {','.join(ids)} fill:{color},color:#fff,stroke:{color}")

    # Legend with class info
    class_info = []
    for r in range(sn.nclasses):
        cname = sn.classnames[r] if r < len(sn.classnames) else f"Class{r}"
        njobs_val = sn.njobs.flat[r] if r < sn.njobs.size else "?"
        if np.isinf(float(njobs_val)):
            class_info.append(f"{cname} (open)")
        else:
            class_info.append(f"{cname} (N={int(njobs_val)})")

    if class_info:
        lines.append("")
        lines.append(f"    %% Classes: {', '.join(class_info)}")

    return "\n".join(lines)


@mcp.tool()
def visualize_model(code: str) -> str:
    """Generate a Mermaid diagram of a queueing network.

    Executes the provided LINE code, finds the Network object, and produces
    a Mermaid graph definition showing nodes, connections, and class info.

    Args:
        code: LINE Python code that creates a Network and calls model.link().

    Returns:
        Mermaid diagram source (paste into any Mermaid renderer).
    """
    try:
        ns = _make_exec_namespace()
    except Exception as e:
        return f"Error importing line_solver: {e}"

    old_stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        _exec_with_timeout(code, ns)
    except Exception:
        sys.stdout = old_stdout
        return f"Error executing code:\n{traceback.format_exc()}"
    finally:
        sys.stdout = old_stdout

    # Find Network object in namespace
    from line_solver import Network
    model = None
    for v in ns.values():
        if isinstance(v, Network):
            model = v
            break

    if model is None:
        return "Error: no Network object found in code. Make sure your code creates a Network."

    try:
        sn = model.get_struct()
        return _generate_mermaid(sn)
    except Exception:
        return f"Error generating diagram:\n{traceback.format_exc()}"


# ---------------------------------------------------------------------------
# 8. Multi-solver comparison
# ---------------------------------------------------------------------------

@mcp.tool()
def compare_solvers(code: str, solvers: str = "") -> str:
    """Run the same model with multiple solvers and compare results.

    Executes the provided LINE code to build a Network, then solves it
    with each requested solver. Returns a merged table with a Solver column.

    Args:
        code: LINE Python code that creates a Network and calls model.link().
        solvers: Comma-separated solver names (default: "MVA,NC,CTMC,SSA,FLD,MAM").

    Returns:
        Combined performance metrics table from all solvers.
    """
    import pandas as pd

    if not solvers.strip():
        solver_list = ["MVA", "NC", "CTMC", "SSA", "FLD", "MAM"]
    else:
        solver_list = [s.strip().upper() for s in solvers.split(",") if s.strip()]

    try:
        ns = _make_exec_namespace()
    except Exception as e:
        return f"Error importing line_solver: {e}"

    old_stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        _exec_with_timeout(code, ns)
    except Exception:
        sys.stdout = old_stdout
        return f"Error executing code:\n{traceback.format_exc()}"
    finally:
        sys.stdout = old_stdout

    # Find Network object
    from line_solver import Network
    model = None
    for v in ns.values():
        if isinstance(v, Network):
            model = v
            break

    if model is None:
        return "Error: no Network object found in code. Make sure your code creates a Network."

    from line_solver.solvers import (
        MVA as SolverMVA_cls, CTMC as SolverCTMC_cls,
        NC as SolverNC_cls, SSA as SolverSSA_cls,
        FLD as SolverFLD_cls, MAM as SolverMAM_cls,
    )
    import numpy as np

    solver_map = {
        "MVA": SolverMVA_cls,
        "CTMC": SolverCTMC_cls,
        "NC": SolverNC_cls,
        "SSA": SolverSSA_cls,
        "FLD": SolverFLD_cls,
        "MAM": SolverMAM_cls,
    }

    all_dfs = []
    errors = []

    # Determine if open model (for CTMC cutoff)
    sn = model.get_struct()
    is_open = bool(np.any(np.isinf(sn.njobs)))

    for solver_name in solver_list:
        solver_cls = solver_map.get(solver_name)
        if solver_cls is None:
            errors.append(f"{solver_name}: unknown solver")
            continue

        try:
            model.reset()
            if solver_name == "CTMC" and is_open:
                s = solver_cls(model, cutoff=50)
            else:
                s = solver_cls(model)

            avg_table = s.avg_table()
            df = avg_table.data if hasattr(avg_table, "data") else pd.DataFrame({"result": [str(avg_table)]})
            df = df.copy()
            df.insert(0, "Solver", solver_name)
            all_dfs.append(df)
        except Exception:
            errors.append(f"{solver_name}: {traceback.format_exc().splitlines()[-1]}")

    if not all_dfs:
        return "No solvers succeeded.\n" + "\n".join(errors)

    merged = pd.concat(all_dfs, ignore_index=True)
    result = merged.to_string(index=False)
    if errors:
        result += "\n\nErrors:\n" + "\n".join(errors)
    return result


# ---------------------------------------------------------------------------
# 9. clear_session tool
# ---------------------------------------------------------------------------

@mcp.tool()
def clear_session(session_id: str) -> str:
    """Remove a persistent code session.

    Args:
        session_id: The session ID to clear.

    Returns:
        Confirmation message.
    """
    if session_id in _sessions:
        del _sessions[session_id]
        return f"Session '{session_id}' cleared."
    return f"Session '{session_id}' not found."


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    mcp.run()
