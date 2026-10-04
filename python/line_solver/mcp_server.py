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
import tempfile
import threading
import time
import traceback
import math
import contextlib
import functools
import itertools

from mcp.server.mcpserver import MCPServer

# ---------------------------------------------------------------------------
# Richer server instructions (enhancement 10)
# ---------------------------------------------------------------------------
_INSTRUCTIONS = """\
LINE is an open-source queueing network solver for performance and capacity modeling.

## Available Tools

| Tool | Purpose |
|------|---------|
| `solve_model_file` | Solve a model DOCUMENT: the whole `line-cli` command line |
| `find_solver` | Which solvers and methods can analyze a given model, and why not |
| `list_solvers` | The vocabulary: solvers, analyses, input and output formats |
| `build_network` | Build ANY network from a description: open, closed or mixed, many stations, many classes, no code |
| `build_topology` | Build a multi-station network from a demand matrix, with no routing |
| `build_from_gallery` | Build one of the ready-made gallery models |
| `list_gallery` | The gallery catalogue, with the arguments each model takes |
| `convert_model` | Convert a model to LINE json, JMT, PNML, .lqnx or Python |
| `sweep_model` | Sweep any field of any model description |
| `solve_line_model` | Execute arbitrary LINE Python code (sandbox-restricted) |
| `analyze_queue` | Build and solve ONE station from arguments: open, closed or mixed |
| `analyze_network` | Build and solve a multi-station network from a demand matrix |
| `sweep_parameter` | Sweep one parameter of a single station and collect metrics |
| `compare_solvers` | Run the same model with multiple solvers side-by-side |
| `visualize_model` | Generate a Mermaid diagram of a queueing network |
| `list_examples` | Browse available example models by category |
| `get_example_code` | Retrieve source code of a named example |
| `check_install` | Environment check: which optional backends are reachable |
| `get_version` | LINE version, interpreter and examples directory |
| `clear_session` | Remove a persistent code session |

## Which Tool

A model that already EXISTS as a file or as text is `solve_model_file`: it is
the `line-cli` surface, so it reads `json`, JMT `jsimg`, layered `lqnx` and
`pnml`, takes any solver or `solver.method`, and serves every analysis LINE
publishes -- transient, response-time CDF and percentiles, state probabilities,
simulated samples, rewards, cache hit/miss, normalizing constants, busy
periods, sensitivities.

A model to be BUILT needs no code. `build_network` takes a description of
stations, classes, service and routing and covers open, closed AND mixed
networks with any number of stations and classes; `build_topology` takes a
demand matrix and supplies the routing itself; `build_from_gallery` returns a
ready-made model (`list_gallery` names them all). Each answers with a
`model_id`, and that id then goes to `solve_model_file` for ANY analysis above,
or to `find_solver`, `visualize_model`, `compare_solvers` or `convert_model`.
Reach for `solve_line_model` only when the model needs Python the builders do
not express, such as a cache, a Petri net or a layered model built from
scratch.

A model wanted ONLY for its numbers skips the build: `analyze_queue` takes the
parameters of one station and `analyze_network` the demand matrix of several,
build the model and solve it in the same call. They are not a reduced world --
any distribution, any scheduling discipline, any number of classes, open,
closed or mixed, and the same solver, method and analysis vocabulary
`solve_model_file` takes. Build the model instead when it is to be kept and
then visualized, converted, or solved several ways.

Ask `find_solver` when it is not obvious which solver fits, and `list_solvers`
for the exact method names either one accepts.

Distributions in a built model are written `name(args)`: `exp(2.0)`,
`erlang(3.0, 2)`, `hyperexp(0.5, 1.0, 4.0)`, `det(0.25)`, `cox2(1.0, 2.0, 0.5)`,
`gamma(2.0, 0.5)`. Every distribution class the library publishes is accepted.

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
| **LDES** | Discrete-event simulation: the widest feature coverage |
| **AUTO** | The default: picks a solver that can run the model as described |

## Output Formats

Most tools accept `format="text"` (default, human-readable) or `format="json"`. \
The json reply carries `results` (and usually `metadata`), or `error` when the \
call failed. `solve_model_file` spells this argument `output` instead, and \
takes `readable`, `json`, `csv`, `pickle` and `mat`; its json reply is keyed by \
analysis name. `list_examples`, `get_example_code`, `visualize_model`, \
`compare_solvers` and `clear_session` answer in text only.

## Code Execution

`solve_line_model` runs user code in a sandboxed namespace with \
`from line_solver import *`, numpy, and pandas pre-loaded.  Dangerous \
builtins (`__import__`, `exec`, `eval`, `open`, `compile`) are blocked.  \
Use `session_id` to persist variables across calls.

All code runs with 'from line_solver import *' pre-loaded.
"""

mcp = MCPServer("LINE Solver", instructions=_INSTRUCTIONS)


def _resolve_examples_dir() -> str:
    """Directory holding the example models, for both supported layouts.

    An installed wheel places them inside the package, at
    ``line_solver/examples``; a source checkout keeps them one level up, at
    ``python/examples``, next to the package rather than inside it.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    installed = os.path.join(here, "examples")
    checkout = os.path.join(os.path.dirname(here), "examples")
    return checkout if not os.path.isdir(installed) and os.path.isdir(checkout) else installed


EXAMPLES_DIR = _resolve_examples_dir()


@contextlib.contextmanager
def _quiet_stdout():
    """Keep anything printed by a tool away from the JSON-RPC channel.

    MCP speaks JSON-RPC over stdout, so a solver banner written there is read
    by the client as a malformed message. LINE's analyzers print progress on
    completion, so every tool runs with stdout swapped for a throwaway buffer.
    The tools that deliberately collect printed output swap in a buffer of
    their own inside this one, which nests without interference.
    """
    saved = sys.stdout
    sys.stdout = io.StringIO()
    try:
        yield
    finally:
        sys.stdout = saved


def _tool(fn):
    """Register ``fn`` as an MCP tool, with its stdout held off the wire."""

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        with _quiet_stdout():
            return fn(*args, **kwargs)

    return mcp.tool()(wrapper)

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


def _error(message: str, format: str = "text") -> str:
    """Render a tool-level error in the caller's requested format.

    A json caller must get json back even when the run never reached a result:
    sweep_parameter reads these replies with json.loads, and a plain-text error
    surfaced there as a parse failure rather than as the reason it failed.
    """
    if format == "json":
        return json.dumps({"error": message}, indent=2, default=str)
    return message


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
# 10. Built models
# ---------------------------------------------------------------------------
# A model built in conversation has to outlive the call that built it, or every
# later analysis would mean rebuilding it. Keyed like a session and GC'd on the
# same clock, so a long-lived server does not accumulate models forever.
_models: dict[str, dict] = {}
_MODEL_TTL = 3600  # 1 hour
_model_seq = itertools.count(1)


def _register_model(model, origin: str) -> str:
    """Store *model* under a fresh id and return the id.

    Args:
        model: Any model object the solve chain accepts.
        origin: Short human-readable note on where it came from, echoed back to
            the caller so a model id is not an opaque handle.
    """
    now = time.time()
    for k in [k for k, v in _models.items() if now - v["last_access"] > _MODEL_TTL]:
        del _models[k]

    model_id = "m%d" % next(_model_seq)
    _models[model_id] = {"model": model, "last_access": now, "origin": origin}
    return model_id


def _resolve_model(model_id: str):
    """Return the model stored under *model_id*, refreshing its TTL.

    Raises:
        ValueError: when the id is unknown, which after an hour is what an
            expired model looks like; the message says so rather than leaving
            the caller to guess at a typo.
    """
    entry = _models.get(model_id)
    if entry is None:
        known = ", ".join(sorted(_models)) or "none"
        raise ValueError(
            "unknown model_id '%s' (models held: %s). A model id is returned by "
            "build_network, build_topology or build_from_gallery, and is dropped "
            "after an hour of disuse." % (model_id, known))
    entry["last_access"] = time.time()
    return entry["model"]


def _model_summary(model, model_id: str, origin: str) -> dict:
    """A compact description of a built model, for the builders to return."""
    summary: dict = {"model_id": model_id, "origin": origin,
                     "kind": type(model).__name__}
    try:
        summary["name"] = model.getName()
    except Exception:
        pass
    try:
        import numpy as np

        sn = model.get_struct()
        summary["stations"] = int(sn.nstations)
        summary["nodes"] = int(sn.nnodes)
        summary["classes"] = int(sn.nclasses)
        summary["chains"] = int(sn.nchains)
        njobs = np.array(sn.njobs).flatten()
        summary["open_classes"] = int(np.sum(np.isinf(njobs)))
        summary["closed_classes"] = int(np.sum(~np.isinf(njobs)))
        summary["population"] = [None if math.isinf(float(n)) else float(n)
                                 for n in njobs]
    except Exception:
        # A LayeredNetwork or an Environment has no NetworkStruct of this
        # shape; the identity above is still worth returning.
        pass
    return summary


def _render_built(model, model_id: str, origin: str, format: str) -> str:
    """Render a builder's reply in the requested format."""
    summary = _model_summary(model, model_id, origin)
    if format == "json":
        return json.dumps({"results": [summary], "metadata": {"model_id": model_id}},
                          indent=2, default=str)
    lines = ["Built %s  (model_id: %s)" % (summary.get("kind", "model"), model_id),
             "  origin: %s" % origin]
    for key in ("name", "nodes", "stations", "classes", "open_classes",
                "closed_classes", "chains"):
        if key in summary:
            lines.append("  %-15s %s" % (key, summary[key]))
    lines.append("")
    lines.append("Solve it with solve_model_file(model_id=\"%s\", analysis=...)." % model_id)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Tools
# ---------------------------------------------------------------------------

@_tool
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
        return _error(f"Error importing line_solver: {e}", format)

    # Capture stdout
    old_stdout = sys.stdout
    sys.stdout = captured = io.StringIO()

    try:
        _exec_with_timeout(code, exec_globals)
    except _TimeoutError as te:
        sys.stdout = old_stdout
        return _error(str(te), format)
    except Exception:
        sys.stdout = old_stdout
        output = captured.getvalue()
        tb = traceback.format_exc()
        parts = []
        if output.strip():
            parts.append(output.strip())
        parts.append(f"Error:\n{tb}")
        return _error("\n".join(parts), format)
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


@_tool
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


@_tool
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

def _distribution_classes() -> dict:
    """Every distribution the library exports, keyed by lowercased name.

    Read off the package rather than listed here, so a distribution added to
    `line_solver.distributions` is accepted by the builder with no edit in this
    file. The abstract bases are dropped; everything else is constructible.
    """
    import line_solver as ls
    from .distributions.base import Distribution

    abstract = {"Distribution", "ContinuousDistribution", "DiscreteDistribution",
                "Markovian"}
    out = {}
    for name in dir(ls):
        obj = getattr(ls, name, None)
        if (isinstance(obj, type) and issubclass(obj, Distribution)
                and name not in abstract):
            out[name.lower()] = obj
    return out


def _parse_distribution(spec, what: str):
    """Resolve a distribution written as `name(args)`.

    Examples: exp(2.0), erlang(3.0, 2), hyperexp(0.5, 1.0, 4.0), det(0.25),
    cox2(1.0, 2.0, 0.5), gamma(2.0, 0.5), replayer('trace.txt'), immediate().

    The name is matched against the library's own distribution classes and the
    arguments are read with ast.literal_eval, so numbers, strings and lists are
    all accepted and nothing is executed.
    """
    import ast

    if isinstance(spec, (int, float)):
        raise ValueError(
            "%s: write the distribution, not a bare number -- a rate and a mean "
            "are not the same thing. Use e.g. \"exp(%s)\"." % (what, spec))
    text = str(spec).strip()
    if not text:
        raise ValueError("%s: no distribution given" % what)
    if "(" not in text or not text.endswith(")"):
        raise ValueError(
            "%s: '%s' is not a distribution; write name(args), e.g. \"exp(2.0)\", "
            "\"erlang(3.0, 2)\", \"det(0.25)\"." % (what, text))

    head, _, rest = text.partition("(")
    known = _distribution_classes()
    cls = known.get(head.strip().lower())
    if cls is None:
        raise ValueError("%s: unknown distribution '%s'. Known: %s"
                         % (what, head.strip(), ", ".join(sorted(known))))

    argtext = rest[:-1].strip()
    try:
        args = ast.literal_eval("(" + argtext + ",)") if argtext else ()
    except (ValueError, SyntaxError) as e:
        raise ValueError("%s: could not read the arguments of '%s': %s"
                         % (what, text, e))
    try:
        return cls(*args)
    except Exception as e:
        import inspect

        try:
            sig = str(inspect.signature(cls.__init__)).replace("self, ", "")
        except (TypeError, ValueError):
            sig = "(...)"
        raise ValueError("%s: %s%s does not take %s: %s"
                         % (what, cls.__name__, sig, list(args), e))


def _model_from_any_source(model_id: str = "", file: str = "", content: str = "",
                           code: str = "", input_format: str = "", verbose: bool = False):
    """Resolve a model from whichever of the four sources was given.

    One place, so every tool that takes a model agrees on what the sources are
    and what it means to give two of them.
    """
    from . import cli as line_cli

    given = [n for n, v in (("model_id", model_id), ("file", file),
                            ("content", content), ("code", code)) if str(v).strip()]
    if not given:
        raise ValueError("no model given: pass `model_id` for a model built here, "
                         "`file` or `content` for a document, or `code` that builds one")
    if len(given) > 1:
        raise ValueError("give the model ONE way, not %d: %s" % (len(given), ", ".join(given)))

    if model_id.strip():
        return _resolve_model(model_id.strip())
    if code.strip():
        return _model_from_code(code)
    return _load_cli_model(file, content,
                           _resolve_input_format(file, content, input_format), verbose)


@_tool
def convert_model(
    target: str,
    model_id: str = "",
    file: str = "",
    content: str = "",
    code: str = "",
    input_format: str = "",
    output_file: str = "",
    name: str = "model",
    format: str = "text",
) -> str:
    """Convert a model between the formats LINE can read and write.

    Takes a model from any source and renders it as another format: the LINE
    json interchange document, a JMT .jsimg, a PNML net, a layered .lqnx, or
    the equivalent Python script. This is how a model built in conversation
    becomes something you can keep, edit, version or open in another tool.

    Args:
        target: json (LINE interchange), jsimg (JMT), pnml (Petri net),
            lqnx (layered), or python (the equivalent script).
        model_id: A model built here by build_network, build_topology or
            build_from_gallery.
        file: Path to a model document to convert instead.
        content: The document itself, when there is no file.
        code: LINE Python code building the model to convert.
        input_format: Format of file/content; detected when not stated.
        output_file: Write the result here instead of returning it inline.
        name: Variable name for the generated model, target="python" only.
        format: text or json. json wraps the document in {"results": [...]}.

    Returns:
        The converted document, or a note naming the file it was written to.
    """
    import tempfile

    try:
        model = _model_from_any_source(model_id, file, content, code, input_format)

        tgt = target.strip().lower()
        targets = ("json", "jsimg", "pnml", "lqnx", "python")
        if tgt not in targets:
            raise ValueError("unknown target '%s'; use one of: %s"
                             % (target, ", ".join(targets)))

        def _via_file(suffix, write):
            """For the writers that only emit to a path."""
            if output_file:
                write(output_file)
                return None
            handle = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
            handle.close()
            try:
                write(handle.name)
                with open(handle.name, "r") as fh:
                    return fh.read()
            finally:
                try:
                    os.unlink(handle.name)
                except OSError:
                    pass

        if tgt == "json":
            from .io.linemodel_io import model_to_dict

            document = json.dumps(model_to_dict(model), indent=2, default=str)
            if output_file:
                with open(output_file, "w") as fh:
                    fh.write(document)
                document = None
        elif tgt == "python":
            from .api.io.code_gen import qn2python

            document = qn2python(model, name)
            if output_file:
                with open(output_file, "w") as fh:
                    fh.write(document)
                document = None
        elif tgt == "jsimg":
            from .api.io.jmt_io import qn2jsimg

            # qn2jsimg answers with the PATH it wrote, not the document.
            written = qn2jsimg(model, output_file or None)
            if output_file:
                document = None
            else:
                with open(written, "r") as fh:
                    document = fh.read()
                try:
                    os.unlink(written)
                except OSError:
                    pass
        elif tgt == "pnml":
            from .io.pnml_io import save_pnml

            document = _via_file(".pnml", lambda path: save_pnml(model, path))
        else:
            from .layered import LayeredNetwork

            if not isinstance(model, LayeredNetwork):
                raise ValueError("target 'lqnx' needs a LayeredNetwork; this model "
                                 "is a %s" % type(model).__name__)
            document = _via_file(".lqnx", lambda path: model.writeXML(path))
    except Exception as e:
        return _error("Error: %s" % e, format)

    if document is None:
        message = "Wrote %s as %s to %s" % (type(model).__name__, tgt, output_file)
        if format == "json":
            return json.dumps({"results": [{"target": tgt, "output_file": output_file}],
                               "metadata": {"written": True}}, indent=2)
        return message

    if format == "json":
        return json.dumps({"results": [{"target": tgt, "document": document}],
                           "metadata": {"written": False}}, indent=2, default=str)
    return document


@_tool
def sweep_model(
    param: str,
    values: str,
    builder: str = "network",
    args: str = "",
    solver: str = "auto",
    analysis: str = "avg",
    format: str = "text",
) -> str:
    """Sweep ANY field of a model description and collect the results.

    The general sweep: where sweep_parameter varies one of five preset knobs of
    the two preset shapes, this varies any field of any builder's arguments --
    a population, a service rate, a server count, one cell of a demand matrix --
    and reports the chosen analysis at each point.

    Args:
        param: Dotted path into `args` naming the field to vary. List elements
            are indexed by number: "classes.1.population", "service.Q1.web",
            "demands.0.1", "params.M".
        values: Comma-separated numbers, start:stop:step, or a JSON list when
            the field is not numeric, e.g. '["exp(2.0)", "exp(4.0)"]' to sweep
            a service distribution. A field that currently holds an integer is
            swept with integers, so a population stays a whole number.
        builder: network (build_network), topology (build_topology), or
            gallery (build_from_gallery).
        args: JSON object of that builder's arguments, with the nested parts as
            REAL json rather than strings. For builder="network":
            '{"stations": [...], "classes": [...], "service": {...},
              "arrival": {...}, "routing": "serial"}'.
        solver: Solver for every point.
        analysis: Analysis for every point; the rows it returns are collected.
        format: text or json.

    Returns:
        One row per point per station-class, with the swept value in the first
        column, plus the reason for any point that failed.
    """
    import copy

    import pandas as pd

    try:
        spec = json.loads(args) if args.strip() else {}
        if not isinstance(spec, dict):
            raise ValueError("`args` is a JSON object of the builder's arguments")

        if values.strip().startswith("["):
            # A JSON list carries whatever the field holds -- a distribution
            # string, a boolean, a nested list -- where the numeric grammar
            # below cannot.
            points = json.loads(values)
            if not isinstance(points, list):
                raise ValueError("`values` as JSON must be a list")
        else:
            try:
                points = _parse_values(values)
            except ValueError as e:
                raise ValueError("could not read `values`: %s" % e)
        if not points:
            raise ValueError("no values to sweep")
        if len(points) > 100:
            raise ValueError("maximum 100 sweep points allowed, got %d" % len(points))

        builders = {"network": build_network, "topology": build_topology,
                    "gallery": build_from_gallery}
        build = builders.get(builder.strip().lower())
        if build is None:
            raise ValueError("unknown builder '%s'; use one of: %s"
                             % (builder, ", ".join(sorted(builders))))

        path = [p for p in str(param).split(".") if p]
        if not path:
            raise ValueError("`param` is a dotted path into `args`, "
                             'e.g. "classes.1.population"')

        def _walk(container, keys):
            """The container holding the last key, and that key."""
            cursor = container
            for step in keys[:-1]:
                cursor = cursor[int(step)] if isinstance(cursor, list) else cursor[step]
            last = keys[-1]
            return cursor, (int(last) if isinstance(cursor, list) else last)

        def _place(container, keys, value):
            cursor, last = _walk(container, keys)
            cursor[last] = value

        # A field holding an int must keep holding an int: _parse_values answers
        # in floats, and a population or a station count of 2.0 is a TypeError
        # inside the constructor rather than a sweep point.
        try:
            cursor, last = _walk(spec, path)
            current = cursor[last]
        except (KeyError, IndexError, TypeError, ValueError) as e:
            raise ValueError("`param` path '%s' does not name a field of `args`: %s"
                             % (param, e))
        if isinstance(current, int) and not isinstance(current, bool):
            points = [int(v) if isinstance(v, float) and v.is_integer() else v
                      for v in points]
    except Exception as e:
        return _error("Error: %s" % e, format)

    rows, errors = [], []
    for value in points:
        try:
            here = copy.deepcopy(spec)
            _place(here, path, value)
            # The builders take their nested arguments as json strings; rebuild
            # that spelling here so the sweep's args can stay real json.
            call = {k: (v if isinstance(v, str) else json.dumps(v))
                    for k, v in here.items()}
            built = build(format="json", **call)
            payload = json.loads(built)
            if "error" in payload:
                raise ValueError(payload["error"])
            mid = payload["metadata"]["model_id"]

            solved = solve_model_file(model_id=mid, solver=solver,
                                      analysis=analysis, output="json")
            data = json.loads(solved)
            if "error" in data:
                raise ValueError(data["error"])
            for _, block in data.items():
                if isinstance(block, list):
                    for record in block:
                        rows.append(dict({param: value}, **record))
        except Exception as e:
            errors.append("%s=%s: %s" % (param, value, e))

    if not rows:
        return _error("No successful runs.\n" + "\n".join(errors), format)

    frame = pd.DataFrame(rows)
    if format == "json":
        return json.dumps({"results": frame.to_dict(orient="records"),
                           "metadata": {"param": param, "points": len(points),
                                        "errors": errors}}, indent=2, default=str)
    text = frame.to_string(index=False)
    if errors:
        text += "\n\nErrors:\n" + "\n".join(errors)
    return text


@_tool
def list_gallery(format: str = "text") -> str:
    """The gallery: ready-made models, each a finished network with no code.

    Every entry is a factory in `line_solver.gallery`, listed with the exact
    arguments it takes, so the catalogue is generated from the module and can
    never drift from it. Between them they cover shapes the parameter-driven
    tools cannot express at all: fork-join, caches, layered networks,
    finite-capacity regions, random environments and the classical
    arrival/service matrix (M/M/1, Erlang, Coxian, hyperexponential, MAP, ...).

    Build one with build_from_gallery, then solve it with solve_model_file.

    Args:
        format: text or json.

    Returns:
        The catalogue: name, signature and one-line description of each entry.
    """
    import inspect

    from . import gallery

    entries = []
    for name in sorted(n for n in dir(gallery) if n.startswith("gallery_")):
        fn = getattr(gallery, name)
        if not callable(fn):
            continue
        try:
            params = [
                {"name": p.name,
                 "default": None if p.default is inspect.Parameter.empty else p.default}
                for p in inspect.signature(fn).parameters.values()
            ]
            signature = str(inspect.signature(fn))
        except (TypeError, ValueError):
            params, signature = [], "()"
        doc = (inspect.getdoc(fn) or "").strip().split("\n")[0]
        entries.append({"name": name, "signature": signature,
                        "params": params, "description": doc})

    if format == "json":
        return json.dumps({"results": entries, "metadata": {"count": len(entries)}},
                          indent=2, default=str)

    lines = ["Gallery models (%d)" % len(entries), "=" * 40, ""]
    for e in entries:
        lines.append("%s%s" % (e["name"], e["signature"]))
        if e["description"]:
            lines.append("    %s" % e["description"])
    lines.append("")
    lines.append('Build one with build_from_gallery(name="gallery_cqn", params=\'{"M": 3}\').')
    return "\n".join(lines)


@_tool
def build_from_gallery(name: str, params: str = "", format: str = "text") -> str:
    """Build a ready-made gallery model and keep it for later analysis.

    The shortest path to a non-trivial model: no stations to declare, no
    routing to write. list_gallery names every entry and its arguments.

    Args:
        name: The factory name, with or without the `gallery_` prefix.
        params: JSON object of arguments for it, e.g. '{"M": 3, "useDelay": true}'.
            Omit for the factory's own defaults.
        format: text or json.

    Returns:
        The model id and a summary of what was built. Pass the id to
        solve_model_file, find_solver, visualize_model or convert_model.
    """
    from . import gallery

    try:
        fname = name if name.startswith("gallery_") else "gallery_" + name
        fn = getattr(gallery, fname, None)
        if fn is None or not callable(fn):
            raise ValueError("unknown gallery model '%s'; list_gallery() names "
                             "every entry" % name)

        kwargs = {}
        if params.strip():
            kwargs = json.loads(params)
            if not isinstance(kwargs, dict):
                raise ValueError('`params` is a JSON object of the factory\'s '
                                 'arguments, e.g. {"M": 3}')

        model = fn(**kwargs)
    except Exception as e:
        return _error("Error: %s" % e, format)

    origin = "%s(%s)" % (fname, params.strip() or "")
    return _render_built(model, _register_model(model, origin), origin, format)


@_tool
def build_network(
    stations: str,
    classes: str,
    service: str = "",
    arrival: str = "",
    routing: str = "serial",
    name: str = "Model",
    format: str = "text",
) -> str:
    """Build ANY ordinary queueing network from a description, with no code.

    Open, closed and mixed; any number of stations; any number of classes; any
    scheduling; any distribution the library publishes. This is the general
    case that analyze_queue (one station) and analyze_network (a demand matrix,
    with the topology supplying the routing) cannot express.

    Args:
        stations: JSON list of station objects, in routing order. Each has a
            `name` and a `type` (Source, Sink, Queue, Delay, Router,
            ClassSwitch, Fork, Join), and may carry `scheduling` (FCFS, PS,
            INF, ... by NAME), `servers` and `capacity`. A Join ALSO needs
            `fork`, naming the Fork station it closes; the pairing is a
            declaration and is not derivable from the routing. Example:
            '[{"name":"Src","type":"Source"},
              {"name":"Q1","type":"Queue","scheduling":"PS","servers":2},
              {"name":"Snk","type":"Sink"}]'
        classes: JSON list of class objects. Each has a `name` and a `type`
            (Open or Closed); a Closed class also needs `population` and
            `refstation`. Open and closed classes MAY appear together, which is
            what makes a mixed network. `priority` is optional. Example:
            '[{"name":"web","type":"Open"},
              {"name":"batch","type":"Closed","population":5,"refstation":"Think"}]'
        service: JSON map station -> class -> distribution, e.g.
            '{"Q1": {"web": "exp(4.0)", "batch": "erlang(2.0, 3)"}}'.
            Distributions are written name(args); list_solvers is unrelated to
            this, the names are the library's own distribution classes.
        arrival: JSON map for the open classes at a Source, either
            '{"Src": {"web": "exp(1.0)"}}' or, when there is exactly one
            Source, the shorter '{"web": "exp(1.0)"}'.
        routing: "serial" chains the stations in the order given, per class
            (closing the cycle unless the last is a Sink). Otherwise a JSON
            list of hops, each [from, to, class, probability] or
            [from, to, from_class, to_class, probability] to switch class on
            the hop.
        name: Model name.
        format: text or json.

    Returns:
        The model id and a summary. Pass the id to solve_model_file for any
        analysis, or to find_solver, visualize_model or convert_model.

    Note:
        Caches, Petri-net places and transitions, finite-capacity regions and
        layered networks are NOT built here. Reach those through
        build_from_gallery, or hand a LINE json document to
        solve_model_file(content=...).
    """
    try:
        from .lang.classes import ClosedClass, OpenClass
        from .lang.network import Network
        from .lang.nodes import (ClassSwitch, Delay, Fork, Join, Queue, Router,
                                 Sink, Source)

        def _load(text, what, default=None):
            if not str(text).strip():
                return default
            try:
                return json.loads(text)
            except json.JSONDecodeError as e:
                raise ValueError("`%s` is not valid JSON: %s" % (what, e))

        station_specs = _load(stations, "stations")
        class_specs = _load(classes, "classes")
        if not isinstance(station_specs, list) or not station_specs:
            raise ValueError("`stations` is a non-empty JSON list of station objects")
        if not isinstance(class_specs, list) or not class_specs:
            raise ValueError("`classes` is a non-empty JSON list of class objects")

        model = Network(name)

        builders = {
            "source": lambda m, n, s: Source(m, n),
            "sink": lambda m, n, s: Sink(m, n),
            "delay": lambda m, n, s: Delay(m, n),
            "router": lambda m, n, s: Router(m, n),
            "classswitch": lambda m, n, s: ClassSwitch(m, n),
            "fork": lambda m, n, s: Fork(m, n),
            "join": lambda m, n, s: Join(m, n),
            "queue": lambda m, n, s: Queue(m, n, _sched(s)),
        }

        def _sched(spec):
            try:
                return _sched_by_name(spec.get("scheduling", "FCFS"))
            except ValueError as e:
                raise ValueError("station '%s': %s" % (spec.get("name"), e))

        nodes = {}
        order = []
        for spec in station_specs:
            if not isinstance(spec, dict) or "name" not in spec or "type" not in spec:
                raise ValueError("each station needs a `name` and a `type`, got %r" % (spec,))
            nname, ntype = str(spec["name"]), str(spec["type"]).strip().lower()
            if nname in nodes:
                raise ValueError("duplicate station name '%s'" % nname)
            build = builders.get(ntype)
            if build is None:
                raise ValueError(
                    "station '%s': unsupported type '%s'. This tool builds %s; for "
                    "caches, Petri nets or layered models use build_from_gallery or "
                    "hand a document to solve_model_file."
                    % (nname, spec["type"], ", ".join(sorted(builders))))
            node = build(model, nname, spec)
            if spec.get("servers") is not None and ntype == "queue":
                node.setNumberOfServers(int(spec["servers"]))
            if spec.get("capacity") is not None:
                node.setCapacity(float(spec["capacity"]))
            nodes[nname] = node
            order.append(node)

        def _station(nname, what):
            if nname not in nodes:
                raise ValueError("%s: unknown station '%s'; declared: %s"
                                 % (what, nname, ", ".join(nodes)))
            return nodes[nname]

        # A Join must name the Fork it closes; the pairing is a declaration and
        # cannot be read off the routing when a model nests two fork-join pairs.
        # Bound after the loop so a Join may be declared before its Fork.
        for spec in station_specs:
            if str(spec["type"]).strip().lower() != "join":
                continue
            jname = str(spec["name"])
            fname = spec.get("fork", spec.get("forkNode"))
            if fname is None:
                raise ValueError(
                    "join '%s' needs a `fork` naming the fork station it closes"
                    % jname)
            nodes[jname]._fork = _station(str(fname), "join '%s' fork" % jname)

        jobclasses = {}
        for spec in class_specs:
            if not isinstance(spec, dict) or "name" not in spec or "type" not in spec:
                raise ValueError("each class needs a `name` and a `type`, got %r" % (spec,))
            cname, ctype = str(spec["name"]), str(spec["type"]).strip().lower()
            if cname in jobclasses:
                raise ValueError("duplicate class name '%s'" % cname)
            prio = int(spec.get("priority", 0))
            if ctype == "open":
                jobclasses[cname] = OpenClass(model, cname, prio)
            elif ctype == "closed":
                if spec.get("population") is None or spec.get("refstation") is None:
                    raise ValueError("closed class '%s' needs `population` and "
                                     "`refstation`" % cname)
                jobclasses[cname] = ClosedClass(
                    model, cname, int(spec["population"]),
                    _station(str(spec["refstation"]), "class '%s' refstation" % cname),
                    prio)
            else:
                raise ValueError("class '%s': type must be Open or Closed, got '%s'"
                                 % (cname, spec["type"]))

        def _jobclass(cname, what):
            if cname not in jobclasses:
                raise ValueError("%s: unknown class '%s'; declared: %s"
                                 % (what, cname, ", ".join(jobclasses)))
            return jobclasses[cname]

        for sname, per_class in (_load(service, "service", {}) or {}).items():
            node = _station(sname, "service")
            if not isinstance(per_class, dict):
                raise ValueError("service['%s'] is a map of class -> distribution" % sname)
            for cname, dist in per_class.items():
                node.setService(_jobclass(cname, "service['%s']" % sname),
                                _parse_distribution(dist, "service['%s']['%s']"
                                                    % (sname, cname)))

        arrivals = _load(arrival, "arrival", {}) or {}
        if arrivals and all(not isinstance(v, dict) for v in arrivals.values()):
            sources = [n for n, o in nodes.items() if isinstance(o, Source)]
            if len(sources) != 1:
                raise ValueError("`arrival` keyed by class needs exactly one Source; "
                                 "found %d. Key it by station instead." % len(sources))
            arrivals = {sources[0]: arrivals}
        for sname, per_class in arrivals.items():
            node = _station(sname, "arrival")
            for cname, dist in per_class.items():
                node.setArrival(_jobclass(cname, "arrival['%s']" % sname),
                                _parse_distribution(dist, "arrival['%s']['%s']"
                                                    % (sname, cname)))

        if str(routing).strip().lower() == "serial":
            model.link(Network.serialRouting(order))
        else:
            hops = _load(routing, "routing")
            if not isinstance(hops, list) or not hops:
                raise ValueError('`routing` is "serial" or a non-empty JSON list of '
                                 'hops [from, to, class, probability]')
            P = model.initRoutingMatrix()
            for hop in hops:
                if len(hop) == 4:
                    frm, to, cls, prob = hop
                    src_cls = dst_cls = str(cls)
                elif len(hop) == 5:
                    frm, to, src_cls, dst_cls, prob = hop
                    src_cls, dst_cls = str(src_cls), str(dst_cls)
                else:
                    raise ValueError("routing hop %r must be [from, to, class, prob] "
                                     "or [from, to, from_class, to_class, prob]" % (hop,))
                P.set(_jobclass(src_cls, "routing"), _jobclass(dst_cls, "routing"),
                      _station(str(frm), "routing"), _station(str(to), "routing"),
                      float(prob))
            model.link(P)
    except Exception as e:
        return _error("Error: %s" % e, format)

    origin = "build_network(%d stations, %d classes)" % (len(nodes), len(jobclasses))
    return _render_built(model, _register_model(model, origin), origin, format)


@_tool
def build_topology(
    kind: str,
    demands: str,
    population: str = "",
    think_time: str = "",
    arrival_rates: str = "",
    scheduling: str = "PS",
    servers: str = "",
    format: str = "text",
) -> str:
    """Build a multi-station network from a DEMAND MATRIX, with no routing.

    This is the textbook (N, Z, D) input: one service demand per station per
    class, and the topology supplies the routing itself. It reaches open,
    closed and mixed multi-station, multi-class networks without a routing
    grammar, which is what the single-station analyze_queue cannot do.

    It REGISTERS the model rather than solving it; analyze_network takes the
    same arguments and answers with the numbers instead.

    Args:
        kind: tandem (open chain), cyclic (closed loop), cluster (open fan-out),
            cluster_closed, or cluster_mixed (open and closed classes together).
        demands: The demand matrix D as JSON, rows = stations, columns =
            classes, e.g. '[[0.5, 0.3], [0.2, 0.4]]'. For cluster_mixed the
            columns are the open classes first, then the closed ones.
        population: JSON list of job counts per CLOSED class, e.g. '[10]'.
            Required by cyclic, cluster_closed and cluster_mixed.
        think_time: JSON list of think times per closed class, e.g. '[5.0]'.
            Used by cluster_closed and cluster_mixed.
        arrival_rates: JSON list of arrival rates per OPEN class, e.g. '[1.0]'.
            Required by tandem, cluster and cluster_mixed.
        scheduling: Scheduling at every queue: FCFS, PS, INF and the rest of
            SchedStrategy, by NAME. A per-station JSON list is accepted too.
        servers: JSON list of server counts per station, e.g. '[1, 2]'.
            Omit for single-server stations.
        format: text or json.

    Returns:
        The model id and a summary of what was built.
    """
    try:
        model, k, D = _topology_model(kind, demands, population, think_time,
                                      arrival_rates, scheduling, servers)
    except Exception as e:
        return _error("Error: %s" % e, format)

    origin = "%s(D=%dx%d)" % (k, D.shape[0], D.shape[1])
    return _render_built(model, _register_model(model, origin), origin, format)


def _topology_model(kind: str, demands: str, population: str = "",
                    think_time: str = "", arrival_rates: str = "",
                    scheduling: str = "PS", servers: str = ""):
    """Build the (N, Z, D) topology both `build_topology` and `analyze_network` take.

    Single-sourced on purpose: the builder registers the model and the analyzer
    solves it, and a second copy of this parsing would be a second grammar for
    the same five arguments.

    Returns:
        (model, resolved kind, the demand matrix D).
    """
    import numpy as np

    from .lang.network import Network

    def _arr(text, what):
        if not str(text).strip():
            return None
        value = json.loads(text)
        return np.array(value, dtype=float)

    D = _arr(demands, "demands")
    if D is None:
        raise ValueError("`demands` is required: the demand matrix D as "
                         "JSON, rows = stations, columns = classes")
    if D.ndim != 2:
        raise ValueError("`demands` must be a matrix (a list of per-station "
                         "rows), got %d dimension(s)" % D.ndim)

    N = _arr(population, "population")
    Z = _arr(think_time, "think_time")
    rates = _arr(arrival_rates, "arrival_rates")
    S = _arr(servers, "servers")

    nstations = D.shape[0]
    sched_names = json.loads(scheduling) if scheduling.strip().startswith("[") \
        else [scheduling] * nstations
    if len(sched_names) != nstations:
        raise ValueError("`scheduling` names %d strategies but `demands` has "
                         "%d stations" % (len(sched_names), nstations))
    strategies = [_sched_by_name(s) for s in sched_names]

    kinds = {"tandem", "cyclic", "cluster", "cluster_closed", "cluster_mixed"}
    k = kind.strip().lower()
    if k not in kinds:
        raise ValueError("unknown kind '%s'; use one of: %s"
                         % (kind, ", ".join(sorted(kinds))))

    def _need(value, argname):
        if value is None:
            raise ValueError("`%s` is required for kind '%s'" % (argname, k))
        return value

    if k == "tandem":
        model = Network.tandem(_need(rates, "arrival_rates"), D, strategies)
    elif k == "cyclic":
        model = Network.cyclic(_need(N, "population"), D, strategies, S)
    elif k == "cluster":
        model = Network.cluster(_need(rates, "arrival_rates"), D, strategies, S)
    elif k == "cluster_closed":
        model = Network.cluster_closed(_need(N, "population"),
                                       _need(Z, "think_time"), D, strategies, S)
    else:
        model = Network.cluster_mixed(_need(rates, "arrival_rates"),
                                      _need(N, "population"),
                                      _need(Z, "think_time"), D, strategies, S)

    return model, k, D


# ---------------------------------------------------------------------------
# 4b. The analyze family: build a model from arguments and solve it in one call
# ---------------------------------------------------------------------------
#
# These two tools are the quick-start surface: they take parameters rather than
# a model, so nothing has to be built, registered and then solved in a second
# call. What they must NOT be is a smaller world than the rest of the server.
# The model side takes any distribution the library publishes, any scheduling
# discipline, any number of classes, and open, closed or mixed populations; the
# solve side is `_run_cli_analysis`, the same dispatch `solve_model_file` uses,
# so every solver, method and analysis is reachable from here too.

_KENDALL = {
    "Exp": "M", "Det": "D", "Erlang": "E", "HyperExp": "H", "APH": "PH",
    "PH": "PH", "Coxian": "Cox", "Cox2": "Cox", "MAP": "MAP", "MMPP2": "MMPP2",
    "Uniform": "U", "Gamma": "Gam", "Pareto": "Par", "Weibull": "W",
    "Lognormal": "LN", "Replayer": "Trace", "Immediate": "I",
}


def _kendall(dist) -> str:
    """The Kendall letter of a distribution, or its own name when it has none."""
    return _KENDALL.get(type(dist).__name__, type(dist).__name__)


def _dist_or_rate(spec, rate, what: str):
    """Resolve a distribution given either as `name(args)` or as a bare RATE.

    The rate shorthand is what keeps `analyze_queue(arrival_rate=2)` a one-line
    call; the spec is what makes every other distribution reachable from the
    same argument. The spec wins when both are given, being the more specific
    statement.

    Returns:
        (distribution, the text to report it by), or (None, "") when neither
        was given.
    """
    if spec is not None and str(spec).strip():
        return _parse_distribution(spec, what), str(spec).strip()
    if rate:
        from .distributions import Exp

        if float(rate) <= 0:
            raise ValueError("%s: a rate must be > 0, got %s" % (what, rate))
        return Exp(float(rate)), "exp(%s)" % rate
    return None, ""


def _dist_or_mean(spec, mean, what: str):
    """As `_dist_or_rate`, for an argument stated as a MEAN rather than a rate.

    Think time is a time, so `think_time=5` is Exp(0.2) and not Exp(5). Getting
    that backwards is silent: the model still solves, at a think time 25 times
    too small.
    """
    if spec is not None and str(spec).strip():
        return _parse_distribution(spec, what), str(spec).strip()
    if mean:
        from .distributions import Exp

        if float(mean) < 0:
            raise ValueError("%s: a mean time must be >= 0, got %s" % (what, mean))
        return Exp(1.0 / float(mean)), "exp(%g)" % (1.0 / float(mean))
    return None, ""


def _class_spec(spec: dict, default_name: str) -> dict:
    """Normalise one class of the analyze family into a single shape."""
    if not isinstance(spec, dict):
        raise ValueError("each class is a JSON object, got %r" % (spec,))
    name = str(spec.get("name", default_name))
    service, service_text = _dist_or_rate(
        spec.get("service"), spec.get("service_rate", 0),
        "class '%s' service" % name)
    if service is None:
        raise ValueError("class '%s' needs a `service` distribution or a "
                         "`service_rate`" % name)
    arrival, arrival_text = _dist_or_rate(
        spec.get("arrival"), spec.get("arrival_rate", 0),
        "class '%s' arrival" % name)
    population = int(spec.get("population", 0) or 0)
    think, think_text = _dist_or_mean(
        spec.get("think"), spec.get("think_time", 0),
        "class '%s' think time" % name)

    if arrival is not None and population:
        raise ValueError("class '%s' is both open (an arrival) and closed (a "
                         "population); give one or the other" % name)
    if arrival is None and not population:
        raise ValueError("class '%s' needs an `arrival` or `arrival_rate` "
                         "(open) or a `population` (closed)" % name)
    if arrival is not None and think is not None:
        raise ValueError("class '%s' is open, so it has no think time; a think "
                         "time belongs to a closed class" % name)

    return {
        "name": name,
        "kind": "open" if arrival is not None else "closed",
        "service": service, "service_text": service_text,
        "arrival": arrival, "arrival_text": arrival_text,
        "population": population,
        "think": think, "think_text": think_text,
        "priority": int(spec.get("priority", 0) or 0),
    }


def _queue_class_specs(classes: str, service: str, service_rate: float,
                       arrival: str, arrival_rate: float,
                       population: int, think_time: float) -> list:
    """The class list of `analyze_queue`, from `classes` or from the scalars."""
    if classes.strip():
        scalars = [n for n, v in (("service", service), ("arrival", arrival),
                                  ("service_rate", service_rate),
                                  ("arrival_rate", arrival_rate),
                                  ("population", population),
                                  ("think_time", think_time)) if v]
        if scalars:
            # Silently ignoring them would answer a different model than the
            # one the caller described.
            raise ValueError("`classes` describes every class, so %s cannot be "
                             "given beside it; move it into the class object"
                             % ", ".join("`%s`" % n for n in scalars))
        specs = json.loads(classes)
        if not isinstance(specs, list) or not specs:
            raise ValueError('`classes` is a non-empty JSON list of class '
                             'objects, e.g. [{"name":"web","arrival_rate":1.0,'
                             '"service_rate":4.0}]')
        return [_class_spec(s, "Class%d" % (i + 1)) for i, s in enumerate(specs)]

    return [_class_spec({"service": service, "service_rate": service_rate,
                         "arrival": arrival, "arrival_rate": arrival_rate,
                         "population": population, "think_time": think_time},
                        "Class1")]


def _queue_utilization(specs: list, servers: int):
    """Offered load of the open classes at the station, or None if unknowable.

    Computed from the MEANS, so it holds for any distribution and not only for
    the exponential pair this tool used to be limited to. A model with no open
    class has no offered load in this sense and answers None.
    """
    if not any(s["kind"] == "open" for s in specs):
        return None
    try:
        rho = 0.0
        for spec in specs:
            if spec["kind"] == "open":
                rho += spec["service"].getMean() / spec["arrival"].getMean()
        return rho / servers
    except Exception:
        return None


def _queue_title(specs: list, servers: int, capacity: int) -> str:
    """Name the model in Kendall notation where that is honest, in words where not."""
    closed = [s for s in specs if s["kind"] == "closed"]
    opened = [s for s in specs if s["kind"] == "open"]
    if len(specs) == 1 and opened:
        return "%s/%s/%d%s" % (_kendall(specs[0]["arrival"]),
                               _kendall(specs[0]["service"]), servers,
                               "/%d" % capacity if capacity > 0 else "")
    if len(specs) == 1:
        return "closed %s/%d, N=%d" % (_kendall(specs[0]["service"]),
                                       servers, closed[0]["population"])
    kind = "mixed" if (opened and closed) else ("open" if opened else "closed")
    return "%s queue, %d classes, k=%d" % (kind, len(specs), servers)


def _sched_by_name(name):
    """Resolve a scheduling discipline BY NAME, over the whole SchedStrategy set.

    Four names used to be accepted by the analyze tools (FCFS, PS, LCFS, INF),
    which left them unable to express disciplines the engines have served for
    years. Resolution is by NAME because SchedStrategy is numbered differently
    in `lang.base` and in `constants`, so an integer here would mean two
    different disciplines depending on which enum was consulted.
    """
    from .lang.base import SchedStrategy

    resolved = getattr(SchedStrategy, str(name).strip().upper(), None)
    if resolved is None:
        published = sorted(n for n in dir(SchedStrategy) if n.isupper())
        raise ValueError("unknown scheduling strategy '%s'; SchedStrategy "
                         "publishes: %s" % (name, ", ".join(published)))
    return resolved


def _sanitize_json(value):
    """Make NaN and Inf legal JSON, which `json.dumps` alone does not."""
    if isinstance(value, float):
        if math.isnan(value):
            return None
        if math.isinf(value):
            return "Inf" if value > 0 else "-Inf"
        return value
    if isinstance(value, dict):
        return {k: _sanitize_json(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_sanitize_json(v) for v in value]
    return value


def _analysis_header(title: str, rows: list, warnings: list) -> list:
    """The text-mode header: the warnings, the model, then one row per line."""
    head = ["WARNING: %s" % w for w in warnings]
    if warnings:
        head.append("")
    head.append("Model: %s" % title)
    width = max([len(str(label)) for label, _ in rows] or [0])
    for label, value in rows:
        head.append("  %-*s = %s" % (width, label, value))
    head.append("")
    return head


def _render_analysis(results: dict, metadata: dict, header: list, format: str) -> str:
    """Render a CLI results dict under the analyze family's own contract.

    One analysis answers with its table alone, under a header naming the model;
    a comma-separated list answers keyed by analysis, the way the CLI does. The
    json reply keeps the `results` plus `metadata` shape the rest of this
    server uses, rather than the CLI's analysis-keyed document.
    """
    from . import cli as line_cli

    values = {k: v for k, v in results.items() if v is not None}

    # `execute_single_analysis` returns a failure AS the result, a string
    # beginning "Error: ". That is the CLI's error channel, so it is rendered
    # as this server's error rather than as a result a json caller would have
    # to sniff for -- and than one sweep_parameter could not read back at all.
    failed = [v for v in values.values() if isinstance(v, str) and v.startswith("Error: ")]
    if values and len(failed) == len(values):
        return _error("\n".join(failed), format)

    single = next(iter(values.values())) if len(values) == 1 else None

    if format == "json":
        if single is not None:
            payload = _sanitize_json(line_cli._convert_to_json_serializable(single))
        else:
            payload = {k: _sanitize_json(line_cli._convert_to_json_serializable(v))
                       for k, v in values.items()}
        return json.dumps({"results": payload, "metadata": metadata},
                          indent=2, default=str)

    if not values:
        body = "(no results: the analysis produced no output)"
    elif single is not None and not isinstance(single, dict):
        body = line_cli._format_single_result(single)
    else:
        body = line_cli.format_readable(values)
    return "\n".join(header) + body


def _selected_solver_name(solver_obj) -> str:
    """Which solver actually RAN, or "" when no choice was recorded.

    Every request goes through SolverAUTO, so the object's own class name is
    always SolverAUTO and says nothing: the selection it made is what a caller
    who wrote `auto`, `fast` or `exact` needs back. `getSelectedSolverName()`
    is NOT that answer, because when nothing has recorded a choice it runs the
    heuristic and returns what AUTO WOULD pick -- which, for an analysis served
    by a family accessor rather than by runAnalyzer, names a solver that never
    ran. The recorded choice is read instead and its absence reported as one.
    """
    name = getattr(solver_obj, "_selected_solver_name", None)
    if name:
        return str(name)
    inner = getattr(solver_obj, "_selected_solver", None)
    if inner is not None:
        return type(inner).__name__.replace("SolverNative", "").replace("Solver", "")
    return ""


def _analyze_and_render(model, title: str, rows: list, warnings: list,
                        metadata: dict, solver: str, analysis: str,
                        format: str, knobs: dict) -> str:
    """Solve a model the analyze family built, and render it or say why not."""
    try:
        results, _, solver_obj = _run_cli_analysis(
            lambda verbose: model, "silent", solver=solver, analysis=analysis,
            **knobs)
    except Exception as e:
        return _error("Error: %s" % e, format)

    # `auto` is a request, not an answer: which solver it resolved to is part
    # of the result and is reported rather than left to be guessed.
    used = _selected_solver_name(solver_obj)
    shown = solver.upper() if not used or used.upper() == solver.upper() \
        else "%s -> %s" % (solver, used)
    metadata["solver"] = solver.upper()
    if used:
        metadata["solver_used"] = used
    metadata["analysis"] = analysis
    return _render_analysis(
        results, metadata,
        _analysis_header(title, rows + [("solver", shown)], warnings), format)


_ANALYZE_KNOBS_DOC = """\
        solver: Any solver or solver.method the CLI takes: auto (the default,
            which picks one that can run the model as described), mva, nc,
            ctmc, mam, fld, ssa, ldes, jmt, lqns, ag, ba, uq, env, and the
            high-level names sim, exact, fast and accurate. list_solvers prints
            the vocabulary; find_solver says which ones fit a given model.
        method: Algorithm within that solver, also spelled solver.method.
        analysis: One analysis or a comma-separated list, exactly as
            solve_model_file takes them: avg (default), sys, stage, chain,
            node, nodechain, all; cdf-respt, cdf-passt, perct-respt; tran-avg;
            prob, prob-aggr, prob-marg; sample, sample-aggr; reward-value;
            normconst; busyperiod; sens; interval.
        node: Node index for the prob and sample analyses (0-based).
        class_idx: Job class index for prob-marg (0-based).
        state: State vector for prob, comma-separated integers.
        events: Events drawn by the sample analyses. Default 1000.
        percentiles: Percentiles of perct-respt. Default "50,90,95,99".
        reward_name: Built-in reward of reward-value: QLen, Tput, Util, RespT,
            WaitT, ArvR, ResidT.
        seed: Random seed of the stochastic solvers (JMT, SSA, LDES).
        samples: Simulation samples / Monte Carlo draws (SSA, LDES, JMT, UQ).
        cutoff: State-space cutoff per open class (CTMC, SSA).
        timespan: Horizon of the transient analyses: "T1", "T0,T1" or "T0:T1".
        options: JSON object of further solver options, e.g. {"stiff": false}.
        format: "text" (default) or "json" (the results plus the metadata)."""


def _knob_doc(fn):
    """Splice the shared solve-side argument docs into a tool's docstring.

    Applied BELOW @_tool, so the substitution happens before the docstring is
    read and registered as the tool's description. The alternative is two
    copies of the same twenty lines, which drift.
    """
    fn.__doc__ = (fn.__doc__ or "").replace("%KNOBS%", _ANALYZE_KNOBS_DOC)
    return fn


def _analyze_knobs(method, node, class_idx, state, events, percentiles,
                   reward_name, seed, samples, cutoff, timespan, options) -> dict:
    """Collect the solve-side arguments the analyze family shares."""
    return {
        "method": method, "node": node, "class_idx": class_idx, "state": state,
        "events": events, "percentiles": percentiles, "reward_name": reward_name,
        "seed": seed, "samples": samples, "cutoff": cutoff,
        "timespan": timespan, "options": options,
    }


@_tool
@_knob_doc
def analyze_queue(
    service: str = "",
    service_rate: float = 0.0,
    arrival: str = "",
    arrival_rate: float = 0.0,
    population: int = 0,
    think_time: float = 0.0,
    classes: str = "",
    servers: int = 1,
    queue_capacity: int = -1,
    scheduling: str = "FCFS",
    solver: str = "auto",
    method: str = "",
    analysis: str = "avg",
    node: int | None = None,
    class_idx: int | None = None,
    state: str = "",
    events: int = 1000,
    percentiles: str = "50,90,95,99",
    reward_name: str = "",
    seed: int | None = None,
    samples: int | None = None,
    cutoff: float | None = None,
    timespan: str = "",
    options: str = "",
    format: str = "text",
) -> str:
    """Build and solve ONE queueing station, without writing code.

    Open, closed or mixed: a class with an arrival is open and is fed by a
    Source, a class with a population is closed and circulates through a Delay
    holding its think time (or straight back into the queue when there is
    none), and both may appear in the same model. Any distribution the library
    publishes, any scheduling discipline, any number of servers, a finite
    capacity, class priorities, and any solver and analysis the CLI serves.

    M/M/1 is analyze_queue(arrival_rate=2, service_rate=5). A closed
    machine-repairman is analyze_queue(population=10, think_time=5,
    service_rate=2). Everything past those two is the arguments below.

    Args:
        service: Service distribution, written name(args): "exp(5.0)",
            "erlang(3.0, 2)", "hyperexp(0.5, 1.0, 4.0)", "det(0.25)",
            "gamma(2.0, 0.5)", "replayer('trace.txt')". Every distribution
            class the library publishes is accepted.
        service_rate: Shorthand for an exponential service, as a RATE (mu).
        arrival: Arrival distribution of an OPEN class, written the same way.
        arrival_rate: Shorthand for a Poisson arrival stream, as a RATE
            (lambda).
        population: Jobs of a CLOSED class circulating in the system (N).
        think_time: MEAN think time of that closed class (Z), spent at a Delay.
            0 sends the class straight back into the queue.
        classes: A JSON list of class objects, for a multiclass or mixed model,
            REPLACING the six arguments above. Each object takes "name",
            "service" or "service_rate", then either "arrival"/"arrival_rate"
            (open) or "population" with an optional "think_time"/"think"
            (closed), plus an optional "priority". Example:
            '[{"name":"web","arrival_rate":1.0,"service_rate":4.0},
              {"name":"batch","population":4,"think_time":2.0,
               "service":"erlang(2.0, 3)"}]'
        servers: Number of servers (k). Default 1.
        queue_capacity: Total capacity including the servers. -1 (default) is
            infinite; a finite value is the /N of the Kendall notation.
        scheduling: Scheduling discipline BY NAME, over the whole SchedStrategy
            set: FCFS, PS, LCFS, INF, SIRO, HOL, SJF, LJF, SEPT, LEPT, SRPT,
            DPS, GPS, FCFSPR, LCFSPR and the rest. Default FCFS.
%KNOBS%

    Returns:
        The requested analysis, under a header naming the model that was built.
    """
    warnings: list = []
    try:
        specs = _queue_class_specs(classes, service, service_rate, arrival,
                                   arrival_rate, population, think_time)
        if servers < 1:
            raise ValueError("`servers` must be >= 1")
        sched = _sched_by_name(scheduling)

        from .lang.classes import ClosedClass, OpenClass
        from .lang.network import Network
        from .lang.nodes import Delay, Queue, Sink, Source

        has_open = any(s["kind"] == "open" for s in specs)
        has_think = any(s["kind"] == "closed" and s["think"] is not None
                        for s in specs)
        rho = _queue_utilization(specs, servers)
        if rho is not None and queue_capacity < 0 and rho >= 1.0 \
                and all(s["kind"] == "open" for s in specs):
            # Only a purely open model is unstable at rho >= 1: a closed class
            # bounds the queue by its own population whatever the open load is.
            return _error(
                "Error: system is unstable (rho = %.4f >= 1). Reduce the "
                "arrival rate, raise the service rate or the server count, or "
                "set a finite queue_capacity." % rho, format)

        model = Network(_queue_title(specs, servers, queue_capacity))
        queue = Queue(model, "Queue", sched)
        if servers > 1:
            queue.setNumberOfServers(servers)
        if queue_capacity > 0:
            queue.setCapacity(queue_capacity)

        source = sink = delay = None
        if has_open:
            source, sink = Source(model, "Source"), Sink(model, "Sink")
        if has_think:
            delay = Delay(model, "Think")

        P = model.initRoutingMatrix()
        for spec in specs:
            if spec["kind"] == "open":
                jobclass = OpenClass(model, spec["name"], spec["priority"])
                source.setArrival(jobclass, spec["arrival"])
                queue.setService(jobclass, spec["service"])
                P.set(jobclass, jobclass, source, queue, 1.0)
                P.set(jobclass, jobclass, queue, sink, 1.0)
            elif spec["think"] is not None:
                jobclass = ClosedClass(model, spec["name"], spec["population"],
                                       delay, spec["priority"])
                delay.setService(jobclass, spec["think"])
                queue.setService(jobclass, spec["service"])
                P.set(jobclass, jobclass, delay, queue, 1.0)
                P.set(jobclass, jobclass, queue, delay, 1.0)
            else:
                # No think time: the class re-enters the queue it just left,
                # rather than passing through a Delay whose mean stands in for
                # zero. The closed tool used Exp(1e9) there, which is a think
                # time of 1e-9 and not the absence of one.
                jobclass = ClosedClass(model, spec["name"], spec["population"],
                                       queue, spec["priority"])
                queue.setService(jobclass, spec["service"])
                P.set(jobclass, jobclass, queue, queue, 1.0)
        model.link(P)
    except Exception as e:
        return _error("Error: %s" % e, format)

    knobs = _analyze_knobs(method, node, class_idx, state, events, percentiles,
                           reward_name, seed, samples, cutoff, timespan, options)
    if solver.lower().startswith("ctmc") and knobs["cutoff"] is None \
            and has_open and queue_capacity < 0:
        # An open CTMC needs a truncation and this model carries none, so the
        # tool supplies the one it has always used rather than refusing.
        closed = sum(s["population"] for s in specs if s["kind"] == "closed")
        knobs["cutoff"] = max(50, closed, int(10 * (rho or 1)))
        warnings.append("CTMC on an infinite-capacity queue truncates the "
                        "state space at cutoff=%d." % knobs["cutoff"])
    if solver.lower().startswith("fld"):
        warnings.append("FLD (fluid) is an approximation; results may differ "
                        "from an exact solver.")

    rows = [("scheduling", scheduling.upper()), ("servers (k)", servers)]
    if queue_capacity > 0:
        rows.append(("capacity", queue_capacity))
    for spec in specs:
        if spec["kind"] == "open":
            rows.append(("class %s" % spec["name"],
                         "open, arrival %s, service %s"
                         % (spec["arrival_text"], spec["service_text"])))
        else:
            rows.append(("class %s" % spec["name"],
                         "closed, N=%d, think %s, service %s"
                         % (spec["population"], spec["think_text"] or "none",
                            spec["service_text"])))
    if rho is not None:
        rows.append(("utilization (rho)", "%.6f" % rho))

    metadata = {
        "model": model.getName(),
        "stations": 1 + (2 if has_open else 0) + (1 if has_think else 0),
        "servers": servers,
        "queue_capacity": queue_capacity,
        "scheduling": scheduling.upper(),
        "classes": [{"name": s["name"], "type": s["kind"],
                     "service": s["service_text"],
                     "arrival": s["arrival_text"] or None,
                     "population": s["population"] or None,
                     "think_time": s["think_text"] or None,
                     "priority": s["priority"]} for s in specs],
    }
    if rho is not None:
        metadata["utilization"] = round(rho, 6)
    if warnings:
        metadata["warnings"] = warnings

    return _analyze_and_render(model, model.getName(), rows, warnings, metadata,
                               solver, analysis, format, knobs)


# ---------------------------------------------------------------------------
# 5. Multi-station quick analysis
# ---------------------------------------------------------------------------

@_tool
@_knob_doc
def analyze_network(
    kind: str,
    demands: str,
    population: str = "",
    think_time: str = "",
    arrival_rates: str = "",
    scheduling: str = "PS",
    servers: str = "",
    solver: str = "auto",
    method: str = "",
    analysis: str = "avg",
    node: int | None = None,
    class_idx: int | None = None,
    state: str = "",
    events: int = 1000,
    percentiles: str = "50,90,95,99",
    reward_name: str = "",
    seed: int | None = None,
    samples: int | None = None,
    cutoff: float | None = None,
    timespan: str = "",
    options: str = "",
    format: str = "text",
) -> str:
    """Build and solve a MULTI-STATION network from a demand matrix, in one call.

    The textbook (N, Z, D) input: one service demand per station per class,
    with the topology supplying the routing. Open, closed and mixed, any
    number of stations and classes, and any solver and analysis the CLI
    serves. It is where analyze_queue stops, which is at one station.

    build_topology plus solve_model_file, for when the model is wanted only
    for its numbers. Build it instead when it is to be kept and then
    visualized, converted, or solved several ways.

    Args:
        kind: tandem (open chain), cyclic (closed loop), cluster (open
            fan-out), cluster_closed, or cluster_mixed (open and closed classes
            together).
        demands: The demand matrix D as JSON, rows = stations, columns =
            classes, e.g. '[[0.5, 0.3], [0.2, 0.4]]'. For cluster_mixed the
            columns are the open classes first, then the closed ones.
        population: JSON list of job counts per CLOSED class, e.g. '[10]'.
            Required by cyclic, cluster_closed and cluster_mixed.
        think_time: JSON list of think times per closed class, e.g. '[5.0]'.
            Used by cluster_closed and cluster_mixed.
        arrival_rates: JSON list of arrival rates per OPEN class, e.g. '[1.0]'.
            Required by tandem, cluster and cluster_mixed.
        scheduling: Scheduling at every queue, BY NAME (FCFS, PS, INF, ...).
            A per-station JSON list is accepted too. Default PS.
        servers: JSON list of server counts per station, e.g. '[1, 2]'.
            Omit for single-server stations.
%KNOBS%

    Returns:
        The requested analysis, under a header naming the model that was built.
    """
    try:
        model, resolved, D = _topology_model(kind, demands, population,
                                             think_time, arrival_rates,
                                             scheduling, servers)
    except Exception as e:
        return _error("Error: %s" % e, format)

    knobs = _analyze_knobs(method, node, class_idx, state, events, percentiles,
                           reward_name, seed, samples, cutoff, timespan, options)
    title = "%s(D=%dx%d)" % (resolved, D.shape[0], D.shape[1])
    rows = [("kind", resolved), ("stations", D.shape[0]),
            ("classes", D.shape[1]), ("scheduling", scheduling)]
    for label, value in (("population (N)", population),
                         ("think time (Z)", think_time),
                         ("arrival rates", arrival_rates),
                         ("servers", servers)):
        if str(value).strip():
            rows.append((label, str(value).strip()))

    metadata = {"model": title, "kind": resolved, "stations": int(D.shape[0]),
                "classes": int(D.shape[1]), "scheduling": scheduling,
                "demands": D.tolist()}
    return _analyze_and_render(model, title, rows, [], metadata, solver,
                               analysis, format, knobs)


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


@_tool
def sweep_parameter(
    param: str,
    values: str,
    service: str = "",
    service_rate: float = 0.0,
    arrival: str = "",
    arrival_rate: float = 0.0,
    population: int = 0,
    think_time: float = 0.0,
    classes: str = "",
    servers: int = 1,
    queue_capacity: int = -1,
    scheduling: str = "FCFS",
    solver: str = "auto",
    method: str = "",
    seed: int | None = None,
    samples: int | None = None,
    format: str = "text",
) -> str:
    """Sweep one parameter of a single station and collect the metrics.

    The model is the one analyze_queue builds, and every argument below the
    first two is that tool's, with the same meaning: open, closed or mixed,
    any distribution, any discipline, any solver. `param` then names which of
    the numeric ones the sweep varies, and the rest stay at the values given.

    To sweep a field of a model built by build_network, build_topology or the
    gallery, use sweep_model instead: it sweeps the DESCRIPTION rather than
    this fixed shape.

    Args:
        param: What to vary: arrival_rate, service_rate, servers, population,
            think_time or queue_capacity.
        values: The values, comma-separated ("0.5,1.0,1.5") or as
            start:stop:step ("0.1:2.0:0.1"). At most 100 points.
        service: Service distribution, as analyze_queue takes it.
        service_rate: Exponential service rate (mu).
        arrival: Arrival distribution of an open class.
        arrival_rate: Poisson arrival rate (lambda).
        population: Closed-class population (N).
        think_time: Mean think time of the closed class (Z).
        classes: JSON list of class objects, for a multiclass or mixed base
            model, as analyze_queue takes it.
        servers: Number of servers (k). Default 1.
        queue_capacity: Total capacity including servers. -1 is infinite.
        scheduling: Scheduling discipline by name. Default FCFS.
        solver: Solver to run at every point. Default auto.
        method: Algorithm within that solver.
        seed: Random seed of the stochastic solvers.
        samples: Simulation samples / Monte Carlo draws.
        format: "text" (default, one table) or "json" (one record per station,
            class and swept value).

    Returns:
        One row per station and class at each swept value, with the swept
        parameter as the leading column.
    """
    import pandas as pd

    try:
        vals = _parse_values(values)
    except ValueError as e:
        return _error("Error parsing values: %s" % e, format)
    if not vals:
        return _error("Error: no values to sweep.", format)
    if len(vals) > 100:
        return _error("Error: maximum 100 sweep points allowed.", format)

    # A parameter is sweepable when it is a NUMBER: the distributions and the
    # discipline are swept by giving several base models, not by interpolating.
    sweepable = {"arrival_rate", "service_rate", "servers", "population",
                 "think_time", "queue_capacity"}
    if param not in sweepable:
        return _error("Unknown parameter '%s'. Use one of: %s"
                      % (param, ", ".join(sorted(sweepable))), format)

    base = {"service": service, "service_rate": service_rate,
            "arrival": arrival, "arrival_rate": arrival_rate,
            "population": population, "think_time": think_time,
            "classes": classes, "servers": servers,
            "queue_capacity": queue_capacity, "scheduling": scheduling,
            "solver": solver, "method": method, "seed": seed,
            "samples": samples, "format": "json"}
    integral = {"servers", "population", "queue_capacity"}

    rows: list = []
    errors: list = []
    for v in vals:
        kwargs = dict(base)
        kwargs[param] = int(v) if param in integral else v
        try:
            data = json.loads(analyze_queue(**kwargs))
        except Exception as exc:
            errors.append("%s=%s: %s" % (param, v, exc))
            continue
        if not isinstance(data, dict) or "results" not in data:
            errors.append("%s=%s: %s" % (param, v, data.get("error", data)
                                         if isinstance(data, dict) else data))
            continue
        records = data["results"]
        if not isinstance(records, list):
            errors.append("%s=%s: the analysis returned no table" % (param, v))
            continue
        for rec in records:
            rec[param] = v
            rows.append(rec)

    if not rows:
        return _error("No successful runs.\n" + "\n".join(errors), format)

    if format == "json":
        payload = {"results": rows,
                   "metadata": {"param": param, "values": vals,
                                "solver": solver.upper(), "points": len(vals)}}
        if errors:
            payload["metadata"]["errors"] = errors
        return json.dumps(_sanitize_json(payload), indent=2, default=str)

    df = pd.DataFrame(rows)
    if param in df.columns:
        df = df[[param] + [c for c in df.columns if c != param]]
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
    # Without these four a Petri net drew as a row of blue queue rectangles:
    # the .get() fallback below is a queue, so an unmapped type was not merely
    # uncoloured, it was mislabelled as a different kind of node.
    9: ("(({label}))", "place"),       # PLACE — circle, as a Petri net draws one
    10: ("[/{label}/]", "transition"),  # TRANSITION — parallelogram (a bar)
    # Rounded rather than a trapezoid on purpose: the trapezoid forms are
    # delimited by backslashes, and every label here is quoted, so `[\"x"\]`
    # would be a new way to lose the whole diagram to a parse error. Colour
    # already tells a logger from a delay.
    11: ("({label})", "logger"),       # LOGGER — rounded
    12: ("[[{label}]]", "region"),     # FINITE_CAPACITY_REGION — subroutine box
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
    "place": "#3F51B5",
    "transition": "#212121",
    "logger": "#9E9E9E",
    "region": "#8BC34A",
}


def _mermaid_label(text) -> str:
    """Quote a node label so Mermaid reads it verbatim.

    A bare label may not contain (), [] or {}, which is exactly what the
    "(IS)" and "(k=N)" server suffixes below append: the unquoted form is a
    parse error, not a cosmetic problem. Quoting is accepted in every shape,
    so it is applied unconditionally rather than sniffed for.
    """
    return '"' + str(text).replace('"', "&quot;") + '"'


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

        node_def = shape_tpl.format(label=_mermaid_label(label))
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


@_tool
def visualize_model(code: str = "", model_id: str = "") -> str:
    """Generate a Mermaid diagram of a queueing network.

    Draws whichever model it is given: one built in conversation (model_id) or
    one described by LINE code.

    Args:
        code: LINE Python code that creates a Network and calls model.link().
        model_id: A model built by build_network, build_topology or
            build_from_gallery, drawn instead of running any code.

    Returns:
        Mermaid diagram source (paste into any Mermaid renderer).
    """
    if model_id.strip():
        try:
            model = _resolve_model(model_id.strip())
        except Exception as e:
            return "Error: %s" % e
    elif code.strip():
        try:
            # _model_from_code resolves a LayeredNetwork and an Environment too,
            # where the loop this used to run saw only a Network.
            model = _model_from_code(code)
        except Exception as e:
            return "Error: %s" % e
    else:
        return ("Error: no model given: pass `model_id` for a model built here, "
                "or `code` that builds one.")

    try:
        sn = model.get_struct()
        return _generate_mermaid(sn)
    except Exception:
        return f"Error generating diagram:\n{traceback.format_exc()}"


# ---------------------------------------------------------------------------
# 8. Multi-solver comparison
# ---------------------------------------------------------------------------

@_tool
def compare_solvers(code: str = "", solvers: str = "", model_id: str = "") -> str:
    """Run the same model with multiple solvers and compare results.

    Solves one model with each requested solver and returns a merged table with
    a Solver column. The model is either one built in conversation (model_id) or
    one described by LINE code.

    Args:
        code: LINE Python code that creates a Network and calls model.link().
        solvers: Comma-separated solver names (default: "MVA,NC,CTMC,SSA,FLD,MAM").
        model_id: A model built by build_network, build_topology or
            build_from_gallery, compared instead of running any code.

    Returns:
        Combined performance metrics table from all solvers.
    """
    import pandas as pd

    if not solvers.strip():
        solver_list = ["MVA", "NC", "CTMC", "SSA", "FLD", "MAM"]
    else:
        solver_list = [s.strip().upper() for s in solvers.split(",") if s.strip()]

    if model_id.strip():
        try:
            model = _resolve_model(model_id.strip())
        except Exception as e:
            return "Error: %s" % e
    elif code.strip():
        try:
            model = _model_from_code(code)
        except Exception as e:
            return "Error: %s" % e
    else:
        return ("Error: no model given: pass `model_id` for a model built here, "
                "or `code` that builds one.")

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

@_tool
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
# 10. The line-cli surface: model files, every solver, every analysis
# ---------------------------------------------------------------------------
# The tools above build a model out of Python code. A command line does not:
# it reads a DOCUMENT, picks a solver by method name, asks for one of forty
# analyses and renders the answer in a chosen format. That whole surface was
# unreachable from this server, so a client holding a `model.json`, a `.lqnx`
# or a JMT `.jsimg` had to paste it into code before it could be solved, and
# every analysis other than the average table had to be reached by writing the
# getter call by hand.
#
# Everything below WRAPS ``line_solver.cli`` -- the native twin of
# ``jline.cli.LineCLI`` and of the C++ ``common/line-cli`` -- rather than
# restating its tables here. The vocabulary a client sees is therefore the
# CLI's own: a solver method name or an analysis added to that module reaches this
# server with no edit, and the two cannot drift into two opinions about what
# LINE can be asked. Keep it that way -- a table copied here is a table that
# will lag.


def _sniff_input_format(content: str) -> str:
    """The format of a model handed over inline, read off its first bytes.

    The CLI reads stdin as ``jsim`` unless told otherwise, which is a JMT
    default a conversational caller rarely means; here the document says what
    it is -- a brace opens LINE's own JSON and the three XML dialects are told
    apart by their root element. An unrecognised document is called json, so
    the failure is the json parser's own message about the caller's text rather
    than a guess reported as a fact. Pass ``input_format`` to settle it.
    """
    head = content.lstrip()[:512].lower()
    if head.startswith("{") or head.startswith("["):
        return "json"
    if "<pnml" in head:
        return "pnml"
    if "<lqn" in head:
        return "lqnx"
    if "<sim" in head or "<archive" in head or "<jsim" in head:
        return "jsimg"
    return "json"


def _resolve_input_format(file: str, content: str, input_format: str) -> str:
    """The format of the model about to be read, WITHOUT reading it.

    Resolved on its own so the per-format solver whitelist can refuse by name
    before the document is parsed: the format follows from the extension, the
    first bytes or the caller's own word, none of which needs the model. Asking
    afterwards, as the CLI does, answers a bad pairing with the parser's
    complaint about the file instead of with the list of solvers that format
    admits.
    """
    from . import cli as line_cli

    fmt = (input_format or "").strip().lower() or None
    if fmt is not None and fmt not in line_cli.SUPPORTED_FORMATS:
        raise ValueError(
            "unknown input format '%s'; supported: %s"
            % (fmt, ", ".join(line_cli.SUPPORTED_FORMATS)))
    if file and content.strip():
        raise ValueError("give the model either as `file` (a path) or as `content` "
                         "(the document itself), not both")
    if file:
        return fmt or line_cli.detect_input_format(file)
    if not content.strip():
        raise ValueError("no model given: pass `file` (a path) or `content` (the document itself)")
    return fmt or _sniff_input_format(content)


def _load_cli_model(file: str, content: str, fmt: str, verbose: bool):
    """Load the model named by a path or given inline, as ``-f`` and stdin do."""
    from . import cli as line_cli

    if file:
        return line_cli.load_model(file, fmt, verbose)
    # `mat` and `pkl` are byte formats and a JSON-RPC string cannot carry them.
    if fmt in ("mat", "pkl"):
        raise ValueError("the '%s' format is binary and must be given as `file`, not inline" % fmt)
    with tempfile.NamedTemporaryFile(mode="w", suffix="." + fmt, delete=False) as fh:
        fh.write(content)
        path = fh.name
    try:
        return line_cli.load_model(path, fmt, verbose)
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass


@contextlib.contextmanager
def _session_verbosity(level: str):
    """Hold ``GlobalConstants.Verbose`` at *level* for one tool call only.

    The CLI sets it once per process and exits. This server outlives any one
    call, so a tool asking for ``debug`` must not leave every later tool
    narrating: the previous level is restored whatever the run does.
    """
    from line_solver import GlobalConstants, VerboseLevel

    wanted = {"silent": VerboseLevel.SILENT,
              "normal": VerboseLevel.STD,
              "debug": VerboseLevel.DEBUG}[level]
    saved = GlobalConstants.getVerbose()
    GlobalConstants.setVerbose(wanted)
    try:
        yield
    finally:
        GlobalConstants.setVerbose(saved)


def _resolve_verbosity(verbosity: str) -> str:
    """`-v` as one of the three levels. `standard` is the C++ CLI's `normal`."""
    level = "normal" if verbosity == "standard" else (verbosity or "silent")
    if level not in ("silent", "normal", "debug"):
        raise ValueError("verbosity is silent, normal (alias standard) or debug "
                         "(got '%s')" % verbosity)
    return level


def _solve_and_format(
    get_model,
    solver: str = "auto",
    analysis: str = "avg",
    output: str = "readable",
    method: str = "",
    node: int | None = None,
    class_idx: int | None = None,
    state: str = "",
    events: int = 1000,
    percentiles: str = "50,90,95,99",
    reward_name: str = "",
    seed: int | None = None,
    samples: int | None = None,
    cutoff: float | None = None,
    timespan: str = "",
    timestep: float | None = None,
    tol: float | None = None,
    iter_tol: float | None = None,
    iter_max: int | None = None,
    multiserver: str = "",
    fork_join: str = "",
    warmupfrac: float | None = None,
    stage_solver: str = "",
    uq_solver: str = "",
    busyperiod: str = "",
    busyperiod_subnet: str = "",
    sens_method: str = "auto",
    sens_scheme: str = "forward",
    sens_step: float | None = None,
    verbosity: str = "silent",
    options: str = "",
    output_file: str = "",
) -> str:
    """The whole `line-cli` solve, for a model from wherever it came.

    Everything past "I have a model" is here: the argument validation, the
    solver token, the knobs, the analysis dispatch and the rendering. It is
    deliberately the ONLY place any of that happens, so a model built in
    conversation reaches exactly the analyses a model read from a file does,
    rather than a second and smaller vocabulary that drifts from the first.

    Args:
        get_model: Callable taking the resolved `verbose` flag and returning
            the model. It runs inside the verbosity context and inside the
            error handler, so a source that has its own validation (a format
            whitelist, a spec check) does it in here and fails like any other
            usage error.
        Everything else: as documented on solve_model_file.

    Returns:
        The analysis rendered in the requested output format.
    """
    from . import cli as line_cli

    err_format = "json" if output == "json" else "text"
    level = "silent"
    try:
        level = _resolve_verbosity(verbosity)

        if output not in line_cli.OUTPUT_FORMATS:
            raise ValueError("unknown output format '%s'; supported: %s"
                             % (output, ", ".join(line_cli.OUTPUT_FORMATS)))

        results, model, _ = _run_cli_analysis(
            get_model, level,
            solver=solver, analysis=analysis, method=method, node=node,
            class_idx=class_idx, state=state, events=events,
            percentiles=percentiles, reward_name=reward_name, seed=seed,
            samples=samples, cutoff=cutoff, timespan=timespan,
            timestep=timestep, tol=tol, iter_tol=iter_tol, iter_max=iter_max,
            multiserver=multiserver, fork_join=fork_join,
            warmupfrac=warmupfrac, stage_solver=stage_solver,
            uq_solver=uq_solver, busyperiod=busyperiod,
            busyperiod_subnet=busyperiod_subnet, sens_method=sens_method,
            sens_scheme=sens_scheme, sens_step=sens_step, options=options)
        formatted = line_cli.format_results(results, output, model,
                                            output_file or None)
    except Exception as e:
        message = "Error: %s" % e
        if level == "debug":
            message += "\n" + traceback.format_exc()
        return _error(message, err_format)

    if isinstance(formatted, bytes):
        formatted = formatted.hex()
    return formatted if formatted.strip() else "(no results: no analysis produced output)"


def _run_cli_analysis(
    get_model,
    level: str = "silent",
    solver: str = "auto",
    analysis: str = "avg",
    method: str = "",
    node: int | None = None,
    class_idx: int | None = None,
    state: str = "",
    events: int = 1000,
    percentiles: str = "50,90,95,99",
    reward_name: str = "",
    seed: int | None = None,
    samples: int | None = None,
    cutoff: float | None = None,
    timespan: str = "",
    timestep: float | None = None,
    tol: float | None = None,
    iter_tol: float | None = None,
    iter_max: int | None = None,
    multiserver: str = "",
    fork_join: str = "",
    warmupfrac: float | None = None,
    stage_solver: str = "",
    uq_solver: str = "",
    busyperiod: str = "",
    busyperiod_subnet: str = "",
    sens_method: str = "auto",
    sens_scheme: str = "forward",
    sens_step: float | None = None,
    options: str = "",
):
    """The `line-cli` solve itself: validation, solver token, knobs, dispatch.

    Returns the results, the model and the solver, RAISING the way the CLI does
    rather than rendering an error. Two callers render it differently and
    neither may grow a second copy of this: `_solve_and_format` hands it to
    `line_cli.format_results` for the document surface, and the analyze_*
    family wraps it in the metadata header those tools have always carried.

    Args:
        get_model: Callable taking the resolved `verbose` flag and returning
            the model, as documented on `_solve_and_format`.
        level: Resolved verbosity, already through `_resolve_verbosity`.
        Everything else: as documented on solve_model_file.

    Returns:
        (results, model, the solver object that produced them).
    """
    from . import cli as line_cli

    verbose = level in ("normal", "debug")

    # Parsed before anything is loaded, so a malformed list is a usage error
    # and not a failure half way through a solve.
    timespan_v = line_cli.parse_timespan(timespan or None)
    busy_orders = line_cli.parse_int_list(busyperiod or None, "busyperiod")
    busy_subnet = line_cli.parse_int_list(busyperiod_subnet or None, "busyperiod_subnet")

    analysis_types = line_cli.validate_analysis_types(analysis)
    line_cli.validate_analysis_solver_compat(analysis_types, solver)
    line_cli.validate_analysis_params(analysis_types, node, class_idx, reward_name or None)

    extra: dict = {}
    if options.strip():
        extra = json.loads(options)
        if not isinstance(extra, dict):
            raise ValueError('`options` is a JSON object of solver options, '
                             'e.g. {"stiff": false}')

    with _session_verbosity(level):
        model = get_model(verbose)

        # The composite families name their inner solver in the method name and
        # not in a separate option, so --uq-solver and --stage-solver fold
        # into the -s token here exactly as they do in the native CLI.
        solver_token = line_cli.resolve_solver_token(
            solver, method or None, uq_solver or None, stage_solver or None)

        knobs = {
            "samples": samples,
            "cutoff": cutoff,
            "timespan": timespan_v,
            "timestep": timestep,
            "tol": tol,
            "iter_tol": iter_tol,
            "iter_max": iter_max,
            "multiserver": multiserver or None,
            "fork_join": fork_join or None,
            "warmupfrac": warmupfrac,
        }
        knobs.update(extra)

        solver_obj = line_cli.solve_model(model, solver_token, seed, verbose, knobs)
        results = line_cli.get_solver_results(
            solver_obj,
            analysis,
            model=model,
            node_idx=node,
            class_idx=class_idx,
            state_str=state or None,
            num_events=events,
            percentiles_str=percentiles,
            reward_name=reward_name or None,
            sens_method=sens_method,
            sens_scheme=sens_scheme,
            sens_step=sens_step,
            busy_subnet=busy_subnet,
            busy_orders=busy_orders,
        )

    return results, model, solver_obj


@_tool
def solve_model_file(
    file: str = "",
    content: str = "",
    model_id: str = "",
    input_format: str = "",
    solver: str = "auto",
    analysis: str = "avg",
    output: str = "readable",
    method: str = "",
    node: int | None = None,
    class_idx: int | None = None,
    state: str = "",
    events: int = 1000,
    percentiles: str = "50,90,95,99",
    reward_name: str = "",
    seed: int | None = None,
    samples: int | None = None,
    cutoff: float | None = None,
    timespan: str = "",
    timestep: float | None = None,
    tol: float | None = None,
    iter_tol: float | None = None,
    iter_max: int | None = None,
    multiserver: str = "",
    fork_join: str = "",
    warmupfrac: float | None = None,
    stage_solver: str = "",
    uq_solver: str = "",
    busyperiod: str = "",
    busyperiod_subnet: str = "",
    sens_method: str = "auto",
    sens_scheme: str = "forward",
    sens_step: float | None = None,
    verbosity: str = "silent",
    options: str = "",
    output_file: str = "",
) -> str:
    """Solve a model DOCUMENT, the whole `line-cli` command line as one tool.

    This is the command-line front end: it reads a model file in any format
    LINE understands, runs any solver, asks for any of the analyses the library
    publishes, and renders the answer. Use it whenever the model already exists
    as a file or as text; use solve_line_model when the model is to be BUILT in
    Python.

    Args:
        file: Path to the model document. Mutually exclusive with content
            and model_id.
        content: The document itself, when there is no file to point at.
        model_id: A model built this session by build_network, build_topology
            or build_from_gallery. Given this, no document is read and every
            analysis below applies to the built model exactly as it would to a
            file.
        input_format: json (LINE's own interchange model), jsim/jsimg/jsimw
            (JMT), lqnx/xml (layered), pnml (place/transition net), mat, pkl.
            Detected from the extension of file, or from the first bytes of
            content, when not stated.
        solver: Solver or solver.method. Network: auto, mva, nc, ctmc, mam,
            fld (alias fluid), ssa, ldes, jmt, lqns (qns methods), ag, ba, uq. Layered: ln,
            ln.mva, ln.nc, ln.comom, lqns. Random environment: env. The
            high-level names auto, sim, exact, fast and accurate pick for you.
        analysis: One analysis or a comma-separated list. avg, sys, chain,
            node, nodechain, stage, all; cache, item; orbit, loss, region-loss,
            deadline; normconst, busyperiod; cdf-respt, cdf-passt, perct-respt;
            tran-avg, tran-cdf-respt, tran-cdf-passt; prob, prob-aggr,
            prob-marg, prob-sys, prob-sys-aggr, prob-sys-marg; sample,
            sample-aggr, sample-sys, sample-sys-aggr; reward, reward-steady,
            reward-value; sens; interval. Several need a particular solver --
            list_solvers reports which.
        output: readable (default), json, csv, or pickle/mat with output_file.
        method: Algorithm within the chosen solver, also spelled solver.method.
        node: Node index for the prob and sample analyses (0-based).
        class_idx: Job class index for prob-marg (0-based).
        state: State vector for prob, comma-separated integers.
        events: Events drawn by the sample analyses. Default 1000.
        percentiles: Percentiles of perct-respt. Default "50,90,95,99".
        reward_name: Built-in reward of reward-value: QLen, Tput, Util, RespT,
            WaitT, ArvR, ResidT.
        seed: Random seed of the stochastic solvers (JMT, SSA, LDES).
        samples: Simulation samples / Monte Carlo draws (SSA, LDES, JMT, UQ).
        cutoff: State-space cutoff per open class (CTMC, SSA).
        timespan: Horizon of the transient analyses: "T1", "T0,T1" or "T0:T1".
        timestep: Fixed transient output step, instead of the adaptive grid.
        tol: General solver tolerance.
        iter_tol: Iteration convergence tolerance.
        iter_max: Maximum iterations.
        multiserver: AMVA multiserver rule (seidmann, softmin, ...).
        fork_join: Fork-join method: default/ht (Heidelberger-Trivedi), mmt.
        warmupfrac: Leading fraction of a simulated path discarded before the
            means are taken, in [0,1).
        stage_solver: Solver run at each stage of an Environment (-s env).
        uq_solver: Engine SolverUQ runs at each design point; -s uq needs it.
        busyperiod: Orders of -a busyperiod, comma-separated. Default 1.
        busyperiod_subnet: 0-based stations forming the busyperiod subnetwork.
        sens_method: Differentiation used by -a sens. Default auto.
        sens_scheme: Finite-difference scheme of -a sens. Default forward.
        sens_step: Finite-difference step of -a sens.
        verbosity: silent (default here), normal (alias standard) or debug.
            debug turns on the solver console, a running progress log, and adds
            the traceback to an error.
        options: JSON object of further solver options, forwarded verbatim the
            way jline.cli.LineCLI forwards a stated option into SolverOptions,
            e.g. {"stiff": false, "keep": true}.
        output_file: Destination of the pickle and mat output formats.

    Returns:
        The analysis rendered in the requested output format.
    """
    from . import cli as line_cli

    def _get(verbose):
        if model_id.strip():
            if file or content.strip():
                raise ValueError("give the model as `model_id` OR as a document "
                                 "(`file`/`content`), not both")
            # A built model has no input format, so there is no per-format
            # solver whitelist to apply: the model object is what the solver
            # will see either way.
            return _resolve_model(model_id.strip())
        # The format, and therefore the whitelist, is settled before the
        # document is read: a solver that format does not admit is refused by
        # name rather than after a parse that was never going to be used.
        fmt = _resolve_input_format(file, content, input_format)
        line_cli.validate_solver_compatibility(fmt, solver)
        return _load_cli_model(file, content, fmt, verbose)

    return _solve_and_format(
        _get,
        solver=solver, analysis=analysis, output=output, method=method,
        node=node, class_idx=class_idx, state=state, events=events,
        percentiles=percentiles, reward_name=reward_name, seed=seed,
        samples=samples, cutoff=cutoff, timespan=timespan, timestep=timestep,
        tol=tol, iter_tol=iter_tol, iter_max=iter_max, multiserver=multiserver,
        fork_join=fork_join, warmupfrac=warmupfrac, stage_solver=stage_solver,
        uq_solver=uq_solver, busyperiod=busyperiod,
        busyperiod_subnet=busyperiod_subnet, sens_method=sens_method,
        sens_scheme=sens_scheme, sens_step=sens_step, verbosity=verbosity,
        options=options, output_file=output_file,
    )


@_tool
def find_solver(
    file: str = "",
    content: str = "",
    input_format: str = "",
    code: str = "",
    model_id: str = "",
    metric: str = "",
    show_all: bool = False,
    format: str = "text",
) -> str:
    """Which solvers and methods can analyze this model (`line-cli --find-solver`).

    Reports rather than solves: one row per (solver, method) pair, saying
    whether it runs, whether it is exact, approximate, a bound or a simulation,
    and which measures it answers. Ask it before choosing a solver, and after a
    solver refuses a model.

    Args:
        file: Path to the model document.
        content: The document itself, when there is no file to point at.
        input_format: As in solve_model_file; detected when not stated.
        code: LINE Python code building the model instead of reading one. The
            first Network, LayeredNetwork or Environment it defines is the one
            reported on.
        model_id: A model built this session by build_network, build_topology
            or build_from_gallery, reported on instead of a document.
        metric: Narrow the report to a measure group ('cdf', 'tran', 'prob',
            'sample') or to an accessor ('getCdfRespT'). Empty keeps every pair.
        show_all: Also list the pairs that are REFUSED, with the reason each
            was refused. This is `--find-solver-all`.
        format: "text" (default) or "json".

    Returns:
        The findSolver table, whose Method column is the method name to pass back
        as solve_model_file's `solver`.
    """
    try:
        if model_id.strip():
            model = _resolve_model(model_id.strip())
        elif code.strip():
            model = _model_from_code(code)
        else:
            model = _load_cli_model(file, content,
                                    _resolve_input_format(file, content, input_format), False)
        table = model.findSolver(metric, show_all)
    except Exception as e:
        return _error("Error: %s" % e, format)

    if format == "json":
        return json.dumps({"results": table.to_dict(orient="records")}, indent=2, default=str)
    return table.to_string(index=False)


def _model_from_code(code: str):
    """The first model object defined by *code*, for the report tools.

    A Network, a LayeredNetwork or an Environment all answer findSolver, and
    which one the code builds is the caller's business, so all three are
    accepted rather than the Network the older tools look for.
    """
    from line_solver import Network
    from line_solver.environment import Environment
    from line_solver.layered import LayeredNetwork

    ns = _make_exec_namespace()
    buffer = io.StringIO()
    saved = sys.stdout
    sys.stdout = buffer
    try:
        _exec_with_timeout(code, ns)
    finally:
        sys.stdout = saved
    for value in ns.values():
        if isinstance(value, (Network, LayeredNetwork, Environment)):
            return value
    raise ValueError("no Network, LayeredNetwork or Environment was defined by the code")


@_tool
def list_solvers(format: str = "text") -> str:
    """The vocabulary: every solver, analysis and format this server accepts.

    The same tables `line-cli --help-all` prints, read off the CLI module so
    they cannot drift: the solver method names and their aliases, the analyses,
    which analyses require which solver, which solvers each input format
    admits, and the output formats.

    Args:
        format: "text" (default) or "json".

    Returns:
        The vocabulary of solve_model_file's `solver`, `analysis`, `input_format`
        and `output` arguments.
    """
    from . import cli as line_cli

    payload = {
        "solvers": list(line_cli.VALID_SOLVERS),
        "solver_aliases": dict(line_cli.SOLVER_ALIASES),
        "solver_method_examples": list(line_cli.SOLVER_METHOD_EXAMPLES),
        "analyses": list(line_cli.VALID_ANALYSIS_TYPES),
        "analysis_requires_solver": {k: list(v)
                                     for k, v in line_cli.ANALYSIS_SOLVER_COMPAT.items()},
        "analysis_requires_node": list(line_cli.ANALYSIS_REQUIRES_NODE),
        "analysis_requires_class": list(line_cli.ANALYSIS_REQUIRES_CLASS),
        "reward_names": list(line_cli.VALID_REWARD_NAMES),
        "input_formats": list(line_cli.SUPPORTED_FORMATS),
        "output_formats": list(line_cli.OUTPUT_FORMATS),
        "solvers_by_input_format": {
            "json": list(line_cli.JSON_COMPATIBLE_SOLVERS),
            "jsim": list(line_cli.JSIM_COMPATIBLE_SOLVERS),
            "lqnx": list(line_cli.LQN_COMPATIBLE_SOLVERS),
            "pnml": list(line_cli.PNML_COMPATIBLE_SOLVERS),
        },
    }

    if format == "json":
        return json.dumps(payload, indent=2, default=str)

    lines = ["LINE solver vocabulary", "=" * 40, ""]
    lines.append("## Solvers")
    lines.append("  " + ", ".join(payload["solvers"]))
    lines.append("  aliases: " + ", ".join("%s -> %s" % kv
                                           for kv in payload["solver_aliases"].items()))
    lines.append("  solver.method examples: " + ", ".join(payload["solver_method_examples"]))
    lines.append("")
    lines.append("## Analyses")
    lines.append("  " + ", ".join(payload["analyses"]))
    lines.append("")
    lines.append("## Analyses that require a particular solver")
    for name, solvers in payload["analysis_requires_solver"].items():
        lines.append("  %-16s %s" % (name, " or ".join(solvers)))
    lines.append("")
    lines.append("## Analyses that require an index")
    lines.append("  node:  " + ", ".join(payload["analysis_requires_node"]))
    lines.append("  class: " + ", ".join(payload["analysis_requires_class"]))
    lines.append("")
    lines.append("## Reward names (reward-value)")
    lines.append("  " + ", ".join(payload["reward_names"]))
    lines.append("")
    lines.append("## Input formats")
    lines.append("  " + ", ".join(payload["input_formats"]))
    for fmt, solvers in payload["solvers_by_input_format"].items():
        lines.append("  %-6s %s" % (fmt, ", ".join(solvers)))
    lines.append("")
    lines.append("## Output formats")
    lines.append("  " + ", ".join(payload["output_formats"]))
    return "\n".join(lines)


@_tool
def check_install(format: str = "text") -> str:
    """Environment check: which optional backends are reachable (`line-cli --install`).

    Reports on the Java runtime the JMT and LDES solvers need and on the
    SageMath symbolic backend the symbolic CTMC and fluid methods start. Never
    fails: a missing dependency costs a subset of the solvers and leaves the
    rest usable.

    Args:
        format: "text" (default) or "json".

    Returns:
        The check's own report, plus each warning it raised.
    """
    import warnings as _warnings

    from .install import line_install

    buffer = io.StringIO()
    saved = sys.stdout
    sys.stdout = buffer
    try:
        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            ok = line_install()
    except Exception as e:
        sys.stdout = saved
        return _error("Error: %s" % e, format)
    finally:
        sys.stdout = saved

    messages = [str(w.message) for w in caught]
    if format == "json":
        return json.dumps({"ok": bool(ok), "log": buffer.getvalue().splitlines(),
                           "warnings": messages}, indent=2, default=str)
    report = buffer.getvalue().rstrip()
    if messages:
        report += "\n\nWarnings:\n" + "\n".join("  - " + m for m in messages)
    report += "\n\n%s" % ("All optional dependencies are in place."
                          if ok else "Some optional dependencies are missing (see above).")
    return report


@_tool
def get_version(format: str = "text") -> str:
    """LINE's version and where this server is reading its parts from.

    The `line-cli -V` report, plus the interpreter and the examples directory,
    which is what a client needs when list_examples comes back empty.

    Args:
        format: "text" (default) or "json".
    """
    from line_solver import GlobalConstants

    payload = {
        "line_version": GlobalConstants.Version,
        "python": sys.version.split()[0],
        "package": os.path.dirname(os.path.abspath(__file__)),
        "examples_dir": EXAMPLES_DIR,
        "examples_dir_exists": os.path.isdir(EXAMPLES_DIR),
    }
    if format == "json":
        return json.dumps(payload, indent=2, default=str)
    return "\n".join("%-20s %s" % (k, v) for k, v in payload.items())

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Run the server over stdio. Entry point of the ``line-mcp`` script."""
    mcp.run()


if __name__ == "__main__":
    main()
