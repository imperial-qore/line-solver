#!/usr/bin/env python3
"""LINE CLI - Command-line interface for LINE queueing network solver.

A single-file script for solving queueing network models using the LINE solver.

Usage:
    python line-cli.py solve model.jsimg -s mva
    python line-cli.py solve model.jsimg -s nc -o json
    python line-cli.py solve model.lqnx -s ln
    python line-cli.py info
    python line-cli.py list solvers
    python line-cli.py server -p 5863     # Start WebSocket server
    python line-cli.py rest -p 8080       # Same server, port 8080 (see cmd_rest)
"""

import argparse
import csv
import io
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple

__version__ = "3.0.7"

# =============================================================================
# Data Models
# =============================================================================

@dataclass
class MetricRow:
    """A single row of metrics from solver output."""
    station: str
    job_class: str
    metric: str
    value: float
    unit: str = ""


@dataclass
class SolveResult:
    """Result of a solver execution."""
    success: bool
    avg_table: Optional[List[MetricRow]] = None
    sys_table: Optional[List[MetricRow]] = None
    raw_output: str = ""
    error_message: Optional[str] = None
    execution_time: float = 0.0


@dataclass
class SolveOptions:
    """Options for solver execution."""
    solver: str = "auto"
    input_format: str = "jsim"
    output_format: str = "readable"
    analysis: str = "all"
    seed: Optional[int] = None
    verbose: bool = False
    # Extended analysis parameters
    node_index: Optional[int] = None
    class_index: Optional[int] = None
    state: Optional[str] = None
    num_events: int = 1000
    percentiles: str = "50,90,95,99"
    reward_name: Optional[str] = None
    # Numeric controls. NONE of these could be expressed through this wrapper,
    # although the JAR has taken every one of them for as long as it has had an
    # SSA branch: a simulation run from here always used the default 10000
    # samples, a CTMC could not be given a cutoff, and a transient analysis
    # answered over a horizon nobody chose.
    samples: Optional[int] = None
    cutoff: Optional[float] = None
    timespan: Optional[str] = None
    timestep: Optional[float] = None
    method: Optional[str] = None
    tol: Optional[float] = None
    iter_tol: Optional[float] = None
    iter_max: Optional[int] = None
    multiserver: Optional[str] = None
    warmupfrac: Optional[float] = None
    stage_solver: Optional[str] = None
    uq_solver: Optional[str] = None
    busyperiod: Optional[str] = None
    busyperiod_subnet: Optional[str] = None
    sens_method: Optional[str] = None
    sens_scheme: Optional[str] = None
    sens_step: Optional[float] = None


# =============================================================================
# Solver and Format Definitions
# =============================================================================

SOLVERS: Dict[str, Dict[str, Any]] = {
    "auto": {
        "name": "Automatic Solver Selection",
        "description": "Automatically selects the best solver for the model",
        "formats": ["jsim", "jsimg", "jsimw", "lqnx", "xml", "json", "pnml"],
    },
    "ctmc": {
        "name": "Continuous-Time Markov Chain",
        "description": "Exact analysis using CTMC state space exploration",
        "formats": ["jsim", "jsimg", "jsimw", "json", "pnml"],
    },
    # THE JAR'S TOKEN IS `ldes`. This entry used to be keyed `des` alone and
    # the wrapper forwarded the key VERBATIM as `-s des`, which the JAR rejects
    # -- so the SSJ engine was unreachable from this CLI by any spelling, while
    # `list solvers` and CLI.md both advertised it. `des` survives as an alias.
    "ldes": {
        "name": "Discrete Event Simulation (LDES)",
        "description": "Discrete event simulation using the SSJ-based LDES engine",
        "formats": ["jsim", "jsimg", "jsimw", "json", "pnml"],
    },
    "ag": {
        "name": "Agent-Based (RCAT/INAP)",
        "description": "Product-form reversed-rate analysis of interacting agents",
        "formats": ["jsim", "jsimg", "jsimw", "json"],
    },
    "ba": {
        "name": "Bound Analysis",
        "description": "Closed-form bounds on throughput and response time",
        "formats": ["jsim", "jsimg", "jsimw", "json"],
    },
    "env": {
        "name": "Random Environment",
        "description": "Blended analysis of a model in a random environment",
        "formats": ["json"],
    },
    "uq": {
        "name": "Uncertainty Quantification",
        "description": "Prior-weighted analysis over a design of models (needs --uq-solver)",
        "formats": ["jsim", "jsimg", "jsimw", "json"],
    },
    "fld": {
        "name": "Fluid/Mean-Field ODE",
        "description": "Approximate analysis using fluid/ODE model",
        "formats": ["jsim", "jsimg", "jsimw", "json"],
    },
    "jmt": {
        "name": "Java Modelling Tools",
        "description": "Discrete event simulation using JMT",
        "formats": ["jsim", "jsimg", "jsimw", "json", "pnml"],
    },
    "ln": {
        "name": "Layered Network",
        "description": "Solver for layered queueing networks (MVA layers)",
        "formats": ["lqnx", "xml", "json"],
    },
    "ln.mva": {
        "name": "Layered Network (MVA layers)",
        "description": "Layered solver with the layer fixed point solved by MVA",
        "formats": ["lqnx", "xml", "json"],
    },
    "ln.nc": {
        "name": "Layered Network (NC layers)",
        "description": "Layered solver with the layer fixed point solved by NC",
        "formats": ["lqnx", "xml", "json"],
    },
    "ln.comom": {
        "name": "Layered Network (CoMoM layers)",
        "description": "Layered solver with NC/CoMoM layers",
        "formats": ["lqnx", "xml", "json"],
    },
    "lqns": {
        "name": "LQN Solver",
        "description": "External LQNS solver integration",
        "formats": ["lqnx", "xml", "json"],
    },
    "mam": {
        "name": "Matrix Analytic Methods",
        "description": "Analysis using matrix analytic methods (supports Fork-Join percentiles)",
        "formats": ["jsim", "jsimg", "jsimw", "json"],
    },
    "mva": {
        "name": "Mean Value Analysis",
        "description": "Analytical solver using Mean Value Analysis algorithm",
        "formats": ["jsim", "jsimg", "jsimw", "json"],
    },
    "nc": {
        "name": "Normalizing Constant",
        "description": "Exact analysis using normalizing constant computation",
        "formats": ["jsim", "jsimg", "jsimw", "json"],
    },
    "qns": {
        "name": "QNS",
        "description": "External QNSolver integration",
        "formats": ["jsim", "jsimg", "jsimw", "json"],
    },
    "ssa": {
        "name": "Stochastic Simulation Algorithm",
        "description": "Stochastic simulation of the model",
        "formats": ["jsim", "jsimg", "jsimw", "json", "pnml"],
    },
}

# Solver aliases (map alias -> canonical name)
SOLVER_ALIASES: Dict[str, str] = {
    "fluid": "fld",
    "qnsolver": "qns",
    # `des` was this wrapper's ONLY spelling of the LDES engine and was
    # forwarded verbatim to a JAR that knows it as `ldes`. Kept as an alias so
    # existing command lines still parse, now resolving to a token that works.
    "des": "ldes",
}


def resolve_solver(solver: str) -> str:
    """Resolve solver alias to canonical name."""
    return SOLVER_ALIASES.get(solver.lower(), solver.lower())


def auto_select_solver(input_format: str, model_file: Optional[str] = None) -> str:
    """Select an appropriate solver based on input format when 'auto' is specified.

    Returns a concrete solver name since the JAR doesn't support 'auto'.
    """
    # For LQN models, use the layered network solver
    if input_format in ("lqnx", "xml"):
        return "ln"
    # A place/transition net needs a solver whose feature set declares
    # Transition; MVA cannot answer for one, so `auto` must not resolve there.
    if input_format == "pnml":
        return "ctmc"
    # A portable JSON model carries a Network, a LayeredNetwork or an
    # Environment; peek at the declared model type so each is routed to the
    # solver that reads it. An Environment accepts `env` and nothing else.
    if input_format == "json" and model_file is not None:
        try:
            with open(model_file, "r") as fh:
                model_type = (json.load(fh).get("model") or {}).get("type")
            if model_type == "LayeredNetwork":
                return "ln"
            if model_type == "Environment":
                return "env"
        except (OSError, ValueError):
            pass
    # For JMT formats and Network JSON, use MVA as the default analytical solver
    return "mva"


INPUT_FORMATS: Dict[str, Dict[str, Any]] = {
    "jsim": {
        "name": "JSIM",
        "description": "JMT simulation model format",
        "extensions": [".jsim"],
    },
    "jsimg": {
        "name": "JSIMG",
        "description": "JMT simulation model with graphics",
        "extensions": [".jsimg"],
    },
    "jsimw": {
        "name": "JSIMW",
        "description": "JMT simulation workspace format",
        "extensions": [".jsimw"],
    },
    "lqnx": {
        "name": "LQNX",
        "description": "Layered Queueing Network XML format",
        "extensions": [".lqnx"],
    },
    "xml": {
        "name": "XML",
        "description": "Generic XML model format",
        "extensions": [".xml"],
    },
    "json": {
        "name": "JSON",
        "description": "LINE portable model format (line-model.schema.json)",
        "extensions": [".json"],
    },
    # The JAR has read PNML since the place/transition import landed; this
    # wrapper simply never listed it, so a .pnml path could not be auto-detected
    # and `-i pnml` was not offered.
    "pnml": {
        "name": "PNML",
        "description": "Place/transition net (ISO/IEC 15909-2)",
        "extensions": [".pnml"],
    },
}

OUTPUT_FORMATS: Dict[str, Dict[str, str]] = {
    "table": {
        "name": "Table",
        "description": "Human-readable table format",
    },
    "json": {
        "name": "JSON",
        "description": "JSON format for programmatic use",
    },
    "csv": {
        "name": "CSV",
        "description": "Comma-separated values format",
    },
    "raw": {
        "name": "Raw",
        "description": "Raw output from JAR (no formatting)",
    },
}

ANALYSIS_TYPES: Dict[str, str] = {
    # Basic
    "all": "Both average and system metrics",
    "avg": "Average performance metrics only",
    "sys": "System-level metrics only",
    "stage": "Stage-based metrics (multi-stage service)",
    "chain": "Chain-level averages",
    "node": "Node-level averages",
    "nodechain": "Node-chain level averages",
    # Cache and the other station-class tables the JAR publishes beside the
    # AvgTable. All were reachable from the JAR and absent from this table, so
    # the wrapper refused analyses its own back end serves.
    "cache": "Cache hit/miss metrics",
    "item": "Per-item cache metrics",
    "orbit": "Retrial orbit metrics",
    "loss": "Class-level loss metrics",
    "region-loss": "Finite-capacity region loss metrics",
    "deadline": "Deadline-miss metrics (EDD/EDF)",
    "normconst": "Log normalizing constant (nc, mva)",
    "busyperiod": "Subnetwork busy period (nc, ldes; see --busyperiod*)",
    "sens": "Sensitivity of the means to the service demands",
    # Distribution
    "cdf-respt": "Response time CDF",
    "cdf-passt": "Passage time CDF",
    "perct-respt": "Response time percentiles (MAM solver)",
    # Transient
    "tran-avg": "Transient average metrics",
    "tran-cdf-respt": "Transient response time CDF",
    "tran-cdf-passt": "Transient passage time CDF",
    # Probability
    "prob": "State probability at node (requires --node)",
    "prob-aggr": "Aggregated state probability (requires --node)",
    "prob-marg": "Marginal state probability (requires --node, --class-idx)",
    "prob-sys": "System state probability",
    "prob-sys-aggr": "Aggregated system state probability",
    "prob-sys-marg": "System marginal probability (requires --state)",
    # Sampling (SSA only)
    "sample": "Sample node state trajectory (requires --node)",
    "sample-aggr": "Sample aggregated node state (requires --node)",
    "sample-sys": "Sample system state trajectory",
    "sample-sys-aggr": "Sample aggregated system state",
    # Reward (CTMC only)
    "reward": "Compute reward metrics",
    "reward-steady": "Steady-state reward",
    "reward-value": "Reward value function (requires --reward-name)",
    # Solver-internal structures
    "generator": "CTMC infinitesimal generator and state space",
    "statevec": "Fluid ODE state vector",
    "moments": "Second-order moment-closure report (fld)",
    "interval": "Design-point envelope (uq)",
}

# Analysis types that require specific solvers
# KEPT IN STEP WITH `LineCLI.ANALYSIS_SOLVER_COMPAT`. A narrower gate here is
# not a conservative one -- it refuses, before the JAR is even started, a solve
# the JAR performs.
ANALYSIS_SOLVER_COMPAT: Dict[str, List[str]] = {
    "sample": ["ssa"],
    "sample-aggr": ["ssa"],
    "sample-sys": ["ssa"],
    "sample-sys-aggr": ["ssa"],
    "reward": ["ctmc"],
    "reward-steady": ["ctmc"],
    "reward-value": ["ctmc"],
    "perct-respt": ["mam"],
    "prob": ["ctmc", "ssa"],
    "prob-aggr": ["ctmc", "ssa", "fld"],
    "prob-marg": ["ctmc", "ssa", "mva", "nc", "mam"],
    "prob-sys": ["ctmc", "ssa"],
    "prob-sys-aggr": ["ctmc", "ssa"],
    "prob-sys-marg": ["ctmc", "ssa", "mva", "nc", "mam"],
    "normconst": ["nc", "mva"],
    "busyperiod": ["nc", "ldes"],
    "generator": ["ctmc"],
    "statevec": ["fld"],
    "moments": ["fld"],
    "interval": ["uq"],
}

# Analysis types that require node index
ANALYSIS_REQUIRES_NODE: List[str] = ["prob", "prob-aggr", "prob-marg", "sample", "sample-aggr"]

# Analysis types that require class index
ANALYSIS_REQUIRES_CLASS: List[str] = ["prob-marg"]

# Valid built-in reward names
VALID_REWARD_NAMES: List[str] = ["QLen", "Tput", "Util", "RespT", "WaitT", "ArvR", "ResidT"]

FORMAT_EXTENSIONS: Dict[str, str] = {
    ".jsimg": "jsimg",
    ".jsimw": "jsimw",
    ".jsim": "jsim",
    ".lqnx": "lqnx",
    ".xml": "xml",
    ".json": "json",
    ".pnml": "pnml",
}

# =============================================================================
# Configuration
# =============================================================================

CONFIG_DIR = Path("~/.config/line-cli").expanduser()
CONFIG_FILE = CONFIG_DIR / "config.yaml"


@dataclass
class ServerConfig:
    """Server configuration settings."""
    host: str = "localhost"
    port: int = 5863


@dataclass
class Config:
    """LINE CLI configuration."""
    jar_path: Optional[Path] = None
    java_path: str = "java"
    default_solver: str = "auto"
    default_output_format: str = "table"
    server: ServerConfig = field(default_factory=ServerConfig)

    def __post_init__(self):
        if self.jar_path is None:
            self.jar_path = self._find_jar()

    def _find_jar(self) -> Optional[Path]:
        """Find jline.jar in common locations."""
        env_jar = os.environ.get("LINE_JAR_PATH")
        if env_jar:
            jar_path = Path(env_jar)
            if jar_path.exists():
                return jar_path

        script_dir = Path(__file__).parent
        possible_paths = [
            script_dir / "jline.jar",  # JAR in same directory as script (dist branch)
            script_dir / "common" / "jline.jar",
            script_dir.parent / "common" / "jline.jar",
            Path.cwd() / "common" / "jline.jar",
            Path.cwd().parent / "common" / "jline.jar",
        ]

        for path in possible_paths:
            if path.exists():
                return path.resolve()

        return None


def load_config() -> Config:
    """Load configuration from file or use defaults."""
    config = Config()

    if CONFIG_FILE.exists():
        try:
            import yaml
            with open(CONFIG_FILE) as f:
                data = yaml.safe_load(f) or {}

            if "jar_path" in data and data["jar_path"]:
                config.jar_path = Path(data["jar_path"])
            if "java_path" in data:
                config.java_path = data["java_path"]
            if "default_solver" in data:
                config.default_solver = data["default_solver"]
            if "default_output_format" in data:
                config.default_output_format = data["default_output_format"]
            if "server" in data:
                server_data = data["server"]
                if "host" in server_data:
                    config.server.host = server_data["host"]
                if "port" in server_data:
                    config.server.port = server_data["port"]
        except ImportError:
            pass  # yaml not installed, use defaults
        except Exception:
            pass  # Use defaults on error

    return config


# =============================================================================
# Output Parsing
# =============================================================================

def parse_readable_output(output: str) -> SolveResult:
    """Parse human-readable table output from JAR."""
    if not output.strip():
        return SolveResult(
            success=False,
            error_message="Empty output from solver",
        )

    error_patterns = [
        r"^Error:",
        r"Exception:",
        r"^Invalid",
        r"Unknown solver",
        r"Unsupported",
    ]
    for pattern in error_patterns:
        if re.search(pattern, output, re.MULTILINE):
            return SolveResult(
                success=False,
                raw_output=output,
                error_message=output.strip(),
            )

    avg_table: List[MetricRow] = []
    sys_table: List[MetricRow] = []

    lines = output.strip().split("\n")

    avg_lines = []
    sys_lines = []
    in_sys_section = False

    for line in lines:
        if line.strip().startswith("Chain") and "JobClasses" in line:
            in_sys_section = True

        if in_sys_section:
            sys_lines.append(line)
        else:
            avg_lines.append(line)

    avg_table = _parse_avg_table(avg_lines)
    sys_table = _parse_sys_table(sys_lines)

    return SolveResult(
        success=True,
        avg_table=avg_table if avg_table else None,
        sys_table=sys_table if sys_table else None,
        raw_output=output,
    )


def _parse_avg_table(lines: List[str]) -> List[MetricRow]:
    """Parse the average metrics table section."""
    result: List[MetricRow] = []

    header_line = None
    header_idx = -1
    for i, line in enumerate(lines):
        if "Station" in line and "JobClass" in line:
            header_line = line
            header_idx = i
            break

    if header_line is None:
        return result

    body_lines = []
    for line in lines[header_idx + 1:]:
        if not line.strip() or line.strip().startswith("==="):
            break
        if re.match(r'^[\-=]+$', line.strip()):
            continue
        body_lines.append(line)

    bounds = _column_bounds([header_line] + body_lines)
    headers = [header_line[a:b].strip() for a, b in bounds]

    if len(headers) < 3:
        return result

    metric_names = headers[2:]

    for line in body_lines:
        parts = [line[a:b].strip() for a, b in bounds]

        if len(parts) >= 3:
            station = parts[0]
            job_class = parts[1]

            if station.lower() == "station" or not station:
                continue

            for i, metric_name in enumerate(metric_names):
                if i + 2 < len(parts):
                    value_str = parts[i + 2]
                    try:
                        value = _parse_float(value_str)
                        result.append(MetricRow(
                            station=station,
                            job_class=job_class,
                            metric=metric_name,
                            value=value,
                            unit="",
                        ))
                    except ValueError:
                        continue

    return result


def _column_bounds(lines: List[str]) -> List[Tuple[int, int]]:
    """Column spans of a whitespace-aligned table.

    The JAR pads every column to a common width and right-aligns the numeric
    ones, header included, so a numeric value can start left of its header.
    Taking the spans from the header text alone therefore cuts values in half;
    the separators are the runs of at least two positions that are blank in
    every line of the table.
    """
    width = max(len(line) for line in lines) if lines else 0
    padded = [line.ljust(width) for line in lines]
    blank = [all(line[i] == ' ' for line in padded) for i in range(width)]

    bounds: List[Tuple[int, int]] = []
    start = 0
    i = 0
    while i < width:
        if not blank[i]:
            i += 1
            continue
        run = i
        while run < width and blank[run]:
            run += 1
        if run - i >= 2 and i > start:
            bounds.append((start, i))
            start = run
        i = run
    if start < width:
        bounds.append((start, width))
    return bounds


def _parse_sys_table(lines: List[str]) -> List[MetricRow]:
    """Parse the system metrics table section."""
    result: List[MetricRow] = []

    if not lines:
        return result

    header_line = None
    header_idx = -1
    for i, line in enumerate(lines):
        if "Chain" in line and ("SysRespT" in line or "SysTput" in line):
            header_line = line
            header_idx = i
            break

    if header_line is None:
        return result

    header_pattern = re.compile(r'(\S+)')
    headers = header_pattern.findall(header_line)

    metric_names = headers[2:] if len(headers) > 2 else []

    for line in lines[header_idx + 1:]:
        if re.match(r'^[\-=]+$', line.strip()):
            continue
        if not line.strip():
            continue

        match = re.match(r'(\S+)\s+\(([^)]+)\)\s+(.*)', line)
        if match:
            chain = match.group(1).strip()
            job_classes = match.group(2).strip()
            values_str = match.group(3).strip()

            values = values_str.split()

            for i, metric_name in enumerate(metric_names):
                if i < len(values):
                    try:
                        value = _parse_float(values[i])
                        result.append(MetricRow(
                            station=chain,
                            job_class=job_classes,
                            metric=metric_name,
                            value=value,
                            unit="",
                        ))
                    except ValueError:
                        continue

    return result


def _parse_fixed_width_line(line: str, positions: List[int], num_cols: int) -> List[str]:
    """Parse a fixed-width line based on column positions."""
    parts = []
    for i in range(len(positions)):
        start = positions[i]
        if i + 1 < len(positions):
            end = positions[i + 1]
        else:
            end = len(line)

        if start < len(line):
            parts.append(line[start:end])
        else:
            parts.append("")

    return parts


def _parse_float(value_str: str) -> float:
    """Parse a float value, handling special cases."""
    value_str = value_str.strip()
    if not value_str:
        raise ValueError("Empty value")

    value_lower = value_str.lower()
    if value_lower == "nan":
        return float("nan")
    elif value_lower == "inf" or value_lower == "infinity":
        return float("inf")
    elif value_lower == "-inf" or value_lower == "-infinity":
        return float("-inf")

    return float(value_str)


def parse_json_output(output: str) -> SolveResult:
    """Parse JSON formatted output from JAR."""
    if not output.strip():
        return SolveResult(
            success=False,
            error_message="Empty output from solver",
        )

    try:
        data = json.loads(output)
    except json.JSONDecodeError:
        return parse_readable_output(output)

    if isinstance(data, list):
        avg_table: List[MetricRow] = []
        sys_table: List[MetricRow] = []

        for item in data:
            if isinstance(item, dict) and "data" in item:
                table_data = item["data"]

                if "Station" in table_data and "JobClass" in table_data:
                    parsed = _parse_avg_table(table_data.split("\n"))
                    avg_table.extend(parsed)
                elif "Chain" in table_data and "JobClasses" in table_data:
                    parsed = _parse_sys_table(table_data.split("\n"))
                    sys_table.extend(parsed)

        return SolveResult(
            success=True,
            avg_table=avg_table if avg_table else None,
            sys_table=sys_table if sys_table else None,
            raw_output=output,
        )

    if isinstance(data, dict):
        avg_table = []
        sys_table = []

        if "avgTable" in data:
            avg_data = data["avgTable"]
            if isinstance(avg_data, list):
                for row in avg_data:
                    if isinstance(row, dict):
                        avg_table.append(MetricRow(
                            station=row.get("station", ""),
                            job_class=row.get("class", row.get("jobClass", "")),
                            metric=row.get("metric", ""),
                            value=float(row.get("value", 0)),
                            unit=row.get("unit", ""),
                        ))

        if "sysTable" in data:
            sys_data = data["sysTable"]
            if isinstance(sys_data, list):
                for row in sys_data:
                    if isinstance(row, dict):
                        sys_table.append(MetricRow(
                            station=row.get("station", ""),
                            job_class=row.get("class", row.get("jobClass", "")),
                            metric=row.get("metric", ""),
                            value=float(row.get("value", 0)),
                            unit=row.get("unit", ""),
                        ))

        return SolveResult(
            success=True,
            avg_table=avg_table if avg_table else None,
            sys_table=sys_table if sys_table else None,
            raw_output=output,
        )

    return SolveResult(
        success=False,
        raw_output=output,
        error_message="Unexpected JSON format",
    )


def extract_error_message(output: str, stderr: str) -> Optional[str]:
    """Extract meaningful error message from JAR output."""
    if stderr.strip():
        match = re.search(r"Exception[^:]*:\s*(.+)", stderr)
        if match:
            return match.group(1).strip()
        return stderr.strip()

    if output.strip():
        match = re.search(r"Error[^:]*:\s*(.+)", output)
        if match:
            return match.group(1).strip()

    return None


# =============================================================================
# Format Utilities
# =============================================================================

def detect_format(file_path: Path) -> Optional[str]:
    """Auto-detect format from file extension."""
    ext = file_path.suffix.lower()
    return FORMAT_EXTENSIONS.get(ext)


def validate_format(format_str: str) -> bool:
    """Check if format string is valid."""
    return format_str.lower() in set(FORMAT_EXTENSIONS.values())


# =============================================================================
# JAR Runner
# =============================================================================

class JarRunnerError(Exception):
    """Exception raised for JAR runner errors."""
    pass


class JarRunner:
    """Runner for LINE JAR execution."""

    def __init__(self, jar_path: Optional[Path] = None, java_path: str = "java"):
        self.jar_path = jar_path
        self.java_path = java_path

    def _validate_jar(self) -> None:
        """Validate that JAR exists and is accessible."""
        if self.jar_path is None:
            raise JarRunnerError(
                "JAR path not configured. Set LINE_JAR_PATH environment variable "
                "or specify jar_path in configuration."
            )
        if not self.jar_path.exists():
            raise JarRunnerError(f"JAR file not found: {self.jar_path}")

    def _validate_java(self) -> None:
        """Validate that Java is available."""
        try:
            result = subprocess.run(
                [self.java_path, "-version"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            if result.returncode != 0:
                raise JarRunnerError(f"Java not available: {result.stderr}")
        except FileNotFoundError:
            raise JarRunnerError(
                f"Java executable not found: {self.java_path}. "
                "Ensure Java is installed and in PATH."
            )
        except subprocess.TimeoutExpired:
            raise JarRunnerError("Java version check timed out.")

    def _build_command(self, model_path: Path, options: SolveOptions) -> List[str]:
        """Build the command line for JAR execution."""
        cmd = [
            self.java_path,
            "-jar",
            str(self.jar_path),
            "-f", str(model_path),
            "-i", options.input_format,
            "-o", options.output_format,
            "-s", options.solver,
            "-a", options.analysis,
        ]

        if options.seed is not None:
            cmd.extend(["-d", str(options.seed)])

        if options.verbose:
            cmd.extend(["-v", "normal"])

        # Extended analysis parameters
        if options.node_index is not None:
            cmd.extend(["-n", str(options.node_index)])

        if options.class_index is not None:
            cmd.extend(["-c", str(options.class_index)])

        if options.state is not None:
            cmd.extend(["--state", options.state])

        if options.num_events != 1000:
            cmd.extend(["--events", str(options.num_events)])

        if options.percentiles != "50,90,95,99":
            cmd.extend(["--percentiles", options.percentiles])

        if options.reward_name is not None:
            cmd.extend(["--reward-name", options.reward_name])

        # Straight pass-through: the flag names are the JAR's own, so this table
        # is the whole of the mapping and a flag added there needs one row here.
        for _flag, _value in (
            ("--samples", options.samples),
            ("--cutoff", options.cutoff),
            ("--timespan", options.timespan),
            ("--timestep", options.timestep),
            ("--method", options.method),
            ("--tol", options.tol),
            ("--iter_tol", options.iter_tol),
            ("--iter_max", options.iter_max),
            ("--multiserver", options.multiserver),
            ("--warmupfrac", options.warmupfrac),
            ("--stage-solver", options.stage_solver),
            ("--uq-solver", options.uq_solver),
            ("--busyperiod", options.busyperiod),
            ("--busyperiod-subnet", options.busyperiod_subnet),
            ("--sens-method", options.sens_method),
            ("--sens-scheme", options.sens_scheme),
            ("--sens-step", options.sens_step),
        ):
            if _value is not None:
                cmd.extend([_flag, str(_value)])

        return cmd

    def solve(self, model_path: Path, options: SolveOptions) -> SolveResult:
        """Execute solver on a model file."""
        self._validate_jar()
        self._validate_java()

        if not model_path.exists():
            return SolveResult(
                success=False,
                error_message=f"Model file not found: {model_path}",
            )

        solver_info = SOLVERS.get(options.solver)
        if solver_info is None:
            return SolveResult(
                success=False,
                error_message=f"Unknown solver: {options.solver}",
            )

        if options.input_format not in solver_info["formats"]:
            return SolveResult(
                success=False,
                error_message=(
                    f"Solver '{options.solver}' does not support format '{options.input_format}'. "
                    f"Supported formats: {', '.join(solver_info['formats'])}"
                ),
            )

        cmd = self._build_command(model_path, options)

        # Warn user if JMT solver may need to download JMT.jar
        if options.solver == "jmt":
            jmt_paths = [
                self.jar_path.parent / "JMT.jar" if self.jar_path else None,
                Path.home() / ".jmt" / "JMT.jar",
                Path("/usr/share/jmt/JMT.jar"),
                Path("/opt/jmt/JMT.jar"),
            ]
            jmt_found = any(p and p.exists() for p in jmt_paths)
            if not jmt_found:
                print("Note: JMT.jar not found. It will be downloaded automatically (~50MB). This may take several minutes, please hold.", file=sys.stderr)

        start_time = time.time()
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
            )
            execution_time = time.time() - start_time

            if result.returncode != 0:
                error_msg = extract_error_message(result.stdout, result.stderr)
                return SolveResult(
                    success=False,
                    raw_output=result.stdout,
                    error_message=error_msg or f"Solver exited with code {result.returncode}",
                    execution_time=execution_time,
                )

            if options.output_format == "json":
                solve_result = parse_json_output(result.stdout)
            else:
                solve_result = parse_readable_output(result.stdout)

            solve_result.execution_time = execution_time
            return solve_result

        except subprocess.TimeoutExpired:
            return SolveResult(
                success=False,
                error_message="Solver execution timed out (5 minutes)",
                execution_time=300.0,
            )
        except Exception as e:
            return SolveResult(
                success=False,
                error_message=f"Execution error: {e}",
            )

    def solve_stdin(self, model_content: str, options: SolveOptions) -> SolveResult:
        """Execute solver with model content from stdin."""
        self._validate_jar()
        self._validate_java()

        cmd = [
            self.java_path,
            "-jar",
            str(self.jar_path),
            "-i", options.input_format,
            "-o", options.output_format,
            "-s", options.solver,
            "-a", options.analysis,
        ]

        if options.seed is not None:
            cmd.extend(["-d", str(options.seed)])

        if options.verbose:
            cmd.extend(["-v", "normal"])

        # Extended analysis parameters
        if options.node_index is not None:
            cmd.extend(["-n", str(options.node_index)])

        if options.class_index is not None:
            cmd.extend(["-c", str(options.class_index)])

        if options.state is not None:
            cmd.extend(["--state", options.state])

        if options.num_events != 1000:
            cmd.extend(["--events", str(options.num_events)])

        if options.percentiles != "50,90,95,99":
            cmd.extend(["--percentiles", options.percentiles])

        if options.reward_name is not None:
            cmd.extend(["--reward-name", options.reward_name])

        # Straight pass-through: the flag names are the JAR's own, so this table
        # is the whole of the mapping and a flag added there needs one row here.
        for _flag, _value in (
            ("--samples", options.samples),
            ("--cutoff", options.cutoff),
            ("--timespan", options.timespan),
            ("--timestep", options.timestep),
            ("--method", options.method),
            ("--tol", options.tol),
            ("--iter_tol", options.iter_tol),
            ("--iter_max", options.iter_max),
            ("--multiserver", options.multiserver),
            ("--warmupfrac", options.warmupfrac),
            ("--stage-solver", options.stage_solver),
            ("--uq-solver", options.uq_solver),
            ("--busyperiod", options.busyperiod),
            ("--busyperiod-subnet", options.busyperiod_subnet),
            ("--sens-method", options.sens_method),
            ("--sens-scheme", options.sens_scheme),
            ("--sens-step", options.sens_step),
        ):
            if _value is not None:
                cmd.extend([_flag, str(_value)])

        # Warn user if JMT solver may need to download JMT.jar
        if options.solver == "jmt":
            jmt_paths = [
                self.jar_path.parent / "JMT.jar" if self.jar_path else None,
                Path.home() / ".jmt" / "JMT.jar",
                Path("/usr/share/jmt/JMT.jar"),
                Path("/opt/jmt/JMT.jar"),
            ]
            jmt_found = any(p and p.exists() for p in jmt_paths)
            if not jmt_found:
                print("Note: JMT.jar not found. It will be downloaded automatically (~50MB). This may take several minutes, please hold.", file=sys.stderr)

        start_time = time.time()
        try:
            result = subprocess.run(
                cmd,
                input=model_content,
                capture_output=True,
                text=True,
                timeout=300,
            )
            execution_time = time.time() - start_time

            if result.returncode != 0:
                error_msg = extract_error_message(result.stdout, result.stderr)
                return SolveResult(
                    success=False,
                    raw_output=result.stdout,
                    error_message=error_msg or f"Solver exited with code {result.returncode}",
                    execution_time=execution_time,
                )

            if options.output_format == "json":
                solve_result = parse_json_output(result.stdout)
            else:
                solve_result = parse_readable_output(result.stdout)

            solve_result.execution_time = execution_time
            return solve_result

        except subprocess.TimeoutExpired:
            return SolveResult(
                success=False,
                error_message="Solver execution timed out (5 minutes)",
                execution_time=300.0,
            )
        except Exception as e:
            return SolveResult(
                success=False,
                error_message=f"Execution error: {e}",
            )

    def start_server(self, host: str = "localhost", port: int = 5863) -> subprocess.Popen:
        """Start LINE in WebSocket server mode."""
        self._validate_jar()
        self._validate_java()

        cmd = [
            self.java_path,
            "-jar",
            str(self.jar_path),
            "-p", str(port),
        ]

        return subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    def start_rest_server(self, port: int = 8080) -> subprocess.Popen:
        """Start the WebSocket server behind the 'rest' subcommand.

        The JAR has no HTTP mode, so this is start_server under another name; the
        HTTP API is the separate rest-api module (jline.rest.LineRestServer).
        """
        self._validate_jar()
        self._validate_java()

        cmd = [
            self.java_path,
            "-jar",
            str(self.jar_path),
            "-p", str(port),
        ]

        return subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    def get_version(self) -> str:
        """Get JAR version."""
        self._validate_jar()
        self._validate_java()

        try:
            result = subprocess.run(
                [self.java_path, "-jar", str(self.jar_path), "-V"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            return result.stdout.strip() or "Unknown"
        except Exception:
            return "Unknown"

    def get_java_version(self) -> str:
        """Get Java version."""
        try:
            result = subprocess.run(
                [self.java_path, "-version"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            output = result.stderr or result.stdout
            if output:
                first_line = output.strip().split("\n")[0]
                return first_line
            return "Unknown"
        except Exception:
            return "Not available"


# =============================================================================
# Display Utilities
# =============================================================================

def print_table(title: str, headers: List[str], rows: List[List[str]], col_widths: Optional[List[int]] = None) -> None:
    """Print a formatted table."""
    if col_widths is None:
        col_widths = [len(h) for h in headers]
        for row in rows:
            for i, cell in enumerate(row):
                if i < len(col_widths):
                    col_widths[i] = max(col_widths[i], len(str(cell)))

    print(f"\n{title}")
    print("=" * len(title))

    header_line = "  ".join(h.ljust(col_widths[i]) for i, h in enumerate(headers))
    print(header_line)
    print("-" * len(header_line))

    for row in rows:
        row_line = "  ".join(str(cell).ljust(col_widths[i]) if i < len(col_widths) else str(cell)
                           for i, cell in enumerate(row))
        print(row_line)


def display_results(result: SolveResult, show_avg: bool = True, show_sys: bool = True) -> None:
    """Display solve results as tables."""
    if not result.success:
        print(f"Error: {result.error_message or 'Unknown error'}", file=sys.stderr)
        return

    printed = False

    if show_avg and result.avg_table:
        headers = ["Station", "Class", "Metric", "Value", "Unit"]
        rows = []
        for row in result.avg_table:
            value_str = f"{row.value:.6g}" if isinstance(row.value, float) else str(row.value)
            rows.append([row.station, row.job_class, row.metric, value_str, row.unit])
        print_table("Average Metrics", headers, rows)
        printed = True

    if show_sys and result.sys_table:
        headers = ["Station", "Class", "Metric", "Value", "Unit"]
        rows = []
        for row in result.sys_table:
            value_str = f"{row.value:.6g}" if isinstance(row.value, float) else str(row.value)
            rows.append([row.station, row.job_class, row.metric, value_str, row.unit])
        print_table("System Metrics", headers, rows)
        printed = True

    # AN ANALYSIS THIS FORMATTER CANNOT RE-PARSE MUST STILL BE SHOWN. The two
    # parsers above know the AvgTable and AvgSysTable layouts and nothing else,
    # so every other analysis -- `cache`, `item`, `orbit`, `normconst`,
    # `busyperiod`, `sens`, the probability and reward families -- rendered as an
    # empty screen while `-o raw` showed the JAR had answered in full. Echoing
    # the solver's own text is the honest fallback: the wrapper does not
    # understand the section, so it does not reformat it.
    #
    # THE CONDITION IS "NOTHING WAS PRINTED", not "the tables are empty". Those
    # are different: `-a sens` produces a table whose header carries Station and
    # JobClass, so `_parse_avg_table` claims it, while `show_avg` is false for
    # any analysis that is not avg/sys/all -- so the rows were parsed, then
    # suppressed, and an emptiness test on the tables would have stayed silent.
    if not printed and result.raw_output:
        text = result.raw_output.strip()
        if text:
            print(text)

    if result.execution_time > 0:
        print(f"\nExecution time: {result.execution_time:.3f}s")


def display_solver_list() -> None:
    """Display available solvers with descriptions."""
    headers = ["Solver", "Name", "Description", "Formats"]
    rows = []
    for solver_id, info in SOLVERS.items():
        rows.append([solver_id, info["name"], info["description"], ", ".join(info["formats"])])
    print_table("Available Solvers", headers, rows)


def display_format_list() -> None:
    """Display supported input and output formats."""
    headers = ["Format", "Name", "Description", "Extensions"]
    rows = []
    for fmt_id, info in INPUT_FORMATS.items():
        rows.append([fmt_id, info["name"], info["description"], ", ".join(info["extensions"])])
    print_table("Input Formats", headers, rows)

    headers = ["Format", "Name", "Description"]
    rows = []
    for fmt_id, info in OUTPUT_FORMATS.items():
        rows.append([fmt_id, info["name"], info["description"]])
    print_table("Output Formats", headers, rows)


def display_analysis_types() -> None:
    """Display analysis types."""
    headers = ["Type", "Description"]
    rows = [[type_id, desc] for type_id, desc in ANALYSIS_TYPES.items()]
    print_table("Analysis Types", headers, rows)


# =============================================================================
# Commands
# =============================================================================

def validate_analysis_types(analysis_str: str) -> List[str]:
    """Validate comma-separated analysis types."""
    types = [t.strip() for t in analysis_str.split(',')]
    for t in types:
        if t not in ANALYSIS_TYPES:
            raise ValueError(f"Invalid analysis type: '{t}'. Valid types: {', '.join(ANALYSIS_TYPES.keys())}")
    return types


def validate_analysis_solver_compat(analysis_types: List[str], solver: str) -> None:
    """Validate that analysis types are compatible with solver."""
    for analysis_type in analysis_types:
        required_solvers = ANALYSIS_SOLVER_COMPAT.get(analysis_type)
        if required_solvers and solver not in required_solvers:
            raise ValueError(
                f"Analysis type '{analysis_type}' requires solver: {' or '.join(required_solvers)}, "
                f"but '{solver}' was specified."
            )


def validate_analysis_params(analysis_types: List[str], node_idx: Optional[int],
                             class_idx: Optional[int], reward_name: Optional[str]) -> None:
    """Validate that required parameters are provided for analysis types."""
    for analysis_type in analysis_types:
        if analysis_type in ANALYSIS_REQUIRES_NODE and node_idx is None:
            raise ValueError(f"Analysis type '{analysis_type}' requires --node parameter.")
        if analysis_type in ANALYSIS_REQUIRES_CLASS and class_idx is None:
            raise ValueError(f"Analysis type '{analysis_type}' requires --class-idx parameter.")
        if analysis_type == 'reward-value' and not reward_name:
            raise ValueError("Analysis type 'reward-value' requires --reward-name parameter.")


def cmd_solve(args: argparse.Namespace) -> int:
    """Execute the solve command."""
    config = load_config()

    # Resolve solver alias
    solver = resolve_solver(args.solver)

    if solver not in SOLVERS:
        print(f"Error: Unknown solver '{args.solver}'. Available: {', '.join(SOLVERS.keys())}", file=sys.stderr)
        return 1

    # Determine input format early for auto solver selection
    input_format = args.input_format
    if args.model_file is not None and input_format is None:
        input_format = detect_format(Path(args.model_file))

    # Handle 'auto' solver - select concrete solver based on input format
    if solver == "auto":
        if input_format is None:
            # Default to mva for unknown formats
            solver = "mva"
        else:
            solver = auto_select_solver(input_format, args.model_file)

    # Validate analysis types
    try:
        analysis_types = validate_analysis_types(args.analysis)
        validate_analysis_solver_compat(analysis_types, solver)
        validate_analysis_params(analysis_types, args.node, args.class_idx, args.reward_name)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # Handle stdin input
    if args.model_file is None:
        if sys.stdin.isatty():
            print("Error: No model file specified and no stdin input", file=sys.stderr)
            print("Provide a model file or pipe input: cat model.jsimg | python line.py solve -i jsimg", file=sys.stderr)
            return 1

        if args.input_format is None:
            print("Error: Input format required when reading from stdin (-i/--input-format)", file=sys.stderr)
            return 1

        model_content = sys.stdin.read()
        if not model_content.strip():
            print("Error: Empty input from stdin", file=sys.stderr)
            return 1

        runner = JarRunner(jar_path=config.jar_path, java_path=config.java_path)
        options = SolveOptions(
            solver=solver,
            input_format=args.input_format,
            output_format="json" if args.output_format == "json" else "readable",
            analysis=args.analysis,
            seed=args.seed,
            verbose=args.verbose,
            node_index=args.node,
            class_index=args.class_idx,
            state=args.state,
            num_events=args.events,
            percentiles=args.percentiles,
            reward_name=args.reward_name,
            samples=args.samples,
            cutoff=args.cutoff,
            timespan=args.timespan,
            timestep=args.timestep,
            method=args.method,
            tol=args.tol,
            iter_tol=args.iter_tol,
            iter_max=args.iter_max,
            multiserver=args.multiserver,
            warmupfrac=args.warmupfrac,
            stage_solver=args.stage_solver,
            uq_solver=args.uq_solver,
            busyperiod=args.busyperiod,
            busyperiod_subnet=args.busyperiod_subnet,
            sens_method=args.sens_method,
            sens_scheme=args.sens_scheme,
            sens_step=args.sens_step,
        )

        if not args.quiet:
            print("Solving model...", file=sys.stderr)

        try:
            result = runner.solve_stdin(model_content, options)
        except JarRunnerError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1

    else:
        model_path = Path(args.model_file)
        if not model_path.exists():
            print(f"Error: Model file not found: {model_path}", file=sys.stderr)
            return 1

        input_format = args.input_format
        if input_format is None:
            input_format = detect_format(model_path)
            if input_format is None:
                print(f"Error: Could not detect format from extension: {model_path.suffix}", file=sys.stderr)
                print("Use -i/--input-format to specify the format", file=sys.stderr)
                return 1
        elif not validate_format(input_format):
            print(f"Error: Invalid input format: {input_format}", file=sys.stderr)
            return 1

        runner = JarRunner(jar_path=config.jar_path, java_path=config.java_path)
        options = SolveOptions(
            solver=solver,
            input_format=input_format,
            output_format="json" if args.output_format == "json" else "readable",
            analysis=args.analysis,
            seed=args.seed,
            verbose=args.verbose,
            node_index=args.node,
            class_index=args.class_idx,
            state=args.state,
            num_events=args.events,
            percentiles=args.percentiles,
            reward_name=args.reward_name,
            samples=args.samples,
            cutoff=args.cutoff,
            timespan=args.timespan,
            timestep=args.timestep,
            method=args.method,
            tol=args.tol,
            iter_tol=args.iter_tol,
            iter_max=args.iter_max,
            multiserver=args.multiserver,
            warmupfrac=args.warmupfrac,
            stage_solver=args.stage_solver,
            uq_solver=args.uq_solver,
            busyperiod=args.busyperiod,
            busyperiod_subnet=args.busyperiod_subnet,
            sens_method=args.sens_method,
            sens_scheme=args.sens_scheme,
            sens_step=args.sens_step,
        )

        if not args.quiet:
            print(f"Solving {model_path.name}...", file=sys.stderr)

        try:
            result = runner.solve(model_path, options)
        except JarRunnerError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1

    if not result.success:
        print(f"Error: {result.error_message or 'Solver failed'}", file=sys.stderr)
        return 1

    show_avg = args.analysis in ("all", "avg")
    show_sys = args.analysis in ("all", "sys")

    output_content = ""

    if args.output_format == "json":
        output_data = {
            "success": result.success,
            "execution_time": result.execution_time,
        }
        if result.avg_table:
            output_data["avgTable"] = [
                {
                    "station": r.station,
                    "class": r.job_class,
                    "metric": r.metric,
                    "value": r.value,
                    "unit": r.unit,
                }
                for r in result.avg_table
            ]
        if result.sys_table:
            output_data["sysTable"] = [
                {
                    "station": r.station,
                    "class": r.job_class,
                    "metric": r.metric,
                    "value": r.value,
                    "unit": r.unit,
                }
                for r in result.sys_table
            ]
        # An analysis with no avg/sys table of its own -- `cache`, `sens`,
        # `busyperiod`, the probability and reward families -- would otherwise
        # be reported as a bare {"success": true} envelope with the answer
        # discarded. The JAR was already asked for -o json, so its own document
        # IS the result and rides under `solver`.
        if "avgTable" not in output_data and "sysTable" not in output_data:
            raw = (result.raw_output or "").strip()
            if raw:
                try:
                    output_data["solver"] = json.loads(raw[raw.index("{"):])
                except (ValueError, json.JSONDecodeError):
                    output_data["solverOutput"] = raw
        output_content = json.dumps(output_data, indent=2)

    elif args.output_format == "csv":
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(["Table", "Station", "Class", "Metric", "Value", "Unit"])

        if show_avg and result.avg_table:
            for r in result.avg_table:
                writer.writerow(["avg", r.station, r.job_class, r.metric, r.value, r.unit])

        if show_sys and result.sys_table:
            for r in result.sys_table:
                writer.writerow(["sys", r.station, r.job_class, r.metric, r.value, r.unit])

        output_content = output.getvalue()

    elif args.output_format == "raw":
        output_content = result.raw_output

    if args.output_file:
        Path(args.output_file).write_text(output_content if output_content else result.raw_output)
        if not args.quiet:
            print(f"Results written to {args.output_file}", file=sys.stderr)
    elif output_content:
        print(output_content)
    else:
        display_results(result, show_avg=show_avg, show_sys=show_sys)

    return 0


def _check_external_tool(name: str) -> bool:
    """Check if an external tool is available in PATH."""
    return shutil.which(name) is not None


def _check_jmt_available(config: Config) -> bool:
    """Check if JMT.jar is available."""
    # Check common JMT locations (same order as JarRunner.solve)
    jmt_paths = []

    # First check in the same directory as jline.jar (common/)
    if config.jar_path:
        jmt_paths.append(config.jar_path.parent / "JMT.jar")

    # Then check other standard locations
    jmt_paths.extend([
        Path.home() / ".jmt" / "JMT.jar",
        Path("/usr/share/jmt/JMT.jar"),
        Path("/opt/jmt/JMT.jar"),
    ])

    env_jmt = os.environ.get("JMT_PATH")
    if env_jmt:
        jmt_paths.insert(0, Path(env_jmt))

    for path in jmt_paths:
        if path.exists():
            return True
    return False


def _check_lqns_available(config: Config) -> bool:
    """Check if LQNS is available."""
    # First check in the same directory as jline.jar (common/)
    if config.jar_path:
        lqns_path = config.jar_path.parent / "lqns"
        if lqns_path.exists():
            return True

    # Then check if lqns is in PATH
    return _check_external_tool("lqns")


def cmd_info(args: argparse.Namespace) -> int:
    """Execute the info command."""
    config = load_config()

    print(f"\nLINE Solver Information")
    print("=" * 40)
    print(f"CLI Version:     {__version__}")
    print(f"Config File:     {CONFIG_FILE if CONFIG_FILE.exists() else 'Not found (using defaults)'}")

    if config.jar_path and config.jar_path.exists():
        print(f"JAR Path:        {config.jar_path}")
        print(f"JAR Status:      Found")

        try:
            runner = JarRunner(jar_path=config.jar_path, java_path=config.java_path)
            jar_version = runner.get_version()
            print(f"JAR Version:     {jar_version}")
        except Exception as e:
            print(f"JAR Version:     Error: {e}")
    else:
        print(f"JAR Path:        Not found")
        print(f"JAR Status:      Missing")

    try:
        runner = JarRunner(jar_path=config.jar_path, java_path=config.java_path)
        java_version = runner.get_java_version()
        print(f"Java:            {java_version}")
    except Exception:
        print(f"Java:            Not available")

    print(f"Default Solver:  {config.default_solver}")
    print(f"Default Output:  {config.default_output_format}")
    print(f"Server Host:     {config.server.host}")
    print(f"Server Port:     {config.server.port}")

    # External tools availability
    print("\nExternal Tools:")
    print("-" * 40)
    jmt_available = _check_jmt_available(config)
    lqns_available = _check_lqns_available(config)
    qnsolver_available = _check_external_tool("qnsolver")

    print(f"JMT:             {'Available' if jmt_available else 'Not found'}")
    print(f"LQNS:            {'Available' if lqns_available else 'Not found'}")
    print(f"QNSolver:        {'Available' if qnsolver_available else 'Not found'}")

    print("\nSolver Compatibility:")
    print("-" * 40)

    jsim_formats = {"jsim", "jsimg", "jsimw"}
    lqn_formats = {"lqnx", "xml"}

    headers = ["Solver", "JSIM/JSIMG/JSIMW", "LQNX/XML"]
    rows = []
    # Sort solvers alphabetically
    for solver_id in sorted(SOLVERS.keys()):
        info = SOLVERS[solver_id]
        supported = set(info["formats"])
        jsim_ok = "Yes" if supported & jsim_formats else "No"
        lqn_ok = "Yes" if supported & lqn_formats else "No"
        rows.append([solver_id, jsim_ok, lqn_ok])

    col_widths = [10, 18, 10]
    header_line = "  ".join(h.ljust(col_widths[i]) for i, h in enumerate(headers))
    print(header_line)
    for row in rows:
        row_line = "  ".join(str(cell).ljust(col_widths[i]) for i, cell in enumerate(row))
        print(row_line)

    print("\nUsage:")
    print("-" * 40)
    print("Specify solver with -s/--solver option:")
    print("  python line-cli.py solve model.jsimg -s mva")
    print("  python line-cli.py solve model.jsimg -s nc")
    print("  python line-cli.py solve model.jsimg -s fld")
    print("  python line-cli.py solve model.lqnx -s ln")
    print()
    print("List all solvers:  python line-cli.py list solvers")
    print("List all options:  python line-cli.py solve --help")

    return 0


def cmd_list(args: argparse.Namespace) -> int:
    """Execute the list command."""
    resource = args.resource

    if resource is None:
        print("\nAvailable resources to list:")
        print("  solvers   - Available solver algorithms")
        print("  formats   - Supported input/output formats")
        print("  analysis  - Analysis types")
        print("\nUsage: python line.py list <resource>")
        return 0

    resource = resource.lower()

    if resource == "solvers":
        display_solver_list()
    elif resource == "formats":
        display_format_list()
    elif resource == "analysis":
        display_analysis_types()
    else:
        print(f"Error: Unknown resource: {resource}", file=sys.stderr)
        print("Available: solvers, formats, analysis", file=sys.stderr)
        return 1

    return 0


def cmd_server(args: argparse.Namespace) -> int:
    """Execute the server command."""
    config = load_config()

    try:
        runner = JarRunner(jar_path=config.jar_path, java_path=config.java_path)
    except JarRunnerError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    print(f"Starting LINE WebSocket server...")
    print(f"  Host: {args.host}")
    print(f"  Port: {args.port}")
    print()

    try:
        process = runner.start_server(host=args.host, port=args.port)
    except JarRunnerError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    print(f"Server started on ws://{args.host}:{args.port}")
    print("Press Ctrl+C to stop the server")
    print()

    def signal_handler(sig, frame):
        print()
        print("Shutting down server...")
        process.terminate()
        try:
            process.wait(timeout=5)
        except Exception:
            process.kill()
        print("Server stopped")
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    try:
        while True:
            line = process.stdout.readline()
            if line:
                print(line.rstrip())
            elif process.poll() is not None:
                break

        returncode = process.returncode
        if returncode != 0:
            stderr = process.stderr.read() if process.stderr else ""
            print(f"Error: Server exited with code {returncode}", file=sys.stderr)
            if stderr:
                print(stderr, file=sys.stderr)
            return returncode

    except KeyboardInterrupt:
        signal_handler(None, None)

    return 0


def cmd_rest(args: argparse.Namespace) -> int:
    """Start the WebSocket server on port 8080; the HTTP API is the rest-api module."""
    config = load_config()

    try:
        runner = JarRunner(jar_path=config.jar_path, java_path=config.java_path)
    except JarRunnerError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    print(f"Starting LINE WebSocket server...")
    print(f"  Port: {args.port}")
    print()

    try:
        process = runner.start_rest_server(port=args.port)
    except JarRunnerError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    print(f"WebSocket server started on ws://localhost:{args.port}")
    print("Press Ctrl+C to stop the server")
    print()
    print("This subcommand is an alias of 'server'; it does not serve HTTP.")
    print("For the HTTP API, build and run the rest-api module:")
    print("  cd io/rest-api && mvn clean package")
    print("  java -cp target/line-rest.jar:common/jline.jar \\")
    print(f"       jline.rest.LineRestServer --port {args.port}")
    print()

    def signal_handler(sig, frame):
        print()
        print("Shutting down WebSocket server...")
        process.terminate()
        try:
            process.wait(timeout=5)
        except Exception:
            process.kill()
        print("Server stopped")
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    try:
        while True:
            line = process.stdout.readline()
            if line:
                print(line.rstrip())
            elif process.poll() is not None:
                break

        returncode = process.returncode
        if returncode != 0:
            stderr = process.stderr.read() if process.stderr else ""
            print(f"Error: WebSocket server exited with code {returncode}", file=sys.stderr)
            if stderr:
                print(stderr, file=sys.stderr)
            return returncode

    except KeyboardInterrupt:
        signal_handler(None, None)

    return 0


# =============================================================================
# Main Entry Point
# =============================================================================

def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        prog="line",
        description="LINE - Queueing network solver command-line interface",
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"LINE CLI version {__version__}",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Solve command
    solve_parser = subparsers.add_parser("solve", help="Solve a queueing network model")
    solve_parser.add_argument(
        "model_file",
        nargs="?",
        help="Path to model file (reads from stdin if not provided)",
    )
    solve_parser.add_argument(
        "-s", "--solver",
        default="auto",
        help=f"Solver algorithm to use (default: auto). Available: {', '.join(sorted(SOLVERS.keys()))}. "
             f"Aliases: {', '.join(f'{k}={v}' for k, v in SOLVER_ALIASES.items())}",
    )
    solve_parser.add_argument(
        "-i", "--input-format",
        help="Input format (auto-detected from extension if not specified)",
    )
    solve_parser.add_argument(
        "-o", "--output-format",
        default="table",
        choices=["table", "json", "csv", "raw"],
        help="Output format (default: table)",
    )
    solve_parser.add_argument(
        "-a", "--analysis",
        default="all",
        help="Analysis type(s), comma-separated for multiple (default: all). "
             f"Available: {', '.join(ANALYSIS_TYPES.keys())}",
    )
    solve_parser.add_argument(
        "-d", "--seed",
        type=int,
        help="Random seed for stochastic solvers",
    )
    solve_parser.add_argument(
        "-n", "--node",
        type=int,
        help="Node index for prob/sample analysis (0-based)",
    )
    solve_parser.add_argument(
        "-c", "--class-idx",
        type=int,
        help="Job class index for prob-marg analysis (0-based)",
    )
    solve_parser.add_argument(
        "--state",
        type=str,
        help="State vector for prob analysis (comma-separated integers)",
    )
    solve_parser.add_argument(
        "--events",
        type=int,
        default=1000,
        help="Number of events for sample analysis (default: 1000)",
    )
    solve_parser.add_argument(
        "--percentiles",
        type=str,
        default="50,90,95,99",
        help="Percentile values for perct-respt (comma-separated, default: 50,90,95,99)",
    )
    solve_parser.add_argument(
        "--reward-name",
        type=str,
        choices=VALID_REWARD_NAMES,
        help=f"Built-in reward name for reward-value analysis: {', '.join(VALID_REWARD_NAMES)}",
    )
    solve_parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output",
    )
    # ---- Numeric controls, forwarded verbatim to the JAR -----------------
    solve_parser.add_argument("--samples", type=int,
        help="Simulation samples / Monte Carlo draws (ssa, ldes, jmt, uq)")
    solve_parser.add_argument("--cutoff", type=float,
        help="State-space cutoff per open class (ctmc, ssa)")
    solve_parser.add_argument("--timespan", "--tspan", dest="timespan", type=str,
        help="Time span of the transient analyses, e.g. --timespan 0,100")
    solve_parser.add_argument("--timestep", type=float,
        help="Fixed transient output step (default: the adaptive ODE grid)")
    solve_parser.add_argument("--method", type=str,
        help="Algorithm within the chosen solver")
    solve_parser.add_argument("--tol", type=float, help="General solver tolerance")
    solve_parser.add_argument("--iter_tol", type=float,
        help="Iteration convergence tolerance")
    solve_parser.add_argument("--iter_max", type=int, help="Maximum iterations")
    solve_parser.add_argument("--multiserver", type=str,
        help="AMVA multiserver rule (seidmann, softmin, ...)")
    solve_parser.add_argument("--warmupfrac", type=float,
        help="Leading fraction of a simulated path discarded before the means")
    solve_parser.add_argument("--stage-solver", type=str,
        help="Solver run at each stage of an Environment ('fluid' or 'ctmc')")
    solve_parser.add_argument("--uq-solver", type=str,
        help="Engine SolverUQ runs at each design point; required by -s uq")
    solve_parser.add_argument("--busyperiod", type=str, metavar="N[,N...]",
        help="Orders of -a busyperiod (default: 1)")
    solve_parser.add_argument("--busyperiod-subnet", type=str, metavar="I[,I...]",
        help="0-based stations forming the -a busyperiod subnetwork (required)")
    solve_parser.add_argument("--sens-method", type=str,
        help="Differentiation of -a sens")
    solve_parser.add_argument("--sens-scheme", type=str,
        help="Finite-difference scheme of -a sens")
    solve_parser.add_argument("--sens-step", type=float,
        help="Finite-difference step of -a sens")
    solve_parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress non-essential output",
    )
    solve_parser.add_argument(
        "-O", "--output-file",
        help="Write results to file instead of stdout",
    )

    # Info command
    subparsers.add_parser("info", help="Display system information and configuration")

    # List command
    list_parser = subparsers.add_parser("list", help="List available resources")
    list_parser.add_argument(
        "resource",
        nargs="?",
        help="Resource to list: solvers, formats, analysis",
    )

    # Server command (WebSocket)
    server_parser = subparsers.add_parser("server", help="Start LINE in WebSocket server mode")
    server_parser.add_argument(
        "-p", "--port",
        type=int,
        default=5863,
        help="WebSocket port to listen on (default: 5863)",
    )
    server_parser.add_argument(
        "-H", "--host",
        default="localhost",
        help="Host address to bind (default: localhost)",
    )

    # REST command (alias of server; the HTTP API lives in io/rest-api/)
    rest_parser = subparsers.add_parser(
        "rest",
        help="Start the WebSocket server on port 8080 (alias of server)",
    )
    rest_parser.add_argument(
        "-p", "--port",
        type=int,
        default=8080,
        help="WebSocket port to listen on (default: 8080)",
    )

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "solve":
        return cmd_solve(args)
    elif args.command == "info":
        return cmd_info(args)
    elif args.command == "list":
        return cmd_list(args)
    elif args.command == "server":
        return cmd_server(args)
    elif args.command == "rest":
        return cmd_rest(args)

    return 0


if __name__ == "__main__":
    sys.exit(main())
