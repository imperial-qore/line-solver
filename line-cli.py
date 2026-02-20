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
    python line-cli.py rest -p 8080       # Start REST API server
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
from typing import Dict, List, Optional, Any

__version__ = "3.0.3"

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


# =============================================================================
# Solver and Format Definitions
# =============================================================================

SOLVERS: Dict[str, Dict[str, Any]] = {
    "auto": {
        "name": "Automatic Solver Selection",
        "description": "Automatically selects the best solver for the model",
        "formats": ["jsim", "jsimg", "jsimw", "lqnx", "xml"],
    },
    "ctmc": {
        "name": "Continuous-Time Markov Chain",
        "description": "Exact analysis using CTMC state space exploration",
        "formats": ["jsim", "jsimg", "jsimw"],
    },
    "des": {
        "name": "Discrete Event Simulation",
        "description": "Discrete event simulation using SSJ library",
        "formats": ["jsim", "jsimg", "jsimw"],
    },
    "fld": {
        "name": "Fluid/Mean-Field ODE",
        "description": "Approximate analysis using fluid/ODE model",
        "formats": ["jsim", "jsimg", "jsimw"],
    },
    "jmt": {
        "name": "Java Modelling Tools",
        "description": "Discrete event simulation using JMT",
        "formats": ["jsim", "jsimg", "jsimw"],
    },
    "ln": {
        "name": "Layered Network",
        "description": "Solver for layered queueing networks",
        "formats": ["lqnx", "xml"],
    },
    "lqns": {
        "name": "LQN Solver",
        "description": "External LQNS solver integration",
        "formats": ["lqnx", "xml"],
    },
    "mam": {
        "name": "Matrix Analytic Methods",
        "description": "Analysis using matrix analytic methods (supports Fork-Join percentiles)",
        "formats": ["jsim", "jsimg", "jsimw"],
    },
    "mva": {
        "name": "Mean Value Analysis",
        "description": "Analytical solver using Mean Value Analysis algorithm",
        "formats": ["jsim", "jsimg", "jsimw"],
    },
    "nc": {
        "name": "Normalizing Constant",
        "description": "Exact analysis using normalizing constant computation",
        "formats": ["jsim", "jsimg", "jsimw"],
    },
    "qns": {
        "name": "QNS",
        "description": "External QNSolver integration",
        "formats": ["jsim", "jsimg", "jsimw"],
    },
    "ssa": {
        "name": "Stochastic Simulation Algorithm",
        "description": "Stochastic simulation of the model",
        "formats": ["jsim", "jsimg", "jsimw"],
    },
}

# Solver aliases (map alias -> canonical name)
SOLVER_ALIASES: Dict[str, str] = {
    "fluid": "fld",
    "qnsolver": "qns",
}


def resolve_solver(solver: str) -> str:
    """Resolve solver alias to canonical name."""
    return SOLVER_ALIASES.get(solver.lower(), solver.lower())


def auto_select_solver(input_format: str) -> str:
    """Select an appropriate solver based on input format when 'auto' is specified.

    Returns a concrete solver name since the JAR doesn't support 'auto'.
    """
    # For LQN models, use the layered network solver
    if input_format in ("lqnx", "xml"):
        return "ln"
    # For JMT formats, use MVA as the default analytical solver
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
    # Sampling (SSA only)
    "sample": "Sample node state trajectory (requires --node)",
    "sample-aggr": "Sample aggregated node state (requires --node)",
    "sample-sys": "Sample system state trajectory",
    "sample-sys-aggr": "Sample aggregated system state",
    # Reward (CTMC only)
    "reward": "Compute reward metrics",
    "reward-steady": "Steady-state reward",
    "reward-value": "Reward value function (requires --reward-name)",
}

# Analysis types that require specific solvers
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
    "prob-aggr": ["ctmc", "ssa"],
    "prob-marg": ["ctmc", "ssa"],
    "prob-sys": ["ctmc", "ssa"],
    "prob-sys-aggr": ["ctmc", "ssa"],
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

    header_pattern = re.compile(r'(\S+)')
    headers = header_pattern.findall(header_line)

    if len(headers) < 3:
        return result

    metric_names = headers[2:]

    col_positions = []
    pos = 0
    for header in headers:
        idx = header_line.find(header, pos)
        col_positions.append(idx)
        pos = idx + len(header)

    for line in lines[header_idx + 1:]:
        if re.match(r'^[\-=]+$', line.strip()):
            continue
        if not line.strip():
            continue

        parts = _parse_fixed_width_line(line, col_positions, len(headers))

        if len(parts) >= 3:
            station = parts[0].strip()
            job_class = parts[1].strip()

            if station.lower() == "station" or not station:
                continue

            for i, metric_name in enumerate(metric_names):
                if i + 2 < len(parts):
                    value_str = parts[i + 2].strip()
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
        """Start LINE in REST API server mode.

        Note: The JAR uses WebSocket server mode (-p) for both WebSocket and REST.
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

    if show_avg and result.avg_table:
        headers = ["Station", "Class", "Metric", "Value", "Unit"]
        rows = []
        for row in result.avg_table:
            value_str = f"{row.value:.6g}" if isinstance(row.value, float) else str(row.value)
            rows.append([row.station, row.job_class, row.metric, value_str, row.unit])
        print_table("Average Metrics", headers, rows)

    if show_sys and result.sys_table:
        headers = ["Station", "Class", "Metric", "Value", "Unit"]
        rows = []
        for row in result.sys_table:
            value_str = f"{row.value:.6g}" if isinstance(row.value, float) else str(row.value)
            rows.append([row.station, row.job_class, row.metric, value_str, row.unit])
        print_table("System Metrics", headers, rows)

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
            solver = auto_select_solver(input_format)

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
    """Execute the rest command to start REST API server."""
    config = load_config()

    try:
        runner = JarRunner(jar_path=config.jar_path, java_path=config.java_path)
    except JarRunnerError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    print(f"Starting LINE REST API server...")
    print(f"  Port: {args.port}")
    print()

    try:
        process = runner.start_rest_server(port=args.port)
    except JarRunnerError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    print(f"REST API server started on http://localhost:{args.port}")
    print("Press Ctrl+C to stop the server")
    print()
    print("Available endpoints:")
    print(f"  POST http://localhost:{args.port}/solve  - Solve a model")
    print(f"  GET  http://localhost:{args.port}/health - Health check")
    print()

    def signal_handler(sig, frame):
        print()
        print("Shutting down REST API server...")
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
            print(f"Error: REST API server exited with code {returncode}", file=sys.stderr)
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

    # REST API command
    rest_parser = subparsers.add_parser("rest", help="Start LINE REST API server")
    rest_parser.add_argument(
        "-p", "--port",
        type=int,
        default=8080,
        help="REST API port to listen on (default: 8080)",
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
