# LINE CLI - Queueing Network Solver

**Version: 3.0.2**

[![Download Sources](https://a.fsdn.com/con/app/sf-download-button)](https://sourceforge.net/projects/line-solver/files/latest/download)

LINE is an open source package for Java, Matlab and Python to analyze queueing models via analytical methods and simulation. The tool features algorithms for the solution of open queueing systems (e.g., M/M/1, M/M/k, M/G/1, ...), open and closed queueing networks, and layered queueing networks.

LINE is developed by the QORE lab at Imperial College London and distributed under the BSD-3 license.

> **Note**: This is a lightweight CLI-only distribution. For full source code with MATLAB, Python, and Java/Kotlin implementations, click the **Download Now** button above.

## Quick Start

```bash
# Solve a queueing model
python line.py solve examples/jsimg/closed_1class_2stat.jsimg

# Solve with specific solver and JSON output
python line.py solve model.jsimg -s fluid -o json

# List available solvers
python line.py list solvers

# Show system info
python line.py info

# Start WebSocket server mode
python line.py server -p 8080
```

## Requirements

- **Python 3.8+** (standard library only, no pip packages needed)
- **Java Runtime Environment (JRE) 8+**

## Installation

1. Download or clone this distribution:
   ```bash
   git clone --single-branch https://github.com/imperial-qore/line-solver.git line-cli
   cd line-cli
   ```

2. Verify setup:
   ```bash
   python line.py info
   ```

Alternatively, set the `LINE_JAR_PATH` environment variable to point to `jline.jar` if placed elsewhere.

## Supported Input Formats

| Format | Extensions | Description |
|--------|------------|-------------|
| JSIMG  | `.jsimg`, `.jsim`, `.jsimw` | JMT queueing network models |
| LQNX   | `.lqnx`, `.xml` | Layered Queueing Network models |

## Available Solvers

| Solver | Description |
|--------|-------------|
| `mva`  | Mean Value Analysis (analytical, fast) |
| `fluid` | Fluid/ODE approximation |
| `ctmc` | Continuous-Time Markov Chain (exact) |
| `jmt`  | Discrete event simulation via JMT |
| `nc`   | Normalizing constant method |
| `ssa`  | Stochastic simulation algorithm |
| `ln`   | Layered network solver |
| `lqns` | External LQN solver |
| `mam`  | Matrix-analytic methods |

## Output Formats

- `table` (default) - Human-readable formatted tables
- `json` - JSON format for programmatic use
- `csv` - Comma-separated values
- `raw` - Raw JAR output

## Examples

The `examples/` directory contains sample models:

- `examples/jsimg/closed_1class_2stat.jsimg` - Closed queueing network
- `examples/jsimg/open_1class_1stat_mg1fcfs.jsimg` - M/G/1 FCFS queue
- `examples/jsimg/open_1class_1stat_mm1k.jsimg` - M/M/1/K finite capacity queue
- `examples/lqn/simple-forwarding.lqnx` - Layered queueing network

## Documentation

- `doc/line-cli-quickstart.md` - CLI quick start guide
- `doc/LINE-java.pdf` - Full Java/Kotlin API documentation
- `doc/LINE-python.pdf` - Python wrapper documentation
- `doc/LINE-cheatsheet.pdf` - Quick reference card

## Full Distribution

For the complete LINE solver with MATLAB, Python, and Java/Kotlin implementations:

- **SourceForge**: https://sourceforge.net/projects/line-solver/
- **PyPI**: `pip install line-solver`
- **Documentation**: https://line-solver.sourceforge.net/

## License

BSD-3-Clause License. See [LICENSE](LICENSE) file.

