# LINE CLI Quick Start

Analyze queueing models with `line.py`. **Requirements:** Java 8+, Python 3.x

## Installation

```bash
# Download and extract
wget https://sf.net/projects/line-solver/files/latest/download -O line-solver.zip
unzip line-solver.zip && cd line-solver

# Verify installation
python line.py info
```

If JAR not found: `export LINE_JAR_PATH=/path/to/jline.jar`

## Basic Usage

```bash
# Analyze a model (uses MVA solver by default)
python line.py solve model.jsimg

# Choose solver: mva, fluid, ctmc, ssa, jmt, nc, mam
python line.py solve model.jsimg -s fluid

# Output formats: table (default), json, csv
python line.py solve model.jsimg -o json > results.json

# Analysis types: all, avg, sys
python line.py solve model.jsimg -a avg,sys

# List available options
python line.py list solvers
python line.py list analysis
```

## Solvers

| Solver | Use Case | Speed |
|--------|----------|-------|
| `mva` | General use (default) | Fast |
| `fluid` | Large models (100+ jobs) | Fast |
| `ctmc` | Exact results, small models | Slow |
| `ssa` | Complex distributions | Medium |
| `jmt` | Simulation validation | Slow |

## Output Metrics

| Metric | Description |
|--------|-------------|
| QLen | Average queue length (jobs) |
| Util | Server utilization (0-1) |
| RespT | Response time (seconds) |
| Tput | Throughput (jobs/sec) |
| ArvR | Arrival rate (jobs/sec) |

## Command Reference

| Option | Description |
|--------|-------------|
| `-s, --solver` | Solver algorithm |
| `-o, --output-format` | table, json, csv, raw |
| `-a, --analysis` | Analysis type(s), comma-separated |
| `-d, --seed` | Random seed for stochastic solvers |
| `-v, --verbose` | Verbose output |
| `-q, --quiet` | Suppress non-essential output |

## Creating Models

Create JSIMG files using the [JSIMgraph GUI](http://jmt.sf.net), or programmatically:

```python
from line_solver import *

model = Network('MyModel')
source = Source(model, 'Source')
queue = Queue(model, 'Queue', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')

job_class = OpenClass(model, 'Jobs')
source.setArrival(job_class, Exp(0.5))
queue.setService(job_class, Exp(1.0))

model.addLink(source, queue)
model.addLink(queue, sink)

solver = SolverMVA(model)
print(solver.avg_table())
```

## More Information

- Documentation: [line-solver.net](http://line-solver.net)
- Issues: [github.com/imperial-qore/line-solver](https://github.com/imperial-qore/line-solver/issues)
