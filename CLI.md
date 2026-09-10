# LINE CLI

`line-cli.py` — a single-file command-line front end for the LINE queueing
network solver. Pure Python 3.9+ standard library, no dependencies; wraps `jline.jar`.

```bash
python line-cli.py <command> [options]   # or ./line-cli.py after chmod +x
```

This is one of five LINE command-line front ends; see "The other front ends"
below for which one to reach for and what each covers.

## solve

```bash
python line-cli.py solve model.json               # default MVA solver
python line-cli.py solve model.json -s fluid      # pick a solver
python line-cli.py solve model.json -o json       # output: table (default) | json | csv
python line-cli.py solve model.json -a avg        # analysis: avg | sys | ...
python line-cli.py solve model.json -s jmt -d 12345     # seed for stochastic solvers
python line-cli.py solve model.json -s ssa --samples 1e6   # simulation run length
cat model.jsimg | python line-cli.py solve -i jsimg -s mva   # read from stdin (external format)
```

Solvers (`list solvers`): `mva` (default), `auto`, `ag`, `ba`, `ctmc`, `env`,
`fld`/`fluid`, `jmt`, `ldes` (alias `des`), `ln`, `ln.mva`, `ln.nc`,
`ln.comom`, `lqns`, `mam`, `nc`, `qns`/`qnsolver`, `ssa`, `uq`.
`env` reads an Environment model and takes no other token; `uq` needs
`--uq-solver` to name the engine it runs at each design point.

Input formats (`list formats`): `json`, `jsim`, `jsimg`, `jsimw`, `lqnx`,
`xml`, `pnml`. `json` is LINE's own portable model format
([`doc/line-model.schema.json`](doc/line-model.schema.json), also read and
written through the library API `save_model`/`load_model`); the others are
external model file types (JMT's JSIMG/JSIMW/JSIM, LQNS's LQNX, ISO/IEC 15909-2
PNML) that can also be passed to the CLI.

Analyses (`list analysis`) cover the average tables (`avg`, `sys`, `chain`,
`node`, `nodechain`, `stage`), the cache and station-class tables (`cache`,
`item`, `orbit`, `loss`, `region-loss`, `deadline`), `normconst`, `busyperiod`,
`sens`, the distribution and transient families (`cdf-respt`, `cdf-passt`,
`perct-respt`, `tran-*`), the probability and sampling families (`prob*`,
`sample*`), the CTMC rewards (`reward*`) and the solver-internal reports
(`generator`, `statevec`, `moments`, `interval`).

Numeric controls are forwarded verbatim to the JAR: `--samples`, `--cutoff`,
`--timespan` (alias `--tspan`), `--timestep`, `--method`, `--tol`,
`--iter_tol`, `--iter_max`, `--multiserver`, `--warmupfrac`, `--stage-solver`,
`--uq-solver`, `--busyperiod`, `--busyperiod-subnet`, `--sens-method`,
`--sens-scheme`, `--sens-step`.

## info

```bash
python line-cli.py info    # system info, JAR status, solver compatibility
```

## server / rest

```bash
python line-cli.py server                  # WebSocket server on port 5863
python line-cli.py server -H 0.0.0.0 -p 8080
```

`rest` is an alias of `server` (execs `java -jar jline.jar -p <port>`): it starts
the WebSocket server (default port 8080) and does **not** serve HTTP. The HTTP
REST API is the separate `io/rest-api/` Maven module, not reachable through this CLI:

```bash
cd io/rest-api && mvn clean package
java -cp target/line-rest.jar:common/jline.jar jline.rest.LineRestServer --port 8080
```

It serves JSON under `/api/v1` (`POST /api/v1/models/solve`, `GET /api/v1/health`,
plus job/analysis/metrics routes); see `io/rest-api/README.md` and `io/rest-api/openapi.yaml`.

## Configuration

Optional `~/.config/line-cli/config.yaml` (needs PyYAML) overrides defaults:
`jar_path`, `java_path`, `default_solver`, `default_output_format`, and a `server:`
block with `host`/`port`. Otherwise the JAR is auto-detected in `common/jline.jar`
or taken from the `LINE_JAR_PATH` env var.

## Requirements

Python 3.9+, Java 8+, `jline.jar`. PyYAML is optional (config file only).

## The other front ends

| Front end | Invocation | Covers |
|---|---|---|
| `line-cli.py` (this one) | `python line-cli.py solve ...` | a thin wrapper over the JAR; everything the JAR serves, plus `csv` output and the `list`/`info` subcommands |
| JAR `LineCLI` | `java -jar common/jline.jar -s ... -a ...` | the reference vocabulary this wrapper forwards to |
| C++ `line-cli` | `common/line-cli -f model.json -s ...` | a superset: `--arith` multiprecision, `--api`, sensitivity/UQ knobs, and a dozen analyses with no JAR spelling. `--help-all` prints the full reference |
| native python | `python3 -m line_solver.cli -f model.json -s ...` | the no-JVM twin; same formats and the same analysis vocabulary, with `pickle`/`mat` output on top |
| MATLAB `linemcr` | `linemcr('-f','model.json','-s','mva')` | the MCR/Docker front end: the average tables (`avg`, `sys`, `chain`, `node`, `nodechain`, `cache`, `item`) only |

Two deliberate differences between the JAR and the C++ CLI, documented because
neither can be removed without breaking a bridge:

- **`-n`/`--node` and `-c`/`--class` are 0-based in the JAR (and in both Python
  CLIs) and 1-based in the C++ `line-cli`**, which indexes stations as MATLAB
  does. `python/line_solver/solvers/cpp_dispatch.py` and `jar_dispatch.py` each
  convert for their own CLI. Both spellings parse in both, so a command line is
  portable; the *number* is not.
- The JAR spells the default verbosity `normal` and the C++ CLI `standard`.
  Both now accept both.
