# line-sage-rest: the LINE symbolic backend

SageMath behind a small JSON/HTTP service, so that all three LINE codebases
share one computer algebra engine and one normal form.

It exists because the symbolic capability was unequal and licence bound:

- the **JAR** has no CAS at all. `SolverCTMC.getSymbolicGenerator` only exploits
  that the generator is *linear* in the event symbols and stores one numeric
  coefficient matrix per event, so it could evaluate and print a symbolic
  generator but never solve, simplify or differentiate with it;
- **MATLAB** needs the Symbolic Math Toolbox, which is licence gated;
- **Python** uses sympy, exact but slow over the rational function field, and
  with no normal form shared with MATLAB.

The service closes all three gaps at once. It is a plain repackaging of
SageMath and adds no algebra of its own: results obtained through it must cite
SageMath, see `/api/v1/citation`.

## Running it

```bash
docker run -d --rm -p 8080:8080 imperialqore/line-sage-rest:latest
curl -s localhost:8080/api/v1/ready
```

Clients also start it themselves from a locally present image, on an ephemeral
port, when nothing is listening; see the resolution order below.

The image keeps the dual-mode entrypoint of the other `line-*-rest` services,
so it is a superset of the stock `sagemath/sagemath` image:

```bash
docker run --rm imperialqore/line-sage-rest sage -c 'print(factor(2^67-1))'
docker run --rm imperialqore/line-sage-rest sage -python /app/selftest.py
```

Build and publish with `../../upload-sage.sh`, which refuses to push an image
whose self-test or readiness check fails.

## Endpoints

| Method | Path | Request | Response |
|---|---|---|---|
| GET | `/api/v1/health` | - | `{status, timestamp}` |
| GET | `/api/v1/ready` | - | `{ready, checks}` (actually solves a chain) |
| GET | `/api/v1/info` | - | `{name, version, sage_version, endpoints, citation}` |
| GET | `/api/v1/citation` | - | upstream attribution |
| POST | `/api/v1/ctmc/solve` | `{Q, symbols, normalize}` | `{pi, num, den, nConnComp, connComp}` |
| POST | `/api/v1/ctmc/sensitivity` | `{Q, symbols, theta, reward}` | `{pi, dpi, Er, S, SS}` |
| POST | `/api/v1/simplify` | `{exprs, form}` | `{results}` |
| POST | `/api/v1/diff` | `{exprs, var, order}` | `{results}` |
| POST | `/api/v1/eval` | `{exprs, values}` | `{exact, values}` |
| POST | `/api/v1/fluid/odes` | `{rhs, vars, want}` | `{jacobian, latex, equilibria}` |

Every response carries `status`, `"ok"` or `"error"`; an error also carries a
machine-readable `code` (`parse`, `input`, `singular`, `absorbing`, `timeout`,
`compute`, `internal`). Errors are reported with HTTP 200 and that field, so a
client distinguishes a refused request from a transport failure.

## Protocol notes

- **Expressions are plain ASCII infix strings**, e.g. `2*x1 - 3*x2`. That is
  exactly what `symbolicGeneratorResult.getSymbolicEntry` (JAR) and
  `SolverCTMC.symbolicEntries` (MATLAB) emit. Both `^` and `**` are accepted;
  `^` is rewritten textually before parsing, because Python reads it as xor and
  would silently misparse `2*x^2` as `(2*x)^2`.
- **Numbers are read exactly.** A decimal literal is converted through its
  source text, so `0.1` is `1/10` and not the binary double nearest to it. This
  is why clients send coefficients as text.
- **Nothing is evaluated.** The parser is an explicit AST walker with a
  whitelist (`sqrt`, `exp`, `log`, `abs`, `min`, `max`), never `eval` or
  `sage_eval`, so the service is safe to expose on a port and cannot resolve a
  name the caller did not declare.
- **Results are in a primitive normal form**: reduced, integer coefficients
  with no common factor, positive leading coefficient in the denominator.
  Without it the same probability comes back scaled differently from state to
  state, because a rational constant is a unit in the fraction field and never
  enters a gcd.
- **Do not compare expression text across codebases.** Symbol numbering x1..xE
  follows event enumeration order, which is matched in practice but not
  structurally guaranteed, and the printed form depends on the Sage version.
  Substitute rates through `/api/v1/eval` and compare numbers.
- **Concurrency**: the server forks per request, so `signal.alarm` can enforce
  `timeout_s` and a runaway solve is killed with its process. Symbolic solves
  grow superpolynomially in the number of states.

## Backend resolution (identical in all three clients)

1. an explicit URL (`options.config.symbolic` in MATLAB and the JAR,
   `set_backend()` in Python);
2. the `LINE_SAGE_URL` environment variable;
3. a service already listening on 8085 or 8080, verified through
   `/api/v1/info` rather than by port: every `imperialqore/line-*-rest`
   service uses 8080, so a health probe alone would accept the LQNS one;
4. a container started by the client from a locally present image, on an
   ephemeral port;
5. nothing, in which case the caller falls back to its native CAS.

The default keeps the native engine (Symbolic Toolbox, sympy) so no existing
result changes silently; the JAR has no native engine and so always uses the
service.

## Clients

| Codebase | Entry point |
|---|---|
| JAR | `jline.api.sym.SymEngines.resolve()`, `SolverCTMC.getSymbolicSolution()` |
| MATLAB | `SAGE.m`, `@SolverCTMC/getSymbolicSolution.m`, `ctmc_solve` with `options.config.symbolic` |
| Python | `line_solver.api.sym`, `ctmc_solve` with `LINE_SYMBOLIC_BACKEND=sage` |

## Tests

- in-image: `sage -python /app/selftest.py` (algebra only, no HTTP)
- JAR: `jline.solvers.ctmc.SolverCTMCSymbolicSageTest` (tagged `remote`)
- MATLAB: `line-test.git/test/testsCTMC/test_symbolic_ctmc.m`, the `test_sage_*` cases
- Python: `line-test.git/test/testsCTMC/test_symbolic_ctmc.py`, `TestSymbolicSageBackend`

All of them skip cleanly when no backend can be resolved.
