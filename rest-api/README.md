# LINE REST API

REST API for solving queueing network models using LINE Solver.

## Prerequisites

Place `jline.jar` in the `common/` directory before building.

## Quick Start

```bash
# Build
mvn clean package

# Run
java -cp target/line-rest.jar:common/jline.jar jline.rest.LineRestServer

# With security options
java -cp target/line-rest.jar:common/jline.jar jline.rest.LineRestServer \
  --port 9090 --api-keys key1,key2 --rate-limit 100
```

## Endpoints

### Health & Info

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/health` | GET | Liveness probe |
| `/api/v1/ready` | GET | Readiness probe |
| `/api/v1/info` | GET | Server info |
| `/api/v1/solvers` | GET | List solvers |

### Models

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/models/solve` | POST | Solve model (sync) |
| `/api/v1/models/solve/async` | POST | Solve model (async) |
| `/api/v1/models/validate` | POST | Validate model |
| `/api/v1/models/convert` | POST | Convert format |
| `/api/v1/models/calibrate` | POST | Calibrate from metrics |

### Jobs

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/jobs` | GET | List jobs |
| `/api/v1/jobs/:id` | GET | Get job status |
| `/api/v1/jobs/:id` | DELETE | Cancel job |
| `/api/v1/jobs/:id/stream` | GET | Stream progress (SSE) |

### Analysis

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/analysis/whatif` | POST | Parameter sweep |
| `/api/v1/analysis/sensitivity` | POST | Sensitivity analysis |
| `/api/v1/analysis/bottleneck` | POST | Bottleneck detection |

### SRE/Metrics

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/metrics` | GET | Prometheus metrics |
| `/api/v1/metrics/json` | GET | JSON metrics |
| `/api/v1/traces/import` | POST | Import traces |

## Request Examples

### Solve Model

```bash
curl -X POST http://localhost:8080/api/v1/models/solve \
  -H "Content-Type: application/json" \
  -d '{
    "model": {"format": "jsimg", "content": "...", "base64": true},
    "solver": "mva",
    "analysis": "all"
  }'
```

### Async Solve

```bash
# Submit job
curl -X POST http://localhost:8080/api/v1/models/solve/async \
  -d '{"model": {...}, "solver": "jmt"}'

# Check status
curl http://localhost:8080/api/v1/jobs/{jobId}
```

## Security

### API Key Authentication

```bash
# Enable via CLI or environment
java -cp ... jline.rest.LineRestServer --api-keys key1,key2
# or
export LINE_API_KEYS="key1,key2"

# Usage
curl -H "X-API-Key: key1" http://localhost:8080/api/v1/solvers
```

Excluded paths: `/health`, `/ready`, `/metrics`

### Rate Limiting

```bash
# 100 requests per 60-second window
java -cp ... jline.rest.LineRestServer --rate-limit 100
# or
export LINE_RATE_LIMIT=100
export LINE_RATE_WINDOW=60
```

## Docker

```bash
# Build
docker build -t line-rest .

# Run
docker run -p 8080:8080 line-rest

# With options
docker run -p 8080:8080 \
  -e LINE_API_KEYS="key1,key2" \
  -e LINE_RATE_LIMIT=100 \
  line-rest
```

## Kubernetes

```bash
cd kubernetes
kubectl apply -f configmap.yaml -f deployment.yaml -f service.yaml

# Create API key secret
kubectl create secret generic line-solver-secrets \
  --from-literal=api-keys="key1,key2"
```

## Solvers

| Solver | Type | Description |
|--------|------|-------------|
| `mva` | Analytical | Mean Value Analysis |
| `jmt` | Simulation | Java Modelling Tools |
| `ssa` | Analytical | Stochastic State-space |
| `ctmc` | Analytical | CTMC |
| `fluid` | Analytical | Fluid approximation |
| `nc` | Analytical | Normalizing constant |
| `des` | Simulation | Discrete Event Simulation |
| `ln` | Analytical | Layered Networks |
| `lqns` | Analytical | LQN Solver |

## OpenAPI

Full API spec: `openapi.yaml`
