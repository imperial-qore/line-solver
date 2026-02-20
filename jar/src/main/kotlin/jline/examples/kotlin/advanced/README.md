# Advanced Examples - Kotlin

This directory contains advanced examples demonstrating specialized LINE solver capabilities using Kotlin notebooks.

## Directory Structure

### CDF Response Time Analysis (`cdfRespT/`)
Examples of cumulative distribution function analysis for response times.
- `cdf_respt_closed.ipynb` - CDF analysis for closed networks with statistical moment calculation

### Cyclic Polling (`cyclicPolling/`)
Examples of polling systems with different service disciplines.
- Coming soon: Exhaustive, gated, and k-limited polling examples

### Initial State Analysis (`initState/`)
Examples analyzing systems starting from specific initial states.
- Coming soon: FCFS and PS systems with non-equilibrium initial conditions

### Layered Cache Queueing (`layeredCQ/`)
Examples of multi-level cache hierarchies.
- Coming soon: Single-host and multi-host cache examples

### Load-Dependent Service (`loadDependent/`)
Examples with service rates depending on queue length or system state.
- Coming soon: Multi-server systems with load-dependent service

### Random Environment (`randomEnv/`)
Examples of queueing systems in randomly changing environments.
- `example_random_env_basics.ipynb` - Basic random environment model with fast/slow server modes
- `renv_fourstages_repairmen.ipynb` - Multi-stage repairmen model with random failures
- `renv_threestages_repairmen.ipynb` - Three-stage repair process with environmental effects
- `renv_twostages_repairmen.ipynb` - Two-stage repair model
- `renv_genqn.ipynb` - General closed queueing network in random environment

### State-Dependent Routing (`stateDepRouting/`)
Examples with routing decisions based on system state.
- Coming soon: Open and closed networks with state-dependent routing

### State Probabilities (`stateProbabilities/`)
Examples for computing steady-state probability distributions.
- Coming soon: Aggregate and detailed state probability analysis

### Switchover Times (`switchoverTimes/`)
Examples of systems with setup/switchover times between different service modes.
- Coming soon: Basic switchover time examples

## Key Advanced Features

### Distribution Analysis
Beyond mean values, LINE can compute:
- **Cumulative Distribution Functions (CDFs)**: Full probability distributions
- **Higher-Order Moments**: Variance, skewness, kurtosis
- **Percentiles**: 95th, 99th percentile response times

### Specialized Queueing Systems
- **Polling Systems**: Cyclic service with various disciplines
- **Cache Models**: Multi-level cache hierarchies with replacement policies
- **Random Environments**: Time-varying service parameters
- **Load-Dependent Service**: Service rates that change with queue length

### Advanced Analysis Techniques
- **Transient Analysis**: Time-dependent behavior from initial states
- **State Space Analysis**: Detailed Markov chain state probabilities
- **Sensitivity Analysis**: Impact of parameter changes
- **Optimization**: Finding optimal system configurations

## CDF Analysis Example

```kotlin
// Get response time CDF from solver
val cdf = solver.cdfRespT

// Calculate mean from CDF by numerical integration
var meanRT = 0.0
for (k in 1 until cdf[i][c].size(0)) {
    val dt = cdf[i][c].get(k, 0) - cdf[i][c].get(k-1, 0)
    val rt = cdf[i][c].get(k, 1)
    meanRT += dt * rt
}
```

## State Probability Analysis

```kotlin
// Get state space and probabilities
val stateSpace = solver.stateSpace
val probabilities = solver.stateProbabilities

// Analyze specific states
for (state in stateSpace) {
    println("State: ${state}, Probability: ${probabilities[state]}")
}
```

## Performance vs. Accuracy Tradeoffs

Different solvers offer various tradeoffs:

| Solver | Accuracy | Speed | Features |
|--------|----------|--------|----------|
| CTMC | Exact | Slow | Full state space |
| JMT | High | Medium | General distributions |
| Fluid | Medium | Fast | Heavy traffic approx |
| SSA | High | Medium | Non-product form |

## Getting Started with Advanced Examples

1. Begin with `cdfRespT/cdf_respt_closed.ipynb` for distribution analysis
2. Explore state-dependent features in specialized directories
3. Use transient analysis for time-dependent behavior
4. Apply optimization techniques for system design

## Common Advanced Patterns

### CDF Collection
```kotlin
val solver = JMT(model, "seed", 23000)
val cdf = solver.cdfRespT
// Process CDF data for statistical analysis
```

### State Space Analysis
```kotlin
val solver = CTMC(model)
val stateSpace = solver.stateSpace.stateSpace
val generator = solver.generator.infGen
```

### Transient Analysis
```kotlin
val solver = FLD(model)
solver.transient(timepoints)
```

These advanced examples demonstrate LINE's capabilities beyond basic steady-state mean value analysis, enabling detailed performance characterization and optimization of complex queueing systems.