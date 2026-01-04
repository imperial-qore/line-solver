# Basic Examples - Kotlin

This directory contains basic examples demonstrating core LINE solver functionality using Kotlin notebooks.

## Directory Structure

### Closed Queueing Networks (`closedQN/`)
Examples of closed queueing networks where jobs circulate within the system.
- `cqn_oneline.ipynb` - One-line creation of cyclic networks using `Network.cyclicPsInf()`
- `cqn_twoqueues.ipynb` - Simple closed network with probabilistic routing
- `cqn_twoclass_erl.ipynb` - Multi-class network with Erlang service distributions

### Open Queueing Networks (`openQN/`)
Examples of open queueing networks with external arrivals and departures.
- `oqn_basic.ipynb` - Basic open network with hyperexponential service and multiple solver comparison

### Mixed Queueing Networks (`mixedQN/`)
Examples combining both open and closed job classes in the same network.
- `mqn_basic.ipynb` - Mixed network with different service distributions per class

### Cache Models (`cacheModel/`)
Examples of cache replacement strategies and performance analysis.
- `cache_replc_routing.ipynb` - Cache with state-dependent routing and hit/miss analysis
- `cache_replc_rr.ipynb` - Cache with Round Robin replacement strategy
- `cache_replc_fifo.ipynb` - Cache with FIFO replacement strategy

### Fork-Join Networks (`forkJoin/`)
Examples of parallel processing with fork-join topologies.
- `fj_basic_open.ipynb` - Basic fork-join with parallel queues and synchronization

### Layered Models (`layeredModel/`)
Examples of Layered Queueing Networks (LQN) for software system modeling.
- `lqn_basic.ipynb` - Multi-tier application with processors, tasks, and synchronous calls

## Common Patterns

### Model Creation
```kotlin
val model = Network("ModelName")
```

### Node Types
- **Source/Sink**: Entry and exit points for open classes
- **Queue**: Service stations with various scheduling strategies
- **Delay**: Infinite server stations (think time, pure delays)
- **Router**: Routing decisions
- **Fork/Join**: Parallel processing synchronization
- **Cache**: State-dependent routing for cache modeling

### Job Classes
- **OpenClass**: Unlimited population from external arrivals
- **ClosedClass**: Fixed population circulating in the system

### Service Processes
- **Exp(rate)**: Exponential distribution
- **Erlang(mean, order)**: Erlang distribution
- **HyperExp(p, λ1, λ2)**: Hyperexponential distribution

### Solvers
- **JMT**: Simulation-based, handles complex models
- **MVA**: Mean Value Analysis, exact for product-form networks
- **CTMC**: Continuous Time Markov Chain, exact but computationally intensive
- **FLD**: Fluid approximation, fast for heavily loaded systems
- **SSA**: Stochastic State-space Analysis
- **NC**: Normalizing Constant method

## Getting Started

1. Start with `closedQN/cqn_oneline.ipynb` for a simple closed network
2. Try `openQN/oqn_basic.ipynb` to see solver comparisons
3. Explore `mixedQN/mqn_basic.ipynb` for combined workloads
4. Advanced features in `forkJoin/` and `layeredModel/`

## Performance Metrics

All examples compute standard performance metrics:
- **QLen**: Average queue length
- **Util**: Station utilization
- **RespT**: Mean response time
- **ResidT**: Mean residence time
- **Tput**: Throughput
- **ArvR**: Arrival rate

Access results using:
```kotlin
val avgTable = solver.avgTable
avgTable.print()
```