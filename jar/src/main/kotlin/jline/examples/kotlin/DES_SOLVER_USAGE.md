# DES Solver Usage in Kotlin Notebooks

## Overview

The DES (Discrete Event Simulation) solver is available in the JAR implementation but is not demonstrated in most Kotlin notebooks. This document explains how to use the DES solver in Kotlin examples.

## Important Notes

1. **DES solver is JAR-only**: The DES solver is implemented in `jar/src/main/kotlin/jline/solvers/des/` using the SSJ library
2. **Not available in MATLAB**: DES solver doesn't have a MATLAB equivalent
3. **Supported features**:
   - Open queueing networks
   - FCFS queues
   - Delay nodes
   - Multiserver (M/M/c)
   - Multiclass workloads
   - Distributions: Exponential, Erlang, HyperExp, PH, APH, Coxian

## Usage Pattern

```kotlin
import jline.solvers.des.*

// Create your model (Network, queues, etc.)
val model = Network("MyModel")
// ... model setup ...

// Create DES solver
val desSolver = DES(
    model,
    "seed", 23000,           // Random seed for reproducibility
    "verbose", true,          // Show progress
    "samples", 10000,         // Number of samples to collect
    "maxEvents", 1000000      // Maximum simulation events
)

// Run solver
val avgTable = desSolver.getAvgTable()
avgTable.print()
```

## Example: M/M/1 Queue with DES

```kotlin
import jline.*
import jline.lang.*
import jline.lang.nodes.*
import jline.lang.processes.*
import jline.lang.constant.*
import jline.solvers.des.*

val model = Network("MM1-DES")

val source = Source(model, "Source")
val queue = Queue(model, "Queue", SchedStrategy.FCFS)
val sink = Sink(model, "Sink")

val jobclass = OpenClass(model, "Class1")
source.setArrival(jobclass, Exp(1.0))
queue.setService(jobclass, Exp(2.0))

model.link(Network.serialRouting(source, queue, sink))

// Use DES solver
val solver = DES(model, "seed", 23000, "samples", 100000)
val results = solver.getAvgTable()
results.print()
```

## Comparison with Other Solvers

DES is particularly useful for:
- **Validation**: Comparing analytical results with simulation
- **Complex scenarios**: When analytical solvers don't support certain features
- **Transient analysis**: Studying time-dependent behavior (future feature)

## Why DES is Not in Most Notebooks

1. **Performance**: DES is slower than analytical solvers for simple models
2. **Stochastic variability**: Results vary between runs (unless seed is fixed)
3. **Analytical focus**: Most examples demonstrate analytical solution methods
4. **MATLAB parity**: Since DES doesn't exist in MATLAB, it's not in ported examples

## Adding DES to Existing Notebooks

To add DES solver to an existing notebook, simply add it to the solver array:

```kotlin
val solvers = arrayOf(
    MVA(model),
    NC(model),
    DES(model, "seed", 23000, "samples", 10000)  // Add this line
)

for (solver in solvers) {
    val avgTable = solver.getAvgTable()
    avgTable.print()
}
```

## See Also

- DES implementation: `jar/src/main/kotlin/jline/solvers/des/`
- DES tests: `jar/src/test/java/jline/examples/basic/OpenExamplesTest.java`
- SSJ library documentation: http://simul.iro.umontreal.ca/ssj/
