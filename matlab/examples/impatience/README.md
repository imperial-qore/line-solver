# Impatience Examples

This directory contains examples demonstrating customer impatience (reneging/abandonment) in queueing networks using LINE's JMT solver.

## What is Impatience?

Customer impatience allows jobs to abandon (renege from) queues after waiting for a certain amount of time. The patience time is drawn from a configurable probability distribution. This is useful for modeling:

- Call centers with customer abandonments
- Web services with timeout mechanisms
- Healthcare systems with patient reneging
- Any system where waiting customers may leave before service

## Example Files

### `example_mm1_exp_impatience.m`
**Purpose:** Simple M/M/1 queue with exponential patience

Demonstrates the basics of configuring impatience at a single queue with a single job class. Uses exponential patience distribution.

**Key concepts:**
- Queue-specific patience configuration using `queue.setPatience(jobclass, distribution)`
- Exponential patience distribution
- Impact of reneging on queue length and response time

### `example_multiclass_mixed_impatience.m`
**Purpose:** Multi-class system with different patience distributions

Shows how different job classes can have different patience behaviors, and how patience can vary across queues in the network.

**Key concepts:**
- Different distributions per class (Exponential, Deterministic, Erlang)
- Multiple queues with different patience parameters
- Comparing reneging rates across classes

### `example_global_patience.m`
**Purpose:** Global class-level patience with per-queue overrides

Demonstrates setting patience at the job class level (applies to all queues) and selectively overriding it at specific queues.

**Key concepts:**
- Global patience: `jobclass.setPatience(distribution)`
- Per-queue override: `queue.setPatience(jobclass, distribution)`
- Precedence: queue-specific settings override global class settings

### `example_no_impatience.m`
**Purpose:** Baseline verification without impatience

A standard queueing model without any impatience configured, serving as a baseline and verifying backward compatibility.

**Key concepts:**
- Backward compatibility verification
- Standard queueing behavior
- XML generation produces `null` for impatience when not configured

## Supported Distributions

The following distributions are supported for patience/impatience:

- **Exponential** (`Exp`) - Memoryless, constant hazard rate
- **Deterministic** (`Det`) - Fixed timeout
- **Erlang** (`Erlang`) - Less variable than exponential
- **Hyperexponential** (`HyperExp`) - More variable than exponential
- **Gamma** (`Gamma`) - Flexible shape
- **Pareto** (`Pareto`) - Heavy-tailed
- **Weibull** (`Weibull`) - Increasing/decreasing hazard rate
- **Lognormal** (`Lognormal`) - Right-skewed
- **Uniform** (`Uniform`) - Constant probability over interval
- **Phase-type** (`PH`, `APH`, `Coxian`) - General distributions

**Not supported:** BMAP, MAP, MMPP2 and other modulated processes

## Usage

To run an example, start LINE and execute the script:

```matlab
lineStart
cd examples/impatience
example_mm1_exp_impatience
```

## Configuration Methods

### Queue-Specific Patience
Set patience for a specific class at a specific queue:
```matlab
queue.setPatience(jobclass, Exp(0.5));  % Mean patience = 2.0
```

### Global Class-Level Patience
Set patience for a class that applies to all queues:
```matlab
jobclass.setPatience(Det(5.0));  % Fixed timeout = 5.0 everywhere
```

### Precedence
If both global and queue-specific patience are set, the queue-specific setting takes precedence.

## Solver Support

**Currently supported:** JMT solver only

The impatience feature leverages JMT's native impatience support (`jmt.engine.NetStrategies.ImpatienceStrategies.Impatience`).

## Additional Resources

- See LINE manual for detailed API documentation
- Consult JMT documentation for underlying simulation mechanics
- Check `doc/latex/manual.tex` for theoretical background

## Copyright

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
