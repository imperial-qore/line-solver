# LINE Solver Python Examples

This directory contains Jupyter notebook examples demonstrating various features of the LINE solver. Examples are organized into subdirectories by category.

## Directory Structure

- **cdfRespT/** - Cumulative distribution function of response times
- **cyclicPolling/** - Cyclic polling system examples
- **initState/** - Initial state configuration examples
- **layeredCacheQueueing/** - Layered networks with caching
- **loadDependent/** - Load-dependent service examples
- **randomEnvironment/** - Random environment models
- **stateDependentRouting/** - State-dependent routing
- **switchoverTimes/** - Switchover time modeling

## Running Examples

Each subdirectory contains Jupyter notebooks that can be run individually:

```bash
# Navigate to advanced directory
cd python/advanced

# Run a specific example
jupyter notebook closedQN/example_closedModel_1.ipynb
```

## Notes

- Examples demonstrate various LINE solver features including different solver algorithms (MVA, NC, LQNS, etc.)
- Some examples may require additional dependencies or external tools (e.g., LQNS solver)
- Fallback implementations are provided where advanced features are not yet fully supported in the Python version