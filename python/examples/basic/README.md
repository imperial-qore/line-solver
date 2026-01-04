# LINE Solver Python Examples

This directory contains Jupyter notebook examples demonstrating various features of the LINE solver. Examples are organized into subdirectories by category.

## Directory Structure

- **cacheModel/** - Cache modeling examples with hit/miss behavior
- **classSwitching/** - Class switching scenarios
- **closedModel/** - Closed queueing network models
- **forkJoin/** - Fork-join queueing models
- **layeredModel/** - Layered queueing network (LQN) examples
- **mixedModel/** - Mixed open/closed queueing models
- **openModel/** - Open queueing network models
- **prioModel/** - Priority queueing examples
- **stochPetriNet/** - Stochastic Petri net examples

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