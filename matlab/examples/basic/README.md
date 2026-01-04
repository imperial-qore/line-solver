# LINE Solver MATLAB Examples - Basic

This directory contains MATLAB examples demonstrating basic features of the LINE solver. Examples are organized into subdirectories by category.

## Directory Structure

- **cacheModel/** - Cache modeling examples with hit/miss behavior
- **classSwitching/** - Class switching scenarios  
- **closedQN/** - Closed queueing network models
- **forkJoin/** - Fork-join queueing models
- **layeredModel/** - Layered queueing network (LQN) examples
- **mixedQN/** - Mixed open/closed queueing models
- **openQN/** - Open queueing network models
- **prioModel/** - Priority queueing examples
- **stochPetriNet/** - Stochastic Petri net examples

## Running Examples

Each subdirectory contains MATLAB scripts that can be run individually:

```matlab
% Navigate to examples directory
cd matlab/examples/basic

% Run a specific example
run('closedQN/example_closedModel_1.m')
```

## Notes

- Examples demonstrate various LINE solver features including different solver algorithms (MVA, NC, LQNS, etc.)
- Some examples may require additional dependencies or external tools (e.g., LQNS solver)
- MATLAB implementation provides the most complete feature set