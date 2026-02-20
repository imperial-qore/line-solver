# LINE Solver MATLAB Examples - Advanced

This directory contains MATLAB examples demonstrating advanced features of the LINE solver. Examples are organized into subdirectories by category.

## Directory Structure

- **agentModels/** - Agent-based modeling examples
- **cdfRespT/** - Cumulative distribution function of response times
- **cyclicPolling/** - Cyclic polling system examples
- **example_custom_solver/** - Custom solver implementation example
- **fcRegion/** - Feasible capacity region analysis
- **initState/** - Initial state configuration examples
- **layeredCQ/** - Layered networks with caching
- **loadDependent/** - Load-dependent service examples
- **randomEnv/** - Random environment models
- **rewardModel/** - Reward models and reinforcement learning routing examples
- **stateDepRouting/** - State-dependent routing
- **stateProbabilities/** - State probability analysis examples
- **switchoverTimes/** - Switchover time modeling

## Running Examples

Each subdirectory contains MATLAB scripts that can be run individually:

```matlab
% Navigate to examples directory
cd matlab/examples/advanced

% Run a specific example
run('cdfRespT/cdf_respt_closed.m')
```

## Notes

- Examples demonstrate various LINE solver features including different solver algorithms (MVA, NC, LQNS, etc.)
- Some examples may require additional dependencies or external tools (e.g., LQNS solver)
- MATLAB implementation provides the most complete feature set
- Custom solver example shows how to extend LINE with new algorithms