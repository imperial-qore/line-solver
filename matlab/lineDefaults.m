function options = lineDefaults(solverName)
if nargin < 1
    solverName = 'Solver'; % global options unless overridden by a solver
end
options = SolverOptions(solverName);
end
