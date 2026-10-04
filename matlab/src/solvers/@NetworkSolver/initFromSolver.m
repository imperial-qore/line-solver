function self = initFromSolver(self, initSolver)
% SELF = INITFROMSOLVER(INITSOLVER)
%
% Warm-start the solver from the steady-state solution of an auxiliary solver.
%
% The auxiliary solver is used to find the steady-state distribution of the
% model, and that distribution decides an integer job placement (see
% NetworkSolver.warmStartPlacement): with SolverCTMC the placement is the mode
% of the exact stationary distribution over the aggregate state space, with
% any other solver the rounded steady-state mean queue lengths conserving each
% closed-class population.
%
% The placement is applied as the model initial state via initFromMarginal,
% which the state-driven solvers honor: SolverFLD starts the ODE integration
% from it, SolverSSA starts the simulated trajectory from it, and SolverJMT
% preloads the stations with it. Note that this modifies the initial state of
% the model object shared with any other solver instance.
%
% Example:
% @code
% solver = SolverSSA(model, SolverMVA(model), 'samples', 50000);
% @endcode

sn = self.model.getStruct(true);
placement = NetworkSolver.warmStartPlacement(initSolver, sn);
self.model.initFromMarginal(placement);
end
