function self = initFromSolver(self, initSolver)
% SELF = INITFROMSOLVER(INITSOLVER)
%
% Warm-start the simulation from the steady-state solution of an auxiliary
% solver.
%
% If the auxiliary solver is a SolverCTMC, the exact stationary distribution
% over the aggregate state space is computed and the initial state is set to
% the mode of that distribution (the most probable aggregate state). For any
% other network solver, the steady-state mean queue lengths are used instead
% and rounded to an integer placement that conserves each closed-class
% population (see NetworkSolver.warmStartPlacement).
%
% Since the simulation starts (approximately) in steady state, the transient
% removal filter is disabled (tranfilter = 'fixed' with warmupfrac = 0), so
% every simulated sample contributes to the estimators.
%
% Example:
% @code
% solver = SolverLDES(model, SolverMVA(model), 'samples', 50000);
% solver.getAvg();  % simulation starts from the MVA steady-state placement
% @endcode

sn = self.model.getStruct(true);
M = sn.nstations;
K = sn.nclasses;

placement = NetworkSolver.warmStartPlacement(initSolver, sn);

% LDES consumes the placement through options.init_sol (station-major vector)
% rather than the model initial state.
self.options.init_sol = reshape(placement', 1, M * K);
self.options.config.tranfilter = 'fixed';
self.options.config.warmupfrac = 0;
end
