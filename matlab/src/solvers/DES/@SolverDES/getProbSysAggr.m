function ProbSysAggr = getProbSysAggr(self)
% PROBSYSAGGR = GETPROBSYSAGGR() Returns joint aggregated system state probability via DES simulation
%
% @brief Estimates joint steady-state probability of the entire aggregated system state via DES simulation
%
% This method estimates the probability of observing the current system state
% (combined per-class job counts across all stateful nodes) using simulation.
% States are aggregated over service phases - only job counts per class matter.
%
% @param self SolverDES instance
%
% @return ProbSysAggr Estimated joint aggregated system state probability value
%         - 0 if state not seen during simulation (rare event)
%
% @note DES estimates probabilities from finite simulation. Accuracy improves with
%       more samples. For complex systems with large state spaces, the target state
%       may be rarely visited.
%
% @warning Simulation-based estimate - results vary between runs unless seed is set.
%
% @see getProbSys - Returns joint phase-detailed system state probabilities
% @see getProbAggr - Returns aggregated state probabilities for a single node
%
% Example:
% @code
% solver = SolverDES(model, 'samples', 100000, 'seed', 42);
%
% % Estimate probability of current joint aggregated system state
% prob_sys_aggr = solver.getProbSysAggr();
% fprintf('Estimated joint aggregated probability = %.6f\n', prob_sys_aggr);
% @endcode

if GlobalConstants.DummyMode
    ProbSysAggr = NaN;
    return
end

% Convert model to Java
jmodel = LINE2JLINE(self.model);

% Create Java DES options
joptions = jline.solvers.des.DESOptions();
joptions.samples = self.options.samples;
joptions.seed = self.options.seed;

% Create Java solver
jsolver = jline.solvers.des.SolverDES(jmodel, joptions);

% Call Java getProbSysAggr method and extract probability value
jresult = jsolver.getProbSysAggr();

% Extract the probability value from the ProbabilityResult object
ProbSysAggr = jresult.probability;

end
