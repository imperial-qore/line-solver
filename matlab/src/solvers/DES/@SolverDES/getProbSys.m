function ProbSys = getProbSys(self)
% PROBSYS = GETPROBSYS() Returns joint system state probability via DES simulation
%
% @brief Estimates joint steady-state probability of the entire system state via DES simulation
%
% This method estimates the probability of observing the current system state
% (combined state across all stateful nodes) using simulation-based estimation.
% States include phase information from service distributions.
%
% @param self SolverDES instance
%
% @return ProbSys Estimated joint system state probability value
%         - 0 if state not seen during simulation (rare event)
%
% @note DES estimates probabilities from finite simulation. Accuracy improves with
%       more samples. For complex systems with large state spaces, the target state
%       may be rarely visited.
%
% @warning Simulation-based estimate - results vary between runs unless seed is set.
%
% @see getProbSysAggr - Returns joint aggregated system state probabilities
% @see getProb - Returns marginal state probabilities for a single node
%
% Example:
% @code
% solver = SolverDES(model, 'samples', 100000, 'seed', 42);
%
% % Estimate probability of current joint system state
% prob_sys = solver.getProbSys();
% fprintf('Estimated joint system probability = %.6f\n', prob_sys);
% @endcode

if GlobalConstants.DummyMode
    ProbSys = NaN;
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

% Call Java getProbSys method and extract probability value
jresult = jsolver.getProbSys();

% Extract the probability value from the ProbabilityResult object
ProbSys = jresult.probability;

end
