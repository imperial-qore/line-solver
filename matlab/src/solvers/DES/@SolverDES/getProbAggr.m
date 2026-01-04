function ProbAggr = getProbAggr(self, node, state)
% PROBAGGR = GETPROBAGGR(NODE, STATE) Returns aggregated state probability via DES simulation
%
% @brief Estimates steady-state aggregated (per-class) state probabilities via DES simulation
%
% This method estimates the probability of observing a specific per-class job distribution
% (e.g., [2 jobs of class 1, 1 job of class 2]) at a station. States are aggregated over
% service phases - only the number of jobs per class matters.
%
% @param self SolverDES instance
% @param node Queue or node object (or node index as integer)
% @param state (optional) Aggregated state specification as vector (per-class job counts).
%             If omitted, uses current network state aggregated over phases.
%
% @return ProbAggr Estimated aggregated state probability value
%         - 0 if state not seen during simulation (rare event)
%
% @note DES estimates probabilities from finite simulation. Accuracy improves with
%       more samples.
%
% @warning Simulation-based estimate - results vary between runs unless seed is set.
%
% @see getProb - Returns phase-detailed state probabilities
% @see getProbSysAggr - Returns joint aggregated system state probabilities
%
% Example:
% @code
% solver = SolverDES(model, 'samples', 50000, 'seed', 42);
% queue1 = model.nodes{1};
%
% % Estimate probability of having 3 class-1 jobs and 2 class-2 jobs
% prob_aggr = solver.getProbAggr(queue1, [3, 2]);
% fprintf('Estimated Pr[n1=3, n2=2] = %.4f\n', prob_aggr);
% @endcode

if GlobalConstants.DummyMode
    ProbAggr = NaN;
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

% Get stateful node from model
sn = self.model.getStruct;
if isa(node, 'Node')
    nodeIdx = node.index;
else
    nodeIdx = node;
end

% Find the corresponding Java node
statefulNodes = jmodel.getStatefulNodes();
isf = sn.nodeToStateful(nodeIdx);
jnode = statefulNodes.get(isf - 1);  % Java is 0-indexed

% Convert state to Java Matrix if provided
if nargin < 3 || isempty(state)
    jstate = [];
else
    jstate = JLINE.to_jline_matrix(state(:)');
end

% Call Java getProbAggr method and extract probability value
if isempty(jstate)
    jresult = jsolver.getProbAggr(jnode);
else
    jresult = jsolver.getProbAggr(jnode, jstate);
end

% Extract the probability value from the ProbabilityResult object
ProbAggr = jresult.probability;

end
