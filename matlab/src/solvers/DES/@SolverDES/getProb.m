function Prob = getProb(self, node, state)
% PROB = GETPROB(NODE, STATE) Returns state probability estimated via DES simulation
%
% @brief Estimates steady-state marginal state probabilities at a node via DES simulation
%
% This method estimates state probabilities by running a discrete-event simulation
% and computing the fraction of total simulation time spent in each state.
%
% @param self SolverDES instance
% @param node Queue or node object (or node index as integer)
% @param state (optional) State specification as vector. If provided, returns
%             probability of this specific state. If omitted, uses current network state.
%
% @return Prob Estimated state probability value
%         - 0 if state not seen during simulation (rare event)
%
% @note DES estimates probabilities from finite simulation. Accuracy improves with
%       more samples. If the specified state is not observed during simulation,
%       probability 0 is returned.
%
% @warning Simulation-based estimate - results vary between runs unless seed is set.
%          Use 'seed' option in constructor for reproducible results.
%
% @see getProbAggr - Returns state probabilities aggregated over phases
% @see getProbSys - Returns joint system state probabilities
%
% Example:
% @code
% solver = SolverDES(model, 'samples', 50000, 'seed', 42);
% queue1 = model.nodes{1};
%
% % Estimate probability of state [2 jobs of class 0, 1 of class 1]
% prob_state = solver.getProb(queue1, [2, 1]);
% fprintf('Estimated Pr[state=[2,1]] = %.4f\n', prob_state);
% @endcode

if GlobalConstants.DummyMode
    Prob = NaN;
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

% Call Java getProb method and extract probability value
if isempty(jstate)
    jresult = jsolver.getProb(jnode);
else
    jresult = jsolver.getProb(jnode, jstate);
end

% Extract the probability value from the ProbabilityResult object
Prob = jresult.probability;

end
