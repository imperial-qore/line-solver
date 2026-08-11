function Prob = getProb(self, node, state)
% PROBSTATE = GETPROBSTATE(NODE, STATE) Returns state probabilities estimated via simulation
%
% @brief Estimates steady-state marginal state probabilities at a node via SSA simulation
%
% This method estimates state probabilities by running a stochastic simulation and
% collecting statistics on how often each state is visited. The SSA solver performs
% CTMC stochastic sampling, recording state transitions and dwell times.
%
% Probabilities are estimated as the fraction of total time spent in each state
% during the simulation. This is an empirical estimate and converges to steady-state
% as the number of samples increases.
%
% @param self SolverSSA instance
% @param node Queue or node object (or node index as integer)
% @param state (optional) State specification as vector. If provided, returns
%             probability of this specific state. If omitted, returns probabilities
%             for all reachable states at the node with phase information.
%
% @return Prob Estimated state probability value
%         - If state parameter specified: scalar probability of that state
%         - If state not specified: probability mass across all states (sum=1)
%         - 0 if state not seen during simulation (rare event)
%         - Includes phase information from service distributions
%
% @note SSA estimates probabilities from finite simulation. Accuracy improves with
%       more samples. If the specified state is not observed during simulation,
%       a warning is issued and probability 0 is returned.
%       For phase-detailed states (CTMC-style), use getProb().
%
% @warning Simulation-based estimate - results vary between runs unless seed is set.
%          Use 'seed' option in constructor for reproducible results.
%          For small models, CTMC provides exact probabilities.
%
% @see getProbAggr - Returns state probabilities aggregated over phases
% @see getProbSys - Returns joint system state probabilities
% @see getAvg - Get average metrics (alternative analysis method)
%
% Example:
% @code
% solver = SolverSSA(model, 'samples', 50000, 'seed', 42);
% queue1 = model.nodes{1};
%
% % Estimate probability of state [2 jobs of class 0, 1 of class 1]
% queue1.setState([2, 1]);
% prob_state = solver.getProb(queue1);
% fprintf('Estimated Pr[state=[2,1]] = %.4f\\n', prob_state);
%
% % Get all state probabilities
% all_probs = solver.getProb(queue1);
% @endcode

if GlobalConstants.DummyMode
    Prob = NaN;
    return
end

switch self.options.method
    case {'default','nrm'}
        self.options.method = 'serial';
end
% we do not use probSysState as that is for joint states
[~, tranSysState] = self.runAnalyzer;
sn = self.getStruct;
isf = sn.nodeToStateful(node.index);
TSS = cell2mat({tranSysState{1},tranSysState{1+isf}});
TSS(:,1)=[TSS(1,1);diff(TSS(:,1))];
if nargin<3 %~exist('state','var')
    if size(sn.state{isf},1)>1
        error('There are multiple station states, choose an initial state as a parameter to getProb.');
    end
    state = sn.state{isf};
end
% add padding of zeros for FCFS stations
state = [zeros(1,size(TSS(:,2:end),2)-size(state,2)),state];
rows = findrows(TSS(:,2:end), state);
if ~isempty(rows)
    Prob = sum(TSS(rows,1))/sum(TSS(:,1));
else
    line_warning(mfilename,'The state was not seen during the simulation.\n');
    Prob = 0;
end
end