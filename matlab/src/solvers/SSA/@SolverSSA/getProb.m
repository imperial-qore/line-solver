function varargout = getProb(self,varargin)
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
% The result recorder captures the scalar this getter returned together
% with the solver that produced it -- see LineResultRecorder. Six of the
% statepr_* goldens hold exactly this number and nothing else, so recording
% it is what makes those goldens attributable instead of "the first bare
% number the example printed". The wrapper exists so that recording happens
% on EVERY exit path of the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getProb_impl(self,varargin{:});
LineResultRecorder.captureScalar(scope, self, 'prob', varargout{1});
end

function Prob = getProb_impl(self, node, state)
% GETPROB_IMPL Implementation of GETPROB; see the wrapper above.


% lang='cpp' answers this from the Prob column of -a prob, the DETAILED marginal
% of the state the model carries; a state named in the CALL stays with the native
% path. See CPPLINE.assertSingleState.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    if nargin >= 3
        CPPLINE.cppUnsupported(self.name, 'getProb(node, state)', ...
            ['line-cli reports the occupancy of the state the model carries and takes no ' ...
            'per-call state, so answering a named one would mean writing it onto the ' ...
            'caller''s model first']);
    end
    CPPLINE.assertSingleState(self.name, 'getProb', self.model);
    if ~isnumeric(node), istc = node.index; else, istc = node; end
    snc = self.getStruct;
    istc = snc.nodeToStation(istc);  % the native reads nodeToStateful(node.index)
    Prob = CPPLINE.probEntry(CPPLINE.probAggr(self.name, self.model, self.options), ...
        'Prob', istc, self.name, 'getProb');
    return
end

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