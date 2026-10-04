function varargout = getProb(self,varargin)
% PNIR = GETPROB(NODE, STATE) Returns state probabilities at equilibrium
%
% @brief Returns state probabilities for a given node, including phase information
%
% This method computes the steady-state marginal state probability distribution
% for jobs at a specified station. Unlike getProb_aggr(), this method returns
% full state-space probabilities including phase information from service distributions.
%
% @param self SolverCTMC instance
% @param node Queue or node object (or node index as integer)
% @param state (optional) State specification as vector. If provided, returns
%             probability of this specific state. If omitted, returns probabilities
%             for all states at the node.
%
% @return Pnir Matrix of state probabilities. If state parameter was provided,
%         returns a scalar probability value. Otherwise returns a matrix where
%         rows represent different states and columns represent job classes.
%         For phase-type distributions, includes phase dimension in state space.
%
% @note Only supported by CTMC and SSA solvers due to their state-space formulation.
%       For aggregated probabilities (without phase info), use getProb_aggr() instead.
%
% @see getProb_aggr - Returns probabilities aggregated over phases
% @see getProbSys - Returns joint system state probabilities
%
% Example:
% @code
% solver = SolverCTMC(model);
% queue1 = model.nodes{1};
%
% % Get all state probabilities for queue1
% prob_matrix = solver.getProb(queue1);
%
% % Get probability for specific state [2 jobs of class 0, 1 job of class 1]
% queue1.setState([2, 1]);
% prob_specific = solver.getProb(queue1);
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

function Pnir = getProb_impl(self, node, state)
% GETPROB_IMPL Implementation of GETPROB; see the wrapper above.


% lang='cpp' answers this from the Prob column of -a prob, the DETAILED marginal
% (phases and buffer arrangement included) of the state the model carries; a
% state prior over several rows is refused by name. See CPPLINE.assertSingleState.
%
% A STATE NAMED IN THE CALL stays with the native path: `-a prob` reports the
% marginal of the state model.json carries and takes no per-call state, so
% serving one here would mean writing it onto the caller's model first.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp') && ~self.isChainSolver()
    if nargin >= 3
        CPPLINE.cppUnsupported(self.name, 'getProb(node, state)', ...
            ['line-cli reports the marginal at the state the model carries and takes no ' ...
            'per-call state, so answering a named one would mean writing it onto the ' ...
            'caller''s model first']);
    end
    CPPLINE.assertSingleState(self.name, 'getProb', self.model);
    if ~isnumeric(node), istc = node.index; else, istc = node; end
    Pnir = CPPLINE.probEntry(CPPLINE.probAggr(self.name, self.model, self.options), ...
        'Prob', istc, self.name, 'getProb');
    return
end

if self.isChainSolver()
    % Chain mode: getProb(state) returns the stationary probability of a
    % single state, identified by its row in the chain state space or, when
    % the chain carries none, by its 1-based state index.
    if nargin<2
        line_error(mfilename,'getProb in chain mode requires a state index or a state vector.');
    end
    pi = self.getProbSys();
    if isempty(self.chainModel.stateSpace)
        idx = node;
        if ~isscalar(idx) || idx~=round(idx) || idx<1 || idx>length(pi)
            line_error(mfilename,'The chain carries no state space, so getProb requires a state index in 1..%d.', length(pi));
        end
    else
        idx = matchrow(self.chainModel.stateSpace, reshape(node,1,[]));
        if idx<1
            line_error(mfilename,'The requested state is not in the chain state space.');
        end
    end
    Pnir = pi(idx);
    return
end

self.assertPhaseTypeStates('getProb');

if GlobalConstants.DummyMode
    Pnir = NaN;
    return
end

if nargin<2 %~exist('node','var')
    line_error(mfilename,'getProb requires to pass a parameter the station of interest.');
end
if ~isfield(self.options,'keep')
    self.options.keep = false;
end
T0 = tic;
sn = self.getStruct;
sn.state = sn.state;
if nargin>=3 %exist('state','var')
    sn.state{node} = state;
end
ind = node.index;
for isf=1:length(sn.state)
    isf_param = sn.nodeToStateful(ind);
    if isf ~= isf_param
        sn.state{isf} = sn.state{isf}*0 -1;
    end
end
Pnir = solver_ctmc_marg(sn, self.options);
self.result.('solver') = getName(self);
self.result.Prob.marginal = Pnir;
runtime = toc(T0);
self.result.runtime = runtime;
Pnir = Pnir(node);
end