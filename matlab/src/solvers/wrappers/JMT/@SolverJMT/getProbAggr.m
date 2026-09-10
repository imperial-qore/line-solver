function varargout = getProbAggr(self,varargin)
% PR = GETPROBAGGR(NODE, STATE_A)
%
% Probability of a SPECIFIC per-class job distribution at a station.
% Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) for given state.
%
% Compare with getProbMarg: returns total queue-length distribution,
% i.e., P(n total jobs) summed over all class combinations.
%
% Input:
%   node    - Node object
%   state_a - Per-class job counts, e.g., [2,1] = 2 class-1, 1 class-2
%
% Output:
%   Pr - Scalar probability in [0,1] (estimated via simulation)
% The result recorder captures the scalar this getter returned together
% with the solver that produced it -- see LineResultRecorder. Six of the
% statepr_* goldens hold exactly this number and nothing else, so recording
% it is what makes those goldens attributable instead of "the first bare
% number the example printed". The wrapper exists so that recording happens
% on EVERY exit path of the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getProbAggr_impl(self,varargin{:});
LineResultRecorder.captureScalar(scope, self, 'probAggr', varargout{1});
end

function Pr = getProbAggr_impl(self, node, state_a)
% GETPROBAGGR_IMPL Implementation of GETPROBAGGR; see the wrapper above.
if GlobalConstants.DummyMode
    Pr = NaN;
    return
end
% lang='cpp' answers this from `-s jmt -a prob`, which weighs ONE instrumented
% run: line-cli logs every station, resamples the traces onto a common grid and
% reports the time each station's declared per-class counts are held for. That is
% the same estimator this file builds out of sampleAggr, and it is the C++ engine
% that ran the simulation, so the number is line-cli's and not MATLAB's. A state
% prior over several rows is refused there by name (CPPLINE.assertSingleState).
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    CPPLINE.assertSingleState(self.name, 'getProbAggr', self.model);
    sn = self.getStruct;
    ist = sn.nodeToStation(node.index);
    flags = {};
    if nargin >= 3 && ~isempty(state_a)
        % `--node` names the station whose declared counts `--state` replaces;
        % every other station keeps its own, as the reference's substitution
        % into sn.state{isf} does.
        flags = [{'--node', num2str(node.index)}, {'--state', ...
            strjoin(arrayfun(@(v) num2str(v), state_a(:)', 'UniformOutput', false), ',')}];
    end
    payload = CPPLINE.analysisViaCpp(self.name, self.model, self.options, 'prob', flags);
    p = struct('ProbSys', [], 'ProbSysAggr', [], 'Prob', [], 'ProbAggr', []);
    if isfield(payload, 'ProbAggr')
        p.ProbAggr = CPPLINE.jsonNumericList(payload.ProbAggr);
    end
    Pr = CPPLINE.probEntry(p, 'ProbAggr', ist, self.name, 'getProbAggr');
    return
end

sn = self.getStruct;
if nargin<3 %~exist('state_a','var')
    state_a = sn.state{sn.nodeToStation(node.index)};
end
stationStateAggr = self.sampleAggr(node);
rows = findrows(stationStateAggr.state, state_a);
t = stationStateAggr.t;
dt = [diff(t);0];
Pr = sum(dt(rows))/sum(dt);
end