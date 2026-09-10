function varargout = getProbAggr(self,varargin)
% PROBAGGR = GETPROBAGGR(NODE, STATE)
% Aggregated state probability at a node, as the residence-time fraction the
% exact joint-state histogram of the LDES run assigns to STATE. Fully
% JSON-mediated (--export-histogram). Mirrors the Python-native getProbAggr().
%
% STATE is a per-class job-count vector, already aggregated over service phases;
% this is exactly the resolution the engine's histogram carries, so unlike
% getProb no conversion is applied. It is therefore NOT a synonym for getProb:
% delegating to it would put a detailed state through State.toMarginal a second
% time. If STATE is omitted the model's current state is aggregated and used.
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

function ProbAggr = getProbAggr_impl(self, node, state)
% GETPROBAGGR_IMPL Implementation of GETPROBAGGR; see the wrapper above.

if nargin < 3
    state = [];
end
if GlobalConstants.DummyMode
    ProbAggr = NaN;
    return
end

sn = self.model.getStruct;
if isa(node, 'Node')
    nodeIdx = node.index;
else
    nodeIdx = node;
end

ist = sn.nodeToStation(nodeIdx);
if isnan(ist) || ist < 1
    line_error(mfilename, sprintf('SolverLDES.getProbAggr: node %s is not a station; the LDES state histogram records station queue lengths only.', sn.nodenames{nodeIdx}));
end

if isempty(state)
    isf = sn.nodeToStateful(nodeIdx);
    if ~iscell(sn.state) || numel(sn.state) < isf || isempty(sn.state{isf})
        line_error(mfilename, 'SolverLDES.getProbAggr: no state was given and the model carries none for this node.');
    end
    [~, state] = State.toMarginal(sn, nodeIdx, sn.state{isf});
end

[space, time] = ldesHistogram(self);
ProbAggr = 0;
for row = 1:size(state, 1)
    ProbAggr = ProbAggr + ldesHistProb(space, time, sn.nclasses, ist, {state(row, :)});
end
end
