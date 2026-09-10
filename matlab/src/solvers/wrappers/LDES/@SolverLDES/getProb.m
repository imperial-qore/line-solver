function varargout = getProb(self,varargin)
% PROB = GETPROB(NODE, STATE)
% Steady-state marginal state probability at a node, as the residence-time
% fraction the exact joint-state histogram of the LDES run assigns to STATE.
% Fully JSON-mediated (--export-histogram). Mirrors the Python-native getProb().
%
% STATE is a detailed node state in the convention State.toMarginal reads; it is
% aggregated here to per-class job counts, which is the resolution the engine's
% histogram carries (it records integer queue lengths, not service phases). If
% STATE is omitted the model's current state for the node is used. Returns 0 if
% the state is never visited during the simulation.
%
% This is NOT computed from sample(): the transient QNt series holds interval
% time-averages of the queue length, so comparing it against an integer state
% matched almost nowhere and reported a near-zero probability for a state the
% chain spends most of its time in.
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

if nargin < 3
    state = [];
end
if GlobalConstants.DummyMode
    Prob = NaN;
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
    line_error(mfilename, sprintf('SolverLDES.getProb: node %s is not a station; the LDES state histogram records station queue lengths only.', sn.nodenames{nodeIdx}));
end

if isempty(state)
    isf = sn.nodeToStateful(nodeIdx);
    if iscell(sn.state) && numel(sn.state) >= isf && ~isempty(sn.state{isf})
        state = sn.state{isf};
    else
        line_error(mfilename, 'SolverLDES.getProb: no state was given and the model carries none for this node.');
    end
end

[~, nir] = State.toMarginal(sn, nodeIdx, state);

% A multi-row state names a SET of states, and its probability is the sum over
% the set, as it is in the CTMC and SSA getProb. The run is made ONCE and every
% row weighed against the same histogram; re-solving per row would draw a fresh
% sample path for each and the sum would not be a probability of anything.
[space, time] = ldesHistogram(self);
Prob = 0;
for row = 1:size(nir, 1)
    Prob = Prob + ldesHistProb(space, time, sn.nclasses, ist, {nir(row, :)});
end
end
