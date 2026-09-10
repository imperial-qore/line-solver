function varargout = getProbSys(self,varargin)
% PROBSYS = GETPROBSYS()
% Joint steady-state probability of the current system state, as the
% residence-time fraction the exact joint-state histogram of the LDES run
% assigns to it. Fully JSON-mediated (--export-histogram). Mirrors the
% Python-native getProbSys().
%
% Every STATION is constrained to the per-class job counts its current state
% aggregates to. A stateful node that is not a station -- a Router, a Cache, a
% Place -- holds no queue length in the engine's histogram and cannot be
% constrained, so it is excluded and the answer is the joint law of the station
% queue lengths alone.
%
% This is NOT computed from sampleSys(): the transient QNt series holds interval
% time-averages of the queue length, so comparing it against integer states
% matched almost nowhere and reported a near-zero probability for a state the
% chain spends most of its time in.
% The result recorder captures the scalar this getter returned together
% with the solver that produced it -- see LineResultRecorder. Six of the
% statepr_* goldens hold exactly this number and nothing else, so recording
% it is what makes those goldens attributable instead of "the first bare
% number the example printed". The wrapper exists so that recording happens
% on EVERY exit path of the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getProbSys_impl(self,varargin{:});
LineResultRecorder.captureScalar(scope, self, 'probSys', varargout{1});
end

function ProbSys = getProbSys_impl(self)
% GETPROBSYS_IMPL Implementation of GETPROBSYS; see the wrapper above.

if GlobalConstants.DummyMode
    ProbSys = NaN;
    return
end

sn = self.model.getStruct;
stations = [];
targets = {};
for ind = 1:sn.nnodes
    ist = sn.nodeToStation(ind);
    if isnan(ist) || ist < 1
        continue
    end
    isf = sn.nodeToStateful(ind);
    if ~iscell(sn.state) || numel(sn.state) < isf || isempty(sn.state{isf})
        line_error(mfilename, sprintf('SolverLDES.getProbSys: the model carries no current state for station %s.', sn.nodenames{ind}));
    end
    [~, nir] = State.toMarginal(sn, ind, sn.state{isf});
    stations(end+1) = ist; %#ok<AGROW>
    targets{end+1} = nir(1, :); %#ok<AGROW>
end

if isempty(stations)
    ProbSys = 0;
    return
end

[space, time] = ldesHistogram(self);
ProbSys = ldesHistProb(space, time, sn.nclasses, stations, targets);
end
