function varargout = getProbSysAggr(self,varargin)
% PROBSYSAGGR = GETPROBSYSAGGR()
% Aggregated joint system-state probability. Fully JSON-mediated. Mirrors the
% Python-native getProbSysAggr().
%
% This equals getProbSys() because the engine's state histogram is aggregated
% already: it records integer queue lengths per station and class, with no phase
% resolution, so the detailed and the aggregate joint question have the same
% answer here. getProbSys aggregates the model's current state before matching,
% so no second aggregation is needed.
% The result recorder captures the scalar this getter returned together
% with the solver that produced it -- see LineResultRecorder. Six of the
% statepr_* goldens hold exactly this number and nothing else, so recording
% it is what makes those goldens attributable instead of "the first bare
% number the example printed". The wrapper exists so that recording happens
% on EVERY exit path of the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getProbSysAggr_impl(self,varargin{:});
LineResultRecorder.captureScalar(scope, self, 'probSysAggr', varargout{1});
end

function ProbSysAggr = getProbSysAggr_impl(self)
% GETPROBSYSAGGR_IMPL Implementation of GETPROBSYSAGGR; see the wrapper above.
ProbSysAggr = self.getProbSys();
end
