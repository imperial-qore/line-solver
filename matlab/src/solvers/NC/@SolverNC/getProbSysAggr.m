function varargout = getProbSysAggr(self,varargin)
% PN = GETPROBSYSSTATEAGGR()
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

function Pn = getProbSysAggr_impl(self)
% GETPROBSYSAGGR_IMPL Implementation of GETPROBSYSAGGR; see the wrapper above.


% lang='cpp' answers this from -a prob; a state prior over several rows is
% refused there by name. See CPPLINE.assertSingleState.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    CPPLINE.assertSingleState(self.name, 'getProbSysAggr', self.model);
    Pn = CPPLINE.probAggr(self.name, self.model, self.options).ProbSysAggr;
    return
end

if GlobalConstants.DummyMode
    Pn = NaN;
    return
end

T0 = tic;
sn = self.getStruct;
% now compute marginal probability
options = self.getOptions;
Solver.resetRandomGeneratorSeed(options.seed);
% solver_nc_jointaggr returns [Pr,G,lG,runtime]: lG is the THIRD output.
[Pn,~,lG] = solver_nc_jointaggr(sn, self.options);
self.result.('solver') = getName(self);
self.result.Prob.logNormConstAggr = lG;
self.result.Prob.joint = Pn;
runtime = toc(T0);
self.result.runtime = runtime;
end