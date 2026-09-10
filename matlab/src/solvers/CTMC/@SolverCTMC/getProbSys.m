function varargout = getProbSys(self,varargin)
% PN = GETPROBSYSSTATE()
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

function Pn = getProbSys_impl(self)
% GETPROBSYS_IMPL Implementation of GETPROBSYS; see the wrapper above.


% lang='cpp' answers this from -a prob; a state prior over several rows is
% refused there by name. See CPPLINE.assertSingleState.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    CPPLINE.assertSingleState(self.name, 'getProbSys', self.model);
    Pn = CPPLINE.probAggr(self.name, self.model, self.options).ProbSys;
    return
end

if self.isChainSolver()
    % Chain mode: the stationary vector of the user-supplied chain.
    if isempty(self.result) || ~isfield(self.result,'Prob')
        self.runAnalyzer();
    end
    Pn = self.result.Prob.joint;
    return
end

self.assertPhaseTypeStates('getProbSys');

if GlobalConstants.DummyMode
    Pn = NaN;
    return
end

if ~isfield(self.options,'keep')
    self.options.keep = false;
end
T0 = tic;
sn = self.getStruct;
if self.model.isStateValid
    Pn = solver_ctmc_joint(sn, self.options);
    self.result.('solver') = getName(self);
    self.result.Prob.joint = Pn;
else
    line_error(mfilename,'The model state is invalid.');
end
runtime = toc(T0);
self.result.runtime = runtime;
end