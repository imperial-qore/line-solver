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


if GlobalConstants.DummyMode
    Pn = NaN;
    return
end

if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    p = CPPLINE.probAggr(self.name, self.model, self.options);
    if isempty(p.ProbSys)
        CPPLINE.cppUnsupported(self.name, 'getProbSys', ...
            'line-cli''s -a prob answer for this solver carries no ProbSys value');
    end
    Pn = p.ProbSys;
    return
end

T0 = tic;
sn = self.getStruct;
% now compute marginal probability
options = self.getOptions;
Solver.resetRandomGeneratorSeed(options.seed);
% solver_nc_joint returns [Pr,G,lG,runtime]: lG is the THIRD output. Taking
% the second stored G into a field every other caller reads as a log, which
% corrupted any later getProb on the same solver object.
[Pn,~,lG] = solver_nc_joint(sn, self.options);
self.result.('solver') = getName(self);
self.result.Prob.logNormConstAggr = lG;
self.result.Prob.joint = Pn;
runtime = toc(T0);
self.result.runtime = runtime;
end