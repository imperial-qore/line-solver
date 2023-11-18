function Pn = getProbSys(self)
% PN = GETPROBSYSSTATE()

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