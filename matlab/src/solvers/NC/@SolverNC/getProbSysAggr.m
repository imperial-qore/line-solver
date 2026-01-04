function Pn = getProbSysAggr(self)
% PN = GETPROBSYSSTATEAGGR()

if GlobalConstants.DummyMode
    Pn = NaN;
    return
end

T0 = tic;
sn = self.getStruct;
% now compute marginal probability
options = self.getOptions;
Solver.resetRandomGeneratorSeed(options.seed);
[Pn,lG] = solver_nc_jointaggr(sn, self.options);
self.result.('solver') = getName(self);
self.result.Prob.logNormConstAggr = lG;
self.result.Prob.joint = Pn;
runtime = toc(T0);
self.result.runtime = runtime;
end