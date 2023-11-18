function Pnir = getProbAggr(self, node, state_a)
% PNIR = GETPROBAGGR(NODE, STATE_A)

if GlobalConstants.DummyMode
    Pnir = NaN;
    return
end

T0 = tic;
sn = self.getStruct;
if nargin<3 %~exist('state_a','var')
    state_a = sn.state{sn.nodeToStateful(node.index)};
end
% now compute marginal probability
ist = sn.nodeToStation(node.index);
sn.state{ist} = state_a;

options = self.getOptions;
Solver.resetRandomGeneratorSeed(options.seed);

self.result.('solver') = getName(self);
if isfield(self.result,'Prob') && isfield(self.result.Prob,'logNormConstAggr') && isfinite(self.result.Prob.logNormConstAggr)
    [Pnir,lG] = solver_nc_margaggr(sn, self.options, self.result.Prob.logNormConstAggr);
else
    [Pnir,lG] = solver_nc_margaggr(sn, self.options);
    self.result.Prob.logNormConstAggr = lG;
end
self.result.Prob.marginal = Pnir;
runtime = toc(T0);
self.result.runtime = runtime;
Pnir = Pnir(ist);
end