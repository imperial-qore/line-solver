function Pnir = getProb(self, node, state)
% PNIR = GETPROB(NODE, STATE)

if GlobalConstants.DummyMode
    Pnir = NaN;
    return
end

if nargin<3 %~exist('state','var')
    state = sn.state{sn.nodeToStateful(node.index)};
end
T0 = tic;
sn = self.getStruct;
% now compute marginal probability
if isa(node,'Node')
    ist = sn.nodeToStation(node.index);
else
    ist = node;    
end
sn.state{ist} = state;

options = self.getOptions;
Solver.resetRandomGeneratorSeed(options.seed);

if ~isempty(self.result) && ~isempty(self.result.Prob) && isfield(self.result.Prob,'logNormConstAggr') && isfinite(self.result.Prob.logNormConstAggr)
    [Pnir,lG] = solver_nc_marg(sn, self.options, self.result.Prob.logNormConstAggr);
else
    [Pnir,lG] = solver_nc_marg(sn, self.options);
    self.result.Prob.logNormConstAggr = lG;
end
self.result.('solver') = getName(self);
self.result.Prob.marginal = Pnir;
runtime = toc(T0);
self.result.runtime = runtime;
Pnir = Pnir(ist);
end