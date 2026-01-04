function Pnir = getProbAggr(self, node, state_a)
% PNIR = GETPROBAGGR(NODE, STATE_A)
%
% Probability of a SPECIFIC per-class job distribution at a station.
% Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) at the node.
%
% Compare with getProbMarg: returns total queue-length distribution,
% i.e., P(n total jobs) summed over all class combinations.
%
% Input:
%   node    - Queue or node object
%   state_a - Per-class job counts, e.g., [2,1] = 2 class-1, 1 class-2
%
% Output:
%   Pnir - Scalar probability in [0,1]

if GlobalConstants.DummyMode
    Pnir = NaN;
    return
end

T0 = tic;
sn = self.getStruct;

% Get unnormalized probability using original NC method
if nargin<3
    state_a = sn.state{sn.nodeToStateful(node.index)};
end

% Store original state and set requested state
ist = sn.nodeToStation(node.index);
original_state = sn.state{ist};
sn.state{ist} = state_a;

options = self.getOptions;
Solver.resetRandomGeneratorSeed(options.seed);

self.result.('solver') = getName(self);
if isfield(self.result,'Prob') && isfield(self.result.Prob,'logNormConstAggr') && isfinite(self.result.Prob.logNormConstAggr)
    [Pnir_vec,lG] = solver_nc_margaggr(sn, self.options, self.result.Prob.logNormConstAggr);
else
    [Pnir_vec,lG] = solver_nc_margaggr(sn, self.options);
    self.result.Prob.logNormConstAggr = lG;
end
self.result.Prob.marginal = Pnir_vec;

% solver_nc_margaggr already returns normalized probabilities
% (it computes P = F_i * G(-i) / G which is properly normalized)
% So we simply extract the probability for this station
Pnir = Pnir_vec(ist);

% Restore original state
sn.state{ist} = original_state;

runtime = toc(T0);
self.result.runtime = runtime;

end
