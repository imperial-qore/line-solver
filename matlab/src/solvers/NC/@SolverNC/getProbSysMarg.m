function [Pn, lPn] = getProbSysMarg(self, nvec, engine)
% [PN, LPN] = GETPROBSYSMARG(NVEC)
% [PN, LPN] = GETPROBSYSMARG(NVEC, ENGINE)
%
% Joint probability that station I holds NVEC(I) jobs IN TOTAL, all classes
% summed out:
%
%   PN = P(n_1 = NVEC(1), ..., n_M = NVEC(M))
%
% Compare with getProbSysAggr, which fixes the PER-CLASS population of every
% station and is a product form; each value returned here is the sum of
% getProbSysAggr over every per-class table with these row sums. Compare also
% with getProbMarg, which is the one-station marginal of this law.
%
% The quantity is a matrix permanent of the demand matrix replicated once per
% job (Ryser 1963 for the evaluation), so it needs no enumeration of that
% fibre.
%
% Input:
%   NVEC   - (1 x nstations) per-station total job counts; must sum to the
%            total closed population
%   ENGINE - permanent engine, one of 'exact' (default), 'spm', 'bethe',
%            'heur', 'huberlaw', 'adapart'. Only 'exact' is exact; the others
%            are refused on a demand matrix with a structural zero rather than
%            having it floored, since they need full support. 'spm' is the
%            saddle-point expansion, the one whose error falls as the per-class
%            populations grow and whose cost does not grow with them; see
%            PFQN_JOINTMARG for the measured accuracy and when it does not
%            apply.
%
% Output:
%   PN  - joint probability
%   LPN - its logarithm, which survives populations PN underflows at

if GlobalConstants.DummyMode
    Pn = NaN;
    lPn = NaN;
    return
end

if nargin < 3
    engine = 'exact';
end

if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    [Pn, lPn] = CPPLINE.probSysMarg(self.name, self.model, self.options, nvec, engine);
    self.lastPermEngine = lower(char(engine));
    return
end

T0 = tic;
sn = self.getStruct;
options = self.getOptions;
Solver.resetRandomGeneratorSeed(options.seed);

% Reuse the constant when a previous getter already paid for it: sweeping the
% whole lattice of total states otherwise recomputes G once per state.
if isfield(self.result,'Prob') && isfield(self.result.Prob,'logNormConstAggr') && isfinite(self.result.Prob.logNormConstAggr)
    lGin = self.result.Prob.logNormConstAggr;
else
    lGin = [];
end

[Pn, lPn, lG] = solver_nc_jointmarg(sn, options, nvec, engine, lGin);

self.lastPermEngine = lower(char(engine));
self.result.('solver') = getName(self);
self.result.Prob.logNormConstAggr = lG;
self.result.Prob.jointMarg = Pn;
self.result.runtime = toc(T0);
end
