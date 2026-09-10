function beta = fes_beta_handle(scalingTable, cutoffs)
% BETA = FES_BETA_HANDLE(SCALINGTABLE, CUTOFFS)
%
% Wrap a flow-equivalent-server (FES) throughput table as a per-class
% class-dependence function beta_{i,r}(n).
%
% SCALINGTABLE is a cell {1 x K} in which SCALINGTABLE{r} is the linearized
% vector of class-r throughputs X_r(n) of the aggregated subnetwork, indexed by
% LJD_LINEARIZE(min(n,CUTOFFS), CUTOFFS). CUTOFFS is the per-class population
% vector the table was tabulated on.
%
% The returned handle takes the per-class population vector n at the station and
% returns the length-K vector of DIMENSIONLESS class-dependence scalings
%   beta_r(n) = X_r(n) * |n| / n_r,
% relative to the nominal rate-1 service of the FES station. The |n|/n_r factor
% cancels the processor-sharing share that the convolution applies (Sauer 1983,
% "Computational Algorithms for State-Dependent Queueing Networks", eq. (40),
% with mu_{r,i}(n) = (n_r/|n|) beta_r(n)), leaving the aggregate completing class
% r at exactly the subnetwork throughput X_r(n). The population is clamped to
% CUTOFFS, so the scaling saturates beyond the tabulated range as the underlying
% table intends. Entries with n_r = 0 are never consulted by the recurrence.
%
% This is the single class-dependence mechanism used across the solvers: the
% exact convolution (PFQN_CONV) reads mu_{r,i}(n) from it, and AMVA-QD reads the
% same handle through PFQN_CDFUN. A table is materialized only at a language
% boundary (see JLINE.handle_to_serializablefun), never in the model.
%
% See also FES_COMPUTE_THROUGHPUTS, PFQN_CDFUN, PFQN_CONV, LJD_LINEARIZE.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

scalingTable = scalingTable(:)';
cutoffs = round(cutoffs(:)');
beta = @(n) fes_beta_eval(n, scalingTable, cutoffs);
end

function v = fes_beta_eval(n, scalingTable, cutoffs)
K = numel(scalingTable);
v = ones(1, K);
n = round(n(:)');
if numel(n) < numel(cutoffs)
    n(end+1:numel(cutoffs)) = 0;
elseif numel(n) > numel(cutoffs)
    n = n(1:numel(cutoffs));
end
nClamped = max(0, min(n, cutoffs));
idx = ljd_linearize(nClamped, cutoffs);
tot = sum(n);
for r = 1:K
    tbl = scalingTable{r};
    if ~isempty(tbl) && idx >= 1 && idx <= numel(tbl) && n(r) > 0
        % see _kb/03-api-layer.md (fes_beta_handle) for rationale
        v(r) = tbl(idx) * tot / n(r);
    end
end
end
