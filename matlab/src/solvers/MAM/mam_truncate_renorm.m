function [meanQ, lossProb, p_norm] = mam_truncate_renorm(D_arr, pie_cell, D0_cell, capK)
% MAM_TRUNCATE_RENORM Finite-buffer marginal for MMAP[K]/PH[K]/1/FCFS.
%
% Solves the infinite-buffer MMAP[K]/PH[K]/1/FCFS queue via BUTools
% MMAPPH1FCFS, truncates the marginal queue length distribution at the
% buffer capacity capK, and renormalizes. For an M/M/1 input the
% renormalized distribution coincides exactly with the M/M/1/K marginal;
% for general MMAP/PH it is an ASTA-style approximation.
%
% INPUT
%   D_arr    - cell array passed as the MMAP argument to MMAPPH1FCFS,
%              e.g. {D0, D_class1, D_class2, ...}
%   pie_cell - cell array of per-class PH initial distributions
%   D0_cell  - cell array of per-class PH subgenerators (T matrices)
%   capK     - buffer capacity (max number of jobs in system, integer >= 1)
%
% OUTPUT
%   meanQ    - mean number of jobs in system (clipped to [0, capK])
%   lossProb - blocking probability = renormalized boundary mass p(N=capK)
%              (PASTA-exact for Poisson arrivals)
%   p_norm   - renormalized truncated marginal as 1x(capK+1) row vector
%
% Multi-class: MMAPPH1FCFS returns the joint (aggregate) marginal at the
% station; callers split per-class with the lambda_k/lambda_total fraction
% (FCFS Little's-law decomposition), mirroring the convention used by the
% closed-class branch of solver_mam_basic.
%
% Multi-server: callers must pre-scale the PH service mean by 1/c (the
% existing "surrogate-delay" trick); the helper then returns the marginal
% of the scaled M/M/1/K system, which is treated as an approximation of
% the M/M/c/K marginal. The c-1/c surrogate-delay term is omitted under
% finite cap because it would double-count physical jobs already bounded
% by capK.
%
% See also MMAPPH1FCFS.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Multi-class needs aggregation: MMAPPH1FCFS('ncDistr') returns the per-class
% marginal P(N_k=n), not the joint P(N_total=n) that the truncation requires.
% Aggregate to a single-class MMAP/PH/1 by summing arrival matrices and
% building a workload-weighted PH mixture for service.
nClasses = numel(pie_cell);
if nClasses > 1
    D0 = D_arr{1};
    Dsum = zeros(size(D0));
    for k=1:nClasses
        Dsum = Dsum + D_arr{k+1};
    end
    e_arr = ones(size(D0,1), 1);
    theta = ctmc_solve(D0 + Dsum);
    lambda_k = zeros(1, nClasses);
    for k=1:nClasses
        lambda_k(k) = theta * D_arr{k+1} * e_arr;
    end
    sumL = sum(lambda_k);
    if sumL > 0
        w = lambda_k / sumL;
    else
        w = ones(1, nClasses) / nClasses;
    end
    n_k = cellfun(@(p) numel(p), pie_cell);
    n_total = sum(n_k);
    alpha_mix = zeros(1, n_total);
    T_mix = zeros(n_total, n_total);
    offset = 0;
    for k=1:nClasses
        alpha_mix(offset+1 : offset+n_k(k)) = w(k) * pie_cell{k}(:)';
        T_mix(offset+1:offset+n_k(k), offset+1:offset+n_k(k)) = D0_cell{k};
        offset = offset + n_k(k);
    end
    D_call = {D0, Dsum};
    pie_call = {alpha_mix};
    D0_call = {T_mix};
else
    D_call = D_arr;
    pie_call = pie_cell;
    D0_call = D0_cell;
end

pdistr = MMAPPH1FCFS(D_call, pie_call, D0_call, 'ncDistr', capK + 1);
pdistr = abs(pdistr(:)');
nLevels = capK + 1;
if numel(pdistr) < nLevels
    pdistr(end+1:nLevels) = 0;
end
p_in = pdistr(1:nLevels);
massIn = sum(p_in);
if massIn <= 0
    p_norm = zeros(1, nLevels);
    p_norm(1) = 1;
else
    p_norm = p_in / massIn;
end
meanQ = max(0, min(capK, (0:capK) * p_norm(:)));
lossProb = p_norm(end);
end
