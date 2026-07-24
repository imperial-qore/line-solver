function svc = mam_svc_mixture(D_arr, pie_cell, D0_cell)
% SVC = MAM_SVC_MIXTURE(D_ARR, PIE_CELL, D0_CELL)
%
% Builds the arrival-weighted phase-type mixture of the per-class service
% laws at a station, as the service descriptor accepted by QSYS_MAPG1K and
% QSYS_MMAPG1K.
%
%   D_ARR    - {D0, D_class1, ..., D_classR} arrival MMAP
%   PIE_CELL - per-class PH initial distributions
%   D0_CELL  - per-class PH subgenerators
%
% The mixture is PH(alpha_mix, T_mix) with alpha_mix = [w_1*pie_1, ...],
% T_mix = blkdiag(D0_1, ...), and w_k = lambda_k/sum_j lambda_j the fraction
% of arrivals belonging to class k. It is therefore the service law of an
% arbitrary packet, and it reduces to the common law exactly (as a
% distribution) when every class shares one.
%
% This is the same construction MAM_TRUNCATE_RENORM already applies when a
% station carries more than one class, and it is factored out here so that
% the exact finite-buffer branch and the truncate-and-renormalize fallback
% rest on identical service assumptions and remain comparable.
%
% See also QSYS_MMAPG1K, MAM_TRUNCATE_RENORM, MAM_DETECT_MMCK.

nClasses = numel(pie_cell);
D0 = D_arr{1};
e_arr = ones(size(D0, 1), 1);
Dsum = zeros(size(D0));
for k = 1:nClasses
    Dsum = Dsum + D_arr{k+1};
end
theta = ctmc_solve(D0 + Dsum);
lambda_k = zeros(1, nClasses);
for k = 1:nClasses
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
for k = 1:nClasses
    alpha_mix(offset+1 : offset+n_k(k)) = w(k) * pie_cell{k}(:)';
    T_mix(offset+1:offset+n_k(k), offset+1:offset+n_k(k)) = D0_cell{k};
    offset = offset + n_k(k);
end

svc = struct('type', 'ph', 'alpha', alpha_mix, 'T', T_mix);
end
