function mu = sn_compat_rate(compat, counts, rates, n)
% MU = SN_COMPAT_RATE(COMPAT, COUNTS, RATES, N)
%
% Total service rate of a station served by heterogeneous server pools with a
% class-compatibility graph. A pool t holds COUNTS(t) identical servers, each
% running at RATES(t), and may serve operand j when COMPAT(t,j) is nonzero. The
% rate the station clears in state N is
%
%   mu(n) = sum_t counts(t)*rates(t)*min(1, sum_{j: compat(t,j) ~= 0} n(j))
%
% the ACTIVATED-SERVER law: a pool contributes its full rate as soon as it is
% compatible with at least one operand PRESENT. This is the order-independent
% reading of a compatibility structure -- at an INTEGER state mu depends on n
% only through its SUPPORT, so it is invariant to the arrival order and to any
% permutation of the microstate, which is exactly the condition an OI station
% has to meet. It is also what PAS_COMPATIBILITY_5CLASS encodes for a flat
% Network, so the layered and flat readings of one compatibility matrix agree.
%
% WHY min(1,.) AND NOT AN INDICATOR. At every integer state the two agree
% exactly -- a pool with at least one compatible job present is fully active,
% one with none is idle -- so nothing about the OI law on the real state lattice
% changes. They part company only at a FRACTIONAL argument, which is what a
% mean-value solver hands this function: AMVA evaluates the rate at a MEAN
% population, and under a hard indicator any operand with a mean above zero,
% however small, activates every pool it touches. A compatibility structure
% would then be invisible to AMVA whenever every operand is a little bit busy --
% which is nearly always. Scaling linearly below one job keeps the structure
% visible at the evaluation point while leaving the integer-state law untouched;
% it is the ordinary continuous relaxation of a step function, and the CTMC and
% simulation paths, which only ever evaluate at integer states, cannot tell the
% difference.
%
% IT IS NOT A MATCHING. A pool of two servers compatible with a class holding
% ONE job contributes both servers here, which over-counts against a
% non-redundant system where one server serves one job. That is deliberate: the
% matching size depends on the counts and not only on the support, so it is NOT
% order independent and would take the station outside the product form the OI
% closure is built on. A model that means the matching wants a different
% station, not a different reading of this one.
%
% Use SN_COMPAT_PEAK for the rate with every pool active, which normalizes the
% utilization of a rate-scaled station as U = T*S/peak.
%
% Reference:
%   Dorsman, Gardner (2024). New directions in pass-and-swap queues. Queueing
%   Systems 107:205-256, Fig. 1.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

npools = numel(counts);
if numel(rates) ~= npools
    line_error(mfilename, 'one rate per pool is required');
end
if size(compat,1) ~= npools
    line_error(mfilename, 'compat must have one row per pool');
end
if numel(n) ~= size(compat,2)
    line_error(mfilename, 'n must have one entry per operand');
end

nrow = n(:)';
nrow(nrow < 0) = 0;
mu = 0;
for t = 1:npools
    % The pool is activated ONCE by the jobs it can reach, not once per operand:
    % its weight is the compatible load, capped at one job.
    load = min(1, sum(nrow(compat(t,:) ~= 0)));
    mu = mu + counts(t) * rates(t) * load;
end
end
