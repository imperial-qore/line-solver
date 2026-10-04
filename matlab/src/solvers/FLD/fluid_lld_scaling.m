function [a, da] = fluid_lld_scaling(lldrow, n)
% [A, DA] = FLUID_LLD_SCALING(LLDROW, N)
%
% Limited load-dependent rate scaling at a CONTINUOUS population.
%
% `sn.lldscaling(i,:)` tabulates the station rate multiplier at integer
% populations 1..lldlimit, and the discrete solvers read it as
% `lldscaling(i, min(n, lldlimit))`. The fluid state is continuous, so the
% table is read here by linear interpolation between consecutive entries,
% clamped to the first entry below n=1 and to the last entry above the table
% end, matching the clamping the CTMC already applies.
%
% Parameters:
%   lldrow - (1 x lldlimit) rate multipliers at populations 1..lldlimit;
%            empty means no load dependence and returns a = 1, da = 0
%   n      - population, scalar or vector, may be non-integer or negative
%
% Returns:
%   a  - interpolated scaling alpha(n)
%   da - its derivative d alpha / d n, zero on the clamped tails
%
% See also FLUID_CAPACITY_CLOSURE, ODE_RATES_CLOSING_FACTORS.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = n(:);
if isempty(lldrow)
    a = ones(size(n));
    da = zeros(size(n));
    return
end

lldrow = lldrow(:)';
L = numel(lldrow);
a = zeros(size(n));
da = zeros(size(n));

lo = n <= 1;
a(lo) = lldrow(1);

hi = n >= L;
a(hi) = lldrow(L);

mid = ~lo & ~hi;
if any(mid)
    k = floor(n(mid));
    frac = n(mid) - k;
    a(mid) = lldrow(k)'.*(1-frac) + lldrow(k+1)'.*frac;
    da(mid) = lldrow(k+1)' - lldrow(k)';
end
end
