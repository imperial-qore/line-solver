function [lambdaFun, isTimeVarying, period] = sn_arrival_rate_fun(sn, ist, r)
% [LAMBDAFUN, ISTIMEVARYING, PERIOD] = SN_ARRIVAL_RATE_FUN(SN, IST, R)
%
% The arrival rate of station IST, class R AS A FUNCTION OF TIME.
%
% The time-varying analyses (Mt/G/Inf, the modified offered load, the
% Gt/Mt/st+GI fluid queue) consume lambda(t) itself, not a mean rate: their
% whole content is the LAG between when work arrives and when it is felt, and a
% time-averaged rate has no lag. LINE carries a time-varying arrival as a MAPt
% or an NHPP, whose sn.proc slot is a piecewise-constant schedule, so lambda(t)
% is read off the segment in force at t.
%
% For any other process the rate is constant and the handle returns it, which
% is what lets a caller ask for the time-varying analysis of a stationary model
% and get the stationary answer rather than an error.
%
% Parameters:
%   sn  - NetworkStruct
%   ist - station index
%   r   - class index
%
% Returns:
%   lambdaFun     - handle lambda(t), vectorized over t
%   isTimeVarying - whether the rate actually depends on t
%   period        - the cycle length when the schedule is cyclic, Inf otherwise
%
% See also SN_SCHEDULE_NOMINAL.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

isTimeVarying = false;
period = Inf;
isNHPP = sn.procid(ist,r) == ProcessType.NHPP;
isSched = isNHPP || sn.procid(ist,r) == ProcessType.MAPT || sn.procid(ist,r) == ProcessType.PHT;
if ~isSched
    lam = sn.rates(ist,r);
    lambdaFun = @(t) lam*ones(size(t));
    return
end
if isNHPP
    % An NHPP slot is {breakpoints, rates, cyclic}: the rates ARE lambda(t),
    % one per interval, so there is no MAP pair to reduce. The layout differs
    % from the MAPt one, so it is read here rather than through
    % SN_SCHEDULE_NOMINAL.
    slot = sn.proc{ist}{r};
    bp = slot{1}(:).';
    segRate = slot{2}(:).';
    cyclic = logical(slot{3});
else
    [~, ~, bp, segD0, segD1, cyclic] = sn_schedule_nominal(sn, ist, r);
    n = numel(segD0);
    % The arrival rate of segment k is pie_k * D1_k * e, the stationary
    % throughput of that segment's MAP.
    segRate = zeros(1, n);
    for k = 1:n
        m = size(segD0{k}, 1);
        if m == 1
            segRate(k) = segD1{k}(1,1);
        else
            pie_k = map_pie({segD0{k}, segD1{k}});
            segRate(k) = pie_k * segD1{k} * ones(m,1);
        end
    end
end
isTimeVarying = any(abs(segRate - segRate(1)) > GlobalConstants.Zero);
if cyclic
    period = bp(end) - bp(1);
end
lambdaFun = @(t) local_piecewise(t, bp, segRate, cyclic);
end

function v = local_piecewise(t, bp, segRate, cyclic)
% Segment k is in force on [bp(k), bp(k+1)). Before the first breakpoint the
% first segment holds and after the last the last one does, so a caller
% integrating over an infinite past (the Mt/G/Inf convolution) gets a defined
% rate everywhere rather than a NaN.
t = t(:).';
v = zeros(size(t));
T0 = bp(1); T1 = bp(end);
u = t;
if cyclic && T1 > T0
    u = T0 + mod(t - T0, T1 - T0);
end
for i = 1:numel(u)
    k = find(u(i) >= bp(1:end-1), 1, 'last');
    if isempty(k)
        k = 1;
    end
    k = min(k, numel(segRate));
    v(i) = segRate(k);
end
v = reshape(v, size(t));
end
