function [D0bar, D1bar, breakpoints, segD0, segD1, cyclic] = sn_schedule_nominal(sn, ist, r)
% [D0BAR, D1BAR, BREAKPOINTS, SEGD0, SEGD1, CYCLIC] = SN_SCHEDULE_NOMINAL(SN, IST, R)
%
% Unpacks the MAPt or PHt slot of sn.proc at station IST, class R and returns
% both the width-weighted time-averaged (D0, D1) pair and the per-segment pairs.
%
% The slot is {breakpoints, segments1, segments2, cyclic}: for a MAPt the two
% cells are the D0 and D1 matrices; for a PHt they are the alpha row vectors and
% the sub-generators S, whose equivalent MAP pair is (S, s*alpha) with
% s = -S*e. The nominal is the stationary carrier of the phase structure that the
% schedule modulates, so its order is the phase count and it is what the fluid
% base rates are built from.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

slot = sn.proc{ist}{r};
if numel(slot) < 4
    line_error(mfilename, 'sn_schedule_nominal: a MAPt/PHt slot of sn.proc must hold {breakpoints, A, B, cyclic}');
end
breakpoints = slot{1}(:).';
first = slot{2};
second = slot{3};
cyclic = logical(slot{4});
n = numel(first);
isMAPt = sn.procid(ist,r) == ProcessType.MAPT;
segD0 = cell(1, n);
segD1 = cell(1, n);
for k = 1:n
    if isMAPt
        segD0{k} = first{k};
        segD1{k} = second{k};
    else
        Sm = second{k};
        segD0{k} = Sm;
        segD1{k} = (-sum(Sm, 2)) * first{k}(:).';
    end
end
widths = diff(breakpoints);
total = sum(widths);
D0bar = zeros(size(segD0{1}));
D1bar = zeros(size(segD1{1}));
for k = 1:n
    D0bar = D0bar + (widths(k)/total) * segD0{k};
    D1bar = D1bar + (widths(k)/total) * segD1{k};
end
end
