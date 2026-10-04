function [D0bar, D1bar, breakpoints, segD0, segD1, cyclic, segMark, segBatch] = sn_schedule_nominal(sn, ist, r)
% [D0BAR, D1BAR, BREAKPOINTS, SEGD0, SEGD1, CYCLIC, SEGMARK, SEGBATCH] = SN_SCHEDULE_NOMINAL(SN, IST, R)
%
% Unpacks the MAPt, PHt, MMAPt, MPHt or BMMAPt slot of sn.proc at station IST, class R
% and returns both the width-weighted time-averaged (D0, D1) pair and the
% per-segment pairs.
%
% The slot is {breakpoints, segments1, segments2, cyclic}: for a MAPt the two
% cells are the D0 and D1 matrices; for a PHt they are the alpha row vectors and
% the sub-generators S, whose equivalent MAP pair is (S, s*alpha) with
% s = -S*e. The nominal is the stationary carrier of the phase structure that the
% schedule modulates, so its order is the phase count and it is what the fluid
% base rates are built from.
%
% A MARKED schedule (MMAPt, and MPHt which is stored lowered to MMAPt form)
% APPENDS a fifth entry holding the per-mark blocks, a 1-by-K cell of 1-by-n
% cells. Its first four entries are exactly a MAPt slot carrying the AGGREGATE
% sum_c D1k{c}, so every caller that wants only the unmarked schedule reads it
% unchanged and needs no knowledge of the marks. SEGMARK is that fifth entry,
% and is empty for an unmarked schedule.
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
% A PHt is the only family whose two cells are NOT already a (D0, D1) pair; the
% marked schedules store theirs lowered, so they take the MAPt arm.
isMAPt = sn.procid(ist,r) ~= ProcessType.PHT;
segMark = {};
if numel(slot) >= 5
    segMark = slot{5};
end
segBatch = {};
if numel(slot) >= 6
    segBatch = slot{6};
end
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
