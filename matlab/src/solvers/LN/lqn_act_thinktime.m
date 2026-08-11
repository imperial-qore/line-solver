function zt = lqn_act_thinktime(lqn, aidx)
% ZT = LQN_ACT_THINKTIME(LQN, AIDX)
%
% Think time of activity AIDX, zero when it has none.
%
% An activity think time is a delay in series with that activity's host demand,
% held at the activity's own task: the task keeps its thread for the whole
% hostdem+thinktime interval, so it serializes against the task multiplicity,
% but the host processor is released for it. This mirrors lqns, whose
% think-time attribute LINE already writes out in writeXML.
%
% lqn.actthink_mean is preallocated as NaN by LayeredNetwork/getStruct and stays
% NaN for an activity that was never given a think time, so the value is
% filtered here rather than added by the callers: a NaN reaching servt or residt
% propagates into the layer solvers, and the NaN-ignoring fallbacks there then
% launder it into a plausible-looking bare-service figure instead of failing.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

zt = 0;
if ~isfield(lqn, 'actthink_mean') && ~isprop(lqn, 'actthink_mean')
    return;
end
if aidx < 1 || aidx > numel(lqn.actthink_mean)
    return;
end
v = lqn.actthink_mean(aidx);
if ~isnan(v) && v > GlobalConstants.FineTol
    zt = v;
end
end
