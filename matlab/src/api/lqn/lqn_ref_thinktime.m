function z = lqn_ref_thinktime(lqn, tidx)
% Z = LQN_REF_THINKTIME(LQN, TIDX)
%
% Declared think time of task TIDX as it enters the thread cycle: the value for
% a REFERENCE task, zero for any other.
%
% A think time is an attribute of the closed customer population a reference
% task stands for, and it is what separates one request of that population from
% the next. On a served task it has no such meaning, and charging it per request
% throttles the task: `lqn_basic`'s T3 has 25 threads and a declared think time
% of 4, and reading it as a per-request delay caps it at 25/(4+0.02) = 6.219
% completions per second. Three independent oracles refuse that reading and put
% the rate at 5 calls per caller request instead -- lqsim 66.5, LDES 66.955,
% lqns 75.6 -- so the think time of a non-reference task does not enter the
% cycle. BOTH SolverLN methods go through this gate: `updateThinkTimes` and
% `updateThinkTimesPH` for the surrogate delay, `getTranAvgCoupled` for its
% trajectory, and `buildLayersRecursive` for the build-time seed. See
% _kb/06-solver-catalog.md (LN section).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

z = 0;
if tidx < 1 || tidx > numel(lqn.think) || isempty(lqn.think{tidx})
    return
end
if ~lqn.isref(tidx)
    return
end
z = lqn.think{tidx}.getMean();
if ~isfinite(z) || z < 0
    z = 0;
end
end
