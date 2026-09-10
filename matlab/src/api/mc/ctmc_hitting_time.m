function h = ctmc_hitting_time(Q, targetStates)
% H = CTMC_HITTING_TIME(Q, TARGETSTATES)
%
% Mean time to reach any state in TARGETSTATES from each state of a CTMC with
% generator Q. Target states have zero hitting time; a state that cannot reach
% the set has an infinite one.
%
% Continuous-time twin of DTMC_HITTING_TIME, and the first-moment special case
% of CTMC_PASSAGE_MOMENTS: (-S) h = 1 on the non-target block, where
% DTMC_HITTING_TIME solves (I - P_NT) h = 1. TARGETSTATES is 1-based here and
% 0-based in the Java, Python and C++ twins.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% MALL does not depend on the initial law, so a uniform one is passed rather
% than paying for the stationary solve that PI0 = [] would trigger.
n = size(Q,1);
mall = ctmc_passage_moments(Q, ones(1,n)/n, targetStates, 1);
h = mall(:,1);
end
