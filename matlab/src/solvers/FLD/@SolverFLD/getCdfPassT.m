function RD = getCdfPassT(self, R)
% RD = GETCDFPASST(R) Returns cumulative distribution function of steady-state
% passage times
%
% @brief Computes the steady-state passage time distribution for each station
% and class
%
% The passage time of a job of class r at station i is the time from its
% arrival at the station to its departure from it. Under the fluid
% approximation this passage is exactly the quantity reported by
% getCdfRespT: both are obtained from solver_fluid_passage_time started from
% the steady-state ODE solution, so this method delegates to it rather than
% duplicating the computation. The two names are kept distinct because
% NetworkSolver declares both, and other solvers may separate them.
%
% For the passage time along a prescribed route, i.e. conditional on a given
% sequence of nodes rather than at a single station, no method is provided:
% the fluid passage time analysis is per station and combining stations would
% require an independence assumption across them.
%
% @param self SolverFLD instance
% @param R (optional) Response time computation handles. If omitted they are
%        obtained from getAvgRespTHandles.
%
% @return RD {nstations x nclasses} cell array. Each non-empty RD{i,r} is an
%         [n x 2] matrix whose first column holds CDF values in [0,1] and
%         whose second column holds the corresponding passage times.
%
% @see getCdfRespT - Steady-state response time distribution (same quantity)
% @see getTranCdfPassT - Passage time distribution during the transient

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    RD = self.getCdfRespT();
else
    RD = self.getCdfRespT(R);
end
end
