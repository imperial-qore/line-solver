function RD = getCdfPT(self, R)
% RD = GETCDFPT(R) Returns cumulative distribution function of steady-state
% passage times
%
% Backward-compatible name for getCdfPassT, to which this delegates. See
% getCdfPassT for the contract.
%
% @param self SolverFLD instance
% @param R (optional) Response time computation handles
%
% @return RD {nstations x nclasses} cell array of [n x 2] matrices holding
%         CDF values and the corresponding passage times
%
% @see getCdfPassT - Steady-state passage time distribution

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    RD = self.getCdfPassT();
else
    RD = self.getCdfPassT(R);
end
end
