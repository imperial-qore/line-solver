function r = solver_ctmc_ratecomplement(D, nonimm, imm, Q12, Q22)
% R = SOLVER_CTMC_RATECOMPLEMENT(D, NONIMM, IMM, Q12, Q22)
%
% Long-run rate of an action, as seen from each tangible (non-vanishing) state,
% given the action's rate filter D over the full state space.
%
% Vanishing states are removed from the generator by stochastic complementation,
% so an action that fires only in vanishing states (a fork firing, a join
% departure, or the firing of an immediate SPN mode) would be lost if its rate
% were read off the tangible rows alone. The rate observed from tangible state s
% is the direct exit rate via the action plus the expected number of firings
% along the vanishing chain entered from s:
%
%   r = D(nonimm,:)*1 + Q12*(-Q22)^(-1)*(D(imm,:)*1)
%
% where Q12 and Q22 are the tangible-to-vanishing and vanishing-to-vanishing
% blocks returned by ctmc_stochcomp.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

r = full(sum(D(nonimm,:),2));
if ~isempty(imm)
    r = r + Q12*((-Q22) \ full(sum(D(imm,:),2)));
end
end
