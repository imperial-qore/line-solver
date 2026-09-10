function [A, phi, esup] = mam_bgchain_envfull(bg, i)
% [A, PHI, ESUP] = MAM_BGCHAIN_ENVFULL(BG, I)
%
% The background modulating chain BG itself, as the environment station I sees,
% WITHOUT the lumping MAM_BGCHAIN_ENV applies.
%
% THIS IS AN ORACLE, NOT AN IMPROVEMENT, AND THE DIFFERENCE WAS MEASURED.
% MAM_BGCHAIN_ENV aggregates BG.Q over the level sets {s : totocc(s,i) = e},
% which is exact just when the partition is lumpable in Kemeny-Snell's sense,
% and nothing in the method checks that it is. It looks like the one
% uncontrolled approximation left in the open half. IT IS NOT ONE: the
% background chain is PRODUCT-FORM BY CONSTRUCTION -- its rates are
% CSHARE(i,e+1) * (n_{i,b}/e) * mu_{i,b}, a load-dependent capacity shared in
% proportion to the class counts, under Markov routing -- so given n_i = e the
% conditional law of the OTHER stations is the product-form law over the
% remaining population, whatever the tagged station's occupancy has been doing.
% The pi-weighted aggregation is then Norton-exact, and the level-dependent
% rescaling (a function of e alone) preserves it. Measured: 6.7e-16 between the
% two on a 3-station single-chain cycle whose up-rates are 0/1/2 within one
% class, 6.7e-16 on a two-background-class chain whose down-rates are 5 and 0.5
% within one class, and no difference at any printed digit on four driver-level
% models checked against SolverCTMC. An INDEPENDENT explicit CTMC of the
% modulated queue agrees to 1.2e-11.
%
% The lumping being exact is a property of the CHAIN, not of the QBD: fed an
% environment that is not a product-form network the QBD does separate the two
% (2.9% on a 4-state counterexample). So keep this function as the check that
% the property still holds if the chain's rate law ever changes -- and do not
% expect it to move a number today. It exists so that the next person to suspect
% the lumping can settle it in one run instead of rebuilding this.
%
% This function returns the chain unlumped, so the QBD's environment axis
% carries the whole closed population vector and the closed dynamics it sees are
% exactly BG.Q. Two things follow:
%
%   - the environment is Markovian by CONSTRUCTION rather than by assumption,
%     and transitions between states that hold the SAME number of closed jobs at
%     station i -- closed jobs moving among the OTHER stations -- exist here and
%     are simply absent from the lumped chain;
%   - the phase count of the QBD grows from (distinct occupancies at i) to
%     (states of the background chain), which is what OPTIONS.CONFIG.QBDPHASES_MAX
%     is there to bound.
%
% What it does NOT remove is the mean-field coupling to the OTHER stations: the
% rates of BG.Q at station j were built at j's capacity share CSHARE(j,.),
% already averaged over j's open occupancy. This is exact in the closed state
% at station I, and mean-field elsewhere.
%
% Unreachable states are dropped, as MAM_BGCHAIN_ENV drops unreachable
% environment states, and the diagonal is rebuilt so what is returned is a
% generator over the states that remain.
%
% Outputs
%   A     (me x me)  generator of the background chain, restricted to its support
%   PHI   (1 x me)   stationary probability of each state
%   ESUP  (1 x me)   closed jobs station i holds in each state; UNLIKE
%                    MAM_BGCHAIN_ENV's, this vector REPEATS -- many chain states
%                    hold the same number of jobs at station i, which is exactly
%                    the information the lumping throws away
%
% See also MAM_BGCHAIN_ENV, MAM_BGCHAIN_CTMC, SOLVER_MAM_BGCHAIN.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

occ = bg.totocc(:, i);
pi0 = bg.pi(:)';

keep = pi0 > GlobalConstants.Zero;
if ~any(keep)
    % degenerate chain: nothing reachable
    A = 0;
    phi = 1;
    esup = 0;
    return;
end

A = full(bg.Q(keep, keep));
esup = occ(keep)';
phi = pi0(keep) / sum(pi0(keep));

me = numel(esup);
if me == 1
    A = 0;
    return;
end

A(1:me+1:end) = 0;
A = A - diag(sum(A, 2));
end
