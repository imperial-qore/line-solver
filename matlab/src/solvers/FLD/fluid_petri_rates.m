function [r, th] = fluid_petri_rates(terms, x, s2, phi, mu, th)
% [R, TH] = FLUID_PETRI_RATES(TERMS, X, S2, PHI, MU)
% [R, TH] = FLUID_PETRI_RATES(TERMS, X, S2, PHI, MU, TH)
%
% The rate of every event column of the fluid Petri drift.
%
%   kind 1  firing of mode j, phase h -> h'
%             single phase : rateBase * theta_j * dep_j
%             multi  phase : rateBase * y(j,h)
%   kind 2  internal phase change of mode j: rateBase * y(j,h)
%   kind 3  exogenous arrival: rateBase, a constant
%   kind 4  firing of an IMMEDIATE mode: PHI(j), an algebraic unknown
%   kind 5  the server latch of a multi-phase mode: MU(j), a free-sign
%           algebraic unknown depositing at the firing process's entry vector
%
% THETA_J IS THE NUMBER OF RUNNING SERVERS a marking supports, closed by
% FLUID_PETRI_THETA. A SINGLE-PHASE mode holds no coordinate of its own, so its
% rate is theta_j directly. A MULTI-PHASE mode carries one coordinate per phase
% holding the servers sitting in it, and its rates are LINEAR in those: the
% marking reaches them through the latch row sum_h y(j,h) = theta_j, not through
% the rate.
%
% THE MARKING-DEPENDENT MULTIPLIER APPLIES TO A FIRING ONLY, and only to a
% single-phase mode: SETFIRINGRATEDEPENDENCE accepts it for exponential timing
% alone, so there is no multi-phase clock for it to speed up.
%
% Parameters:
%   terms - FLUID_PETRI_TERMS output
%   x     - state vector
%   s2    - the closure covariance entries
%   phi   - (nimm x 1) immediate firing flows, in TERMS.immIdx order
%   mu    - (nlatch x 1) server-latch flows, in TERMS.latchMode order
%   th    - optional precomputed FLUID_PETRI_THETA result
%
% Returns:
%   r  - (nev x 1) rates
%   th - the enabling terms, so a caller taking the Jacobian next reuses them
%
% See also FLUID_PETRI_THETA, FLUID_PETRI_JACOBIAN, SOLVER_FLUID_PETRI.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 6 || isempty(th)
    th = fluid_petri_theta(terms, x, s2);
end
if nargin < 4 || isempty(phi)
    phi = zeros(numel(terms.immIdx),1);
end
if nargin < 5 || isempty(mu)
    mu = zeros(numel(terms.latchMode),1);
end

r = zeros(terms.nev,1);
immPos = zeros(numel(terms.modes),1);
immPos(terms.immIdx) = 1:numel(terms.immIdx);
latchPos = zeros(numel(terms.modes),1);
latchPos(terms.latchMode) = 1:numel(terms.latchMode);

for e = 1:terms.nev
    k = terms.evKind(e);
    switch k
        case 3
            r(e) = terms.rateBase(e);
        case 4
            r(e) = phi(immPos(terms.evMode(e)));
        case 5
            r(e) = mu(latchPos(terms.evMode(e)));
        otherwise
            j = terms.evMode(e);
            md = terms.modes(j);
            if md.nph == 1
                r(e) = terms.rateBase(e) * th.theta(j) * th.dep(j);
            else
                r(e) = terms.rateBase(e) * x(md.zblk(terms.evPhase(e)));
            end
    end
end
end
