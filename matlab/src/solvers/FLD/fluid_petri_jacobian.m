function [A, Jr] = fluid_petri_jacobian(terms, x, s2, phi, th)
% [A, JR] = FLUID_PETRI_JACOBIAN(TERMS, X, S2, PHI)
% [A, JR] = FLUID_PETRI_JACOBIAN(TERMS, X, S2, PHI, TH)
%
% Drift Jacobian of the fluid Petri net, A = D * dR/dX.
%
% This is what the Lyapunov equation of the linear noise approximation is
% written about, so it has to be the derivative of the SAME rate vector
% FLUID_PETRI_RATES returns, not a finite-difference stand-in: a covariance
% solved about an inconsistent Jacobian is not the covariance of anything.
%
%   kind 1, single phase : rateBase * ( dTHETA*dep + THETA*dDEP )
%   kind 1/2, multi phase: rateBase * e_y(j,h), since the rate is LINEAR in the
%                          servers sitting in that phase and reads no marking
%   kind 3               : zero, the arrival rate is a constant
%   kind 4/5             : zero, an immediate flow and a server latch are
%                          independent unknowns of the DAE, not functions of x
%
% THE VARIANCES ARE HELD, as they are in FLUID_PETRI_THETA and in the queueing
% twin FLUID_DRIFT_JACOBIAN: the closure covariance is pinned by its own
% consistency row, so it does not move along a derivative with respect to x.
%
% Parameters:
%   terms - FLUID_PETRI_TERMS output
%   x     - state vector
%   s2    - the closure covariance entries
%   phi   - immediate firing flows (unused in the derivative, kept for symmetry)
%   th    - optional precomputed FLUID_PETRI_THETA result
%
% Returns:
%   A  - (nstate x nstate) drift Jacobian
%   Jr - (nev x nstate) rate Jacobian
%
% See also FLUID_PETRI_RATES, FLUID_LYAPUNOV, FLUID_DRIFT_JACOBIAN.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5 || isempty(th)
    th = fluid_petri_theta(terms, x, s2);
end

n = terms.nstate;
Jr = zeros(terms.nev, n);
for e = 1:terms.nev
    k = terms.evKind(e);
    if k == 3 || k == 4 || k == 5
        continue
    end
    j = terms.evMode(e);
    md = terms.modes(j);
    base = terms.rateBase(e);
    if md.nph == 1
        Jr(e, th.dslot{j}) = Jr(e, th.dslot{j}) + base * th.dep(j) * th.dval{j}.';
        if ~isempty(th.depslot{j})
            Jr(e, th.depslot{j}) = Jr(e, th.depslot{j}) + base * th.theta(j) * th.depval{j}.';
        end
    else
        zc = md.zblk(terms.evPhase(e));
        Jr(e, zc) = Jr(e, zc) + base;
    end
end
A = terms.D * Jr;
end
