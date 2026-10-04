function [u, it, converged, resnorm] = fluid_dae_newton(resid, u, tol, maxit, lb)
% [U, IT, CONVERGED, RESNORM] = FLUID_DAE_NEWTON(RESID, U, TOL, MAXIT, LB)
%
% Damped PROJECTED Newton with a finite-difference Jacobian and an Armijo
% backtrack on the residual norm, the solver behind the DAE formulation of the
% fluid closures.
%
% The step is solved in least squares rather than by a square factorisation: the
% drift block is rank deficient by exactly the number of conserved quantities,
% and the conservation rows restore that rank, so the stacked system is
% consistent and overdetermined rather than square. The same holds for the
% Petri route, where the conserved quantities are the net's P-invariants.
%
% Parameters:
%   resid - handle resid(u, quiet) returning the stacked residual, and [] when
%           it cannot be evaluated at that iterate so the line search backs off
%   u     - starting iterate
%   tol   - convergence tolerance on the infinity norm of the residual
%   maxit - maximum Newton steps
%   lb    - lower bound per unknown (-Inf where free), or the scalar NFREE for
%           the layout "the first NFREE unknowns are free, the rest are
%           non-negative"; see FLUID_DAE_PROJECT
%
% Returns:
%   u         - the last iterate
%   it        - Newton steps taken
%   converged - whether RESNORM fell below TOL
%   resnorm   - infinity norm of the residual at U
%
% See also FLUID_DAE_FDJAC, FLUID_DAE_PROJECT, SOLVER_FLUID_DAE, SOLVER_FLUID_PETRI.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5
    lb = [];
end
u = fluid_dae_project(u, lb);
G = resid(u, false);
resnorm = norm(G, inf);
converged = resnorm < tol;
it = 0;
while ~converged && it < maxit
    it = it + 1;
    J = fluid_dae_fdjac(resid, u, G);
    warnstate = warning('off','MATLAB:rankDeficientMatrix');
    du = -(J \ G);
    warning(warnstate);
    if any(~isfinite(du))
        du = -(pinv(J) * G);
    end
    lam = 1;
    stepped = false;
    for ls = 1:25
        % project BEFORE evaluating, so the accepted point and the residual
        % that measured it are the same feasible point
        un = fluid_dae_project(u + lam*du, lb);
        Gn = resid(un, true);
        if ~isempty(Gn) && all(isfinite(Gn)) && norm(Gn,inf) < resnorm*(1 - 1e-4*lam)
            u = un; G = Gn; resnorm = norm(Gn,inf);
            stepped = true;
            break
        end
        lam = lam/2;
    end
    if ~stepped
        break % no descent along this direction; report the last iterate
    end
    converged = resnorm < tol;
end
converged = resnorm < tol;
end
