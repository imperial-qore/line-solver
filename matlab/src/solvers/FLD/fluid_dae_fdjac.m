function J = fluid_dae_fdjac(resid, u, G)
% J = FLUID_DAE_FDJAC(RESID, U, G)
%
% Finite-difference Jacobian of a DAE residual at U, given its value G there.
%
% Forward differences by default, and a BACKWARD difference for any column whose
% forward step lands where the residual cannot be evaluated -- the closure
% raises on a non-hyperbolic fixed point, so a step that crosses into that
% region returns empty rather than a number. A column that fails both ways is
% left at zero and the least-squares step absorbs it.
%
% Parameters:
%   resid - handle resid(u, quiet) returning the residual, [] where it cannot
%           be evaluated
%   u     - the iterate
%   G     - the residual at U
%
% Returns:
%   J - (numel(G) x numel(u)) Jacobian
%
% See also FLUID_DAE_NEWTON, SOLVER_FLUID_DAE, SOLVER_FLUID_PETRI.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = numel(u);
m = numel(G);
J = zeros(m,n);
for j = 1:n
    h = max(1e-7*abs(u(j)), 1e-9);
    up = u; up(j) = up(j) + h;
    Gp = resid(up, true);
    if isempty(Gp) || any(~isfinite(Gp))
        um = u; um(j) = um(j) - h;
        Gm = resid(um, true);
        if isempty(Gm) || any(~isfinite(Gm))
            continue % column left at zero; the least-squares step absorbs it
        end
        J(:,j) = (G - Gm)/h;
    else
        J(:,j) = (Gp - G)/h;
    end
end
end
