function u = fluid_dae_project(u, lb)
% U = FLUID_DAE_PROJECT(U, LB)
%
% An iterate projected onto its feasible box, the lower bound only.
%
% WHY PROJECTED, AND NOT MERELY CLAMPED INSIDE THE RESIDUAL. The bounded
% unknowns of the DAE -- a variance, an admission throttle, an immediate firing
% flow -- are read through max(0,.) by the residual that uses them. Clamping
% there ALONE is a trap: once an iterate goes negative the residual stops
% depending on it, so the finite-difference column is exactly zero, the solver
% has no derivative to climb back on, and the unknown is pinned at the boundary
% for good. That is not hypothetical -- it pinned the throttle at zero in every
% configuration of the capacity sweep, and the region then settled BELOW its cap
% instead of on it, while the residual reported the constraint as the only unmet
% equation.
%
% Projecting the ITERATE instead keeps every evaluation inside the feasible box,
% where the forward difference across max(0,.) is live even exactly at zero, so
% the boundary can be left again.
%
% Parameters:
%   u  - the iterate
%   lb - lower bound per unknown, one entry per unknown; -Inf for a free one.
%        It is a VECTOR, never a count. The earlier form also accepted a scalar
%        NFREE meaning "the first NFREE unknowns are free", which is ambiguous
%        the moment U holds a single unknown: a one-place, one-class net with no
%        closure pair produced lb = -Inf, that was read as NFREE = -Inf, and the
%        projection indexed U(-Inf+1:end).
%
% Returns:
%   u - the projected iterate
%
% See also FLUID_DAE_NEWTON, SOLVER_FLUID_DAE, SOLVER_FLUID_PETRI.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(lb)
    return
end
lb = lb(:);
if numel(lb) ~= numel(u)
    line_error(mfilename, sprintf('The bound vector has %d entries for %d unknowns.', numel(lb), numel(u)));
end
fin = isfinite(lb);
u(fin) = max(lb(fin), u(fin));
end
