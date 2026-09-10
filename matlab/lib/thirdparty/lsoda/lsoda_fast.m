function [T, Y] = lsoda_fast(odefun, tspan, y0, options)
% LSODA_FAST  Fast LSODA variant, limited Adams/BDF order
%
%   The non-stiff fast slot of options.odesolvers, i.e. the LSODA counterpart
%   of @ode23. Equivalent to the JAR's fastODESolver order budget: Adams and
%   BDF order 3. Auto-switching is left on, as this slot is chosen where the
%   caller has not asked for a stiff integrator.
%
%   See also: lsoda_odesolve, lsoda_accurate, lsoda_fast_stiff

    if nargin < 4, options = struct(); end
    [T, Y] = lsoda_odesolve(odefun, tspan, y0, options, 3, 3, false);
end
