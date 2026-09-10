function [T, Y] = lsoda_accurate(odefun, tspan, y0, options)
% LSODA_ACCURATE  Accurate LSODA variant, full Adams/BDF order
%
%   The non-stiff accurate slot of options.odesolvers, i.e. the LSODA
%   counterpart of @ode113: full Adams order 12, BDF order 5 for the stiff
%   regions the auto-switcher may enter.
%
%   See also: lsoda_odesolve, lsoda_fast, lsoda_accurate_stiff

    if nargin < 4, options = struct(); end
    [T, Y] = lsoda_odesolve(odefun, tspan, y0, options, 12, 5, false);
end
