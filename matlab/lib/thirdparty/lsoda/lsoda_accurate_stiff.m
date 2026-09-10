function [T, Y] = lsoda_accurate_stiff(odefun, tspan, y0, options)
% LSODA_ACCURATE_STIFF  Accurate LSODA variant for stiff problems, BDF pinned
%
%   The stiff accurate slot of options.odesolvers, i.e. the LSODA counterpart
%   of @ode15s, and the slot the fluid solver reaches by default
%   (options.stiff is true). Order budget 12/5, as the JAR's
%   accurateStiffODESolver, with the BDF half PINNED (forceStiff).
%
%   See also: lsoda_odesolve, lsoda_fast_stiff, lsoda_accurate

    if nargin < 4, options = struct(); end
    [T, Y] = lsoda_odesolve(odefun, tspan, y0, options, 12, 5, true);
end
