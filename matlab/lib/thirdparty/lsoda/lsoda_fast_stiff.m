function [T, Y] = lsoda_fast_stiff(odefun, tspan, y0, options)
% LSODA_FAST_STIFF  Fast LSODA variant for stiff problems, BDF pinned
%
%   The stiff fast slot of options.odesolvers, i.e. the LSODA counterpart of
%   @ode23s. Order budget 3/3, as the JAR's fastStiffODESolver. The BDF half is
%   PINNED (forceStiff): a slot the caller reached by asking for a stiff
%   integrator must not be answered by the Adams half, which loses its
%   stability bound at a fixed point. See lsoda_matlab for the mechanism.
%
%   See also: lsoda_odesolve, lsoda_accurate_stiff, lsoda_fast

    if nargin < 4, options = struct(); end
    [T, Y] = lsoda_odesolve(odefun, tspan, y0, options, 3, 3, true);
end
