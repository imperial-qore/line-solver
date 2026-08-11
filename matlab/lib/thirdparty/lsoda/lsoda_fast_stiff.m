function [T, Y] = lsoda_fast_stiff(odefun, tspan, y0, options)
% LSODA_FAST_STIFF  Fast LSODA variant optimized for stiff problems
%
%   Equivalent to JAR's fastStiffODESolver: LSODA(minstep, maxstep, tol,
%   tol, 3, 3). Limited order for both Adams (3) and BDF (3) for speed.
%
%   See also: lsoda_solve, lsoda_accurate_stiff, lsoda_fast

    if nargin < 4, options = struct(); end
    if ~isfield(options, 'MaxOrdNonStiff'), options.MaxOrdNonStiff = 3; end
    if ~isfield(options, 'MaxOrdStiff'), options.MaxOrdStiff = 3; end
    [T, Y] = lsoda_solve(odefun, tspan, y0, options);
end
