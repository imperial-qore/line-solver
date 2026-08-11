function [T, Y] = lsoda_accurate_stiff(odefun, tspan, y0, options)
% LSODA_ACCURATE_STIFF  Accurate LSODA variant for stiff problems
%
%   Equivalent to JAR's accurateStiffODESolver: LSODA(minstep, maxstep,
%   tol, tol, 12, 5). Full Adams order (12) and BDF order (5) for
%   maximum accuracy.
%
%   See also: lsoda_solve, lsoda_fast_stiff, lsoda_accurate

    if nargin < 4, options = struct(); end
    if ~isfield(options, 'MaxOrdNonStiff'), options.MaxOrdNonStiff = 12; end
    if ~isfield(options, 'MaxOrdStiff'), options.MaxOrdStiff = 5; end
    [T, Y] = lsoda_solve(odefun, tspan, y0, options);
end
