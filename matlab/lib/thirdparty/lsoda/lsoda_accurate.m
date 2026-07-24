function [T, Y] = lsoda_accurate(odefun, tspan, y0, options)
% LSODA_ACCURATE  Accurate LSODA variant with full-order Adams/BDF methods
%
%   Equivalent to JAR's accurateOdeSolver: full Adams order (12) for
%   maximum accuracy on nonstiff problems, BDF order 5 for stiff regions.
%
%   See also: lsoda_solve, lsoda_fast, lsoda_accurate_stiff

    if nargin < 4, options = struct(); end
    if ~isfield(options, 'MaxOrdNonStiff'), options.MaxOrdNonStiff = 12; end
    if ~isfield(options, 'MaxOrdStiff'), options.MaxOrdStiff = 5; end
    [T, Y] = lsoda_solve(odefun, tspan, y0, options);
end
