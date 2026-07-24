function [T, Y] = lsoda_fast(odefun, tspan, y0, options)
% LSODA_FAST  Fast LSODA variant with low-order Adams/BDF methods
%
%   Equivalent to JAR's fastOdeSolver: limited Adams order for speed.
%   Uses MaxOrdNonStiff=3 to limit nonstiff method order (similar to
%   DormandPrince54 used in JAR for non-stiff problems).
%
%   See also: lsoda_solve, lsoda_accurate, lsoda_fast_stiff

    if nargin < 4, options = struct(); end
    if ~isfield(options, 'MaxOrdNonStiff'), options.MaxOrdNonStiff = 3; end
    if ~isfield(options, 'MaxOrdStiff'), options.MaxOrdStiff = 3; end
    [T, Y] = lsoda_solve(odefun, tspan, y0, options);
end
