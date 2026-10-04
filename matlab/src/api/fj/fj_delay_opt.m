%{ @file fj_delay_opt.m
 %  @brief Deterministic subtask delays that minimise mean dispersion
 %
 %  @author LINE Development Team
%}

%{
 % @brief Deterministic subtask delays that minimise mean dispersion
 %
 % @details
 % Chooses the vector of deterministic delays d = (d_1,...,d_N) that minimises
 % the mean subtask dispersion E[D_d] of fj_dispersion. Holding back a fast
 % branch costs little at the last completion and buys a great deal at the
 % first, so the minimiser is generally interior and strictly positive on every
 % branch but the slowest.
 %
 % The objective is minimised by cyclic coordinate descent with a golden
 % section line search on each coordinate, which is deterministic and needs no
 % derivative; the delay of the branch with the largest mean is pinned at zero,
 % because adding a constant to every delay shifts both order statistics
 % equally and leaves the dispersion unchanged. The search therefore returns
 % the representative with min(d) = 0.
 %
 % @par Syntax:
 % @code
 % d = fj_delay_opt(shape, rate)
 % [d, Edisp, Emax] = fj_delay_opt(shape, rate, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>shape<td>Vector of N Erlang stage counts (positive integers)
 % <tr><td>rate<td>Vector of N Erlang stage rates (positive)
 % <tr><td>options<td>Optional struct with fields tol, npoints, maxsweeps and dtol
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>d<td>Vector of N optimal non-negative delays with min(d) = 0
 % <tr><td>Edisp<td>Mean dispersion attained at d
 % <tr><td>Emax<td>Mean completion time of the last subtask at d, the price paid
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Section 6.2.
 %
 % Original: I. Tsimashenka, W. J. Knottenbelt, "Reduction of Subtask
 % Dispersion in Fork-Join Systems", EPEW, 2013.
%}
function [d, Edisp, Emax] = fj_delay_opt(shape, rate, options)

shape = shape(:)';
rate = rate(:)';
N = numel(shape);

if nargin < 3 || isempty(options)
    options = struct();
end
if ~isfield(options, 'tol')
    options.tol = 1e-10;
end
if ~isfield(options, 'npoints')
    options.npoints = 2001;
end
if ~isfield(options, 'maxsweeps')
    options.maxsweeps = 40;
end
if ~isfield(options, 'dtol')
    options.dtol = 1e-8;
end

if numel(rate) ~= N
    line_error(mfilename, 'shape and rate must have the same length. Got %d and %d.', N, numel(rate));
end
if N < 2
    d = zeros(1, N);
    [Edisp, Emax] = fj_dispersion(shape, rate, d, options);
    return
end

means = shape ./ rate;
% Delaying past the slowest branch never helps, so that is the search ceiling
ub = max(means) + 8 * max(sqrt(shape) ./ rate);

d = zeros(1, N);
fcur = objective(shape, rate, d, options);
for sweep = 1:options.maxsweeps
    fprev = fcur;
    for i = 1:N
        [d(i), fcur] = golden(shape, rate, d, i, 0, ub, options);
    end
    % Normalise so that the smallest delay is zero
    d = d - min(d);
    fcur = objective(shape, rate, d, options);
    if abs(fprev - fcur) <= options.dtol * max(1, abs(fprev))
        break
    end
end

[Edisp, Emax] = fj_dispersion(shape, rate, d, options);

end

function f = objective(shape, rate, d, options)
f = fj_dispersion(shape, rate, d, options);
end

function [xbest, fbest] = golden(shape, rate, d, i, lo, hi, options)
% Golden section search on coordinate i, holding the other delays fixed
invphi = (sqrt(5) - 1) / 2;
a = lo;
b = hi;
c = b - invphi * (b - a);
dd = a + invphi * (b - a);
dc = d; dc(i) = c; fc = objective(shape, rate, dc, options);
dv = d; dv(i) = dd; fd = objective(shape, rate, dv, options);
for it = 1:60
    if fc < fd
        b = dd; dd = c; fd = fc;
        c = b - invphi * (b - a);
        dc = d; dc(i) = c; fc = objective(shape, rate, dc, options);
    else
        a = c; c = dd; fc = fd;
        dd = a + invphi * (b - a);
        dv = d; dv(i) = dd; fd = objective(shape, rate, dv, options);
    end
    if (b - a) <= options.dtol * max(1, hi)
        break
    end
end
if fc < fd
    xbest = c; fbest = fc;
else
    xbest = dd; fbest = fd;
end
end
