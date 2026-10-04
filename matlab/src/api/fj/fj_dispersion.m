%{ @file fj_dispersion.m
 %  @brief Mean subtask dispersion of a split-merge system with Erlang branches
 %
 %  @author LINE Development Team
%}

%{
 % @brief Mean subtask dispersion of a split-merge system with Erlang branches
 %
 % @details
 % Subtask dispersion is the interval between the first and the last subtask
 % completion of a request. For heterogeneous branches with distribution
 % functions F_i shifted by deterministic delays d_i, the heterogeneous order
 % statistics of Eq. (24) give
 %
 %   E[X_(N)] = integral_0^inf [ 1 - prod_i F_i(x - d_i) ] dx,
 %   E[X_(1)] = integral_0^inf prod_i [ 1 - F_i(x - d_i) ] dx,
 %   E[D_d]   = E[X_(N)] - E[X_(1)]
 %            = integral_0^inf [ 1 - prod_i F_i(x-d_i) - prod_i (1-F_i(x-d_i)) ] dx.
 %
 % The integrand of the last line is what is evaluated; it is non-negative and
 % vanishes at both ends, whereas the difference of the two products printed in
 % the survey is not the dispersion and can go negative.
 %
 % Branch i is an Erlang with shape(i) stages of rate rate(i), which is the
 % split-merge equivalent used in the delay-scheduling construction: a subtask
 % with q other subtasks ahead of it in its parallel queue behaves as an
 % Erlang(q+1, mu). Quadrature is composite Simpson on a horizon widened until
 % the completion probability is within TOL of one.
 %
 % @par Syntax:
 % @code
 % Edisp = fj_dispersion(shape, rate)
 % [Edisp, Emax, Emin] = fj_dispersion(shape, rate, d)
 % [Edisp, Emax, Emin] = fj_dispersion(shape, rate, d, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>shape<td>Vector of N Erlang stage counts (positive integers)
 % <tr><td>rate<td>Vector of N Erlang stage rates (positive)
 % <tr><td>d<td>Vector of N non-negative deterministic delays (optional, default zeros)
 % <tr><td>options<td>Optional struct with fields tol and npoints
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Edisp<td>Mean subtask dispersion
 % <tr><td>Emax<td>Mean completion time of the last subtask
 % <tr><td>Emin<td>Mean completion time of the first subtask
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Eqs. (24)-(25)
 % and Section 6.2.
 %
 % Original: I. Tsimashenka, W. J. Knottenbelt, "Reduction of Subtask
 % Dispersion in Fork-Join Systems", EPEW, 2013.
%}
function [Edisp, Emax, Emin] = fj_dispersion(shape, rate, d, options)

shape = shape(:)';
rate = rate(:)';
N = numel(shape);

if nargin < 3 || isempty(d)
    d = zeros(1, N);
end
d = d(:)';
if nargin < 4 || isempty(options)
    options = struct();
end
if ~isfield(options, 'tol')
    options.tol = 1e-10;
end
if ~isfield(options, 'npoints')
    options.npoints = 4001;
end

if numel(rate) ~= N || numel(d) ~= N
    line_error(mfilename, 'shape, rate and d must have the same length. Got %d, %d, %d.', N, numel(rate), numel(d));
end
if any(shape < 1) || any(shape ~= round(shape))
    line_error(mfilename, 'Erlang stage counts must be positive integers.');
end
if any(rate <= 0)
    line_error(mfilename, 'Erlang stage rates must be positive.');
end
if any(d < 0)
    line_error(mfilename, 'Delays must be non-negative.');
end

% Horizon: widen until every shifted branch is essentially complete
U = max(d) + 8 * max(shape ./ rate);
for it = 1:60
    prodF = 1;
    for i = 1:N
        prodF = prodF * erlang_cdf(U - d(i), shape(i), rate(i));
    end
    if 1 - prodF < options.tol
        break
    end
    U = 2 * U;
end

npoints = options.npoints;
if mod(npoints, 2) == 0
    npoints = npoints + 1;
end
x = linspace(0, U, npoints);
h = x(2) - x(1);
w = ones(1, npoints);
w(2:2:end-1) = 4;
w(3:2:end-2) = 2;

Fprod = ones(1, npoints);
Sprod = ones(1, npoints);
for i = 1:N
    Fi = erlang_cdf(x - d(i), shape(i), rate(i));
    Fprod = Fprod .* Fi;
    Sprod = Sprod .* (1 - Fi);
end

Emax = (h / 3) * sum(w .* (1 - Fprod));
Emin = (h / 3) * sum(w .* Sprod);
Edisp = Emax - Emin;

end

function F = erlang_cdf(t, k, mu)
% Erlang-k distribution function, zero on the negative half line
F = zeros(size(t));
pos = t > 0;
tp = t(pos);
acc = zeros(size(tp));
term = ones(size(tp));
for j = 0:(k - 1)
    if j > 0
        term = term .* (mu * tp) / j;
    end
    acc = acc + term;
end
F(pos) = 1 - exp(-mu * tp) .* acc;
end
