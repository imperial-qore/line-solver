%{ @file fj_xmax_hz_het.m
 %  @brief Harrison-Zertal approximation of the maximum of general variables
 %
 %  @author LINE Development Team
%}

%{
 % @brief Harrison-Zertal approximation of the maximum of general variables
 %
 % @details
 % Expected maximum of K independent but not identically distributed
 % non-negative variables, each supplied through its first two moments and its
 % distribution function. Writing S for a sub-collection and alpha_i = 1/m1_i,
 % the recurrence averages, over which branch is singled out, the expected
 % maximum of the remaining branches plus the residual life still owed by the
 % branch singled out:
 %
 %   I(S) = (1/|S|) * sum_{i in S} [ I(S \ i)
 %                                   + (m2_i / (2*m1_i)) * L*_{S\i}(alpha_i) ],
 %
 % anchored at I({i}) = m1_i. The transform of the maximum over a
 % sub-collection is recovered from the product of the distribution functions,
 %
 %   L*_T(s) = s * integral_0^inf exp(-s t) * prod_{j in T} F_j(t) dt,
 %
 % evaluated by composite Simpson quadrature on a horizon widened until the
 % product of the distribution functions is within TOL of one.
 %
 % For identically distributed branches the recurrence collapses onto the
 % closed form of fj_xmax_hz, and for identical exponential branches it is
 % exact at H_K/lambda.
 %
 % @par Syntax:
 % @code
 % Xmax = fj_xmax_hz_het(m1, m2, cdf)
 % Xmax = fj_xmax_hz_het(m1, m2, cdf, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>m1<td>Vector of K branch means (positive)
 % <tr><td>m2<td>Vector of K branch second moments, m2(i) >= m1(i)^2
 % <tr><td>cdf<td>Cell array of K function handles, cdf{i}(t) = P(X_i <= t)
 % <tr><td>options<td>Optional struct with fields tol and npoints
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Xmax<td>Approximate expected maximum
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Eq. (46).
 %
 % Original: P. G. Harrison, S. Zertal, "Queueing Models of RAID Systems with
 % Maxima of Waiting Times", Performance Evaluation 64(7-8), 2007.
%}
function Xmax = fj_xmax_hz_het(m1, m2, cdf, options)

if nargin < 4 || isempty(options)
    options = struct();
end
if ~isfield(options, 'tol')
    options.tol = 1e-10;
end
if ~isfield(options, 'npoints')
    options.npoints = 2001;
end

m1 = m1(:)';
m2 = m2(:)';
K = numel(m1);

if numel(m2) ~= K || numel(cdf) ~= K
    line_error(mfilename, 'm1, m2 and cdf must have the same length. Got %d, %d, %d.', K, numel(m2), numel(cdf));
end
if any(m1 <= 0)
    line_error(mfilename, 'All branch means must be positive.');
end
if any(m2 < m1.^2)
    line_error(mfilename, 'Some second moment is below the square of its mean.');
end
if K > 14
    line_error(mfilename, 'The recurrence enumerates 2^K sub-collections with a quadrature each; K=%d is too large.', K);
end

alpha = 1 ./ m1;
resid = m2 ./ (2 * m1);

% Quadrature horizon: widen until every branch is essentially complete
U = 8 * max(m1);
for it = 1:60
    prodF = 1;
    for j = 1:K
        prodF = prodF * cdf{j}(U);
    end
    if 1 - prodF < options.tol
        break
    end
    U = 2 * U;
end

nmask = 2^K;
Ival = zeros(1, nmask);
for mask = 1:(nmask - 1)
    bits = bitget(mask, 1:K);
    members = find(bits == 1);
    card = numel(members);
    if card == 1
        Ival(mask + 1) = m1(members);
        continue
    end
    acc = 0;
    for i = members
        rest = bitxor(mask, bitshift(1, i - 1));
        acc = acc + Ival(rest + 1) + resid(i) * lst_max(cdf, rest, K, alpha(i), U, options.npoints);
    end
    Ival(mask + 1) = acc / card;
end

Xmax = Ival(nmask);

end

function L = lst_max(cdf, mask, K, s, U, npoints)
% L*_T(s) = s * integral_0^U exp(-s t) * prod_{j in T} F_j(t) dt by Simpson
if mask == 0
    L = 1;
    return
end
if mod(npoints, 2) == 0
    npoints = npoints + 1;
end
t = linspace(0, U, npoints);
h = t(2) - t(1);
g = exp(-s * t);
bits = bitget(mask, 1:K);
for j = find(bits == 1)
    g = g .* cdf{j}(t);
end
w = ones(1, npoints);
w(2:2:end-1) = 4;
w(3:2:end-2) = 2;
L = s * (h / 3) * sum(w .* g);
end
