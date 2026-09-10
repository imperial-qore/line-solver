%{ @file fj_xmax_moments_het.m
 %  @brief Moments of the maximum of heterogeneous exponentials by recurrence
 %
 %  @author LINE Development Team
%}

%{
 % @brief Moments of the maximum of heterogeneous exponentials by recurrence
 %
 % @details
 % Exact moments of Y = max(X_1,...,X_K) with independent X_i ~ Exp(lambda_i),
 % obtained by differentiating the Harrison and Zertal transform recurrence n
 % times and evaluating at the origin:
 %
 %   M_m(lambda, n) = [ n * M_m(lambda, n-1)
 %                      + sum_{j=1..m} lambda_j * M_{m-1}(lambda \ j, n) ]
 %                    / sum_{j=1..m} lambda_j,
 %
 % with M_m(lambda, 0) = 1 and M_0(., n) = 0 for n >= 1. Every moment up to
 % order n is produced by the same sweep, so all of them are returned.
 %
 % Eq. (30) of the survey prints the second sum WITHOUT the lambda_j weight.
 % That form is not the derivative of Eq. (29) and does not reproduce the
 % textbook two-variable answer 1/l1 + 1/l2 - 1/(l1+l2); the weight is
 % restored here, and fj_xmax_het gives the independent inclusion-exclusion
 % check.
 %
 % @par Syntax:
 % @code
 % M = fj_xmax_moments_het(lambda)
 % M = fj_xmax_moments_het(lambda, n)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>Vector of K positive exponential rates
 % <tr><td>n<td>Highest moment order (optional, default 1)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>M<td>Row vector of moments of orders 1..n of the maximum
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Eq. (30).
 %
 % Original: P. G. Harrison, S. Zertal, "Queueing Models of RAID Systems with
 % Maxima of Waiting Times", Performance Evaluation 64(7-8), 2007.
%}
function M = fj_xmax_moments_het(lambda, n)

if nargin < 2 || isempty(n)
    n = 1;
end

lambda = lambda(:)';
K = numel(lambda);

if K < 1
    line_error(mfilename, 'lambda must contain at least one rate.');
end
if any(lambda <= 0)
    line_error(mfilename, 'All exponential rates must be positive.');
end
if n < 1 || n ~= round(n)
    line_error(mfilename, 'Moment order n must be a positive integer. Got n=%g.', n);
end
if K > 22
    line_error(mfilename, 'The recurrence enumerates 2^K sub-collections; K=%d is too large.', K);
end

nmask = 2^K;
% tab(mask+1, k+1) holds the k-th moment of the maximum over the sub-collection
% selected by mask; order 0 is one for every sub-collection, empty or not
tab = zeros(nmask, n + 1);
tab(:, 1) = 1;

for order = 1:n
    % The empty sub-collection has a zero maximum, so all its moments vanish
    tab(1, order + 1) = 0;
    for mask = 1:(nmask - 1)
        bits = bitget(mask, 1:K);
        members = find(bits == 1);
        tot = sum(lambda(members));
        acc = order * tab(mask + 1, order);
        for j = members
            acc = acc + lambda(j) * tab(bitxor(mask, bitshift(1, j - 1)) + 1, order + 1);
        end
        tab(mask + 1, order + 1) = acc / tot;
    end
end

M = tab(nmask, 2:(n + 1));

end
