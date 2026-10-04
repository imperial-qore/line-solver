%{ @file fj_xmax_het.m
 %  @brief Exact moments of the maximum of heterogeneous exponentials
 %
 %  @author LINE Development Team
%}

%{
 % @brief Exact moments of the maximum of heterogeneous exponentials
 %
 % @details
 % Exact n-th moment of Y = max(X_1,...,X_K) for independent but not
 % identically distributed exponential variables, X_i ~ Exp(lambda_i), by
 % inclusion-exclusion on the survival function:
 %
 %   E[Y^n] = sum over the nonempty subsets S of {1..K} of
 %              (-1)^(|S|+1) * n! / ( sum_{i in S} lambda_i )^n .
 %
 % For n = 1 and K = 2 this collapses to the textbook
 % 1/l1 + 1/l2 - 1/(l1+l2), and for equal rates to H_K/lambda.
 %
 % The cost is 2^K - 1 terms, so the enumeration is refused beyond K = 24;
 % use fj_xmax_moments_het for the O(2^K) but numerically safer recursion, or
 % fj_xmax_exp when the rates are equal.
 %
 % @par Syntax:
 % @code
 % Xmax = fj_xmax_het(lambda)
 % Mn = fj_xmax_het(lambda, n)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>Vector of K positive exponential rates
 % <tr><td>n<td>Moment order (optional, default 1)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Mn<td>n-th moment of the maximum
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Eq. (27) and
 % the inclusion-exclusion expansion of Section 7.3 on page 17:47.
%}
function Mn = fj_xmax_het(lambda, n)

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
if K > 24
    line_error(mfilename, 'Inclusion-exclusion over K=%d rates needs 2^K terms; use fj_xmax_moments_het instead.', K);
end

nfact = factorial(n);
Mn = 0;
% Enumerate the 2^K - 1 nonempty subsets through their bit patterns
for mask = 1:(2^K - 1)
    bits = bitget(mask, 1:K);
    rate = sum(lambda(bits == 1));
    card = sum(bits);
    Mn = Mn + (-1)^(card + 1) * nfact / rate^n;
end

end
