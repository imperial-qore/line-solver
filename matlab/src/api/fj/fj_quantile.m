%{ @file fj_quantile.m
 %  @brief Quantile approximation for maximum of K random variables
 %
 %  @author LINE Development Team
%}

%{
 % @brief Quantile approximation for maximum of K random variables
 %
 % @details
 % Approximates the q-th quantile of the maximum of K i.i.d. random
 % variables from the standard Gumbel (Type I extreme value) distribution.
 %
 % x(K,q) ≈ ln(K) - ln(ln(1/q))
 %
 % where P(M_K <= x(K,q)) = q is the q-th quantile.
 %
 % Note: This approximation is inaccurate for small values of K.
 %
 % For general distributions, the quantile of the maximum can be computed
 % using the inverse CDF:
 %   x(K,q) = F^{-1}(q^{1/K})
 %
 % @par Syntax:
 % @code
 % x_Kq = fj_quantile(K, q)                  % Standard Gumbel approx
 % x_Kq = fj_quantile(K, q, F_inv)           % General distribution
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of random variables (positive integer)
 % <tr><td>q<td>Quantile probability (0 < q < 1)
 % <tr><td>F_inv<td>Optional: inverse CDF function handle for general distributions
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>x_Kq<td>q-th quantile of the maximum of K random variables
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Page 17:22.
%}
function x_Kq = fj_quantile(K, q, F_inv)

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

if any(q <= 0) || any(q >= 1)
    line_error(mfilename, 'Quantile q must satisfy 0 < q < 1.');
end

if nargin < 3 || isempty(F_inv)
    % Standard Gumbel approximation
    % x(K,q) ≈ ln(K) - ln(ln(1/q))
    x_Kq = log(K) - log(log(1 ./ q));
else
    % General distribution using inverse CDF
    % The CDF of the maximum is F_max(x) = [F(x)]^K
    % So the q-th quantile satisfies: [F(x)]^K = q
    % Therefore: F(x) = q^(1/K) and x = F^{-1}(q^{1/K})
    x_Kq = F_inv(q.^(1 / K));
end

end
