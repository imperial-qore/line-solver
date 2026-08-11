%{ @file fj_xmax_approx.m
 %  @brief General approximation for expected maximum of K random variables
 %
 %  @author LINE Development Team
%}

%{
 % @brief General approximation for expected maximum of K random variables
 %
 % @details
 % Approximates the expected value of the maximum of K i.i.d. random
 % variables using the mean-variance approximation from David [1970].
 %
 % X_K^max ≈ mu_X + sigma_X * G(K)
 %
 % where G(K) depends on the distribution type:
 %   - Exponential: G(K) = H_K - 1 (where H_K is the K-th harmonic number)
 %   - Uniform:     G(K) = sqrt(3) * (K-1) / (K+1)
 %   - EVD:         G(K) = sqrt(6) * ln(K) / pi
 %   - General:     G(K) <= (K-1) / sqrt(2K-1)  (upper bound)
 %
 % @par Syntax:
 % @code
 % Xmax = fj_xmax_approx(K, mu_X, sigma_X, dist_type)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of random variables (positive integer)
 % <tr><td>mu_X<td>Mean of the distribution
 % <tr><td>sigma_X<td>Standard deviation of the distribution
 % <tr><td>dist_type<td>Distribution type: 'exp', 'uniform', 'evd', 'bound' (default: 'exp')
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Xmax<td>Approximate expected maximum X_K^max
 % <tr><td>GK<td>The G(K) factor used
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Eq. (22) and (23) on page 17:16.
 %
 % Original: H.A. David, "Order Statistics", Wiley, 1970.
%}
function [Xmax, GK] = fj_xmax_approx(K, mu_X, sigma_X, dist_type)

if nargin < 4
    dist_type = 'exp';
end

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

if sigma_X < 0
    line_error(mfilename, 'Standard deviation sigma_X must be non-negative.');
end

% Compute G(K) based on distribution type
switch lower(dist_type)
    case 'exp'
        % Exponential distribution: G(K) = H_K - 1
        GK = fj_harmonic(K) - 1;

    case 'uniform'
        % Uniform distribution: G(K) = sqrt(3) * (K-1) / (K+1)
        GK = sqrt(3) * (K - 1) / (K + 1);

    case 'evd'
        % Extreme Value Distribution: G(K) = sqrt(6) * ln(K) / pi
        GK = sqrt(6) * log(K) / pi;

    case 'bound'
        % Upper bound for any distribution: G(K) <= (K-1) / sqrt(2K-1)
        GK = (K - 1) / sqrt(2 * K - 1);

    otherwise
        line_error(mfilename, 'Unknown distribution type: %s. Valid: exp, uniform, evd, bound.', dist_type);
end

% X_K^max ≈ mu_X + sigma_X * G(K)
Xmax = mu_X + sigma_X * GK;

end
