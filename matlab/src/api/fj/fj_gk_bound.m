%{ @file fj_gk_bound.m
 %  @brief Compute G(K) factors for expected maximum approximation
 %
 %  @author LINE Development Team
%}

%{
 % @brief Compute G(K) factors for expected maximum approximation
 %
 % @details
 % Computes the G(K) factor used in the approximation:
 %   X_K^max ≈ mu_X + sigma_X * G(K)
 %
 % For standardized distributions (mean=0, variance=1), G(K) = X_K^max.
 %
 % Available G(K) formulas:
 %   - Exponential: G(K) = H_K - 1
 %   - Uniform:     G(K) = sqrt(3) * (K-1) / (K+1)
 %   - EVD:         G(K) = sqrt(6) * ln(K) / pi
 %   - Upper bound: G(K) <= (K-1) / sqrt(2K-1)  (tight bound from David [1970])
 %
 % @par Syntax:
 % @code
 % GK = fj_gk_bound(K)                    % Returns struct with all
 % GK = fj_gk_bound(K, 'exp')             % Exponential G(K)
 % GK = fj_gk_bound(K, 'uniform')         % Uniform G(K)
 % GK = fj_gk_bound(K, 'evd')             % EVD G(K)
 % GK = fj_gk_bound(K, 'bound')           % Upper bound
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of random variables (positive integer)
 % <tr><td>type<td>Optional: 'exp', 'uniform', 'evd', 'bound', or 'all' (default)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>GK<td>G(K) value(s) - scalar or struct depending on type
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Eq. (22) and (23) on page 17:16.
 %
 % Original: H.A. David, "Order Statistics", Wiley, 1970.
%}
function GK = fj_gk_bound(K, type)

if nargin < 2
    type = 'all';
end

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

switch lower(type)
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

    case 'all'
        % Return struct with all G(K) values
        GK = struct();
        GK.K = K;
        GK.exponential = fj_harmonic(K) - 1;
        GK.uniform = sqrt(3) * (K - 1) / (K + 1);
        GK.evd = sqrt(6) * log(K) / pi;
        GK.upper_bound = (K - 1) / sqrt(2 * K - 1);

    otherwise
        line_error(mfilename, 'Unknown type: %s. Valid: exp, uniform, evd, bound, all.', type);
end

end
