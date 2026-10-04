%{ @file fj_xmax_emma.m
 %  @brief Expected maximum using EMMA method
 %
 %  @author LINE Development Team
%}

%{
 % @brief Expected maximum using EMMA method
 %
 % @details
 % Computes the expected value of the maximum of K i.i.d. random variables
 % using the EMMA (Expected Maximum from Marginal Approximation) method.
 %
 % The method is based on the property that:
 %   [F_X(E[Y_K])]^K ≈ φ = 0.570376
 %
 % Therefore:
 %   E[Y_K] ≈ F^{-1}(φ^{1/K})
 %
 % For the exponential distribution with rate μ:
 %   E[Y_K] = -(1/μ) * ln(1 - φ^{1/K})
 %
 % The constant φ = exp(-exp(-γ)) ≈ 0.570376, where γ ≈ 0.5772 is the
 % Euler-Mascheroni constant.
 %
 % @par Syntax:
 % @code
 % Xmax = fj_xmax_emma(K, mu)                    % Exponential with rate mu
 % Xmax = fj_xmax_emma(K, F_inv)                 % General CDF inverse
 % Xmax = fj_xmax_emma(K, mu, 'exp')             % Explicit exponential
 % Xmax = fj_xmax_emma(K, F_inv, 'general')      % Explicit general
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of random variables (positive integer)
 % <tr><td>param<td>Rate mu for exponential, or inverse CDF function handle
 % <tr><td>dist_type<td>Optional: 'exp' (default) or 'general'
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Xmax<td>Approximate expected maximum E[Y_K]
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Eq. (44) on page 17:22.
 %
 % Original: Y. Sun, K.L. Peterson, "Computing Extreme Order Statistics from
 % Large Data Sets", 2012.
%}
function Xmax = fj_xmax_emma(K, param, dist_type)

% EMMA constant: φ = exp(-exp(-γ)) where γ is Euler-Mascheroni constant
phi = 0.570376;

if nargin < 3
    % Determine type from param
    if isa(param, 'function_handle')
        dist_type = 'general';
    else
        dist_type = 'exp';
    end
end

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

switch lower(dist_type)
    case 'exp'
        % Exponential distribution with rate mu
        mu = param;
        if mu <= 0
            line_error(mfilename, 'Rate mu must be positive.');
        end
        % E[Y_K] = -(1/mu) * ln(1 - phi^(1/K))
        Xmax = -(1 / mu) * log(1 - phi^(1 / K));

    case 'general'
        % General distribution: param is the inverse CDF (quantile function)
        F_inv = param;
        if ~isa(F_inv, 'function_handle')
            line_error(mfilename, 'For general distribution, param must be a function handle for the inverse CDF.');
        end
        % E[Y_K] = F^{-1}(phi^{1/K})
        Xmax = F_inv(phi^(1 / K));

    otherwise
        line_error(mfilename, 'Unknown distribution type: %s. Valid: exp, general.', dist_type);
end

end
