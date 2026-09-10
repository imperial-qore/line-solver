%{ @file fj_xmax_normal.m
 %  @brief Expected maximum for normal distribution
 %
 %  @author LINE Development Team
%}

%{
 % @brief Expected maximum for normal distribution
 %
 % @details
 % Computes the expected value of the maximum of K i.i.d. normal random
 % variables using the approximation from Johnson et al. [1995].
 %
 % E[Y_K] ≈ μ + σ * [sqrt(2*ln(K)) - (ln(ln(K)) - ln(4π) + 2γ) / (2*sqrt(2*ln(K)))]
 %
 % where γ ≈ 0.5772 is the Euler-Mascheroni constant.
 %
 % A bias correction from Petzold [2000] can be applied:
 %   δ(K) = 0.1727 * K^(-0.2750)
 %
 % The variance of the maximum is approximately:
 %   Var[Y_K] ≈ 1.64492 * σ^2 / (2*ln(K))
 %
 % A simpler approximation from Arnold [1980] is also available:
 %   E[Y_K] ≈ μ + σ * sqrt(2*ln(K))
 %
 % @par Syntax:
 % @code
 % Xmax = fj_xmax_normal(K, mu, sigma)
 % Xmax = fj_xmax_normal(K, mu, sigma, 'johnson')     % Johnson et al. (default)
 % Xmax = fj_xmax_normal(K, mu, sigma, 'arnold')      % Arnold approximation
 % Xmax = fj_xmax_normal(K, mu, sigma, 'corrected')   % With Petzold correction
 % [Xmax, Vmax] = fj_xmax_normal(...)                 % Also return variance
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of random variables (positive integer >= 2)
 % <tr><td>mu<td>Mean of the normal distribution
 % <tr><td>sigma<td>Standard deviation of the normal distribution
 % <tr><td>method<td>Optional: 'johnson' (default), 'arnold', or 'corrected'
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Xmax<td>Approximate expected maximum E[Y_K]
 % <tr><td>Vmax<td>Approximate variance of maximum Var[Y_K]
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Eq. (45) on page 17:23.
 %
 % N.L. Johnson, S. Kotz, N. Balakrishnan, "Continuous Univariate
 % Distributions", Vol. 2, Wiley, 1995.
%}
function [Xmax, Vmax] = fj_xmax_normal(K, mu, sigma, method)

if nargin < 4
    method = 'johnson';
end

if K < 2
    line_error(mfilename, 'K must be at least 2 for normal approximation. Got K=%d.', K);
end

if sigma < 0
    line_error(mfilename, 'Standard deviation sigma must be non-negative.');
end

% Euler-Mascheroni constant
gamma_em = 0.5772156649015329;

% Common term: sqrt(2*ln(K))
sqrt_2lnK = sqrt(2 * log(K));

switch lower(method)
    case 'arnold'
        % Simple approximation: E[Y_K] ≈ μ + σ * sqrt(2*ln(K))
        GK = sqrt_2lnK;
        Xmax = mu + sigma * GK;

    case 'johnson'
        % Johnson et al. approximation (Eq. 45)
        % E[Y_K] ≈ μ + σ * [sqrt(2*ln(K)) - (ln(ln(K)) - ln(4π) + 2γ) / (2*sqrt(2*ln(K)))]
        correction = (log(log(K)) - log(4 * pi) + 2 * gamma_em) / (2 * sqrt_2lnK);
        GK = sqrt_2lnK - correction;
        Xmax = mu + sigma * GK;

    case 'corrected'
        % Johnson with Petzold bias correction
        correction = (log(log(K)) - log(4 * pi) + 2 * gamma_em) / (2 * sqrt_2lnK);
        GK = sqrt_2lnK - correction;
        % Subtract bias: δ(K) = 0.1727 * K^(-0.2750)
        delta_K = 0.1727 * K^(-0.2750);
        GK = GK - delta_K;
        Xmax = mu + sigma * GK;

    otherwise
        line_error(mfilename, 'Unknown method: %s. Valid: johnson, arnold, corrected.', method);
end

% Compute variance if requested
if nargout > 1
    % Var[Y_K] ≈ 1.64492 * σ^2 / (2*ln(K))
    Vmax = 1.64492 * sigma^2 / (2 * log(K));
end

end
