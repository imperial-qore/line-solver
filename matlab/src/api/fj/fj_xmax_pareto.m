%{ @file fj_xmax_pareto.m
 %  @brief Expected maximum for Pareto distribution
 %
 %  @author LINE Development Team
%}

%{
 % @brief Expected maximum for Pareto distribution
 %
 % @details
 % Computes the expected value of the maximum of K i.i.d. Pareto random
 % variables using order statistics formulas.
 %
 % Pareto distribution CDF:
 %   F(x) = 1 - α(x + γ)^{-β}, x >= 0, β > 2
 %
 % For the standardized Pareto with mean 1:
 %   α = γ^β, γ = β - 1
 %   Mean = 1, Second moment M = 2 + 2/(β - 2)
 %
 % The j-th moment of the i-th order statistic is:
 %   m_{i,j} = Γ(K+1) * Γ(K - j + 1 - i/α) / [Γ(K + 1 - i/α) * Γ(K - j + 1)] * k^i
 %
 % For the maximum (i = K), the expected value can be computed numerically
 % or approximated using the characteristic maximum method.
 %
 % @par Syntax:
 % @code
 % Xmax = fj_xmax_pareto(K, beta)               % Standard Pareto (mean=1)
 % Xmax = fj_xmax_pareto(K, beta, k)            % Pareto with scale k
 % [Xmax, MK] = fj_xmax_pareto(K, beta, k)      % Also return characteristic max
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of random variables (positive integer)
 % <tr><td>beta<td>Shape parameter (beta > 2 required for finite moments)
 % <tr><td>k<td>Optional: scale parameter (default: beta - 1 for mean = 1)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Xmax<td>Expected maximum E[Y_K]
 % <tr><td>MK<td>Characteristic maximum M_K
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Page 17:24 and 17:26.
 %
 % N.L. Johnson, S. Kotz, N. Balakrishnan, "Continuous Univariate
 % Distributions", Vol. 1, Wiley, 1995.
%}
function [Xmax, MK] = fj_xmax_pareto(K, beta, k)

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

if beta <= 2
    line_error(mfilename, 'Shape parameter beta must be > 2 for finite moments. Got beta=%.4f.', beta);
end

% Default scale parameter for mean = 1
if nargin < 3 || isempty(k)
    k = beta - 1;
end

if k <= 0
    line_error(mfilename, 'Scale parameter k must be positive.');
end

% Pareto survival function: S(x) = (k / (k + x))^beta for x >= 0
S_pareto = @(x) (k ./ (k + max(x, 0))).^beta;

% Find m_K such that S(m_K) = 1/K
% (k / (k + m_K))^beta = 1/K
% k + m_K = k * K^(1/beta)
% m_K = k * (K^(1/beta) - 1)
mK = k * (K^(1/beta) - 1);

% Compute expected maximum using numerical integration
% E[Y_K] = integral_0^inf [1 - F(x)^K] dx = integral_0^inf [1 - (1-S(x))^K] dx
% This can be split into:
% E[Y_K] = integral_0^inf P(max > x) dx = integral_0^inf [1 - (1-S(x))^K] dx

% For Pareto, we can use the order statistic formula
% E[Y_(K)] for the K-th order statistic (maximum)
Xmax = compute_pareto_max_expectation(K, beta, k);

% Characteristic maximum
if nargout > 1
    % M_K = m_K + K * integral_{m_K}^{inf} S(x) dx
    % For Pareto: integral_{m_K}^{inf} (k/(k+x))^beta dx
    %           = k^beta * integral_{m_K}^{inf} (k+x)^{-beta} dx
    %           = k^beta * [(k+x)^{1-beta} / (1-beta)]_{m_K}^{inf}
    %           = k^beta * (k+m_K)^{1-beta} / (beta-1)

    if beta > 1
        tail_integral = k^beta * (k + mK)^(1 - beta) / (beta - 1);
        MK = mK + K * tail_integral;
    else
        % For beta <= 1, the integral diverges
        MK = Inf;
    end
end

end

function Xmax = compute_pareto_max_expectation(K, beta, k)
% Compute E[max(X_1, ..., X_K)] for Pareto distribution
% Using numerical integration: E[Y_K] = integral_0^inf [1 - F(x)^K] dx

% CDF: F(x) = 1 - (k/(k+x))^beta
F_pareto = @(x) 1 - (k ./ (k + max(x, 0))).^beta;

% Integrand: 1 - F(x)^K = 1 - [1 - (k/(k+x))^beta]^K
integrand = @(x) 1 - F_pareto(x).^K;

% Need to integrate from 0 to infinity
% For numerical stability, split the integral and use appropriate bounds
% The tail decays as (k/(k+x))^beta ~ k^beta / x^beta for large x

% Find a reasonable upper limit where integrand is negligible
upper_limit = k * K^(2/beta) * 10;  % Conservative estimate

Xmax = integral(integrand, 0, upper_limit, 'RelTol', 1e-8, 'AbsTol', 1e-10);
end
