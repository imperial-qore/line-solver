%{ @file fj_rmax_erlang.m
 %  @brief Maximum response time R_K^max for Erlang service times
 %
 %  @author LINE Development Team
%}

%{
 % @brief Maximum response time R_K^max for Erlang service times
 %
 % @details
 % Computes the expected value of the maximum response time for K
 % M/E_k/1 queues (Poisson arrivals, Erlang-k service times).
 %
 % For K=2 (two servers), uses the closed-form formula:
 %   R_2^max = R_1 + R_2 - sum_{m=0}^{k1-1} sum_{n=0}^{k2-1} C(m+n,m) * mu1^m * mu2^n / (mu1+mu2)^{m+n+1}
 %
 % For general K, uses numerical integration:
 %   R_K^max = integral_0^inf [1 - prod_{i=1}^{K} F_Erlang(t)] dt
 %
 % where F_Erlang(t) = 1 - exp(-mu*t) * sum_{j=0}^{k-1} (mu*t)^j / j!
 %
 % @par Syntax:
 % @code
 % Rmax = fj_rmax_erlang(K, k, lambda, mu)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of parallel servers (positive integer)
 % <tr><td>k<td>Number of Erlang stages (positive integer)
 % <tr><td>lambda<td>Arrival rate
 % <tr><td>mu<td>Service rate per Erlang stage (mean service = k/mu)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Rmax<td>Expected maximum response time R_K^max
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Eq. (33) and (34) on page 17:18.
%}
function Rmax = fj_rmax_erlang(K, k, lambda, mu)

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

if k < 1
    line_error(mfilename, 'k (number of stages) must be a positive integer. Got k=%d.', k);
end

% Compute mean service time and utilization
mean_service = k / mu;
rho = lambda * mean_service;

if rho >= 1
    line_error(mfilename, 'System is unstable: rho = %.4f >= 1.', rho);
end

% Compute response time for M/E_k/1 queue using Pollaczek-Khinchin formula
% R = mean_service / (1 - rho) * [1 + rho * (1 + 1/k) / 2]
% Simplified: R = mean_service * (1 + rho*(k+1)/(2*k*(1-rho)))
cv2 = 1 / k;  % Squared CV of Erlang-k
R_single = mean_service * (1 + rho * (1 + cv2) / (2 * (1 - rho)));

if K == 2
    % Use closed-form formula (Eq. 34)
    % For identical servers: k1 = k2 = k, mu1 = mu2 = k/R
    mu_resp = k / R_single;  % Effective rate for response time

    % R_2^max = 2*R - sum_{m=0}^{k-1} sum_{n=0}^{k-1} C(m+n,m) * mu^m * mu^n / (2*mu)^{m+n+1}
    correction = 0;
    for m = 0:(k-1)
        for n = 0:(k-1)
            correction = correction + nchoosek(m + n, m) * (mu_resp^m) * (mu_resp^n) / ((2 * mu_resp)^(m + n + 1));
        end
    end
    Rmax = 2 * R_single - correction;
else
    % Use numerical integration for general K
    % The response time distribution for M/E_k/1 can be approximated by Erlang
    % with parameters matched to the response time moments

    % For M/E_k/1, use the approximation that response time is approximately
    % Erlang with adjusted parameters
    k_resp = ceil(1 / cv2_response(k, rho));  % Effective stages for response time
    mu_resp = k_resp / R_single;

    % Erlang CDF for response time
    erlang_cdf = @(t) erlang_cdf_func(t, k_resp, mu_resp);

    % Expected maximum: integral_0^inf [1 - F(t)^K] dt
    integrand = @(t) 1 - erlang_cdf(t).^K;

    % Numerical integration
    upper_limit = R_single * 20;  % Well beyond the tail
    Rmax = integral(integrand, 0, upper_limit, 'RelTol', 1e-10, 'AbsTol', 1e-12);
end

end

function cv2 = cv2_response(k, rho)
% Squared CV of response time for M/E_k/1 queue (approximation)
% Based on: c2_R ≈ (c2_X + rho*(1+c2_X)/2) / (1 + rho*(1+c2_X)/2)
cv2_X = 1 / k;
cv2 = (cv2_X + rho) / (1 + rho);  % Simplified approximation
cv2 = max(cv2, 1/20);  % Bound to avoid numerical issues
end

function F = erlang_cdf_func(t, k, mu)
% Erlang-k CDF: F(t) = 1 - exp(-mu*t) * sum_{j=0}^{k-1} (mu*t)^j / j!
F = zeros(size(t));
for i = 1:numel(t)
    if t(i) <= 0
        F(i) = 0;
    else
        S = 0;
        for j = 0:(k-1)
            S = S + (mu * t(i))^j / factorial(j);
        end
        F(i) = 1 - exp(-mu * t(i)) * S;
    end
end
end
