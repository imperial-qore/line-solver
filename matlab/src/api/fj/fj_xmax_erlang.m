%{ @file fj_xmax_erlang.m
 %  @brief Expected maximum of K i.i.d. Erlang-k random variables
 %
 %  @author LINE Development Team
%}

%{
 % @brief Expected maximum of K i.i.d. Erlang-k random variables
 %
 % @details
 % Computes the expected value of the maximum of K independent and
 % identically distributed Erlang random variables with k stages
 % and rate mu per stage.
 %
 % For k=2 (Erlang-2):
 %   X_K^max = (1/mu) * sum_{n=1}^{K} C(K,n) * (-1)^{n-1} * sum_{m=1}^{n} C(n,m) * m! / (2*n^{m+1})
 %
 % For general k-stage Erlang, the formula involves k-fold convolution
 % of exponential distributions.
 %
 % The Erlang distribution has mean k/mu and variance k/mu^2, giving
 % CV = 1/sqrt(k) < 1.
 %
 % @par Syntax:
 % @code
 % Xmax = fj_xmax_erlang(K, k, mu)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of parallel servers (positive integer)
 % <tr><td>k<td>Number of Erlang stages (positive integer, default=2)
 % <tr><td>mu<td>Rate parameter per stage (mean service time = k/mu)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Xmax<td>Expected maximum service time X_K^max
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Page 17:13.
%}
function Xmax = fj_xmax_erlang(K, k, mu)

if nargin < 2
    k = 2;  % Default to Erlang-2
end

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

if k < 1
    line_error(mfilename, 'k (number of stages) must be a positive integer. Got k=%d.', k);
end

if mu <= 0
    line_error(mfilename, 'Rate mu must be positive. Got mu=%.4f.', mu);
end

% For Erlang-2 (k=2), use the closed-form formula from the paper
if k == 2
    % X_K^max = (1/mu) * sum_{n=1}^{K} C(K,n) * (-1)^{n-1} * sum_{m=1}^{n} C(n,m) * m! / (2*n^{m+1})
    outer_sum = 0;
    for n = 1:K
        inner_sum = 0;
        for m = 1:n
            inner_sum = inner_sum + nchoosek(n, m) * factorial(m) / (2 * n^(m + 1));
        end
        outer_sum = outer_sum + nchoosek(K, n) * ((-1)^(n - 1)) * inner_sum;
    end
    Xmax = outer_sum / mu;
else
    % For general k-stage Erlang, use numerical integration
    % The CDF of Erlang-k is: F(t) = 1 - exp(-mu*t) * sum_{j=0}^{k-1} (mu*t)^j / j!
    % The expected max is: X_K^max = integral_0^inf [1 - F(t)^K] dt

    % Erlang CDF function
    erlang_cdf = @(t) 1 - exp(-mu * t) .* sum_erlang_terms(mu * t, k);

    % Integrand: 1 - F(t)^K
    integrand = @(t) 1 - erlang_cdf(t).^K;

    % Numerical integration (upper limit chosen based on Erlang mean and variance)
    mean_erlang = k / mu;
    upper_limit = mean_erlang * 10 + 10 * sqrt(k) / mu;  % Well beyond the tail
    Xmax = integral(integrand, 0, upper_limit, 'RelTol', 1e-10, 'AbsTol', 1e-12);
end

end

function S = sum_erlang_terms(x, k)
% Helper function to compute sum_{j=0}^{k-1} x^j / j!
% Handles both scalar and vector inputs
S = zeros(size(x));
for j = 0:(k-1)
    S = S + (x.^j) / factorial(j);
end
end
