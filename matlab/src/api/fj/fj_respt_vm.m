%{ @file fj_respt_vm.m
 %  @brief Varma-Makowski approximation for K-way F/J response time
 %
 %  @author LINE Development Team
%}

%{
 % @brief Varma-Makowski approximation for K-way F/J response time
 %
 % @details
 % Computes an approximate mean response time for a K-way Fork-Join
 % queueing system with Poisson arrivals and exponential service times.
 %
 % R_K^{F/J}(rho) ≈ [H_K + (A_K - H_K) * rho] * (mu - lambda)^{-1}
 %
 % where A_K = sum_{i=1}^{K} C(K,i) * (-1)^{i-1} * sum_{m=1}^{i} C(i,m) * (m-1)! / i^{m+1}
 %
 % and C(n,k) denotes the binomial coefficient "n choose k".
 %
 % @par Syntax:
 % @code
 % R = fj_respt_vm(K, lambda, mu)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of parallel servers (positive integer)
 % <tr><td>lambda<td>Arrival rate
 % <tr><td>mu<td>Service rate (mu > lambda for stability)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>R<td>Approximate K-way F/J response time R_K^{F/J}(rho)
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Table I, Eq. (4) on page 17:10.
 %
 % Original: S. Varma and A.M. Makowski, "Interpolation Approximations for
 % Symmetric Fork-Join Queues", Performance Evaluation, 20, 1994.
%}
function R = fj_respt_vm(K, lambda, mu)

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

% Compute utilization
rho = lambda / mu;

if rho >= 1
    line_error(mfilename, 'System is unstable: rho = lambda/mu = %.4f >= 1. Require lambda < mu.', rho);
end

% Compute harmonic number H_K
H_K = fj_harmonic(K);

% Compute A_K = sum_{i=1}^{K} C(K,i) * (-1)^{i-1} * sum_{m=1}^{i} C(i,m) * (m-1)! / i^{m+1}
A_K = 0;
for i = 1:K
    inner_sum = 0;
    for m = 1:i
        inner_sum = inner_sum + nchoosek(i, m) * factorial(m - 1) / (i^(m + 1));
    end
    A_K = A_K + nchoosek(K, i) * ((-1)^(i - 1)) * inner_sum;
end

% Varma-Makowski approximation (Eq. 4):
% R_K^{F/J}(rho) ≈ [H_K + (A_K - H_K) * rho] * (mu - lambda)^{-1}
R_rho = 1 / (mu - lambda);
R = (H_K + (A_K - H_K) * rho) * R_rho;

end
