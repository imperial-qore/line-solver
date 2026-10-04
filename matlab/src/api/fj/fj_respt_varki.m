%{ @file fj_respt_varki.m
 %  @brief Varki et al. approximation for K-way F/J response time
 %
 %  @author LINE Development Team
%}

%{
 % @brief Varki et al. approximation for K-way F/J response time
 %
 % @details
 % Computes an approximate mean response time for a K-way Fork-Join
 % queueing system with Poisson arrivals and exponential service times.
 % This is the mean of pessimistic (upper) and optimistic (lower) bounds.
 %
 % R_K^{F/J}(rho) ≈ (1/mu) * [H_K + (rho/(2*(1-rho))) * (S1 + (1-2*rho)*S2)]
 %
 % where:
 %   S1 = sum_{i=1}^{K} 1/(i - rho)
 %   S2 = sum_{i=1}^{K} 1/(i*(i - rho))
 %
 % @par Syntax:
 % @code
 % R = fj_respt_varki(K, lambda, mu)
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
 % Table I, Eq. (5) on page 17:10.
 %
 % Original: E. Varki, A. Merchant, J. Xu, and X. Qiu, "Issues and Challenges
 % in the Performance Analysis of Real Disk Arrays", IEEE Trans. Parallel
 % Distrib. Syst., 25(6), 2014 (published as Varki et al. [2013] in draft).
%}
function R = fj_respt_varki(K, lambda, mu)

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

% Compute S1 = sum_{i=1}^{K} 1/(i - rho)
S1 = sum(1 ./ ((1:K) - rho));

% Compute S2 = sum_{i=1}^{K} 1/(i*(i - rho))
S2 = sum(1 ./ ((1:K) .* ((1:K) - rho)));

% Varki et al. approximation (Eq. 5):
% R_K^{F/J}(rho) ≈ (1/mu) * [H_K + (rho/(2*(1-rho))) * (S1 + (1-2*rho)*S2)]
R = (1 / mu) * (H_K + (rho / (2 * (1 - rho))) * (S1 + (1 - 2 * rho) * S2));

end
