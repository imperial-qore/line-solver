%{ @file fj_bounds.m
 %  @brief Upper and lower bounds for K-way F/J response time
 %
 %  @author LINE Development Team
%}

%{
 % @brief Upper and lower bounds for K-way F/J response time
 %
 % @details
 % Computes pessimistic (upper) and optimistic (lower) bounds for the
 % mean response time of a K-way Fork-Join queueing system.
 %
 % Upper bound (pessimistic):
 %   R_K^max(rho) = H_K / (mu * (1 - rho)) = H_K * R(rho)
 %
 % Lower bound (optimistic) from Varki et al.:
 %   R_K^{F/J(opt)}(rho) = (1/mu) * [H_K + S_{K(K-rho)}]
 %   where S_{K(K-rho)} = sum_{j=1}^{K} (1/j) * (rho/(j - rho))
 %
 % The bounds satisfy: R_K^{F/J(opt)}(rho) <= R_K^{F/J}(rho) <= R_K^max(rho)
 %
 % @par Syntax:
 % @code
 % [Rmax, Rmin] = fj_bounds(K, lambda, mu)
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
 % <tr><td>Rmax<td>Upper bound (pessimistic): R_K^max(rho)
 % <tr><td>Rmin<td>Lower bound (optimistic): R_K^{F/J(opt)}(rho)
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Eq. (1) and Eq. (2) on page 17:9.
%}
function [Rmax, Rmin] = fj_bounds(K, lambda, mu)

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

% Upper bound (Eq. 1): R_K^max(rho) = H_K / (mu * (1 - rho))
Rmax = H_K / (mu * (1 - rho));

% Lower bound (Eq. 2): R_K^{F/J(opt)}(rho) = (1/mu) * [H_K + S_{K(K-rho)}]
% where S_{K(K-rho)} = sum_{j=1}^{K} (1/j) * (rho/(j - rho))
S_K = sum((1 ./ (1:K)) .* (rho ./ ((1:K) - rho)));
Rmin = (1 / mu) * (H_K + S_K);

end
