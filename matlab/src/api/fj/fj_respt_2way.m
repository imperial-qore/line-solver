%{ @file fj_respt_2way.m
 %  @brief Exact two-way Fork-Join response time R_2^{F/J}(rho)
 %
 %  @author LINE Development Team
%}

%{
 % @brief Exact two-way Fork-Join response time R_2^{F/J}(rho)
 %
 % @details
 % Computes the exact mean response time for a 2-way (K=2) Fork-Join
 % queueing system with Poisson arrivals and exponential service times.
 % This is derived from the analysis in Flatto and Hahn [1984].
 %
 % R_2^{F/J}(rho) = (H_2 - rho/8) * R(rho) = (12 - rho)/8 * R(rho)
 %
 % where H_2 = 1.5 and R(rho) = (mu - lambda)^{-1} is the M/M/1 mean
 % response time.
 %
 % Note that R_2^{F/J}(rho) < R_2^max(rho), with the difference being
 % rho/8 * R(rho).
 %
 % @par Syntax:
 % @code
 % R = fj_respt_2way(lambda, mu)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>Arrival rate
 % <tr><td>mu<td>Service rate (mu > lambda for stability)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>R<td>Exact 2-way F/J response time R_2^{F/J}(rho)
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Eq. (6) on page 17:10.
 %
 % Original derivation: L. Flatto and S. Hahn, "Two Parallel Queues Created
 % by Arrivals with Two Demands I", SIAM J. Appl. Math., 44(5), 1984.
%}
function R = fj_respt_2way(lambda, mu)

% Compute utilization
rho = lambda / mu;

if rho >= 1
    line_error(mfilename, 'System is unstable: rho = lambda/mu = %.4f >= 1. Require lambda < mu.', rho);
end

% Compute M/M/1 mean response time
R_rho = qsys_mm1(lambda, mu);

% Exact 2-way F/J response time: R_2^{F/J}(rho) = (H_2 - rho/8) * R(rho)
% H_2 = 1 + 1/2 = 1.5
% (H_2 - rho/8) = (1.5 - rho/8) = (12 - rho)/8
R = (12 - rho) / 8 * R_rho;

end
