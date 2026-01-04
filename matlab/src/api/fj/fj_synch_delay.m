%{ @file fj_synch_delay.m
 %  @brief Synchronization delay S_2(rho) for two-way Fork-Join
 %
 %  @author LINE Development Team
%}

%{
 % @brief Synchronization delay S_2(rho) for two-way Fork-Join
 %
 % @details
 % Computes the mean synchronization delay for a 2-way Fork-Join
 % queueing system. This is the average time tasks spend waiting
 % at synchronization queues for their siblings to complete.
 %
 % S_2(rho) = (1/2) * (1 - rho/4) * R(rho)
 %
 % The F/J response time decomposes as:
 %   R_2^{F/J}(rho) = R(rho) + S_2(rho)
 %
 % where R(rho) is the delay at the server (M/M/1 response time).
 %
 % @par Syntax:
 % @code
 % S = fj_synch_delay(lambda, mu)
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
 % <tr><td>S<td>Mean synchronization delay S_2(rho)
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Eq. (8) on page 17:11.
%}
function S = fj_synch_delay(lambda, mu)

% Compute utilization
rho = lambda / mu;

if rho >= 1
    line_error(mfilename, 'System is unstable: rho = lambda/mu = %.4f >= 1. Require lambda < mu.', rho);
end

% Compute M/M/1 mean response time
R_rho = qsys_mm1(lambda, mu);

% Synchronization delay (Eq. 8): S_2(rho) = (1/2) * (1 - rho/4) * R(rho)
S = 0.5 * (1 - rho / 4) * R_rho;

end
