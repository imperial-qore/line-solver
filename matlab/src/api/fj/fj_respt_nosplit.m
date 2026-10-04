%{ @file fj_respt_nosplit.m
 %  @brief Mean response time of the distributed no-splitting parallel system
 %
 %  @author LINE Development Team
%}

%{
 % @brief Mean response time of the distributed no-splitting parallel system
 %
 % @details
 % In the distributed no-splitting policy a job is a sequence of K tasks that
 % is routed in one piece to a single server chosen uniformly among the K
 % servers, so each server sees a Poisson stream of rate lambda/K whose service
 % time is the sum of K exponential stages of rate mu. That is an M/E_K/1
 % queue, and the Pollaczek-Khinchine mean response time reduces to
 %
 %   R(rho) = 1/(mu - lambda),   rho = lambda/mu,
 %   R_{D/NS} = [ K - (K-1)*rho/2 ] * R(rho),
 %
 % which is the reference against which the splitting policies of the same
 % section are judged: splitting wins because it leaves fewer servers idle.
 % At K = 1 the expression collapses to the M/M/1 response time.
 %
 % @par Syntax:
 % @code
 % R = fj_respt_nosplit(K, lambda, mu)
 % [R, rho] = fj_respt_nosplit(K, lambda, mu)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of servers, equal to the number of tasks per job
 % <tr><td>lambda<td>Total job arrival rate
 % <tr><td>mu<td>Task service rate at each server, mu > lambda for stability
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>R<td>Mean job response time under distributed no splitting
 % <tr><td>rho<td>Utilization of each server, lambda/mu
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Section 6.1
 % on page 17:33.
 %
 % Original: R. Nelson, D. Towsley, A. N. Tantawi, "Performance Analysis of
 % Parallel Processing Systems", IEEE Trans. Software Eng. 14(4), 1988.
%}
function [R, rho] = fj_respt_nosplit(K, lambda, mu)

if K < 1 || K ~= round(K)
    line_error(mfilename, 'K must be a positive integer. Got K=%g.', K);
end
if lambda <= 0
    line_error(mfilename, 'The arrival rate lambda must be positive. Got lambda=%g.', lambda);
end
if mu <= 0
    line_error(mfilename, 'The service rate mu must be positive. Got mu=%g.', mu);
end

rho = lambda / mu;
if rho >= 1
    line_error(mfilename, 'System is unstable: rho = lambda/mu = %.4f >= 1. Require lambda < mu.', rho);
end

% M/M/1 response time at the same utilization
Rmm1 = 1 / (mu - lambda);

R = (K - (K - 1) * rho / 2) * Rmm1;

end
