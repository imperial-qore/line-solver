%{ @file fj_mg1_respt_moments.m
 %  @brief Mean and variance of the M/G/1 response time, as ForkTail inputs
 %
 %  @author LINE Development Team
%}

%{
 % @brief Mean and variance of the M/G/1 response time, as ForkTail inputs
 %
 % @details
 % Closes the white-box route of fj_tail_forktail for a fork branch that is
 % an M/G/1 FCFS queue, from the first three moments of its service time:
 %
 %   E[T] = E[S]*(1 + rho/(1-rho) * (1+SCV_S)/2)
 %   V[T] = E[W]^2 + lambda*E[S^3]/(3*(1-rho)) + E[S^2] - E[S]^2
 %
 % with rho = lambda*E[S] and E[W] = lambda*E[S^2]/(2*(1-rho)), i.e. the
 % Pollaczek-Khinchine mean waiting time. The third moment enters only the
 % variance, so an exponential or phase-type branch needs no extra input
 % beyond what the service distribution already reports (getMoments(3) on a
 % Markovian law, or 3*ES*VS + ES^3 + skewness*VS^(3/2) in general).
 %
 % @par Syntax:
 % @code
 % [ET, VT] = fj_mg1_respt_moments(lambda, ES, ES2, ES3)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>Arrival rate at the branch
 % <tr><td>ES<td>First moment of the service time
 % <tr><td>ES2<td>Second moment of the service time
 % <tr><td>ES3<td>Third moment of the service time
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>ET<td>Mean response time
 % <tr><td>VT<td>Variance of the response time
 % </table>
 %
 % @par References:
 % M. Nguyen, S. Alesawi, N. Li, H. Che, H. Jiang, "ForkTail: A Black-Box
 % Fork-Join Tail Latency Prediction Model for User-Facing Datacenter
 % Workloads", ACM HPDC 2018, Eqs. (10) and (11).
%}
function [ET, VT] = fj_mg1_respt_moments(lambda, ES, ES2, ES3)

rho = lambda * ES;
if rho >= 1
    line_error(mfilename, 'The branch is unstable (rho = %g >= 1); the response time moments do not exist.', rho);
end
if ~isfinite(ES3)
    line_error(mfilename, ['The service law has no finite third moment, so the ForkTail response ' ...
        'time variance is undefined (a Pareto branch with shape <= 3, for instance). ' ...
        'Use a service law with three finite moments.']);
end
scvS = (ES2 - ES^2) / ES^2;
ET = ES * (1 + rho/(1-rho) * (1 + scvS)/2);
EW = lambda * ES2 / (2*(1-rho));
VT = EW^2 + lambda*ES3/(3*(1-rho)) + ES2 - ES^2;
end
