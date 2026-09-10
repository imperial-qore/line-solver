%{ @file snc_mean_backlog.m
 %  @brief Upper bound on the mean backlog from the tail bound
 %
 %  @author LINE Development Team
%}

%{
 % @brief Upper bound on the mean backlog from the tail bound
 %
 % @details
 % The counterpart of snc_mean_delay for the backlog: integrating the tail
 % bound of snc_bound_backlog gives E[B] = int_0^inf P{B>b} db and hence an
 % upper bound on the mean. At fixed theta the bound is K*exp(-theta*b) with
 %
 %   K = exp(theta*(sigmaA+sigmaS)) / (1-exp(-theta*(rhoS-rhoA))),
 %
 % so the clipped integral is (log(K)+1)/theta when K >= 1 and K/theta
 % otherwise, in closed form. The unit of the answer is the unit of the
 % envelopes: jobs when the pair is snc_env_poisson with snc_srv_exp, units of
 % work when it is snc_env_cpoisson with snc_srv_rate.
 %
 % SolverBA does NOT use this for its queue-length column: it applies Little's
 % law to the response-time bound instead, so that Q and R stay consistent with
 % the exact open-network throughput. The two are close but not identical,
 % since each optimizes its own theta.
 %
 % @par Syntax:
 % @code
 % EB = snc_mean_backlog(arv,srv)
 % [EB,theta] = snc_mean_backlog(arv,srv,thetamax)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>arv<td>Handle theta -> [sigmaA,rhoA], e.g. an snc_env_* function
 % <tr><td>srv<td>Handle theta -> [sigmaS,rhoS], e.g. snc_srv_exp
 % <tr><td>thetamax<td>Optional upper end of the theta search, default 1e3
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>EB<td>Upper bound on E[B], Inf if no feasible theta exists
 % <tr><td>theta<td>Minimizing theta, NaN if no feasible theta exists
 % </table>
 %
 % @par Reference:
 % M. Fidler and A. Rizk, "A Guide to the Stochastic Network Calculus",
 % IEEE Communications Surveys and Tutorials, 17(1):92-105, 2015, Sec. IV-B.
%}
function [EB, theta] = snc_mean_backlog(arv, srv, thetamax)

if nargin < 2
    line_error(mfilename, 'Usage: EB = snc_mean_backlog(arv,srv).');
end
if ~isa(arv, 'function_handle') || ~isa(srv, 'function_handle')
    line_error(mfilename, 'arv and srv must be function handles of theta.');
end
if nargin < 3
    thetamax = [];
end

[EB, theta] = snc_thetaopt(@(t) obj(t, arv, srv), thetamax);
end

function x = obj(theta, arv, srv)
[sA, rA] = arv(theta);
[sS, rS] = srv(theta);
if ~isfinite(sA) || ~isfinite(sS) || ~isfinite(rA) || ~isfinite(rS) || rS <= rA
    x = Inf;
    return
end
logK = theta * (sA + sS) - log(1 - exp(-theta * (rS - rA)));
if logK >= 0
    x = (logK + 1) / theta;
else
    x = exp(logK) / theta;
end
end
