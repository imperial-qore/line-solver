%{ @file snc_perc_backlog.m
 %  @brief Backlog quantile at a prescribed violation probability
 %
 %  @author LINE Development Team
%}

%{
 % @brief Backlog quantile at a prescribed violation probability
 %
 % @details
 % Inverts the backlog bound of snc_bound_backlog in b: the smallest level for
 % which the bound certifies P{B > b} <= EPS is, at fixed theta,
 %
 %   b(theta) = sigmaA + sigmaS - log(eps*(1-exp(-theta*(rhoS-rhoA))))/theta,
 %
 % and the reported quantile is the minimum of b(theta) over the feasible
 % thetas. Note that the minimizing theta differs from the one of the forward
 % bound at a given level, which is why the inversion is done in closed form and
 % re-optimized rather than by a search on snc_bound_backlog.
 %
 % @par Syntax:
 % @code
 % b = snc_perc_backlog(arv,srv,eps)
 % [b,theta] = snc_perc_backlog(arv,srv,eps,thetamax)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>arv<td>Handle theta -> [sigmaA,rhoA], e.g. an snc_env_* function
 % <tr><td>srv<td>Handle theta -> [sigmaS,rhoS], e.g. snc_srv_rate
 % <tr><td>eps<td>Violation probability, 0 < eps < 1
 % <tr><td>thetamax<td>Optional upper end of the theta search, default 1e3
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>b<td>Backlog quantile, units of work, Inf if the station is unstable
 % <tr><td>theta<td>Minimizing theta, NaN if no feasible theta exists
 % </table>
 %
 % @par Reference:
 % M. Fidler and A. Rizk, "A Guide to the Stochastic Network Calculus",
 % IEEE Communications Surveys and Tutorials, 17(1):92-105, 2015, Sec. IV-B.
%}
function [b, theta] = snc_perc_backlog(arv, srv, eps, thetamax)

if nargin < 3
    line_error(mfilename, 'Usage: b = snc_perc_backlog(arv,srv,eps).');
end
if ~isa(arv, 'function_handle') || ~isa(srv, 'function_handle')
    line_error(mfilename, 'arv and srv must be function handles of theta.');
end
if eps <= 0 || eps >= 1
    line_error(mfilename, 'eps must lie in (0,1). Got %g.', eps);
end
if nargin < 4
    thetamax = [];
end

[b, theta] = snc_thetaopt(@(t) obj(t, arv, srv, eps), thetamax);
b = max(b, 0);
end

function x = obj(theta, arv, srv, eps)
[sA, rA] = arv(theta);
[sS, rS] = srv(theta);
if ~isfinite(sA) || ~isfinite(sS) || ~isfinite(rA) || ~isfinite(rS) || rS <= rA
    x = Inf;
    return
end
x = sA + sS - log(eps * (1 - exp(-theta * (rS - rA)))) / theta;
end
