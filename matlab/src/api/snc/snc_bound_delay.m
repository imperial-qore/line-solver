%{ @file snc_bound_delay.m
 %  @brief Violation probability of a delay target (stochastic network calculus)
 %
 %  @author LINE Development Team
%}

%{
 % @brief Violation probability of a delay target (stochastic network calculus)
 %
 % @details
 % For a flow with arrival envelope (sigmaA,rhoA) served by an element with
 % service envelope (sigmaS,rhoS), the virtual delay D(t) of the stable station
 % rhoA < rhoS obeys, for every theta > 0,
 %
 %   P{D(t) > d} <= exp(-theta*(rhoS*d-sigmaA-sigmaS)) /
 %                  (1-exp(-theta*(rhoS-rhoA))),
 %
 % the horizontal rather than vertical deviation between the arrival and
 % service envelopes. The returned value is the infimum over theta, clipped at
 % 1. On the M/M/1 parameterization (snc_env_cpoisson feeding snc_srv_rate with
 % C = mu) the optimal theta approaches mu-lambda, so the bound reproduces the
 % exact asymptotic decay rate of the waiting-time tail.
 %
 % @par Syntax:
 % @code
 % eps = snc_bound_delay(arv,srv,d)
 % [eps,theta] = snc_bound_delay(arv,srv,d,thetamax)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>arv<td>Handle theta -> [sigmaA,rhoA], e.g. an snc_env_* function
 % <tr><td>srv<td>Handle theta -> [sigmaS,rhoS], e.g. snc_srv_rate
 % <tr><td>d<td>Delay target, slots
 % <tr><td>thetamax<td>Optional upper end of the theta search, default 1e3
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>eps<td>Upper bound on P{D > d}, in [0,1]
 % <tr><td>theta<td>Minimizing theta, NaN if no feasible theta exists
 % </table>
 %
 % @par Reference:
 % M. Fidler and A. Rizk, "A Guide to the Stochastic Network Calculus",
 % IEEE Communications Surveys and Tutorials, 17(1):92-105, 2015, Sec. IV-B.
%}
function [eps, theta] = snc_bound_delay(arv, srv, d, thetamax)

if nargin < 3
    line_error(mfilename, 'Usage: eps = snc_bound_delay(arv,srv,d).');
end
if ~isa(arv, 'function_handle') || ~isa(srv, 'function_handle')
    line_error(mfilename, 'arv and srv must be function handles of theta.');
end
if d < 0
    line_error(mfilename, 'd must be nonnegative. Got %g.', d);
end
if nargin < 4
    thetamax = [];
end

[eps, theta] = snc_thetaopt(@(t) obj(t, arv, srv, d), thetamax);
if ~isfinite(eps) || eps > 1
    eps = 1; % no feasible theta, or the bound is vacuous at this target
end
end

function e = obj(theta, arv, srv, d)
[sA, rA] = arv(theta);
[sS, rS] = srv(theta);
if ~isfinite(sA) || ~isfinite(sS) || ~isfinite(rA) || ~isfinite(rS) || rS <= rA
    e = Inf;
    return
end
e = exp(-theta * (rS * d - sA - sS)) / (1 - exp(-theta * (rS - rA)));
end
