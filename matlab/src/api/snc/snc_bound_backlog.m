%{ @file snc_bound_backlog.m
 %  @brief Violation probability of a backlog level (stochastic network calculus)
 %
 %  @author LINE Development Team
%}

%{
 % @brief Violation probability of a backlog level (stochastic network calculus)
 %
 % @details
 % For a flow with arrival envelope (sigmaA,rhoA) served by an element with
 % service envelope (sigmaS,rhoS), the backlog B(t) of the stable station
 % rhoA < rhoS obeys, for every theta > 0,
 %
 %   P{B(t) > b} <= exp(-theta*(b-sigmaA-sigmaS)) / (1-exp(-theta*(rhoS-rhoA))),
 %
 % the union bound over the start of the backlogged period summed as a geometric
 % series on the unit-slot time axis. The returned value is the infimum over
 % theta, clipped at 1, and is an upper bound on the tail, never an estimate of
 % it: the decay rate is asymptotically exact and the prefactor is loose.
 %
 % @par Syntax:
 % @code
 % eps = snc_bound_backlog(arv,srv,b)
 % [eps,theta] = snc_bound_backlog(arv,srv,b,thetamax)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>arv<td>Handle theta -> [sigmaA,rhoA], e.g. an snc_env_* function
 % <tr><td>srv<td>Handle theta -> [sigmaS,rhoS], e.g. snc_srv_rate
 % <tr><td>b<td>Backlog level, units of work
 % <tr><td>thetamax<td>Optional upper end of the theta search, default 1e3
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>eps<td>Upper bound on P{B > b}, in [0,1]
 % <tr><td>theta<td>Minimizing theta, NaN if no feasible theta exists
 % </table>
 %
 % @par Reference:
 % M. Fidler and A. Rizk, "A Guide to the Stochastic Network Calculus",
 % IEEE Communications Surveys and Tutorials, 17(1):92-105, 2015, Sec. IV-B.
%}
function [eps, theta] = snc_bound_backlog(arv, srv, b, thetamax)

if nargin < 3
    line_error(mfilename, 'Usage: eps = snc_bound_backlog(arv,srv,b).');
end
if ~isa(arv, 'function_handle') || ~isa(srv, 'function_handle')
    line_error(mfilename, 'arv and srv must be function handles of theta.');
end
if b < 0
    line_error(mfilename, 'b must be nonnegative. Got %g.', b);
end
if nargin < 4
    thetamax = [];
end

[eps, theta] = snc_thetaopt(@(t) obj(t, arv, srv, b), thetamax);
if ~isfinite(eps) || eps > 1
    eps = 1; % no feasible theta, or the bound is vacuous at this level
end
end

function e = obj(theta, arv, srv, b)
[sA, rA] = arv(theta);
[sS, rS] = srv(theta);
if ~isfinite(sA) || ~isfinite(sS) || ~isfinite(rA) || ~isfinite(rS) || rS <= rA
    e = Inf;
    return
end
e = exp(-theta * (b - sA - sS)) / (1 - exp(-theta * (rS - rA)));
end
