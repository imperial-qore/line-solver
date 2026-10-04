%{ @file snc_mean_delay.m
 %  @brief Upper bound on the mean delay from the tail bound
 %
 %  @author LINE Development Team
%}

%{
 % @brief Upper bound on the mean delay from the tail bound
 %
 % @details
 % Integrates the delay tail bound of snc_bound_delay over the whole positive
 % axis, which for a nonnegative delay gives E[D] = int_0^inf P{D>d} dd and
 % hence an upper bound on the MEAN. At fixed theta the bound is K*exp(-a*d)
 % with a = theta*rhoS and
 %
 %   K = exp(theta*(sigmaA+sigmaS)) / (1-exp(-theta*(rhoS-rhoA))),
 %
 % so, clipping the bound at 1 where it exceeds it, the integral is available in
 % CLOSED FORM: (log(K)+1)/a when K >= 1 and K/a otherwise. No quadrature is
 % involved, so the result is a bound and not a bound plus a discretization
 % error; the only numerical step is the minimization over theta.
 %
 % This is the entry point SolverBA calls for the 'snc.upper' response-time
 % column. IT IS A LOOSE MEAN BOUND AND THAT IS INHERENT: on the M/M/1 read in
 % job units (snc_env_poisson feeding snc_srv_exp) it returns 2.4x the exact
 % 1/(mu-lambda) at rho = 0.1 and 10.4x at rho = 0.95, because the prefactor of
 % the tail bound, not its decay rate, is what dominates an integral over the
 % whole axis. The tail bound it integrates is the sharp object; use
 % snc_perc_delay when the quantile is what matters.
 %
 % @par Syntax:
 % @code
 % ED = snc_mean_delay(arv,srv)
 % [ED,theta] = snc_mean_delay(arv,srv,thetamax)
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
 % <tr><td>ED<td>Upper bound on E[D], Inf if no feasible theta exists
 % <tr><td>theta<td>Minimizing theta, NaN if no feasible theta exists
 % </table>
 %
 % @par Reference:
 % M. Fidler and A. Rizk, "A Guide to the Stochastic Network Calculus",
 % IEEE Communications Surveys and Tutorials, 17(1):92-105, 2015, Sec. IV-B.
%}
function [ED, theta] = snc_mean_delay(arv, srv, thetamax)

if nargin < 2
    line_error(mfilename, 'Usage: ED = snc_mean_delay(arv,srv).');
end
if ~isa(arv, 'function_handle') || ~isa(srv, 'function_handle')
    line_error(mfilename, 'arv and srv must be function handles of theta.');
end
if nargin < 3
    thetamax = [];
end

[ED, theta] = snc_thetaopt(@(t) obj(t, arv, srv), thetamax);
end

function x = obj(theta, arv, srv)
[sA, rA] = arv(theta);
[sS, rS] = srv(theta);
if ~isfinite(sA) || ~isfinite(sS) || ~isfinite(rA) || ~isfinite(rS) || rS <= rA
    x = Inf;
    return
end
logK = theta * (sA + sS) - log(1 - exp(-theta * (rS - rA)));
a = theta * rS;
if logK >= 0
    x = (logK + 1) / a;
else
    x = exp(logK) / a;
end
end
