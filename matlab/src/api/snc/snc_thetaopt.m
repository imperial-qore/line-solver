%{ @file snc_thetaopt.m
 %  @brief Minimize a Chernoff bound over the free parameter theta
 %
 %  @author LINE Development Team
%}

%{
 % @brief Minimize a Chernoff bound over the free parameter theta
 %
 % @details
 % Every bound in the snc family holds for each theta > 0 for which the arrival
 % MGF is finite and the station is stable, so the reported bound is the
 % infimum over theta. FUN is evaluated on a logarithmic grid, non-finite
 % values (a diverging MGF, an unstable leftover rate) are discarded, and the
 % best grid point is refined by fminbnd in log10(theta). The two-stage search
 % is used because the feasible set is an interval whose endpoints are not
 % known in closed form once envelopes are composed, and a plain fminbnd over
 % the whole range would step into the infeasible region.
 %
 % @par Syntax:
 % @code
 % [val,theta] = snc_thetaopt(fun)
 % [val,theta] = snc_thetaopt(fun,thetamax)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>fun<td>Handle theta -> scalar objective to minimize
 % <tr><td>thetamax<td>Optional upper end of the search range, default 1e3
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>val<td>Minimum found, Inf if no feasible theta exists
 % <tr><td>theta<td>Minimizing theta, NaN if no feasible theta exists
 % </table>
%}
function [val, theta] = snc_thetaopt(fun, thetamax)

if nargin < 1
    line_error(mfilename, 'Usage: [val,theta] = snc_thetaopt(fun).');
end
if nargin < 2 || isempty(thetamax)
    thetamax = 1e3;
end
if thetamax <= 0
    line_error(mfilename, 'thetamax must be positive. Got %g.', thetamax);
end

grid = logspace(-6, log10(thetamax), 600);
fval = Inf(1, numel(grid));
for i = 1:numel(grid)
    fval(i) = safeval(fun, grid(i));
end

[val, imin] = min(fval);
if val >= 1e299
    val = Inf;
    theta = NaN;
    return
end
theta = grid(imin);

lo = grid(max(imin - 1, 1));
hi = grid(min(imin + 1, numel(grid)));
if hi > lo
    obj = @(x) safeval(fun, 10^x);
    [xopt, vopt] = fminbnd(obj, log10(lo), log10(hi), optimset('TolX', 1e-10));
    if vopt < val
        val = vopt;
        theta = 10^xopt;
    end
end
end

function v = safeval(fun, theta)
v = fun(theta);
if ~isscalar(v) || ~isreal(v) || ~isfinite(v)
    v = 1e300; % infeasible theta, kept finite so that fminbnd can compare it
end
end
