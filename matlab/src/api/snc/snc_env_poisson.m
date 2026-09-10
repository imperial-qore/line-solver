%{ @file snc_env_poisson.m
 %  @brief MGF arrival envelope of a Poisson flow with unit-size jobs
 %
 %  @author LINE Development Team
%}

%{
 % @brief MGF arrival envelope of a Poisson flow with unit-size jobs
 %
 % @details
 % Returns the pair (sigma,rho) such that the cumulative arrivals A(s,t) of a
 % Poisson process of rate LAMBDA, each job carrying one unit of work, satisfy
 % the (sigma(theta),rho(theta))-constrained bound
 %
 %   E[exp(theta*A(s,t))] <= exp(theta*(rho*(t-s) + sigma)),   theta > 0.
 %
 % Since log E[exp(theta*A(0,t))] = lambda*t*(exp(theta)-1) exactly, the
 % envelope is tight with a zero burst term:
 %
 %   rho(theta) = lambda*(exp(theta)-1)/theta,   sigma(theta) = 0.
 %
 % Time is slotted with unit slot length, the convention of the whole snc
 % family; rho is therefore work per slot.
 %
 % @par Syntax:
 % @code
 % [sigma,rho] = snc_env_poisson(lambda,theta)
 % arv = @(theta) snc_env_poisson(lambda,theta);   % handle form for the bounds
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>Arrival rate, jobs per slot
 % <tr><td>theta<td>Chernoff parameter, theta > 0
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sigma<td>Burst term of the envelope (0 here)
 % <tr><td>rho<td>Rate term of the envelope, work per slot
 % </table>
 %
 % @par Reference:
 % M. Fidler and A. Rizk, "A Guide to the Stochastic Network Calculus",
 % IEEE Communications Surveys and Tutorials, 17(1):92-105, 2015, Sec. IV.
%}
function [sigma, rho] = snc_env_poisson(lambda, theta)

if nargin < 2
    line_error(mfilename, 'Usage: [sigma,rho] = snc_env_poisson(lambda,theta).');
end
if lambda < 0
    line_error(mfilename, 'lambda must be nonnegative. Got %g.', lambda);
end
if theta <= 0
    line_error(mfilename, 'theta must be positive. Got %g.', theta);
end

sigma = 0;
rho = lambda * (exp(theta) - 1) / theta;
end
