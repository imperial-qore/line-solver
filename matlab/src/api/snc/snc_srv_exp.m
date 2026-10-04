%{ @file snc_srv_exp.m
 %  @brief MGF service envelope of an exponential server, in job units
 %
 %  @author LINE Development Team
%}

%{
 % @brief MGF service envelope of an exponential server, in job units
 %
 % @details
 % A single server with Exp(MU) service times completes jobs at the epochs of a
 % Poisson process of rate MU while it is busy, so its cumulative service S(s,t)
 % counted in JOBS is Poisson with mean mu*(t-s) and
 %
 %   log E[exp(-theta*S(s,t))] = mu*(t-s)*(exp(-theta)-1),
 %
 % giving the tight envelope rho(theta) = mu*(1-exp(-theta))/theta, sigma = 0.
 %
 % THIS IS THE SERVICE ELEMENT TO USE WHENEVER THE WORK UNIT IS THE JOB. Pairing
 % snc_srv_rate with a job-counting arrival envelope would model a server that
 % completes jobs at deterministic intervals, an M/D/1, and would UNDERSTATE the
 % delay of an exponential server rather than bound it. The M/M/1 read with this
 % element instead reproduces both exact decay rates: the backlog bound decays as
 % (lambda/mu)^n in jobs and the delay bound as exp(-(mu-lambda)*d) in time,
 % since the optimal theta tends to log(mu/lambda).
 %
 % Job units also compose across hops: a departure envelope from snc_output is a
 % job count and is directly the arrival envelope of the next station, whereas
 % service-time work units differ from station to station.
 %
 % @par Syntax:
 % @code
 % [sigma,rho] = snc_srv_exp(mu,theta)
 % srv = @(theta) snc_srv_exp(mu,theta);
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>mu<td>Service rate, jobs per slot
 % <tr><td>theta<td>Chernoff parameter, theta > 0
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sigma<td>Burst term of the service envelope (0 here)
 % <tr><td>rho<td>Rate term of the service envelope, jobs per slot
 % </table>
 %
 % @par Reference:
 % M. Fidler and A. Rizk, "A Guide to the Stochastic Network Calculus",
 % IEEE Communications Surveys and Tutorials, 17(1):92-105, 2015, Sec. IV.
%}
function [sigma, rho] = snc_srv_exp(mu, theta)

if nargin < 2
    line_error(mfilename, 'Usage: [sigma,rho] = snc_srv_exp(mu,theta).');
end
if mu <= 0
    line_error(mfilename, 'mu must be positive. Got %g.', mu);
end
if theta <= 0
    line_error(mfilename, 'theta must be positive. Got %g.', theta);
end

sigma = 0;
rho = mu * (1 - exp(-theta)) / theta;
end
