%{ @file snc_env_cpoisson.m
 %  @brief MGF arrival envelope of a compound Poisson flow with Exp job sizes
 %
 %  @author LINE Development Team
%}

%{
 % @brief MGF arrival envelope of a compound Poisson flow with Exp job sizes
 %
 % @details
 % Jobs arrive Poisson at rate LAMBDA and each carries an Exp(MU) amount of
 % work, so the cumulative work A(s,t) is a compound Poisson process with
 %
 %   log E[exp(theta*A(0,t))] = lambda*t*theta/(mu-theta),   0 < theta < mu,
 %
 % giving the tight envelope rho(theta) = lambda/(mu-theta), sigma(theta) = 0.
 % The MGF diverges at theta >= mu, where rho is returned as Inf so that the
 % theta search in snc_thetaopt discards the point.
 %
 % Fed to a constant-rate server of rate MU (snc_srv_rate) this is the network
 % calculus model of the M/M/1 queue: the resulting delay bound decays at rate
 % mu-lambda, the exact asymptotic decay rate of the M/M/1 waiting time.
 %
 % @par Syntax:
 % @code
 % [sigma,rho] = snc_env_cpoisson(lambda,mu,theta)
 % arv = @(theta) snc_env_cpoisson(lambda,mu,theta);
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>Job arrival rate, jobs per slot
 % <tr><td>mu<td>Rate of the Exp job size, so mean work per job is 1/mu
 % <tr><td>theta<td>Chernoff parameter, theta > 0
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sigma<td>Burst term of the envelope (0 here)
 % <tr><td>rho<td>Rate term of the envelope, Inf when theta >= mu
 % </table>
 %
 % @par Reference:
 % M. Fidler and A. Rizk, "A Guide to the Stochastic Network Calculus",
 % IEEE Communications Surveys and Tutorials, 17(1):92-105, 2015, Sec. IV.
%}
function [sigma, rho] = snc_env_cpoisson(lambda, mu, theta)

if nargin < 3
    line_error(mfilename, 'Usage: [sigma,rho] = snc_env_cpoisson(lambda,mu,theta).');
end
if lambda < 0 || mu <= 0
    line_error(mfilename, 'lambda must be nonnegative and mu positive. Got %g, %g.', lambda, mu);
end
if theta <= 0
    line_error(mfilename, 'theta must be positive. Got %g.', theta);
end

sigma = 0;
if theta >= mu
    rho = Inf; % the job-size MGF diverges, no envelope at this theta
else
    rho = lambda / (mu - theta);
end
end
