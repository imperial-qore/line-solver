%{ @file snc_leftover.m
 %  @brief Leftover service envelope under blind multiplexing
 %
 %  @author LINE Development Team
%}

%{
 % @brief Leftover service envelope under blind multiplexing
 %
 % @details
 % A server with envelope (sigmaS,rhoS) shared with a cross flow of arrival
 % envelope (sigmaX,rhoX) leaves the flow of interest, under blind (arbitrary)
 % multiplexing, the service process S(s,t)-X(s,t), whose envelope is
 %
 %   rho = rhoS - rhoX,   sigma = sigmaS + sigmaX.
 %
 % The subtraction of the two exponential forms is exact when the two processes
 % are independent; otherwise the pair must be split by Hoelder's inequality,
 % which this elementary version does not do. A nonpositive rho means the cross
 % traffic can exhaust the server, and the bound functions return a violation
 % probability of 1 in that case rather than a number.
 %
 % @par Syntax:
 % @code
 % [sigma,rho] = snc_leftover(sigmaS,rhoS,sigmaX,rhoX)
 %
 % % leftover of a rate-C server shared with a Poisson cross flow, as the
 % % service handle expected by snc_bound_delay:
 % function [s,r] = srv(theta)
 %   [sX,rX] = snc_env_poisson(lambdaX,theta);
 %   [sS,rS] = snc_srv_rate(C,theta);
 %   [s,r] = snc_leftover(sS,rS,sX,rX);
 % end
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sigmaS<td>Burst term of the service envelope
 % <tr><td>rhoS<td>Rate term of the service envelope
 % <tr><td>sigmaX<td>Burst term of the cross-flow arrival envelope
 % <tr><td>rhoX<td>Rate term of the cross-flow arrival envelope
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sigma<td>Burst term of the leftover service envelope
 % <tr><td>rho<td>Rate term of the leftover service envelope
 % </table>
 %
 % @par Reference:
 % M. Fidler and A. Rizk, "A Guide to the Stochastic Network Calculus",
 % IEEE Communications Surveys and Tutorials, 17(1):92-105, 2015, Sec. V-B.
%}
function [sigma, rho] = snc_leftover(sigmaS, rhoS, sigmaX, rhoX)

if nargin < 4
    line_error(mfilename, 'Usage: [sigma,rho] = snc_leftover(sigmaS,rhoS,sigmaX,rhoX).');
end

sigma = sigmaS + sigmaX;
rho = rhoS - rhoX;
end
