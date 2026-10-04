%{ @file snc_output.m
 %  @brief Output arrival envelope of a flow leaving a server
 %
 %  @author LINE Development Team
%}

%{
 % @brief Output arrival envelope of a flow leaving a server
 %
 % @details
 % The departures of a flow with arrival envelope (sigmaA,rhoA) from a server
 % with service envelope (sigmaS,rhoS) admit, for independent processes and a
 % stable station rhoA < rhoS, the arrival envelope
 %
 %   rho = rhoA,
 %   sigma = sigmaA + sigmaS - log(1-exp(-theta*(rhoS-rhoA)))/theta.
 %
 % The rate is conserved and the server adds burstiness. This is what carries a
 % flow across a feed-forward network one hop at a time; for a tandem traversed
 % by the same flow, snc_conv gives the tighter end-to-end result.
 %
 % @par Syntax:
 % @code
 % [sigma,rho] = snc_output(sigmaA,rhoA,sigmaS,rhoS,theta)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sigmaA,rhoA<td>Arrival envelope of the flow entering the server
 % <tr><td>sigmaS,rhoS<td>Service envelope of the server
 % <tr><td>theta<td>Chernoff parameter, theta > 0
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sigma<td>Burst term of the departure envelope, Inf if unstable
 % <tr><td>rho<td>Rate term of the departure envelope, equal to rhoA
 % </table>
 %
 % @par Reference:
 % M. Fidler and A. Rizk, "A Guide to the Stochastic Network Calculus",
 % IEEE Communications Surveys and Tutorials, 17(1):92-105, 2015, Sec. V-A.
%}
function [sigma, rho] = snc_output(sigmaA, rhoA, sigmaS, rhoS, theta)

if nargin < 5
    line_error(mfilename, 'Usage: [sigma,rho] = snc_output(sigmaA,rhoA,sigmaS,rhoS,theta).');
end
if theta <= 0
    line_error(mfilename, 'theta must be positive. Got %g.', theta);
end

rho = rhoA;
if ~isfinite(rhoA) || ~isfinite(rhoS) || rhoS <= rhoA
    sigma = Inf; % unstable station, no exponential-form departure envelope
    return
end

sigma = sigmaA + sigmaS - log(1 - exp(-theta * (rhoS - rhoA))) / theta;
end
