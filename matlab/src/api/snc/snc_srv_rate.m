%{ @file snc_srv_rate.m
 %  @brief MGF service envelope of a constant-rate work-conserving server
 %
 %  @author LINE Development Team
%}

%{
 % @brief MGF service envelope of a constant-rate work-conserving server
 %
 % @details
 % A work-conserving server of capacity C offers the deterministic service
 % process S(s,t) = C*(t-s), whose exponential-form envelope
 %
 %   E[exp(-theta*S(s,t))] <= exp(-theta*(rho*(t-s) - sigma))
 %
 % holds with rho(theta) = C and sigma(theta) = 0 for every theta > 0. This is
 % the elementary service element of the snc family; a station shared by cross
 % traffic is obtained from it through snc_leftover, and a tandem of stations
 % through snc_conv.
 %
 % @par Syntax:
 % @code
 % [sigma,rho] = snc_srv_rate(C)
 % [sigma,rho] = snc_srv_rate(C,theta)   % theta ignored, kept for handle
 %                                       % compatibility
 % srv = @(theta) snc_srv_rate(C,theta);
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>C<td>Server capacity, work per slot
 % <tr><td>theta<td>Optional Chernoff parameter, unused
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sigma<td>Burst term of the service envelope (0 here)
 % <tr><td>rho<td>Rate term of the service envelope, equal to C
 % </table>
 %
 % @par Reference:
 % M. Fidler and A. Rizk, "A Guide to the Stochastic Network Calculus",
 % IEEE Communications Surveys and Tutorials, 17(1):92-105, 2015, Sec. IV.
%}
function [sigma, rho] = snc_srv_rate(C, theta)

if nargin < 1
    line_error(mfilename, 'Usage: [sigma,rho] = snc_srv_rate(C).');
end
if C <= 0
    line_error(mfilename, 'C must be positive. Got %g.', C);
end
if nargin >= 2 && theta <= 0
    line_error(mfilename, 'theta must be positive. Got %g.', theta);
end

sigma = 0;
rho = C;
end
