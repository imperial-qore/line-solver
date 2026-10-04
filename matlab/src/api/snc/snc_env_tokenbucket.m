%{ @file snc_env_tokenbucket.m
 %  @brief Deterministic token-bucket arrival envelope
 %
 %  @author LINE Development Team
%}

%{
 % @brief Deterministic token-bucket arrival envelope
 %
 % @details
 % A flow policed by a (b,r) token bucket satisfies A(s,t) <= b + r*(t-s) with
 % probability one, hence E[exp(theta*A(s,t))] <= exp(theta*(r*(t-s)+b)) for
 % every theta > 0. The envelope is therefore constant in theta:
 %
 %   sigma(theta) = b,   rho(theta) = r.
 %
 % This is the deterministic network calculus arrival curve read as a degenerate
 % MGF envelope, so it can be mixed freely with the stochastic ones.
 %
 % @par Syntax:
 % @code
 % [sigma,rho] = snc_env_tokenbucket(b,r)
 % [sigma,rho] = snc_env_tokenbucket(b,r,theta)   % theta ignored, kept for
 %                                                % handle compatibility
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>b<td>Bucket depth, units of work
 % <tr><td>r<td>Token rate, work per slot
 % <tr><td>theta<td>Optional Chernoff parameter, unused
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sigma<td>Burst term of the envelope, equal to b
 % <tr><td>rho<td>Rate term of the envelope, equal to r
 % </table>
 %
 % @par Reference:
 % J.-Y. Le Boudec and P. Thiran, "Network Calculus", LNCS 2050, Springer, 2001.
%}
function [sigma, rho] = snc_env_tokenbucket(b, r, theta)

if nargin < 2
    line_error(mfilename, 'Usage: [sigma,rho] = snc_env_tokenbucket(b,r).');
end
if b < 0 || r < 0
    line_error(mfilename, 'b and r must be nonnegative. Got %g, %g.', b, r);
end
if nargin >= 3 && theta <= 0
    line_error(mfilename, 'theta must be positive. Got %g.', theta);
end

sigma = b;
rho = r;
end
