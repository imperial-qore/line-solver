%{ @file snc_conv.m
 %  @brief Min-plus convolution of two service envelopes (tandem concatenation)
 %
 %  @author LINE Development Team
%}

%{
 % @brief Min-plus convolution of two service envelopes (tandem concatenation)
 %
 % @details
 % Two stations traversed in series offer the flow their min-plus convolution.
 % For independent servers with exponential-form envelopes (sigma1,rho1) and
 % (sigma2,rho2), summing the geometric series over the intermediate epoch gives
 %
 %   rho = min(rho1,rho2),
 %   sigma = sigma1 + sigma2 - log(1-exp(-theta*|rho1-rho2|))/theta.
 %
 % This is the pay-bursts-only-once result: the end-to-end burst term grows
 % additively rather than the delay bounds of the two stations being summed.
 %
 % The series diverges when rho1 = rho2, so equal rates are handled by shifting
 % the slower server down by DELTA, which is the usual regularization; DELTA
 % then trades rate against burst and can be optimized jointly with theta.
 %
 % @par Syntax:
 % @code
 % [sigma,rho] = snc_conv(sigma1,rho1,sigma2,rho2,theta)
 % [sigma,rho] = snc_conv(sigma1,rho1,sigma2,rho2,theta,delta)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sigma1,rho1<td>Envelope of the first station
 % <tr><td>sigma2,rho2<td>Envelope of the second station
 % <tr><td>theta<td>Chernoff parameter, theta > 0
 % <tr><td>delta<td>Optional rate separation used when rho1 = rho2, default
 %                  1e-2*min(rho1,rho2)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sigma<td>Burst term of the concatenated service envelope
 % <tr><td>rho<td>Rate term of the concatenated service envelope
 % </table>
 %
 % @par Reference:
 % M. Fidler and A. Rizk, "A Guide to the Stochastic Network Calculus",
 % IEEE Communications Surveys and Tutorials, 17(1):92-105, 2015, Sec. V-C.
 %
 % Original: F. Ciucu, A. Burchard, J. Liebeherr, "Scaling Properties of
 % Statistical End-to-End Bounds in the Network Calculus", IEEE Trans. on
 % Information Theory, 52(6):2300-2312, 2006.
%}
function [sigma, rho] = snc_conv(sigma1, rho1, sigma2, rho2, theta, delta)

if nargin < 5
    line_error(mfilename, 'Usage: [sigma,rho] = snc_conv(sigma1,rho1,sigma2,rho2,theta).');
end
if theta <= 0
    line_error(mfilename, 'theta must be positive. Got %g.', theta);
end
if nargin < 6 || isempty(delta)
    delta = 1e-2 * min(rho1, rho2);
end
if delta <= 0
    line_error(mfilename, 'delta must be positive. Got %g.', delta);
end

if ~isfinite(rho1) || ~isfinite(rho2)
    sigma = Inf;
    rho = min(rho1, rho2);
    return
end

gap = abs(rho1 - rho2);
if gap <= delta
    gap = delta; % equal rates: shift the slower server down to close the series
    rho = min(rho1, rho2) - delta;
else
    rho = min(rho1, rho2);
end

if rho <= 0
    sigma = Inf;
    return
end

sigma = sigma1 + sigma2 - log(1 - exp(-theta * gap)) / theta;
end
