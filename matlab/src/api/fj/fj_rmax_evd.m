%{ @file fj_rmax_evd.m
 %  @brief Maximum response time using Extreme Value Distribution approximation
 %
 %  @author LINE Development Team
%}

%{
 % @brief Maximum response time using Extreme Value Distribution approximation
 %
 % @details
 % Approximates the expected value of the maximum response time for K
 % parallel servers using the Extreme Value Distribution (EVD).
 %
 % R_K^max(rho) = R(rho) + (sqrt(6)*ln(K)/pi) * sigma_R(rho)
 %
 % where R(rho) is the mean response time and sigma_R(rho) is the standard
 % deviation of response time.
 %
 % This approximation is based on fitting the response time distribution
 % to a Gumbel (Type I extreme value) distribution with:
 %   - Location parameter: a = R - gamma*b
 %   - Scale parameter: b = sqrt(6)*sigma_R/pi
 %
 % where gamma ≈ 0.5772 is the Euler-Mascheroni constant.
 %
 % Note: Thomasian et al. [2007] found that dividing the correction term
 % by 1.27 improves accuracy.
 %
 % @par Syntax:
 % @code
 % Rmax = fj_rmax_evd(K, R, sigma_R)
 % Rmax = fj_rmax_evd(K, R, sigma_R, calibrated)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of parallel servers (positive integer)
 % <tr><td>R<td>Mean response time
 % <tr><td>sigma_R<td>Standard deviation of response time
 % <tr><td>calibrated<td>Optional: use calibrated formula (default: false)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Rmax<td>Approximate maximum response time R_K^max
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Eq. (43) on page 17:21.
%}
function Rmax = fj_rmax_evd(K, R, sigma_R, calibrated)

if nargin < 4
    calibrated = false;
end

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

if R <= 0
    line_error(mfilename, 'Mean response time R must be positive.');
end

if sigma_R < 0
    line_error(mfilename, 'Standard deviation sigma_R must be non-negative.');
end

% EVD approximation: R_K^max = R + (sqrt(6)*ln(K)/pi) * sigma_R
correction_factor = sqrt(6) * log(K) / pi;

if calibrated
    % Calibrated version from Thomasian et al. [2007]
    % Divides the correction term by 1.27 for improved accuracy
    correction_factor = correction_factor / 1.27;
end

Rmax = R + correction_factor * sigma_R;

end
