%{ @file fj_xmax_exp.m
 %  @brief Expected maximum of K i.i.d. exponential random variables
 %
 %  @author LINE Development Team
%}

%{
 % @brief Expected maximum of K i.i.d. exponential random variables
 %
 % @details
 % Computes the expected value of the maximum of K independent and
 % identically distributed exponential random variables with rate mu
 % (mean service time x = 1/mu).
 %
 % X_K^max = H_K / mu = H_K * x
 %
 % where H_K is the K-th Harmonic number.
 %
 % This is a fundamental quantity in Fork-Join analysis. For a split-merge
 % (SM) queueing system, the maximum throughput is lambda_K^SM = 1 / X_K^max.
 %
 % For large K, X_K^max ≈ (ln(K) + gamma) / mu, where gamma ≈ 0.57721 is
 % the Euler-Mascheroni constant.
 %
 % @par Syntax:
 % @code
 % Xmax = fj_xmax_exp(K, mu)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of parallel servers (positive integer)
 % <tr><td>mu<td>Service rate (mean service time is 1/mu)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Xmax<td>Expected maximum service time X_K^max = H_K / mu
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Page 17:2 and Eq. (10) on page 17:11.
%}
function Xmax = fj_xmax_exp(K, mu)

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

if mu <= 0
    line_error(mfilename, 'Service rate mu must be positive. Got mu=%.4f.', mu);
end

% Compute harmonic number H_K
H_K = fj_harmonic(K);

% Expected maximum: X_K^max = H_K / mu
Xmax = H_K / mu;

end
