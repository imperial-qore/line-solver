%{ @file fj_harmonic.m
 %  @brief Compute Harmonic sum H_K = sum(1/k) for k=1 to K
 %
 %  @author LINE Development Team
%}

%{
 % @brief Compute Harmonic sum H_K = sum(1/k) for k=1 to K
 %
 % @details
 % Computes the K-th Harmonic number, which is the sum of reciprocals
 % from 1 to K. This is a fundamental quantity in Fork-Join analysis.
 %
 % For large K, H_K ≈ ln(K) + γ, where γ ≈ 0.57721 is the Euler-Mascheroni
 % constant.
 %
 % @par Syntax:
 % @code
 % H = fj_harmonic(K)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of parallel servers (positive integer)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>H<td>Harmonic sum H_K = 1 + 1/2 + 1/3 + ... + 1/K
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
%}
function H = fj_harmonic(K)

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

% Compute harmonic sum: H_K = sum_{k=1}^{K} 1/k
H = sum(1 ./ (1:K));

end
