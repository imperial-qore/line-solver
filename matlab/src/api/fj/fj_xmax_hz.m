%{ @file fj_xmax_hz.m
 %  @brief Harrison-Zertal approximation of the maximum of i.i.d. variables
 %
 %  @author LINE Development Team
%}

%{
 % @brief Harrison-Zertal approximation of the maximum of i.i.d. variables
 %
 % @details
 % Closed form obtained by collapsing the Harrison and Zertal recurrence onto
 % K identically distributed branches described by their first two moments:
 %
 %   X_K^max ~ m1 + ( m2 / (2*m1) ) * ( H_K - 1 ).
 %
 % The correction is the equilibrium (residual life) mean of the branch
 % distribution scaled by H_K - 1, so it reads as "one branch, plus the
 % expected residual work still owed by the branches that finish later".
 % Writing m2/(2*m1) = m1*(1+SCV)/2 gives the equivalent
 % m1 * [ 1 + (1+SCV)/2 * (H_K - 1) ], which is exact for the exponential
 % distribution (SCV = 1) and reduces to m1 at K = 1 for every branch law.
 %
 % @par Syntax:
 % @code
 % Xmax = fj_xmax_hz(m1, m2, K)
 % [Xmax, resid] = fj_xmax_hz(m1, m2, K)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>m1<td>Mean of the branch distribution (positive)
 % <tr><td>m2<td>Second moment of the branch distribution, m2 >= m1^2
 % <tr><td>K<td>Number of branches (positive integer)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Xmax<td>Approximate expected maximum
 % <tr><td>resid<td>Equilibrium mean m2/(2*m1) used as the per-branch inflation
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Eq. (46) and
 % the identically distributed specialization on page 17:23.
 %
 % Original: P. G. Harrison, S. Zertal, "Queueing Models of RAID Systems with
 % Maxima of Waiting Times", Performance Evaluation 64(7-8), 2007.
%}
function [Xmax, resid] = fj_xmax_hz(m1, m2, K)

if m1 <= 0
    line_error(mfilename, 'The branch mean m1 must be positive. Got m1=%g.', m1);
end
if m2 < m1^2
    line_error(mfilename, 'The second moment m2=%g is below m1^2=%g, which no distribution attains.', m2, m1^2);
end
if K < 1 || K ~= round(K)
    line_error(mfilename, 'K must be a positive integer. Got K=%g.', K);
end

% Equilibrium mean of the branch distribution
resid = m2 / (2 * m1);

Xmax = m1 + resid * (fj_harmonic(K) - 1);

end
