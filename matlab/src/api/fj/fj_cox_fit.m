%{ @file fj_cox_fit.m
 %  @brief Two-stage Coxian fit of a mean and a squared coefficient of variation
 %
 %  @author LINE Development Team
%}

%{
 % @brief Two-stage Coxian fit of a mean and a squared coefficient of variation
 %
 % @details
 % Marie's balanced-stage fit of a two-stage Coxian law to a target mean 1/mu
 % and squared coefficient of variation c2. The two stages are required to
 % contribute equally to the mean, 1/mu1 = q/mu2, which closes the system of
 % two moment equations in three unknowns and yields
 %
 %   mu1 = 2*mu,   q = 1/(2*c2),   mu2 = 2*mu*q = mu/c2.
 %
 % The fit needs q <= 1, hence c2 >= 0.5; below that the balanced-stage
 % condition is infeasible and an Erlang stage count is the natural choice
 % instead. The admissible number of exponential stages of an Erlang or
 % same-rate Coxian representation of the same target is bracketed by
 %
 %   ceil(1/c2) <= k <= floor(1/c2) + 1,
 %
 % which is also returned.
 %
 % @par Syntax:
 % @code
 % [mu1, mu2, q] = fj_cox_fit(m1, c2)
 % [mu1, mu2, q, kmin, kmax] = fj_cox_fit(m1, c2)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>m1<td>Target mean (positive)
 % <tr><td>c2<td>Target squared coefficient of variation, c2 >= 0.5
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>mu1<td>Rate of the first stage
 % <tr><td>mu2<td>Rate of the second stage
 % <tr><td>q<td>Probability that the second stage is visited
 % <tr><td>kmin<td>Smallest admissible Erlang stage count for the same c2
 % <tr><td>kmax<td>Largest admissible Erlang stage count for the same c2
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Section 4.5,
 % Eqs. (39)-(40) and the discussion of Algorithm 3.2.5.
 %
 % Original: R. Marie, "Methodes iteratives de resolution de modeles
 % mathematiques de systemes informatiques", RAIRO Informatique 12(2), 1978;
 % A. O. Allen, "Probability, Statistics and Queueing Theory", Academic Press,
 % 2nd ed., 1990, Algorithm 3.2.5.
%}
function [mu1, mu2, q, kmin, kmax] = fj_cox_fit(m1, c2)

if m1 <= 0
    line_error(mfilename, 'The target mean m1 must be positive. Got m1=%g.', m1);
end
if c2 < 0.5
    line_error(mfilename, 'The balanced-stage Coxian fit needs c2 >= 0.5. Got c2=%g; use an Erlang with %d stages instead.', c2, max(1, ceil(1 / c2)));
end

mu = 1 / m1;

% Balanced-stage conditions: both stages contribute half of the mean
mu1 = 2 * mu;
q = 1 / (2 * c2);
mu2 = mu / c2;

% Erlang stage-count bracket for the same squared coefficient of variation
kmin = ceil(1 / c2);
kmax = floor(1 / c2) + 1;
kmin = max(kmin, 1);
kmax = max(kmax, kmin);

end
