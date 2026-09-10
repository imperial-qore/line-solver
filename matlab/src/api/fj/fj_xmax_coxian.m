%{ @file fj_xmax_coxian.m
 %  @brief Expected maximum of K i.i.d. two-stage Coxian variables
 %
 %  @author LINE Development Team
%}

%{
 % @brief Expected maximum of K i.i.d. two-stage Coxian variables
 %
 % @details
 % Exact expected maximum of K independent copies of X = T1 + B*T2, where
 % T1 ~ Exp(mu1), T2 ~ Exp(mu2) and B is Bernoulli(q). The survival function
 % is a two-term exponential mixture,
 %
 %   S(t) = A*exp(-mu1 t) + B*exp(-mu2 t),
 %   A = (1-q) + q*mu2/(mu2-mu1),   B = -q*mu1/(mu2-mu1),
 %
 % so expanding 1 - (1-S)^K binomially and integrating term by term gives
 %
 %   E[Y_K] = sum_{j=1..K} (-1)^(j+1) binom(K,j)
 %              sum_{i=0..j} binom(j,i) A^(j-i) B^i / ((j-i)*mu1 + i*mu2).
 %
 % When the two stage rates coincide the mixture degenerates into
 % S(t) = (1 + q*mu*t)*exp(-mu*t) and the same expansion is carried out with
 % the polynomial integrals integral t^m exp(-j mu t) dt = m!/(j mu)^(m+1),
 % which the implementation selects automatically.
 %
 % @par Syntax:
 % @code
 % Xmax = fj_xmax_coxian(K, mu1, mu2, q)
 % [Xmax, m1, c2] = fj_xmax_coxian(K, mu1, mu2, q)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of branches (positive integer)
 % <tr><td>mu1<td>Rate of the first stage (positive)
 % <tr><td>mu2<td>Rate of the second stage (positive)
 % <tr><td>q<td>Probability that the second stage is visited, in [0,1]
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Xmax<td>Exact expected maximum of the K branches
 % <tr><td>m1<td>Mean of a single branch, 1/mu1 + q/mu2
 % <tr><td>c2<td>Squared coefficient of variation of a single branch
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Section 4.5,
 % Eqs. (38)-(40).
 %
 % Original: P. M. Chen, D. Towsley, "A Performance Evaluation of RAID
 % Architectures", IEEE Trans. Computers 45(10), 1996.
%}
function [Xmax, m1, c2] = fj_xmax_coxian(K, mu1, mu2, q)

if K < 1 || K ~= round(K)
    line_error(mfilename, 'K must be a positive integer. Got K=%g.', K);
end
if mu1 <= 0 || mu2 <= 0
    line_error(mfilename, 'Both stage rates must be positive. Got mu1=%g, mu2=%g.', mu1, mu2);
end
if q < 0 || q > 1
    line_error(mfilename, 'The branching probability q must lie in [0,1]. Got q=%g.', q);
end
if K > 60
    line_error(mfilename, 'The binomial expansion loses precision beyond K=60. Got K=%d.', K);
end

% Moments of a single two-stage Coxian branch
m1 = 1 / mu1 + q / mu2;
var1 = 1 / mu1^2 + q * (2 - q) / mu2^2;
c2 = var1 / m1^2;

Xmax = 0;
if abs(mu2 - mu1) > 1e-12 * max(mu1, mu2)
    % Distinct stage rates: the survival function is a two-exponential mixture
    A = (1 - q) + q * mu2 / (mu2 - mu1);
    B = -q * mu1 / (mu2 - mu1);
    for j = 1:K
        cj = nchoosek(K, j) * (-1)^(j + 1);
        inner = 0;
        for i = 0:j
            rate = (j - i) * mu1 + i * mu2;
            inner = inner + nchoosek(j, i) * A^(j - i) * B^i / rate;
        end
        Xmax = Xmax + cj * inner;
    end
else
    % Coincident stage rates: S(t) = (1 + q*mu*t)*exp(-mu*t)
    mu = mu1;
    for j = 1:K
        cj = nchoosek(K, j) * (-1)^(j + 1);
        inner = 0;
        for i = 0:j
            % Coefficient of (q*mu*t)^i inside S(t)^j, integrated against exp(-j mu t)
            inner = inner + nchoosek(j, i) * (q * mu)^i * factorial(i) / (j * mu)^(i + 1);
        end
        Xmax = Xmax + cj * inner;
    end
end

end
