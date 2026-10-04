%{ @file fj_xmax_hyperexp.m
 %  @brief Expected maximum of K i.i.d. Hyperexponential-2 random variables
 %
 %  @author LINE Development Team
%}

%{
 % @brief Expected maximum of K i.i.d. Hyperexponential-2 random variables
 %
 % @details
 % Computes the expected value of the maximum of K independent and
 % identically distributed Hyperexponential random variables with
 % two branches.
 %
 % The pdf of H2 is: f(t) = p1*mu1*exp(-mu1*t) + p2*mu2*exp(-mu2*t)
 % where p1 + p2 = 1 and p1, p2 > 0.
 %
 % The expected maximum is:
 %   X_K^max = sum_{n=1}^{K} (-1)^{n+1} * sum_{m=0}^{n} C(n,m) * p1^m * p2^{n-m} / (m*mu1 + (n-m)*mu2)
 %
 % The H2 distribution has CV > 1, making it useful for modeling
 % high-variability service times.
 %
 % @par Syntax:
 % @code
 % Xmax = fj_xmax_hyperexp(K, p1, mu1, mu2)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of parallel servers (positive integer)
 % <tr><td>p1<td>Probability of branch 1 (0 < p1 < 1)
 % <tr><td>mu1<td>Rate of branch 1
 % <tr><td>mu2<td>Rate of branch 2
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Xmax<td>Expected maximum service time X_K^max
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
 % Page 17:13, Eq. (15).
%}
function Xmax = fj_xmax_hyperexp(K, p1, mu1, mu2)

if K < 1
    line_error(mfilename, 'K must be a positive integer. Got K=%d.', K);
end

if p1 <= 0 || p1 >= 1
    line_error(mfilename, 'p1 must be in (0,1). Got p1=%.4f.', p1);
end

if mu1 <= 0 || mu2 <= 0
    line_error(mfilename, 'Rates mu1 and mu2 must be positive.');
end

p2 = 1 - p1;

% X_K^max = sum_{n=1}^{K} (-1)^{n+1} * sum_{m=0}^{n} C(n,m) * p1^m * p2^{n-m} / (m*mu1 + (n-m)*mu2)
Xmax = 0;
for n = 1:K
    inner_sum = 0;
    for m = 0:n
        denominator = m * mu1 + (n - m) * mu2;
        if denominator > 0
            inner_sum = inner_sum + nchoosek(n, m) * (p1^m) * (p2^(n - m)) / denominator;
        end
    end
    Xmax = Xmax + ((-1)^(n + 1)) * inner_sum;
end

end
