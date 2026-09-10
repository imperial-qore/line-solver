%{ @file fj_char_max_blom.m
 %  @brief Blom-corrected plotting position for the characteristic maximum
 %
 %  @author LINE Development Team
%}

%{
 % @brief Blom-corrected plotting position for the characteristic maximum
 %
 % @details
 % The characteristic maximum m_K is the quantile at which the survival
 % function of a branch drops to 1/K. The naive plotting position
 % m_K = F^-1(K/(K+1)) is biased; Blom's correction replaces it with
 %
 %   m_K = F^-1( (K - alpha) / (K - alpha - beta + 1) ),
 %
 % which for alpha = beta = 0 falls back on the naive position. The survey
 % quotes alpha = 0.4886 and beta = 0.3140, which are the defaults here.
 %
 % For the standard normal branch the position is bracketed without any
 % inversion, for K >= 5, by
 %
 %   sqrt(2 ln K - ln ln K - 3) < m_K < sqrt(2 ln K - ln ln K),
 %
 % and the leading term alone gives the Kruskal-Weiss estimate
 % m_K ~ mu + sigma sqrt(2 ln K).
 %
 % @par Syntax:
 % @code
 % mK = fj_char_max_blom(K)
 % [mK, lo, hi] = fj_char_max_blom(K)
 % [mK, lo, hi] = fj_char_max_blom(K, Finv, alpha, beta)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>K<td>Number of i.i.d. copies (positive integer)
 % <tr><td>Finv<td>Quantile function handle (optional; the standard normal by default)
 % <tr><td>alpha<td>Blom numerator offset (optional, default 0.4886)
 % <tr><td>beta<td>Blom denominator offset (optional, default 0.3140)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>mK<td>Blom-corrected characteristic maximum
 % <tr><td>lo<td>Kruskal-Weiss lower bracket, valid for the standard normal and K >= 5
 % <tr><td>hi<td>Kruskal-Weiss upper bracket, valid for the standard normal and K >= 5
 % </table>
 %
 % @par Reference:
 % A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
 % ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014, Section 4.9
 % on page 17:26.
 %
 % Original: G. Blom, "Statistical Estimates and Transformed Beta Variables",
 % Wiley, 1958; W. Kruskal, G. Weiss, "Allocation of Observations in
 % Comparison of Treatments", 1985.
%}
function [mK, lo, hi] = fj_char_max_blom(K, Finv, alpha, beta)

if K < 1 || K ~= round(K)
    line_error(mfilename, 'K must be a positive integer. Got K=%g.', K);
end
if nargin < 3 || isempty(alpha)
    alpha = 0.4886;
end
if nargin < 4 || isempty(beta)
    beta = 0.3140;
end

den = K - alpha - beta + 1;
if den <= 0
    line_error(mfilename, 'The Blom offsets leave a non-positive denominator at K=%d.', K);
end
q = (K - alpha) / den;
if q <= 0 || q >= 1
    line_error(mfilename, 'The Blom plotting position %g fell outside (0,1).', q);
end

if nargin < 2 || isempty(Finv)
    % Standard normal quantile through the inverse error function
    mK = sqrt(2) * erfinv(2 * q - 1);
else
    mK = Finv(q);
end

% Kruskal-Weiss bracket for the standard normal, meaningful from K = 5 on
lo = NaN;
hi = NaN;
if K >= 5
    z = 2 * log(K) - log(log(K));
    if z > 3
        lo = sqrt(z - 3);
    end
    if z > 0
        hi = sqrt(z);
    end
end

end
