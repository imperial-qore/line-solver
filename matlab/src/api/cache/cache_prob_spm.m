%{ @file cache_prob_spm.m
 %  @brief Computes cache hit probabilities using saddle-point method
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes cache hit probabilities using saddle-point approximation
 %
 % @details
 % This function computes cache hit probability distribution using the
 % saddle-point method (SPM) for approximation.
 %
 % @par Syntax:
 % @code
 % prob = cache_prob_spm(gamma, m)
 % prob = cache_prob_spm(gamma, m, lE)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities
 % <tr><td>m<td>Cache capacity vector
 % <tr><td>lE<td>(Optional) Pre-computed log of normalizing constant
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>prob<td>Cache hit probability distribution
 % </table>
%}
function prob = cache_prob_spm(gamma,m, lE)
[n,h]=size(gamma);
if nargin < 3
[~, lE] = cache_spm(gamma,m);
end
for i=1:n
    for j=1:h
        [~, lEi] = cache_spm(gamma(setdiff(1:n,i),:),oner(m,j));
        prob(i,1+j) = m(j) * gamma(i,j) * exp(lEi-lE);
    end
    prob(i,1) = abs(1 - sum(prob(i,2:end)));
end
end
