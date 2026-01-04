%{ @file cache_prob_erec.m
 %  @brief Computes cache hit probabilities using recursive method
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes cache hit probabilities recursively
 %
 % @details
 % This function computes cache hit probability distribution using a
 % recursive method based on normalizing constants.
 %
 % @par Syntax:
 % @code
 % prob = cache_prob_erec(gamma, m)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities
 % <tr><td>m<td>Cache capacity vector
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>prob<td>Cache hit probability distribution
 % </table>
%}
function prob = cache_prob_erec(gamma,m)
[n,h]=size(gamma);
E = cache_erec(gamma, m);
for i=1:n
    for j=1:h
        Ei = cache_erec(gamma(setdiff(1:n,i),:),oner(m,j));
        prob(i,1+j) = m(j) * gamma(i,j) * Ei / E;
    end
    prob(i,1) = abs(1 - sum(prob(i,2:end)));
end
end
