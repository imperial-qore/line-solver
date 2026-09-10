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
 % prob = cache_prob_erec(gamma, m, sigma, k)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities
 % <tr><td>m<td>Cache capacity vector
 % <tr><td>sigma<td>(Optional) item storage costs (sizes), 1 x n
 % <tr><td>k<td>(Optional) per-list storage cost caps, 1 x h
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>prob<td>Cache hit probability distribution
 % </table>
%}
function prob = cache_prob_erec(gamma,m,sigma,k)
[n,h]=size(gamma);
if nargin<3 || isempty(sigma) || isempty(k)
    sigma = []; k = [];
else
    sigma = sigma(:).'; k = k(:).';
end
E = cache_erec(gamma, m, sigma, k);
prob = zeros(n,h+1);
for i=1:n
    others = setdiff(1:n,i);
    for j=1:h
        if isempty(sigma)
            Ei = cache_erec(gamma(others,:),oner(m,j));
        else
            kij = k; kij(j) = kij(j) - sigma(i);
            if kij(j) < 0
                Ei = 0;
            else
                Ei = cache_erec(gamma(others,:),oner(m,j),sigma(others),kij);
            end
        end
        prob(i,1+j) = m(j) * gamma(i,j) * Ei / E;
    end
    prob(i,1) = abs(1 - sum(prob(i,2:end)));
end
end
