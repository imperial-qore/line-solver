%{ @file cache_erec.m
 %  @brief Recursive computation of the normalizing constant for cache models
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes the normalizing constant for cache models recursively
 %
 % @details
 % This function recursively computes the normalizing constant E for
 % cache models with given item popularity probabilities and cache capacity.
 %
 % @par Syntax:
 % @code
 % E = cache_erec(gamma, m)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities
 % <tr><td>m<td>Cache capacity
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>E<td>Normalizing constant
 % </table>
%}
function E=cache_erec(gamma,m)
E = sub_cache_erec(gamma,m,length(gamma));
end

function E=sub_cache_erec(gamma,m,k)
h=length(m);
if sum(m)==0
    E=1;
    return
end
if sum(m)>k || min(m)<0
    E=0;
    return
end

if k==1 && sum(m)==1
    j = find(m);
    E=gamma(1,j);
    return
end
E = sub_cache_erec(gamma,m,k-1);
for j=1:h
    if m(j)>0
        E = E + gamma(k,j)*m(j)*sub_cache_erec(gamma,oner(m,j),k-1);
    end
end
end
