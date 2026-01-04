%{ @file cache_t_lrum.m
 %  @brief Computes characteristic times for LRU(m) cache levels
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes characteristic times for LRU(m) cache
 %
 % @details
 % This function computes the characteristic time for each level of an
 % LRU(m) cache using fixed-point iteration.
 %
 % @par Syntax:
 % @code
 % t = cache_t_lrum(gamma, m)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities or arrival rates
 % <tr><td>m<td>Cache capacity vector
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>t<td>Characteristic time for each cache level
 % </table>
%}
function t = cache_t_lrum(gamma,m)
% fsolve solution to to T1,T2,...,Th, i.e., the characteristic time of
% each list
[n,h]=size(gamma);
fun = @time;
x = ones(1,h);
options = optimoptions('fsolve','MaxIter',1e5,'MaxFunEvals',1e6,'Display','off');
t = fsolve(fun,x,options);

function F = time(x)
F = zeros(1,h);

% the probability of each item at each list
trans = zeros(n,h);
logtrans = zeros(n,h);
denom = zeros(1,n);
capa = zeros(1,h);
stablecapa = zeros(1,h);
for k =1:n
    for j = 1:h
        trans(k,j) = exp(gamma(k,j)*x(j))-1;  %birth and death
        logtrans(k,j) = log(trans(k,j));
    end
    denom(k) = sum(exp(cumsum(logtrans(k,:))));
end
for l = 1:h
    for k = 1:n
        capa(l) = capa(l) + prod(trans(k,1:l))/(1+denom(k));
        stablecapa(l) = stablecapa(l) + exp(sum(log(trans(k,1:l)))-log(1+denom(k)));
    end
    F(l) = m(l)-stablecapa(l);
end
end
end
