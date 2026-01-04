%{ @file cache_t_hlru.m
 %  @brief Computes characteristic times for hierarchical LRU cache levels
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes characteristic times for hierarchical LRU cache
 %
 % @details
 % This function computes the characteristic time for each level of a
 % hierarchical LRU cache using fixed-point iteration.
 %
 % @par Syntax:
 % @code
 % t = cache_t_hlru(gamma, m)
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
function t = cache_t_hlru(gamma,m)
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
for a = 1:h
    temp1 = ones(n,1);
    temp2 = zeros(n,1);
    probh = zeros(n,1);
for k = 1:n
    for s = 1:a
        temp1(k) = temp1(k)*(1-exp(-gamma(k,s)*x(s)));
    end

    middtemp = 1;
    middtemp2 = 0;
    for l = 1:a-1
        for s = 1:l
            middtemp = middtemp*(1-exp(-gamma(k,s)*x(s)));
        end
        middtemp2 = middtemp2+middtemp;
    end
    temp2(k) = exp(-gamma(k,a)*x(a))*(1+middtemp2);

    probh(k) = temp1(k)/(temp1(k)+temp2(k));
end

F(a) = m(a)-sum(probh);
end
end

end
