%{ @file cache_xi_iter.m
 %  @brief Computes Lagrange multipliers using Gast-van Houdt fixed-point iteration
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes Lagrange multipliers using Gast-van Houdt fixed-point iteration
 %
 % @details
 % This function computes Lagrange multipliers for cache models using
 % the Gast-van Houdt algorithm fixed-point iteration method.
 %
 % @par Syntax:
 % @code
 % z = cache_xi_iter(gamma, m)
 % z = cache_xi_iter(gamma, m, tmax)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities
 % <tr><td>m<td>Cache capacity
 % <tr><td>tmax<td>(Optional) Maximum time limit
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>z<td>Converged Lagrange multipliers
 % </table>
%}
function z=cache_xi_iter(gamma,m,tmax)
if nargin<3
    tmax=Inf;
end

n=size(gamma,1);
f=m./n;
h=size(f,2);

pp=zeros(h+1,n);
pp(1,:)=ones(1,n);
for i=1:h
     pp(i+1,:)=gamma(:,i);
end

z_old=zeros(1,h+1);
z=ones(1,h+1);

while (max(abs(z-z_old)) > 10^(-12)*max(abs(z_old)))
    z_old=z;
    temp=n*z*pp;
    for i=1:h % update z(i+1)
        a=temp-n*z(i+1)*pp(i+1,:);
        Fi=sum(pp(i+1,:)./(n*pp(i+1,:)+a)); %Fi(1)
        if (Fi > f(i))
            zi_min=0;
            zi_max=1;
        else
            zi_min=1;
            zi_max=2;
            while (sum(zi_max*pp(i+1,:)./(n*zi_max*pp(i+1,:)+a)) < f(i)) %Fi(zi_max)
                zi_min=zi_max;
                zi_max=zi_max*2;
            end
        end
        for x=1:50
            z(i+1)=(zi_min+zi_max)/2;
            if (sum(z(i+1)*pp(i+1,:)./(n*z(i+1)*pp(i+1,:)+a)) < f(i)) %Fi(zi_max)
                zi_min=z(i+1);
            else
                zi_max=z(i+1);
            end
        end
    end
end

z=z(2:end);
end
