%{ @file cache_spm.m
 %  @brief Computes normalizing constant using saddle-point method
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes normalizing constant using saddle-point approximation
 %
 % @details
 % This function computes the normalizing constant for cache models using
 % the saddle-point method (SPM) for approximation.
 %
 % @par Syntax:
 % @code
 % [Z, lZ, xi] = cache_spm(gamma, m)
 % [Z, lZ, xi] = cache_spm(gamma, m, xi0)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities
 % <tr><td>m<td>Cache capacity vector
 % <tr><td>xi0<td>(Optional) Initial guess for Lagrange multipliers
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Z<td>Normalizing constant
 % <tr><td>lZ<td>Log of normalizing constant
 % <tr><td>xi<td>Lagrange multipliers
 % </table>
%}
function [Z,lZ,xi]=cache_spm(gamma,m,xi0)
gamma=gamma(find(sum(gamma,2)>0),:); %#ok<FNDSB>
h=length(m);
n=length(gamma);
mt=sum(m);
if n==mt
    line_warning(mfilename,'The number of items equals the cache capacity.\n');
    Z=cache_erec(gamma,m);
    lZ=log(Z);
    xi = cache_xi_iter(gamma,m);
    return
end

if nargin<3
    xi = cache_xi_iter(gamma,m);
else
    xi = cache_xi_iter(gamma,m,xi0);
end

for k=1:n
    S(k) = 0;
    for l=1:h
        S(k) = S(k) + gamma(k,l) * xi(l);
    end
end

%% phi
phi = 0;
for k=1:n
    phi = phi + log(1+S(k)) ;
end
phi = phi - log(xi) * m';

%% A
delta=eye(h);
for j=1:h
    for l=1:h
        C1=0;
        for k=1:n
            C1=C1+gamma(k,j)/(1+S(k));
        end
        C2=0;
        for k=1:n
            C2=C2+gamma(k,j)*gamma(k,l)/(1+S(k))^2;
        end
        C(j,l) = delta(j,l) * C1 - xi(j) * C2;
    end
end

%%
Z = exp(phi) * sqrt(2*pi)^(-h) * prod(factorial(m)) / prod(sqrt(xi)) / sqrt(det(C));
lZ = (-h) * log(sqrt(2*pi)) + (phi) + sum(factln((m))) - sum(log(sqrt(xi))) - log(sqrt(det(C)));
lZ=real(lZ); % remove small imaginary part roundoffs

end
