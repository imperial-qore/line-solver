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
n=size(gamma,1); % item count is the row count; length() is max(n,h) and reads past gamma when h>n
mt=sum(m);
if n==mt
    % Degenerate saddle: every item is cached, so the capacity equations force
    % every multiplier to infinity and cache_xi_iter cannot converge. Take Z
    % from the exact recursion and report the limit, rather than iterating.
    line_warning(mfilename,'The number of items equals the cache capacity.\n');
    Z=cache_erec(gamma,m);
    lZ=log(Z);
    xi=inf(1,h);
    return
end

% A list with no capacity has xi=0, which is a boundary of the Laplace integral
% rather than a direction of it, so it must leave the expansion: kept, its
% -sum_l log(sqrt(xi_l)) prefactor diverges and Z comes out far too large.
% Dropping it is exact, since setting z_l=0 in the generating function removes
% list l from E(m) and prod_l m_l! is unchanged because 0!=1.
keep=find(m>0);
hk=length(keep);
xi=zeros(1,h);
if hk==0
    Z=1; lZ=0; % E(0)=1 and prod_l m_l!=1
    return
end
gk=gamma(:,keep);
mk=m(keep);

if nargin<3
    xik = cache_xi_iter(gk,mk);
else
    xik = cache_xi_iter(gk,mk,xi0);
end
xi(keep)=xik;

S = zeros(1, n);
for k=1:n
    for l=1:hk
        S(k) = S(k) + gk(k,l) * xik(l);
    end
end

%% phi
phi = 0;
for k=1:n
    phi = phi + log(1+S(k)) ;
end
phi = phi - log(xik) * mk';

%% A
delta=eye(hk);
C = zeros(hk, hk);
for j=1:hk
    for l=1:hk
        C1=0;
        for k=1:n
            C1=C1+gk(k,j)/(1+S(k));
        end
        C2=0;
        for k=1:n
            C2=C2+gk(k,j)*gk(k,l)/(1+S(k))^2;
        end
        C(j,l) = delta(j,l) * C1 - xik(j) * C2;
    end
end

%%
Z = exp(phi) * sqrt(2*pi)^(-hk) * prod(factorial(mk)) / prod(sqrt(xik)) / sqrt(det(C));
lZ = (-hk) * log(sqrt(2*pi)) + (phi) + sum(factln((mk))) - sum(log(sqrt(xik))) - log(sqrt(det(C)));
lZ=real(lZ); % remove small imaginary part roundoffs

end
