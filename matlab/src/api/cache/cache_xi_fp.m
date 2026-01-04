%{ @file cache_xi_fp.m
 %  @brief Fixed-point iteration for computing cache performance metrics
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes cache performance metrics using fixed-point iteration
 %
 % @details
 % This function uses fixed-point iteration to compute cache performance
 % metrics including Lagrange multipliers, miss probabilities, and hit
 % probabilities for given item popularity and cache capacity.
 %
 % @par Syntax:
 % @code
 % [xi, pi0, pij, it] = cache_xi_fp(gamma, m)
 % [xi, pi0, pij, it] = cache_xi_fp(gamma, m, xi)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities
 % <tr><td>m<td>Cache capacity
 % <tr><td>xi<td>(Optional) Initial guess for Lagrange multipliers
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>xi<td>Converged Lagrange multipliers
 % <tr><td>pi0<td>Miss probability per item
 % <tr><td>pij<td>Hit probability per item per list
 % <tr><td>it<td>Number of iterations
 % </table>
%}
function [xi,pi0,pij,it] = cache_xi_fp(gamma,m,xi)
[n,h]=size(gamma);
tol=1e-14;
pi0=ones(1,n)/(h+1);
pij=zeros(n,h);
xi=zeros(1,h);
if nargin<3
    for l=1:h
        xi(l) = m(l)/mean(gamma(:,l))/(n+sum(m)-1);
    end
end
for it=1:1e4
    pi0_1=pi0;
    xi = m ./ (pi0_1*gamma);
    pij = abs(gamma .* repmat(xi,n,1)) ./ abs(1+gamma*xi');
    pi0=max(tol,1-sum(pij,2)');
    DELTA=norm(abs(1-pi0(:)./pi0_1(:)),1);
    if DELTA<tol
        xi(xi<0)=tol;
        return
    end
end
xi(xi<0)=tol;
end
