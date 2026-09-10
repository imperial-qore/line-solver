%{ @file cache_miss_fpi.m
 %  @brief Computes miss rates using fixed-point iteration
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes cache miss rates using FPI method
 %
 % @details
 % This function computes global, per-user, and per-item miss rates for
 % cache models using the fixed-point iteration (FPI) method.
 %
 % @par Syntax:
 % @code
 % [M, MU, MI, pi0] = cache_miss_fpi(gamma, m, lambda)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities
 % <tr><td>m<td>Cache capacity vector
 % <tr><td>lambda<td>(Optional) Arrival rates per user per item
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>M<td>Global miss rate
 % <tr><td>MU<td>Per-user miss rate
 % <tr><td>MI<td>Per-item miss rate
 % <tr><td>pi0<td>Per-item miss probability
 % </table>
%}
function [M,MU,MI,pi0]=cache_miss_fpi(gamma,m,lambda)
% FPI method
[n,~]=size(gamma);
xi = cache_xi_fp(gamma,m);
MI=zeros(n,1);
for i=1:n
    MI(i) = sum(lambda(:,i,1))/(1+gamma(i,:)*xi(:));
end
M=sum(MI);
if nargin>2
    u=size(lambda,1);
    n=size(lambda,2);
    MU=zeros(u,1);
    pi0=zeros(1,n);
    for i=1:n
        pi0(i) = 1/(1+gamma(i,:)*xi(:));
        for v=1:u
            MU(v) = MU(v) + lambda(v,i,1)*pi0(i);
        end
    end
end
end
