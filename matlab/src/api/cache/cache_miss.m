%{ @file cache_miss.m
 %  @brief Cache miss rate computation
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes miss rates for a cache system
 %
 % @details
 % This function computes various miss rate metrics for a cache system
 % with given item popularity distribution and cache capacity.
 %
 % @par Syntax:
 % @code
 % [M,MU,MI,pi0] = cache_miss(gamma, m)
 % [M,MU,MI,pi0] = cache_miss(gamma, m, lambda)
 % [M,MU,MI,pi0] = cache_miss(gamma, m, lambda, sigma, k)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities
 % <tr><td>m<td>Cache capacity
 % <tr><td>lambda<td>(Optional) Arrival rates per user per item
 % <tr><td>sigma<td>(Optional) item storage costs (sizes), 1 x n
 % <tr><td>k<td>(Optional) per-list storage cost caps, 1 x h
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
 %
 % @see cache_erec
%}
function [M,MU,MI,pi0]=cache_miss(gamma,m,lambda,sigma,cap)
% M: global miss rate
% MU: per-user miss rate
% MI: per-item miss rate
% pi0: per-item miss probability
if nargin<4 || isempty(sigma) || isempty(cap)
    sigma = []; cap = [];
else
    sigma = sigma(:).'; cap = cap(:).';
end
ma=m; ma(1)=ma(1)+1;
E = cache_erec(gamma,m,sigma,cap);
M = cache_erec(gamma,ma,sigma,cap) / E;
if nargin>2
    u=size(lambda,1);
    n=size(lambda,2);
    MU=zeros(u,1);
    pi0=zeros(1,n);
    for k=1:n
        others = setdiff(1:n,k);
        if isempty(sigma)
            pi0(k) = cache_erec(gamma(others,:),m) / E;
        else
            pi0(k) = cache_erec(gamma(others,:),m,sigma(others),cap) / E;
        end
    end
    for v=1:u
        MU(v) = sum(lambda(v,:,1).*pi0);
    end
    MI=zeros(n,1);
    for k=1:n
        MI(k) = MI(k) + sum(lambda(:,k,1))*pi0(k);
    end
end
end