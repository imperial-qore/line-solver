%{ @file cache_miss_spm.m
 %  @brief Computes miss rates using saddle-point method
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes cache miss rates using saddle-point approximation
 %
 % @details
 % This function computes global, per-user, and per-item miss rates for
 % cache models using the saddle-point method (SPM).
 %
 % @par Syntax:
 % @code
 % [M, MU, MI, pi0, lE] = cache_miss_spm(gamma, m, lambda)
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
 % <tr><td>lE<td>Log of normalizing constant
 % </table>
%}
function [M,MU,MI,pi0,lE]=cache_miss_spm(gamma,m,lambda)
% M: global miss rate
% MU: per-user miss rate
% MI: per-item miss rate
ma=m; ma(1)=ma(1)+1;
[~,lE,xi] = cache_spm(gamma,m);
[~,lEa,~] = cache_spm(gamma,ma);
M =  exp(lEa - lE);
if nargin>2
    u=size(lambda,1);
    n=size(lambda,2);
    %h=size(lambda,3)-1;
    pi0=zeros(1,n);
    % compute MU
    MU=zeros(u,1);
    if nargout>1
        for k=1:n
            if sum(gamma(k,:))>0
                [~,lE1(k)]=cache_spm(gamma(setdiff(1:n,k),:),m,xi);
                pi0(k) = exp(lE1(k) - lE);
                if pi0(k)>1 || pi0(k)<0 % try to recompute xi
                    [~,lE1(k)]=cache_spm(gamma(setdiff(1:n,k),:),m);
                    pi0(k) = exp(lE1(k) - lE);
                end
                for v=1:u
                    MU(v) = MU(v) + lambda(v,k,1)*pi0(k);
                end
            end
        end
    end
    % compute MI
    if nargout>2
        MI=zeros(n,1);
        for k=1:n
            if sum(gamma(k,:))>0
                MI(k) = MI(k) + sum(lambda(:,k,1))*exp(lE1(k)-lE);
            else
                MI(k)=0;
            end
        end
    end
end
end
