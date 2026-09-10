%{ @file cache_mva.m
 %  @brief Exact Mean Value Analysis for caches
 %
 %  @author LINE Development Team
%}

%{
 % @brief Performs exact Mean Value Analysis for cache systems
 %
 % @details
 % This function computes exact performance metrics for cache systems
 % using Mean Value Analysis with given item popularity probabilities
 % or arrival rates and cache capacity.
 %
 % @par Syntax:
 % @code
 % [pi, pi0, pij, x, u, E] = cache_mva(gamma, m)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities or arrival rates
 % <tr><td>m<td>Cache capacity vector (size h for h levels)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>pi<td>State probability vector
 % <tr><td>pi0<td>Miss probability per item
 % <tr><td>pij<td>Hit probability per item per list
 % <tr><td>x<td>Throughput vector
 % <tr><td>u<td>Utilization per item
 % <tr><td>E<td>Normalizing constant
 % </table>
%}
function [pi,pi0,pij,x,u,E] = cache_mva(gamma,m)
[n,h]=size(gamma);
SS=[];
for l=1:h
    SS = State.cartesian(SS, [1:(m(l)+1)]');
end
SS=SS-1;
pi=zeros(size(SS,1),n);
pij=zeros(size(SS,1),n,h);
x=zeros(1,h);
E=1;
%Ecur=SS(1,:);
for s=1:size(SS,1)
    mcur = SS(s,:);
    for l=1:h
        mcur_l = oner(mcur,l);
        s_l = matchrow(SS,mcur_l);
        if s_l > 0
            x(l) = mcur(l)/(gamma(:,l)'*(1-pi(s_l,:))');
            pij(s,:,l) = gamma(:,l)'.*(1-pi(s_l,:))*x(l);
            pi(s,:) = pi(s,:) + pij(s,:,l);
        end
    end
end
s = matchrow(SS,m);
pi=pi(s,:)';
pij=reshape(pij(s,:,:),n,h);
pi0=1-pi;
if nargout>2
     for l=1:h
         for k=1:n
             u(k,l)=x(l)*gamma(k,l);
         end
     end
end
end
