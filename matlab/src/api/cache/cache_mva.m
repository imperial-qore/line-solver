function [pi,pi0,pij,x,u,E] = cache_mva(gamma,m)
[n,h]=size(gamma);
SS=[];
for l=1:h
    SS = ssg_decorate(SS, [1:(m(l)+1)]');
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

function SS = ssg_decorate(SS, SS2)
% SS = ssg_decorate(SS1, SS2)
% INPUT:
% SS1 : a state space (n1,...,nk)
% SS2 : a state space (m1,...,mk)
% OUTPUT:
% SS  : a state space (n1,...,nk,m1,...,mk)
% EXAMPLE:
% SS = ssg_closed_single(3,2)
% RR=ssg_renv(10)
% ssg_decorate(SS,RR)
%
% Copyright (c) 2012-2014, Imperial College London 
% All rights reserved.


if isempty(SS)
    SS = SS2;
    return
else if isempty(SS2)
    return
end

n1 = size(SS,1); m1 = size(SS,2);
n2 = size(SS2,1); m2 = size(SS2,2);
SS = repmat(SS, n2, 1);

curStates = 1:n1;
for s=1:n2    
    SS(curStates,(m1+1):(m1+m2)) = repmat(SS2(s,:),length(curStates),1);
    curStates = (curStates) + n1;
end

end
end