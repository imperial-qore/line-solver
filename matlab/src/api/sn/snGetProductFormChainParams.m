function [lambda,D,N,Z,mu,S,V]= snGetProductFormChainParams(sn)
% [LAMBDA,D,N,Z,MU,S,V]= NETWORKSTRUCT_GETPRODUCTFORMCHAINPARAMS(SN)

% mu also returns max(S) elements after population |N| as this is
% required by MVALDMX

[lambda,~,~,~,mu,~] = snGetProductFormParams(sn);
queueIndex = find(sn.nodetype == NodeType.Queue);
delayIndex = find(sn.nodetype == NodeType.Delay);
ignoreIndex = find(sn.nodetype == NodeType.Source | sn.nodetype == NodeType.Join);
[Dchain,~,Vchain,~,Nchain,~,~] = snGetDemandsChain(sn);
lambda_chains = zeros(1,sn.nchains);
for c=1:sn.nchains
    lambda_chains(c) = nansum(lambda(sn.inchain{c}));
    %if sn.refclass(c)>0
        %D_chains(:,c) = Dchain(find(isfinite(sn.nservers)),c)/alpha(sn.refstat(c),sn.refclass(c));
        %Z_chains(:,c) = Dchain(find(isinf(sn.nservers)),c)/alpha(sn.refstat(c),sn.refclass(c));
    %else
        D_chains(:,c) = Dchain(sn.nodeToStation(queueIndex),c);
        Z_chains(:,c) = Dchain(sn.nodeToStation(delayIndex),c);
    %end
end
S = sn.nservers(sn.nodeToStation(queueIndex));
lambda = lambda_chains;
N = Nchain;
%D_chains(sn.nodeToStation(ignoreIndex),:) =[];
Vchain(sn.nodeToStation(ignoreIndex),:) =[];
D = D_chains;
Z = Z_chains;

V = Vchain;
if isempty(Z)
    Z = 0*N;
end
end
    
    
