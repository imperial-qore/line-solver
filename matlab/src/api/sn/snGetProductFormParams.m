function [lambda,D,N,Z,mu,S,V]= snGetProductFormParams(sn)
% [LAMBDA,D,N,Z,MU,S,V]= SNGETPRODUCTFORMPARAMS()

% mu also returns max(S) elements after population |N| as this is
% required by MVALDMX

R = sn.nclasses;
N = sn.njobs;
queueIndices = find(sn.nodetype == NodeType.Queue);
delayIndices = find(sn.nodetype == NodeType.Delay);
sourceIndex = find(sn.nodetype == NodeType.Source);
Mq = length(queueIndices); % number of queues
Mz = length(delayIndices); % number of delays
lambda = zeros(1,R);
S = sn.nservers(sn.nodeToStation(queueIndices));
for r=1:R
    if isinf(N(r))
        lambda(r) = sn.rates(sn.nodeToStation(sourceIndex),r);
    end
end

D = zeros(Mq,R);
Nct = sum(N(isfinite(N)));
mu = ones(Mq, ceil(Nct)+max(S(isfinite(S))));
for i=1:Mq
    for r=1:R
        c = find(sn.chains(:,r),1);
        if sn.refclass(c)>0
            D(i,r) = sn.visits{c}(sn.nodeToStateful(queueIndices(i)),r) / sn.rates(sn.nodeToStation(queueIndices(i)),r) / sn.visits{c}(sn.stationToStateful(sn.refstat(r)),sn.refclass(c));
        else
            D(i,r) = sn.visits{c}(sn.nodeToStateful(queueIndices(i)),r) / sn.rates(sn.nodeToStation(queueIndices(i)),r);
        end
    end
    mu(i,1:size(mu,2)) = min(1:size(mu,2), sn.nservers(sn.nodeToStation(queueIndices(i))));
end
Z = zeros(max(1,Mz),R);
for i=1:Mz
    for r=1:R
        c = find(sn.chains(:,r),1);
        if sn.refclass(c)>0
            Z(r) = sn.visits{c}(sn.nodeToStateful(delayIndices(i)),r) / sn.rates(sn.nodeToStation(delayIndices(i)),r) / sn.visits{c}(sn.stationToStateful(sn.refstat(r)),sn.refclass(c));
        else
            Z(r) = sn.visits{c}(sn.nodeToStateful(delayIndices(i)),r) / sn.rates(sn.nodeToStation(delayIndices(i)),r);
        end
    end
end
V=cellsum(sn.visits);
end
