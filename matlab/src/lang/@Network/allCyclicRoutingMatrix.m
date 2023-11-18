function P = allCyclicRoutingMatrix(self)
% P = ALLCYCLICROUTINGMATRIX()

M = self.getNumberOfNodes;
K = self.getNumberOfClasses;
P = cellzeros(K,K,M,M);
for k=1:K
    P{k} = circul(M);
end
end