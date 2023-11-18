function P = initRoutingMatrix(self)
% P = INITROUTINGMATRIX()

M = self.getNumberOfNodes;
K = self.getNumberOfClasses;
P = cellzeros(K,K,M,M);
end