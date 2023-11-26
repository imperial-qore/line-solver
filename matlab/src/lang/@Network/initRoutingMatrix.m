function P = initRoutingMatrix(self)
% P = INITROUTINGMATRIX()

M = self.getNumberOfNodes;
K = self.getNumberOfClasses;
P = RoutingMatrix(cellzeros(K,K,M,M));
end