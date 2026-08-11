function P = initRoutingMatrix(self)
% P = INITROUTINGMATRIX()
if ~isempty(self.obj)
    P = self.obj.initRoutingMatrix();
    return
end
M = self.getNumberOfNodes;
K = self.getNumberOfClasses;
P = RoutingMatrix(cellzeros(K,K,M,M));
end