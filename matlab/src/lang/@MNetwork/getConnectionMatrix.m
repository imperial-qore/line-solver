function connections = getConnectionMatrix(self)
% [CONNMATRIX] = GETCONNECTIONMATRIX()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% connection matrix entries are setup by the addLink method
connections = self.connections;
if size(connections,1)< self.getNumberOfNodes
    connections(self.getNumberOfNodes,1)=0;
end
if size(connections,2)< self.getNumberOfNodes
    connections(1,self.getNumberOfNodes)=0;
end
end
