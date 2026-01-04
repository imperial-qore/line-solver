function addLink(self, sourceNode, destNode)
% ADDLINK(SOURCENODE, DESTNODE)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
self.connections(sourceNode.index, destNode.index) = 1;
end
