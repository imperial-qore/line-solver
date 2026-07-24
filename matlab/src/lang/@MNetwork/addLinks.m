function addLinks(self, nodesList)
% ADDLINKS(NODESLIST)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
for i=1:size(nodesList,1)
    self.addLink(self.nodes{nodesList(i,1)}, self.nodes{nodesList(i,2)});
end
end
