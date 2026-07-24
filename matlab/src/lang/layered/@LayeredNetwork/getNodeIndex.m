function idx = getNodeIndex(self,node)
% IDX = GETNODEINDEX(SELF,NODE)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

idx = find(cellfun(@any,strfind(self.getNodeNames,node.name)));
end
