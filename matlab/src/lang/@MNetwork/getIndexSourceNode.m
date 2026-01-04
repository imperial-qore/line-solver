function index = getIndexSourceNode(self)
% INDEX = GETINDEXSOURCENODE()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

index = find(cellisa(self.nodes,'Source'));
end
