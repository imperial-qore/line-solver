function sinkidx = getIndexSinkNode(self)
% INDEX = GETINDEXSINKNODE()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(self.sinkidx)
    self.sinkidx = find(cellisa(self.nodes,'Sink'));
end
sinkidx = self.sinkidx;
end
