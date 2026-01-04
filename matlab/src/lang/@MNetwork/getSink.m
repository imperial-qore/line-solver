function node = getSink(self)
% NODE = GETSINK()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

idx = self.getIndexSinkNode;
if isempty(idx)
    %    line_warning(mfilename,'The model does not have a Sink station.');
    node = [];
    return;
end
node = self.nodes{idx};
end
