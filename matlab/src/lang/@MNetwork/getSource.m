function node = getSource(self)
% NODE = GETSOURCE()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

idx = self.getIndexSourceNode;
if isempty(idx)
    %    line_warning(mfilename,'The model does not have a Source station.');
    node = [];
    return
end
node = self.nodes{idx};
end
