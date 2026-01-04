function fcr = addRegion(self, nodes)
% ADDREGION(NODES)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%
fcr = Region(nodes, self.classes);
regionIndex = length(self.regions) + 1;
fcr.setName(sprintf('FCR%d', regionIndex));
self.regions{end+1,1} = fcr;
end
