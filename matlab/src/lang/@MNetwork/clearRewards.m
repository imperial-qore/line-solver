function self = clearRewards(self)
% CLEARREWARDS Remove all reward definitions from the network
%
% SELF = CLEARREWARDS(SELF)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isempty(self.sn)
    self.sn.reward = {};
end

end
