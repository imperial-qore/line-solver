function rewards = getReward(self, name)
% GETREWARD Get reward definition(s) from the network
%
% REWARDS = GETREWARD(SELF) returns all reward definitions as a cell array
%
% REWARDS = GETREWARD(SELF, NAME) returns the reward definition with the
% specified name, or empty if not found
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(self.sn) || isempty(self.sn.reward)
    rewards = {};
    return;
end

if nargin < 2 || isempty(name)
    % Return all rewards
    rewards = self.sn.reward;
else
    % Find specific reward by name
    name = char(name);
    rewards = {};
    for i = 1:length(self.sn.reward)
        if strcmp(self.sn.reward{i}.name, name)
            rewards = self.sn.reward{i};
            return;
        end
    end
end

end
