function self = setReward(self, name, rewardFn)
% SETREWARD Define a reward function on the network for CTMC analysis
%
% SELF = SETREWARD(SELF, NAME, REWARDFN)
%
% NAME     - String identifier for the reward
% REWARDFN - Function handle with signature:
%            @(state) -> double              (recommended - new API)
%            @(state, sn) -> double          (backward compatible)
%
%            state: RewardState object providing intuitive state access:
%              - state.at(node)              - Jobs at node (all classes)
%              - state.at(node, class)       - Jobs at node for specific class
%              - state.forClass(class)       - Jobs of class (all nodes)
%
%            RewardState also supports aggregation via RewardStateView:
%              - state.at(node).total()      - Sum jobs across classes at node
%              - state.at(node).max()        - Max jobs across classes at node
%              - state.at(node).min()        - Min jobs across classes at node
%              - state.at(node).count()      - Count non-zero entries
%
%            sn: NetworkStruct (optional, for advanced users to access
%                rates, capacities, etc.)
%
% Multiple rewards can be defined by calling setReward multiple times
% with different names. Use clearRewards() to remove all rewards.
%
% Recommended: Use Reward templates for common metrics
%   Reward.queueLength(node, [class])
%   Reward.utilization(node, [class])
%   Reward.blocking(node)
%
% Examples (New API - Recommended):
%   % Basic: reference nodes and classes directly
%   model.setReward('QLen_Q1_C1', @(state) state.at(queue1, class1));
%
%   % Aggregation: sum across classes
%   model.setReward('TotalJobs_Q1', @(state) state.at(queue1).total());
%
%   % Aggregation: sum across stations
%   model.setReward('ClassTotal_C1', @(state) state.forClass(class1).total());
%
%   % Using templates
%   model.setReward('Util_Q1', Reward.utilization(queue1));
%   model.setReward('Block_Q1', Reward.blocking(queue1));
%
%   % Custom composite rewards
%   model.setReward('Cost', @(state) ...
%       2.0 * state.at(queue1).total() + ...
%       5.0 * state.at(queue2).total());
%
%   % Advanced: using NetworkStruct parameter (for power-users)
%   model.setReward('AdvancedMetric', @(state, sn) ...
%       state.at(queue1, class1) * sn.rates(sn.nodeToStation(queue1.index), 1));
%
% See also: RewardState, RewardStateView, Reward, getAvgReward
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~ischar(name) && ~isstring(name)
    error('Reward name must be a string');
end
if ~isa(rewardFn, 'function_handle')
    error('Reward function must be a function handle @(state, sn) -> double');
end

% Initialize reward cell array if needed
if isempty(self.sn)
    self.sn = NetworkStruct();
end
if ~isfield(self.sn, 'reward') || isempty(self.sn.reward)
    self.sn.reward = {};
end

% Check if reward with this name already exists and update it
name = char(name);
found = false;
for i = 1:length(self.sn.reward)
    if strcmp(self.sn.reward{i}.name, name)
        self.sn.reward{i}.fn = rewardFn;
        found = true;
        break;
    end
end

% Add new reward if not found
if ~found
    reward.name = name;
    reward.fn = rewardFn;
    reward.type = 'state';
    self.sn.reward{end+1} = reward;
end

end
