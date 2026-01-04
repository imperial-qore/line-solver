function [t, V, names, stateSpace] = getReward(self, rewardName)
% GETREWARD Get reward value function over time from CTMC analysis
%
% [T, V, NAMES, STATESPACE] = GETREWARD(SELF) returns all rewards
%
% [T, V, NAMES, STATESPACE] = GETREWARD(SELF, REWARDNAME) returns the
% specific reward with the given name
%
% OUTPUTS:
%   t          - Time vector [1 x Tmax+1]
%   V          - Value functions: cell array {nRewards x 1} where each
%                V{r} is [Tmax+1 x nStates] matrix, or single matrix if
%                rewardName is specified
%   names      - Cell array of reward names, or single name string
%   stateSpace - State space matrix [nStates x nDims]
%
% The value function V{r}(k,s) represents the cumulative reward accumulated
% starting from state s after k time steps (uniformized).
%
% Example:
%   model.setReward('qlen', @(state, sn) state(2));
%   solver = SolverCTMC(model);
%   [t, V, names] = solver.getReward();
%   plot(t, mean(V{1}, 2));  % Plot mean reward over time
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    rewardName = [];
end

% Check if results are already cached
if isfield(self.result, 'Reward') && ~isempty(self.result.Reward) && ...
        isfield(self.result.Reward, 'V') && ~isempty(self.result.Reward.V)
    V = self.result.Reward.V;
    t = self.result.Reward.t;
    names = self.result.Reward.names;
    stateSpace = self.result.Reward.stateSpace;
else
    % Run the reward analyzer
    [V, t, names, stateSpace] = self.runRewardAnalyzer();
end

% Filter to specific reward if requested
if ~isempty(rewardName)
    rewardName = char(rewardName);
    idx = find(strcmp(names, rewardName));
    if isempty(idx)
        error('Reward "%s" not found. Available rewards: %s', ...
            rewardName, strjoin(names, ', '));
    end
    V = V{idx};
    names = names{idx};
end

end
