function [R, names] = getAvgReward(self)
% GETAVGREWARD Get steady-state expected reward for all defined rewards
%
% [R, NAMES] = GETAVGREWARD(SELF)
%
% OUTPUTS:
%   R     - Vector of steady-state expected rewards [nRewards x 1]
%   names - Cell array of reward names {nRewards x 1}
%
% The steady-state expected reward is computed as the average reward rate
% from the value iteration, which converges to E[r] = pi * r where pi is
% the steady-state distribution.
%
% Example:
%   model.setReward('qlen', @(state, sn) state(2));
%   model.setReward('util', @(state, sn) min(state(2), 1));
%   solver = SolverCTMC(model);
%   [R, names] = solver.getAvgReward();
%   fprintf('Average queue length: %f\n', R(1));
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Check if results are already cached
if isfield(self.result, 'Reward') && ~isempty(self.result.Reward) && ...
        isfield(self.result.Reward, 'steadyState') && ~isempty(self.result.Reward.steadyState)
    R = self.result.Reward.steadyState;
    names = self.result.Reward.names;
else
    % Run the reward analyzer
    self.runRewardAnalyzer();
    R = self.result.Reward.steadyState;
    names = self.result.Reward.names;
end

end
