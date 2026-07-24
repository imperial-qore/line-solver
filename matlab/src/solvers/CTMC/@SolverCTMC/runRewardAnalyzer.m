function [V, t, names, stateSpace] = runRewardAnalyzer(self)
% RUNREWARDANALYZER Compute steady-state and transient rewards
%
% [V, T, NAMES, STATESPACE] = RUNREWARDANALYZER(SELF)
%
% Computes both:
% - Steady-state expected rewards using the CTMC stationary distribution:
%     E[r] = sum_s pi(s) * r(s)
% - Transient value functions using uniformization with value iteration
%
% Results are cached in self.result.Reward for subsequent calls.
%
% OUTPUTS:
%   V          - Cell array of value functions {nRewards x 1}
%                Each V{r} is [Tmax+1 x nStates] matrix
%   t          - Time vector [1 x Tmax+1]
%   names      - Cell array of reward names {nRewards x 1}
%   stateSpace - State space matrix [nStates x nDims]
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Get network structure with rewards
sn = self.model.getStruct(true);

% Validate that rewards are defined
if isempty(sn.reward)
    error('No rewards defined. Use model.setReward(name, @(state) ...) before calling reward analysis.');
end

% Run the reward solver
T0 = tic;
[steadyState, names, stateSpace, pi, V, t] = solver_ctmc_reward(sn, self.options);
runtime = toc(T0);

% Cache results
self.result.Reward.V = V;
self.result.Reward.t = t;
self.result.Reward.names = names;
self.result.Reward.stateSpace = stateSpace;
self.result.Reward.steadyState = steadyState;
self.result.Reward.pi = pi;
self.result.Reward.runtime = runtime;

end
