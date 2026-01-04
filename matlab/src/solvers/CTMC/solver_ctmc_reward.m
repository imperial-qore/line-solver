function [steadyStateRewards, rewardNames, stateSpace, pi, V, t] = solver_ctmc_reward(sn, options)
% SOLVER_CTMC_REWARD Compute steady-state and transient rewards
%
% [STEADYSTATE, NAMES, STATESPACE, PI, V, T] = SOLVER_CTMC_REWARD(SN, OPTIONS)
%
% Computes both steady-state expected rewards using the stationary distribution:
%   E[r] = sum_s pi(s) * r(s)
%
% and transient value functions using uniformization with value iteration:
%   V^{k+1}(s) = r(s) + sum_{s'} P(s,s') * V^k(s')
%
% INPUTS:
%   sn      - NetworkStruct with .reward field populated
%   options - Solver options with:
%             .rewardIterations - Number of iterations for transient (default 1000)
%
% OUTPUTS:
%   steadyStateRewards - Vector of steady-state expected rewards [nRewards x 1]
%   rewardNames        - Cell array of reward names
%   stateSpace         - State space matrix [nStates x nDims]
%   pi                 - Stationary distribution [1 x nStates]
%   V                  - Cell array of value functions {nRewards x 1}
%                        Each V{r} is [Tmax+1 x nStates] matrix
%   t                  - Time vector (scaled by uniformization rate)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Validate rewards are defined
if isempty(sn.reward)
    error('No rewards defined. Use model.setReward() before calling this solver.');
end

if isinf(options.cutoff)
    options.cutoff = 100;
    line_printf(mfilename,'Cutoff option missing, setting cutoff=100');
end

% Get default iterations for transient analysis
if isfield(options, 'rewardIterations') && ~isempty(options.rewardIterations)
    Tmax = options.rewardIterations;
else
    Tmax = 1000;
end

% Get generator and state space
[Q, ~, stateSpaceAggr, ~, ~, ~, sn] = solver_ctmc(sn, options);
stateSpace = stateSpaceAggr;  % Return aggregated state space to user

nstates = size(Q, 1);
nRewards = length(sn.reward);

% Compute stationary distribution for steady-state rewards
pi = ctmc_solve(Q, options);
pi(pi < GlobalConstants.Zero) = 0;
pi = pi / sum(pi);  % Normalize

% Build index maps for RewardState
% nodeToStationMap: node.index -> station index
nodeToStationMap = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
classToIndexMap = containers.Map('KeyType', 'int32', 'ValueType', 'int32');

for ind = 1:sn.nnodes
    if sn.isstation(ind)
        nodeToStationMap(int32(ind)) = sn.nodeToStation(ind);
    end
end

for r = 1:sn.nclasses
    classToIndexMap(int32(r)) = r;
end

% Build reward vectors (one per reward definition)
R = zeros(nRewards, nstates);
rewardNames = cell(nRewards, 1);
for r = 1:nRewards
    rewardNames{r} = sn.reward{r}.name;
    rewardFn = sn.reward{r}.fn;

    for s = 1:nstates
        % Create RewardState for this state row
        stateRow = stateSpaceAggr(s, :);
        rewardState = RewardState(stateRow, sn, nodeToStationMap, classToIndexMap);

        % Try calling with RewardState (new API, single argument)
        try
            R(r, s) = rewardFn(rewardState);
        catch ME
            % If that fails, try backward compatibility with @(state, sn) signature
            try
                R(r, s) = rewardFn(stateSpaceAggr(s,:), sn);
            catch ME2
                % If both fail, report the original error from the new API
                rethrow(ME);
            end
        end
    end
end

% Compute steady-state expected rewards: E[r] = pi * r'
steadyStateRewards = zeros(nRewards, 1);
for r = 1:nRewards
    steadyStateRewards(r) = pi * R(r, :)';
end

% Compute transient value functions via uniformization
% Uniformization: find maximum exit rate
q = max(abs(diag(Q)));
if q == 0
    q = 1;  % Handle absorbing states edge case
end

% Build transition probability matrix P = Q/q + I
P = Q/q + speye(nstates);

% Value iteration
V = cell(nRewards, 1);
for r = 1:nRewards
    V{r} = zeros(Tmax+1, nstates);
    v_prev = zeros(1, nstates);
    for k = 1:Tmax
        % V^{k+1}(s) = r(s) + sum_{s'} P(s,s') * V^k(s')
        v_new = R(r, :) + v_prev * P';
        V{r}(k+1, :) = v_new;
        v_prev = v_new;
    end
end

% Time scaling: convert iteration index to continuous time
t = (0:Tmax) / q;

end
