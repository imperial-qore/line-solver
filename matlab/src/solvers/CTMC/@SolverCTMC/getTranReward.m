function [Rt, t, names] = getTranReward(self, rewardName)
% GETTRANREWARD Get transient expected reward metrics over time
%
% [RT, T, NAMES] = GETTRANREWARD(SELF) returns transient expected rewards
% for all defined rewards
%
% [RT, T, NAMES] = GETTRANREWARD(SELF, REWARDNAME) returns the transient
% expected reward for a specific named reward
%
% This method computes the transient expected reward E[r(X(t))] for each
% time point t, where X(t) is the system state at time t and r is the
% reward function. The computation uses the transient probability
% distribution pi_t(s) from CTMC uniformization:
%
%   E[r(X(t))] = sum_s pi_t(s) * r(s)
%
% OUTPUTS:
%   Rt    - Cell array {nRewards x 1} where each Rt{r} is a struct with:
%           .t      - Time vector [nTimePoints x 1]
%           .metric - Expected reward at each time [nTimePoints x 1]
%           .name   - Reward name string
%           Or a single struct if rewardName is specified
%   t     - Time vector [nTimePoints x 1]
%   names - Cell array of reward names, or single string if rewardName specified
%
% Example:
%   model = Network('example');
%   % ... model setup ...
%   model.setReward('QLen', @(state) state.at(queue, class1));
%   solver = SolverCTMC(model, 'timespan', [0, 10]);
%   [Rt, t, names] = solver.getTranReward();
%   plot(t, Rt{1}.metric);
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    rewardName = [];
end

options = self.getOptions;

% Validate timespan
if ~isfield(options, 'timespan') || ~isfinite(options.timespan(2))
    line_error(mfilename, 'getTranReward requires a finite timespan. Use SolverCTMC(model, ''timespan'', [0, T]).');
end

sn = self.model.getStruct(true);

% Validate that rewards are defined
if isempty(sn.reward)
    line_error(mfilename, 'No rewards defined. Use model.setReward(name, @(state) ...) before calling getTranReward.');
end

nRewards = length(sn.reward);

% Get transient probabilities and state space
[InfGen, StateSpace, StateSpaceAggr, ~, ~, ~, sn] = solver_ctmc(sn, options);

% Compute initial state probability
state = [];
for ist = 1:sn.nnodes
    if sn.isstateful(ist)
        isf = sn.nodeToStateful(ist);
        state = [state, zeros(1, size(sn.space{isf}, 2) - length(sn.state{isf})), sn.state{isf}];
    end
end

nstates = size(InfGen, 1);
pi0 = zeros(1, nstates);

state0 = matchrow(StateSpace, state);
if state0 == -1
    line_error(mfilename, 'Initial state not contained in the state space.');
end
pi0(state0) = 1;

% Compute transient probabilities
[pit, t] = ctmc_transient(InfGen, pi0, options.timespan(1), options.timespan(2), options.stiff, [], options.timestep);
pit(pit < GlobalConstants.Zero) = 0;

% Build index maps for RewardState
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
names = cell(nRewards, 1);
for r = 1:nRewards
    names{r} = sn.reward{r}.name;
    rewardFn = sn.reward{r}.fn;

    for s = 1:nstates
        % Create RewardState for this state row
        stateRow = StateSpaceAggr(s, :);
        rewardState = RewardState(stateRow, sn, nodeToStationMap, classToIndexMap);

        % Try calling with RewardState (new API, single argument)
        try
            R(r, s) = rewardFn(rewardState);
        catch ME
            % If that fails, try backward compatibility with @(state, sn) signature
            try
                R(r, s) = rewardFn(StateSpaceAggr(s, :), sn);
            catch
                rethrow(ME);
            end
        end
    end
end

% Compute transient expected rewards: E[r(X(t))] = pit * R'
Rt = cell(nRewards, 1);
for r = 1:nRewards
    reward_t = pit * R(r, :)';

    metricVal = struct();
    metricVal.t = t;
    metricVal.metric = reward_t;
    metricVal.name = names{r};
    Rt{r} = metricVal;
end

% Filter to specific reward if requested
if ~isempty(rewardName)
    rewardName = char(rewardName);
    idx = find(strcmp(names, rewardName));
    if isempty(idx)
        line_error(mfilename, 'Reward "%s" not found. Available rewards: %s', ...
            rewardName, strjoin(names, ', '));
    end
    Rt = Rt{idx};
    names = names{idx};
end

end
