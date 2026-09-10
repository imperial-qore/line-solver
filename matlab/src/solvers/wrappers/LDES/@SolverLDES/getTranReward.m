function [Rt, t, names] = getTranReward(self, rewardName)
% [RT, T, NAMES] = GETTRANREWARD([REWARDNAME]) Transient reward r(X(t)) via LDES
%
% Returns the transient reward trajectory r(X(t)) for each reward function defined
% on the model via setReward, sampled along the simulated path. The LDES simulator
% exports the integer joint-state trajectory; the reward functions are evaluated here
% on each state, so the result is correct also for nonlinear rewards.
%
% OUTPUTS:
%   Rt    - Cell array {nRewards x 1} of structs with fields .t, .metric, .name
%           (or a single struct if REWARDNAME is given)
%   t     - Time vector [nPoints x 1]
%   names - Cell array of reward names (or a single name if REWARDNAME is given)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    rewardName = [];
end

if GlobalConstants.DummyMode
    Rt = {}; t = []; names = {};
    return
end

sn = self.model.getStruct(true);
if isempty(sn.reward)
    line_error(mfilename, 'No rewards defined. Use model.setReward(name, @(state) ...) before calling getTranReward.');
end

% Run the LDES simulation (fully JSON-mediated) requesting the time-ordered
% joint-state trajectory via --export-histogram --trajectory.
data = self.solveCli(self.getOptions, {'--export-histogram', '--trajectory'});
space = [];
t = [];
if isstruct(data) && isfield(data, 'stateHistogram') && ~isempty(data.stateHistogram)
    space = ldesJson2mat(ldesGetField(data.stateHistogram, 'trajSpace', []), [], []);
    t = ldesJson2mat(ldesGetField(data.stateHistogram, 'trajTime', []), [], []);
end
t = t(:);

nRewardsAll = length(sn.reward);
allNames = cell(nRewardsAll, 1);
for r = 1:nRewardsAll
    allNames{r} = sn.reward{r}.name;
end

% Build index maps for RewardState (station-major aggregated layout)
nodeToStationMap = configureDictionary('int32', 'int32');
classToIndexMap = configureDictionary('int32', 'int32');
for ind = 1:sn.nnodes
    if sn.isstation(ind)
        nodeToStationMap(int32(ind)) = sn.nodeToStation(ind);
    end
end
for r = 1:sn.nclasses
    classToIndexMap(int32(r)) = r;
end

nPoints = size(space, 1);
Rt = cell(nRewardsAll, 1);
names = allNames;
for r = 1:nRewardsAll
    rewardFn = sn.reward{r}.fn;
    metric = zeros(nPoints, 1);
    for s = 1:nPoints
        stateRow = space(s, :);
        rewardState = RewardState(stateRow, sn, nodeToStationMap, classToIndexMap);
        try
            metric(s) = rewardFn(rewardState);
        catch ME
            try
                metric(s) = rewardFn(stateRow, sn);
            catch
                rethrow(ME);
            end
        end
    end
    Rt{r} = struct('t', t, 'metric', metric, 'name', allNames{r});
end

% Filter to a single named reward if requested
if ~isempty(rewardName)
    idx = find(strcmp(allNames, rewardName), 1);
    if isempty(idx)
        line_error(mfilename, sprintf('Reward ''%s'' not found.', rewardName));
    end
    Rt = Rt{idx};
    names = allNames{idx};
end

end
