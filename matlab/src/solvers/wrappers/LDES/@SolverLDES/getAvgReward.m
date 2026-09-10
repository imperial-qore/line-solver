function [R, names] = getAvgReward(self)
% [R, NAMES] = GETAVGREWARD() Steady-state expected Markov reward via LDES simulation
%
% Computes the steady-state expected value E[r]=sum_s pi(s) r(s) of each reward
% function defined on the model via setReward. The LDES simulator exports the exact
% joint-state residence-time histogram; the reward functions are evaluated here on
% each visited state, so the result is correct also for nonlinear rewards (e.g. E[n^2]).
%
% OUTPUTS:
%   R     - Vector of steady-state expected rewards [nRewards x 1]
%   names - Cell array of reward names {nRewards x 1}
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if GlobalConstants.DummyMode
    R = []; names = {};
    return
end

sn = self.model.getStruct(true);
if isempty(sn.reward)
    line_error(mfilename, 'No rewards defined. Use model.setReward(name, @(state) ...) before calling getAvgReward.');
end

% Run the LDES simulation (fully JSON-mediated) requesting the exact
% joint-state residence-time histogram via --export-histogram.
data = self.solveCli(self.getOptions, {'--export-histogram'});
space = [];
time = [];
if isstruct(data) && isfield(data, 'stateHistogram') && ~isempty(data.stateHistogram)
    space = ldesJson2mat(ldesGetField(data.stateHistogram, 'space', []), [], []);
    time = ldesJson2mat(ldesGetField(data.stateHistogram, 'time', []), [], []);
end

nRewards = length(sn.reward);
names = cell(nRewards, 1);
R = zeros(nRewards, 1);

if isempty(space) || isempty(time) || sum(time) <= 0
    for r = 1:nRewards
        names{r} = sn.reward{r}.name;
    end
    return
end

w = time(:) / sum(time);

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

nstates = size(space, 1);
for r = 1:nRewards
    names{r} = sn.reward{r}.name;
    rewardFn = sn.reward{r}.fn;
    acc = 0;
    for s = 1:nstates
        stateRow = space(s, :);
        rewardState = RewardState(stateRow, sn, nodeToStationMap, classToIndexMap);
        try
            val = rewardFn(rewardState);
        catch ME
            try
                val = rewardFn(stateRow, sn);
            catch
                rethrow(ME);
            end
        end
        acc = acc + w(s) * val;
    end
    R(r) = acc;
end

end
