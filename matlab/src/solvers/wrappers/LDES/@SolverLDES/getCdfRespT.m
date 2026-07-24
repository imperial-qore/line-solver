function RD = getCdfRespT(self, R)
% RD = GETCDFRESPT(R)
%
% Empirical response time CDF, measured on the simulated sample path.
%
% A simulator must report what it observed: the base NetworkSolver
% implementation fabricates an exponential law with the right mean, which is
% an analytical fallback and says nothing about the tail. The LDES engine
% records every per-job response time (Solver_ssj.responseTimeSamples ->
% LDESResult.respTimeSamples -> the JSON "transient.respTimeSamples" block,
% exported under --trajectory), so the CDF here is the ecdf of those samples.
%
% RD is an (nstations x nclasses) cell of [F(t), t] matrices, the convention
% every other getCdfRespT follows. A (station,class) pair with no observation
% is left empty rather than filled with a guess.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
RD = cell(sn.nstations, sn.nclasses);
if GlobalConstants.DummyMode
    return
end
if nargin < 2
    R = self.getAvgRespTHandles; %#ok<NASGU>
end

T0 = tic;
% Ask the engine for the per-job response time samples. They are a
% distributional output, not a trajectory, so they have their own CLI flag and
% their own top-level JSON block; the flag also forces the jar runner, since
% the AOT native binary predates it (see solveCli).
data = self.solveCli(self.getOptions, {'--respt-samples'});
samplesCell = ldesResptSamples(data, sn.nstations, sn.nclasses);
if isempty(samplesCell)
    line_error(mfilename, ['The LDES run returned no response time samples, so an empirical CDF ' ...
        'cannot be built. Increase options.samples, or use getPerctRespT(...,''forktail'') ' ...
        'for the analytical fork-join tail.']);
end

for i = 1:sn.nstations
    for r = 1:sn.nclasses
        if i <= size(samplesCell,1) && r <= size(samplesCell,2)
            samples = samplesCell{i, r};
            if ~isempty(samples)
                x = sort(samples(:));
                F = (1:numel(x))' / numel(x);
                % collapse repeated observations, keeping the largest CDF value
                [x, lastIdx] = unique(x, 'last');
                RD{i, r} = [F(lastIdx), x];
            end
        end
    end
end

runtime = toc(T0);
self.setDistribResults(RD, runtime);
end

function out = ldesResptSamples(data, nstations, nclasses)
% OUT = LDESRESPTSAMPLES(DATA, NSTATIONS, NCLASSES)
% Normalize the engine's respTimeSamples block into an (nstations x nclasses)
% cell of column vectors. jsondecode yields a cell of cells when the per-class
% sample counts differ and a numeric array when they happen to match, so both
% shapes have to be accepted.
out = {};
if ~isstruct(data) || ~isfield(data, 'respTimeSamples') || isempty(data.respTimeSamples)
    return
end
rts = data.respTimeSamples;
if ~iscell(rts)
    rts = num2cell(rts, [2 3]);
end
out = cell(nstations, nclasses);
for i = 1:min(numel(rts), nstations)
    stationEntry = rts{i};
    if ~iscell(stationEntry)
        stationEntry = num2cell(stationEntry, 2);
    end
    for r = 1:min(numel(stationEntry), nclasses)
        v = stationEntry{r};
        if ~isempty(v)
            out{i, r} = double(v(:));
        end
    end
end
end
