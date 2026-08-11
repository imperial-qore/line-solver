function res = runTransientJson(self, numEvents)
% RES = RUNTRANSIENTJSON(NUMEVENTS)
% Run a fully JSON-mediated transient LDES simulation over the horizon
% [0, samples] (or [0, numEvents] when numEvents>0) with trajectory export,
% and return the parsed trajectory. Mirrors the Python-native _run_transient.
%
% RES has fields:
%   .t   - (nTimePoints x 1) time vector, or [] if none
%   .QNt - 1 x Nouter cell; RES.QNt{i} is a 1 x R cell whose {r} entry is the
%          (nTimePoints x 2) [value, time] trajectory of class r at the i-th
%          stateful node, or [] if absent. Nouter matches the engine's outer
%          transient dimension (stateful-node major, as consumed by sample()).
%   .respTimeSamples - (nstations x nclasses) cell of per-job response time
%          samples recorded by the engine (empty when the run produced none).
%          These are the sample-path observations that getCdfRespT turns into
%          an empirical CDF, so that a simulator reports measured percentiles
%          instead of the exponential approximation of the base class.

options = self.getOptions;
% The horizon is passed explicitly via --timespan; options.samples (the
% steady-state event budget) is left untouched since the engine ignores it
% in transient mode.
S = options.samples;
if nargin >= 2 && ~isempty(numEvents) && numEvents > 0
    S = numEvents;
end
extraFlags = {'--timespan', sprintf('0,%.10g', S), '--trajectory'};

% see _kb/06-solver-catalog.md (Wrappers: LDES ensemble transient needs --replications)
% --replications and --numthreads are emitted centrally by solveCli for every
% analysis, steady-state included, so nothing is added here.

data = self.solveCli(options, extraFlags);

res = struct('t', [], 'QNt', {{}}, 'respTimeSamples', {{}});
if ~isstruct(data) || ~isfield(data, 'transient') || isempty(data.transient)
    return;
end
tran = data.transient;
if isfield(tran, 't') && ~isempty(tran.t)
    tvec = ldesJson2mat(tran.t, [], []);
    res.t = tvec(:);
end
if isfield(tran, 'respTimeSamples') && ~isempty(tran.respTimeSamples)
    sn = self.model.getStruct;
    res.respTimeSamples = cell(sn.nstations, sn.nclasses);
    rts = tran.respTimeSamples;
    if ~iscell(rts)
        rts = num2cell(rts, [2 3]);
    end
    for i = 1:min(numel(rts), sn.nstations)
        stationEntry = rts{i};
        if ~iscell(stationEntry)
            stationEntry = num2cell(stationEntry, 2);
        end
        for r = 1:min(numel(stationEntry), sn.nclasses)
            v = stationEntry{r};
            if ~isempty(v)
                res.respTimeSamples{i, r} = double(v(:));
            end
        end
    end
end
if isfield(tran, 'QNt') && ~isempty(tran.QNt)
    sn = self.model.getStruct;
    Nouter = sn.nstateful;
    R = sn.nclasses;
    % Normalize into an Nouter x R cell of [nRows x 2] matrices, then repackage
    % as a 1 x Nouter cell of 1 x R class-cells (the layout sample() consumes).
    grid = ldesTrajCell(tran.QNt, Nouter, R);
    res.QNt = cell(1, Nouter);
    for i = 1:Nouter
        classCells = cell(1, R);
        for r = 1:R
            classCells{r} = grid{i, r};
        end
        res.QNt{i} = classCells;
    end
end
end
