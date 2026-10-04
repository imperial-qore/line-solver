function [b, count, names] = getAvgBusyPeriod(self, stations, jobclass, n)
% [B, COUNT, NAMES] = GETAVGBUSYPERIOD(STATIONS, JOBCLASS, N)
%
% Mean busy period of order N for a set of stations, measured along the
% simulated sample path. A busy period of order n runs from the instant an
% arrival raises the jobs held by the set to n up to the instant the set falls
% back below n (H. Daduna, "Busy Periods for Subnetworks in Stochastic
% Networks: Mean Value Analysis", J. ACM 35(3), 1988).
%
% Input:
%   stations - stations forming the subnetwork, as objects, names or station
%              indexes; omitted or empty returns the whole measured table
%   jobclass - job class object, name or index; -1 or empty counts every class
%   n        - busy period order(s), 1 for the ordinary busy period
%
% Output:
%   b     - mean busy period duration(s); NaN where no period completed
%   count - number of completed busy periods behind each mean
%   names - target names, returned with the whole table
%
% Periods still in progress when the warmup ends are discarded, since their
% start is not observable. A period of zero length is an artifact of two
% simultaneous events and is not an observation.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if GlobalConstants.DummyMode
    b = NaN; count = 0; names = {};
    return
end

if nargin < 2
    stations = [];
end
if nargin < 3 || isempty(jobclass)
    jobclass = -1;
end
if nargin < 4 || isempty(n)
    n = 1;
end

sn = self.model.getStruct(false);
subnet = local_stations(self, sn, stations);
cls = local_class(self, sn, jobclass);

args = {'--busyperiod', num2str(max(n))};
if numel(subnet) > 1
    % station indexes cross the wire zero-based
    args = [args, {'--busyperiod-subnet', strjoin(arrayfun(@(s) num2str(s-1), ...
        subnet, 'UniformOutput', false), ',')}];
end
data = self.solveCli(self.getOptions, args);

if ~isstruct(data) || ~isfield(data, 'busyPeriods') || isempty(data.busyPeriods)
    line_error(mfilename, 'The LDES engine returned no busy period measurement.');
end
targets = data.busyPeriods.targets;
if ~iscell(targets)
    targets = num2cell(targets);
end

names = cell(numel(targets), 1);
for t = 1:numel(targets)
    names{t} = targets{t}.name;
end

if isempty(stations)
    % whole table: one row per target, one column per order
    b = zeros(numel(targets), max(n));
    count = zeros(numel(targets), max(n));
    for t = 1:numel(targets)
        b(t, :) = local_vec(targets{t}.mean);
        count(t, :) = local_vec(targets{t}.count);
    end
    b(count == 0) = NaN;
    return
end

wanted = sort(subnet(:))';
for t = 1:numel(targets)
    tgt = targets{t};
    st = sort(local_vec(tgt.stations)) + 1;
    if numel(st) == numel(wanted) && all(st == wanted) && tgt.class == cls
        allMean = local_vec(tgt.mean);
        allCount = local_vec(tgt.count);
        b = allMean(n);
        count = allCount(n);
        b(count == 0) = NaN;
        names = tgt.name;
        return
    end
end
line_error(mfilename, 'No busy period target matches the requested stations and class.');
end

function v = local_vec(x)
if iscell(x)
    v = cell2mat(x);
else
    v = x;
end
v = v(:)';
end

function subnet = local_stations(self, sn, stations)
if isempty(stations)
    subnet = [];
    return
end
subnet = zeros(1, numel(stations));
for t = 1:numel(stations)
    if isnumeric(stations)
        subnet(t) = stations(t);
    else
        if iscell(stations)
            st = stations{t};
        else
            st = stations(t);
        end
        if ischar(st) || isstring(st)
            subnet(t) = sn.nodeToStation(self.model.getNodeIndex(char(st)));
        else
            subnet(t) = sn.nodeToStation(st.index);
        end
    end
end
end

function cls = local_class(self, sn, jobclass) %#ok<INUSD>
if isnumeric(jobclass)
    if jobclass < 0
        cls = -1;
    else
        cls = jobclass - 1;  % zero-based on the wire
    end
elseif ischar(jobclass) || isstring(jobclass)
    cls = find(strcmp(sn.classnames, char(jobclass)), 1) - 1;
else
    cls = jobclass.index - 1;
end
end
