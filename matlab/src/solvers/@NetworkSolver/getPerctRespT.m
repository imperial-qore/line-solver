function [PercRT, PercTable] = getPerctRespT(self, percentiles, jobclass, method)
% [PERCRT, PERCTABLE] = GETPERCTRESPT(SELF, PERCENTILES, JOBCLASS, METHOD)
% Extract response time percentiles from CDF or solver-specific results
%
% Parameters:
%   self - NetworkSolver instance
%   percentiles - Array of percentile values (e.g., [0.90, 0.95, 0.99] or [90, 95, 99])
%   jobclass - (optional) specific job class to retrieve percentiles for
%   method - (optional) 'default' extracts the percentiles from getCdfRespT;
%            'forktail' predicts the fork-join request tail latency with the
%            ForkTail approximation (fj_tail_forktail), which needs only the
%            mean and variance of the per-branch task response times and so
%            applies to heterogeneous branches and mixed service laws
%
% Returns:
%   PercRT - Struct array with fields for each class:
%            .class - class name
%            .percentiles - percentile levels (as percentages)
%            .values - percentile values (response times)
%            .method - extraction method ('cdf' or solver-specific)
%   PercTable - Formatted table for display
%
% This method provides a generic interface for extracting response time
% percentiles. Solvers that compute CDFs via getCdfRespT can use this
% method to extract percentile values. Solvers may override this method
% to provide more efficient or accurate percentile computation.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4 || isempty(method)
    method = 'default';
end

% Normalize percentiles to [0, 1] range
if any(percentiles > 1)
    percentiles = percentiles / 100;
end

sn = self.getStruct();

self.lastPerctMethod = lower(char(method));
if strcmpi(method, 'forktail')
    [PercRT, PercTable] = forktail_percentiles(self, sn, percentiles, ...
        resolve_classes(sn, jobclass, nargin >= 3 && ~isempty(jobclass)));
    return
end

% Determine which classes to return
if nargin < 3 || isempty(jobclass)
    % Return all classes
    classes = 1:sn.nclasses;
else
    % Find class index
    if ischar(jobclass) || isstring(jobclass)
        % Class name provided
        classIdx = find(strcmp(sn.classnames, jobclass));
        if isempty(classIdx)
            line_error(mfilename, 'Job class "%s" not found in model.', jobclass);
        end
        classes = classIdx;
    else
        % Class index provided
        classes = jobclass;
    end
end

% Try to get CDF for percentile extraction
try
    RD = self.getCdfRespT();
catch ME
    line_error(mfilename, ['Unable to compute percentiles. getCdfRespT not available for this solver.\n', ...
        'Error: %s'], ME.message);
end

% Build PercRT per (station, class) from the [F(t), t] CDF cell (CDF value
% first, time second). see _kb/06-solver-catalog.md for rationale
PercRT = struct([]);
idx = 0;

for ci = 1:length(classes)
    r = classes(ci);
    for ist = 1:sn.nstations
        cdfData = [];
        if iscell(RD)
            if size(RD,2) >= r && size(RD,1) >= ist
                cdfData = RD{ist, r};
            elseif size(RD,1) >= r && size(RD,2) == 1
                cdfData = RD{r};       % per-class cell, as some overrides return
            end
        elseif isstruct(RD) && isfield(RD, 'C') && r <= length(RD.C)
            cdfData = RD.C{r};
        elseif isstruct(RD) && isfield(RD, 'CDF') && r <= length(RD.CDF)
            cdfData = RD.CDF{r};
        end
        if isempty(cdfData)
            continue    % e.g. the Source station, which has no response time
        end

        if isstruct(cdfData)
            if isfield(cdfData, 't') && isfield(cdfData, 'p')
                times = cdfData.t(:);
                probs = cdfData.p(:);
            elseif isfield(cdfData, 'x') && isfield(cdfData, 'f')
                times = cdfData.x(:);
                probs = cumsum(cdfData.f(:)); % convert a PDF to a CDF
            else
                line_warning(mfilename, 'Unrecognized CDF structure at station %d for class %d.', ist, r);
                continue;
            end
        elseif size(cdfData, 2) >= 2
            probs = cdfData(:, 1);
            times = cdfData(:, 2);
        else
            line_warning(mfilename, 'Unable to parse CDF data at station %d for class %d.', ist, r);
            continue;
        end

        % Ensure the CDF is monotonic and free of repeated levels
        [probs, sortIdx] = sort(probs);
        times = times(sortIdx);
        [probs, uniqueIdx] = unique(probs);
        times = times(uniqueIdx);

        if length(times) > 1 && length(probs) > 1
            percValues = max(interp1(probs, times, percentiles, 'linear', 'extrap'), 0);
            idx = idx + 1;
            PercRT(idx).station = sn.nodenames{sn.stationToNode(ist)};
            PercRT(idx).class = sn.classnames{r};
            PercRT(idx).percentiles = percentiles;  % kept in fractional form [0,1]
            PercRT(idx).values = percValues;
            PercRT(idx).method = 'cdf';
        end
    end
end

if isempty(PercRT)
    line_warning(mfilename, 'No usable response time CDF was returned by %s.', class(self));
end

% Build PercTable for display
if nargout > 1
    Station = {};
    JobClass = {};
    Percentile = [];
    ResponseTime = [];

    for idx = 1:length(PercRT)
        nPercentiles = length(PercRT(idx).percentiles);
        for p = 1:nPercentiles
            Station{end+1,1} = PercRT(idx).station;
            JobClass{end+1,1} = PercRT(idx).class;
            Percentile(end+1,1) = PercRT(idx).percentiles(p);
            ResponseTime(end+1,1) = PercRT(idx).values(p);
        end
    end

    Station = label(Station);
    JobClass = label(JobClass);
    PercTable = Table(Station, JobClass, Percentile, ResponseTime);
end

end

function classes = resolve_classes(sn, jobclass, given)
% CLASSES = RESOLVE_CLASSES(SN, JOBCLASS, GIVEN)
% Class indexes requested by the caller, by name or by index.
if ~given
    classes = 1:sn.nclasses;
    return
end
if ischar(jobclass) || isstring(jobclass)
    classes = find(strcmp(sn.classnames, jobclass));
    if isempty(classes)
        line_error(mfilename, 'Job class "%s" not found in model.', jobclass);
    end
elseif isobject(jobclass) && isprop(jobclass, 'index')
    classes = jobclass.index;
else
    classes = jobclass;
end
end

function [PercRT, PercTable] = forktail_percentiles(self, sn, percentiles, classes)
% [PERCRT, PERCTABLE] = FORKTAIL_PERCENTILES(SELF, SN, PERCENTILES, CLASSES)
% Fork-join request tail latency by the ForkTail approximation. Each branch
% between a fork and its join must be a single station: the approximation is
% defined on the per-branch task response time, and a multi-station branch has
% no M/G/1 moments of its own (see fj_tail_forktail).
forks = find(sn.nodetype == NodeType.Fork)';
if isempty(forks)
    line_error(mfilename, 'The forktail method requires a model with a Fork node.');
end
if numel(forks) > 1
    line_error(mfilename, 'The forktail method supports a single fork-join pair; this model has %d forks.', numel(forks));
end
f = forks(1);
joinIdx = find(sn.fj(f,:));
if isempty(joinIdx)
    line_error(mfilename, 'The fork node has no matching join; the request response time is undefined.');
end

branches = find(sn.connmatrix(f,:));
for b = branches
    if ~sn.isstation(b)
        line_error(mfilename, 'Branch node %s is not a station; the forktail method needs one queueing station per branch.', sn.nodenames{b});
    end
    if ~sn.connmatrix(b, joinIdx)
        line_error(mfilename, 'Branch station %s does not feed the join directly; the forktail method needs one station per branch.', sn.nodenames{b});
    end
end

[~, UN, ~, TN] = self.getAvg();
PercRT = struct([]);
idx = 0;
for r = classes(:)'
    ET = zeros(1, numel(branches));
    VT = ET;
    rho = ET;
    ok = true;
    for bi = 1:numel(branches)
        ist = sn.nodeToStation(branches(bi));
        lambda = TN(ist, r);
        if lambda <= GlobalConstants.FineTol
            ok = false;    % the class does not traverse this fork
            break
        end
        svc = self.model.nodes{branches(bi)}.getService(self.model.classes{r});
        ES = svc.getMean();
        VS = svc.getVar();
        ES2 = VS + ES^2;
        ES3 = svc.getSkewness()*VS^1.5 + 3*ES*ES2 - 2*ES^3;
        [ET(bi), VT(bi)] = fj_mg1_respt_moments(lambda, ES, ES2, ES3);
        rho(bi) = UN(ist, r);
    end
    if ~ok
        continue
    end
    % ForkTail is a heavy-traffic result; under-shoots at low load.
    % see _kb/06-solver-catalog.md for rationale
    if max(rho) < 0.5
        line_warning(mfilename, 'ForkTail is a heavy-traffic approximation; the busiest branch is at utilization %.2f, so the tail is likely under-predicted.\n', max(rho));
    end
    values = zeros(1, numel(percentiles));
    for pi = 1:numel(percentiles)
        values(pi) = fj_tail_forktail(ET, VT, [], percentiles(pi));
    end
    idx = idx + 1;
    PercRT(idx).class = sn.classnames{r};
    PercRT(idx).percentiles = percentiles;
    PercRT(idx).values = values;
    PercRT(idx).method = 'forktail';
end

if nargout > 1
    JobClass = {}; Percentile = []; ResponseTime = [];
    for i = 1:length(PercRT)
        for p = 1:length(PercRT(i).percentiles)
            JobClass{end+1,1} = PercRT(i).class; %#ok<AGROW>
            Percentile(end+1,1) = PercRT(i).percentiles(p); %#ok<AGROW>
            ResponseTime(end+1,1) = PercRT(i).values(p); %#ok<AGROW>
        end
    end
    JobClass = label(JobClass);
    PercTable = Table(JobClass, Percentile, ResponseTime);
end
end
