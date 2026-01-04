function [PercRT, PercTable] = getPerctRespT(self, percentiles, jobclass)
% [PERCRT, PERCTABLE] = GETPERCTRESPT(SELF, PERCENTILES, JOBCLASS)
% Extract response time percentiles from CDF or solver-specific results
%
% Parameters:
%   self - NetworkSolver instance
%   percentiles - Array of percentile values (e.g., [0.90, 0.95, 0.99] or [90, 95, 99])
%   jobclass - (optional) specific job class to retrieve percentiles for
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

% Normalize percentiles to [0, 1] range
if any(percentiles > 1)
    percentiles = percentiles / 100;
end

sn = self.getStruct();

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

% Build PercRT struct array by extracting from CDF
PercRT = struct([]);

for idx = 1:length(classes)
    r = classes(idx);

    % CDF structure may vary by solver, try common patterns
    cdfData = [];

    % Pattern 1: RD is a struct with .C field (common in CTMC, Fluid, MAM)
    if isfield(RD, 'C') && r <= length(RD.C)
        cdfData = RD.C{r};
    % Pattern 2: RD is a struct with .CDF field
    elseif isfield(RD, 'CDF') && r <= length(RD.CDF)
        cdfData = RD.CDF{r};
    % Pattern 3: RD is a cell array
    elseif iscell(RD) && r <= length(RD)
        cdfData = RD{r};
    end

    if isempty(cdfData)
        line_warning(mfilename, 'No CDF data available for class %d.', r);
        continue;
    end

    % Extract time and probability arrays from CDF structure
    if isstruct(cdfData)
        if isfield(cdfData, 't') && isfield(cdfData, 'p')
            times = cdfData.t(:);
            probs = cdfData.p(:);
        elseif isfield(cdfData, 'x') && isfield(cdfData, 'f')
            times = cdfData.x(:);
            probs = cumsum(cdfData.f(:)); % Convert PDF to CDF
        else
            line_warning(mfilename, 'Unrecognized CDF structure for class %d.', r);
            continue;
        end
    elseif size(cdfData, 2) >= 2
        % Matrix format: [time, prob]
        times = cdfData(:, 1);
        probs = cdfData(:, 2);
    else
        line_warning(mfilename, 'Unable to parse CDF data for class %d.', r);
        continue;
    end

    % Ensure CDF is monotonic and normalized
    [probs, sortIdx] = sort(probs);
    times = times(sortIdx);

    % Remove duplicates
    [probs, uniqueIdx] = unique(probs);
    times = times(uniqueIdx);

    % Interpolate to find response times at requested percentiles
    if length(times) > 1 && length(probs) > 1
        percValues = interp1(probs, times, percentiles, 'linear', 'extrap');

        % Ensure non-negative response times
        percValues = max(percValues, 0);

        PercRT(idx).class = sn.classnames{r};
        PercRT(idx).percentiles = percentiles;  % Keep in fractional form [0,1]
        PercRT(idx).values = percValues;
        PercRT(idx).method = 'cdf';
    else
        line_warning(mfilename, 'Insufficient CDF data points for class %d.', r);
    end
end

% Build PercTable for display
if nargout > 1
    JobClass = {};
    Percentile = [];
    ResponseTime = [];

    for idx = 1:length(PercRT)
        nPercentiles = length(PercRT(idx).percentiles);
        for p = 1:nPercentiles
            JobClass{end+1,1} = PercRT(idx).class;
            Percentile(end+1,1) = PercRT(idx).percentiles(p);
            ResponseTime(end+1,1) = PercRT(idx).values(p);
        end
    end

    JobClass = label(JobClass);
    PercTable = Table(JobClass, Percentile, ResponseTime);
end

end
