function [PercRT, PercTable] = getPerctRespT(self, percentiles, jobclass)
% [PERCRT, PERCTABLE] = GETPERCTRESPT(SELF, PERCENTILES, JOBCLASS)
% Retrieve response time percentile results from FJ_codes analysis or CDF
%
% Parameters:
%   self - SolverMAM instance
%   percentiles - Array of percentile values to compute (e.g., [0.90, 0.95, 0.99] or [90, 95, 99])
%   jobclass - (optional) specific job class to retrieve percentiles for
%
% Returns:
%   PercRT - Struct array with fields for each class:
%            .class - class name
%            .percentiles - percentile levels (e.g., [90, 95, 99])
%            .values - percentile values (response times)
%            .K - number of parallel queues in Fork-Join
%   PercTable - Formatted table for display with columns:
%               JobClass, Percentile, ResponseTime
%
% Reference:
%   Z. Qiu, J.F. Pérez, and P. Harrison, "Beyond the Mean in Fork-Join Queues:
%   Efficient Approximation for Response-Time Tails", IFIP Performance 2015.
%   Copyright 2015 Imperial College London
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Normalize percentiles to [0, 1] range
if any(percentiles > 1)
    percentiles = percentiles / 100;
end

sn = self.getStruct();

% Check if FJ_codes percentile results are available
if isfield(self.result, 'Percentile') && ~isempty(self.result.Percentile)
    % Use FJ_codes results (pre-computed)
    percResults = self.result.Percentile;
    useFJCodes = true;
else
    % Fall back to CDF-based percentile extraction
    useFJCodes = false;
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

% Build PercRT struct array
PercRT = struct([]);

if useFJCodes
    % Use pre-computed FJ_codes results
    for idx = 1:length(classes)
        r = classes(idx);

        if r > length(percResults.RT) || isempty(percResults.RT{r})
            line_warning(mfilename, 'No percentile results for class %d.', r);
            continue;
        end

        % Extract percentile data for this class
        classPercResults = percResults.RT{r};

        % Interpolate to get requested percentiles (if different from stored)
        storedPercentiles = classPercResults.percentiles;  % In percentage form from mainFJ
        storedValues = classPercResults.RTp;

        % Convert requested percentiles to percentage form for interpolation
        % (percentiles are already normalized to [0,1] on line 28)
        requestedPercentiles = percentiles * 100;

        % Interpolate to requested percentiles
        interpValues = interp1(storedPercentiles, storedValues, requestedPercentiles, 'linear', 'extrap');

        PercRT(idx).class = sn.classnames{r};
        PercRT(idx).percentiles = percentiles;  % Keep in fractional form [0,1]
        PercRT(idx).values = interpValues;
        PercRT(idx).K = percResults.K;
        PercRT(idx).method = percResults.method;
    end
else
    % Extract percentiles from CDF
    try
        RD = self.getCdfRespT();
    catch
        line_error(mfilename, 'Unable to compute percentiles. getCdfRespT not available for this model.');
    end

    for idx = 1:length(classes)
        r = classes(idx);

        % Find CDF data for this class
        % RD is a cell array with format RD{station, class} = [F, X]
        % where F are CDF probabilities and X are time values
        % Search through stations to find non-empty CDF data
        cdfData = [];
        for i = 1:size(RD, 1)
            if r <= size(RD, 2) && ~isempty(RD{i, r})
                cdfData = RD{i, r};
                break;
            end
        end

        if ~isempty(cdfData)
            % CDF format is [F, X] where F = probabilities, X = times
            probs = cdfData(:, 1);
            times = cdfData(:, 2);

            % Remove duplicate probability values for interpolation
            [probs_unique, idx_unique] = unique(probs, 'last');
            times_unique = times(idx_unique);

            % Interpolate to find times at requested percentiles
            percValues = interp1(probs_unique, times_unique, percentiles, 'linear', 'extrap');

            PercRT(idx).class = sn.classnames{r};
            PercRT(idx).percentiles = percentiles;  % Keep in fractional form [0,1]
            PercRT(idx).values = percValues;
            PercRT(idx).method = 'cdf';
        end
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
