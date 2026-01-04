function [result, parsed] = getResults(self)
% [RESULT, PARSED] = GETRESULTS()

options = self.getOptions;

% Check if result file exists, if not run the analyzer first
% Need to handle case where filePath/fileName haven't been set yet
needsAnalyzer = false;
if isempty(self.filePath) || isempty(self.fileName)
    needsAnalyzer = true;
else
    switch options.method
        case {'jsim','default'}
            fileName = [getFileName(self),'.jsim-result.jsim'];
        otherwise
            fileName = [getFileName(self),'.jmva-result.jmva'];
    end
    filePath = [getFilePath(self),filesep,fileName];
    if ~exist(filePath,'file')
        needsAnalyzer = true;
    end
end
if needsAnalyzer
    % Result file doesn't exist, need to run the simulation first
    % runAnalyzer() calls getResults() internally and stores result in self.result
    % It may also clean up the temp directory, so we return the cached result
    runAnalyzer(self, options);
    result = self.result;
    parsed = struct(); % parsed data not available after cleanup
    return;
end

switch options.method
    case {'jsim','default'}
        [result, parsed] = self.getResultsJSIM;
    otherwise
        [result, parsed] = self.getResultsJMVA;
end

sn = self.model.getStruct;

% Extend matrices to include FCR rows at indices nstations+1 to nstations+nregions
% This allows FCR metrics to be stored directly without custom fields
totalRows = sn.nstations + sn.nregions;
result.Avg.Q = zeros(totalRows, sn.nclasses);
result.Avg.U = zeros(totalRows, sn.nclasses);
result.Avg.R = zeros(totalRows, sn.nclasses);
result.Avg.T = zeros(totalRows, sn.nclasses);
result.Avg.A = zeros(totalRows, sn.nclasses);
result.Avg.W = zeros(totalRows, sn.nclasses);
result.Avg.Tard = zeros(totalRows, sn.nclasses);
result.Avg.SysTard = zeros(1, sn.nclasses);

% Set FCR U and A to NaN (JMT doesn't provide these for regions)
if sn.nregions > 0
    fcrRowStart = sn.nstations + 1;
    fcrRowEnd = sn.nstations + sn.nregions;
    result.Avg.U(fcrRowStart:fcrRowEnd, :) = NaN;
    result.Avg.A(fcrRowStart:fcrRowEnd, :) = NaN;
end

% Initialize CI storage (half-widths)
[confintEnabled, ~] = Solver.parseConfInt(options.confint);
if confintEnabled
    result.Avg.QCI = zeros(sn.nstations, sn.nclasses);
    result.Avg.UCI = zeros(sn.nstations, sn.nclasses);
    result.Avg.RCI = zeros(sn.nstations, sn.nclasses);
    result.Avg.TCI = zeros(sn.nstations, sn.nclasses);
    result.Avg.ACI = zeros(sn.nstations, sn.nclasses);
    result.Avg.WCI = zeros(sn.nstations, sn.nclasses);
    result.Avg.TardCI = zeros(sn.nstations, sn.nclasses);
    result.Avg.SysTardCI = zeros(1, sn.nclasses);
end

for m=1:length(result.metric)
    metric = result.metric{m};
    % Compute CI half-width from JMT bounds
    if confintEnabled && isfield(metric, 'upperLimit') && isfield(metric, 'lowerLimit')
        ciHalfWidth = (metric.upperLimit - metric.lowerLimit) / 2;
    else
        ciHalfWidth = 0;
    end

    % Check if this is a region (FCR) metric
    % FCR metrics can be identified by either nodeType='region' or station name starting with 'FCRegion'
    isFCRMetric = false;
    fcrStationName = '';
    if isfield(metric, 'nodeType') && strcmp(metric.nodeType, 'region')
        isFCRMetric = true;
        if isfield(metric, 'station')
            fcrStationName = metric.station;
        end
    elseif isfield(metric, 'station') && startsWith(metric.station, 'FCRegion')
        % JMT may not echo nodeType, detect by station name pattern
        isFCRMetric = true;
        fcrStationName = metric.station;
    end

    if isFCRMetric && ~isempty(fcrStationName) && startsWith(fcrStationName, 'FCRegion')
        fcrIndexStr = fcrStationName(9:end);  % Extract number after "FCRegion"
        fcrIndex = str2double(fcrIndexStr);
        if ~isnan(fcrIndex) && fcrIndex >= 1 && fcrIndex <= sn.nregions
            % FCR row index in result matrices: nstations + fcrIndex
            fcrRowIdx = sn.nstations + fcrIndex;
            % FCR metrics from JMT are aggregate (not per-class), distribute across classes
            switch metric.measureType
                case 'Number of Customers'
                    for r = 1:sn.nclasses
                        result.Avg.Q(fcrRowIdx, r) = metric.meanValue / sn.nclasses;
                    end
                case 'Response Time'
                    for r = 1:sn.nclasses
                        result.Avg.R(fcrRowIdx, r) = metric.meanValue;
                    end
                case 'Residence Time'
                    for r = 1:sn.nclasses
                        result.Avg.W(fcrRowIdx, r) = metric.meanValue;
                    end
                case 'Throughput'
                    for r = 1:sn.nclasses
                        result.Avg.T(fcrRowIdx, r) = metric.meanValue / sn.nclasses;
                    end
                case 'Arrival Rate'
                    for r = 1:sn.nclasses
                        result.Avg.A(fcrRowIdx, r) = metric.meanValue / sn.nclasses;
                    end
            end
        end
        continue;  % Skip to next metric
    end

    switch metric.measureType
        case MetricType.toText(MetricType.QLen)
            i = sn.nodeToStation(find(sn.nodenames == metric.station));
            r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
            if isinf(sn.njobs(r))
                result.Avg.Q(i,r) = metric.meanValue;
                if confintEnabled
                    result.Avg.QCI(i,r) = ciHalfWidth;
                end
            else % 'closed'
                N = sn.njobs;
                chainIdx = find(sn.classnames == metric.class);
                if metric.analyzedSamples > sum(sn.njobs(chainIdx))  % for a class to be considered recurrent we ask more samples than jobs in the corresponding closed chain
                    result.Avg.Q(i,r) = metric.meanValue;
                    if confintEnabled
                        result.Avg.QCI(i,r) = ciHalfWidth;
                    end
                else
                    result.Avg.Q(i,r) = 0;
                end
            end
        case MetricType.toText(MetricType.Util)
            i = sn.nodeToStation(find(sn.nodenames == metric.station));
            r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
            if isinf(sn.njobs(r))
                result.Avg.U(i,r) = metric.meanValue;
                if confintEnabled
                    result.Avg.UCI(i,r) = ciHalfWidth;
                end
            else % 'closed'
                N = sn.njobs;
                chainIdx = find(sn.classnames == metric.class);
                if metric.analyzedSamples > sum(sn.njobs(chainIdx))  % for a class to be considered recurrent we ask more samples than jobs in the corresponding closed chain
                    result.Avg.U(i,r) = metric.meanValue;
                    if confintEnabled
                        result.Avg.UCI(i,r) = ciHalfWidth;
                    end
                else
                    result.Avg.U(i,r) = 0;
                end
            end
        case MetricType.toText(MetricType.RespT)
            i = sn.nodeToStation(find(sn.nodenames == metric.station));
            r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
            if isinf(sn.njobs(r))
                result.Avg.R(i,r) = metric.meanValue;
                if confintEnabled
                    result.Avg.RCI(i,r) = ciHalfWidth;
                end
            else % 'closed'
                N = sn.njobs;
                chainIdx = find(sn.classnames == metric.class);
                if metric.analyzedSamples > sum(sn.njobs(chainIdx))  % for a class to be considered recurrent we ask more samples than jobs in the corresponding closed chain
                    result.Avg.R(i,r) = metric.meanValue;
                    if confintEnabled
                        result.Avg.RCI(i,r) = ciHalfWidth;
                    end
                else
                    result.Avg.R(i,r) = 0;
                end
            end
        case MetricType.toText(MetricType.ResidT)
            % JMT ResidT is inconsistently defined with LINE's on some
            % difficult class switching cases, hence we recompute it at the
            % level of the NetworkSolver class to preserve consistency

%             i = sn.nodeToStation(find(sn.nodenames == metric.station));
%             r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
%             if isinf(sn.njobs(r))
%                 result.Avg.W(i,r) = metric.meanValue;
%             else % 'closed'
%                 N = sn.njobs;
%                 chainIdx = find(sn.classnames == metric.class);
%                 if metric.analyzedSamples > sum(sn.njobs(chainIdx))  % for a class to be considered recurrent we ask more samples than jobs in the corresponding closed chain
%                     result.Avg.W(i,r) = metric.meanValue;
%                 else
%                     result.Avg.W(i,r) = 0;
%                 end
%             end
        case MetricType.toText(MetricType.ArvR)
            i = sn.nodeToStation(find(sn.nodenames == metric.station));
            r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
            if isinf(sn.njobs(r))
                result.Avg.A(i,r) = metric.meanValue;
                if confintEnabled
                    result.Avg.ACI(i,r) = ciHalfWidth;
                end
            else % 'closed'
                N = sn.njobs;
                chainIdx = find(sn.classnames == metric.class);
                if metric.analyzedSamples > sum(sn.njobs(chainIdx))  % for a class to be considered recurrent we ask more samples than jobs in the corresponding closed chain
                    result.Avg.A(i,r) = metric.meanValue;
                    if confintEnabled
                        result.Avg.ACI(i,r) = ciHalfWidth;
                    end
                else
                    result.Avg.A(i,r) = 0;
                end
            end
        case MetricType.toText(MetricType.Tput)
            i = sn.nodeToStation(find(sn.nodenames == metric.station));
            r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
            if isinf(sn.njobs(r))
                result.Avg.T(i,r) = metric.meanValue;
                if confintEnabled
                    result.Avg.TCI(i,r) = ciHalfWidth;
                end
            else % 'closed'
                N = sn.njobs;
                chainIdx = find(sn.classnames == metric.class);
                if metric.analyzedSamples > sum(sn.njobs(chainIdx))  % for a class to be considered recurrent we ask more samples than jobs in the corresponding closed chain
                    result.Avg.T(i,r) = metric.meanValue;
                    if confintEnabled
                        result.Avg.TCI(i,r) = ciHalfWidth;
                    end
                else
                    result.Avg.T(i,r) = 0;
                end
            end
        case MetricType.toText(MetricType.Tard)
            i = sn.nodeToStation(find(sn.nodenames == metric.station));
            r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
            if isinf(sn.njobs(r))
                result.Avg.Tard(i,r) = metric.meanValue;
                if confintEnabled
                    result.Avg.TardCI(i,r) = ciHalfWidth;
                end
            else % 'closed'
                N = sn.njobs;
                chainIdx = find(sn.classnames == metric.class);
                if metric.analyzedSamples > sum(sn.njobs(chainIdx))
                    result.Avg.Tard(i,r) = metric.meanValue;
                    if confintEnabled
                        result.Avg.TardCI(i,r) = ciHalfWidth;
                    end
                else
                    result.Avg.Tard(i,r) = 0;
                end
            end
        case MetricType.toText(MetricType.SysTard)
            r = find(cellfun(@(c) strcmp(c,metric.class), sn.classnames));
            if isinf(sn.njobs(r))
                result.Avg.SysTard(1,r) = metric.meanValue;
                if confintEnabled
                    result.Avg.SysTardCI(1,r) = ciHalfWidth;
                end
            else % 'closed'
                N = sn.njobs;
                chainIdx = find(sn.classnames == metric.class);
                if metric.analyzedSamples > sum(sn.njobs(chainIdx))
                    result.Avg.SysTard(1,r) = metric.meanValue;
                    if confintEnabled
                        result.Avg.SysTardCI(1,r) = ciHalfWidth;
                    end
                else
                    result.Avg.SysTard(1,r) = 0;
                end
            end

        case 'Cache Hit Rate'
            % Store cache hit rate for later processing
            % Find the cache node by station name
            cacheNodeIdx = find(cellfun(@(c) strcmp(c, metric.station), sn.nodenames));
            hitClassIdx = find(cellfun(@(c) strcmp(c, metric.class), sn.classnames));
            if ~isempty(cacheNodeIdx) && ~isempty(hitClassIdx) && sn.nodetype(cacheNodeIdx) == NodeType.Cache
                % Initialize cache hit rate storage if not exists
                if ~isfield(result, 'CacheHitRate')
                    result.CacheHitRate = struct();
                end
                % Store hit rate keyed by cache node index
                fieldName = sprintf('node%d', cacheNodeIdx);
                if ~isfield(result.CacheHitRate, fieldName)
                    result.CacheHitRate.(fieldName) = zeros(1, sn.nclasses);
                end
                % Find the original class (the one that maps to this hit class)
                hitclass = sn.nodeparam{cacheNodeIdx}.hitclass;
                for r = 1:length(hitclass)
                    if hitclass(r) == hitClassIdx
                        result.CacheHitRate.(fieldName)(r) = metric.meanValue;
                        break;
                    end
                end
            end

    end
end
% Set self.result AFTER the for loop completes, so FCR metrics are included
self.result = result;
end
