function [exportClasses, cacheClasses] = getExportableClasses(self)
% [EXPORTCLASSES, CACHECLASSES] = GETEXPORTABLECLASSES()
%
% Returns a logical array indicating which classes should be exported to JMT.
% Classes with 0 customers are normally skipped unless they are used as
% cache hit/miss classes (since jobs will switch into them).
%
% Returns:
%   exportClasses - logical array, true for classes to export
%   cacheClasses  - array of class indices used as cache hit/miss classes

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
numOfClasses = sn.nclasses;

% Build set of classes used as cache hit/miss classes
cacheClasses = [];
for i = 1:length(sn.nodeparam)
    if ~isempty(sn.nodeparam{i}) && isfield(sn.nodeparam{i}, 'hitclass')
        cacheClasses = [cacheClasses, sn.nodeparam{i}.hitclass(sn.nodeparam{i}.hitclass > 0)];
    end
    if ~isempty(sn.nodeparam{i}) && isfield(sn.nodeparam{i}, 'missclass')
        cacheClasses = [cacheClasses, sn.nodeparam{i}.missclass(sn.nodeparam{i}.missclass > 0)];
    end
end
cacheClasses = unique(cacheClasses);

% Build set of classes that can receive jobs via class-switching
% A class can receive jobs if it's in the same chain as another class with jobs
classSwitchClasses = [];
for c = 1:sn.nchains
    classesInChain = find(sn.chains(c,:) > 0);
    if length(classesInChain) > 1
        % Check if any class in this chain has jobs
        chainHasJobs = false;
        for r = classesInChain
            if sn.njobs(r) > 0
                chainHasJobs = true;
                break;
            end
        end
        if chainHasJobs
            % All classes in this chain can receive jobs via class-switching
            classSwitchClasses = [classSwitchClasses, classesInChain];
        end
    end
end
classSwitchClasses = unique(classSwitchClasses);

% Determine which classes to export
exportClasses = true(1, numOfClasses);
for r = 1:numOfClasses
    % Skip closed classes with 0 customers unless they are cache classes
    % or can receive jobs via class-switching
    if isfinite(sn.njobs(r)) && sn.njobs(r) == 0 && ...
            ~ismember(r, cacheClasses) && ~ismember(r, classSwitchClasses)
        exportClasses(r) = false;
    end
end

end
