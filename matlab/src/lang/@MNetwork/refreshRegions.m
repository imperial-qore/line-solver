function sn = refreshRegions(self)
% SN = REFRESHREGIONS() Populate finite capacity region information in sn struct
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Extract finite capacity region information
% region is a cell array of size F (number of regions)
% region{f} is Matrix(M, K+1) where:
%   entry (i,r) = max jobs of class r at station i in region f
%   entry (i,K+1) = global max jobs at station i in region f
%   -1 = infinite capacity
sn = self.sn;
if ~isempty(self.regions)
    F = length(self.regions);
    sn.nregions = F;
    sn.region = cell(F, 1);

    M = sn.nstations;
    K = sn.nclasses;

    % regionrule(f, r) = DropStrategy for class r in region f
    sn.regionrule = DropStrategy.DROP * ones(F, K);  % Default to drop
    % regionweight(f, r) = class weight for class r in region f
    sn.regionweight = ones(F, K);  % Default weight = 1.0
    % regionsz(f, r) = class size/memory for class r in region f
    sn.regionsz = ones(F, K);  % Default size = 1

    for f = 1:F
        fcr = self.regions{f};
        % Matrix with M rows (stations) and K+1 columns (K classes + 1 global)
        regionMatrix = -1 * ones(M, K + 1);  % Initialize all to infinite (-1)

        % Find which stations are in this region and set their capacities
        for n = 1:length(fcr.nodes)
            node = fcr.nodes{n};
            for i = 1:M
                if self.stations{i} == node
                    % Set per-class max jobs for this station in this region
                    % fcr.classMaxJobs is an array indexed by class index
                    for r = 1:K
                        regionMatrix(i, r) = fcr.classMaxJobs(r);
                    end
                    % Set global max jobs for this station in this region (column K+1)
                    regionMatrix(i, K + 1) = fcr.globalMaxJobs;
                    break;
                end
            end
        end

        % Extract drop rule for each class in this region
        for r = 1:K
            % fcr.dropRule is a DropStrategy array indexed by class index
            sn.regionrule(f, r) = fcr.dropRule(r);
        end

        % Extract class weights and sizes for this region
        for r = 1:K
            sn.regionweight(f, r) = fcr.classWeight(r);
            sn.regionsz(f, r) = fcr.classSize(r);
        end

        sn.region{f} = regionMatrix;
    end
else
    sn.nregions = 0;
    sn.region = {};
    sn.regionrule = [];
    sn.regionweight = [];
    sn.regionsz = [];
end
self.sn = sn;
end
