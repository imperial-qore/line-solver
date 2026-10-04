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
    % regionlincon{f,1}/regionlincon{f,2} = linear constraint pair (A,b) for region f
    sn.regionlincon = cell(F, 2);
    % regionmaxmem{f} = Matrix(M,1): global memory budget of region f replicated
    % on each member station row, -1 = unbounded. Kept separate from the
    % job-count column so a region can cap memory while leaving jobs unbounded;
    % mirrors the Java NetworkStruct.regionmaxmem field.
    sn.regionmaxmem = cell(F, 1);
% see _kb/04-networkstruct.md (refreshRegions.m) for rationale
    sn.regionmembers = cell(F, 1);

    for f = 1:F
        fcr = self.regions{f};
        % Matrix with M rows (stations) and K+1 columns (K classes + 1 global)
        regionMatrix = -1 * ones(M, K + 1);  % Initialize all to infinite (-1)
        regionMemMatrix = -1 * ones(M, 1);   % Initialize all to unbounded (-1)
        regionMemberMask = false(M, 1);      % membership, independent of the caps

        % Find which stations are in this region and set their capacities
        for n = 1:length(fcr.nodes)
            node = fcr.nodes{n};
            for i = 1:M
                if self.stations{i} == node
                    regionMemberMask(i) = true;
                    % see _kb/04-networkstruct.md (refreshRegions.m) for rationale
                    for r = 1:K
                        if r > length(fcr.classMaxJobs)
                            regionMatrix(i, r) = FiniteCapacityRegion.UNBOUNDED;
                            continue;
                        end
                        cap_r = fcr.classMaxJobs(r);
                        if fcr.classMaxMemory(r) ~= FiniteCapacityRegion.UNBOUNDED && fcr.classSize(r) > 0
                            memjobs = floor(fcr.classMaxMemory(r) / fcr.classSize(r));
                            if cap_r == FiniteCapacityRegion.UNBOUNDED
                                cap_r = memjobs;
                            else
                                cap_r = min(cap_r, memjobs);
                            end
                        end
                        regionMatrix(i, r) = cap_r;
                    end
                    % Set global max jobs for this station in this region (column K+1)
                    regionMatrix(i, K + 1) = fcr.globalMaxJobs;
                    % Replicate the region-global memory budget on this member row
                    regionMemMatrix(i, 1) = fcr.globalMaxMemory;
                    break;
                end
            end
        end

        % Extract drop rule for each class in this region. Classes beyond
        % the region's vectors default to WAITQ, weight 1 and size 1.
        for r = 1:K
            % fcr.dropRule is a DropStrategy array indexed by class index
            if r > length(fcr.dropRule)
                sn.regionrule(f, r) = DropStrategy.WAITQ;
            else
                sn.regionrule(f, r) = fcr.dropRule(r);
            end
        end

        % Extract class weights and sizes for this region
        for r = 1:K
            if r <= length(fcr.classWeight)
                sn.regionweight(f, r) = fcr.classWeight(r);
            end
            if r <= length(fcr.classSize)
                sn.regionsz(f, r) = fcr.classSize(r);
            end
        end

        sn.region{f} = regionMatrix;
        sn.regionmaxmem{f} = regionMemMatrix;
        sn.regionmembers{f} = regionMemberMask;

        % Capture the linear-constraint pair (A,b) on the same cell row if
        % set. A is padded with zero columns up to K, so classes appended
        % after the region was created stay unconstrained and engines can
        % index A by class without going out of bounds.
        if ismethod(fcr, 'hasLinearConstraints') && fcr.hasLinearConstraints()
            [linConA, linConB] = fcr.getLinearConstraints();
            if size(linConA, 2) < K
                linConA(:, end+1:K) = 0;
            end
            sn.regionlincon{f, 1} = linConA;
            sn.regionlincon{f, 2} = linConB;
        end
    end
else
    sn.nregions = 0;
    sn.region = {};
    sn.regionrule = [];
    sn.regionweight = [];
    sn.regionsz = [];
    sn.regionlincon = {};
    sn.regionmaxmem = {};
    sn.regionmembers = {};
end
self.sn = sn;
end
