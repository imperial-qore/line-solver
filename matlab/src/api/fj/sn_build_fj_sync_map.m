function fjSyncMap = sn_build_fj_sync_map(sn)
% FJSYNCMAP = SN_BUILD_FJ_SYNC_MAP(SN)
%
% Builds a fork-join synchronization map from LINE's sn structure.
% For each (Fork, Join) pair, identifies which source nodes feed into the
% join and must be synchronized using mmap_max.
%
% Output:
%   fjSyncMap.nodeSync(joinIdx, srcIdx) = groupId
%     - groupId > 0 means srcIdx belongs to sync group groupId at joinIdx
%     - groupId == 0 means srcIdx is an independent (non-synced) flow
%   fjSyncMap.forkOfGroup(groupId) = forkIdx
%   fjSyncMap.joinOfGroup(groupId) = joinIdx
%   fjSyncMap.nGroups = total number of sync groups
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

I = sn.nnodes;
K = sn.nclasses;

% nodeSync is indexed by (destination node, source node)
nodeSync = zeros(I, I);

groupId = 0;
forkOfGroup = [];
joinOfGroup = [];

% Iterate over all (Fork, Join) pairs from sn.fj
forkIndices = find(any(sn.fj, 2))';
for forkIdx = forkIndices
    joinIdx = find(sn.fj(forkIdx, :));
    if isempty(joinIdx)
        continue;
    end
    for ji = 1:length(joinIdx)
        jnd = joinIdx(ji);
        groupId = groupId + 1;
        forkOfGroup(groupId) = forkIdx;
        joinOfGroup(groupId) = jnd;

        % Find nodes on the parallel path between this fork and join.
        % A node is on the parallel path if:
        %   1. The fork routes to it (for at least one class), AND
        %   2. It routes to the join (for at least one class)
        for ind = 1:I
            if ind == forkIdx || ind == jnd
                continue;
            end
            forkRoutesToNode = false;
            nodeRoutesToJoin = false;
            for k = 1:K
                if sn.rtnodes((forkIdx-1)*K+k, (ind-1)*K+k) > 0
                    forkRoutesToNode = true;
                end
                if sn.rtnodes((ind-1)*K+k, (jnd-1)*K+k) > 0
                    nodeRoutesToJoin = true;
                end
            end
            if forkRoutesToNode && nodeRoutesToJoin
                nodeSync(jnd, ind) = groupId;
            end
        end
    end
end

fjSyncMap.nodeSync = nodeSync;
fjSyncMap.forkOfGroup = forkOfGroup;
fjSyncMap.joinOfGroup = joinOfGroup;
fjSyncMap.nGroups = groupId;
end
