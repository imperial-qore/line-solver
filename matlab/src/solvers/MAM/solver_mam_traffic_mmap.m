function ARV = solver_mam_traffic_mmap(sn, DEP, config, fjSyncMap)
% ARV = SOLVER_MAM_TRAFFIC_MMAP(SN, DEP, CONFIG, FJSYNCMAP)
% FJ-aware traffic solver extending solver_mam_traffic with mmap_max
% synchronization at join points.
%
% DEP{i,r} is the departure process of class r from i in (D0,D1) format
% fjSyncMap is built by sn_build_fj_sync_map
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

I = sn.nnodes;
R = sn.nclasses;

if ~isfield(config, 'fj_sync_q_len')
    config.fj_sync_q_len = 2;
end

% In this function we use indexing over all non-ClassSwitch nodes
non_cs_classes = [];
isNCS = zeros(1,I);
nodeToNCS = zeros(1,I);
NCStoNode = [];
for ind=1:I
    if sn.nodetype(ind) ~= NodeType.ClassSwitch
        non_cs_classes(end+1:end+R)= ((ind-1)*R+1):(ind*R);
        isNCS(ind) = true;
        nodeToNCS(ind) = sum(isNCS);
        NCStoNode(end+1) = ind;
    else
        isNCS(ind) = false;
    end
end

% Hide the nodes that are class switches
rtncs = dtmc_stochcomp(sn.rtnodes, non_cs_classes);
Inc = I - sum(sn.nodetype == NodeType.ClassSwitch);

% DEP is indexed by node (ind, r) - convert to MMAP format
MMAP = cell(I, R);
for ind=1:I
    for r=1:R
        MMAP{ind,r} = DEP{ind,r};
        if isempty(MMAP{ind,r}) || any(any(isnan(MMAP{ind,r}{1})))
            MMAP{ind,r} = {[0],[0],[0]}; % no arrivals from this class
        else
            MMAP{ind,r}{3} = MMAP{ind,r}{2};
        end
    end
end

ARV = cell(Inc,1);
DEP_NCS = cell(Inc,R);
LINKS = cell(Inc,Inc);

% Build the nodeSync matrix in NCS indexing
nodeSyncNCS = zeros(Inc, Inc);
for ind=1:I
    if isNCS(ind)
        inc = nodeToNCS(ind);
        for jnd=1:I
            if isNCS(jnd)
                jnc = nodeToNCS(jnd);
                if fjSyncMap.nodeSync(ind, jnd) > 0
                    nodeSyncNCS(inc, jnc) = fjSyncMap.nodeSync(ind, jnd);
                end
            end
        end
    end
end

% First determine all outgoing flows from all nodes
for ind=1:I
    if isNCS(ind)
        inc = nodeToNCS(ind);
        switch sn.nodetype(ind)
            case {NodeType.Source, NodeType.Delay, NodeType.Queue, NodeType.Fork, NodeType.Join}
                % obtain departure maps (MMAP is now node-indexed)
                if R>1
                    % Order-preserving bounded class-by-class superposition;
                    % see _kb/06-solver-catalog.md for rationale
                    DEP_NCS{inc} = MMAP{ind,1};
                    for rr=2:R
                        DEP_NCS{inc} = mmap_super(DEP_NCS{inc}, MMAP{ind,rr});
                        if length(DEP_NCS{inc}{1}) > config.space_max
                            DEP_NCS{inc} = mmap_compress(DEP_NCS{inc}, config);
                        end
                    end
                else
                    DEP_NCS{inc} = MMAP{ind,1};
                end
                Psplit = zeros(R,Inc*R);
                for r=1:R
                    for jnd = 1:I
                        if isNCS(jnd)
                            jnc = nodeToNCS(jnd);
                            for s=1:R
                                Psplit(r,(jnc-1)*R+s) = rtncs((inc-1)*R+r, (jnc-1)*R+s);
                            end
                        end
                    end
                end

                [Fsplit{1:Inc}] = npfqn_traffic_split_cs(DEP_NCS{inc}, Psplit, config);
                for jnc=1:Inc
                    LINKS{inc,jnc} = Fsplit{jnc};
                    LINKS{inc,jnc} = mmap_normalize(LINKS{inc,jnc});
                end
        end
    end
end

% Then determine all incoming flows, with FJ synchronization
for ind=1:I
    if isNCS(ind) && sn.nodetype(ind) ~= NodeType.Source
        inc = nodeToNCS(ind);

        % Partition incoming links into sync groups and independent flows
        syncGroupsAtNode = unique(nodeSyncNCS(inc, :));
        syncGroupsAtNode = syncGroupsAtNode(syncGroupsAtNode > 0);

        independentFlows = {};
        syncFlows = struct();

        for jnc=1:Inc
            if isempty(LINKS{jnc,inc}) || sum(mmap_lambda(LINKS{jnc,inc})) <= GlobalConstants.FineTol
                continue;
            end
            gid = nodeSyncNCS(inc, jnc);
            if gid == 0
                % Independent flow
                independentFlows{end+1} = LINKS{jnc,inc};
            else
                % Synchronized flow — group by sync group ID
                fname = ['g' num2str(gid)];
                if ~isfield(syncFlows, fname)
                    syncFlows.(fname) = {};
                end
                syncFlows.(fname){end+1} = LINKS{jnc,inc};
            end
        end

        % Process synchronized flows: apply mmap_max iteratively within each group
        syncResults = {};
        for g = 1:length(syncGroupsAtNode)
            gid = syncGroupsAtNode(g);
            fname = ['g' num2str(gid)];
            if ~isfield(syncFlows, fname)
                continue;
            end
            groupFlows = syncFlows.(fname);
            if isempty(groupFlows)
                continue;
            end

            % Start with first flow, iteratively apply mmap_max with remaining
            syncedFlow = groupFlows{1};
            for f = 2:length(groupFlows)
                syncedFlow = mmap_max(syncedFlow, groupFlows{f}, config.fj_sync_q_len);
                % Normalize mmap_max output before use; see _kb/06-solver-catalog.md for rationale
                syncedFlow = mmap_normalize(syncedFlow);

                % Compress after each mmap_max step if state space is too large
                if length(syncedFlow{1}) > config.space_max
                    syncedFlow = mmap_compress(syncedFlow, struct('method', config.compress));
                end
            end
            syncResults{end+1} = syncedFlow;
        end

        % Merge synced flows with independent flows
        allFlows = [syncResults, independentFlows];

        if length(allFlows) > 1
            ARV{ind} = npfqn_traffic_merge(allFlows, config);
        elseif length(allFlows) == 1
            ARV{ind} = allFlows{1};
        else
            % No flows — take a zero-rate link
            for jnc=1:Inc
                if ~isempty(LINKS{jnc,inc})
                    ARV{ind} = LINKS{jnc,inc};
                    break;
                end
            end
            if isempty(ARV{ind})
                ARV{ind} = {[0],[0],[0]};
            end
        end
    else
        ARV{ind} = [];
    end
end
end
