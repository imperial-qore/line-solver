function initInterlock(self)
% INITINTERLOCK Build the interlock path table and locate common parents
%
% Interlocking arises when requests issued by one client reach a common
% lower-level server along two or more independent paths, so that arrivals
% a layer decomposition treats as independent are in fact correlated. The
% static half of the correction is built here, once per solve:
%
%   Phase A: the path table path(a,b) of Franks (1999), Sec. 4.2, holding
%            the calls to entry b caused by one invocation of entry a, with
%            a unit diagonal. LINE keeps a second table restricted to the
%            phase-1 flow, so that deferred (phase-2) work can be scored
%            separately when the correction is applied
%   Phase B: the common-parent finder of Fig. 4.2, retaining only the
%            entries at which the flow genuinely splits
%   Phase C: the source tasks and the source count n_s of Eq. (4.7)
%
% Only the flow computation in UPDATEPOPULATIONS depends on the iterate,
% so the tables built here are reused across the fixed-point iteration.
%
% The phase-aware tables, the branch-point test and the source count are
% refinements beyond the published algorithm, which assumes one path table
% and counts source tasks directly.
%
% Reference: G. Franks, "Performance Analysis of Distributed Server
% Systems", PhD thesis, Carleton University, 1999, Ch. 4; published as
% G. Franks, "Traffic dependencies in client-server systems and their
% effect on performance prediction", IEEE IPDS, 1995, pp. 24-33.

lqn = self.lqn;

% Phase A: path table of Sec. 4.2, il_all(a,b) = calls to entry b per
% invocation of entry a; il_ph1 restricts the same count to phase-1 flow
nentries = lqn.nentries;
il_all = zeros(nentries, nentries);
il_ph1 = zeros(nentries, nentries);

for e = 1:nentries
    eidx = lqn.eshift + e;
    visited = false(nentries, 1);
    [il_all, il_ph1] = traceInterlockPaths(lqn, eidx, e, 1.0, 1.0, visited, il_all, il_ph1, 0);
end

self.il_table_all = il_all;
self.il_table_ph1 = il_ph1;

% Phase B+C: common parents of Fig. 4.2 and the source count n_s
nslots = lqn.nhosts + lqn.ntasks;
self.il_common_entries = cell(lqn.tshift + lqn.ntasks, 1);
self.il_source_tasks_all = cell(lqn.tshift + lqn.ntasks, 1);
self.il_source_tasks_ph2 = cell(lqn.tshift + lqn.ntasks, 1);
self.il_num_sources = zeros(lqn.tshift + lqn.ntasks, 1);

% Process task servers
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if lqn.isref(tidx) || lqn.sched(tidx) == SchedStrategy.INF
        continue;
    end
    [ce, sa, sp, ns] = findInterlockForServer(lqn, tidx, il_all, il_ph1);
    self.il_common_entries{tidx} = ce;
    self.il_source_tasks_all{tidx} = sa;
    self.il_source_tasks_ph2{tidx} = sp;
    self.il_num_sources(tidx) = ns;
end

% Process host servers
for h = 1:lqn.nhosts
    hidx = h;
    if lqn.sched(hidx) == SchedStrategy.INF
        continue;
    end
    [ce, sa, sp, ns] = findInterlockForServer(lqn, hidx, il_all, il_ph1);
    self.il_common_entries{hidx} = ce;
    self.il_source_tasks_all{hidx} = sa;
    self.il_source_tasks_ph2{hidx} = sp;
    self.il_num_sources(hidx) = ns;
end

end

%% Phase A: depth-first path tracing that fills the path table
function [il_all, il_ph1] = traceInterlockPaths(lqn, eidx, root_e, prob_all, prob_ph1, visited, il_all, il_ph1, depth)
e = eidx - lqn.eshift;
if e < 1 || e > lqn.nentries
    return;
end
if visited(e)
    return;
end
visited(e) = true;

% Accumulate path(root,e); the diagonal is unity by definition
il_all(root_e, e) = il_all(root_e, e) + prob_all;
il_ph1(root_e, e) = il_ph1(root_e, e) + prob_ph1;

% Follow synchronous calls from activities of this entry
acts = lqn.actsof{eidx};
for aidx = acts(:)'
    if aidx <= lqn.ashift || aidx > lqn.ashift + lqn.nacts
        continue;
    end
    a = aidx - lqn.ashift;

    % Deferred work is only traced at the root, so a nested phase 2 ends the path
    if depth > 0 && isfield(lqn, 'actphase') && lqn.actphase(a) > 1
        continue;
    end

    is_ph1 = true;
    if isfield(lqn, 'actphase') && lqn.actphase(a) > 1
        is_ph1 = false;
    end

    % Follow calls from this activity
    if aidx <= length(lqn.callsof)
        calls_from_act = lqn.callsof{aidx};
    else
        calls_from_act = [];
    end
    for cidx = calls_from_act(:)'
        if cidx < 1 || cidx > lqn.ncalls
            continue;
        end
        if lqn.calltype(cidx) ~= CallType.SYNC
            continue;
        end
        call_mean = lqn.callproc_mean(cidx);
        if call_mean <= 0
            continue;
        end
        dst_eidx = lqn.callpair(cidx, 2);
        dst_e = dst_eidx - lqn.eshift;
        if dst_e < 1 || dst_e > lqn.nentries
            continue;
        end

        next_all = prob_all * call_mean;
        if is_ph1
            next_ph1 = prob_ph1 * call_mean;
        else
            next_ph1 = 0;
        end

        [il_all, il_ph1] = traceInterlockPaths(lqn, dst_eidx, root_e, next_all, next_ph1, visited, il_all, il_ph1, depth + 1);
    end
end

visited(e) = false;
end

%% Phase B+C: common parents and sources for one server
function [commonEntries, srcAll, srcPh2, numSources] = findInterlockForServer(lqn, serverIdx, il_all, il_ph1)
commonEntries = [];
srcAll = [];
srcPh2 = [];
numSources = 0;

% Get server entry numbers
serverEntryNums = getServerEntryNums(lqn, serverIdx);
if isempty(serverEntryNums)
    return;
end

% Get client tasks
clientTasks = getClientTasks(lqn, serverIdx);
if length(clientTasks) < 1
    return;
end

% Get client entries that reach the server
clientEntryPairs = []; % [taskIdx, entryNum]
for ct = clientTasks(:)'
    for ce = lqn.entriesof{ct}(:)'
        ce_num = ce - lqn.eshift;
        if ce_num < 1 || ce_num > lqn.nentries
            continue;
        end
        for se_num = serverEntryNums(:)'
            if il_all(ce_num, se_num) > 0
                clientEntryPairs(end+1, :) = [ct, ce_num]; %#ok<AGROW>
                break;
            end
        end
    end
end

if size(clientEntryPairs, 1) < 2
    return;
end

% Common-parent finder of Fig. 4.2, kept only where the flow splits
commonEntriesSet = [];
nPairs = size(clientEntryPairs, 1);
for i = 1:nPairs
    for j = (i+1):nPairs
        if clientEntryPairs(i, 1) == clientEntryPairs(j, 1)
            continue; % Same task
        end
        entryA_num = clientEntryPairs(i, 2);
        entryC_num = clientEntryPairs(j, 2);

        % Search all tasks for common parents
        for t = 1:lqn.ntasks
            tidx = lqn.tshift + t;
            entries_of_task = lqn.entriesof{tidx};
            for ex = entries_of_task(:)'
                for ey = entries_of_task(:)'
                    ex_num = ex - lqn.eshift;
                    ey_num = ey - lqn.eshift;
                    if ex_num < 1 || ey_num < 1 || ex_num > lqn.nentries || ey_num > lqn.nentries
                        continue;
                    end
                    if il_all(ex_num, entryA_num) > 0 && il_all(ey_num, entryC_num) > 0
                        if isBranchPointCheck(lqn, ex, entryA_num + lqn.eshift, ey, entryC_num + lqn.eshift, il_all)
                            commonEntriesSet(end+1) = ex; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end
end
commonEntriesSet = unique(commonEntriesSet);

if isempty(commonEntriesSet)
    return;
end

% A second pruning pass, dropping a common entry whose own owner already
% lists it as a common entry, needs the tables of every server at once
% and is not applied here

commonEntries = commonEntriesSet;

% Phase C: source tasks feeding the interlocked paths
interlockedTasks = [];
for ce_eidx = commonEntries(:)'
    itasks = findInterlockedTasks(lqn, ce_eidx, serverIdx, il_all);
    interlockedTasks = union(interlockedTasks, itasks);
end

% All source tasks = tasks owning common entries
allSrcTasks = [];
ph2SrcTasks = [];
for ce_eidx = commonEntries(:)'
    owner_tidx = lqn.parent(ce_eidx);
    allSrcTasks(end+1) = owner_tidx; %#ok<AGROW>
    if hasPhase2Activities(lqn, ce_eidx)
        ph2SrcTasks(end+1) = owner_tidx; %#ok<AGROW>
    end
end
allSrcTasks = unique(allSrcTasks);
ph2SrcTasks = unique(ph2SrcTasks);

% Remove interlocked tasks from allSrcTasks
allSrcTasks = setdiff(allSrcTasks, interlockedTasks);

% Ph2 sources: interlocked tasks with phase-2 activities reaching server
ph2SrcTasks = [];
for it = interlockedTasks(:)'
    for ie = lqn.entriesof{it}(:)'
        if hasPhase2Activities(lqn, ie)
            ie_num = ie - lqn.eshift;
            if ie_num >= 1 && ie_num <= lqn.nentries
                for se_num = serverEntryNums(:)'
                    if il_all(ie_num, se_num) - il_ph1(ie_num, se_num) > 0
                        ph2SrcTasks(end+1) = it; %#ok<AGROW>
                        break;
                    end
                end
            end
        end
    end
end
ph2SrcTasks = unique(ph2SrcTasks);

% Add external sources (tasks calling into interlocked paths from outside)
for it = interlockedTasks(:)'
    for ie = lqn.entriesof{it}(:)'
        [calling_idx, ~] = find(lqn.iscaller(:, ie));
        for ci = calling_idx(:)'
            if ci > lqn.tshift && ci <= lqn.tshift + lqn.ntasks
                if ~ismember(ci, interlockedTasks)
                    allSrcTasks(end+1) = ci; %#ok<AGROW>
                end
            end
        end
    end
end
allSrcTasks = unique(allSrcTasks);

% n_s of Eq. (4.7), counting task copies rather than tasks
nsrc = 0;
for st = allSrcTasks(:)'
    nsrc = nsrc + lqn.mult(st);
end

% Deferred flow from an interlocked task would add one further source per
% task, but SolverLN does not model phase-2 concurrency, so it is not counted
dstEntryNums = serverEntryNums;
if false %#ok<BDLOG> % disabled until LN supports phase-2
for it = interlockedTasks(:)'
    for ie = lqn.entriesof{it}(:)'
        ie_num = ie - lqn.eshift;
        if ie_num < 1 || ie_num > lqn.nentries
            continue;
        end
        for de = dstEntryNums(:)'
            ph2_flow = il_all(ie_num, de) - il_ph1(ie_num, de);
            if ph2_flow > 0
                nsrc = nsrc + 1;
                break;
            end
        end
    end
end
end % if false (phase-2 source counting disabled)

srcAll = allSrcTasks;
srcPh2 = ph2SrcTasks;
numSources = nsrc;
end

%% Get server entry numbers
function nums = getServerEntryNums(lqn, serverIdx)
nums = [];
if serverIdx <= lqn.nhosts
    for tidx = lqn.tasksof{serverIdx}(:)'
        for se = lqn.entriesof{tidx}(:)'
            nums(end+1) = se - lqn.eshift; %#ok<AGROW>
        end
    end
else
    for se = lqn.entriesof{serverIdx}(:)'
        nums(end+1) = se - lqn.eshift; %#ok<AGROW>
    end
end
end

%% Get client tasks for a server
function clientTasks = getClientTasks(lqn, serverIdx)
if serverIdx <= lqn.nhosts
    clientTasks = lqn.tasksof{serverIdx};
else
    server_entries = lqn.entriesof{serverIdx};
    clientTasks = [];
    for se = server_entries(:)'
        [calling_idx, ~] = find(lqn.iscaller(:, se));
        for ci = calling_idx(:)'
            if ci > lqn.tshift && ci <= lqn.tshift + lqn.ntasks
                clientTasks(end+1) = ci; %#ok<AGROW>
            end
        end
    end
    clientTasks = unique(clientTasks);
end
end

%% Branch-point test: do the two paths diverge before the common server
function result = isBranchPointCheck(lqn, srcX_eidx, entryA_eidx, srcY_eidx, entryB_eidx, il_all)
entryA_num = entryA_eidx - lqn.eshift;
entryB_num = entryB_eidx - lqn.eshift;
taskA = lqn.parent(entryA_eidx);
taskB = lqn.parent(entryB_eidx);
taskX = lqn.parent(srcX_eidx);

% Multiserver client: if X, A, B same task => not branch point
if taskX == taskA && taskX == taskB
    result = false;
    return;
end

% Quick check: direct call
if srcX_eidx == entryA_eidx || srcY_eidx == entryB_eidx
    result = true;
    return;
end

% Check downstream calls diverge to different tasks
dstTasks_X = getCallDstTasks(lqn, srcX_eidx, entryA_num, il_all);
dstTasks_Y = getCallDstTasks(lqn, srcY_eidx, entryB_num, il_all);

for dx = dstTasks_X(:)'
    for dy = dstTasks_Y(:)'
        if dx ~= dy
            result = true;
            return;
        end
    end
end
result = false;
end

%% Get destination tasks of sync calls from an entry reaching a target
function dstTasks = getCallDstTasks(lqn, src_eidx, target_e_num, il_all)
dstTasks = [];
acts = lqn.actsof{src_eidx};
for aidx = acts(:)'
    if aidx <= lqn.ashift || aidx > lqn.ashift + lqn.nacts
        continue;
    end
    if aidx <= length(lqn.callsof)
        calls = lqn.callsof{aidx};
    else
        calls = [];
    end
    for cidx = calls(:)'
        if cidx < 1 || cidx > lqn.ncalls
            continue;
        end
        if lqn.calltype(cidx) ~= CallType.SYNC
            continue;
        end
        dst_eidx = lqn.callpair(cidx, 2);
        dst_e = dst_eidx - lqn.eshift;
        if dst_e >= 1 && dst_e <= lqn.nentries && il_all(dst_e, target_e_num) > 0
            dstTasks(end+1) = lqn.parent(dst_eidx); %#ok<AGROW>
        end
    end
end
dstTasks = unique(dstTasks);
end

%% Get interlocked tasks on paths from an entry to a server
function itasks = findInterlockedTasks(lqn, src_eidx, serverIdx, il_all)
visited = false(lqn.nentries, 1);
itasks = traceToServerRec(lqn, src_eidx, serverIdx, il_all, visited, [], true);
end

function itasks = traceToServerRec(lqn, eidx, serverIdx, il_all, visited, itasks, isHead)
e = eidx - lqn.eshift;
if e < 1 || e > lqn.nentries || visited(e)
    return;
end

ownerTask = lqn.parent(eidx);

% Check if we reached the server
if ownerTask == serverIdx
    return;
end
if serverIdx <= lqn.nhosts && lqn.parent(ownerTask) == serverIdx
    return;
end

visited(e) = true;

% Follow synchronous calls from ALL phases (unlike traceInterlockPaths,
% interlocked task detection needs to follow phase-2 paths to find
% intermediate tasks that are sources of non-interlocked traffic)
acts = lqn.actsof{eidx};
found = false;
for aidx = acts(:)'
    if aidx <= lqn.ashift || aidx > lqn.ashift + lqn.nacts
        continue;
    end

    if aidx <= length(lqn.callsof)
        calls = lqn.callsof{aidx};
    else
        calls = [];
    end
    for cidx = calls(:)'
        if cidx < 1 || cidx > lqn.ncalls || lqn.calltype(cidx) ~= CallType.SYNC
            continue;
        end
        dst_eidx = lqn.callpair(cidx, 2);
        dst_task = lqn.parent(dst_eidx);

        % Check if destination reaches server
        reachesServer = false;
        if dst_task == serverIdx
            reachesServer = true;
        elseif serverIdx <= lqn.nhosts && lqn.parent(dst_task) == serverIdx
            reachesServer = true;
        else
            dst_e = dst_eidx - lqn.eshift;
            serverEntryNums = getServerEntryNums(lqn, serverIdx);
            for se_num = serverEntryNums(:)'
                if il_all(dst_e, se_num) > 0
                    reachesServer = true;
                    break;
                end
            end
        end

        if reachesServer
            itasks = traceToServerRec(lqn, dst_eidx, serverIdx, il_all, visited, itasks, false);
            found = true;
        end
    end
end

if found && ~isHead
    itasks(end+1) = ownerTask;
    itasks = unique(itasks);
end

visited(e) = false;
end

%% Check if entry has phase-2 activities
function result = hasPhase2Activities(lqn, eidx)
result = false;
if ~isfield(lqn, 'actphase')
    return;
end
for aidx = lqn.actsof{eidx}(:)'
    a = aidx - lqn.ashift;
    if a > 0 && a <= lqn.nacts && lqn.actphase(a) > 1
        result = true;
        return;
    end
end
end
