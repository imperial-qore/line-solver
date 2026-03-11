function updatePopulations(self, it)
% UPDATEPOPULATIONS Apply interlock correction to call residence times
%
% Uses LQNS V5-style interlock analysis: computes interlock probability
% for each (client, server) pair and reduces the waiting time component
% of call residence times proportionally.
%
% LQNS applies interlock as arrival-rate scaling inside MVA:
%   L_k = (1 - prIL_k) * X_k * R_k
% affecting only waiting time, not utilization. LINE approximates this
% by adjusting callresidt post-MVA to remove interlocked waiting:
%   R_adj = S + (1 - prIL) * W
% where S = service time, W = waiting time = R - S.
%
% This function is called AFTER updateMetrics (which computes raw
% callresidt from layer MVA results) and BEFORE updateThinkTimes.

lqn = self.lqn;

% Save originals for proportional entry_servt update
callresidt_orig = self.callresidt(:);
residt_orig = self.residt(:);
adjusted = false;

% For each sync call, check if destination server has interlock
for cidx = 1:lqn.ncalls
    if lqn.calltype(cidx) ~= CallType.SYNC
        continue;
    end

    dst_eidx = lqn.callpair(cidx, 2);
    server_tidx = lqn.parent(dst_eidx);

    % Find the server entity with interlock data
    server_for_il = [];
    if server_tidx <= length(self.il_common_entries) && ~isempty(self.il_common_entries{server_tidx})
        server_for_il = server_tidx;
    else
        % Check host server
        if server_tidx > lqn.tshift
            host_idx = lqn.parent(server_tidx);
            if host_idx >= 1 && host_idx <= length(self.il_common_entries) && ~isempty(self.il_common_entries{host_idx})
                server_for_il = host_idx;
            end
        end
    end
    if isempty(server_for_il)
        continue;
    end

    % Get client task (activity -> task via parent)
    src_aidx = lqn.callpair(cidx, 1);
    client_tidx = lqn.parent(src_aidx);

    % Compute prIL using interlockedFlow formula
    prIL = computeInterlockProb(self, lqn, client_tidx, server_for_il);

    if prIL <= GlobalConstants.FineTol
        continue;
    end

    % Compute waiting time reduction
    S = self.servt(dst_eidx);  % service time at destination entry
    call_mean = lqn.callproc_mean(cidx);

    if call_mean <= 0 || self.callservt(cidx) <= 0
        continue;
    end

    RN = self.callservt(cidx) / call_mean;  % response time per visit
    W = max(0, RN - S);  % waiting time per visit

    if W > GlobalConstants.FineTol
        RN_adj = S + (1 - prIL) * W;
        scale = RN_adj / RN;
        self.callservt(cidx) = self.callservt(cidx) * scale;
        self.callresidt(cidx) = self.callresidt(cidx) * scale;
        if self.callservt(cidx) > 0
            self.callservtproc{cidx} = Exp.fitMean(self.callservt(cidx));
        end
        adjusted = true;
    end
end

% Pass 2: Host-level interlock — reduce processor queueing in residt
for h = 1:lqn.nhosts
    hidx = h;
    if isempty(self.il_common_entries{hidx})
        continue;
    end

    % Compute prIL and processor utilization for each task on this host
    host_tasks = lqn.tasksof{hidx}(:)';
    task_prIL = zeros(length(host_tasks), 1);
    task_util = zeros(length(host_tasks), 1);
    for ti = 1:length(host_tasks)
        tidx = host_tasks(ti);
        task_prIL(ti) = computeInterlockProb(self, lqn, tidx, hidx);
        % Compute task's processor utilization
        for eidx = lqn.entriesof{tidx}(:)'
            for aidx = lqn.actsof{eidx}(:)'
                task_util(ti) = task_util(ti) + self.tput(aidx) * lqn.hostdem_mean(aidx);
            end
        end
    end

    U_total = sum(task_util);
    U_interlocked = sum(task_util(task_prIL > GlobalConstants.FineTol));
    if U_total <= GlobalConstants.FineTol || U_interlocked <= GlobalConstants.FineTol
        continue;
    end
    il_fraction = U_interlocked / U_total;

    for ti = 1:length(host_tasks)
        if task_prIL(ti) <= GlobalConstants.FineTol
            continue;
        end
        tidx = host_tasks(ti);
        % Scale prIL by fraction of utilization that is interlocked
        effective_prIL = task_prIL(ti) * il_fraction;
        for eidx = lqn.entriesof{tidx}(:)'
            for aidx = lqn.actsof{eidx}(:)'
                D = lqn.hostdem_mean(aidx);
                if D > 0 && self.residt(aidx) > D + GlobalConstants.FineTol
                    W_proc = self.residt(aidx) - D;
                    self.residt(aidx) = D + (1 - effective_prIL) * W_proc;
                    adjusted = true;
                end
            end
        end
    end
end

if ~adjusted
    return;
end

% Recompute entry service times from adjusted callresidt/residt
% Use proportional scaling to preserve visit ratio adjustments
% applied in updateMetricsDefault
entry_servt_old = self.servtmatrix * [residt_orig; callresidt_orig];
entry_servt_new = self.servtmatrix * [self.residt; self.callresidt(:)];

for eidx = (lqn.eshift+1):(lqn.eshift+lqn.nentries)
    if entry_servt_old(eidx) > GlobalConstants.FineTol
        ratio = entry_servt_new(eidx) / entry_servt_old(eidx);
        self.servt(eidx) = self.servt(eidx) * ratio;
        self.residt(eidx) = self.residt(eidx) * ratio;
        if self.servt(eidx) > 0
            self.servtproc{eidx} = Exp.fitMean(self.servt(eidx));
        end
    end
end
end

%% Compute interlock probability for a (client, server) pair
function prIL = computeInterlockProb(self, lqn, client_tidx, server_idx)
commonEntries = self.il_common_entries{server_idx};
numSources = self.il_num_sources(server_idx);
allSrcTasks = self.il_source_tasks_all{server_idx};
ph2SrcTasks = self.il_source_tasks_ph2{server_idx};

if numSources == 0 || isempty(commonEntries)
    prIL = 0;
    return;
end

% Get client entries
client_entries = lqn.entriesof{client_tidx};

% Compute interlocked flow (LQNS interlockedFlow formula)
sum_flow = 0;
for ce_eidx = commonEntries(:)'
    srcTask = lqn.parent(ce_eidx);
    ce_num = ce_eidx - lqn.eshift;

    for dstA_eidx = client_entries(:)'
        dstA_num = dstA_eidx - lqn.eshift;

        if dstA_num < 1 || dstA_num > lqn.nentries
            continue;
        end
        if self.il_table_all(ce_num, dstA_num) <= 0
            continue;
        end

        % Get source entry throughput
        ce_tput = getEntryTput(self, lqn, ce_eidx, srcTask);

        if ce_tput <= GlobalConstants.FineTol
            continue;
        end

        % Check maxPhase for the source entry
        hasP2 = hasPhase2Check(lqn, ce_eidx);

        if ~hasP2 && ismember(srcTask, allSrcTasks)
            sum_flow = sum_flow + ce_tput * self.il_table_all(ce_num, dstA_num);
        elseif hasP2 && ismember(srcTask, allSrcTasks)
            sum_flow = sum_flow + ce_tput * self.il_table_ph1(ce_num, dstA_num);
        end

        ph2 = self.il_table_all(ce_num, dstA_num) - self.il_table_ph1(ce_num, dstA_num);
        if ph2 > 0 && ismember(srcTask, ph2SrcTasks)
            sum_flow = sum_flow + ce_tput * ph2;
        end
    end
end

% Get client throughput
client_tput = getTaskTput(self, lqn, client_tidx);

if client_tput <= GlobalConstants.FineTol
    prIL = 0;
    return;
end

client_threads = lqn.mult(client_tidx);
IL = min(sum_flow, client_tput) / (client_tput * client_threads * numSources);
prIL = IL / lqn.mult(server_idx);
prIL = min(prIL, 1.0);
end

%% Helper: get entry throughput
function tput = getEntryTput(self, lqn, eidx, taskIdx)
tput = self.tput(eidx);
if tput <= GlobalConstants.FineTol
    % Try first activity
    acts = lqn.actsof{eidx};
    if ~isempty(acts)
        tput = self.tput(acts(1));
    end
end
if tput <= GlobalConstants.FineTol
    tput = self.tput(taskIdx);
end
end

%% Helper: get task throughput
function tput = getTaskTput(self, lqn, tidx)
tput = self.tput(tidx);
if tput <= GlobalConstants.FineTol
    for eidx = lqn.entriesof{tidx}(:)'
        et = self.tput(eidx);
        if et <= GlobalConstants.FineTol
            acts = lqn.actsof{eidx};
            if ~isempty(acts)
                et = self.tput(acts(1));
            end
        end
        tput = tput + et;
    end
end
end

%% Helper: check if entry has phase-2 activities
function result = hasPhase2Check(lqn, eidx)
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
