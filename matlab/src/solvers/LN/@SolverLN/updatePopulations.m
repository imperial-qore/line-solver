function updatePopulations(self, it)
% UPDATEPOPULATIONS Apply the interlock correction to call residence times
%
% The path tables built by INITINTERLOCK are combined here with the current
% iterate to obtain, for each (client, server) pair, the interlocked flow
% of Eq. (4.3) of Franks (1999),
%
%   lambda_IL = sum over common parents of throughput times path(e,a),
%
% and from it the interlock probability, the share of that flow that the
% layer decomposition would otherwise count twice. Eq. (4.7) removes one
% source in n_s from the queue length inside MVA. SolverLN instead applies
% the equivalent correction to the residence times returned by the layer,
%
%   R_adj = S + (1 - prIL) * W,   W = R - S,
%
% which leaves service and utilization untouched and removes only the
% interlocked share of the waiting time. A second pass does the same for
% the processor queueing seen by the tasks of each host, weighted by the
% fraction of host utilization that is itself interlocked.
%
% Called after UPDATEMETRICS, which produces the raw callresidt from the
% layer solutions, and before UPDATETHINKTIMES.
%
% Reference: G. Franks, "Performance Analysis of Distributed Server
% Systems", PhD thesis, Carleton University, 1999, Ch. 4; published as
% G. Franks, "Traffic dependencies in client-server systems and their
% effect on performance prediction", IEEE IPDS, 1995, pp. 24-33.

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

    % Interlock probability for this client and server. This path serves a TASK,
    % not a processor, so the m' rule of Li/lqns leaves the source population
    % alone; the product IR*Pr(IL) reproduces the scalar this branch used before.
    [IRc, PrILc] = computeInterlockProb(self, lqn, client_tidx, server_for_il, false);
    prIL = IRc * PrILc;

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

% Pass 2: same correction on the processor queueing seen at each host.
% Every layer starts the pass without a matrix, so a host that stops being
% interlocked does not keep last iteration's correction alive.
if isempty(self.ilscaling)
    self.ilscaling = cell(length(self.ensemble),1);
end
for e = 1:length(self.ensemble)
    self.ilscaling{e} = [];
    if isfield(self.solvers{e}.options.config,'interlock')
        self.solvers{e}.options.config.interlock = [];
    end
end
for h = 1:lqn.nhosts
    hidx = h;
    if isempty(self.il_common_entries{hidx})
        continue;
    end

    % Li and Franks (2015): IR of Eq. (4) and Pr(IL) per task on this host. The
    % host of a task layer is a PROCESSOR, which is what selects the m' rule.
    host_tasks = lqn.tasksof{hidx}(:)';
    task_prIL = zeros(length(host_tasks), 1);   % IR, Eq. (4)
    task_PrIL = zeros(length(host_tasks), 1);   % Pr(IL)
    task_util = zeros(length(host_tasks), 1);
    for ti = 1:length(host_tasks)
        tidx = host_tasks(ti);
        [task_prIL(ti), task_PrIL(ti)] = computeInterlockProb(self, lqn, tidx, hidx, true);
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

    % When the layer solver carries Eq. (4.7) inside its own MVA, the interlock
    % goes to the layer as a chain-level matrix and the residence times are left
    % untouched. Scaling them here as well would remove the same waiting twice,
    % and would still leave the layer's own THROUGHPUT uncorrected, which is
    % what breaks flow balance across a call: the reported task rate then comes
    % from a cycle time the correction has already shortened elsewhere.
    e = self.idxhash(hidx);
    if ~isnan(e) && layerTakesInterlock(self, e)
        ILmat = buildLayerInterlock(self, e, host_tasks, task_prIL, task_PrIL);
        self.ilscaling{e} = ILmat;
        self.solvers{e}.options.config.interlock = ILmat;
        continue;
    end

    for ti = 1:length(host_tasks)
        if task_prIL(ti) <= GlobalConstants.FineTol
            continue;
        end
        tidx = host_tasks(ti);
        % Weight by the share of host utilization that is interlocked. The rate
        % is the SAME Eq. (5) product IR*Pr(IL) that pass 1 applies to a call and
        % that buildLayerInterlock puts in the layer matrix -- IR alone is a flow
        % SHARE, ~1 whenever a layer has a single common source, and using it
        % here removed the whole processor queueing rather than the interlocked
        % part of it, which broke flow balance across a call: X(callee)/X(caller)
        % came out 0.68 instead of 1 on a unit call.
        effective_prIL = task_prIL(ti) * task_PrIL(ti) * il_fraction;
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

% The entry servt is rescaled only when it was itself assembled from these
% residence times, which is the default path. After the moment3 pass the entry
% servt is the MEAN OF AN APH CONVOLUTION of the activities' own response
% distributions, and a ratio of residence-time sums is not a correction to it:
% applying it multiplied the entry law by the entry's visit ratio and reported a
% service time BELOW that of the single activity the entry contains. The
% residence times keep their correction either way. See BUGS.md BUG-97.
momentLaws = strcmp(self.options.method,'moment3') ...
    && ~isempty(self.momentPassDone) && self.momentPassDone;
for eidx = (lqn.eshift+1):(lqn.eshift+lqn.nentries)
    if entry_servt_old(eidx) > GlobalConstants.FineTol
        ratio = entry_servt_new(eidx) / entry_servt_old(eidx);
        if ~momentLaws
            self.servt(eidx) = self.servt(eidx) * ratio;
            if self.servt(eidx) > 0
                self.servtproc{eidx} = Exp.fitMean(self.servt(eidx));
            end
        end
        self.residt(eidx) = self.residt(eidx) * ratio;
    end
end
end

%% True when the layer solver applies Eq. (4.7) inside its own MVA
function tf = layerTakesInterlock(self, e)
% Only the MVA layer solver reads options.config.interlock, and only a layer
% whose sole queueing stations are the host's own tasks can take a matrix that
% was built for that host: under flat layering one layer holds every server, so
% the correction stays on the residence times there.
% A layer whose MVA path has no interlock term would be moved to another
% algorithm by the matrix alone: exact multiserver MVA would become AMVA, the
% linearizer would become the load-dependent forward step. That swap is worth
% far more than the correction it carries, and on a layer sitting near a
% bifurcation it turns the LN iteration into a limit cycle. Such a layer keeps
% the residt scaling instead.
tf = isa(self.solvers{e}, 'SolverMVA') && ~isFlatLayering(self) && ...
    mva_carries_interlock(self.ensemble{e}.getStruct(), self.solvers{e}.options);
end

function tf = isFlatLayering(self)
tf = isfield(self.options.config,'layering') && ...
    any(strcmpi(self.options.config.layering, {'flat','squashed'}));
end

%% Class-level interlock matrix of one host layer
function IL = buildLayerInterlock(self, e, host_tasks, task_prIL, task_PrIL)
% IL(r,s) is the share of the class-s queue that a class-r arrival must not see
% at the host. The matrix is CLASS-indexed, not chain-indexed, so that a later
% refreshChains cannot leave it stale; the layer solver aggregates it to chains
% against the struct it is about to solve. Two classes are interlocked only if
% BOTH their tasks are, which is the 0/1 relation ir_mkj of Eq. (5); the diagonal
% stays zero, since a request always sees its own class in full. The entry is the
% Eq. (5) product Pr(IL_ms)*IR_ms*IR_mr, asymmetric in (r,s) because Pr(IL) is
% taken from the QUEUED class s, so that the layer's ILw(r,s) = 1-IL(r,s) is the
% lower-level adjustment rate r_lower.
nclasses = length(self.ensemble{e}.classes);
IL = zeros(nclasses,nclasses);

class_prIL = zeros(nclasses,1);   % IR
class_PrIL = zeros(nclasses,1);   % Pr(IL)
for r = 1:nclasses
    tidx = clientTaskOfClass(self, e, r);
    if isnan(tidx)
        continue
    end
    ti = find(host_tasks == tidx, 1);
    if ~isempty(ti)
        class_prIL(r) = task_prIL(ti);
        class_PrIL(r) = task_PrIL(ti);
    end
end

for r = 1:nclasses
    if class_prIL(r) <= GlobalConstants.FineTol
        continue
    end
    for s = 1:nclasses
        if s == r || class_prIL(s) <= GlobalConstants.FineTol
            continue
        end
        IL(r,s) = class_PrIL(s) * class_prIL(s) * class_prIL(r);
    end
end

if ~any(IL(:) > GlobalConstants.FineTol)
    IL = []; % nothing interlocked, keep the layer on the plain MVA path
end
end

%% Task that a layer class belongs to, NaN when the class names no task
function tidx = clientTaskOfClass(self, e, r)
lqn = self.lqn;
attr = self.ensemble{e}.classes{r}.attribute;
tidx = NaN;
if isempty(attr)
    return
end
switch attr(1)
    case LayeredNetworkElement.TASK
        tidx = attr(2);
    case LayeredNetworkElement.ENTRY
        tidx = lqn.parent(attr(2));
    case LayeredNetworkElement.ACTIVITY
        tidx = lqn.parent(attr(2));
    case LayeredNetworkElement.CALL
        tidx = lqn.parent(lqn.callpair(attr(2),1));
end
if ~isnan(tidx) && (tidx <= lqn.tshift || tidx > lqn.tshift + lqn.ntasks)
    tidx = NaN;
end
end

%% Interlock probability for one (client, server) pair
function [IR, prIL] = computeInterlockProb(self, lqn, client_tidx, server_idx, isProcessorHost)
commonEntries = self.il_common_entries{server_idx};
numSources = self.il_num_sources(server_idx);
allSrcTasks = self.il_source_tasks_all{server_idx};
ph2SrcTasks = self.il_source_tasks_ph2{server_idx};

if numSources == 0 || isempty(commonEntries)
    IR = 0;
    prIL = 0;
    return;
end

% Get client entries
client_entries = lqn.entriesof{client_tidx};

% Interlocked flow lambda^IL of Li and Franks (2015), Eq. (4), and alongside it
% the flow weighted by 1/m', which gives Pr(IL) of Eq. (3). m' is derived from
% the COMMON SOURCE population m the way lqns does it in
% Interlock::ilrate_pril_flow: at a PROCESSOR the source population is doubled
% above 3 customers and squared at or below it, which is what turns m = 4 into
% the pril = 1/8 its trace reports. The rationale is that two customers of the
% same source contend for the processor through TWO intermediate tasks, so the
% pairwise chance of self-contention falls faster than 1/m.
sum_flow = 0;
sum_pril = 0;
for ce_eidx = commonEntries(:)'
    srcTask = lqn.parent(ce_eidx);
    ce_num = ce_eidx - lqn.eshift;
    % population of this common source, in customer copies
    m_src = lqn.mult(srcTask);
    if ~isfinite(m_src) || m_src < 1
        m_src = 1;
    end
    if isProcessorHost
        if m_src > 3
            m_eff = m_src + m_src;
        else
            m_eff = m_src * m_src;
        end
    else
        m_eff = m_src;
    end

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

        % Deferred flow is scored separately from the phase-1 flow
        hasP2 = hasPhase2Check(lqn, ce_eidx);

        if ~hasP2 && ismember(srcTask, allSrcTasks)
            contrib = ce_tput * self.il_table_all(ce_num, dstA_num);
            sum_flow = sum_flow + contrib;
            sum_pril = sum_pril + contrib / m_eff;
        elseif hasP2 && ismember(srcTask, allSrcTasks)
            contrib = ce_tput * self.il_table_ph1(ce_num, dstA_num);
            sum_flow = sum_flow + contrib;
            sum_pril = sum_pril + contrib / m_eff;
        end

        ph2 = self.il_table_all(ce_num, dstA_num) - self.il_table_ph1(ce_num, dstA_num);
        if ph2 > 0 && ismember(srcTask, ph2SrcTasks)
            contrib = ce_tput * ph2;
            sum_flow = sum_flow + contrib;
            sum_pril = sum_pril + contrib / m_eff;
        end
    end
end

% Get client throughput
client_tput = getTaskTput(self, lqn, client_tidx);

if client_tput <= GlobalConstants.FineTol
    IR = 0;
    prIL = 0;
    return;
end

% Li and Franks (2015): IR is Eq. (4), the share of this chain's flow that is
% interlocked; Pr(IL) is the 1/m'-weighted mean of that same flow. The two are
% multiplied into the Eq. (5) rate by buildLayerInterlock, so neither carries the
% source count on its own -- that lives in m'. This replaces the superseded
% (n_s-1)/n_s discount of Franks (1999), Eq. (4.7).
IR = min(sum_flow, client_tput) / client_tput;
IR = min(1.0, max(0.0, IR));
if sum_flow <= GlobalConstants.FineTol
    prIL = 0;
else
    prIL = sum_pril / sum_flow;
end
prIL = min(1.0, max(0.0, prIL));
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
