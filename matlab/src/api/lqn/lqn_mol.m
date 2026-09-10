%{
 % @brief Method of Layers on the SRVN decomposition of a layered queueing
 %        network whose entries carry no activity graph.
 %
 % @details
 % A compact, self-contained reimplementation of the layered fixed point that
 % @SolverLN runs, restricted to LQNs in which every entry binds exactly one
 % activity and there are no activity precedences. It decomposes the model the
 % way lqns --srvn-layering does -- one submodel per processor and one per
 % called task -- and sweeps them in the two phases of Rolia-Sevcik's Method of
 % Layers: all software (task) submodels, then all hardware (processor) ones.
 %
 % Every submodel is a closed multiclass queueing network with ONE station and
 % one class per client task, so it is solved by pfqn_qdamva rather than by
 % building a Network object. The surrogate client delay of @SolverLN collapses
 % into the think-time vector Z of that call.
 %
 % Where this differs from @SolverLN's srvn.cs: a submodel here carries one
 % class per client TASK with visit-weighted demands, where srvn.cs carries one
 % class per activity and encodes the call multiplicities as routing. On an
 % entry-only model the two agree on the structure and differ only in the
 % aggregation, so the throughputs and processor utilizations track closely
 % while entry response times spread more.
 %
 % @par Scope:
 % Entry-only models. Activity graphs (fork/join, OR-branches, loops, second
 % phases, forwarding), asynchronous calls, caches, setup tasks, admission
 % constraints, replication and open arrivals are REFUSED, not approximated.
 %
 % @par Syntax:
 % @code
 % [QN,UN,RN,TN] = lqn_mol(lsn)
 % [QN,UN,RN,TN,info] = lqn_mol(lsn, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lsn<td>LayeredNetworkStruct (from LayeredNetwork.getStruct())
 % <tr><td>options<td>optional struct: iter_max (200), iter_tol (1e-6),
 %                    relax_factor (0.5)
 % </table>
 %
 % @par Returns:
 % Four (nidx x 1) vectors in the column convention @SolverLN and LQNS report,
 % so they line up with LN(model).getAvgTable cell for cell:
 % <table>
 % <tr><th>Index<th>QN (QLen)<th>UN (Util)<th>RN (RespT)<th>TN (Tput)
 % <tr><td>host<td>NaN<td>processor utilization<td>NaN<td>NaN
 % <tr><td>task<td>sum of entry T*S<td>sum of entry proc util<td>NaN<td>cycle rate
 % <tr><td>entry<td>T*S<td>proc util<td>service time<td>invocation rate
 % <tr><td>activity<td>T*S<td>proc util<td>service time<td>execution rate
 % </table>
 % info carries iter, resid, and the converged servt/residt/callservt/thinkt.
%}
function [QN,UN,RN,TN,info] = lqn_mol(lsn, options)
% [QN,UN,RN,TN,INFO] = LQN_MOL(LSN, OPTIONS)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(options)
    options = struct();
end
if ~isfield(options,'iter_max'),     options.iter_max = 200;    end
if ~isfield(options,'iter_tol'),     options.iter_tol = 1e-6;   end
if ~isfield(options,'relax_factor'), options.relax_factor = 0.5; end

lqn_mol_assert(lsn);

nidx = lsn.nidx;
eidxs = lsn.eshift + (1:lsn.nentries);
tidxs = lsn.tshift + (1:lsn.ntasks);

%% static per-entry data: bound activity, host demand, owning task and host
actof  = zeros(nidx,1);   % entry -> its single activity
dem    = zeros(nidx,1);   % entry -> host demand of that activity
taskof = zeros(nidx,1);   % entry -> task
hostof = zeros(nidx,1);   % task  -> host
for tidx = tidxs
    hostof(tidx) = lsn.parent(tidx);
end
for eidx = eidxs
    aidx = lsn.actsof{eidx};
    actof(eidx) = aidx(1);
    d = lsn.hostdem_mean(actof(eidx));
    if isnan(d), d = 0; end
    dem(eidx) = d;
    taskof(eidx) = lsn.parent(eidx);
end

%% static per-call data
ncalls = lsn.ncalls;
callsrc = zeros(ncalls,1);   % call -> calling ENTRY (its activity's entry)
calldst = zeros(ncalls,1);   % call -> called entry
cally   = zeros(ncalls,1);   % call -> mean number of calls per invocation
entryOfAct = zeros(nidx,1);
for eidx = eidxs
    entryOfAct(actof(eidx)) = eidx;
end
for cidx = 1:ncalls
    callsrc(cidx) = entryOfAct(lsn.callpair(cidx,1));
    calldst(cidx) = lsn.callpair(cidx,2);
    y = lsn.callproc_mean(cidx);
    if isnan(y), y = 0; end
    cally(cidx) = y;
end
% calls issued by each entry, and calls arriving at each entry
callsFrom = cell(nidx,1);
callsTo   = cell(nidx,1);
for idx = 1:nidx, callsFrom{idx} = []; callsTo{idx} = []; end
for cidx = 1:ncalls
    callsFrom{callsrc(cidx)}(end+1) = cidx;
    callsTo{calldst(cidx)}(end+1) = cidx;
end

%% populations and server counts, from maxmult (mult is wrong for INF tasks)
npop = ones(nidx,1);
for idx = 1:(lsn.tshift + lsn.ntasks)
    m = lsn.maxmult(idx);
    if ~isfinite(m) || m < 1, m = 1; end
    npop(idx) = m;
end

%% layer sets: one hardware layer per populated host, one software layer per
%% called non-reference task (buildLayers.m:87-118)
hostLayers = [];
for hidx = 1:lsn.nhosts
    if ~isempty(lsn.tasksof{hidx})
        hostLayers(end+1) = hidx; %#ok<AGROW>
    end
end
taskLayers = [];
isCalled = false(nidx,1);
for cidx = 1:ncalls
    isCalled(taskof(calldst(cidx))) = true;
end
for tidx = tidxs
    if ~lsn.isref(tidx) && isCalled(tidx)
        taskLayers(end+1) = tidx; %#ok<AGROW>
    end
end

%% fixed-point state
residt    = dem;              % entry -> processor residence of its activity
servt     = zeros(nidx,1);    % entry -> response time seen by a caller
callservt = zeros(ncalls,1);  % call  -> blocking time per call
thinkt    = zeros(nidx,1);    % task  -> surrogate idle time
share     = zeros(nidx,1);    % entry -> fraction of its task's invocations
Xtask     = zeros(nidx,1);
Xentry    = zeros(nidx,1);
busyth    = zeros(nidx,1);    % task  -> mean number of busy threads
zref      = zeros(nidx,1);
for tidx = tidxs
    zref(tidx) = lqn_ref_thinktime(lsn, tidx);
    thinkt(tidx) = zref(tidx);
    es = lsn.entriesof{tidx};
    if ~isempty(es)
        share(es) = 1 / numel(es);
    end
end
% seed servt bottom-up over the call graph so a callee is priced before its
% caller; a cycle in the call graph just leaves the residual demand seeded at 0
servt(eidxs) = dem(eidxs);
for pass = 1:max(1,lsn.nentries)
    for eidx = eidxs
        s = residt(eidx);
        for cidx = callsFrom{eidx}
            s = s + cally(cidx) * servt(calldst(cidx));
        end
        servt(eidx) = s;
    end
end
for cidx = 1:ncalls
    callservt(cidx) = servt(calldst(cidx));
end

om = options.relax_factor;
resid = Inf;
iter = 0;
while iter < options.iter_max
    iter = iter + 1;
    servt_prev = servt;
    thinkt_prev = thinkt;

    % ---- phase 1: software layers (thread contention at each called task)
    for tidx = taskLayers
        [gcl, callers] = mol_task_layer(tidx);
        for k = 1:numel(callers)
            ctask = callers(k);
            for eidx = lsn.entriesof{ctask}
                for cidx = callsFrom{eidx}
                    if taskof(calldst(cidx)) == tidx
                        newv = gcl(k) * servt(calldst(cidx));
                        callservt(cidx) = om*newv + (1-om)*callservt(cidx);
                    end
                end
            end
        end
    end

    % ---- phase 2: hardware layers (processor contention at each host)
    for hidx = hostLayers
        [f, tsks] = mol_host_layer(hidx);
        for k = 1:numel(tsks)
            for eidx = lsn.entriesof{tsks(k)}
                residt(eidx) = f(k) * dem(eidx);
            end
        end
    end

    % ---- recompose entry service times
    for eidx = eidxs
        s = residt(eidx);
        for cidx = callsFrom{eidx}
            s = s + cally(cidx) * callservt(cidx);
        end
        servt(eidx) = om*s + (1-om)*servt(eidx);
    end

    % ---- throughputs, entry shares, think-time closure
    mol_throughputs();
    for tidx = tidxs
        if lsn.isref(tidx)
            thinkt(tidx) = zref(tidx);
            continue
        end
        if Xtask(tidx) <= GlobalConstants.FineTol
            continue
        end
        % Idle time of a thread per cycle. updateThinkTimes splits this into an
        % INF arm (njobs - util) and a finite arm (njobs*abs(1-util)) only
        % because LINE reports busy SERVERS at an infinite server and a busy
        % FRACTION at a finite one; carrying the count in both cases makes the
        % two arms the same expression.
        newz = max(0, abs(npop(tidx) - busyth(tidx)) / Xtask(tidx) - zref(tidx));
        thinkt(tidx) = om*newz + (1-om)*thinkt(tidx);
    end

    % both halves of the state must settle: servt alone can sit still for an
    % iteration while the think times are still moving
    resid = max([abs(servt(eidxs) - servt_prev(eidxs)) ./ max(1, abs(servt(eidxs))); ...
                 abs(thinkt(tidxs) - thinkt_prev(tidxs)) ./ max(1, abs(thinkt(tidxs)))]);
    if resid < options.iter_tol
        break
    end
end
mol_throughputs();

%% assemble the reported vectors
QN = nan(nidx,1); UN = nan(nidx,1); RN = nan(nidx,1); TN = nan(nidx,1);
procutil = zeros(nidx,1);
for eidx = eidxs
    hidx = hostof(taskof(eidx));
    procutil(eidx) = Xentry(eidx) * dem(eidx) / mol_host_servers(hidx);
end
for eidx = eidxs
    QN(eidx) = Xentry(eidx) * servt(eidx);
    UN(eidx) = procutil(eidx);
    RN(eidx) = servt(eidx);
    TN(eidx) = Xentry(eidx);
    aidx = actof(eidx);
    QN(aidx) = QN(eidx);
    UN(aidx) = UN(eidx);
    RN(aidx) = RN(eidx);
    TN(aidx) = TN(eidx);
end
for tidx = tidxs
    es = lsn.entriesof{tidx};
    QN(tidx) = sum(QN(es));
    UN(tidx) = sum(UN(es));
    RN(tidx) = NaN;
    TN(tidx) = Xtask(tidx);
end
for hidx = 1:lsn.nhosts
    u = 0;
    for tidx = lsn.tasksof{hidx}
        u = u + UN(tidx);
    end
    QN(hidx) = NaN;
    UN(hidx) = u;
    RN(hidx) = NaN;
    TN(hidx) = NaN;   % no throughput is defined at a processor, as in LQNS
end

info = struct('iter', iter, 'resid', resid, 'servt', servt, 'residt', residt, ...
    'callservt', callservt, 'thinkt', thinkt, 'share', share, ...
    'hostLayers', hostLayers, 'taskLayers', taskLayers);

%% ---------------------------------------------------------------------
    function c = mol_host_servers(hidx)
        % LINE scales a station utilization into [0,1] whatever its
        % multiplicity, and reports busy servers at an infinite server
        if lsn.sched(hidx) == SchedStrategy.INF
            c = 1;
        else
            c = npop(hidx);
        end
    end

    function z = mol_cycle_outside(tidx, excludeTask)
        % Time a thread of TIDX spends away from EXCLUDETASK in one cycle:
        % its think time, its own processor residence, and its blocking at
        % every callee other than EXCLUDETASK
        z = thinkt(tidx);
        for eidx = lsn.entriesof{tidx}
            w = share(eidx) * residt(eidx);
            for cidx = callsFrom{eidx}
                if taskof(calldst(cidx)) ~= excludeTask
                    w = w + share(eidx) * cally(cidx) * callservt(cidx);
                end
            end
            z = z + w;
        end
    end

    function [f, tsks] = mol_host_layer(hidx)
        % Hardware submodel: the processor is the station, its tasks the
        % classes. Returns the queueing inflation f = R/L per task.
        tsks = lsn.tasksof{hidx};
        K = numel(tsks);
        L = zeros(1,K); Z = zeros(1,K); N = zeros(1,K);
        for k = 1:K
            tidx = tsks(k);
            N(k) = npop(tidx);
            d = 0; z = thinkt(tidx);
            for eidx = lsn.entriesof{tidx}
                d = d + share(eidx) * dem(eidx);
                for cidx = callsFrom{eidx}
                    z = z + share(eidx) * cally(cidx) * callservt(cidx);
                end
            end
            L(k) = d; Z(k) = z;
        end
        f = ones(1,K);
        if lsn.sched(hidx) == SchedStrategy.INF
            return % a delay processor never queues
        end
        mu = mol_mu(N, npop(hidx));
        [~,~,~,~,R] = pfqn_qdamva(L, N, Z, mu);
        for k = 1:K
            if L(k) > GlobalConstants.FineTol
                f(k) = R(1,k) / L(k);
            end
        end
    end

    function [g, callers] = mol_task_layer(tidx)
        % Software submodel: the task is the station, its caller tasks the
        % classes. Returns the inflation g = R/L per caller.
        callers = [];
        for eidx = lsn.entriesof{tidx}
            for cidx = callsTo{eidx}
                callers(end+1) = taskof(callsrc(cidx)); %#ok<AGROW>
            end
        end
        callers = unique(callers);
        K = numel(callers);
        g = ones(1,K);
        if K == 0
            return
        end
        L = zeros(1,K); Z = zeros(1,K); N = zeros(1,K);
        for k = 1:K
            ctask = callers(k);
            N(k) = npop(ctask);
            d = 0;
            for eidx = lsn.entriesof{ctask}
                for cidx = callsFrom{eidx}
                    if taskof(calldst(cidx)) == tidx
                        d = d + share(eidx) * cally(cidx) * servt(calldst(cidx));
                    end
                end
            end
            L(k) = d;
            Z(k) = mol_cycle_outside(ctask, tidx);
        end
        if lsn.sched(tidx) == SchedStrategy.INF
            return % an infinite-thread task never queues for a thread
        end
        % The AMVA U output is X*L*g, where g is the RECIPROCAL RATE MULTIPLIER
        % at the current congestion, not 1/c -- it is not a busy-server count,
        % so nothing here reads it. Occupancy comes from Little's law in
        % mol_throughputs.
        [~,~,~,~,R] = pfqn_qdamva(L, N, Z, mol_mu(N, npop(tidx)));
        for k = 1:K
            if L(k) > GlobalConstants.FineTol
                g(k) = R(1,k) / L(k);
            end
        end
    end

    function mol_throughputs()
        % Reference tasks set the pace; every other rate follows from the call
        % rates, so the entries are visited in call-graph order until stable
        for t = tidxs
            if lsn.isref(t)
                cyc = thinkt(t);
                for eidx = lsn.entriesof{t}
                    cyc = cyc + share(eidx) * servt(eidx);
                end
                if cyc > GlobalConstants.FineTol
                    Xtask(t) = npop(t) / cyc;
                else
                    Xtask(t) = 0;
                end
                for eidx = lsn.entriesof{t}
                    Xentry(eidx) = Xtask(t) * share(eidx);
                end
            else
                Xtask(t) = 0;
                Xentry(lsn.entriesof{t}) = 0;
            end
        end
        for pass = 1:max(1,lsn.ntasks)
            for eidx = eidxs
                if lsn.isref(taskof(eidx))
                    continue
                end
                x = 0;
                for cidx = callsTo{eidx}
                    x = x + cally(cidx) * Xentry(callsrc(cidx));
                end
                Xentry(eidx) = x;
            end
            for t = tidxs
                if lsn.isref(t), continue, end
                es = lsn.entriesof{t};
                Xtask(t) = sum(Xentry(es));
                if Xtask(t) > GlobalConstants.FineTol
                    share(es) = Xentry(es) / Xtask(t);
                end
            end
        end
        % Mean busy threads, by Little's law over the entries the task serves.
        % This is the occupancy the think-time closure needs, and it is exact
        % given the throughputs -- unlike the layer AMVA's own U.
        for t = tidxs
            u = 0;
            for eidx = lsn.entriesof{t}
                u = u + Xentry(eidx) * servt(eidx);
            end
            busyth(t) = u;
        end
    end
end

%% -------------------------------------------------------------------------
function mu = mol_mu(N, c)
% MU = MOL_MU(N, C) is the queue-dependent rate multiplier row of a C-server
% station over a population of SUM(N). pfqn_lldfun SKIPS a constant row, so a
% single server must come back as ones and not as a scalar 1.
smax = max(2, ceil(sum(N)));
if ~isfinite(c) || c <= 1
    mu = ones(1, smax);
else
    mu = min(1:smax, c);
end
end

%% -------------------------------------------------------------------------
function lqn_mol_assert(lsn)
% Refuse every feature this decomposition does not represent, naming the
% element, rather than returning a number that quietly ignores it.
for eidx = lsn.eshift + (1:lsn.nentries)
    acts = lsn.actsof{eidx};
    if numel(acts) ~= 1
        line_error(mfilename, sprintf(['Entry %s binds %d activities. lqn_mol solves ' ...
            'entry-only models; use SolverLN for an activity graph.'], ...
            lsn.hashnames{eidx}, numel(acts)));
    end
end
ashift = lsn.ashift;
for a = 1:lsn.nacts
    aidx = ashift + a;
    succ = find(lsn.graph(aidx, ashift + (1:lsn.nacts)));
    if ~isempty(succ)
        line_error(mfilename, sprintf(['Activity %s has an activity precedence. lqn_mol ' ...
            'solves entry-only models; use SolverLN for an activity graph.'], ...
            lsn.hashnames{aidx}));
    end
    if ~isempty(lsn.actphase) && a <= numel(lsn.actphase) && lsn.actphase(a) ~= 1
        line_error(mfilename, sprintf('Activity %s is in phase %d. lqn_mol supports phase 1 only.', ...
            lsn.hashnames{aidx}, lsn.actphase(a)));
    end
end
for cidx = 1:lsn.ncalls
    if lsn.calltype(cidx) ~= CallType.SYNC
        line_error(mfilename, sprintf('Call %s is %s. lqn_mol supports synchronous calls only.', ...
            lsn.callhashnames{cidx}, CallType.toText(lsn.calltype(cidx))));
    end
end
for idx = 1:(lsn.tshift + lsn.ntasks)
    if ~isempty(lsn.iscache) && idx <= numel(lsn.iscache) && lsn.iscache(idx)
        line_error(mfilename, sprintf('%s is a cache task, which lqn_mol does not model.', lsn.hashnames{idx}));
    end
    if ~isempty(lsn.hassetup) && idx <= numel(lsn.hassetup) && lsn.hassetup(idx)
        line_error(mfilename, sprintf('%s has a setup time, which lqn_mol does not model.', lsn.hashnames{idx}));
    end
    if lsn.repl(idx) ~= 1
        line_error(mfilename, sprintf('%s is replicated %d times, which lqn_mol does not model.', ...
            lsn.hashnames{idx}, lsn.repl(idx)));
    end
    if ~isempty(lsn.lincon) && idx <= size(lsn.lincon,1) && ~isempty(lsn.lincon{idx,1})
        line_error(mfilename, sprintf('%s carries an admission constraint, which lqn_mol does not model.', lsn.hashnames{idx}));
    end
    s = lsn.sched(idx);
    if idx <= lsn.nhosts
        ok = ismember(s, [SchedStrategy.PS, SchedStrategy.FCFS, SchedStrategy.INF]);
    else
        ok = ismember(s, [SchedStrategy.PS, SchedStrategy.FCFS, SchedStrategy.INF, SchedStrategy.REF]);
    end
    if ~ok
        line_error(mfilename, sprintf('%s is scheduled %s, which lqn_mol does not model.', ...
            lsn.hashnames{idx}, SchedStrategy.toText(s)));
    end
end
if isfield(lsn,'callgroups') && ~isempty(lsn.callgroups)
    line_error(mfilename, 'This model uses routed call groups, which lqn_mol does not model.');
end
for eidx = lsn.eshift + (1:lsn.nentries)
    if ~isempty(lsn.arrival) && eidx <= numel(lsn.arrival) && ~isempty(lsn.arrival{eidx})
        line_error(mfilename, sprintf('Entry %s has an open arrival. lqn_mol solves closed models only.', ...
            lsn.hashnames{eidx}));
    end
end
end
