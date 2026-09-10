function [outspace, outrate, outprob, eventCache, outstart, outpreempt] = afterEventStation(sn, ind, inspace, event, class, isSimulation, eventCache, ...
    M, R, S, phasessz, phaseshift, pie, isf, ismkvmod, ismkvmodclass, lldscaling, lldlimit, cdscaling, ...
    hasOnlyExp, ist, K, Ks, mu, phi, proc, capacity, classcap, V, space_buf, space_srv, space_var, key, noPromote)
% NOPROMOTE (optional, default false): when true, a DEP at an FCFS-family
% station does not promote the head-of-buffer job into the vacated server.
% It is set only for the active (departure) half of an immediate-feedback
% self-loop (sn.immfeed), so the job that self-loops holds the server and the
% subsequent passive arrival re-enters service instead of re-queueing behind
% the waiting jobs. See State.afterEvent and solver_ssa.
if nargin < 34 || isempty(noPromote)
    noPromote = false;
end
% A retrial station keeps blocked arrivals in an orbit (represented by the
% buffer slots) instead of an ordered waiting line: on service completion the
% freed server is NOT filled from the orbit (no promotion); orbiting jobs
% re-enter only through RETRY events at the memoryless retrial rate.
isRetrialStation = isfield(sn,'retrialProc') && ~isempty(sn.retrialProc) ...
    && ist > 0 && any(~cellfun(@isempty, sn.retrialProc(ist,:)));
outspace = [];
outrate = [];
outprob = 1;
% START/PREEMPT annotations, one row per successor and one column per class.
% They are instantaneous tags on the arcs built below, never events in their
% own right: nothing here reads them back, so no rate, probability or state
% depends on them. Sites that start or preempt a job write their rows with
% State.tagArc; every other arc stays a zero row, filled in by State.tagPad
% in the caller once the successor list is final.
outstart = [];
outpreempt = [];
% Pass-and-swap / order-independent stations use a dedicated ordered-list
% representation and rate function mu(c); handle them separately.
if sn.sched(ist) == SchedStrategy.PAS
    [outspace, outrate, outprob, eventCache, outstart, outpreempt] = State.afterEventStationPAS(sn, ind, ist, inspace, event, class, isSimulation, eventCache, R, V, key);
    return;
end
% Server breakdown. The status is the trailing local-variable column (0 = down,
% 1 = up), exclusive with the BAS marker and the polling controller. A down
% server does not serve, so DEP and PHASE are suppressed unless a degraded
% down-server rate was configured, in which case the completion rate is rescaled
% by downRate/upRate at the end of the handler. Gating here keeps every
% scheduling branch below unaware of the server status.
isBreakdownStation = isfield(sn,'hasbreakdown') && ~isempty(sn.hasbreakdown) ...
    && numel(sn.hasbreakdown) >= ind && sn.hasbreakdown(ind) == 1;
downRateScale = 1;
if isBreakdownStation && ~isempty(inspace) && (event == EventType.DEP || event == EventType.PHASE)
    if inspace(1,end) == 0 % server down
        downRate = 0;
        if ~isempty(sn.downServiceRates) && size(sn.downServiceRates,1) >= ist ...
                && size(sn.downServiceRates,2) >= class
            downRate = sn.downServiceRates(ist,class);
        end
        if downRate <= 0
            outspace = [];
            outrate = [];
            outprob = [];
            return
        end
        upRate = sn.rates(ist,class);
        if ~isfinite(upRate) || upRate <= 0
            line_error(mfilename, sprintf(['Station ''%s'' declares a down-server service rate for class ''%s'' but has ' ...
                'no finite up-server service rate to rescale.'], sn.nodenames{ind}, sn.classnames{class}));
        end
        downRateScale = downRate / upRate;
    end
end
switch event
    case EventType.FAILURE
        % An up server fails at the memoryless rate breakdownMu, whether or not
        % it is serving. Only the status column changes: jobs in service are not
        % lost and, service being memoryless here, they resume on repair. The
        % passive half of the synchronization is LOCAL, so no job moves.
        outspace = [];
        outrate = [];
        outprob = [];
        if isBreakdownStation && ~isempty(inspace) && inspace(1,end) == 1
            outspace = inspace;
            outspace(:,end) = 0;
            outrate = sn.breakdownMu(ist);
            outprob = 1;
        end
    case EventType.REPAIR
        % A down server is restored at the memoryless rate repairMu.
        outspace = [];
        outrate = [];
        outprob = [];
        if isBreakdownStation && ~isempty(inspace) && inspace(1,end) == 0
            outspace = inspace;
            outspace(:,end) = 1;
            outrate = sn.repairMu(ist);
            outprob = 1;
        end
    case EventType.ARV %% passive
        % A Place holds a marking, not a service facility: it has no servers and
        % no service phases, so an arriving token only increments the class
        % marking. The scheduling branches below would instead write an entering
        % job into a phase slot that a Place state does not carry; MATLAB grows
        % the row to fit, and the widened state can no longer be matched against
        % the marking state space, so the successor hash fails and the arrival is
        % dropped. Closed SPNs never notice, because there tokens reach a Place
        % through FIRE (State.afterGlobalEvent) rather than through ARV.
        if sn.nodetype(ind) == NodeType.Place
            outspace = inspace;
            outspace(:,class) = outspace(:,class) + 1;
            % passive action: the rate is set by the active node
            outrate = -1*ones(size(outspace,1),1);
            outprob = ones(size(outspace,1),1);
            return
        end
        % Signal / catastrophe arrival: a signal class is a G-network negative
        % customer. Instead of joining, it removes job(s) of its target class
        % from this station and is itself annihilated. CATASTROPHE empties the
        % station of all jobs. The event stays passive (rate set by the active
        % signal source); when there is no target job the signal simply vanishes
        % (destination state unchanged).
        if isfield(sn,'issignal') && ~isempty(sn.issignal) && sn.issignal(class)
            % A REPLY signal is not a negative customer: it completes a
            % synchronous call, releasing the server this station holds for the
            % caller, and then joins as an ordinary job. Only stations that
            % actually hold a block for it take this path; elsewhere a REPLY
            % class is a plain job class and falls through to the normal
            % arrival handling below.
            if isfield(sn,'signaltype') && ~isempty(sn.signaltype) && numel(sn.signaltype) >= class ...
                    && ~isempty(sn.signaltype{class}) && ~any(isnan(sn.signaltype{class})) ...
                    && sn.signaltype{class} == SignalType.REPLY
                rinfoR = State.replyBlockInfo(sn, ind);
                if rinfoR.width > 0
                    [outspace, outrate, outprob, outstart] = State.afterEventStationReply(sn, ind, ist, class, K, Ks, S, pie, space_buf, space_srv, space_var);
                    if isSimulation && size(outprob,1) > 1
                        cum_prob = cumsum(outprob) / sum(outprob);
                        firing_ctr = 1 + max([0,find( rand > cum_prob' )]);
                        outspace = outspace(firing_ctr,:);
                        outrate = outrate(firing_ctr,:);
                        outstart = outstart(firing_ctr,:);
                        outprob = 1;
                    end
                    return
                end
            else
            [outspace, outrate, outprob, outstart] = State.afterEventStationSignal(sn, ind, ist, inspace, class, K, Ks, S, pie, space_buf, space_srv, space_var);
            % Which job a signal removes is in general a random choice, so the
            % generator needs every destination and its probability. A simulation
            % instead walks one sample path, so it draws a single successor from
            % outprob, as the balking branch below does for its passive action.
            if isSimulation && size(outprob,1) > 1
                cum_prob = cumsum(outprob) / sum(outprob);
                firing_ctr = 1 + max([0,find( rand > cum_prob' )]); % select action
                outspace = outspace(firing_ctr,:);
                outstart = outstart(firing_ctr,:);
                outrate = -1;
                outprob = 1;
            end
            return;
            end
        end
        % Ordinary SPN Place arrival: see _kb/11-conventions-and-gotchas.md
        % ("An ordinary Place must be special-cased before the generic
        % scheduling switch") for the rationale.
        isOrdinaryPlace = isfield(sn,'nodetype') && ~isempty(sn.nodetype) ...
            && sn.nodetype(ind) == NodeType.Place;
        if isOrdinaryPlace && isfield(sn,'isqueueingplace') && ~isempty(sn.isqueueingplace) ...
                && ist <= numel(sn.isqueueingplace) && sn.isqueueingplace(ist)
            isOrdinaryPlace = false;
        end
        if isOrdinaryPlace
            outspace = inspace;
            nrows = size(outspace,1);
            if inspace(1,class) < classcap(ist,class)
                outspace(:,class) = outspace(:,class) + 1;
                outrate = -ones(nrows,1);   % passive: rate set by active source
                outprob = ones(nrows,1);
            else
                % place full: arrival is blocked and lost (no state change)
                outrate = -ones(nrows,1);
                outprob = zeros(nrows,1);
            end
            return;
        end
        % return if there is no space to accept the arrival
        [ni,nir] = State.toMarginalAggr(sn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
        % otherwise check scheduling strategy
        pentry = pie{ist}{class};
        % For Place nodes (INF scheduling with NaN service), use uniform entry probability
        if all(isnan(pentry))
            pentry = ones(size(pentry)) / length(pentry);
        end
        outprob = [];
        outprob_k = [];
        for kentry = 1:K(class)
            space_var_k = space_var;
            space_srv_k = space_srv;
            space_buf_k = space_buf;
            % Per-row tag of the arrival arc: the class that takes a server in
            % this row (start_k) and the class it displaces (preempt_k), 0 for
            % neither. They follow every row selection applied below, so that
            % they can be filtered with en_o at the single append site.
            start_k = zeros(size(space_srv_k,1),1);
            preempt_k = zeros(size(space_srv_k,1),1);
            switch sn.sched(ist)
                case SchedStrategy.EXT % source, can receive any "virtual" arrival from the sink as long as it is from an open class
                    if isinf(sn.njobs(class))
                        outspace = inspace;
                        outrate = -1*zeros(size(outspace,1)); % passive action, rate is unspecified
                        outprob = ones(size(outspace,1));
                        break
                    end
                case {SchedStrategy.PS, SchedStrategy.INF, SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, SchedStrategy.LPS}
                    % job enters service immediately
                    if space_srv_k(:,Ks(class)+kentry) < classcap(ist,class)
                        space_srv_k(:,Ks(class)+kentry) = space_srv_k(:,Ks(class)+kentry) + 1;
                        start_k(:) = class;
                        outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                    else
                        outprob_k = pentry(kentry)*zeros(size(space_srv_k,1));
                    end
                case {SchedStrategy.SIRO, SchedStrategy.SEPT, SchedStrategy.LEPT}
                    % Idle-server test on server occupancy, not total count ni: these
                    % agree for work-conserving states but an immediate-feedback
                    % self-loop transiently yields an idle server with a non-empty
                    % buffer, where the fed-back job must re-enter the vacated server
                    % (mirrors the FCFS case below).
                    if sum(space_srv_k,2)<S(ist)
                        space_srv_k(:,Ks(class)+kentry) = space_srv_k(:,Ks(class)+kentry) + 1;
                        start_k(:) = class;
                        outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                    else
                        space_buf_k(:,class) = space_buf_k(:,class) + 1;
                        outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                    end
                case SchedStrategy.POLLING
                    % The controller, not the arrival, decides who is served:
                    % an arriving job joins its class buffer and waits for the
                    % server to walk to it, even when the facility is idle,
                    % because the server is then in a switchover. The exception
                    % is a parked server, which State.pollingNext only ever
                    % produces with an empty station and immediate switchovers:
                    % it therefore reaches the arriving job in zero time and
                    % starts a visit on it at once.
                    pinfoA = State.pollingInfo(sn, ind);
                    srvclassA = 0;
                    for rA=1:R
                        if sum(space_srv_k(1,(Ks(rA)+1):(Ks(rA)+K(rA)))) > 0
                            srvclassA = rA;
                            break
                        end
                    end
                    [posA, swkA] = State.pollingGet(pinfoA, space_var_k, srvclassA);
                    if srvclassA == 0 && swkA == 0
                        % Parked. The arriving job is the only work present, so
                        % the walk necessarily resolves to a visit on its own
                        % class: it enters service here, in phase kentry, and
                        % never occupies the buffer.
                        nbufA = space_buf_k(1,1:R);
                        nbufA(class) = nbufA(class) + 1;
                        [qA, ~, budgetA] = State.pollingNext(pinfoA, posA, nbufA, R, true);
                        space_srv_k(:,Ks(class)+kentry) = space_srv_k(:,Ks(class)+kentry) + 1;
                        start_k(:) = class;
                        space_var_k = State.pollingSet(pinfoA, space_var_k, qA, 0, budgetA);
                    else
                        space_buf_k(:,class) = space_buf_k(:,class) + 1;
                    end
                    outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS, SchedStrategy.LCFSPRIO}
                    % find states with all servers busy - this
                    % needs not to be moved

                    % if MAP service, when empty restart from the phase
                    % stored in space_var for this class
                    if ~ismkvmodclass(class) || (ismkvmodclass(class) && kentry == space_var(sum(sn.nvars(ind,1:class))))
                        if ismkvmodclass(class)
                            pentry = zeros(size(pentry));
                            pentry(kentry) = 1.0;
                        end
                        % Servers held by synchronous calls awaiting a REPLY are
                        % NOT available to an arriving job: subtract them from the
                        % server count. Zero for every model without reply signals.
                        [~, nbA] = State.replyBlocked(sn, ind, space_var_k);
                        SeffA = S(ist) - nbA;
                        all_busy_srv = sum(space_srv_k,2) >= SeffA;

                        % find and modify states with an idle server
                        idle_srv = sum(space_srv_k,2) < SeffA;
                        space_srv_k(idle_srv, end-sum(K)+Ks(class)+kentry) = space_srv_k(idle_srv,end-sum(K)+Ks(class)+kentry) + 1; % job enters service
                        start_k(idle_srv) = class;

                        % this section dynamically grows the number of
                        % elements in the buffer

                        if any(ni < capacity(ist))
                            if any(nir(:,class) < classcap(ist,class)) % if there is room
                                if ~any(space_buf_k(:)==0) % but the buffer has no empty slots
                                    % append job slot
                                    space_buf_k = [zeros(size(space_buf_k,1),1),space_buf_k];
                                end
                            end
                        end
                        %end
                        %get position of first empty slot
                        empty_slots = -1*ones(size(all_busy_srv,1),1);
                        if size(space_buf_k,2) == 0
                            empty_slots(all_busy_srv) = false;
                        elseif size(space_buf_k,2) == 1
                            empty_slots(all_busy_srv) = space_buf_k(all_busy_srv,:)==0;
                        else
                            empty_slots(all_busy_srv) = max(bsxfun(@times, space_buf_k(all_busy_srv,:)==0, [1:size(space_buf_k,2)]),[],2);
                        end

                        % ignore states where the buffer has no empty slots.
                        % A structurally free buffer column is NOT enough: the job may
                        % only be placed when the CAPACITY also permits it. Without this
                        % the job was written into the slot and the whole row was then
                        % vetoed by the capacity filter below (en_o), which returned an
                        % EMPTY outspace -- i.e. the arrival event never fired, so the
                        % upstream departure was never counted and the reported arrival
                        % rate collapsed from the OFFERED rate to the CARRIED one
                        % (BUG-85: M/M/1/1 reported ArvR 0.444444 instead of 0.800000,
                        % hiding the loss entirely, while M/M/1/2 -- whose buffer is
                        % full rather than absent -- correctly reported 0.800000).
                        % Leaving the state unchanged instead makes the arrival a
                        % self-loop: the event fires, ArvR counts the offered job, and
                        % the job is lost. A self-loop cancels on the generator diagonal,
                        % so the stationary distribution -- and hence QLen/Util/Tput --
                        % cannot move.
                        % Gate the placement on the capacity ONLY for a PHYSICAL
                        % finite capacity. When the bound is a state-space cutoff
                        % (open class, no physical cap; solver_ssa folds the cutoff
                        % into the capacity/classcap arguments), keep the pre-change
                        % behaviour: place the job structurally and let the en_o
                        % capacity filter below delete the beyond-cutoff row
                        % (truncation). Firing the capacity gate at a cutoff would
                        % send the arrival to the loss/block branch, turning a
                        % truncation into a self-loop. See State.isPhysicalCapacity.
                        if State.isPhysicalCapacity(sn, ist, class)
                            hasRoom = (ni < capacity(ist)) & (nir(:,class) < classcap(ist,class));
                        else
                            hasRoom = true(size(empty_slots));
                        end
                        wbuf_empty = empty_slots>0 & hasRoom;
                        if any(wbuf_empty)
                            space_srv_k = space_srv_k(wbuf_empty,:);
                            space_buf_k = space_buf_k(wbuf_empty,:);
                            space_var_k = space_var_k(wbuf_empty,:);
                            start_k = start_k(wbuf_empty);
                            preempt_k = preempt_k(wbuf_empty);
                            empty_slots = empty_slots(wbuf_empty);
                            space_buf_k(sub2ind(size(space_buf_k),1:size(space_buf_k,1),empty_slots')) = class;
                            %outspace(all_busy_srv(wbuf_empty),:) = [space_buf, space_srv, space_var];
                        elseif any(all_busy_srv) && ~State.arrivalIsLost(sn, ist, class)
                            % The arrival cannot be placed (all servers busy, no room).
                            % Whether that is a LOSS or a BLOCK is decided by the class
                            % type, not by the drop rule: a CLOSED job has nowhere to go,
                            % so it cannot be lost, and an explicit blocking rule
                            % (BAS/BBS/RSRD) asks for blocking too. Removing the rows
                            % returns an EMPTY outspace, which disables the upstream
                            % departure until space frees (and is what the CTMC
                            % become-blocked edge tests for, so true BAS fires). For a
                            % LOST arrival the rows are kept unchanged instead, so the
                            % event still fires and ArvR counts the offered job -- see
                            % State.arrivalIsLost and the hasRoom comment above.
                            space_srv_k = space_srv_k(idle_srv,:);
                            space_buf_k = space_buf_k(idle_srv,:);
                            space_var_k = space_var_k(idle_srv,:);
                            start_k = start_k(idle_srv);
                            preempt_k = preempt_k(idle_srv);
                        end
                        outprob_k = pentry(kentry)*ones(size(space_srv_k,1),1);
                    else
                        outprob_k = 0*ones(size(space_srv_k,1),1); % zero probability event
                    end
                case {SchedStrategy.FCFSPR,SchedStrategy.FCFSPI,SchedStrategy.FCFSPRPRIO,SchedStrategy.FCFSPIPRIO,SchedStrategy.LCFSPR,SchedStrategy.LCFSPI,SchedStrategy.LCFSPRPRIO,SchedStrategy.LCFSPIPRIO}
                    % find states with all servers busy - this
                    % must not be moved
                    all_busy_srv = sum(space_srv_k,2) >= S(ist);
                    % find states with an idle server
                    idle_srv = sum(space_srv_k,2) < S(ist);

                    % reorder states so that idle ones come first
                    space_buf_k_reord = space_buf_k(idle_srv,:);
                    space_srv_k_reord = space_srv_k(idle_srv,:);
                    space_var_k_reord = space_var_k(idle_srv,:);
                    start_k_reord = start_k(idle_srv);
                    preempt_k_reord = preempt_k(idle_srv);

                    % if idle, the job enters service in phase kentry
                    if any(idle_srv)
                        space_srv_k_reord(:, end-sum(K)+Ks(class)+kentry) = space_srv_k_reord(:,end-sum(K)+Ks(class)+kentry) + 1;
                        start_k_reord(:) = class;
                        outprob_k = pentry(kentry);
                    else
                        % if all busy, expand output states for all possible choices of job class to preempt
                        psentry = ones(size(space_buf_k_reord,1),1); % probability scaling due to preemption
                        isPrioAware = (sn.sched(ist) == SchedStrategy.FCFSPRPRIO || sn.sched(ist) == SchedStrategy.FCFSPIPRIO || sn.sched(ist) == SchedStrategy.LCFSPRPRIO || sn.sched(ist) == SchedStrategy.LCFSPIPRIO);
                        isLcfsPrioFamily = (sn.sched(ist) == SchedStrategy.LCFSPRPRIO || sn.sched(ist) == SchedStrategy.LCFSPIPRIO);
                        for classpreempt = 1:R
                            % For priority variants (or LCFSPR/FCFSPR when priorities differ),
                            % only higher-priority jobs can preempt
                            isPrioSched = (sn.sched(ist) == SchedStrategy.FCFSPRPRIO || sn.sched(ist) == SchedStrategy.FCFSPIPRIO || sn.sched(ist) == SchedStrategy.LCFSPRPRIO || sn.sched(ist) == SchedStrategy.LCFSPIPRIO);
                            % Priority-awareness is a property of the DECLARED policy, never of
                            % the data. LCFSPRPRIO/FCFSPRPRIO are the priority-aware variants;
                            % inferring the discipline from ~all(classprio==classprio(1)) turned
                            % plain LCFSPR/FCFSPR into something that is neither the base policy
                            % nor the PRIO variant. The FCFSPR case below states the rule outright
                            % ("FCFS preempt-resume (no priority)").
                            isPrioAware = isPrioSched;
                            if isPrioAware
                                % Across priority groups a strictly higher-priority arrival preempts.
                                % WITHIN a group the base discipline decides: LCFS-PR keeps the NEWEST
                                % job in service, so an equal-priority arrival preempts; FCFS-PR never
                                % lets an arrival preempt. Refusing on equality for both left
                                % LCFSPRPRIO unable to preempt anything when the priorities were flat,
                                % so no arrival could enter a busy station at all.
                                isLcfsPrioFamily = (sn.sched(ist) == SchedStrategy.LCFSPRPRIO || sn.sched(ist) == SchedStrategy.LCFSPIPRIO);
                                if isLcfsPrioFamily
                                    cannotPreempt = sn.classprio(class) > sn.classprio(classpreempt);
                                else
                                    cannotPreempt = sn.classprio(class) >= sn.classprio(classpreempt);
                                end
                                if cannotPreempt
                                    continue;
                                end
                            end
                            for phasepreempt = 1:K(classpreempt) % phase of job to preempt
                                si_preempt = space_srv_k(:, (end-sum(K)+Ks(classpreempt)+phasepreempt));
                                busy_preempt = si_preempt > 0; % states where there is at least on class-r job in execution
                                if any(busy_preempt)
                                    psentry = [psentry; si_preempt(busy_preempt) ./ sum(space_srv_k,2)];
                                    space_srv_k_preempt = space_srv_k(busy_preempt,:);
                                    space_buf_k_preempt = space_buf_k(busy_preempt,:);
                                    space_var_k_preempt = space_var_k(busy_preempt,:);
                                    space_srv_k_preempt(:, end-sum(K)+Ks(classpreempt)+phasepreempt) = space_srv_k_preempt(:,end-sum(K)+Ks(classpreempt)+phasepreempt) - 1; % remove preempted job
                                    space_srv_k_preempt(:, end-sum(K)+Ks(class)+kentry) = space_srv_k_preempt(:,end-sum(K)+Ks(class)+kentry) + 1;

                                    % dynamically grow buffer lenght in
                                    % simulation
                                    if isSimulation
                                        if ni < capacity(ist) && nir(class) < classcap(ist,class) % if there is room
                                            if ~any(space_buf_k_preempt(:)==0) % but the buffer has no empty slots
                                                % append job slot
                                                space_buf_k_preempt = [zeros(size(space_buf_k_preempt,1),2),space_buf_k]; % append two columns for (class, preempt-phase)
                                            end
                                        end
                                    end

                                    %get position of first empty slot
                                    empty_slots = -1*ones(sum(busy_preempt),1);
                                    if size(space_buf_k_preempt,2) == 0
                                        empty_slots(busy_preempt) = false;
                                    elseif size(space_buf_k_preempt,2) == 2 % 2 due to (class, preempt-phase) pairs
                                        empty_slots(busy_preempt) = space_buf_k_preempt(busy_preempt,1:2:end)==0;
                                    else
                                        empty_slots(busy_preempt) = max(bsxfun(@times, space_buf_k_preempt(busy_preempt,:)==0, [1:size(space_buf_k_preempt,2)]),[],2)-1; %-1 due to (class, preempt-phase) pairs
                                    end

                                    % ignore states where the buffer has no empty slots
                                    wbuf_empty = empty_slots>0;
                                    if any(wbuf_empty)
                                        space_srv_k_preempt = space_srv_k_preempt(wbuf_empty,:);
                                        space_buf_k_preempt = space_buf_k_preempt(wbuf_empty,:);
                                        space_var_k_preempt = space_var_k_preempt(wbuf_empty,:);
                                        empty_slots = empty_slots(wbuf_empty);
                                        if sn.sched(ist) == SchedStrategy.LCFSPR || sn.sched(ist) == SchedStrategy.LCFSPRPRIO || sn.sched(ist) == SchedStrategy.FCFSPR || sn.sched(ist) == SchedStrategy.FCFSPRPRIO % preempt-resume
                                            space_buf_k_preempt(sub2ind(size(space_buf_k_preempt),1:size(space_buf_k_preempt,1),empty_slots')+1) = phasepreempt;
                                        elseif sn.sched(ist) == SchedStrategy.LCFSPI || sn.sched(ist) == SchedStrategy.LCFSPIPRIO || sn.sched(ist) == SchedStrategy.FCFSPI || sn.sched(ist) == SchedStrategy.FCFSPIPRIO % preempt-independent
                                            space_buf_k_preempt(sub2ind(size(space_buf_k_preempt),1:size(space_buf_k_preempt,1),empty_slots')+1) = 1;
                                        end
                                        space_buf_k_preempt(sub2ind(size(space_buf_k_preempt),1:size(space_buf_k_preempt,1),empty_slots')) = classpreempt;
                                        %outspace(all_busy_srv(wbuf_empty),:) = [space_buf, space_srv, space_var];
                                    end
                                    space_srv_k_reord = [space_srv_k_reord; space_srv_k_preempt];
                                    space_buf_k_reord = [space_buf_k_reord; space_buf_k_preempt];
                                    space_var_k_reord = [space_var_k_reord; space_var_k_preempt];
                                    % the displaced job leaves the server and the
                                    % arriving one takes it, on the same arc
                                    nprm = size(space_srv_k_preempt,1);
                                    preempt_k_reord = [preempt_k_reord; classpreempt*ones(nprm,1)];
                                    start_k_reord = [start_k_reord; class*ones(nprm,1)];
                                end
                            end
                        end

                        % Rows where the arrival can preempt nothing: every busy server holds a
                        % job it may not displace. The job cannot seize a server, so it WAITS in
                        % the buffer. Without this branch the loop above emits no state at all
                        % for those rows -- the arrival transition simply does not exist, so the
                        % class can never enter a busy station and its queue is silently
                        % understated. Python already carries this fallback (its BUG-70 fix);
                        % MATLAB never received it. The job is stored as a (class, entry-phase)
                        % pair exactly as a preempted job is, so on promotion it resumes from
                        % that phase (for exponential service resume == restart).
                        if isPrioAware
                            canPreempt = false(size(space_srv_k,1),1);
                            for cp = 1:R
                                if isLcfsPrioFamily
                                    preemptable = sn.classprio(class) <= sn.classprio(cp);
                                else
                                    preemptable = sn.classprio(class) < sn.classprio(cp);
                                end
                                if preemptable
                                    for ph = 1:K(cp)
                                        canPreempt = canPreempt | (space_srv_k(:, end-sum(K)+Ks(cp)+ph) > 0);
                                    end
                                end
                            end
                            waitRows = all_busy_srv & ~canPreempt;
                            if any(waitRows)
                                space_srv_k_wait = space_srv_k(waitRows,:);
                                space_buf_k_wait = space_buf_k(waitRows,:);
                                space_var_k_wait = space_var_k(waitRows,:);
                                % Grow the buffer in simulation, exactly as the preemption branch
                                % above does. The SSA state vector starts ONE (class,phase) pair
                                % wide and only ever widens at these two sites, so without this a
                                % WAITING arrival was dropped as soon as that pair was taken: the
                                % arrival transition vanished, the station saturated at S+1 jobs
                                % and the queue length fell well below the exact answer. It bit
                                % FCFSPIPRIO hardest, the one preempt-family discipline whose
                                % arrivals both wait and cannot preempt (LCFSPI/FCFSPI preempt
                                % unconditionally and LCFSPIPRIO on equal priority, so all three
                                % reach the preemption grow instead). Measured on a Delay+Queue
                                % closed model, 3+1 jobs: QLen 0.4902 against the exact 0.8028.
                                if isSimulation && ~any(space_buf_k_wait(:)==0) ...
                                        && ni(1) < capacity(ist) && nir(1,class) < classcap(ist,class)
                                    space_buf_k_wait = [zeros(size(space_buf_k_wait,1),2), space_buf_k_wait];
                                end
                                % rightmost empty (class,phase) pair, matching the preemption store
                                empty_wait = -1*ones(sum(waitRows),1);
                                if size(space_buf_k_wait,2) > 0
                                    empty_wait = max(bsxfun(@times, space_buf_k_wait==0, 1:size(space_buf_k_wait,2)),[],2)-1;
                                end
                                keep_wait = empty_wait > 0;
                                if any(keep_wait)
                                    space_srv_k_wait = space_srv_k_wait(keep_wait,:);
                                    space_buf_k_wait = space_buf_k_wait(keep_wait,:);
                                    space_var_k_wait = space_var_k_wait(keep_wait,:);
                                    ew = empty_wait(keep_wait);
                                    nw = size(space_buf_k_wait,1);
                                    space_buf_k_wait(sub2ind(size(space_buf_k_wait),(1:nw)',ew)) = class;
                                    space_buf_k_wait(sub2ind(size(space_buf_k_wait),(1:nw)',ew+1)) = kentry;
                                    psentry = [psentry; ones(nw,1)];
                                    space_srv_k_reord = [space_srv_k_reord; space_srv_k_wait];
                                    space_buf_k_reord = [space_buf_k_reord; space_buf_k_wait];
                                    space_var_k_reord = [space_var_k_reord; space_var_k_wait];
                                    % the arrival preempts nothing and waits: no tag
                                    start_k_reord = [start_k_reord; zeros(nw,1)];
                                    preempt_k_reord = [preempt_k_reord; zeros(nw,1)];
                                end
                            end
                        end
                        outprob_k = pentry(kentry) * psentry .* ones(size(space_srv_k_reord,1),1);
                    end
                    space_buf_k = space_buf_k_reord; % save reordered output states
                    space_srv_k = space_srv_k_reord; % save reordered output states
                    space_var_k = space_var_k_reord; % save reordered output states
                    start_k = start_k_reord; % tags follow the same reordering
                    preempt_k = preempt_k_reord;
            end
            % form the new state
            outspace_k = [space_buf_k, space_srv_k, space_var_k];
            % remove states where new arrival violates capacity or cutoff constraints
            [oi,oir] = State.toMarginalAggr(sn,ind,outspace_k,K,Ks,space_buf_k,space_srv_k,space_var_k);
            en_o = classcap(ist,class)>= oir(:,class) & capacity(ist)*ones(size(oi,1),1) >= oi;

            if size(outspace,2)>size(outspace_k(en_o,:),2)
                outspace = [outspace; zeros(1,size(outspace,2)-size(outspace_k(en_o,:),2)),outspace_k(en_o,:)];
            elseif size(outspace,2)<size(outspace_k(en_o,:),2)
                outspace = [zeros(size(outspace,1),size(outspace_k(en_o,:),2)-size(outspace,2)), outspace; outspace_k(en_o,:)];
            else
                outspace = [outspace; outspace_k(en_o,:)];
            end
            outrate = [outrate; -1*ones(size(outspace_k(en_o,:),1),1)]; % passive action, rate is unspecified
            outprob = [outprob; outprob_k(en_o,:)];
            % tag the arcs just appended; en_o drops the rows the capacity
            % filter deleted, so the tags stay aligned with outspace
            nblk_a = size(outspace_k(en_o,:),1);
            if numel(start_k) ~= size(outspace_k,1)
                line_error(mfilename, sprintf(['Arrival tag vector holds %d rows against %d successor rows at station ''%s'': ' ...
                    'a scheduling branch reselected rows without carrying start_k/preempt_k with them.'], ...
                    numel(start_k), size(outspace_k,1), sn.nodenames{ind}));
            end
            outstart = State.tagArc(outstart, size(outspace,1), nblk_a, R, start_k(en_o));
            outpreempt = State.tagArc(outpreempt, size(outspace,1), nblk_a, R, preempt_k(en_o));
        end
        % Balking (QUEUE_LENGTH strategy): with probability balkProb the
        % arriving class-r job refuses to join, based on the pre-arrival total
        % station population ni. A balked job is lost (destination state left
        % unchanged); the admitted branches are scaled by (1-balkProb). Only
        % the QUEUE_LENGTH strategy is a pure function of the state vector and
        % is therefore admissible in state-space solvers; EXPECTED_WAIT /
        % COMBINED depend on the mean wait and are rejected in the analyzers.
        if isfield(sn,'balkingStrategy') && ~isempty(sn.balkingStrategy) ...
                && sn.balkingStrategy(ist,class) == BalkingStrategy.QUEUE_LENGTH ...
                && ~isempty(outspace)
            balkProb = 0;
            thresholds = sn.balkingThresholds{ist,class};
            qlen = ni(1); % pre-arrival total station population
            for ti = 1:numel(thresholds)
                th = thresholds{ti};
                if qlen >= th{1} && qlen <= th{2}
                    balkProb = th{3};
                    break
                end
            end
            if balkProb > 0
                outprob = outprob * (1 - balkProb);
                % balked branch: job lost, destination state unchanged
                inrow = inspace;
                if size(outspace,2) > size(inrow,2)
                    inrow = [zeros(1,size(outspace,2)-size(inrow,2)), inrow];
                end
                outspace = [outspace; inrow];
                outrate  = [outrate; -1];
                outprob  = [outprob; balkProb];
            end
        end
        if isSimulation
            if size(outprob,1) > 1
                cum_prob = cumsum(outprob) / sum(outprob);
                firing_ctr = 1 + max([0,find( rand > cum_prob' )]); % select action
                outstart = State.tagPad(outstart, size(outspace,1), R);
                outpreempt = State.tagPad(outpreempt, size(outspace,1), R);
                outspace = outspace(firing_ctr,:);
                outstart = outstart(firing_ctr,:); % the tags of the sampled arc
                outpreempt = outpreempt(firing_ctr,:);
                outrate = -1;
                outprob = 1;
            end
        end
    case EventType.DEP
        % Marked (MMAP) source class: the shared modulating chain lives in the
        % carrier's phase block (mark index 1); this class's departures fire
        % from there using its per-mark matrix D1k (M3A cell index 2+mark).
        markofclass = -1;
        phclass = class;
        if sn.sched(ist) == SchedStrategy.EXT && isfield(sn,'markidx') ...
                && ~isempty(sn.markidx) && ist <= size(sn.markidx,1) ...
                && sn.markidx(ist,class) > 0
            markofclass = sn.markidx(ist,class);
            phclass = find(sn.markidx(ist,:) == 1, 1);
        end
        if any(any(space_srv(:,(Ks(phclass)+1):(Ks(phclass)+K(phclass))))) % something is busy
            if hasOnlyExp && (sn.sched(ist) == SchedStrategy.PS || sn.sched(ist) == SchedStrategy.DPS || sn.sched(ist) == SchedStrategy.GPS || sn.sched(ist) == SchedStrategy.INF || sn.sched(ist) == SchedStrategy.PSPRIO || sn.sched(ist) == SchedStrategy.DPSPRIO || sn.sched(ist) == SchedStrategy.GPSPRIO || sn.sched(ist) == SchedStrategy.LPS)
                nir = space_srv;
                ni = sum(nir,2);
                sir = nir;
                kir = sir;
            else
                [ni,nir,sir,kir] = State.toMarginal(sn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
            end
            switch sn.routing(ind,class)
                case RoutingStrategy.RROBIN
                    idx = find(space_var(sum(sn.nvars(ind,1:(R+class)))) == sn.nodeparam{ind}{class}.outlinks);
                    if idx < length(sn.nodeparam{ind}{class}.outlinks)
                        space_var(sum(sn.nvars(ind,1:(R+class)))) = sn.nodeparam{ind}{class}.outlinks(idx+1);
                    else
                        space_var(sum(sn.nvars(ind,1:(R+class)))) = sn.nodeparam{ind}{class}.outlinks(1);
                    end
                case RoutingStrategy.WRROBIN
                    % WRR slot holds a POSITION in weighted_outlinks; advance it
                    % cyclically (mirrors afterEventRouter and sub_wrr). Without
                    % this advance the position stays fixed and only the initial
                    % destination is ever selected.
                    slot = sum(sn.nvars(ind,1:(R+class)));
                    if isfield(sn.nodeparam{ind}{class}, 'weighted_outlinks') ...
                            && ~isempty(sn.nodeparam{ind}{class}.weighted_outlinks)
                        cycle_len = length(sn.nodeparam{ind}{class}.weighted_outlinks);
                    else
                        cycle_len = length(sn.nodeparam{ind}{class}.outlinks);
                    end
                    pos = space_var(slot);
                    if pos < 1 || pos >= cycle_len
                        space_var(slot) = 1;
                    else
                        space_var(slot) = pos + 1;
                    end
            end
            if sir(phclass)>0 % is a job of class is in service
                outprob = [];
                for k=1:K(phclass)
                    space_srv = inspace(:,(end-sum(K)-V+1):(end-V)); % server state
                    space_buf = inspace(:,1:(end-sum(K)-V)); % buffer state
                    rate = zeros(size(space_srv,1),1);
                    en =  space_srv(:,Ks(phclass)+k) > 0;
                    if any(en)
                        switch sn.sched(ist)
                            case SchedStrategy.EXT % source, can produce an arrival from phase-k as long as it is from an open class
                                if isinf(sn.njobs(class))
                                    if markofclass > 0
                                        % per-mark arrival matrix over the shared chain
                                        D1_srv = proc{ist}{class}{2+markofclass};
                                    else
                                        D1_srv = proc{ist}{class}{2};
                                    end
                                    for kentry = 1:K(phclass)
                                        arv_rate = D1_srv(k, kentry);
                                        if arv_rate <= 0
                                            continue
                                        end
                                        space_srv = inspace(:,(end-sum(K)-V+1):(end-V));
                                        space_srv(en,Ks(phclass)+k) = space_srv(en,Ks(phclass)+k) - 1;
                                        space_srv(en,Ks(phclass)+kentry) = space_srv(en,Ks(phclass)+kentry) + 1;
                                        outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                        if isinf(ni)
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*arv_rate*ones(size(inspace(en,:),1),1)];
                                        else
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*arv_rate*ones(size(inspace(en,:),1),1)];
                                        end
                                        outprob = [outprob; ones(size(space_buf(en,:),1),1)];
                                    end
                                end
                            case SchedStrategy.INF % move first job in service
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(en,class,k); % assume active
                                % if state is unchanged, still add with rate 0
                                outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                else
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                end
                                outprob = [outprob; ones(size(rate(en,:),1),1)];
                            case {SchedStrategy.PS, SchedStrategy.LPS} % move first job in service
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*(kir(en,class,k)./ni(en)).*min(ni(en),S(ist)); % assume active
                                % if state is unchanged, still add with rate 0
                                outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                else
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                end
                                outprob = [outprob; ones(size(rate(en,:),1),1)];
                                %                                    end
                            case SchedStrategy.PSPRIO
                                % unclear if LD scaling should be with
                                % ni or with niprio, for now left as ni
                                % for consistency with HOL multiserver
                                if all(ni(en) <= S(ist))
                                    % n <= c: all jobs get service, priority doesn't matter
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*(kir(en,class,k)./ni(en)).*min(ni(en),S(ist)); % assume active
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                elseif sn.classprio(class) == min(sn.classprio(nir>0)) % if this class is in the most urgent priority group (lower value = higher priority in LINE)
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    niprio(en) = sum(nir(sn.classprio==sn.classprio(class))); % jobs at the same priority class
                                    rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*(kir(en,class,k)./niprio(en)).*min(niprio(en),S(ist)); % assume active
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(niprio(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                else % n > c and not highest priority: set rate to zero
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    outrate = [outrate; zeros(size(rate(en,:),1),1)];
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                end
                            case SchedStrategy.DPS
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                if S(ist) > 1
                                    line_error(mfilename,'Multi-server DPS stations are not supported yet.');
                                end
                                % in GPS, the scheduling parameter are the weights
                                w_i = sn.schedparam(ist,:);
                                w_i = w_i / sum(w_i);
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k))*(kir(en,class,k)/nir(class))*w_i(class)*nir(class)./(sum(repmat(w_i,sum(en),1)*nir',2));
                                % if state is unchanged, still add with rate 0
                                outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                else
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                end
                                outprob = [outprob; ones(size(rate(en,:),1),1)];
                            case SchedStrategy.DPSPRIO
                                if all(ni(en) <= S(ist))
                                    % n <= c: all jobs get service, behave like regular DPS
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server DPS stations are not supported yet.');
                                    end
                                    w_i = sn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k))*(kir(en,class,k)/nir(class))*w_i(class)*nir(class)./(sum(repmat(w_i,sum(en),1)*nir',2));
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                elseif sn.classprio(class) == min(sn.classprio(nir>0)) % if this class is in the most urgent priority group (lower value = higher priority in LINE)
                                    nirprio = nir;
                                    nirprio(sn.classprio~=sn.classprio(class)) = 0; % ignore jobs of lower priority
                                    niprio = sum(nirprio);
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server DPS stations are not supported yet.');
                                    end
                                    % in GPS, the scheduling parameter are the weights
                                    w_i = sn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k))*(kir(en,class,k)/nirprio(class))*w_i(class)*nirprio(class)./(sum(repmat(w_i,sum(en),1)*nirprio',2));
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nirprio,en,class).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nirprio,en,class).*lldscaling(ist,min(niprio(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                else % n > c and not most urgent priority: set rate to zero
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    outrate = [outrate; zeros(size(rate(en,:),1),1)];
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                end
                            case SchedStrategy.GPS
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                if S(ist) > 1
                                    line_error(mfilename,'Multi-server GPS stations are not supported yet.');
                                end
                                % in GPS, the scheduling parameter are the weights
                                w_i = sn.schedparam(ist,:);
                                w_i = w_i / sum(w_i);
                                cir = min(nir,ones(size(nir)));
                                rate = mu{ist}{class}(k)*(phi{ist}{class}(k))*(kir(en,class,k)/nir(class))*w_i(class)/(w_i*cir(:)); % assume active
                                % if state is unchanged, still add with rate 0
                                outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                else
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                end
                                outprob = [outprob; ones(size(rate(en,:),1),1)];
                            case SchedStrategy.GPSPRIO
                                if all(ni(en) <= S(ist))
                                    % n <= c: all jobs get service, behave like regular GPS
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server GPS stations are not supported yet.');
                                    end
                                    w_i = sn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    cir = min(nir,ones(size(nir)));
                                    rate = mu{ist}{class}(k)*(phi{ist}{class}(k))*(kir(en,class,k)/nir(class))*w_i(class)/(w_i*cir(:)); % assume active
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                elseif sn.classprio(class) == min(sn.classprio(nir>0)) % if this class is in the most urgent priority group (lower value = higher priority in LINE)
                                    nirprio = nir;
                                    nirprio(sn.classprio~=sn.classprio(class)) = 0; % ignore jobs of lower priority
                                    niprio = sum(nirprio);
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server DPS stations are not supported yet.');
                                    end
                                    % in GPS, the scheduling parameter are the weights
                                    w_i = sn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    cir = min(nirprio,ones(size(nirprio)));
                                    rate = mu{ist}{class}(k)*(phi{ist}{class}(k))*(kir(en,class,k)/nirprio(class))*w_i(class)/(w_i*cir(:)); % assume active
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nirprio,en,class).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nirprio,en,class).*lldscaling(ist,min(niprio(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                else % n > c and not most urgent priority: set rate to zero
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    outrate = [outrate; zeros(size(rate(en,:),1),1)];
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                end
                            case SchedStrategy.FCFS % move first job in service
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                % Servers held for a pending REPLY are unavailable, so a
                                % job waits in the buffer already when ni exceeds the
                                % REMAINING servers. Zero for models without replies.
                                [~, nbD] = State.replyBlocked(sn, ind, space_var);
                                en_wbuf = en & ni>(S(ist)-nbD); %states with jobs in buffer
                                if noPromote || isRetrialStation % immediate feedback / retrial orbit: hold server, do not promote a waiting job
                                    en_wbuf(:) = false;
                                end
                                % Synchronous call: this departing job keeps its server
                                % until its REPLY signal returns here, so the server is
                                % NOT handed to a waiting job; it is recorded as held in
                                % the reply block instead. Mirrors LDES, which omits the
                                % markServerIdle call and records a pendingReply.
                                if isfield(sn,'replyblock') && ~isempty(sn.replyblock) ...
                                        && size(sn.replyblock,1) >= ind && sn.replyblock(ind,class) > 0
                                    rinfoD = State.replyBlockInfo(sn, ind);
                                    en_wbuf(:) = false;
                                    space_var(en, rinfoD.slot(class)) = space_var(en, rinfoD.slot(class)) + 1;
                                end
                                for kdest=1:K(class) % new phase
                                    space_buf_kd = space_buf;
                                    space_var_kd = space_var;
                                    if ismkvmodclass(class)
                                        space_var_kd(en,sum(sn.nvars(ind,1:class))) = kdest;
                                    end
                                    rate_kd = rate;
                                    rate_kd(en) = proc{ist}{class}{2}(k,kdest).*kir(en,class,k); % assume active
                                    % first process states without jobs in buffer
                                    en_wobuf = ~en_wbuf;
                                    if any(en_wobuf) %any state without jobs in buffer
                                        outspace = [outspace; space_buf_kd(en_wobuf,:), space_srv(en_wobuf,:), space_var_kd(en_wobuf,:)];
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,end).*rate_kd(en_wobuf,:)];
                                        else
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate_kd(en_wobuf,:)];
                                        end
                                    end
                                    % now process states with jobs in buffer
                                    outprob = [outprob; ones(size(rate_kd(en_wobuf,:),1),1)];
                                    if any(en_wbuf) %any state with jobs in buffer
                                        % get class of job at head
                                        start_svc_class = space_buf_kd(en_wbuf,end);
                                        if start_svc_class > 0 % redunant if?
                                            % update input buffer
                                            space_buf_kd(en_wbuf,:) = [zeros(sum(en_wbuf),1),space_buf_kd(en_wbuf,1:end-1)];
                                            % probability vector for the next job of starting in phase kentry
                                            if ismkvmodclass(start_svc_class) % if markov-modulated
                                                if start_svc_class==class % if successive service from the same class
                                                    kentry_range = kdest; % new job enters in phase left by departing job
                                                else % resume phase from local variables
                                                    kentry_range = space_var_kd(en,sum(sn.nvars(ind,1:start_svc_class)));
                                                end
                                                pentry_svc_class = 0*pie{ist}{start_svc_class};
                                                pentry_svc_class(kentry_range) = 1.0;
                                            else % if i.i.d.
                                                pentry_svc_class = pie{ist}{start_svc_class};
                                                kentry_range = 1:K(start_svc_class);
                                            end
                                            for kentry = kentry_range
                                                space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                                outspace = [outspace; space_buf_kd(en,:), space_srv(en,:), space_var_kd(en,:)];
                                                % the head of the buffer takes the vacated server
                                                cls_d = zeros(sum(en),1);
                                                cls_d(en_wbuf(en)) = start_svc_class;
                                                outstart = State.tagArc(outstart, size(outspace,1), sum(en), R, cls_d);
                                                rate_k = rate_kd;
                                                rate_k(en_wbuf,:) = rate_kd(en_wbuf,:)*pentry_svc_class(kentry);
                                                if isinf(ni) % use limited load-dependence at the latest user-provided level
                                                    newrate = State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate_k(en,:);
                                                else
                                                    newrate = State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en,:);
                                                end
                                                outrate = [outrate; newrate];
                                                % Zero the probability of the branches added *by this kentry* when
                                                % their rate is zero. The mask must index newrate, not the whole
                                                % accumulated outrate: outrate grows by one entry per kentry while
                                                % outprob_cur has one entry per enabled state, so masking on outrate
                                                % silently GREW outprob_cur (MATLAB expands on out-of-bound
                                                % assignment) once any earlier branch had rate 0 -- which happens
                                                % whenever pentry_svc_class has zeros, i.e. a PH whose entry vector
                                                % pie does not reach every phase. That misaligned outprob against
                                                % outspace/outrate, so the sampled branch read a bogus (often 0)
                                                % probability and depRatesSamples under-counted departures.
                                                outprob_cur = ones(size(rate_kd(en,:),1),1);
                                                outprob_cur(newrate==0.0) = 0;
                                                outprob = [outprob; outprob_cur(:)];
                                                space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) - 1;
                                            end
                                        end
                                    end
                                end
                                % if state is unchanged, still add with rate 0
                            case SchedStrategy.HOL % FCFS priority
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                if noPromote || isRetrialStation % immediate feedback / retrial orbit: hold server, do not promote a waiting job
                                    en_wbuf(:) = false;
                                end
                                en_wobuf = ~en_wbuf;
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                priogroup = [Inf,sn.classprio]; % Inf for empty positions (lower value = higher priority)
                                space_buf_groupg = arrayfun(@(x) priogroup(1+x), space_buf);
                                start_classprio = min(space_buf_groupg(en_wbuf,:),[],2); % min finds highest priority
                                isrowmax = space_buf_groupg == repmat(start_classprio, 1, size(space_buf_groupg,2));
                                [~,rightmostMaxPosFlipped]=max(fliplr(isrowmax),[],2);
                                rightmostMaxPos = size(isrowmax,2) - rightmostMaxPosFlipped + 1;
                                start_svc_class = space_buf(en_wbuf, rightmostMaxPos);
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];
                                if start_svc_class > 0
                                    pentry_svc_class = pie{ist}{start_svc_class};
                                    for kentry = 1:K(start_svc_class)
                                        space_srv_k = space_srv;
                                        space_buf_k = space_buf;
                                        space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                        for j=find(en_wbuf)'
                                            space_buf_k(j,:) = [0, space_buf_k(j,1:rightmostMaxPos(j)-1), space_buf_k(j,(rightmostMaxPos(j)+1):end)];
                                        end
                                        % if state is unchanged, still add with rate 0
                                        outspace = [outspace; space_buf_k(en_wbuf,:), space_srv_k(en_wbuf,:), space_var(en_wbuf,:)];
                                        outstart = State.tagArc(outstart, size(outspace,1), sum(en_wbuf), R, start_svc_class(:));
                                        rate_k = rate;
                                        rate_k(en_wbuf,:) = rate(en_wbuf,:) * pentry_svc_class(kentry);
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,end).*rate_k(en_wbuf,:)];
                                        else
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en_wbuf,:)];
                                        end
                                        outprob = [outprob; ones(size(rate_k(en_wbuf,:),1),1)];
                                    end
                                end
                            case SchedStrategy.LCFSPRIO % LCFS priority - like HOL but LCFS order within priority groups
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                if noPromote || isRetrialStation % immediate feedback / retrial orbit: hold server, do not promote a waiting job
                                    en_wbuf(:) = false;
                                end
                                en_wobuf = ~en_wbuf;
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                priogroup = [Inf,sn.classprio]; % Inf for empty positions (lower value = higher priority)
                                space_buf_groupg = arrayfun(@(x) priogroup(1+x), space_buf);
                                start_classprio = min(space_buf_groupg(en_wbuf,:),[],2); % min finds highest priority
                                isrowmax = space_buf_groupg == repmat(start_classprio, 1, size(space_buf_groupg,2));
                                % LCFS: Find leftmost (first) position instead of rightmost for LCFS order
                                [~,leftmostMaxPos]=max(isrowmax,[],2);
                                start_svc_class = space_buf(en_wbuf, leftmostMaxPos);
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];
                                if start_svc_class > 0
                                    pentry_svc_class = pie{ist}{start_svc_class};
                                    for kentry = 1:K(start_svc_class)
                                        space_srv_k = space_srv;
                                        space_buf_k = space_buf;
                                        space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                        for j=find(en_wbuf)'
                                            % LCFS: Remove from leftmost position instead of rightmost
                                            space_buf_k(j,:) = [0, space_buf_k(j,1:leftmostMaxPos(j)-1), space_buf_k(j,(leftmostMaxPos(j)+1):end)];
                                        end
                                        % if state is unchanged, still add with rate 0
                                        outspace = [outspace; space_buf_k(en_wbuf,:), space_srv_k(en_wbuf,:), space_var(en_wbuf,:)];
                                        outstart = State.tagArc(outstart, size(outspace,1), sum(en_wbuf), R, start_svc_class(:));
                                        rate_k = rate;
                                        rate_k(en_wbuf,:) = rate_k(en_wbuf,:)*pentry_svc_class(kentry);
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,end).*rate_k(en_wbuf,:)];
                                        else
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en_wbuf,:)];
                                        end
                                        outprob = [outprob; ones(size(rate_k(en_wbuf,:),1),1)];
                                        space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) - 1;
                                    end
                                end
                            case SchedStrategy.LCFS % move last job in service
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                if noPromote || isRetrialStation % immediate feedback / retrial orbit: hold server, do not promote a waiting job
                                    en_wbuf(:) = false;
                                end
                                % Plain LCFS is not priority-aware: it always promotes the most
                                % recent arrival. Priorities are honored by SchedStrategy.LCFSPRIO,
                                % which has its own case below ("LCFS order within priority
                                % groups"), exactly as FCFS relates to FCFSPRIO. Branching here on
                                % ~all(classprio == classprio(1)) silently turned every LCFS
                                % station with distinct class priorities into an LCFSPRIO one,
                                % contradicting the "Priorities will be ignored" warning raised by
                                % refreshStruct and, for three or more classes, yielding an
                                % all-zero queue length.
                                [~, colfirstnnz] = max( space_buf(en_wbuf,:) ~=0, [], 2 ); % find first nnz column
                                start_svc_class = space_buf(en_wbuf,colfirstnnz); % job entering service
                                space_buf(en_wbuf,colfirstnnz)=0;
                                if isempty(start_svc_class)
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    if isSimulation && nargin>=7 && isobject(eventCache)
                                        eventCache{key} = {outprob, outspace, outrate, outstart, outpreempt};
                                    end
                                    return
                                end
                                for kentry = 1:K(start_svc_class)
                                    pentry_svc_class = pie{ist}{start_svc_class};
                                    space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    cls_d = zeros(sum(en),1);
                                    cls_d(en_wbuf(en)) = start_svc_class;
                                    outstart = State.tagArc(outstart, size(outspace,1), sum(en), R, cls_d);
                                    rate_k = rate;
                                    rate_k(en_wbuf,:) = rate(en_wbuf,:)*pentry_svc_class(kentry);
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate_k(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) - 1;
                                end
                            case SchedStrategy.LCFSPR % move last job in service (preempt-resume)
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                % Plain LCFSPR resumes the most recently preempted job. It is not
                                % priority-aware: LCFSPRPRIO is, and it carries its own handling.
                                % Selecting by priority whenever the class priorities happened to
                                % differ made this case answer for a policy the user never declared
                                % (see the arrival path above for the same rule).
                                [~, colfirstnnz] = max( space_buf(en_wbuf,:) ~=0, [], 2 ); % find first nnz column
                                start_svc_class = space_buf(en_wbuf,colfirstnnz); % job entering service
                                kentry = space_buf(en_wbuf,colfirstnnz+1); % entry phase of job resuming service
                                space_buf(en_wbuf,colfirstnnz)=0;% zero popped job
                                space_buf(en_wbuf,colfirstnnz+1)=0; % zero popped phase
                                if isempty(start_svc_class)
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    if isSimulation && nargin>=7 && isobject(eventCache)
                                        eventCache{key} = {outprob, outspace, outrate, outstart, outpreempt};
                                    end
                                    return
                                end
                                space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                % if state is unchanged, still add with rate 0
                                outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                % the most recently preempted job resumes on the freed server
                                cls_d = zeros(sum(en),1);
                                cls_d(en_wbuf(en)) = start_svc_class;
                                outstart = State.tagArc(outstart, size(outspace,1), sum(en), R, cls_d);
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                else
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                end
                                outprob = [outprob; ones(size(rate(en,:),1),1)];
                            case SchedStrategy.LCFSPI % move last job in service (preempt-independent)
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                [~, colfirstnnz] = max( space_buf(en_wbuf,:) ~=0, [], 2 ); % find first nnz column
                                start_svc_class = space_buf(en_wbuf,colfirstnnz); % job entering service
                                space_buf(en_wbuf,colfirstnnz)=0;% zero popped job
                                space_buf(en_wbuf,colfirstnnz+1)=0; % zero popped phase (ignored for LCFSPI)
                                if isempty(start_svc_class)
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    if isSimulation && nargin>=7 && isobject(eventCache)
                                        eventCache{key} = {outprob, outspace, outrate, outstart, outpreempt};
                                    end
                                    return
                                end
                                % For LCFSPI, jobs restart from pie distribution instead of stored phase
                                pentry_svc_class = pie{ist}{start_svc_class};
                                for kentry = 1:K(start_svc_class)
                                    space_srv_k = space_srv;
                                    space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv_k(en,:), space_var(en,:)];
                                    cls_d = zeros(sum(en),1);
                                    cls_d(en_wbuf(en)) = start_svc_class;
                                    outstart = State.tagArc(outstart, size(outspace,1), sum(en), R, cls_d);
                                    rate_k = rate;
                                    rate_k(en_wbuf,:) = rate_k(en_wbuf,:) * pentry_svc_class(kentry);
                                    if isinf(ni) % hit limited load-dependence
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate_k(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate_k(en,:),1),1)];
                                end
                            case SchedStrategy.FCFSPR % FCFS preempt-resume (no priority)
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                % FCFS: Find rightmost (last) non-zero position (longest waiting job)
                                [~, colLastNnz] = max(fliplr(space_buf(en_wbuf,:) ~= 0), [], 2);
                                colLastNnz = size(space_buf,2) - colLastNnz; % convert to actual position (odd index for class)
                                start_svc_class = space_buf(en_wbuf, colLastNnz); % job entering service
                                kentry = space_buf(en_wbuf, colLastNnz+1); % entry phase of job resuming service
                                % The (class,phase) buffer is RIGHT-aligned, so removing the oldest
                                % pair -- which sits at the right -- must pad a whole empty PAIR on
                                % the left rather than leave a hole in place. A hole is a layout
                                % fromMarginal never enumerates, so the successor was unreachable and
                                % every full-buffer state came out absorbing (ctmc_solve then reported
                                % "no recurrent state"). Same rule as the LCFSPRPRIO arm below.
                                jbuf = 0;
                                for j = find(en_wbuf)'
                                    jbuf = jbuf + 1;
                                    cpop = colLastNnz(jbuf);
                                    space_buf(j,:) = [0, 0, space_buf(j,1:cpop-1), space_buf(j,(cpop+2):end)];
                                end
                                if isempty(start_svc_class)
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni)
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    if isSimulation && nargin>=7 && isobject(eventCache)
                                        eventCache{key} = {outprob, outspace, outrate, outstart, outpreempt};
                                    end
                                    return
                                end
                                space_srv(en_wbuf, Ks(start_svc_class)+kentry) = space_srv(en_wbuf, Ks(start_svc_class)+kentry) + 1;
                                outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                % the longest-waiting preempted job resumes on the freed server
                                cls_d = zeros(sum(en),1);
                                cls_d(en_wbuf(en)) = start_svc_class;
                                outstart = State.tagArc(outstart, size(outspace,1), sum(en), R, cls_d);
                                if isinf(ni)
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                else
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                end
                                outprob = [outprob; ones(size(rate(en,:),1),1)];
                            case SchedStrategy.FCFSPI % FCFS preempt-independent (no priority)
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                % FCFS: Find rightmost (last) non-zero position (longest waiting job)
                                [~, colLastNnz] = max(fliplr(space_buf(en_wbuf,:) ~= 0), [], 2);
                                colLastNnz = size(space_buf,2) - colLastNnz; % convert to actual position (odd index for class)
                                start_svc_class = space_buf(en_wbuf, colLastNnz); % job entering service
                                % The (class,phase) buffer is RIGHT-aligned, so removing the oldest
                                % pair -- which sits at the right -- must pad a whole empty PAIR on
                                % the left rather than leave a hole in place. A hole is a layout
                                % fromMarginal never enumerates, so the successor was unreachable and
                                % every full-buffer state came out absorbing (ctmc_solve then reported
                                % "no recurrent state"). Same rule as the LCFSPRPRIO arm below.
                                jbuf = 0;
                                for j = find(en_wbuf)'
                                    jbuf = jbuf + 1;
                                    cpop = colLastNnz(jbuf);
                                    space_buf(j,:) = [0, 0, space_buf(j,1:cpop-1), space_buf(j,(cpop+2):end)];
                                end
                                if isempty(start_svc_class)
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    if isinf(ni)
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    if isSimulation && nargin>=7 && isobject(eventCache)
                                        eventCache{key} = {outprob, outspace, outrate, outstart, outpreempt};
                                    end
                                    return
                                end
                                % For FCFSPI, jobs restart from pie distribution instead of stored phase
                                pentry_svc_class = pie{ist}{start_svc_class};
                                for kentry = 1:K(start_svc_class)
                                    space_srv_k = space_srv;
                                    space_srv_k(en_wbuf, Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf, Ks(start_svc_class)+kentry) + 1;
                                    outspace = [outspace; space_buf(en,:), space_srv_k(en,:), space_var(en,:)];
                                    cls_d = zeros(sum(en),1);
                                    cls_d(en_wbuf(en)) = start_svc_class;
                                    outstart = State.tagArc(outstart, size(outspace,1), sum(en), R, cls_d);
                                    rate_k = rate;
                                    rate_k(en_wbuf,:) = rate_k(en_wbuf,:) * pentry_svc_class(kentry);
                                    if isinf(ni)
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate_k(en,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en,:)];
                                    end
                                    outprob = [outprob; ones(size(rate_k(en,:),1),1)];
                                end
                            case SchedStrategy.FCFSPRPRIO % FCFS preempt-resume with priority groups
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                en_wobuf = ~en_wbuf;
                                priogroup = [Inf,sn.classprio]; % Inf for empty positions (lower value = higher priority)
                                % Buffer stores [class,phase,class,phase,...] pairs;
                                % only inspect class columns (odd: 1,3,5,...) for priority
                                class_cols = 1:2:size(space_buf,2);
                                space_buf_class = space_buf(:, class_cols);
                                space_buf_class_groupg = arrayfun(@(x) priogroup(1+x), space_buf_class);
                                start_classprio = min(space_buf_class_groupg(en_wbuf,:),[],2); % min finds highest priority
                                isrowmax = space_buf_class_groupg == repmat(start_classprio, 1, size(space_buf_class_groupg,2));
                                % FCFS: Find rightmost (last) class position for highest priority class (longest waiting)
                                [~,rightmostClassPos]=max(fliplr(isrowmax),[],2);
                                rightmostClassPos = size(space_buf_class_groupg,2) - rightmostClassPos + 1;
                                rightmostMaxPos = 2*rightmostClassPos - 1; % convert to buffer column index
                                start_svc_class = space_buf(en_wbuf, rightmostMaxPos); % job entering service
                                kentry = space_buf(en_wbuf, rightmostMaxPos+1); % entry phase of job resuming service (preempt-resume)

                                % Handle states without buffer jobs
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni)
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];

                                % Handle states with buffer jobs
                                if any(en_wbuf) && start_svc_class > 0
                                    space_srv_k = space_srv;
                                    space_buf_k = space_buf;
                                    space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                    for j=find(en_wbuf)'
                                        % Remove both class and phase from rightmost position (preempt-resume)
                                        space_buf_k(j,:) = [0, 0, space_buf_k(j,1:rightmostMaxPos(j)-1), space_buf_k(j,(rightmostMaxPos(j)+2):end)];
                                    end
                                    outspace = [outspace; space_buf_k(en_wbuf,:), space_srv_k(en_wbuf,:), space_var(en_wbuf,:)];
                                    outstart = State.tagArc(outstart, size(outspace,1), sum(en_wbuf), R, start_svc_class(:));
                                    if isinf(ni)
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,end).*rate(en_wbuf,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wbuf,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en_wbuf,:),1),1)];
                                end
                            case SchedStrategy.LCFSPRPRIO % LCFS preempt-resume with priority groups
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                en_wobuf = ~en_wbuf;
                                priogroup = [Inf,sn.classprio]; % Inf for empty positions (lower value = higher priority)
                                % Buffer stores [class,phase,class,phase,...] pairs;
                                % only inspect class columns (odd: 1,3,5,...) for priority
                                class_cols = 1:2:size(space_buf,2);
                                space_buf_class = space_buf(:, class_cols);
                                space_buf_class_groupg = arrayfun(@(x) priogroup(1+x), space_buf_class);
                                start_classprio = min(space_buf_class_groupg(en_wbuf,:),[],2); % min finds highest priority
                                isrowmax = space_buf_class_groupg == repmat(start_classprio, 1, size(space_buf_class_groupg,2));
                                % LCFS: Find leftmost (first) class position for highest priority class (most recently preempted)
                                [~,leftmostClassPos]=max(isrowmax,[],2);
                                leftmostMaxPos = 2*leftmostClassPos - 1; % convert to buffer column index
                                start_svc_class = space_buf(en_wbuf, leftmostMaxPos); % job entering service
                                kentry = space_buf(en_wbuf, leftmostMaxPos+1); % entry phase of job resuming service (preempt-resume)
                                
                                % Handle states without buffer jobs
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni)
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];
                                
                                % Handle states with buffer jobs
                                if any(en_wbuf) && start_svc_class > 0
                                    space_srv_k = space_srv;
                                    space_buf_k = space_buf;
                                    space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                    for j=find(en_wbuf)'
                                        % Remove both class and phase from leftmost position (preempt-resume)
                                        % The (class,phase) pair buffer is RIGHT-aligned: occupied pairs sit at the
                                        % right, empty pairs pad the left. Priority can promote a pair from the
                                        % middle, so closing the hole must pad a whole empty PAIR on the left --
                                        % as State.afterEventStationSignal's dropWaiting does. Padding one slot
                                        % left and one right instead left a trailing empty slot, e.g.
                                        % [2 1 3 1] -> [0 3 1 0], a layout fromMarginal never enumerates: the
                                        % successor was unreachable, which is what made the generator reducible.
                                        space_buf_k(j,:) = [0, 0, space_buf_k(j,1:leftmostMaxPos(j)-1), space_buf_k(j,(leftmostMaxPos(j)+2):end)];
                                    end
                                    outspace = [outspace; space_buf_k(en_wbuf,:), space_srv_k(en_wbuf,:), space_var(en_wbuf,:)];
                                    outstart = State.tagArc(outstart, size(outspace,1), sum(en_wbuf), R, start_svc_class(:));
                                    if isinf(ni)
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,end).*rate(en_wbuf,:)];
                                    else
                                        outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wbuf,:)];
                                    end
                                    outprob = [outprob; ones(size(rate(en_wbuf,:),1),1)];
                                end
                            case SchedStrategy.FCFSPIPRIO % FCFS preempt-independent with priority groups
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                en_wobuf = ~en_wbuf;
                                priogroup = [Inf,sn.classprio]; % Inf for empty positions (lower value = higher priority)
                                % Buffer stores [class,phase,class,phase,...] pairs;
                                % only inspect class columns (odd: 1,3,5,...) for priority.
                                % Scanning every column instead mapped a PHASE index p
                                % through priogroup(1+p), so a phase masqueraded as class p
                                % and could win the min: the promoted job was then read off
                                % a phase column, i.e. a class with nothing waiting, and the
                                % pair removal straddled two entries. The PR-PRIO arms above
                                % already restrict to class_cols for this reason.
                                class_cols = 1:2:size(space_buf,2);
                                space_buf_class = space_buf(:, class_cols);
                                space_buf_class_groupg = arrayfun(@(x) priogroup(1+x), space_buf_class);
                                start_classprio = min(space_buf_class_groupg(en_wbuf,:),[],2); % min finds highest priority
                                isrowmax = space_buf_class_groupg == repmat(start_classprio, 1, size(space_buf_class_groupg,2));
                                % FCFS: Find rightmost (last) class position for highest priority class
                                [~,rightmostClassPos]=max(fliplr(isrowmax),[],2);
                                rightmostClassPos = size(space_buf_class_groupg,2) - rightmostClassPos + 1;
                                rightmostMaxPos = 2*rightmostClassPos - 1; % convert to buffer column index
                                start_svc_class = space_buf(en_wbuf, rightmostMaxPos); % job entering service

                                % Handle states without buffer jobs
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni)
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];

                                % Handle states with buffer jobs
                                if any(en_wbuf) && start_svc_class > 0
                                    % For FCFSPIPRIO, jobs restart from pie distribution instead of stored phase
                                    pentry_svc_class = pie{ist}{start_svc_class};
                                    for kentry = 1:K(start_svc_class)
                                        space_srv_k = space_srv;
                                        space_buf_k = space_buf;
                                        space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                        for j=find(en_wbuf)'
                                            % Remove both class and phase from rightmost position (preempt-independent ignores stored phase)
                                            space_buf_k(j,:) = [0, 0, space_buf_k(j,1:rightmostMaxPos(j)-1), space_buf_k(j,(rightmostMaxPos(j)+2):end)];
                                        end
                                        outspace = [outspace; space_buf_k(en_wbuf,:), space_srv_k(en_wbuf,:), space_var(en_wbuf,:)];
                                        outstart = State.tagArc(outstart, size(outspace,1), sum(en_wbuf), R, start_svc_class(:));
                                        rate_k = rate;
                                        rate_k(en_wbuf,:) = rate_k(en_wbuf,:) * pentry_svc_class(kentry);
                                        if isinf(ni)
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,end).*rate_k(en_wbuf,:)];
                                        else
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en_wbuf,:)];
                                        end
                                        outprob = [outprob; ones(size(rate_k(en_wbuf,:),1),1)];
                                    end
                                end
                            case SchedStrategy.LCFSPIPRIO % LCFS preempt-independent with priority groups
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % assume active
                                en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                en_wobuf = ~en_wbuf;
                                priogroup = [Inf,sn.classprio]; % Inf for empty positions (lower value = higher priority)
                                % Buffer stores [class,phase,class,phase,...] pairs;
                                % only inspect class columns (odd: 1,3,5,...) for priority.
                                % Scanning every column instead mapped a PHASE index p
                                % through priogroup(1+p), so a phase masqueraded as class p
                                % and could win the min: the promoted job was then read off
                                % a phase column, i.e. a class with nothing waiting, and the
                                % pair removal straddled two entries. The PR-PRIO arms above
                                % already restrict to class_cols for this reason.
                                class_cols = 1:2:size(space_buf,2);
                                space_buf_class = space_buf(:, class_cols);
                                space_buf_class_groupg = arrayfun(@(x) priogroup(1+x), space_buf_class);
                                start_classprio = min(space_buf_class_groupg(en_wbuf,:),[],2); % min finds highest priority
                                isrowmax = space_buf_class_groupg == repmat(start_classprio, 1, size(space_buf_class_groupg,2));
                                % LCFS: Find leftmost (first) class position for highest priority class
                                [~,leftmostClassPos]=max(isrowmax,[],2);
                                leftmostMaxPos = 2*leftmostClassPos - 1; % convert to buffer column index
                                start_svc_class = space_buf(en_wbuf, leftmostMaxPos); % job entering service
                                
                                % Handle states without buffer jobs
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni)
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];
                                
                                % Handle states with buffer jobs
                                if any(en_wbuf) && start_svc_class > 0
                                    % For LCFSPIPRIO, jobs restart from pie distribution instead of stored phase
                                    pentry_svc_class = pie{ist}{start_svc_class};
                                    for kentry = 1:K(start_svc_class)
                                        space_srv_k = space_srv;
                                        space_buf_k = space_buf;
                                        space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                        for j=find(en_wbuf)'
                                            % Remove both class and phase from leftmost position (preempt-independent ignores stored phase)
                                            % The (class,phase) pair buffer is RIGHT-aligned: occupied pairs sit at the
                                        % right, empty pairs pad the left. Priority can promote a pair from the
                                        % middle, so closing the hole must pad a whole empty PAIR on the left --
                                        % as State.afterEventStationSignal's dropWaiting does. Padding one slot
                                        % left and one right instead left a trailing empty slot, e.g.
                                        % [2 1 3 1] -> [0 3 1 0], a layout fromMarginal never enumerates: the
                                        % successor was unreachable, which is what made the generator reducible.
                                        space_buf_k(j,:) = [0, 0, space_buf_k(j,1:leftmostMaxPos(j)-1), space_buf_k(j,(leftmostMaxPos(j)+2):end)];
                                        end
                                        outspace = [outspace; space_buf_k(en_wbuf,:), space_srv_k(en_wbuf,:), space_var(en_wbuf,:)];
                                        outstart = State.tagArc(outstart, size(outspace,1), sum(en_wbuf), R, start_svc_class(:));
                                        rate_k = rate;
                                        rate_k(en_wbuf,:) = rate_k(en_wbuf,:) * pentry_svc_class(kentry);
                                        if isinf(ni)
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,end).*rate_k(en_wbuf,:)];
                                        else
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en_wbuf,:)];
                                        end
                                        outprob = [outprob; ones(size(rate_k(en_wbuf,:),1),1)];
                                    end
                                end
                            case SchedStrategy.SIRO
                                rate = zeros(size(space_srv,1),1);
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % this is for states not in en_buf
                                space_srv = inspace(:,(end-sum(K)-V+1):(end-V)); % server state (clear of the local vars)
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                % first record departure in states where the buffer is empty
                                en_wobuf = en & sum(space_buf(en,:),2) == 0;
                                if noPromote % immediate feedback: hold server, do not promote a waiting job
                                    en_wobuf = en;
                                end
                                outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                if isinf(ni) % hit limited load-dependence
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,end).*rate(en_wobuf,:)];
                                else
                                    outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wobuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate(en_wobuf,:)];
                                end
                                outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];
                                % let's go now to states where the buffer is non-empty
                                % (promotion suppressed under immediate feedback: R*(~noPromote)=0)
                                for r=1:(R*(~noPromote))
                                    rate_r = rate;
                                    space_buf = inspace(:,1:(end-sum(K)-V)); % buffer state (clear of the local vars)
                                    en_wbuf = en & space_buf(en,r) > 0; % states where the buffer is non-empty
                                    space_buf(en_wbuf,r) = space_buf(en_wbuf,r) - 1; % remove from buffer
                                    space_srv_r = space_srv;
                                    pentry_svc_class = pie{ist}{r};
                                    pick_prob = (nir(r)-sir(r)) / (ni-sum(sir));
                                    if pick_prob >= 0
                                        rate_r(en_wbuf,:) = rate_r(en_wbuf,:) * pick_prob;
                                    end
                                    for kentry=1:K(r)
                                        space_srv_r(en_wbuf,Ks(r)+kentry) = space_srv_r(en_wbuf,Ks(r)+kentry) + 1; % bring job in service
                                        outspace = [outspace; space_buf(en_wbuf,:), space_srv_r(en_wbuf,:), space_var(en_wbuf,:)];
                                        outstart = State.tagArc(outstart, size(outspace,1), sum(en_wbuf), R, r);
                                        rate_k = rate_r;
                                        rate_k(en_wbuf,:) = rate_k(en_wbuf,:) * pentry_svc_class(kentry);
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,end).*rate_k(en_wbuf,:)];
                                        else
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en_wbuf,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en_wbuf,:)];
                                        end
                                        outprob = [outprob; ones(size(rate(en_wbuf,:),1),1)];
                                        space_srv_r(en_wbuf,Ks(r)+kentry) = space_srv_r(en_wbuf,Ks(r)+kentry) - 1; % bring job in service
                                    end
                                end
                            case SchedStrategy.POLLING
                                % A completion ends the visit unless the
                                % discipline still allows another job of the
                                % same class to be taken; when it ends, the
                                % server walks the cyclic order to wherever the
                                % next tangible controller state lies. Which of
                                % the two happens depends on the buffer of each
                                % state row, so the rows are resolved one by one
                                % rather than as a single vectorized promotion.
                                % space_buf/space_srv are re-sliced from inspace by the
                                % enclosing k-loop; space_var deliberately is NOT, since
                                % it already carries the RROBIN/WRROBIN pointer advance
                                % applied above and re-slicing would discard it.
                                pinfoD = State.pollingInfo(sn, ind);
                                for rowD = find(en(:))'
                                    rateD = mu{ist}{class}(k)*phi{ist}{class}(k)*kir(rowD,class,k);
                                    if rateD <= 0
                                        continue
                                    end
                                    [~, swkD, ctrD] = State.pollingGet(pinfoD, space_var(rowD,:), class);
                                    if swkD ~= 0
                                        continue % no job can complete while the server is walking
                                    end
                                    bufD = space_buf(rowD,:);
                                    nbufD = bufD(1,1:R);
                                    srvD = space_srv(rowD,:);
                                    srvD(1,Ks(class)+k) = srvD(1,Ks(class)+k) - 1; % record departure
                                    switch pinfoD.ptype
                                        case PollingType.EXHAUSTIVE
                                            ctrnextD = 0;
                                            goonD = nbufD(class) > 0;
                                        case PollingType.GATED
                                            ctrnextD = ctrD - 1; % one of the gated jobs completed
                                            goonD = ctrnextD > 0;
                                        case PollingType.KLIMITED
                                            ctrnextD = ctrD - 1; % one of the K permitted services used
                                            goonD = ctrnextD > 0 && nbufD(class) > 0;
                                        case PollingType.DECREMENTING
                                            ctrnextD = ctrD; % the target level is fixed for the visit
                                            goonD = nbufD(class) > ctrD;
                                    end
                                    if goonD
                                        qD = class; modeD = 1; budgetD = ctrnextD;
                                    else
                                        [qD, modeD, budgetD] = State.pollingNext(pinfoD, class, nbufD, R, false);
                                    end
                                    [rowsD, probsD] = State.pollingLand(pinfoD, qD, modeD, budgetD, ...
                                        bufD, srvD, space_var(rowD,:), K, Ks, pie{ist}, R);
                                    for jD = 1:size(rowsD,1)
                                        outspace = [outspace; rowsD(jD,:)];
                                        % mode 1 pulls a waiting class-qD job into the server;
                                        % a switchover (2) or a park (0) starts nobody
                                        if modeD == 1
                                            outstart = State.tagArc(outstart, size(outspace,1), 1, R, qD);
                                        end
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,rowD,class).*lldscaling(ist,end).*rateD.*probsD(jD)];
                                        else
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,rowD,class).*lldscaling(ist,min(ni(rowD),lldlimit)).*rateD.*probsD(jD)];
                                        end
                                        outprob = [outprob; 1];
                                    end
                                end
                            case {SchedStrategy.SEPT,SchedStrategy.LEPT} % move last job in service
                                rate = zeros(size(space_srv,1),1);
                                rate(en) = mu{ist}{class}(k)*(phi{ist}{class}(k)).*kir(:,class,k); % this is for states not in en_buf
                                space_srv = inspace(:,(end-sum(K)-V+1):(end-V)); % server state (clear of the local vars)
                                space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                space_buf = inspace(:,1:(end-sum(K)-V)); % buffer state (clear of the local vars)
                                % in SEPT, the scheduling parameter is the priority order of the class means
                                % en_wbuf: states where the buffer is non-empty
                                % sept_class: class to pick in service
                                % sn.schedparam(ist,r) holds the RANK of class r's mean service
                                % time (ascending for SEPT, descending for LEPT, see
                                % MNetwork.refreshScheduling). Scanning the buffer columns in
                                % schedparam order visits the classes as P(1),P(2),..., which
                                % coincides with rank order only when P is an involution -- true
                                % of every 2-class model, but false in general, in which case the
                                % wrong class is promoted. The discipline needs the inverse map
                                % rank -> class, so invert the permutation first.
                                [~, classByRank] = sort(sn.schedparam(ist,1:R));
                                [en_wbuf, first_class_inrow] = max(space_buf(:,classByRank)~=0, [], 2);
                                sept_class = classByRank(first_class_inrow); % this is different for sept and lept
                                if noPromote % immediate feedback: hold server, do not promote a waiting job
                                    en_wbuf(:) = false;
                                end

                                space_buf(en_wbuf,sept_class) = space_buf(en_wbuf,sept_class) - 1; % remove from buffer
                                pentry = pie{ist}{sept_class};
                                for kentry=1:K(sept_class)
                                    space_srv(en_wbuf,Ks(sept_class)+kentry) = space_srv(en_wbuf,Ks(sept_class)+kentry) + 1; % bring job in service
                                    cls_d = zeros(sum(en),1);
                                    cls_d(en_wbuf(en)) = sept_class;
                                    if isSimulation
                                        % break the tie
                                        outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                        outstart = State.tagArc(outstart, size(outspace,1), sum(en), R, cls_d);
                                        rate_k = rate;
                                        rate_k(en,:) = rate_k(en,:) * pentry(kentry);
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate_k(en,:)];
                                        else
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en,:)];
                                        end
                                        outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    else
                                        outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                        outstart = State.tagArc(outstart, size(outspace,1), sum(en), R, cls_d);
                                        rate_k = rate;
                                        rate_k(en,:) = rate_k(en,:) * pentry(kentry);
                                        if isinf(ni) % hit limited load-dependence
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,end).*rate_k(en,:)];
                                        else
                                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,en,class).*lldscaling(ist,min(ni(en),lldlimit)).*rate_k(en,:)];
                                        end
                                        outprob = [outprob; ones(size(rate(en,:),1),1)];
                                    end
                                    space_srv(en_wbuf,Ks(sept_class)+kentry) = space_srv(en_wbuf,Ks(sept_class)+kentry) - 1; % bring job in service
                                end
                            otherwise
                                line_error(mfilename,sprintf('Scheduling strategy %s is not supported.', SchedStrategy.toText(sn.sched(ist))));
                        end
                    end
                end
                if isSimulation
                    outstart = State.tagPad(outstart, size(outspace,1), R);
                    outpreempt = State.tagPad(outpreempt, size(outspace,1), R);
                    if nargin>=7 && isobject(eventCache)
                        eventCache{key} = {outprob, outspace, outrate, outstart, outpreempt};
                    end

                    if size(outspace,1) > 1
                        tot_rate = sum(outrate);
                        cum_rate = cumsum(outrate) / tot_rate;
                        firing_ctr = 1 + max([0,find( rand > cum_rate' )]); % select action
                        outspace = outspace(firing_ctr,:);
                        outstart = outstart(firing_ctr,:); % the tags of the sampled arc
                        outpreempt = outpreempt(firing_ctr,:);
                        outrate = sum(outrate);
                        outprob = outprob(firing_ctr,:);
                    end
                end
            end
        end
    case EventType.PHASE
        outspace = [];
        outrate = [];
        outprob = [];
        [ni,nir,~,kir] = State.toMarginal(sn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
        if nir(class)>0
            for k=1:K(class)
                en = space_srv(:,Ks(class)+k) > 0;
                if any(en)
                    for kdest=setdiff(1:K(class),k) % new phase
                        rate = 0;
                        space_srv_k = space_srv(en,:);
                        space_buf_k = space_buf(en,:);
                        space_var_k = space_var(en,:);
                        if ismkvmodclass(class) && ~isempty(space_var_k)
                            space_var_k(sum(sn.nvars(ind,1:class))) = kdest;
                        end
                        space_srv_k(:,Ks(class)+k) = space_srv_k(:,Ks(class)+k) - 1;
                        space_srv_k(:,Ks(class)+kdest) = space_srv_k(:,Ks(class)+kdest) + 1;
                        switch sn.sched(ist)
                            case SchedStrategy.EXT
                                rate = proc{ist}{class}{1}(k,kdest); % move next job forward
                            case SchedStrategy.INF
                                rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k); % assume active
                            case {SchedStrategy.PS, SchedStrategy.LPS}
                                rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)./ni(:).*min(ni(:),S(ist)); % assume active
                            case SchedStrategy.PSPRIO
                                if all(ni <= S(ist)) || sn.classprio(class) == min(sn.classprio(nir>0))
                                    rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)./ni(:).*min(ni(:),S(ist)); % assume active
                                else
                                    rate = 0; % not in most urgent priority group
                                end
                            case SchedStrategy.DPSPRIO
                                if all(ni <= S(ist)) || sn.classprio(class) == min(sn.classprio(nir>0))
                                    w_i = sn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    nirprio = nir;
                                    if ~all(ni <= S(ist))
                                        nirprio(sn.classprio~=sn.classprio(class)) = 0;
                                    end
                                    rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)*w_i(class)./(sum(repmat(w_i,size(nirprio,1),1)*nirprio',2)); % assume active
                                else
                                    rate = 0; % not in most urgent priority group
                                end
                            case SchedStrategy.GPSPRIO
                                if all(ni <= S(ist)) || sn.classprio(class) == min(sn.classprio(nir>0))
                                    w_i = sn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    nirprio = nir;
                                    if ~all(ni <= S(ist))
                                        nirprio(sn.classprio~=sn.classprio(class)) = 0;
                                    end
                                    cir = min(nirprio,ones(size(nirprio)));
                                    rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)/nirprio(class)*w_i(class)/(w_i*cir(:)); % assume active
                                else
                                    rate = 0; % not in most urgent priority group
                                end
                            case SchedStrategy.DPS
                                if S(ist) > 1
                                    line_error(mfilename,'Multi-server DPS not supported yet');
                                end
                                w_i = sn.schedparam(ist,:);
                                w_i = w_i / sum(w_i);
                                rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)*w_i(class)./(sum(repmat(w_i,size(nir,1),1)*nir',2)); % assume active
                            case SchedStrategy.GPS
                                if S(ist) > 1
                                    line_error(mfilename,'Multi-server GPS not supported yet');
                                end
                                cir = min(nir,ones(size(nir)));
                                w_i = sn.schedparam(ist,:); w_i = w_i / sum(w_i);
                                rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k)/nir(class)*w_i(class)/(w_i*cir(:)); % assume active

                            % 2026-07-29: the preempt family completed here. Only
                            % LCFSPR/LCFSPI were listed, so a job in service at
                            % FCFSPR, FCFSPI and PRIO variants could never advance its service
                            % phase: no PHASE transition was emitted, the states
                            % with a busy server were unreachable, and the
                            % generator collapsed onto an absorbing state
                            % ("no recurrent state" from ctmc_solve). Buffered
                            % jobs are frozen under both PR and PI, so the rate
                            % is the same one every non-sharing discipline uses.
                            case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS, SchedStrategy.LCFSPRIO, SchedStrategy.SIRO, SchedStrategy.SEPT, SchedStrategy.LEPT, SchedStrategy.POLLING, ...
                                    SchedStrategy.LCFSPR, SchedStrategy.LCFSPI, SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO, ...
                                    SchedStrategy.FCFSPR, SchedStrategy.FCFSPI, SchedStrategy.FCFSPRPRIO, SchedStrategy.FCFSPIPRIO}
                                rate = proc{ist}{class}{1}(k,kdest)*kir(:,class,k); % assume active
                        end
                        % if the class cannot be served locally,
                        % then rate = NaN since mu{i,class}=NaN
                        if isinf(ni) % hit limited load-dependence
                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,1:size(rate,1),class).*lldscaling(ist,end).*rate];
                        else
                            outrate = [outrate; State.cdclassfactor(cdscaling{ist},nir,1:size(rate,1),class).*lldscaling(ist,min(ni(en),lldlimit)).*rate];
                        end
                        outspace = [outspace; space_buf_k, space_srv_k, space_var_k];
                        outprob = [outprob; ones(size(rate,1),1)];
                    end
                end
            end
            if isSimulation
                if nargin>=7 && isobject(eventCache)
                    eventCache{key} = {outprob, outspace, outrate, outstart, outpreempt};
                end

                if size(outspace,1) > 1
                    tot_rate = sum(outrate);
                    cum_rate = cumsum(outrate) / tot_rate;
                    firing_ctr = 1 + max([0,find( rand > cum_rate' )]); % select action
                    outspace = outspace(firing_ctr,:);
                    outrate = sum(outrate);
                    outprob = outprob(firing_ctr,:);
                end
            end
        end
    case EventType.SWITCH
        % Switchover of a polling server walking towards buffer `class`. Unlike
        % PHASE, which carries only the internal transitions of a phase-type
        % and leaves the absorption to DEP, this event carries both: a
        % completed switchover moves no job and so has no departure to attach
        % the absorption to (see refreshSync). It is therefore also emitted for
        % a single-phase switchover, where it consists of the absorption alone.
        outspace = [];
        outrate = [];
        outprob = [];
        pinfoS = State.pollingInfo(sn, ind);
        if isempty(pinfoS) || ~pinfoS.hasSw(class)
            return
        end
        for rowS = 1:size(inspace,1)
            [posS, swkS] = State.pollingGet(pinfoS, space_var(rowS,:), 0);
            if posS ~= class || swkS == 0
                continue % the server is not inside the switchover into `class`
            end
            D0S = pinfoS.swD0{class};
            D1S = pinfoS.swD1{class};
            % internal transitions of the switchover phase-type
            for kdestS = setdiff(1:pinfoS.Ksw(class), swkS)
                if D0S(swkS,kdestS) <= 0
                    continue
                end
                varS = State.pollingSet(pinfoS, space_var(rowS,:), class, kdestS, 0);
                outspace = [outspace; space_buf(rowS,:), space_srv(rowS,:), varS];
                outrate = [outrate; D0S(swkS,kdestS)];
                outprob = [outprob; 1];
            end
            % absorption: the server arrives at buffer `class` and either opens
            % a visit there or walks on
            rateS = sum(D1S(swkS,:));
            if rateS <= 0
                continue
            end
            nbufS = space_buf(rowS,1:R);
            [qS, modeS, budgetS] = State.pollingNext(pinfoS, class, nbufS, R, true);
            [rowsS, probsS] = State.pollingLand(pinfoS, qS, modeS, budgetS, ...
                space_buf(rowS,:), space_srv(rowS,:), space_var(rowS,:), K, Ks, pie{ist}, R);
            for jS = 1:size(rowsS,1)
                % A switchover that completes over an empty buffer starts the
                % next leg at once, and when that leg re-enters the same phase
                % of the same switchover (the single-buffer case, or any
                % memoryless switchover) the landing state is the departure
                % state. Such a self-loop is not a transition: emitting it would
                % add a spurious rate to the row and inflate the exit rate.
                if isequal(rowsS(jS,:), inspace(rowS,:))
                    continue
                end
                outspace = [outspace; rowsS(jS,:)];
                outrate = [outrate; rateS*probsS(jS)];
                outprob = [outprob; 1];
                % A completed switchover that opens a visit (mode 1) pulls a
                % waiting class-qS job into the server, so it starts service
                % just as an ARV or a DEP promotion does. This is the one
                % service start a polling station reaches through neither, and
                % leaving it untagged would break startRate == TN + preemptRate
                % there for no reason other than the name of the carrier event.
                if modeS == 1
                    outstart = State.tagArc(outstart, size(outspace,1), 1, R, qS);
                end
            end
        end
    case EventType.RENEGE
        % Exponential-patience reneging: each waiting (queued, not-in-service)
        % class-r job abandons the queue at a memoryless rate impatienceMu, so
        % the aggregate rate out of this state is (waiting count) * mu. One
        % waiting job is removed from the buffer and leaves the system (the
        % passive half of the sync is LOCAL). Removing the newest waiting
        % class-r job and re-padding a zero keeps the buffer in the canonical
        % right-aligned form used by the arrival handler; for memoryless
        % patience all waiting jobs are exchangeable, so the choice does not
        % affect the marginal distribution.
        outspace = [];
        outrate = [];
        outprob = [];
        [~,nir,sir] = State.toMarginal(sn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
        waiting_r = nir(class) - sir(class);
        if waiting_r > 0
            slot = find(space_buf == class, 1, 'first');
            if ~isempty(slot)
                space_buf_k = space_buf;
                space_buf_k(slot) = [];
                space_buf_k = [0, space_buf_k];
                outspace = [space_buf_k, space_srv, space_var];
                outrate = waiting_r * sn.impatienceMu(ist,class);
                outprob = 1;
            end
        end
    case EventType.RETRY
        % Exponential retrial: an orbiting (buffered) class-r job retries entry
        % into the station. The retry succeeds only when a server is free, in
        % which case one orbiting job enters service (in an entry phase drawn
        % from pie); otherwise the job stays in orbit and the event is a no-op
        % (not generated). The aggregate rate is (orbit size) * retrialMu.
        outspace = [];
        outrate = [];
        outprob = [];
        [~,nir,sir] = State.toMarginal(sn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
        orbit_r = nir(class) - sir(class); % orbiting class-r jobs (held in buffer)
        if orbit_r > 0 && sum(space_srv,2) < S(ist)
            slot = find(space_buf == class, 1, 'first');
            if ~isempty(slot)
                space_buf_k = space_buf;
                space_buf_k(slot) = [];
                space_buf_k = [0, space_buf_k];
                pentry = pie{ist}{class};
                if all(isnan(pentry))
                    pentry = ones(size(pentry)) / length(pentry);
                end
                for kentry = 1:K(class)
                    if pentry(kentry) <= 0
                        continue;
                    end
                    space_srv_k = space_srv;
                    space_srv_k(:,Ks(class)+kentry) = space_srv_k(:,Ks(class)+kentry) + 1;
                    % LINEAR policy: every orbiting job carries its own timer, so
                    % the aggregate retrial rate scales with the orbit size.
                    % CONSTANT policy: one controller retries on behalf of the
                    % whole orbit, so the rate does not depend on the orbit size.
                    retrialRate = orbit_r * sn.retrialMu(ist,class);
                    if isfield(sn,'retrialPolicy') && ~isempty(sn.retrialPolicy) ...
                            && size(sn.retrialPolicy,1) >= ist && size(sn.retrialPolicy,2) >= class ...
                            && sn.retrialPolicy(ist,class) == RetrialPolicy.CONSTANT
                        retrialRate = sn.retrialMu(ist,class);
                    end
                    outspace = [outspace; space_buf_k, space_srv_k, space_var];
                    outrate  = [outrate; retrialRate * pentry(kentry)];
                    outprob  = [outprob; 1];
                    % A successful retry is the only way into the server at a
                    % retrial station (a DEP there never promotes from the
                    % orbit), so it carries the START the invariant needs.
                    outstart = State.tagArc(outstart, size(outspace,1), 1, R, class);
                end
                if isSimulation && size(outspace,1) > 1
                    cr = cumsum(outrate) / sum(outrate);
                    fc = 1 + max([0, find(rand > cr')]);
                    outstart = State.tagPad(outstart, size(outspace,1), R);
                    outpreempt = State.tagPad(outpreempt, size(outspace,1), R);
                    outspace = outspace(fc,:);
                    outstart = outstart(fc,:);
                    outpreempt = outpreempt(fc,:);
                    outrate = sum(outrate);
                    outprob = 1;
                end
            end
        end
end

% True BAS: when the front job is already blocked (completed, held at the server), this
% DEP is the *instant transfer* of that job downstream — fire at rate 1e7 (effectively
% instant) and clear the blocked marker in the successor. The complementary become-blocked
% transition (b:0->1 when the destination is full) is added by the CTMC generator / SSA
% engine, since only they can see the destination's occupancy. Blocked marker = last state
% column when this station declares it (nvars col 2*R+1 == 1).
% The blocked marker shares nvars col 2*R+1 with the polling controller, so gate on
% the dedicated sn.isbasblocking field (set for the blocking station under BOTH the
% upstream and destination declaration forms) rather than the station's own drop
% rule, which fails for a destination-declared BAS. See BUG-83.
if event == EventType.DEP && ~isempty(sn.isbasblocking) && numel(sn.isbasblocking) >= ind ...
        && sn.isbasblocking(ind) == 1 ...
        && ~isempty(inspace) && inspace(1,end) == 1 && ~isempty(outspace)
    outspace(:,end) = 0;
    outrate(:) = 1e7;
end

% Degraded service while the server is down: the scheduling handlers computed the
% completion rate from the up-server service process, so rescale it to the
% configured down-server rate. downRateScale is 1 whenever the server is up or no
% degraded rate was configured, so this is a no-op in every other model.
if downRateScale ~= 1 && ~isempty(outrate)
    outrate = outrate * downRateScale;
end

end