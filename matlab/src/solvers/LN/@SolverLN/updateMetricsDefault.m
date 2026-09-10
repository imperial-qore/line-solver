function updateMetricsDefault(self, it)
ensemble = self.ensemble;
lqn = self.lqn;

% obtain the activity service times
self.servt = zeros(lqn.nidx,1);
self.residt = zeros(lqn.nidx,1);
for r=1:size(self.servt_classes_updmap,1)
    idx = self.servt_classes_updmap(r,1);
    aidx = self.servt_classes_updmap(r,2);
    nodeidx = self.servt_classes_updmap(r,3);
    classidx = self.servt_classes_updmap(r,4);

    % store the residence times and tput at this layer to become
    % the servt / tputs of aidx in another layer, as needed
    iter_min = min(30,ceil(self.options.iter_max/4));
    wnd_size = (it-self.averagingstart+1);

    % Compute residt from QN/TN_ref instead of WN to avoid
    % fork+loop visit distortion (WN uses visits from DTMC solve
    % which are distorted when Fork non-stochastic rows coexist
    % with loop back-edges in the routing matrix)
    layerIdx = self.idxhash(idx);
    layerSn = ensemble{layerIdx}.getStruct();
    % refclass is 0 on a chain with no reference class (open chain): no TN_ref
    [hasRef, refstat_k, refclass_c] = ln_layer_refcell(layerSn, classidx);

    if ~isempty(self.averagingstart) && it>=iter_min % assume steady-state
        self.servt(aidx) = 0;
        self.residt(aidx) = 0;
        self.tput(aidx) = 0;
        for w=0:(wnd_size-1)
            self.servt(aidx) = self.servt(aidx) + self.results{end-w,layerIdx}.RN(nodeidx,classidx) / wnd_size;
            if hasRef
                TN_ref = self.results{end-w,layerIdx}.TN(refstat_k, refclass_c);
            else
                TN_ref = 0;
            end
            if TN_ref > GlobalConstants.FineTol
                self.residt(aidx) = self.residt(aidx) + self.results{end-w,layerIdx}.QN(nodeidx,classidx) / TN_ref / wnd_size;
            else
                self.residt(aidx) = self.residt(aidx) + self.results{end-w,layerIdx}.WN(nodeidx,classidx) / wnd_size;
            end
            self.tput(aidx) = self.tput(aidx) + self.results{end-w,layerIdx}.TN(nodeidx,classidx) / wnd_size;
        end
    else
        self.servt(aidx) = self.results{end,layerIdx}.RN(nodeidx,classidx);
        if hasRef
            TN_ref = self.results{end,layerIdx}.TN(refstat_k, refclass_c);
        else
            TN_ref = 0;
        end
        QN_val = self.results{end,layerIdx}.QN(nodeidx,classidx);
        if TN_ref > GlobalConstants.FineTol
            self.residt(aidx) = QN_val / TN_ref;
        else
            self.residt(aidx) = self.results{end,layerIdx}.WN(nodeidx,classidx);
        end
        self.tput(aidx) = self.results{end,layerIdx}.TN(nodeidx,classidx);
    end

    % see _kb/06-solver-catalog.md (LN section) for rationale
    zt_act = lqn_act_thinktime(lqn, aidx);
    if zt_act > 0
        self.servt(aidx) = self.servt(aidx) + zt_act;
        self.residt(aidx) = self.residt(aidx) + zt_act;
    end

    % see _kb/06-solver-catalog.md (LN section) for rationale
    if aidx > lqn.ashift && aidx <= lqn.ashift + lqn.nacts
        % This is an activity - find its bound entry
        for eidx = (lqn.eshift+1):(lqn.eshift+lqn.nentries)
            if full(lqn.graph(eidx, aidx)) > 0
                % Found bound entry - check if async-only
                hasSyncCallers = full(any(lqn.issynccaller(:, eidx)));
                hasAsyncCallers = full(any(lqn.isasynccaller(:, eidx)));
                if hasAsyncCallers && ~hasSyncCallers
                    % Async-only target: use RN (response time per visit)
                    % instead of WN (residence time with visit ratio)
                    self.residt(aidx) = self.servt(aidx);  % servt already has RN
                end
                break;
            end
        end
    end

    % Recover from Inf/NaN: snap back to previous iteration's value
    if it > 1
        if (isinf(self.servt(aidx)) || isnan(self.servt(aidx))) && ~isnan(self.servt_prev(aidx))
            self.servt(aidx) = self.servt_prev(aidx);
        end
        if (isinf(self.residt(aidx)) || isnan(self.residt(aidx))) && ~isnan(self.residt_prev(aidx))
            self.residt(aidx) = self.residt_prev(aidx);
        end
        if (isinf(self.tput(aidx)) || isnan(self.tput(aidx))) && ~isnan(self.tput_prev(aidx))
            self.tput(aidx) = self.tput_prev(aidx);
        end
    end

    % Apply under-relaxation if enabled and not first iteration
    omega = self.relax_omega;
    if omega < 1.0 && it > 1
        if ~isnan(self.servt_prev(aidx))
            self.servt(aidx) = omega * self.servt(aidx) + (1 - omega) * self.servt_prev(aidx);
        end
        if ~isnan(self.residt_prev(aidx))
            self.residt(aidx) = omega * self.residt(aidx) + (1 - omega) * self.residt_prev(aidx);
        end
        if ~isnan(self.tput_prev(aidx))
            self.tput(aidx) = omega * self.tput(aidx) + (1 - omega) * self.tput_prev(aidx);
        end
    end
    % Store current values for next iteration
    self.servt_prev(aidx) = self.servt(aidx);
    self.residt_prev(aidx) = self.residt(aidx);
    self.tput_prev(aidx) = self.tput(aidx);

    % Safeguard against MVA numerical instability producing extreme values
    % (matches Python _update_metrics_default max_servt guard)
    max_servt = 1e10;
    if self.servt(aidx) > 0 && self.servt(aidx) <= max_servt
        self.servtproc{aidx} = Exp.fitMean(self.servt(aidx));
    end
    self.tputproc{aidx} = Exp.fitRate(self.tput(aidx));
end

% Phase-2 support: split activity service times by phase
% Note: Overtaking probability is computed later after entry throughput is available
if self.hasPhase2
    % Reset phase-specific arrays
    self.servt_ph1 = zeros(lqn.nidx, 1);
    self.servt_ph2 = zeros(lqn.nidx, 1);

    % Split activity service times by phase
    for a = 1:lqn.nacts
        aidx = lqn.ashift + a;
        if lqn.actphase(a) == 1
            self.servt_ph1(aidx) = self.servt(aidx);
        else
            self.servt_ph2(aidx) = self.servt(aidx);
        end
    end

    % Aggregate phase service times to entry level
    for e = 1:lqn.nentries
        eidx = lqn.eshift + e;
        acts = lqn.actsof{eidx};
        for aidx = acts
            a = aidx - lqn.ashift;
            if a > 0 && a <= lqn.nacts
                if lqn.actphase(a) == 1
                    self.servt_ph1(eidx) = self.servt_ph1(eidx) + self.servt_ph1(aidx);
                else
                    self.servt_ph2(eidx) = self.servt_ph2(eidx) + self.servt_ph2(aidx);
                end
            end
        end
    end
end

% obtain throughput for activities in thinkt_classes_updmap (needed for async calls)
% this ensures tputproc is set for activities that make async calls from client nodes
for r=1:size(self.thinkt_classes_updmap,1)
    idx = self.thinkt_classes_updmap(r,1);
    aidx = self.thinkt_classes_updmap(r,2);
    nodeidx = self.thinkt_classes_updmap(r,3);
    classidx = self.thinkt_classes_updmap(r,4);

    % only update if not already set by servt_classes_updmap processing
    if isempty(self.tputproc) || length(self.tputproc) < aidx || isempty(self.tputproc{aidx})
        iter_min = min(30,ceil(self.options.iter_max/4));
        wnd_size = (it-self.averagingstart+1);
        if ~isempty(self.averagingstart) && it>=iter_min % assume steady-state
            self.tput(aidx) = 0;
            for w=0:(wnd_size-1)
                self.tput(aidx) = self.tput(aidx) + self.results{end-w,self.idxhash(idx)}.TN(nodeidx,classidx) / wnd_size;
            end
        else
            self.tput(aidx) = self.results{end,self.idxhash(idx)}.TN(nodeidx,classidx);
        end
        self.tputproc{aidx} = Exp.fitRate(self.tput(aidx));
    end
end

% obtain the call residence time
self.callservt = zeros(lqn.ncalls,1);
self.callresidt = zeros(lqn.ncalls,1);
% A call into a replicated callee registers ONE ROW PER REPLICA STATION, and
% every row rewrites the same entry, so a residence built by Little from that
% station's queue keeps one replica's share of it. The replicas are symmetric --
% same law, same 1/n routing -- so the whole call is that share times the number
% of rows. Without this, an entry whose activity calls a task replicated r times
% reported a response time r times short of the activity's own: 0.0942105 against
% 0.182631 at r=3, where LDES measures 0.182575. One row leaves it unchanged.
% Count DISTINCT stations, not rows: the same call is registered once per
% position it is reached from in the activity graph, and counting those would
% inflate every forked or looped call in a model with no replication at all.
callRowsPerLayer = ones(size(self.call_classes_updmap,1),1);
for r=1:size(self.call_classes_updmap,1)
    sib = self.call_classes_updmap(:,1) == self.call_classes_updmap(r,1) ...
        & self.call_classes_updmap(:,2) == self.call_classes_updmap(r,2) ...
        & self.call_classes_updmap(:,3) > 1;
    callRowsPerLayer(r) = numel(unique(self.call_classes_updmap(sib,3)));
end
for r=1:size(self.call_classes_updmap,1)
    idx = self.call_classes_updmap(r,1);
    cidx = self.call_classes_updmap(r,2);
    nodeidx = self.call_classes_updmap(r,3);
    classidx = self.call_classes_updmap(r,4);
    if self.call_classes_updmap(r,3) > 1
        if nodeidx == 1
            self.callservt(cidx) = 0;
        else
            % A job blocked at an admission constraint is counted at no station
            % (JMT WAITQ convention), so its wait is absent from RN. Recover it
            % by Little from the layer population deficit, otherwise the caller
            % never sees the blocking and the fixed point loses flow balance.
            % The population is conserved per chain, not per class: a job in a
            % layer switches class along the activity graph, so the call class
            % itself carries population 0. Split the chain deficit by
            % throughput, which gives every region-visiting class the same wait.
            fcrWait = 0;
            eidxLayer = self.idxhash(idx);
            if ~isempty(self.layerHasRegion) && self.layerHasRegion(eidxLayer)
                chains = self.layerChains{eidxLayer};
                cix = find(chains(:,classidx), 1);
                if ~isempty(cix)
                    chainClasses = find(chains(cix,:));
                    QNlayer = self.results{end, eidxLayer}.QN;
                    TNlayer = self.results{end, eidxLayer}.TN;
                    chainPop = 0;
                    for k = chainClasses
                        chainPop = chainPop + self.ensemble{eidxLayer}.classes{k}.population;
                    end
                    deficit = chainPop - sum(sum(QNlayer(:,chainClasses)));
                    xregion = sum(TNlayer(nodeidx,chainClasses));
                    if isfinite(deficit) && deficit > 0 && xregion > GlobalConstants.FineTol
                        fcrWait = deficit / xregion;
                    end
                end
            end
            self.callservt(cidx) = (self.results{end, self.idxhash(idx)}.RN(nodeidx,classidx) + fcrWait) * self.lqn.callproc{cidx}.getMean;
            % Normalise per chain-reference visit, as residt does above. WN
            % divides by the class's own reference rate when the layer is open
            % (an INF client task), which is per-ENTRY visit, and the entry
            % rescaling below would then count the call once per entry.
            callSn = ensemble{eidxLayer}.getStruct();
            % refclass is 0 on a chain with no reference class (open chain)
            [hasCallRef, callRefstat, callRefclass] = ln_layer_refcell(callSn, classidx);
            TN_ref = 0;
            if hasCallRef
                TN_ref = self.results{end, eidxLayer}.TN(callRefstat, callRefclass);
            end
            if TN_ref > GlobalConstants.FineTol
                self.callresidt(cidx) = max(1,callRowsPerLayer(r)) * ...
                    self.results{end, eidxLayer}.QN(nodeidx,classidx) / TN_ref + fcrWait;
            else
                self.callresidt(cidx) = self.results{end, eidxLayer}.WN(nodeidx,classidx) + fcrWait;
            end
        end
        % Recover from Inf/NaN (e.g. a transiently unstable open chain in a
        % layer): snap back to the previous iteration's value, otherwise the
        % under-relaxation below makes Inf absorbing.
        if (isinf(self.callservt(cidx)) || isnan(self.callservt(cidx)))
            if it > 1 && isfinite(self.callservt_prev(cidx))
                self.callservt(cidx) = self.callservt_prev(cidx);
            else
                self.callservt(cidx) = 0;
            end
        end
        if (isinf(self.callresidt(cidx)) || isnan(self.callresidt(cidx)))
            if it > 1 && isfinite(self.callresidt_prev(cidx))
                self.callresidt(cidx) = self.callresidt_prev(cidx);
            else
                self.callresidt(cidx) = 0;
            end
        end
        % Growth rate capping removed - it prevents callservt from converging
        % to the correct value when initial values are near-zero (Immediate)
        % Apply under-relaxation to call service times
        omega = self.relax_omega;
        if omega < 1.0 && it > 1 && ~isnan(self.callservt_prev(cidx))
            self.callservt(cidx) = omega * self.callservt(cidx) + (1 - omega) * self.callservt_prev(cidx);
        end
        self.callservt_prev(cidx) = self.callservt(cidx);
        self.callresidt_prev(cidx) = self.callresidt(cidx);
    end
end

% Obtain the join times for AND-Join activities. This runs after callresidt is
% known: a branch activity with an Immediate host demand spends all its time in
% its synchronous calls, so a branch time read from residt alone would be zero
% and the join correction would silently vanish.
% see _kb/06-solver-catalog.md (LN section) for rationale
self.joint = zeros(lqn.nidx,1);
joint_excess = zeros(lqn.nidx,1);
% PRE_AND marks the branch tails, not the join target, so the joins are the activities
% whose predecessors carry that mark.
branchtails = find(lqn.actpretype == ActivityPrecedenceType.PRE_AND)';
joinedacts = [];
for tailidx = branchtails
    succs = find(lqn.graph(tailidx, :) > 0);
    joinedacts = [joinedacts, succs]; %#ok<AGROW>
end
joinedacts = unique(joinedacts);
joinedacts = joinedacts(joinedacts > lqn.ashift & joinedacts <= lqn.ashift + lqn.nacts);
for aidx = joinedacts
    branches = fj_branch_members(lqn, aidx);
    nbranches = numel(branches);
    if nbranches == 0
        continue;
    end
    branch_times = zeros(1, nbranches);
    for bi = 1:nbranches
        branch_times(bi) = sum(self.residt(branches{bi}));
        for baidx = branches{bi}(:)'
            for cidx = lqn.callsof{baidx}
                if lqn.calltype(cidx) == CallType.SYNC
                    branch_times(bi) = branch_times(bi) + self.callresidt(cidx);
                end
            end
        end
    end
    if nbranches == 1
        self.joint(aidx) = branch_times(1);
        continue;
    end
    quorum = nbranches;
    if isfield(lqn, 'actquorum') && aidx <= length(lqn.actquorum)
        q = full(lqn.actquorum(aidx));
        if q >= 1 && q <= nbranches
            quorum = q;
        end
    end
    % Branch times are taken as exponential, so the variance is the square of the mean.
    self.joint(aidx) = fj_quorum_moments(branch_times, branch_times.^2, quorum);
    joint_excess(aidx) = self.joint(aidx) - sum(branch_times);
end

% then resolve the entry servt summing up these contributions
entry_servt = self.servtmatrix*[self.residt;self.callresidt(:)];
entry_servt(1:lqn.eshift) = 0;

% see _kb/06-solver-catalog.md (LN section) for rationale
joins_with_excess = find(joint_excess ~= 0)';
for eidx = (lqn.eshift+1):(lqn.eshift+lqn.nentries)
    for aidx = joins_with_excess
        if self.servtmatrix(eidx, aidx) > 0
            entry_servt(eidx) = entry_servt(eidx) + joint_excess(aidx);
        end
    end
    entry_servt(eidx) = max(entry_servt(eidx), 0);
end

% A SetupTask's cold start is charged HERE, to the entry, and with the
% probability that the thread was actually found powered down. It is not host
% demand, so it does not belong to any activity's residence -- reporting it there
% put RespT(A2) at 1.29479 against 0.333178 from LDES on lqn_setup, the bare
% demand. The probability is the same closure method 'srvn.ph' uses (phSetupProb):
% a thread released at a reply powers off only if its delay-off countdown D
% expires before the next request, so p = E[I]/(E[I]+d), and admission takes an
% ACTIVE thread before waking a sleeping one, so the pool that cycles is
% max(1,b) threads at offered load b. Charging the FULL setup mean every time,
% which the open QBD decomposition did, ran lqn_setup 12.67% below LDES.
for eidx = (lqn.eshift+1):(lqn.eshift+lqn.nentries)
    entry_servt(eidx) = entry_servt(eidx) + lqn_setup_charge(self, lqn, lqn.parent(eidx));
end

% see _kb/06-solver-catalog.md (LN section) for rationale



% this block fixes the problem that ResidT is scaled so that the
% task has Vtask=1, but in call servt the entries need to have Ventry=1
for eidx=(lqn.eshift+1):(lqn.eshift+lqn.nentries)
    tidx = lqn.parent(eidx); % task of entry
    hidx = lqn.parent(tidx); %host of entry
    if ~self.ignore(tidx) && ~self.ignore(hidx)
        % Check if this entry has sync callers (which create closed classes)
        hasSyncCallers = full(any(lqn.issynccaller(:, eidx)));

        if hasSyncCallers
            % Original logic for entries with sync callers
            % get class in host layer of task and entry
            tidxclass = ensemble{self.idxhash(hidx)}.attribute.tasks(find(ensemble{self.idxhash(hidx)}.attribute.tasks(:,2) == tidx),1);
            eidxclass = ensemble{self.idxhash(hidx)}.attribute.entries(find(ensemble{self.idxhash(hidx)}.attribute.entries(:,2) == eidx),1);
            task_tput  = sum(self.results{end,self.idxhash(hidx)}.TN(ensemble{self.idxhash(hidx)}.attribute.clientIdx,tidxclass));
            entry_tput = sum(self.results{end,self.idxhash(hidx)}.TN(ensemble{self.idxhash(hidx)}.attribute.clientIdx,eidxclass));
            if entry_tput > GlobalConstants.Zero
                self.servt(eidx) = entry_servt(eidx) * task_tput / entry_tput;
                self.residt(eidx) = entry_servt(eidx) * task_tput / entry_tput;
            else
                self.servt(eidx) = entry_servt(eidx);
                self.residt(eidx) = entry_servt(eidx);
            end
        else
            % For async-only targets, use entry_servt directly
            % No throughput ratio scaling needed since there are no closed classes
            self.servt(eidx) = entry_servt(eidx);
            self.residt(eidx) = entry_servt(eidx);
        end
    end
end

% Phase-2 support: compute overtaking probability and apply correction
% This must happen AFTER entry throughput is available (computed above)
if self.hasPhase2
    for e = 1:lqn.nentries
        eidx = lqn.eshift + e;
        tidx = lqn.parent(eidx);

        if self.servt_ph2(eidx) > GlobalConstants.FineTol
            if lqn.isref(tidx) || ~full(any(lqn.issynccaller(:, eidx)))
                self.residt(eidx) = self.servt(eidx);
                continue;
            end
            % Get entry throughput (use task throughput as approximation if entry not available)
            if self.tput(eidx) > GlobalConstants.FineTol
                entry_tput = self.tput(eidx);
            elseif self.tput(tidx) > GlobalConstants.FineTol
                entry_tput = self.tput(tidx);
            else
                entry_tput = 0;
            end

            % Compute overtaking probability now that throughput is available.
            % Layer-1 LQNS phased-server Markov chain (aligns with lqns -t
            % overtaking); see overtake_prob_markov / lqn_overtake_markov.
            if entry_tput > GlobalConstants.FineTol
                self.prOvertake(e) = self.overtake_prob_markov(eidx);
            else
                self.prOvertake(e) = 0;
            end

            % Caller's response time = phase-1 only + P(overtake) * phase-2
            % Phase-2 only delays if overtaking occurs
            overtake_delay = self.prOvertake(e) * self.servt_ph2(eidx);

            % The caller sees phase-1 + overtaking correction (not full phase-2)
            % self.servt(eidx) remains unchanged (phase-1 + phase-2) for utilization calculation
            self.residt(eidx) = self.servt_ph1(eidx) + overtake_delay;
        end
    end
end

%self.servt(lqn.eshift+1:lqn.eshift+lqn.nentries) = entry_servt(lqn.eshift+1:lqn.eshift+lqn.nentries);
%entry_servt((lqn.ashift+1):end) = 0;
for r=1:size(self.call_classes_updmap,1)
    cidx = self.call_classes_updmap(r,2);
    eidx = lqn.callpair(cidx,2);
    if self.call_classes_updmap(r,3) > 1
        if self.servt(eidx) > 0
            self.servtproc{eidx} = Exp.fitMean(self.servt(eidx));
        end
    end
end

% determine call response times processes
for r=1:size(self.call_classes_updmap,1)
    cidx = self.call_classes_updmap(r,2);
    eidx = lqn.callpair(cidx,2);
    if self.call_classes_updmap(r,3) > 1
        if it==1
            % note that respt is per visit, so number of calls is 1
            self.callservt(cidx) = self.servt(eidx);
            self.callservtproc{cidx} = self.servtproc{eidx};
        else
            % note that respt is per visit, so number of calls is 1
            if self.callservt(cidx) > 0
                self.callservtproc{cidx} = Exp.fitMean(self.callservt(cidx));
            end
        end
    end
end

% see _kb/06-solver-catalog.md (LN section) for rationale
if self.hasPhase2
    for cidx = 1:lqn.ncalls
        if lqn.calltype(cidx) ~= CallType.SYNC
            continue;
        end
        target_eidx = lqn.callpair(cidx, 2);
        e_tgt = target_eidx - lqn.eshift;
        if e_tgt >= 1 && e_tgt <= lqn.nentries ...
                && self.servt_ph2(target_eidx) > GlobalConstants.FineTol
            callmean = lqn.callproc{cidx}.getMean();
            eff = self.residt(target_eidx);   % servt_ph1 + prOt * servt_ph2
            if eff > 0
                self.callservt(cidx)     = eff * callmean;
                self.callresidt(cidx)    = eff * callmean;
                self.callservtproc{cidx} = Exp.fitMean(eff * callmean);
            end
        end
    end
end

self.ensemble = ensemble;
end
