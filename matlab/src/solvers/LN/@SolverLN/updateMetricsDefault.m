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
    c = find(layerSn.chains(:, classidx), 1);
    refclass_c = layerSn.refclass(c);
    refstat_k = layerSn.refstat(classidx);

    if ~isempty(self.averagingstart) && it>=iter_min % assume steady-state
        self.servt(aidx) = 0;
        self.residt(aidx) = 0;
        self.tput(aidx) = 0;
        for w=1:(wnd_size-1)
            self.servt(aidx) = self.servt(aidx) + self.results{end-w,layerIdx}.RN(nodeidx,classidx) / wnd_size;
            TN_ref = self.results{end-w,layerIdx}.TN(refstat_k, refclass_c);
            if TN_ref > GlobalConstants.FineTol
                self.residt(aidx) = self.residt(aidx) + self.results{end-w,layerIdx}.QN(nodeidx,classidx) / TN_ref / wnd_size;
            else
                self.residt(aidx) = self.residt(aidx) + self.results{end-w,layerIdx}.WN(nodeidx,classidx) / wnd_size;
            end
            self.tput(aidx) = self.tput(aidx) + self.results{end-w,layerIdx}.TN(nodeidx,classidx) / wnd_size;
        end
    else
        self.servt(aidx) = self.results{end,layerIdx}.RN(nodeidx,classidx);
        TN_ref = self.results{end,layerIdx}.TN(refstat_k, refclass_c);
        if TN_ref > GlobalConstants.FineTol
            self.residt(aidx) = self.results{end,layerIdx}.QN(nodeidx,classidx) / TN_ref;
        else
            self.residt(aidx) = self.results{end,layerIdx}.WN(nodeidx,classidx);
        end
        self.tput(aidx) = self.results{end,layerIdx}.TN(nodeidx,classidx);
    end

    % Fix for async-only entry targets: use RN (response time per visit) for residt
    % The host layer closed model incorrectly splits residence time (WN) between
    % activities when an entry only receives async calls (no sync callers).
    % For async-only entries, use RN instead of WN since the async arrivals
    % don't share the closed chain's visit ratio - each async arrival gets
    % the full response time per visit.
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

    self.servtproc{aidx} = Exp.fitMean(self.servt(aidx));
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
            for w=1:(wnd_size-1)
                self.tput(aidx) = self.tput(aidx) + self.results{end-w,self.idxhash(idx)}.TN(nodeidx,classidx) / wnd_size;
            end
        else
            self.tput(aidx) = self.results{end,self.idxhash(idx)}.TN(nodeidx,classidx);
        end
        self.tputproc{aidx} = Exp.fitRate(self.tput(aidx));
    end
end

% TODO: obtain the join times
%self.joint = zeros(lqn.nidx,1);
%joinedacts = find(lqn.actpretype == ActivityPrecedenceType.PRE_AND)';

% obtain the call residence time
self.callservt = zeros(lqn.ncalls,1);
self.callresidt = zeros(lqn.ncalls,1);
for r=1:size(self.call_classes_updmap,1)
    idx = self.call_classes_updmap(r,1);
    cidx = self.call_classes_updmap(r,2);
    nodeidx = self.call_classes_updmap(r,3);
    classidx = self.call_classes_updmap(r,4);
    if self.call_classes_updmap(r,3) > 1
        if nodeidx == 1
            self.callservt(cidx) = 0;
        else
            self.callservt(cidx) = self.results{end, self.idxhash(idx)}.RN(nodeidx,classidx) * self.lqn.callproc{cidx}.getMean;
            self.callresidt(cidx) = self.results{end, self.idxhash(idx)}.WN(nodeidx,classidx);
        end
        % Apply under-relaxation to call service times
        omega = self.relax_omega;
        if omega < 1.0 && it > 1 && ~isnan(self.callservt_prev(cidx))
            self.callservt(cidx) = omega * self.callservt(cidx) + (1 - omega) * self.callservt_prev(cidx);
        end
        self.callservt_prev(cidx) = self.callservt(cidx);
    end
end

% then resolve the entry servt summing up these contributions
%entry_servt = zeros(lqn.nidx,1);
entry_servt = self.servtmatrix*[self.residt;self.callresidt(:)]; % Sum the residT of all the activities connected to this entry
entry_servt(1:lqn.eshift) = 0;


% Propagate forwarding calls: add target entry's service time to source entry
% When e0 forwards to e1 with probability p, callers of e0 see:
% e0's service + p * e1's service
% Process in topological order to handle forwarding chains correctly
for cidx = 1:lqn.ncalls
    if lqn.calltype(cidx) == CallType.FWD
        source_eidx = lqn.callpair(cidx, 1);
        target_eidx = lqn.callpair(cidx, 2);
        fwd_prob = lqn.callproc{cidx}.getMean();

        % Get target entry's service time
        target_servt = entry_servt(target_eidx);

        % If target entry doesn't have computed service time (forwarding-only target),
        % compute it directly from its activities' host demands
        if target_servt == 0 || isnan(target_servt)
            target_servt = 0;
            % Sum host demands of all activities bound to this entry
            for aidx = lqn.actsof{target_eidx}
                if ~isempty(lqn.hostdem{aidx})
                    target_servt = target_servt + lqn.hostdem{aidx}.getMean();
                end
            end
        end

        entry_servt(source_eidx) = entry_servt(source_eidx) + fwd_prob * target_servt;
    end
end

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
            %entry_servt_refstat = self.ensemble{self.idxhash(hidx)}.classes{tidxclass}.refstat;
            %entry_servt_z = entry_servt_refstat.serviceProcess{self.ensemble{self.idxhash(hidx)}.classes{tidxclass}.index}.getMean();
            %entry_servt(eidx) = self.ensemble{self.idxhash(hidx)}.classes{tidxclass}.population / entry_tput - entry_servt_z;
            self.servt(eidx) = entry_servt(eidx) * task_tput / max(GlobalConstants.Zero, entry_tput);
            self.residt(eidx) = entry_servt(eidx) * task_tput / max(GlobalConstants.Zero, entry_tput);
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
            % Get entry throughput (use task throughput as approximation if entry not available)
            if self.tput(eidx) > GlobalConstants.FineTol
                entry_tput = self.tput(eidx);
            elseif self.tput(tidx) > GlobalConstants.FineTol
                entry_tput = self.tput(tidx);
            else
                entry_tput = 0;
            end

            % Compute overtaking probability now that throughput is available
            if entry_tput > GlobalConstants.FineTol
                self.prOvertake(e) = self.overtake_prob(eidx);
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
        self.servtproc{eidx} = Exp.fitMean(self.servt(eidx));
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
            self.callservtproc{cidx} = Exp.fitMean(self.callservt(cidx));
        end
    end
end

self.ptaskcallers = zeros(size(self.ptaskcallers));
% determine ptaskcallers for direct callers to tasks
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if ~lqn.isref(tidx)
        [calling_idx, ~] = find(lqn.iscaller(:, lqn.entriesof{tidx})); %#ok<ASGLU>
        callers = intersect(lqn.tshift+(1:lqn.ntasks), unique(calling_idx)');
        caller_tput = zeros(1,lqn.ntasks);
        for caller_idx=callers(:)'
            caller_idxclass = self.ensemble{self.idxhash(tidx)}.attribute.tasks(1+find(self.ensemble{self.idxhash(tidx)}.attribute.tasks(2:end,2) == caller_idx),1);
            caller_tput(caller_idx-lqn.tshift)  = sum(self.results{end,self.idxhash(tidx)}.TN(self.ensemble{self.idxhash(tidx)}.attribute.clientIdx,caller_idxclass));
        end
        self.ptaskcallers(tidx,(lqn.tshift+1):(lqn.tshift+lqn.ntasks))= caller_tput / max(GlobalConstants.Zero, sum(caller_tput));
    end
end

% determine ptaskcallers for direct callers to hosts
for hidx = 1:lqn.nhosts
    if ~self.ignore(tidx) && ~self.ignore(hidx)
        caller_tput = zeros(1,lqn.ntasks);
        callers = lqn.tasksof{hidx};
        for caller_idx=callers
            caller_idxclass = self.ensemble{self.idxhash(hidx)}.attribute.tasks(find(self.ensemble{self.idxhash(hidx)}.attribute.tasks(:,2) == caller_idx),1);
            caller_tput(caller_idx-lqn.tshift)  = caller_tput(caller_idx-lqn.tshift) + sum(self.results{end,self.idxhash(hidx)}.TN(self.ensemble{self.idxhash(hidx)}.attribute.clientIdx,caller_idxclass));
        end
        self.ptaskcallers(hidx,(lqn.tshift+1):(lqn.tshift+lqn.ntasks)) = caller_tput / max(GlobalConstants.Zero, sum(caller_tput));
    end
end

% impute call probability using a DTMC random walk on the taskcaller graph
P = self.ptaskcallers;
P = dtmc_makestochastic(P); % hold mass at reference stations when there
self.ptaskcallers_step{1} = P; % configure step 1
for h = 1:lqn.nhosts
    hidx=h;
    for tidx = lqn.tasksof{hidx}
        % initialize the probability mass on tidx
        x0 = zeros(length(self.ptaskcallers),1);
        x0(hidx) = 1;
        x0=x0(:)';
        % start the walk backward to impute probability of indirect callers
        x = x0*P; % skip since pcallers already calculated in this case
        for step=2:self.nlayers % upper bound on maximum dag height
            x = x*P;
            % here self.ptaskcallers_step{step}(tidx,remidx) is the
            % probability that the request to tidx comes from remidx
            self.ptaskcallers_step{step}(tidx,:) = x(:);
            % here self.ptaskcallers_step{step}(hidx,remidx) is the
            % probability that the request to hidx comes from remidx
            % through the call that tidx puts on hidx, which is
            % weighted by the relative tput of tidx on the tasks running on
            % hidx
            self.ptaskcallers_step{step}(hidx,:) = self.ptaskcallers(hidx,tidx)*x(:)';
            if sum(x(find(lqn.isref)))>1.0-self.options.tol %#ok<FNDSB>
                % if all the probability mass has reached backwards the
                % reference stations, then stop
                break;
            end
            self.ptaskcallers(:,tidx) = max([self.ptaskcallers(:,tidx), x(:)],[],2);
        end
    end
end
self.ensemble = ensemble;
end
