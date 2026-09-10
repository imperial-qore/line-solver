function updateThinkTimes(self, it)
% Update the think times of all callers at iteration it. The method handles
% differently the case where a caller is a ref task than the case where the
% caller is a queueing station. A coarse heuristic is used when one or more
% callers are themselves infinite servers.

% Under a PH encoding a caller reaches the server once per invocation, so the
% station rate is not the task's invocation rate -- see updateThinkTimesPH
if self.isPHEncoding()
    updateThinkTimesPH(self, it);
    return
end

% create local variable due to MATLAB's slow access to self properties
lqn = self.lqn;
idxhash = self.idxhash;
results = self.results;

% main code starts here
if size(lqn.iscaller,2) > 0 % ignore models without callers
    torder = 1:(lqn.ntasks); % set sequential order to update the tasks
    % solve all task models
    for t = torder
        tidx = lqn.tshift + t;
        % Only a REFERENCE task's think time separates one request from the
        % next. Charging it on a served task made it a per-request delay and
        % throttled the task: lqn_basic's T3, 25 threads and a declared think
        % time of 4, was capped at 25/(4+0.02)=6.219 completions per second
        % where lqsim reports 66.7, LDES 66.955 and lqns 75.6. See
        % lqn_ref_thinktime and _kb/06-solver-catalog.md (LN section).
        tidx_thinktime = lqn_ref_thinktime(lqn, tidx);
        %if ~lqn.isref(tidx) && ~isnan(idxhash(tidx)) % update tasks ignore ref tasks and empty tasks
        if ~isnan(self.idxhash(tidx)) % this skips all REF tasks
            % obtain total self.tput of task t
            % mean throughput of task t in the model where it is a server, summed across replicas
            njobs = max(self.njobs(tidx,:)); % we use njobs to adapt to interlocking corrections
            [layer_t, station_t] = self.layerOf(tidx);
            self.tput(tidx) = lqn.repl(tidx)*sum(results{end,layer_t}.TN(station_t,:),2);
            % The closure below is a CYCLE of ONE thread, so both sides of it
            % must be per replica. `tput` is the element's total and `station_t`
            % is one representative station, so the rate is divided back down
            % and the thread pool is taken per replica -- `njobs` is the total
            % in whichever layer used that convention for this caller. Mixing
            % them halved the think time of the two-replica probe and reported
            % its task at 5.86443 where LDES measures 8.30737. Inert at repl 1.
            nrep = max(1, lqn.repl(tidx));
            xcyc = self.tput(tidx) / nrep;
            if nrep > 1 && isfinite(lqn.maxmult(tidx)) && lqn.maxmult(tidx) > 0
                njobs = lqn.maxmult(tidx);
            end
            if lqn.sched(tidx) == SchedStrategy.INF % first we consider the update where t is an infinite server
                % obtain total self.utilization of task t
                self.util(tidx) = sum(results{end,layer_t}.UN(station_t,:),2);
                % key LQN think-time update: in LINE an infinite-server self.utilization is dimensionally a mean number of jobs
                self.thinkt(tidx) = max(GlobalConstants.Zero, (njobs-self.util(tidx)) / xcyc - tidx_thinktime);
            else % otherwise we consider the case where t is a regular queueing station (other than an infinite server)
                self.util(tidx) = sum(results{end,layer_t}.UN(station_t,:),2); % self.utilization of t as a server
                % key LQN think-time update: in LINE self.utilization is scaled to [0,1] for all queueing stations regardless of server count
                self.thinkt(tidx) = max(GlobalConstants.Zero, njobs*abs(1-self.util(tidx)) / xcyc - tidx_thinktime);
            end
            % Phase-2 tail: the entry replies after phase 1, so the station
            % serves the caller for residt, but the thread stays busy for
            % servt. What is left is busy time, not think time -- without this
            % the task cycles slower than its callers call it, and throughput
            % is not conserved across the call.
            if ~isFlatLayering(self)
                self.thinkt(tidx) = max(GlobalConstants.Zero, self.thinkt(tidx) - phase2Tail(self, lqn, tidx));
            end
            % The cold start goes the OTHER way from the phase-2 tail. A caller
            % class cycles as delay + station service, and the station serves only
            % the host demand: the charge is on the entry, not on any activity's
            % demand, so the station never sees it and the delay has to carry it.
            % Without this the callee layer cycled faster than its callers drive
            % it, 0.529412 against 0.5 on lqn_setup, with the caller conserved and
            % the callee not. Zero for every task without a setup.
            self.thinkt(tidx) = max(GlobalConstants.Zero, ...
                self.thinkt(tidx) + lqn_setup_charge(self, lqn, tidx));
            % Recover from Inf/NaN: snap back to previous iteration's value
            if it > 1 && ~isnan(self.thinkt_prev(tidx))
                if isinf(self.thinkt(tidx)) || isnan(self.thinkt(tidx))
                    self.thinkt(tidx) = self.thinkt_prev(tidx);
                end
            end
            % Apply under-relaxation to think time if enabled
            omega = self.relax_omega;
            if omega < 1.0 && it > 1 && ~isnan(self.thinkt_prev(tidx))
                rawT = self.thinkt(tidx);
                prevT = self.thinkt_prev(tidx);
                % If recovering from crash (prev much larger than raw), snap to raw
                if prevT > 10 * rawT && rawT > GlobalConstants.FineTol
                    self.thinkt_prev(tidx) = rawT; % reset prev to allow recovery
                end
                self.thinkt(tidx) = omega * self.thinkt(tidx) + (1 - omega) * self.thinkt_prev(tidx);
            end
            self.thinkt_prev(tidx) = self.thinkt(tidx);
            self.thinktproc{tidx} = Exp.fitMean(self.thinkt(tidx) + tidx_thinktime);
        else % ref task, forwarding target or open-arrival target (no task layer)
            % A task reached only by an entry arrival has no caller, so no task
            % layer, so nothing above would ever set its surrogate delay and the
            % thread pool on its host layer would cycle against an Immediate one.
            % The stream is known, so the cycle is closed on it directly: the same
            % construction as a forwarding target below. See lqn_open_arrival.
            % Only when the arrival is the sole way in: with a caller or a
            % forwarding source the host layer still carries the open stream as a
            % class of its own, and closing the chain on the rate as well would
            % load the layer twice. The predicate matches buildLayersRecursive.
            arvrate = 0;
            if ~lqn.isref(tidx) && isfield(lqn,'arrival') && ~isempty(lqn.arrival) && iscell(lqn.arrival) ...
                    && ~any(any(full(lqn.issynccaller(:, lqn.entriesof{tidx})))) ...
                    && ~any(any(full(lqn.isasynccaller(:, lqn.entriesof{tidx})))) ...
                    && ~isFwdTargetTask(lqn, tidx)
                for eidx_arv = lqn.entriesof{tidx}
                    if eidx_arv <= numel(lqn.arrival) && ~isempty(lqn.arrival{eidx_arv})
                        m = lqn.arrival{eidx_arv}.getMean();
                        if isfinite(m) && m > GlobalConstants.FineTol
                            arvrate = arvrate + 1/m;
                        end
                    end
                end
            end
            if arvrate > GlobalConstants.FineTol
                njobs = max(self.njobs(tidx,:));
                if ~(njobs > 0)
                    njobs = lqn.maxmult(tidx);
                end
                tidx_thinktime = lqn_ref_thinktime(lqn, tidx);
                self.tput(tidx) = arvrate;
                host_residt = 0;
                for eidx_arv = lqn.entriesof{tidx}
                    for aidx_arv = lqn.actsof{eidx_arv}
                        if ~isnan(self.residt(aidx_arv))
                            host_residt = host_residt + self.residt(aidx_arv);
                        end
                    end
                end
                self.thinkt(tidx) = max(GlobalConstants.Zero, njobs/arvrate - host_residt - tidx_thinktime);
                omega = self.relax_omega;
                if omega < 1.0 && it > 1 && ~isnan(self.thinkt_prev(tidx))
                    self.thinkt(tidx) = omega * self.thinkt(tidx) + (1 - omega) * self.thinkt_prev(tidx);
                end
                self.thinkt_prev(tidx) = self.thinkt(tidx);
                self.thinktproc{tidx} = Exp.fitMean(self.thinkt(tidx) + tidx_thinktime);
                continue
            end
            % Check if this is a forwarding target task
            isFwdTarget = false;
            if ~lqn.isref(tidx)
                for eidx_fwd = lqn.entriesof{tidx}
                    for cidx_fwd = 1:lqn.ncalls
                        if lqn.calltype(cidx_fwd) == CallType.FWD && lqn.callpair(cidx_fwd, 2) == eidx_fwd
                            isFwdTarget = true;
                            source_eidx = lqn.callpair(cidx_fwd, 1);
                            source_tidx = lqn.parent(source_eidx);
                            fwd_prob = lqn.callproc{cidx_fwd}.getMean();
                            break;
                        end
                    end
                    if isFwdTarget; break; end
                end
            end
            if isFwdTarget
                % see _kb/06-solver-catalog.md (LN section) for rationale
                njobs = max(self.njobs(tidx,:));
                % a forwarding target is never a reference task, so its declared
                % think time does not enter the cycle either
                tidx_thinktime = lqn_ref_thinktime(lqn, tidx);
                arrival_rate = self.tput(source_tidx) * fwd_prob;
                if arrival_rate > GlobalConstants.FineTol && njobs > 0
                    self.tput(tidx) = arrival_rate;
                    % Subtract the processor response time for the target's
                    % activities (already computed by updateMetricsDefault)
                    target_eidx = lqn.callpair(cidx_fwd, 2);
                    host_residt = 0;
                    for aidx_fwd = lqn.actsof{target_eidx}
                        host_residt = host_residt + self.residt(aidx_fwd);
                    end
                    self.thinkt(tidx) = max(GlobalConstants.Zero, njobs / arrival_rate - host_residt - tidx_thinktime);
                else
                    % Source throughput not yet available; use large think time
                    self.thinkt(tidx) = 1000;
                end
                % Apply under-relaxation
                omega = self.relax_omega;
                if omega < 1.0 && it > 1 && ~isnan(self.thinkt_prev(tidx))
                    self.thinkt(tidx) = omega * self.thinkt(tidx) + (1 - omega) * self.thinkt_prev(tidx);
                end
                self.thinkt_prev(tidx) = self.thinkt(tidx);
                self.thinktproc{tidx} = Exp.fitMean(self.thinkt(tidx) + tidx_thinktime);
            else
                self.thinkt(tidx) = GlobalConstants.FineTol;
                self.thinktproc{tidx} = Immediate();
            end
        end
    end
end
end

function tail = phase2Tail(self, lqn, tidx)
% Mean time per request that a thread of task TIDX stays busy after replying,
% i.e. servt - residt averaged over the entries by their share of the requests
tail = 0;
if ~self.hasPhase2
    return
end
wtot = 0;
for eidx = lqn.entriesof{tidx}
    if self.servt_ph2(eidx) <= GlobalConstants.FineTol
        continue
    end
    w = self.tput(eidx);
    if w <= GlobalConstants.FineTol
        w = 1;
    end
    tail = tail + w * max(0, self.servt(eidx) - self.residt(eidx));
    wtot = wtot + w;
end
if wtot > GlobalConstants.FineTol
    tail = tail / wtot;
end
end

function tf = isFwdTargetTask(lqn, tidx)
% True when some entry of task TIDX is the target of a forwarding call
tf = false;
if lqn.ncalls == 0
    return
end
ents = lqn.entriesof{tidx};
for cidx = 1:lqn.ncalls
    if full(lqn.calltype(cidx)) == CallType.FWD && any(full(lqn.callpair(cidx,2)) == ents)
        tf = true;
        return
    end
end
end

function tf = isFlatLayering(self)
% True when the ensemble is a single flat layer
tf = isfield(self.options.config,'layering') && ...
    any(strcmpi(self.options.config.layering, {'flat','squashed'}));
end
