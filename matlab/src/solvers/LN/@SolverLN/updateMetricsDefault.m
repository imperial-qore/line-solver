function updateMetricsDefault(self, it)
ensemble = self.ensemble;
lqn = self.lqn;

% obtain the activity service times
self.servt = zeros(lqn.nidx,1);
for r=1:size(self.servt_classes_updmap,1)
    idx = self.servt_classes_updmap(r,1);
    aidx = self.servt_classes_updmap(r,2);
    nodeidx = self.servt_classes_updmap(r,3);
    classidx = self.servt_classes_updmap(r,4);

    % store the residence times and tput at this layer to become
    % the servt / tputs of aidx in another layer, as needed
    iter_min = min(30,ceil(self.options.iter_max/4));
    wnd_size = (it-self.averagingstart+1);
    if ~isempty(self.averagingstart) && it>=iter_min % assume steady-state
        self.servt(aidx) = 0;
        self.tput(aidx) = 0;
        for w=1:(wnd_size-1)
            self.servt(aidx) = self.servt(aidx) + self.results{end-w,self.idxhash(idx)}.WN(nodeidx,classidx) / wnd_size;
            self.tput(aidx) = self.tput(aidx) + self.results{end-w,self.idxhash(idx)}.TN(nodeidx,classidx) / wnd_size;
        end
    else
        self.servt(aidx) = self.results{end,self.idxhash(idx)}.WN(nodeidx,classidx);
        self.tput(aidx) = self.results{end,self.idxhash(idx)}.TN(nodeidx,classidx);
    end
    self.servtproc{aidx} = Exp.fitMean(self.servt(aidx));
    self.tputproc{aidx} = Exp.fitRate(self.tput(aidx));
end

% obtain the call residence time
self.callresidt = zeros(lqn.ncalls,1);
for r=1:size(self.call_classes_updmap,1)
    idx = self.call_classes_updmap(r,1);
    cidx = self.call_classes_updmap(r,2);
    nodeidx = self.call_classes_updmap(r,3);
    classidx = self.call_classes_updmap(r,4);
    if self.call_classes_updmap(r,3) > 1
        if nodeidx == 1
            self.callresidt(cidx) = 0;
        else
            self.callresidt(cidx) = self.results{end, self.idxhash(idx)}.WN(nodeidx,classidx);
        end
    end
end

% then resolve the entry servt summing up these contributions
entry_servt = self.servtmatrix*[self.servt;self.callresidt]; % Sum the residT of all the activities connected to this entry
entry_servt(1:lqn.eshift) = 0;

% this block fixes the problem that ResidT is scaled so that the
% task as Vtask=1, but in call servt the entries need to have Ventry=1
for eidx=(lqn.eshift+1):(lqn.eshift+lqn.nentries)
    tidx = lqn.parent(eidx); % task of entry
    hidx = lqn.parent(tidx); %host of entry
    % eidxname = lqn.names{eidx};
    % get class in host layer of task and entry
    tidxclass = ensemble{self.idxhash(hidx)}.attribute.tasks(find(ensemble{self.idxhash(hidx)}.attribute.tasks(:,2) == tidx),1);
    eidxclass = ensemble{self.idxhash(hidx)}.attribute.entries(find(ensemble{self.idxhash(hidx)}.attribute.entries(:,2) == eidx),1);
    task_tput  = sum(self.results{end,self.idxhash(hidx)}.TN(ensemble{self.idxhash(hidx)}.attribute.clientIdx,tidxclass));
    entry_tput = sum(self.results{end,self.idxhash(hidx)}.TN(ensemble{self.idxhash(hidx)}.attribute.clientIdx,eidxclass));
    self.servt(eidx) = entry_servt(eidx) * task_tput / max(GlobalConstants.Zero, entry_tput);
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
            self.callresidt(cidx) = self.servt(eidx);
            self.callresidtproc{cidx} = self.servtproc{eidx};
        else
            % note that respt is per visit, so number of calls is 1
            self.callresidtproc{cidx} = Exp.fitMean(self.callresidt(cidx));
        end
    end
end

self.ptaskcallers = zeros(size(self.ptaskcallers));
% determine ptaskcallers for direct callers to tasks
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if ~lqn.isref(tidx)
        [calling_idx, called_entries] = find(lqn.iscaller(:, lqn.entriesof{tidx})); %#ok<ASGLU>
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
    caller_tput = zeros(1,lqn.ntasks);
    callers = lqn.tasksof{hidx};
    for caller_idx=callers
        caller_idxclass = self.ensemble{self.idxhash(hidx)}.attribute.tasks(find(self.ensemble{self.idxhash(hidx)}.attribute.tasks(:,2) == caller_idx),1);
        caller_tput(caller_idx-lqn.tshift)  = caller_tput(caller_idx-lqn.tshift) + sum(self.results{end,self.idxhash(hidx)}.TN(self.ensemble{self.idxhash(hidx)}.attribute.clientIdx,caller_idxclass));
    end
    self.ptaskcallers(hidx,(lqn.tshift+1):(lqn.tshift+lqn.ntasks)) = caller_tput / max(GlobalConstants.Zero, sum(caller_tput));
end

% impute call probability using a DTMC random walk on the taskcaller graph
P = self.ptaskcallers;
P = dtmc_makestochastic(P); % hold mass at reference stations when there
self.ptaskcallers_step{1} = P;
for h = 1:lqn.nhosts
    hidx=h;
    for tidx = lqn.tasksof{hidx}
        % initialize the probability mass on tidx
        x0 = zeros(length(self.ptaskcallers),1);
        x0(hidx) = 1;
        x0=x0(:)';
        step = 1;
        % start the walk backward to impute probability of indirect callers
        x = x0*P; % skip since pcallers already calculated in this case
        for e=1:self.nlayers % upper bound on maximum dag height
            step = step + 1;
            x = x*P;
            if sum(x(find(lqn.isref)))>1.0-self.options.tol %#ok<FNDSB>
                % if all the probability mass has reached backwards the
                % reference stations, then stop
                break;
            end
            self.ptaskcallers_step{step}(:,tidx) = x(:);
            self.ptaskcallers(:,tidx) = max([self.ptaskcallers(:,tidx), x(:)],[],2);
        end
    end
end
self.ensemble = ensemble;
end