function sync = refreshSync(self)
% SYNC = REFRESHSYNC()

sn = self.sn;
local = self.getNumberOfNodes+1;
nclasses = sn.nclasses;
sync = {};
emptystate = cellzeros(sn.nnodes,1,0,0);
if any(sn.isstatedep(:))
    rtmask = self.sn.rtfun(emptystate, emptystate);
else
    rtmask = ceil(self.sn.rt);
end

for ind=1:sn.nnodes
    for r=1:sn.nclasses
        if sn.isstation(ind) && sn.phases(sn.nodeToStation(ind),r)> 1
            % Phase-change action
            sync{end+1,1} = struct('active',cell(1),'passive',cell(1));
            sync{end,1}.active{1} = Event(EventType.PHASE, ind, r);
            sync{end,1}.passive{1} = Event(EventType.LOCAL, local, r, 1.0);
        end
        if sn.isstation(ind) && isfield(sn,'impatienceClass') && ~isempty(sn.impatienceClass) ...
                && sn.impatienceClass(sn.nodeToStation(ind),r) == ImpatienceType.RENEGING ...
                && sn.impatienceType(sn.nodeToStation(ind),r) == ProcessType.EXP
            % Reneging action (exponential patience): a waiting class-r job
            % abandons the queue at a memoryless per-job rate; the job leaves
            % the system (passive LOCAL, mirroring the phase-change action).
            sync{end+1,1} = struct('active',cell(1),'passive',cell(1));
            sync{end,1}.active{1} = Event(EventType.RENEGE, ind, r);
            sync{end,1}.passive{1} = Event(EventType.LOCAL, local, r, 1.0);
        end
        if sn.isstation(ind) && isfield(sn,'retrialProc') && ~isempty(sn.retrialProc) ...
                && ~isempty(sn.retrialProc{sn.nodeToStation(ind),r}) ...
                && sn.retrialType(sn.nodeToStation(ind),r) == ProcessType.EXP
            % see _kb/04-networkstruct.md (refreshGlobalSync/refreshSync) for rationale
            sync{end+1,1} = struct('active',cell(1),'passive',cell(1));
            sync{end,1}.active{1} = Event(EventType.RETRY, ind, r);
            sync{end,1}.passive{1} = Event(EventType.LOCAL, local, r, 1.0);
        end
        if sn.isstation(ind) && r == 1 && isfield(sn,'hasbreakdown') && ~isempty(sn.hasbreakdown) ...
                && numel(sn.hasbreakdown) >= ind && sn.hasbreakdown(ind) == 1
            % Server failure and repair actions. Both are properties of the
            % SERVER, not of a class, so exactly one pair is emitted per station
            % (guarded on r == 1) rather than one pair per class. Like the
            % phase-change action they move no job, so the passive half is LOCAL.
            sync{end+1,1} = struct('active',cell(1),'passive',cell(1));
            sync{end,1}.active{1} = Event(EventType.FAILURE, ind, r);
            sync{end,1}.passive{1} = Event(EventType.LOCAL, local, r, 1.0);
            sync{end+1,1} = struct('active',cell(1),'passive',cell(1));
            sync{end,1}.active{1} = Event(EventType.REPAIR, ind, r);
            sync{end,1}.passive{1} = Event(EventType.LOCAL, local, r, 1.0);
        end
        if sn.isstation(ind) && sn.sched(sn.nodeToStation(ind)) == SchedStrategy.POLLING
            % see _kb/04-networkstruct.md (refreshGlobalSync/refreshSync) for rationale
            pinfoSync = State.pollingInfo(sn, ind);
            if ~isempty(pinfoSync) && pinfoSync.hasSw(r)
                sync{end+1,1} = struct('active',cell(1),'passive',cell(1));
                sync{end,1}.active{1} = Event(EventType.SWITCH, ind, r);
                sync{end,1}.passive{1} = Event(EventType.LOCAL, local, r, 1.0);
            end
        end
        if sn.isstateful(ind)
            if sn.nodetype(ind) == NodeType.Cache
                if ~isnan(sn.nodeparam{ind}.pread{r}) % class can read
                    sync{end+1,1}.active{1} = Event(EventType.READ, ind, r);
                    sync{end,1}.passive{1} = Event(EventType.READ, local, r, 1.0);
                end
            elseif sn.nodetype(ind) == NodeType.Transition
                for m=1:sn.nodeparam{ind}.nmodes
                    % server phase change
                    sync{end+1,1}.active{1} = Event(EventType.PHASE, ind, m);
                    sync{end,1}.passive{1} = Event(EventType.LOCAL, local, m, 1.0);
                end
            end
            if sn.nodetype(ind) == NodeType.Fork
                % Stateful Fork nodes (FJ-augmented copies) do not emit
                % departure syncs: the atomic multi-branch emission is
                % handled by the fork firing synchronizations (sn.fjsync)
                continue
            end
            isf = sn.nodeToStateful(ind);
            for jnd=1:sn.nnodes
                if sn.isstateful(jnd)
                    jsf = sn.nodeToStateful(jnd);
                    for s=1:nclasses
                        p = rtmask((isf-1)*nclasses+r,(jsf-1)*nclasses+s);
                        if p > 0
                            new_sync = struct('active',cell(1),'passive',cell(1));
                            new_sync.active{1} = Event(EventType.DEP, ind, r);
                            switch sn.routing(ind,s)
                                case {RoutingStrategy.RROBIN, RoutingStrategy.WRROBIN, RoutingStrategy.JSQ, RoutingStrategy.SQ, RoutingStrategy.SDR}
                                    new_sync.passive{1} = Event(EventType.ARV, jnd, s, @(state_before, state_after) at(self.sn.rtfun(state_before, state_after), (isf-1)*nclasses+r, (jsf-1)*nclasses+s));
                                otherwise
                                    new_sync.passive{1} = Event(EventType.ARV, jnd, s, sn.rt((isf-1)*nclasses+r, (jsf-1)*nclasses+s));
                            end
                            sync{end+1,1} = new_sync;
                        end
                    end
                end
            end
        end
    end
end
if ~isempty(self.sn) %&& isprop(self.sn,'nvars')
    self.sn.sync = sync;
end
end
