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
                                case {RoutingStrategy.RROBIN, RoutingStrategy.WRROBIN, RoutingStrategy.JSQ, RoutingStrategy.RL}
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
