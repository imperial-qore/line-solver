function updateLayers(self, it)
lqn = self.lqn;
ensemble = self.ensemble;
idxhash = self.idxhash;
thinktproc = self.thinktproc;
servtproc = self.servtproc;
thinkt_classes_updmap = self.thinkt_classes_updmap;
arvproc_classes_updmap = self.arvproc_classes_updmap;
tputproc = self.tputproc;
call_classes_updmap = self.call_classes_updmap;
callservtproc = self.callservtproc;

% reassign service times
for r=1:size(thinkt_classes_updmap,1)
    if mod(it, 2)==1 % elevator
        ri = size(thinkt_classes_updmap,1) - r + 1;
    else
        ri = r;
    end
    idx = thinkt_classes_updmap(ri,1);
    aidx = thinkt_classes_updmap(ri,2);
    nodeidx = thinkt_classes_updmap(ri,3);
    classidx = thinkt_classes_updmap(ri,4);
    class = ensemble{idxhash(thinkt_classes_updmap(ri,1))}.classes{classidx};
    % here update the number of jobs in the task chain
    if aidx <= lqn.tshift + lqn.ntasks
        % aidx here is actually set to tidx in buildLayersRecursive
        switch class.type
            case JobClassType.CLOSED
                if self.options.config.interlocking
                    class.population = self.njobs(aidx,idx);
                end
        end
    end
    node = ensemble{idxhash(idx)}.nodes{nodeidx};
    switch nodeidx
        case  ensemble{idxhash(idx)}.attribute.clientIdx
            if lqn.type(aidx) == LayeredNetworkElement.TASK
                if lqn.sched(aidx) ~= SchedStrategy.REF
                    if ~isempty(thinktproc{aidx}) % this is empty for isolated components, which can be ignored
                        %                        if it==1
                        node.setService(class, thinktproc{aidx});
                        %                         else
                        %                             if ~node.serviceProcess{class.index}.isImmediate()
                        %                                 node.serviceProcess{class.index}.setMean(thinktproc{aidx}.getMean);
                        %                                 node.server.serviceProcess{class.index}{end}.setMean(thinktproc{aidx}.getMean);
                        %                             else
                        %                                 node.setService(class, thinktproc{aidx});
                        %                             end
                        %                         end
                    end
                else
                    %                    if it==1
                    node.setService(class, servtproc{aidx});
                    %                     else
                    %                         if ~node.serviceProcess{class.index}.isImmediate()
                    %                             node.serviceProcess{class.index}.setMean(servtproc{aidx}.getMean);
                    %                             node.server.serviceProcess{class.index}{end}.setMean(servtproc{aidx}.getMean);
                    %                         else
                    %                             node.setService(class, servtproc{aidx});
                    %                         end
                    %                     end
                end
            else
                %                if it==1
                node.setService(class, servtproc{aidx});
                %                 else
                %                     if ~node.serviceProcess{class.index}.isImmediate()
                %                         node.serviceProcess{class.index}.setMean(servtproc{aidx}.getMean);
                %                         node.server.serviceProcess{class.index}{end}.setMean(servtproc{aidx}.getMean);
                %                     else
                %                         node.setService(class, servtproc{aidx});
                %                     end
                %                 end
            end
        case ensemble{idxhash(idx)}.attribute.serverIdx
            node.setService(class, servtproc{aidx});
    end
end

% reassign arrival rates
for r=1:size(arvproc_classes_updmap,1)
    %for r=1:0
    if  mod(it, 2)==1 % elevator
        ri = size(arvproc_classes_updmap,1) - r + 1;
    else
        ri = r;
    end
    idx = arvproc_classes_updmap(ri,1);
    eidx_or_cidx = arvproc_classes_updmap(ri,2);
    nodeidx = arvproc_classes_updmap(ri,3);
    classidx = arvproc_classes_updmap(ri,4);
    class = ensemble{idxhash(arvproc_classes_updmap(ri,1))}.classes{classidx};
    node = ensemble{idxhash(idx)}.nodes{nodeidx};

    if eidx_or_cidx < 0  % Entry-level arrival (negative index)
        eidx = -eidx_or_cidx;
        node.setArrival(class, lqn.arrival{eidx});
    else  % Async call arrival (positive index)
        cidx = eidx_or_cidx;
        node.setArrival(class, tputproc{lqn.callpair(cidx,1)});
    end
end

% reassign call service time / response time
for c=1:size(call_classes_updmap,1)
    if  mod(it, 2)==1 % elevator
        ci = size(call_classes_updmap,1) - c + 1;
    else
        ci = c;
    end
    idx = call_classes_updmap(ci,1);
    cidx = call_classes_updmap(ci,2);
    nodeidx = call_classes_updmap(ci,3);
    class = ensemble{idxhash(call_classes_updmap(ci,1))}.classes{call_classes_updmap(ci,4)};
    node = ensemble{idxhash(idx)}.nodes{nodeidx};
    switch nodeidx
        case ensemble{idxhash(idx)}.attribute.clientIdx % client
            node.setService(class, callservtproc{cidx});
        case ensemble{idxhash(idx)}.attribute.serverIdx % the call is processed by the server, then replace with the svc time
            eidx = lqn.callpair(cidx,2);
            %tidx = lqn.parent(eidx);
            %eidxclass = self.ensemble{self.idxhash(tidx)}.attribute.calls(find(self.ensemble{self.idxhash(tidx)}.attribute.calls(:,4) == eidx),1);
            %eidxchain = find(self.ensemble{self.idxhash(tidx)}.getStruct.chains(:,eidxclass)>0);
            %qn = self.ensemble{self.idxhash(tidx)}.getStruct;
            %servtproc{eidx}.setMean(servtproc{eidx}.getMean * qn.visits{eidxchain}(1,qn.refclass(eidxchain)) / qn.visits{eidxchain}(2,eidxclass))
            %%task_tput = sum(self.results{end,self.idxhash(tidx)}.TN(self.ensemble{self.idxhash(tidx)}.attribute.serverIdx,eidxclass))
            %%entry_tput = sum(self.results{end,self.idxhash(tidx)}.TN(self.ensemble{self.idxhash(tidx)}.attribute.serverIdx,eidxclass))
            %            if it==1

            node.setService(class, servtproc{eidx});
            %             else
            %                 if ~node.serviceProcess{class.index}.isImmediate()
            %                     node.serviceProcess{class.index}.setMean(servtproc{eidx}.getMean);
            %                     node.server.serviceProcess{class.index}{end}.setMean(servtproc{eidx}.getMean);
            %                 else
            %                     node.setService(class, servtproc{eidx});
            %                 end
            %             end

    end
end

% reassign activity think times
actthinkt_classes_updmap = self.actthinkt_classes_updmap;
for r=1:size(actthinkt_classes_updmap,1)
    if mod(it, 2)==1 % elevator
        ri = size(actthinkt_classes_updmap,1) - r + 1;
    else
        ri = r;
    end
    idx = actthinkt_classes_updmap(ri,1);
    aidx = actthinkt_classes_updmap(ri,2);
    nodeidx = actthinkt_classes_updmap(ri,3);
    classidx = actthinkt_classes_updmap(ri,4);
    class = ensemble{idxhash(actthinkt_classes_updmap(ri,1))}.classes{classidx};
    node = ensemble{idxhash(idx)}.nodes{nodeidx};
    % Activity think-time is at client delay node
    if ~isempty(lqn.actthink{aidx})
        node.setService(class, lqn.actthink{aidx});
    end
end
end
