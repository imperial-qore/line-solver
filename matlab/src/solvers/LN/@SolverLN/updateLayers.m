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
callresidtproc = self.callresidtproc;

% reassign service times
for r=1:size(thinkt_classes_updmap,1)
    if mod(it, 0)
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
    if aidx < lqn.tshift + lqn.ntasks
        % aidx here is actually set to tidx in buildLayersRecursive
        switch class.type
            case 'closed'
                class.population = self.njobs(aidx,idx);
        end
    end
    node = ensemble{idxhash(idx)}.nodes{nodeidx};
    switch nodeidx
        case  ensemble{idxhash(idx)}.attribute.clientIdx
            if lqn.type(aidx) == LayeredNetworkElement.TASK
                if lqn.schedid(aidx) ~= SchedStrategy.ID_REF
                    if ~isempty(thinktproc{aidx}) % this is empty for isolated components, which can be ignored
                        %                        if it==1
                        node.setService(class, thinktproc{aidx});
                        %                         else
                        %                             if ~node.serviceProcess{class.index}.isImmediate()
                        %                                 node.serviceProcess{class.index}.updateMean(thinktproc{aidx}.getMean);
                        %                                 node.server.serviceProcess{class.index}{end}.updateMean(thinktproc{aidx}.getMean);
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
                    %                             node.serviceProcess{class.index}.updateMean(servtproc{aidx}.getMean);
                    %                             node.server.serviceProcess{class.index}{end}.updateMean(servtproc{aidx}.getMean);
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
                %                         node.serviceProcess{class.index}.updateMean(servtproc{aidx}.getMean);
                %                         node.server.serviceProcess{class.index}{end}.updateMean(servtproc{aidx}.getMean);
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
    if mod(it, 0)
        ri = size(arvproc_classes_updmap,1) - r + 1;
    else
        ri = r;
    end
    idx = arvproc_classes_updmap(ri,1);
    cidx = arvproc_classes_updmap(ri,2);
    nodeidx = arvproc_classes_updmap(ri,3);
    classidx = arvproc_classes_updmap(ri,4);
    class = ensemble{idxhash(arvproc_classes_updmap(ri,1))}.classes{classidx};
    node = ensemble{idxhash(idx)}.nodes{nodeidx};
    node.setArrival(class, tputproc{lqn.callpair(cidx,1)});
end

% reassign call service time / response time
for c=1:size(call_classes_updmap,1)
    if mod(it, 0)
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
            node.setService(class, callresidtproc{cidx});
        case ensemble{idxhash(idx)}.attribute.serverIdx % the call is processed by the server, then replace with the svc time
            eidx = lqn.callpair(cidx,2);
            %tidx = lqn.parent(eidx);
            %eidxclass = self.ensemble{self.idxhash(tidx)}.attribute.calls(find(self.ensemble{self.idxhash(tidx)}.attribute.calls(:,4) == eidx),1);
            %eidxchain = find(self.ensemble{self.idxhash(tidx)}.getStruct.chains(:,eidxclass)>0);
            %qn = self.ensemble{self.idxhash(tidx)}.getStruct;
            %servtproc{eidx}.updateMean(servtproc{eidx}.getMean * qn.visits{eidxchain}(1,qn.refclass(eidxchain)) / qn.visits{eidxchain}(2,eidxclass))
            %%task_tput = sum(self.results{end,self.idxhash(tidx)}.TN(self.ensemble{self.idxhash(tidx)}.attribute.serverIdx,eidxclass))
            %%entry_tput = sum(self.results{end,self.idxhash(tidx)}.TN(self.ensemble{self.idxhash(tidx)}.attribute.serverIdx,eidxclass))
            %            if it==1
            
            node.setService(class, servtproc{eidx});
            %             else
            %                 if ~node.serviceProcess{class.index}.isImmediate()
            %                     node.serviceProcess{class.index}.updateMean(servtproc{eidx}.getMean);
            %                     node.server.serviceProcess{class.index}{end}.updateMean(servtproc{eidx}.getMean);
            %                 else
            %                     node.setService(class, servtproc{eidx});
            %                 end
            %             end

    end
end
end