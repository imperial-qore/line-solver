function updateRoutingProbabilities(self, it)
for u = 1:length(self.unique_route_prob_updmap)
    if mod(it,0) % do elevator
        idx = self.unique_route_prob_updmap(u);
    else
        idx = self.unique_route_prob_updmap(length(self.unique_route_prob_updmap)-u+1);
    end
    idx_updated = false;
    P = self.ensemble{self.idxhash(idx)}.getLinkedRoutingMatrix;
    for r = find(self.route_prob_updmap(:,1) == idx)'
        host = self.route_prob_updmap(r,1);
        tidx_caller = self.route_prob_updmap(r,2);
        eidx = self.route_prob_updmap(r,3);
        nodefrom = self.route_prob_updmap(r,4);
        nodeto = self.route_prob_updmap(r,5);
        classidxfrom = self.route_prob_updmap(r,6);
        classidxto = self.route_prob_updmap(r,7);
        if ~isempty(self.ensemble{self.idxhash(idx)}.items) % if idx is a cache node
            Xtot = sum(self.results{end,self.idxhash(host)}.TN(self.ensemble{self.idxhash(host)}.attribute.serverIdx,:));
            if Xtot > 0
                hm_tput = sum(self.results{end,self.idxhash(host)}.TN(self.ensemble{self.idxhash(host)}.attribute.serverIdx,classidxto));
                P{classidxfrom,classidxto}(nodefrom, nodeto) = hm_tput / Xtot;
                idx_updated = true;
            end
        else % if idx is not a cache
            Xtot = sum(self.results{end,self.idxhash(tidx_caller)}.TN(self.ensemble{self.idxhash(tidx_caller)}.attribute.serverIdx,:));
            if Xtot > 0
                eidxclass = self.ensemble{self.idxhash(tidx_caller)}.attribute.calls(find(self.ensemble{self.idxhash(tidx_caller)}.attribute.calls(:,4) == eidx),1); %#ok<FNDSB>
                entry_tput = sum(self.results{end,self.idxhash(tidx_caller)}.TN(self.ensemble{self.idxhash(tidx_caller)}.attribute.serverIdx,eidxclass));
                P{classidxfrom,classidxto}(nodefrom, nodeto) = entry_tput / Xtot;
                idx_updated = true;
            end
        end
    end
    if idx_updated
        self.ensemble{self.idxhash(idx)}.relink(P);
    end
end
end

