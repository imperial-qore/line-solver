function [outspace, outrate, outprob, eventCache] = afterEventCache(sn, ind, event, class, isSimulation, eventCache, R, space_buf, space_srv, space_var, key)
% job arrives in class, then reads and moves into hit or miss
% class, then departs

outspace = [];
outrate = [];
outprob = 1;
switch event
    case EventType.ARV
        space_srv(:,class) = space_srv(:,class) + 1;
        outspace = [space_srv, space_var]; % buf is empty
        outrate = -1*ones(size(outspace,1)); % passive action, rate is unspecified
    case EventType.DEP
        if space_srv(class)>0
            space_srv(:,class) = space_srv(:,class) - 1;
            switch sn.routing(ind,class)
                case RoutingStrategy.RROBIN
                    idx = find(space_var(sum(sn.nvars(ind,1:(R+class)))) == sn.nodeparam{ind}{class}.outlinks);
                    if idx < length(sn.nodeparam{ind}{class}.outlinks)
                        space_var(sum(sn.nvars(ind,1:(R+class)))) = sn.nodeparam{ind}{class}.outlinks(idx+1);
                    else
                        space_var(sum(sn.nvars(ind,1:(R+class)))) = sn.nodeparam{ind}{class}.outlinks(1);
                    end
            end
            outspace = [space_srv, space_var]; % buf is empty
            outrate = GlobalConstants.Immediate*ones(size(outspace,1)); % immediate action
        end
    case EventType.READ
        n = sn.nodeparam{ind}.nitems; % n items
        m = sn.nodeparam{ind}.itemcap; % capacity
        ac = sn.nodeparam{ind}.accost; % access cost
        hitclass = sn.nodeparam{ind}.hitclass;
        missclass = sn.nodeparam{ind}.missclass;
        h = length(m);
        replacement_id = sn.nodeparam{ind}.replacestrat;
        if space_srv(class)>0 && sum(space_srv)==1 %  a job of class is in
            p = sn.nodeparam{ind}.pread{class};
            en =  space_srv(:,class) > 0;
            space_srv_k = [];
            space_var_k = [];
            outrate = [];
            if any(en)
                for e=find(en)'
                    if isSimulation
                        % pick one item
                        kset = 1 + max([0,find( rand > cumsum(p) )]);
                        % pick one entry list for cache miss
                        % do not move this entry
                        l = 1 + max([0,find( rand > cumsum(ac{class,kset}(1,:)) )]);
                    else
                        kset = 1:n;
                    end
                    for k=kset % request to item k
                        space_srv_e = space_srv(e,:);
                        space_srv_e(class) = space_srv_e(class) - 1;
                        var = space_var(e,:);
                        posk = find(k==var,1,'first');

                        if isempty(posk) % CACHE MISS, can enter any list based on accessCost
                            space_srv_e(missclass(class)) = space_srv_e(missclass(class)) + 1;
                            switch replacement_id
                                case {ReplacementStrategy.FIFO, ReplacementStrategy.LRU, ReplacementStrategy.SFIFO}
                                    if isSimulation
                                        listidx = l - 1; % l is accessCost column index, listidx is actual list (1-indexed)
                                        if listidx > 0 % only cache if listidx is valid (l >= 2)
                                            varp = var;
                                            varp(cpos(listidx,2):cpos(listidx,m(listidx))) = var(cpos(listidx,1):cpos(listidx,m(listidx)-1));
                                            varp(cpos(listidx,1)) = k; % head of list listidx
                                            space_srv_k = [space_srv_k; space_srv_e];
                                            space_var_k = [space_var_k; varp];
                                            %% no p(k) weighting since that goes in the outprob vec
                                            outrate(end+1,1) = GlobalConstants.Immediate;
                                            outprob(end+1,1) = ac{class,k}(1,l) * p(k);
                                        else
                                            % Cache reject (l=1): pass through without caching
                                            space_srv_k = [space_srv_k; space_srv_e];
                                            space_var_k = [space_var_k; var];
                                            outrate(end+1,1) = GlobalConstants.Immediate;
                                            outprob(end+1,1) = ac{class,k}(1,l) * p(k);
                                        end
                                    else
                                        for l=2:(h+1) % iterate over all possible target lists (columns 2 to h+1)
                                            listidx = l - 1; % l is accessCost column index, listidx is actual list (1-indexed)
                                            varp = var;
                                            varp(cpos(listidx,2):cpos(listidx,m(listidx))) = var(cpos(listidx,1):cpos(listidx,m(listidx)-1));
                                            varp(cpos(listidx,1)) = k; % head of list listidx
                                            space_srv_k = [space_srv_k; space_srv_e];
                                            space_var_k = [space_var_k; varp];
                                            outrate(end+1,1) = ac{class,k}(1,l) * p(k) * GlobalConstants.Immediate;
                                            outprob(end+1,1) = 1;
                                        end
                                    end
                                case ReplacementStrategy.RR
                                    if isSimulation
                                        listidx = l - 1; % l is accessCost column index, listidx is actual list (1-indexed)
                                        if listidx > 0 % only cache if listidx is valid (l >= 2)
                                            varp = var; % var'
                                            r = randi(m(listidx),1,1);
                                            varp(cpos(listidx,r)) = k;
                                            space_srv_k = [space_srv_k; space_srv_e];
                                            space_var_k = [space_var_k; (varp)];
                                            outrate(end+1,1) = GlobalConstants.Immediate;
                                            outprob(end+1,1) = ac{class,k}(1,l) * p(k);
                                        else
                                            % Cache reject (l=1): pass through without caching
                                            space_srv_k = [space_srv_k; space_srv_e];
                                            space_var_k = [space_var_k; var];
                                            outrate(end+1,1) = GlobalConstants.Immediate;
                                            outprob(end+1,1) = ac{class,k}(1,l) * p(k);
                                        end
                                    else
                                        for l=2:(h+1) % iterate over all possible target lists
                                            listidx = l - 1; % l is accessCost column index, listidx is actual list (1-indexed)
                                            for r=1:m(listidx) % random position in list listidx
                                                varp = var;
                                                varp(cpos(listidx,r)) = k;
                                                space_srv_k = [space_srv_k; space_srv_e];
                                                space_var_k = [space_var_k; (varp)];
                                                outrate(end+1,1) = ac{class,k}(1,l) * p(k)/m(listidx) * GlobalConstants.Immediate;
                                            end
                                        end
                                    end
                            end
                        elseif posk <= sum(m(1:h-1)) % CACHE HIT in list i < h, move to list i+1
                            space_srv_e(hitclass(class)) = space_srv_e(hitclass(class)) + 1;
                            i = min(find(posk <= cumsum(m)));
                            j = posk - sum(m(1:i-1));

                            switch replacement_id
                                case ReplacementStrategy.FIFO
                                    if isSimulation
                                        varp = var;
                                        inew = i+probchoose(ac{class,k}(1+i,(1+i):end)/sum(ac{class,k}(1+i,(1+i):end)))-1; % can choose i
                                        if inew~=i
                                            varp(cpos(i,j)) = var(cpos(inew,m(inew)));
                                            varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                                            varp(cpos(inew,1)) = k;
                                        end
                                        %varp(cpos(i,j)) = var(cpos(i+1,m(i+1)));
                                        %varp(cpos(i+1,2):cpos(i+1,m(i+1))) = var(cpos(i+1,1):cpos(i+1,m(i+1)-1));
                                        %varp(cpos(i+1,1)) = k;

                                        space_srv_k = [space_srv_k; space_srv_e];
                                        space_var_k = [space_var_k; varp];
                                        outrate(end+1,1) = GlobalConstants.Immediate;
                                        outprob(end+1,1) = ac{class,k}(1+i,1+inew) * p(k);
                                    else
                                        for inew = i:h
                                            varp = var;
                                            varp(cpos(i,j)) = var(cpos(inew,m(inew)));
                                            varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                                            varp(cpos(inew,1)) = k;
                                            space_srv_k = [space_srv_k; space_srv_e];
                                            space_var_k = [space_var_k; varp];
                                            outrate(end+1,1) = ac{class,k}(1+i,1+inew) * p(k) * GlobalConstants.Immediate;
                                        end
                                    end
                                case ReplacementStrategy.RR
                                    if isSimulation
                                        inew = i+probchoose(ac{class,k}(1+i,(1+i):end)/sum(ac{class,k}(1+i,(1+i):end)))-1; % can choose i
                                        varp = var;
                                        r = randi(m(inew),1,1);
                                        varp(cpos(i,j)) = var(cpos(inew,r));
                                        varp(cpos(inew,r)) = k;
                                        space_srv_k = [space_srv_k; space_srv_e];
                                        space_var_k = [space_var_k; varp];
                                        outrate(end+1,1) = GlobalConstants.Immediate;
                                        outprob(end+1,1) = ac{class,k}(1+i,1+inew) * p(k)/m(inew);
                                    else
                                        for inew = i:h
                                            for r=1:m(inew) % random position in new list
                                                varp = var;
                                                varp(cpos(i,j)) = var(cpos(inew,r));
                                                varp(cpos(inew,r)) = k;
                                                space_srv_k = [space_srv_k; space_srv_e];
                                                space_var_k = [space_var_k; varp];
                                                outrate(end+1,1) = ac{class,k}(1+i,1+inew) * p(k)/m(inew) * GlobalConstants.Immediate;
                                            end
                                        end
                                    end
                                case {ReplacementStrategy.LRU, ReplacementStrategy.SFIFO}
                                    if isSimulation
                                        varp = var;
                                        inew = i+probchoose(ac{class,k}(1+i,(1+i):end)/sum(ac{class,k}(1+i,(1+i):end)))-1; % can choose i
                                        varp(cpos(i,2):cpos(i,j)) = var(cpos(i,1):cpos(i,j-1));
                                        varp(cpos(i,1)) = var(cpos(inew,m(inew)));
                                        varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                                        varp(cpos(inew,1)) = k;
                                        space_srv_k = [space_srv_k; space_srv_e];
                                        space_var_k = [space_var_k; varp];
                                        outrate(end+1,1) = GlobalConstants.Immediate;
                                        outprob(end+1,1) = ac{class,k}(1+i,1+inew) * p(k);
                                    else
                                        for inew = i:h
                                            varp = var;
                                            varp(cpos(i,2):cpos(i,j)) = var(cpos(i,1):cpos(i,j-1));
                                            varp(cpos(i,1)) = var(cpos(inew,m(inew)));
                                            varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                                            varp(cpos(inew,1)) = k;
                                            space_srv_k = [space_srv_k; space_srv_e];
                                            space_var_k = [space_var_k; varp];
                                            outrate(end+1,1) = ac{class,k}(1+i,1+inew) * p(k) * GlobalConstants.Immediate;
                                        end
                                    end
                            end
                        else % CACHE HIT in list h
                            space_srv_e(hitclass(class)) = space_srv_e(hitclass(class)) + 1;
                            i=h;
                            j = posk - sum(m(1:i-1));
                            switch replacement_id
                                case {ReplacementStrategy.RR, ReplacementStrategy.FIFO, ReplacementStrategy.SFIFO}
                                    space_srv_k = [space_srv_k; space_srv_e];
                                    space_var_k = [space_var_k; var];
                                    if isSimulation
                                        outrate(end+1,1) = GlobalConstants.Immediate;
                                        outprob(end+1,1) = p(k);
                                    else
                                        outrate(end+1,1) = p(k) * GlobalConstants.Immediate;
                                    end
                                case ReplacementStrategy.LRU
                                    varp = var;
                                    varp(cpos(h,2):cpos(h,j)) = var(cpos(h,1):cpos(h,j-1));
                                    varp(cpos(h,1)) = var(cpos(h,j));
                                    space_srv_k = [space_srv_k; space_srv_e];
                                    space_var_k = [space_var_k; varp];
                                    if isSimulation
                                        outrate(end+1,1) = GlobalConstants.Immediate;
                                        outprob(end+1,1) = p(k);
                                    else
                                        outrate(end+1,1) = p(k) * GlobalConstants.Immediate;
                                    end
                            end
                        end
                    end
                end
                % if state is unchanged, still add with rate 0
                outspace = [space_srv_k, space_var_k];
            end
        end
end
    function pos = cpos(i,j)
        % POS = CPOS(I,J)

        pos = sum(m(1:i-1)) + j;
    end

end

