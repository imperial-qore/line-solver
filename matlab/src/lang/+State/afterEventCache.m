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
            % A retrieval-class job departs the cache only to BEGIN a retrieval
            % (cache -> queue); record the item in the per-item occupancy bitmap.
            beginBlocked = false;
            if isfield(sn.nodeparam{ind}, 'retrievalClassIndices')
                rciDep = sn.nodeparam{ind}.retrievalClassIndices;
            else
                rciDep = [];
            end
            if any(rciDep == class)
                if isfield(sn.nodeparam{ind}, 'totalCacheCapacity')
                    tccDep = sn.nodeparam{ind}.totalCacheCapacity;
                else
                    tccDep = sum(sn.nodeparam{ind}.itemcap);
                end
                pDep = sn.nodeparam{ind}.pread{class};
                item = find(pDep == 1, 1, 'first'); % retrieval class reads its one-hot item
                if ~isempty(item) && (tccDep + item <= size(space_var,2))
                    if any(space_var(:, tccDep + item) ~= 0)
                        beginBlocked = true; % already retrieving this item -> no transition
                    else
                        space_var(:, tccDep + item) = 1;
                    end
                end
            end
            if ~beginBlocked
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
        end
    case EventType.READ
        n = sn.nodeparam{ind}.nitems; % n items
        m = sn.nodeparam{ind}.itemcap; % capacity
        ac = sn.nodeparam{ind}.accost; % access cost
        hitclass = sn.nodeparam{ind}.hitclass;
        missclass = sn.nodeparam{ind}.missclass;
        h = length(m);
        replacement_id = sn.nodeparam{ind}.replacestrat;
        % retrieval-system parameters (set by Cache.setRetrievalSystem; defaults below)
        if isfield(sn.nodeparam{ind}, 'totalCacheCapacity')
            totalCacheCapacity = sn.nodeparam{ind}.totalCacheCapacity;
        else
            totalCacheCapacity = sum(m);
        end
        if isfield(sn.nodeparam{ind}, 'retrievalClasses')
            retrievalClasses = sn.nodeparam{ind}.retrievalClasses;
        else
            retrievalClasses = [];
        end
        if isfield(sn.nodeparam{ind}, 'retrievalClassIndices')
            retrievalClassIndices = sn.nodeparam{ind}.retrievalClassIndices;
        else
            retrievalClassIndices = [];
        end
        if space_srv(class)>0 && sum(space_srv)==1 %  a job of class is in
            p = sn.nodeparam{ind}.pread{class};
            en =  space_srv(:,class) > 0;
            space_srv_k = [];
            space_var_k = [];
            outrate = [];
            if any(en)
                for e=find(en)'
                    isFromRetrieval = any(retrievalClassIndices == class);
                    if isSimulation || isFromRetrieval
                        % pick one item (a returning retrieval reads exactly its own item)
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
                        % posk is searched only in the cache-contents region (1..totalCacheCapacity).
                        % The trailing retrieval-system slots (if any) are out of scope here.
                        posk = find(k==var(1:totalCacheCapacity),1,'first');
                        % A retrieval-complete READ always completes the miss that started
                        % the retrieval, so it must take the cache-miss branch even if the
                        % state-space enumeration produced a (then unreachable) state in
                        % which item k is already cached. This also avoids indexing the
                        % undefined hitClass of a retrieval-complete class.
                        if isFromRetrieval
                            posk = [];
                        end

                        if isempty(posk) % CACHE MISS or RETRIEVAL begin/return
                            % Retrieval-system occupancy bitmap: column (totalCacheCapacity+k)
                            % is non-zero iff item k is currently being retrieved.
                            inRetrieval = (totalCacheCapacity + k <= size(var,2)) && var(totalCacheCapacity + k) ~= 0;
                            rClass = -1;
                            if ~isempty(retrievalClasses) && ...
                                    class <= size(retrievalClasses,2) && ...
                                    k <= size(retrievalClasses,1)
                                rClass = retrievalClasses(k, class);
                            end

                            % A returning retrieval that is not recorded in the bitmap is an
                            % unreachable event-loop artifact; do not continue it.
                            if isFromRetrieval && ~inRetrieval
                                continue
                            end

                            % Begin a retrieval: the job switches to the retrieval class for
                            % item k. The occupancy bitmap is set on the subsequent DEP (when
                            % the job departs the cache for the queue), so the cache contents
                            % are unchanged here. A concurrent request for an item already
                            % being retrieved is served by the in-flight retrieval and exits
                            % (its class is not incremented, so the job is absorbed).
                            if ~isFromRetrieval && rClass ~= -1
                                if ~inRetrieval
                                    space_srv_e(rClass) = space_srv_e(rClass) + 1;
                                end
                                space_srv_k = [space_srv_k; space_srv_e];
                                space_var_k = [space_var_k; var];
                                if isSimulation
                                    outrate(end+1,1) = GlobalConstants.Immediate;
                                    outprob(end+1,1) = p(k);
                                else
                                    outrate(end+1,1) = p(k) * GlobalConstants.Immediate;
                                end
                                continue
                            end

                            % Item has now been retrieved (or there is no retrieval system):
                            % mark it as a miss and clear its retrieval-system bit.
                            space_srv_e(missclass(class)) = space_srv_e(missclass(class)) + 1;
                            if totalCacheCapacity + k <= size(var,2)
                                var(totalCacheCapacity + k) = 0;
                            end
                            switch replacement_id
                                case {ReplacementStrategy.FIFO, ReplacementStrategy.LRU, ReplacementStrategy.SFIFO, ReplacementStrategy.HLRU}
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
                                        % Cache reject (l=1): pass through without caching
                                        if ac{class,k}(1,1) > 0
                                            space_srv_k = [space_srv_k; space_srv_e];
                                            space_var_k = [space_var_k; var];
                                            outrate(end+1,1) = ac{class,k}(1,1) * p(k) * GlobalConstants.Immediate;
                                            outprob(end+1,1) = 1;
                                        end
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
                                        % Cache reject (l=1): pass through without caching
                                        if ac{class,k}(1,1) > 0
                                            space_srv_k = [space_srv_k; space_srv_e];
                                            space_var_k = [space_var_k; var];
                                            outrate(end+1,1) = ac{class,k}(1,1) * p(k) * GlobalConstants.Immediate;
                                        end
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
                                case ReplacementStrategy.QLRU
                                    % q-LRU: on a miss the item is admitted (LRU head insert)
                                    % with probability q, else it passes through uncached.
                                    if isfield(sn.nodeparam{ind},'qlru')
                                        qadm = sn.nodeparam{ind}.qlru;
                                    else
                                        qadm = 1.0;
                                    end
                                    if isSimulation
                                        listidx = l - 1;
                                        if listidx > 0 && rand <= qadm
                                            varp = var;
                                            varp(cpos(listidx,2):cpos(listidx,m(listidx))) = var(cpos(listidx,1):cpos(listidx,m(listidx)-1));
                                            varp(cpos(listidx,1)) = k;
                                            space_srv_k = [space_srv_k; space_srv_e];
                                            space_var_k = [space_var_k; varp];
                                            outrate(end+1,1) = GlobalConstants.Immediate;
                                            outprob(end+1,1) = ac{class,k}(1,l) * p(k);
                                        else
                                            space_srv_k = [space_srv_k; space_srv_e];
                                            space_var_k = [space_var_k; var];
                                            outrate(end+1,1) = GlobalConstants.Immediate;
                                            outprob(end+1,1) = ac{class,k}(1,l) * p(k);
                                        end
                                    else
                                        % pass-through mass: structural reject plus (1-q) non-admission
                                        rejw = ac{class,k}(1,1) + (1-qadm)*(1 - ac{class,k}(1,1));
                                        if rejw > 0
                                            space_srv_k = [space_srv_k; space_srv_e];
                                            space_var_k = [space_var_k; var];
                                            outrate(end+1,1) = rejw * p(k) * GlobalConstants.Immediate;
                                            outprob(end+1,1) = 1;
                                        end
                                        for l=2:(h+1)
                                            listidx = l - 1;
                                            varp = var;
                                            varp(cpos(listidx,2):cpos(listidx,m(listidx))) = var(cpos(listidx,1):cpos(listidx,m(listidx)-1));
                                            varp(cpos(listidx,1)) = k;
                                            space_srv_k = [space_srv_k; space_srv_e];
                                            space_var_k = [space_var_k; varp];
                                            outrate(end+1,1) = qadm * ac{class,k}(1,l) * p(k) * GlobalConstants.Immediate;
                                            outprob(end+1,1) = 1;
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
                                case {ReplacementStrategy.LRU, ReplacementStrategy.SFIFO, ReplacementStrategy.HLRU, ReplacementStrategy.QLRU}
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
                                case {ReplacementStrategy.LRU, ReplacementStrategy.HLRU, ReplacementStrategy.QLRU}
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

