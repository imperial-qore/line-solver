function [rt, rtfun, rtnodes, sn] = refreshRoutingMatrix(self, rates)
% [RT, RTFUN, CSMASK, RTNODES, SN] = REFRESHROUTINGMATRIX(RATES)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.sn;
if nargin == 1
    if isempty(sn)
        line_error(mfilename,'refreshRoutingMatrix cannot retrieve station rates, pass them as an input parameters.');
    else
        rates = sn.rates;
    end
end
M = sn.nnodes;
K = sn.nclasses;
arvRates = zeros(1,K);
stateful = find(sn.isstateful)';

indSource = find(sn.nodetype == NodeType.Source);
indOpenClasses = find(sn.njobs == Inf);
for r = indOpenClasses
    arvRates(r) = rates(sn.nodeToStation(indSource),r);
end

[rt, rtnodes, linksmat, chains] = self.getRoutingMatrix(arvRates);
sn = self.sn;
sn.chains = chains;

if self.enableChecks
    for r=1:K
        if all(sn.routing(:,r) == -1)
            line_error(mfilename,sprintf('Routing strategy in class %d is unspecified at all nodes.',r));
        end
    end
end

isStateDep = any(sn.isstatedep(:,3));

rnodefuncell = cell(M*K,M*K);

if isStateDep
    for ind=1:M % from
        for jnd=1:M % to
            for r=1:K
                for s=1:K
                    if sn.isstatedep(ind,3)
                        switch sn.routing(ind,r)
                            case RoutingStrategy.RROBIN
                                rnodefuncell{(ind-1)*K+r, (jnd-1)*K+s} = @(state_before, state_after) sub_rr(ind, jnd, r, s, linksmat, state_before, state_after);
                            case RoutingStrategy.WRROBIN
                                rnodefuncell{(ind-1)*K+r, (jnd-1)*K+s} = @(state_before, state_after) sub_wrr(ind, jnd, r, s, linksmat, state_before, state_after);
                            case RoutingStrategy.JSQ
                                rnodefuncell{(ind-1)*K+r, (jnd-1)*K+s} = @(state_before, state_after) sub_jsq(ind, jnd, r, s, linksmat, state_before, state_after);
                            case RoutingStrategy.SQ
                                rnodefuncell{(ind-1)*K+r, (jnd-1)*K+s} = @(state_before, state_after) sub_sq(ind, jnd, r, s, linksmat, state_before, state_after);
                            case RoutingStrategy.SDR
                                rnodefuncell{(ind-1)*K+r, (jnd-1)*K+s} = @(state_before, state_after) sub_sdr(ind, jnd, r, s, linksmat, state_before, state_after);
                            otherwise
                                rnodefuncell{(ind-1)*K+r, (jnd-1)*K+s} = @(~,~) rtnodes((ind-1)*K+r, (jnd-1)*K+s);
                        end
                    else
                        rnodefuncell{(ind-1)*K+r, (jnd-1)*K+s} = @(~,~) rtnodes((ind-1)*K+r, (jnd-1)*K+s);
                    end
                end
            end
        end
    end
end

statefulNodesClasses = [];
for ind=getIndexStatefulNodes(self)
    statefulNodesClasses(end+1:end+K)= ((ind-1)*K+1):(ind*K);
end

% we now generate the node routing matrix for the given state and then
% lump the states for non-stateful nodes so that run gives the routing
% table for stateful nodes only
statefulNodesClasses = [];
for ind=stateful
    statefulNodesClasses(end+1:end+K)= ((ind-1)*K+1):(ind*K);
end

if isStateDep
    rtfunraw = @(state_before, state_after) dtmc_stochcomp(cell2mat(cellfun(@(f) f(state_before, state_after), rnodefuncell,'UniformOutput',false)), statefulNodesClasses);
    rtfun = rtfunraw;
    %rtfun = memoize(rtfunraw); % memoize to reduce the number of stoch comp calls
    %rtfun.CacheSize = 6000^2;
else
    rtfun = @(state_before, state_after) dtmc_stochcomp(rtnodes, statefulNodesClasses);
end

nchains = size(chains,1);
inchain = cell(1,nchains);
for c=1:nchains
    inchain{c} = find(chains(c,:));
end

sn.rt = rt;
sn.rtnodes = rtnodes;
%sn.rtorig = RoutingMatrix.rtnodes2rtorig(sn); %% causes issues to
%JLINE.from_line_links in convering example_closedModel_6
sn.rtfun = rtfun;
sn.chains = chains;
sn.nchains = nchains;
sn.inchain = inchain;
for c=1:sn.nchains
    if range(sn.refstat(inchain{c}))>0
        line_error(mfilename,sprintf('Classes within chain %d (classes: %s) have different reference stations.',c,mat2str(find(sn.chains(c,:)))));
    end
end
self.sn = sn;

    function p = sub_rr(ind, jnd, r, s, linksmat, state_before, state_after)
        % P = SUB_RR(IND, JND, R, S, LINKSMAT, STATE_BEFORE, STATE_AFTER)

        isf = self.sn.nodeToStateful(ind);
        if isempty(state_before{isf})
            p = min(linksmat(ind,jnd),1);
        else
            if r==s
                p = double(state_after{isf}(sub_routeslot(ind, r, state_after{isf}))==jnd);
            else
                p = 0;
            end
        end
    end

    function slot = sub_routeslot(ind, r, state_i)
        % SLOT = SUB_ROUTESLOT(IND, R, STATE_I)
        % Column of STATE_I holding the round-robin pointer of class R at node
        % IND.
        %
        % The per-class routing pointers are appended in class order, and only
        % for the classes that actually route round-robin, followed by the
        % node-level block (nvars column 2R+1) if the node has one. The pointer
        % is therefore located by counting back from the end of the state: over
        % the node-level block, then over the pointers of the classes after r.
        %
        % Reading it as state_i(end-R+r) instead, as this did, assumes that the
        % pointers are exactly the last R columns. That silently reads the wrong
        % column whenever the node carries a trailing block of its own (the
        % polling controller of a POLLING station), and also whenever only some
        % of the classes route round-robin, in which case a pointer is confused
        % with a service phase.
        R = self.sn.nclasses;
        nodeblock = self.sn.nvars(ind, 2*R+1);
        after = sum(self.sn.nvars(ind, (R+r+1):(2*R)));
        slot = numel(state_i) - nodeblock - after;
    end

    function p = sub_wrr(ind, jnd, r, s, linksmat, state_before, state_after)
        % P = SUB_WRR(IND, JND, R, S, LINKSMAT, STATE_BEFORE, STATE_AFTER)
        % WRR slot holds a POSITION in weighted_outlinks; map it to the
        % destination node index before comparing.

        isf = self.sn.nodeToStateful(ind);
        if isempty(state_before{isf})
            p = min(linksmat(ind,jnd),1);
        else
            if r==s
                pos = state_after{isf}(sub_routeslot(ind, r, state_after{isf}));
                np = self.sn.nodeparam{ind}{r};
                if isfield(np,'weighted_outlinks') && ~isempty(np.weighted_outlinks) ...
                        && pos >= 1 && pos <= length(np.weighted_outlinks)
                    dest = np.weighted_outlinks(pos);
                elseif pos >= 1 && pos <= length(np.outlinks)
                    dest = np.outlinks(pos);
                else
                    p = 0; return;
                end
                p = double(dest == jnd);
            else
                p = 0;
            end
        end
    end

    function p = sub_jsq(ind, jnd, r, s, linksmat, state_before, state_after) %#ok<INUSD>
        % P = SUB_JSQ(IND, JND, R, S, LINKSMAT, STATE_BEFORE, STATE_AFTER) %#OK<INUSD>

        isf = self.sn.nodeToStateful(ind);
        if isempty(state_before{isf})
            p = min(linksmat(ind,jnd),1);
        else
            if r==s
                n = Inf*ones(1,self.sn.nnodes);
                for knd=1:self.sn.nnodes
                    if linksmat(ind,knd)
                        ksf = self.sn.nodeToStateful(knd);
                        n(knd) = State.toMarginal(self.sn, knd, state_before{ksf});
                    end
                end
                if n(jnd) == min(n)
                    p = 1 / sum(n == min(n));
                else
                    p = 0;
                end
            else
                p = 0;
            end
        end
    end


    function p = sub_sdr(ind, jnd, r, s, linksmat, state_before, state_after) %#ok<INUSD>
        % P = SUB_SDR(IND, JND, R, S, LINKSMAT, STATE_BEFORE, STATE_AFTER)
        % Krzesinski (1987) product-form state-dependent routing, eq. (10).
        % IND is the entry center e of Q(V,V); the probability of proceeding to
        % the entry center of a branch is a function of the total branch and
        % subnetwork populations, and the residual mass returns the customer to
        % the departure center d, which is the busy form of waiting of Sec. 2.5.

        isf = self.sn.nodeToStateful(ind);
        if isempty(state_before{isf})
            p = min(linksmat(ind,jnd),1);
            return
        end
        if r ~= s
            p = 0; return;
        end
        sdr = self.sn.nodeparam{ind}{r}.sdr;
        n = zeros(1, self.sn.nnodes);
        for knd = 1:self.sn.nnodes
            ksf = self.sn.nodeToStateful(knd);
            if ~isnan(ksf) && ksf >= 1 && ~isempty(state_before{ksf})
                n(knd) = State.toMarginal(self.sn, knd, state_before{ksf});
            end
        end
        [Pb, Ped] = pfqn_sdrprob(sdr, n);
        p = 0;
        for b = 2:numel(sdr.branch)
            if sdr.entryOf(b) == jnd
                p = p + Pb(b);
            end
        end
        if sdr.departure == jnd
            p = p + Ped;
        end
    end

    function p = sub_sq(ind, jnd, r, s, linksmat, state_before, state_after) %#ok<INUSD>
        % P = SUB_SQ(IND, JND, R, S, LINKSMAT, STATE_BEFORE, STATE_AFTER)
        % SQ(d), shortest queue of d: matches LDES semantics in
        % Solver_ssj.kt:selectSQDestination.
        %
        % Enumerate the ndest^d ordered tuples obtained by sampling d
        % destinations WITH replacement, break ties by first occurrence in the
        % tuple, and return the marginal probability that jnd is chosen.
        % Dispatcher memory is not supported, so the closure is a function of
        % the queue lengths alone and carries no auxiliary state.
        %
        % Memoization: the marginal-probability vector p_dest(jnd) depends only
        % on (n[1..ndest], d). Across (ind, r) we cache by the per-call
        % queue-length signature so each unique (n, d) is enumerated once and
        % the per-jnd probabilities are reused.
        persistent sqCache
        if isempty(sqCache)
            sqCache = configureDictionary('string','cell');
        end

        isf = self.sn.nodeToStateful(ind);
        if isempty(state_before{isf})
            p = min(linksmat(ind,jnd),1);
            return;
        end
        if r ~= s
            p = 0; return;
        end
        eligible = find(linksmat(ind,:));
        ndest = numel(eligible);
        if ndest == 0 || ~any(eligible == jnd)
            p = 0; return;
        end
        np = self.sn.nodeparam{ind}{r};
        if isfield(np,'d') && ~isempty(np.d)
            d = np.d;
        else
            d = 2;
        end
        d = max(1, min(d, ndest));
        n = zeros(1, ndest);
        for i = 1:ndest
            knd = eligible(i);
            ksf = self.sn.nodeToStateful(knd);
            n(i) = State.toMarginal(self.sn, knd, state_before{ksf});
        end
        jnd_pos = find(eligible == jnd, 1);

        % Cache key: (per-class) d, ndest, n vector. Same vector reused across
        % all jnd values for this (ind, r) pair.
        cacheKey = sprintf('d%d|nd%d|n%s', d, ndest, mat2str(n));
        if isKey(sqCache, cacheKey)
            pVec = sqCache{cacheKey};
        else
            pVec = zeros(1, ndest);
            n_tuples = ndest^d;
            for ti = 0:(n_tuples-1)
                tuple = zeros(1, d);
                rem = ti;
                for c = 1:d
                    tuple(c) = mod(rem, ndest) + 1;
                    rem = floor(rem / ndest);
                end
                sub_n = n(tuple);
                minval = min(sub_n);
                winner_pos_in_tuple = find(sub_n == minval, 1);
                winner = tuple(winner_pos_in_tuple);
                pVec(winner) = pVec(winner) + 1;
            end
            pVec = pVec / n_tuples;
            % Bound cache size to avoid unbounded growth in long simulations.
            if numEntries(sqCache) > 50000
                sqCache = remove(sqCache, keys(sqCache));
            end
            sqCache{cacheKey} = pVec;
        end
        p = pVec(jnd_pos);
    end

end


































% function [rt, rtfun, rtnodes, sn] = refreshRoutingMatrix(self, rates)
% % [RT, RTFUN, CSMASK, RTNODES, SN] = REFRESHROUTINGMATRIX(RATES)
% %
% % Copyright (c) 2012-2026, Imperial College London
% % All rights reserved.
% 
% sn = self.sn;
% if nargin == 1
%     if isempty(sn)
%         line_error(mfilename,'refreshRoutingMatrix cannot retrieve station rates, pass them as an input parameters.');
%     else
%         rates = sn.rates;
%     end
% end
% M = sn.nnodes;
% K = sn.nclasses;
% arvRates = zeros(1,K);
% stateful = find(sn.isstateful)';
% 
% indSource = find(sn.nodetype == NodeType.Source);
% indOpenClasses = find(sn.njobs == Inf);
% for r = indOpenClasses
%     arvRates(r) = rates(sn.nodeToStation(indSource),r);
% end
% 
% [rt, rtnodes, linksmat, chains] = self.getRoutingMatrix(arvRates);
% sn = self.sn;
% sn.chains = chains;
% 
% if self.enableChecks
%     for r=1:K
%         if all(sn.routing(:,r) == -1)
%             line_error(mfilename,sprintf('Routing strategy in class %d is unspecified at all nodes.',r));
%         end
%     end
% end
% 
% isStateDep = any(sn.isstatedep(:,3));
% 
% rnodefuncell = cell(M*K,M*K);
% 
% if isStateDep
%     for ind=1:M % from
%         for jnd=1:M % to
%             for r=1:K
%                 for s=1:K
%                     if sn.isstatedep(ind,3)
%                         switch sn.routing(ind,r)
%                             case RoutingStrategy.RROBIN
%                                 rnodefuncell{(ind-1)*K+r, (jnd-1)*K+s} = @(state_before, state_after) sub_rr(ind, jnd, r, s, linksmat, state_before, state_after);
%                             case RoutingStrategy.WRROBIN
%                                 rnodefuncell{(ind-1)*K+r, (jnd-1)*K+s} = @(state_before, state_after) sub_wrr(ind, jnd, r, s, linksmat, state_before, state_after);
%                             case RoutingStrategy.JSQ
%                                 rnodefuncell{(ind-1)*K+r, (jnd-1)*K+s} = @(state_before, state_after) sub_jsq(ind, jnd, r, s, linksmat, state_before, state_after);
%                             otherwise
%                                 rnodefuncell{(ind-1)*K+r, (jnd-1)*K+s} = @(~,~) rtnodes((ind-1)*K+r, (jnd-1)*K+s);
%                         end
%                     else
%                         rnodefuncell{(ind-1)*K+r, (jnd-1)*K+s} = @(~,~) rtnodes((ind-1)*K+r, (jnd-1)*K+s);
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% statefulNodesClasses = [];
% for ind=getIndexStatefulNodes(self)
%     statefulNodesClasses(end+1:end+K)= ((ind-1)*K+1):(ind*K);
% end
% 
% % we now generate the node routing matrix for the given state and then
% % lump the states for non-stateful nodes so that run gives the routing
% % table for stateful nodes only
% statefulNodesClasses = [];
% for ind=stateful
%     statefulNodesClasses(end+1:end+K)= ((ind-1)*K+1):(ind*K);
% end
% 
% if isStateDep
%     rtfunraw = @(state_before, state_after) dtmc_stochcomp(cell2mat(cellfun(@(f) f(state_before, state_after), rnodefuncell,'UniformOutput',false)), statefulNodesClasses);
%     rtfun = rtfunraw;
%     %rtfun = memoize(rtfunraw); % memoize to reduce the number of stoch comp calls
%     %rtfun.CacheSize = 6000^2;
% else
%     rtfun = @(state_before, state_after) dtmc_stochcomp(rtnodes, statefulNodesClasses);
% end
% 
% nchains = size(chains,1);
% inchain = cell(1,nchains);
% for c=1:nchains
%     inchain{c} = find(chains(c,:));
% end
% 
% sn.rt = rt;
% sn.rtnodes = rtnodes;
% sn.rtfun = rtfun;
% sn.chains = chains;
% sn.nchains = nchains;
% sn.inchain = inchain;
% for c=1:sn.nchains
%     if range(sn.refstat(inchain{c}))>0
%         line_error(mfilename,sprintf('Classes within chain %d (classes: %s) have different reference stations.',c,mat2str(find(sn.chains(c,:)))));
%     end
% end
% self.sn = sn;
% 
%     function p = sub_rr(ind, jnd, r, s, linksmat, state_before, state_after)
%         % P = SUB_RR(IND, JND, R, S, LINKSMAT, STATE_BEFORE, STATE_AFTER)
% 
%         R = sn.nclasses;
%         isf = sn.nodeToStateful(ind);
%         if isempty(state_before{isf})
%             p = min(linksmat(ind,jnd),1);
%         else
%             if r==s
%                 p = double(state_after{isf}(end-R+r)==jnd);
%             else
%                 p = 0;
%             end
%         end
%     end
% 
%     function p = sub_wrr(ind, jnd, r, s, linksmat, state_before, state_after)
%         % P = SUB_WRR(IND, JND, R, S, LINKSMAT, STATE_BEFORE, STATE_AFTER)
% 
%         R = sn.nclasses;
%         isf = sn.nodeToStateful(ind);
%         if isempty(state_before{isf})
%             p = min(linksmat(ind,jnd),1);
%         else
%             if r==s
%                 p = double(state_after{isf}(end-R+r)==jnd);
%             else
%                 p = 0;
%             end
%         end
%     end
% 
%     function p = sub_jsq(ind, jnd, r, s, linksmat, state_before, state_after) %#ok<INUSD>
%         % P = SUB_JSQ(IND, JND, R, S, LINKSMAT, STATE_BEFORE, STATE_AFTER) %#OK<INUSD>
% 
%         isf = sn.nodeToStateful(ind);
%         if isempty(state_before{isf})
%             p = min(linksmat(ind,jnd),1);
%         else
%             if r==s
%                 n = Inf*ones(1,sn.nnodes);
%                 for knd=1:sn.nnodes
%                     if linksmat(ind,knd)
%                         ksf = sn.nodeToStateful(knd);
%                         n(knd) = State.toMarginal(sn, knd, state_before{ksf});
%                     end
%                 end
%                 if n(jnd) == min(n)
%                     p = 1 / sum(n == min(n));
%                 else
%                     p = 0;
%                 end
%             else
%                 p = 0;
%             end
%         end
%     end
% 
% end
