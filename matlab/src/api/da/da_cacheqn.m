function [res, hitprob, missprob, it, sn, cacheinfo] = da_cacheqn(sn, missfun, netfun, options)
% [RES,HITPROB,MISSPROB,IT,SN,CACHEINFO] = DA_CACHEQN(SN, MISSFUN, NETFUN, OPTIONS)
%
% Shared decomposition-aggregation driver for integrated cache-queueing
% models. Alternates between (i) the solution of each cache in isolation,
% given the current per-class arrival rates, and (ii) the solution of the
% surrounding queueing network with the caches replaced by class switches
% routing according to the hit/miss probabilities, until the cache arrival
% rates reach a fixed point (driven by da_fpi with the 1-norm).
%
% SN:      NetworkStruct; mutated in place (cache nodes relabeled as
%          ClassSwitch, routing and visits refreshed) and returned.
% MISSFUN: MISSRATE = MISSFUN(GAMMA, M, LAMBDA_CACHE, CH) solves the
%          isolated cache and returns the per-class miss rates
%          (1 x nclasses); CH is the cache node parameter struct
%          (sn.nodeparam), e.g. to select the algorithm by replacement
%          strategy.
% NETFUN:  RES = NETFUN(SN) solves the surrounding queueing network; the
%          returned struct must include the field XN (1 x nclasses system
%          throughputs); all other fields are passed through to the caller.
% OPTIONS: solver options struct (iter_max, iter_tol).
%
% Returns the last NETFUN result RES, the per-cache hit/miss probabilities
% (length(caches) x nclasses, rows ordered as find(sn.nodetype==Cache)),
% the iteration count IT, the mutated SN, and CACHEINFO with the converged
% isolated-cache inputs per cache (fields node, gamma, m, lambda_cache,
% Rcost, strat) for per-item occupancy reporting.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

I = sn.nnodes;
K = sn.nclasses;

statefulNodes = find(sn.isstateful)';
statefulNodesClasses = [];
for ind = statefulNodes
    statefulNodesClasses(end+1:end+K) = ((ind-1)*K+1):(ind*K);
end
caches = find(sn.nodetype == NodeType.Cache);
ncaches = length(caches);

hitprob = zeros(ncaches, K);
missprob = zeros(ncaches, K);
res = struct();

cacheinfo = struct();
cacheinfo.node = caches;
cacheinfo.gamma = cell(1, ncaches);
cacheinfo.m = cell(1, ncaches);
cacheinfo.lambda_cache = cell(1, ncaches);
cacheinfo.Rcost = cell(1, ncaches);
cacheinfo.strat = cell(1, ncaches);

% initialization sweep: seed the cache arrival rates with random values and
% relabel the caches as class switches
lambda0 = zeros(1, K);
for cIdx = 1:ncaches
    ch = sn.nodeparam{caches(cIdx)};
    inputClass = find(ch.hitclass);
    lambda0(inputClass) = rand(1, length(inputClass));
    sn.nodetype(caches(cIdx)) = NodeType.ClassSwitch;
end

fpopts = options;
fpopts.config.da_norm = @(d) norm(d, 1);
[~, it] = da_fpi(@da_sweep, lambda0, fpopts);

    function [xnew, xref] = da_sweep(x, itnum) %#ok<INUSD>
        lambda = x;
        for cIdx = 1:ncaches
            ind = caches(cIdx);
            ch = sn.nodeparam{ind};
            hitClass = ch.hitclass;
            missClass = ch.missclass;
            inputClass = find(hitClass);

            % solution of isolated cache
            [gamma, lambda_cache, Rcost] = da_cache_isolate(ch, lambda);
            cacheinfo.gamma{cIdx} = gamma;
            cacheinfo.m{cIdx} = ch.itemcap;
            cacheinfo.lambda_cache{cIdx} = lambda_cache;
            cacheinfo.Rcost{cIdx} = Rcost;
            cacheinfo.strat{cIdx} = ch.replacestrat;

            missrate = missfun(gamma, ch.itemcap, lambda_cache, ch);
            missrate = reshape(missrate, 1, []); % row, as legacy indexed assignment did
            missprob(cIdx,:) = missrate ./ lambda; % NaN if no arrivals
            hitprob(cIdx,:) = 1 - missprob(cIdx,:);
            hitprob(isnan(hitprob)) = 0;
            missprob(isnan(missprob)) = 0;

            % bring back the isolated model results into the queueing model
            for r = inputClass
                sn.rtnodes((ind-1)*K+r,:) = 0;
                for jnd = 1:I
                    if sn.connmatrix(ind, jnd)
                        sn.rtnodes((ind-1)*K+r, (jnd-1)*K+hitClass(r)) = hitprob(cIdx, r);
                        sn.rtnodes((ind-1)*K+r, (jnd-1)*K+missClass(r)) = missprob(cIdx, r);
                    end
                end
            end
            sn.rt = dtmc_stochcomp(sn.rtnodes, statefulNodesClasses);
        end
        [visits, nodevisits, sn] = sn_refresh_visits(sn, sn.chains, sn.rt, sn.rtnodes);
        sn.visits = visits;
        sn.nodevisits = nodevisits;

        % aggregation step: solve the surrounding queueing network
        res = netfun(sn);

        % update the cache arrival rates from the network solution
        nv = cellsum(nodevisits);
        for cIdx = 1:ncaches
            ind = caches(cIdx);
            inputClass = find(sn.nodeparam{ind}.hitclass);
            for r = inputClass
                c = find(sn.chains(:,r));
                inchain = find(sn.chains(c,:));
                if sn.refclass(c) > 0
                    lambda(r) = sum(res.XN(inchain)) * nv(ind,r) / nv(sn.stationToNode(sn.refstat(r)), sn.refclass(c));
                else
                    lambda(r) = sum(res.XN(inchain)) * nv(ind,r) / nv(sn.stationToNode(sn.refstat(r)), r);
                end
            end
        end
        xnew = lambda;
        xref = x;
    end
end
