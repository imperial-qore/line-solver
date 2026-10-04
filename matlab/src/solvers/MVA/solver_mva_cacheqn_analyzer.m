function [Q,U,R,T,C,X,lG,hitprob,missprob,runtime,it] = solver_mva_cacheqn_analyzer(self, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER] = SOLVER_MVA_CACHEQN_ANALYZER(SELF, OPTIONS)
%
% Integrated cache-queueing analyzer: delegates the decomposition-
% aggregation alternation between the isolated caches and the queueing
% network to da_cacheqn, supplying the MVA-specific isolated-cache miss
% algorithm (exact cache_mva or FPI approximation) and network solver.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

snorig = self.model.getStruct;
sn = snorig;
K = sn.nclasses;

line_debug('MVA cacheqn analyzer starting: method=%s, nclasses=%d', options.method, K);

caches = find(sn.nodetype == NodeType.Cache);

switch options.method
    case 'exact'
        missfun = @miss_exact;
    otherwise
        missfun = @miss_fpi;
end

[res, hitprob_pc, missprob_pc, it, ~, cacheinfo] = da_cacheqn(sn, missfun, @netsolve, options);
Q = res.Q; U = res.U; R = res.R; T = res.T; C = res.C; X = res.X;
lG = res.lG; runtime = res.runtime;

% legacy contract: hit/miss probabilities indexed by node row
hitprob = zeros(length(caches), K);
missprob = zeros(length(caches), K);
for ci = 1:length(caches)
    hitprob(caches(ci),:) = hitprob_pc(ci,:);
    missprob(caches(ci),:) = missprob_pc(ci,:);
end

% per-item occupancy [nitems x (lists+1)] from converged access factors (RR/FIFO
% exact recursion skipped, NaN, for >10 items); see _kb/09-ldes-and-cache.md
for ci = 1:length(caches)
    itemprob = da_cacheqn_itemprob(cacheinfo, ci);
    if ~isempty(itemprob)
        self.model.nodes{caches(ci)}.setResultItemProb(itemprob);
    end
end

    function missrate = miss_exact(gamma, m, lambda_cache, ~)
        u = size(lambda_cache, 1);
        [~,~,pij] = cache_mva(gamma, m);
        pij = [abs(1-sum(pij,2)), pij];
        missrate = zeros(1, u);
        for v = 1:u
            missrate(v) = lambda_cache(v,:,1) * pij(:,1);
        end
    end

    function missrate = miss_fpi(gamma, m, lambda_cache, ~)
        line_debug('Default method: using FPI approximation for cache\n');
        line_debug('Using FPI approximation, calling cache_miss_fpi');
        [~, missrate] = cache_miss_fpi(gamma, m, lambda_cache);
    end

    function res = netsolve(snit)
        res = struct();
        switch options.method
            case {'aba.upper', 'aba.lower', 'bjb.upper', 'bjb.lower', 'pb.upper', 'pb.lower', 'gb.upper', 'gb.lower', 'sb.upper', 'sb.lower'}
                [res.Q,res.U,res.R,res.T,res.C,res.X,res.lG,res.runtime] = solver_ba_analyzer(snit, options);
            otherwise
                if ~isempty(snit.lldscaling) || ~isempty(snit.cdscaling) || ~isempty(snit.jdscaling)
                    [res.Q,res.U,res.R,res.T,res.C,res.X,res.lG,res.runtime] = solver_mvald_analyzer(snit, options);
                else
                    [res.Q,res.U,res.R,res.T,res.C,res.X,res.lG,res.runtime] = solver_mva_analyzer(snit, options);
                end
        end
        res.XN = res.X;
    end
end
