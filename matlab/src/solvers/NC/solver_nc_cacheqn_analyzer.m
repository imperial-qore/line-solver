function [QN,UN,RN,TN,CN,XN,lG,hitprob,missprob,runtime,it,method] = solver_nc_cacheqn_analyzer(self, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER] = SOLVER_NC_CACHEQN_ANALYZER(SELF, OPTIONS)
%
% Integrated cache-queueing analyzer: delegates the decomposition-
% aggregation alternation between the isolated caches and the queueing
% network to da_cacheqn, supplying the NC-specific isolated-cache miss
% algorithm (exact recursion or SPM approximation) and network solver.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

snorig = self.model.getStruct;
sn = snorig;
K = sn.nclasses;

line_debug('NC cacheqn analyzer starting: method=%s, nclasses=%d', options.method, K);

caches = find(sn.nodetype == NodeType.Cache);
for ci = 1:length(caches)
    ch = sn.nodeparam{caches(ci)};
    if ch.nitems < sum(ch.itemcap) + 2
        line_error(mfilename,'NC requires the number of items to exceed the cache capacity at least by 2.');
    end
end

switch options.method
    case 'exact'
        missfun = @miss_exact;
        method = 'exact';
    otherwise
        missfun = @miss_spm;
        method = 'spm';
end

[res, hitprob_pc, missprob_pc, it] = da_cacheqn(sn, missfun, @netsolve, options);
QN = res.Q; UN = res.U; RN = res.R; TN = res.T; CN = res.C; XN = res.X;
lG = res.lG; runtime = res.runtime;

% legacy contract: hit/miss probabilities indexed by node row
hitprob = zeros(length(caches), K);
missprob = zeros(length(caches), K);
for ci = 1:length(caches)
    hitprob(caches(ci),:) = hitprob_pc(ci,:);
    missprob(caches(ci),:) = missprob_pc(ci,:);
end

    function missrate = miss_exact(gamma, m, lambda_cache, ~)
        u = size(lambda_cache, 1);
        pij = cache_prob_erec(gamma, m);
        missrate = zeros(1, u);
        for v = 1:u
            missrate(v) = lambda_cache(v,:,1) * pij(:,1);
        end
    end

    function missrate = miss_spm(gamma, m, lambda_cache, ~)
        line_debug('Default method: using SPM approximation for cache\n');
        line_debug('Using SPM approximation, calling cache_miss_spm');
        [~, missrate] = cache_miss_spm(gamma, m, lambda_cache);
    end

    function res = netsolve(snit)
        res = struct();
        if ~isempty(snit.lldscaling) || ~isempty(snit.cdscaling) || ~isempty(snit.jdscaling)
            [res.Q,res.U,res.R,res.T,res.C,res.X,res.lG,res.runtime] = solver_ncld_analyzer(snit, options);
        else
            [res.Q,res.U,res.R,res.T,res.C,res.X,res.lG,res.runtime] = solver_nc_analyzer(snit, options);
        end
        res.XN = res.X;
    end
end
