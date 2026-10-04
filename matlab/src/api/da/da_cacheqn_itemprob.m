function itemprob = da_cacheqn_itemprob(cacheinfo, ci)
% ITEMPROB = DA_CACHEQN_ITEMPROB(CACHEINFO, CI)
%
% Per-item occupancy [nitems x (lists+1)] of cache CI from the converged access
% factors returned by da_cacheqn; column 1 is the miss probability, columns
% 2..end the per-list ones. Returns [] when the cache carries no access factors.
%
% This is the EMBEDDED (per-request) occupancy: it is the stationary law of the
% cache-content chain seen at request instants, which coincides with the
% time-stationary one only under PASTA. SolverCTMC reports the time-weighted
% counterpart instead.
%
% The RR/FIFO exact recursion is skipped past 10 items and NaN reported, since
% an approximation of it is not a distribution; see _kb/09-ldes-and-cache.md.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

itemprob = [];
if isempty(cacheinfo.gamma{ci})
    return
end
ni = size(cacheinfo.gamma{ci},1);
hi = numel(cacheinfo.m{ci});
if cacheinfo.strat{ci} == ReplacementStrategy.LRU
    itemprob = cache_ttl_lrua(cacheinfo.lambda_cache{ci}, cacheinfo.Rcost{ci}, cacheinfo.m{ci});
elseif ni > 10
    line_warning(mfilename, 'Per-item cache occupancy (getAvgItemTable) requires the exact algorithm for RR/FIFO and is skipped for caches with more than 10 items (%d items); reporting NaN.', ni);
    itemprob = NaN(ni, hi+1);
else
    itemprob = cache_prob_erec(cacheinfo.gamma{ci}, cacheinfo.m{ci});
end
end
