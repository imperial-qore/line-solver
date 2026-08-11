function [QN,UN,RN,TN,CN,XN,lG,pij,runtime,method,hitproblist,itemprob] = solver_nc_cache_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,PIJ,RUNTIME,METHOD,HITPROBLIST] = SOLVER_NC_CACHE_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0=tic;
QN = []; UN = [];
RN = []; TN = [];
CN = [];
XN = zeros(1,sn.nclasses);
lG = NaN;
iter = NaN;

line_debug('NC cache analyzer starting: method=%s, nclasses=%d', options.method, sn.nclasses);

source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
sourceRate = sn.rates(source_ist,:);
sourceRate(isnan(sourceRate)) = 0;
TN(source_ist,:) = sourceRate;

ch = sn.nodeparam{sn.nodetype == NodeType.Cache};

m = ch.itemcap;
n = ch.nitems;

if n<m+2
    line_error(mfilename,'NC requires the number of items to exceed the cache capacity at least by 2.');
end

h = length(m);
u = sn.nclasses;
lambda = zeros(u,n,h);

for v=1:u
    for k=1:n
        for l=1:(h+1)
            if ~isnan(ch.pread{v})
                lambda(v,k,l) = sourceRate(v) * ch.pread{v}(k);
            end
        end
    end
end

R = ch.accost;
if isempty(R)
    % Default linear cache routing: items flow from list l to list l+1
    R = cell(u, n);
    for v = 1:u
        for k = 1:n
            Rmat = diag(ones(1, h), 1);
            Rmat(h+1, h+1) = 1;
            R{v, k} = Rmat;
        end
    end
end
gamma = cache_gamma_lp(lambda,R);
% per-list hit probabilities are genuine only on the exact branch; see
% _kb/09-ldes-and-cache.md on cache-analyzer per-list/per-item reporting
pijlist = []; % genuine per-list occupancy (n x h), set in the exact branch
switch options.method
    case 'exact'
        % cache_prob_erec is exact only for the exchangeable (RR/FIFO/RANDOM)
        % family; see _kb/09-ldes-and-cache.md on cache-analyzer reporting
        switch sn.nodeparam{sn.nodetype == NodeType.Cache}.replacestrat
            case {ReplacementStrategy.RR, ReplacementStrategy.FIFO}
                % supported
            otherwise
                line_error(mfilename,'NC does not support exact solution of the specified cache replacement policy; use the default (approximate) method or SolverCTMC.');
        end
        line_debug('Using exact method, calling cache_prob_erec');
        [pij] = cache_prob_erec(gamma, m);
        missRate = zeros(1,u);
        for v=1:u
            missRate(v) = lambda(v,:,1)*pij(:,1);
        end
        pijlist = pij(:, 2:end);
        method='exact';
    case 'sampling'
        line_debug('Using sampling method, calling cache_miss_is');
        [~,missRate,~,~,lE] = cache_miss_is(gamma, m, lambda, options.samples);
        pij = cache_prob_is(gamma, m, options.samples);
        method='sampling';
    otherwise
        line_debug('Default method: using SPM approximation method\n');
        line_debug('Using SPM approximation method, calling cache_miss_spm');
        [~,missRate,~,~,lE] = cache_miss_spm(gamma, m, lambda);
        pij = cache_prob_spm(gamma, m, lE);
        method='spm';
end

for r = 1:sn.nclasses
    if length(ch.hitclass)>=r && ch.missclass(r)>0 && ch.hitclass(r)>0
        XN(ch.missclass(r)) = XN(ch.missclass(r)) + missRate(r);
        XN(ch.hitclass(r)) = XN(ch.hitclass(r)) + (sourceRate(r) - missRate(r));
    end
end

% per-list (per-level) hit probabilities (access-weighted), only where the
% exact algorithm produced a genuine per-list occupancy matrix.
hitproblist = NaN(u, h);
if ~isempty(pijlist)
    for v=1:u
        if any(~isnan(ch.pread{v}))
            pread = ch.pread{v}(:).';
            for l=1:h
                hitproblist(v,l) = pread * pijlist(:,l);
            end
        end
    end
end
% per-item occupancy [nitems x (lists+1)] (col 1 = miss) from the exact
% cache_prob_erec recursion (skipped, NaN, for >10 items); see
% _kb/09-ldes-and-cache.md on cache-analyzer per-list/per-item reporting
if n > 10
    line_warning(mfilename, 'Per-item cache occupancy (getAvgItemTable) requires the exact algorithm and is skipped for caches with more than 10 items (%d items); reporting NaN.', n);
    itemprob = NaN(n, h+1);
elseif ~isempty(pijlist)
    itemprob = pij; % exact branch: already [n x (h+1)] with col 1 = miss
else
    itemprob = cache_prob_erec(gamma, m); % [n x (h+1)] with col 1 = miss
end
runtime=toc(T0);
end
