function [QN,UN,RN,TN,CN,XN,lG,pij,runtime,method,hitproblist,itemprob,listcost] = solver_nc_cache_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,PIJ,RUNTIME,METHOD,HITPROBLIST,ITEMPROB,LISTCOST] = SOLVER_NC_CACHE_ANALYZER(QN, OPTIONS)

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
[gamma,~,~,~,parent] = cache_gamma_lp(lambda,R);

% per-item storage costs and per-list cost caps (ton21cache Sec. IX)
sigma = []; costcap = [];
if isfield(ch,'itemsize'), sigma = ch.itemsize(:).'; end
if isfield(ch,'costcap'), costcap = ch.costcap(:).'; end
if ~isempty(costcap)
    if isempty(sigma)
        line_error(mfilename,'Storage cost caps require per-item sizes; call Cache.setItemSizes first.');
    end
    if numel(sigma)~=n
        line_error(mfilename,'The item size vector must have one entry per item.');
    end
    if numel(costcap)~=h
        line_error(mfilename,'The cost cap vector must have one entry per cache list.');
    end
    viol = cache_cost_pathcheck(gamma, sigma, costcap, parent);
    if ~isempty(viol)
        line_warning(mfilename, 'Storage cost caps block the promotion path of item %d into list %d at list %d (and %d further pairs). The exact recursion normalizes over all size-feasible states, which is then a strict superset of the states the cache can reach; cross-check with SolverLDES.', viol(1,1), viol(1,2), viol(1,3), size(viol,1)-1);
    end
end
% sizes without caps place no constraint on the state space, but still feed
% the mean per-list storage cost reported below
% per-list hit probabilities are genuine only on the exact branch; see
% _kb/09-ldes-and-cache.md on cache-analyzer per-list/per-item reporting
pijlist = []; % genuine per-list occupancy (n x h), set in the exact branch
cacheMethod = options.method;
% 'rayint' is an alias of 'spm': on a cache both name the SPM saddle point, and
% the method name stays live for solver_nc_retrieval_analyzer's delayed-hit expansion.
if any(strcmp(cacheMethod,{'rayint','spm'}))
    cacheMethod = 'default';
end
% The SPM family serves its size-tilted form (cache_spm_size) once the items
% carry storage costs. The saddle escapes to infinity at sum(m) = n, so the
% size-free saddle point takes over there rather than the exact recursion,
% which would refuse every replacement policy outside RR/FIFO.
useSpmSize = strcmp(cacheMethod,'default') && ~isempty(sigma) && sum(m) < n;
if ~isempty(costcap) && ~useSpmSize && ~any(strcmp(cacheMethod,{'exact','sampling'}))
    % The size-free SPM, FPI and mean-field methods have no cost-capped
    % counterpart: the k - sigma_i 1_j argument couples item sizes into the
    % recursion graph, which only cache_spm_size and the exact path carry.
    if n*prod(m+1)*prod(costcap+1) <= 1e6
        cacheMethod = 'exact';
    else
        cacheMethod = 'sampling';
    end
    line_warning(mfilename,'Method ''%s'' does not support storage cost caps; using ''%s'' instead.', options.method, cacheMethod);
end
switch cacheMethod
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
        [pij] = cache_prob_erec(gamma, m, sigma, costcap);
        missRate = zeros(1,u);
        for v=1:u
            missRate(v) = lambda(v,:,1)*pij(:,1);
        end
        pijlist = pij(:, 2:end);
        method='exact';
    case 'sampling'
        line_debug('Using sampling method, calling cache_miss_is');
        [~,missRate,~,~,lE] = cache_miss_is(gamma, m, lambda, options.samples, sigma, costcap);
        pij = cache_prob_is(gamma, m, options.samples, sigma, costcap);
        method='sampling';
    otherwise
        if useSpmSize
            % Size-tilted SPM: a 2h Newton solve whose cost does not grow with
            % the (m,k) lattice the exact recursion walks. O(1/n), so it wants
            % room between the occupancies and n.
            line_debug('Using size-tilted SPM expansion, calling cache_spm_size');
            raycap = costcap;
            if isempty(raycap)
                % Sizes but no caps: cap each list at the dearest load it can
                % hold, which is exactly slack, so the cost coordinate leaves
                % the saddle and the expansion degenerates to the size-free one.
                srt = sort(sigma,'descend');
                raycap = zeros(1,h);
                for j=1:h
                    raycap(j) = sum(srt(1:m(j)));
                end
            end
            [~,lE,rayout] = cache_spm_size(gamma, m, sigma, raycap);
            pij = rayout.pij;
            missRate = zeros(1,u);
            for v=1:u
                missRate(v) = lambda(v,:,1)*pij(:,1);
            end
            method='spm.size';
        else
            line_debug('Default method: using SPM approximation method\n');
            line_debug('Using SPM approximation method, calling cache_miss_spm');
            [~,missRate,~,~,lE] = cache_miss_spm(gamma, m, lambda);
            pij = cache_prob_spm(gamma, m, lE);
            method='spm';
        end
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
    itemprob = cache_prob_erec(gamma, m, sigma, costcap); % [n x (h+1)] with col 1 = miss
end
% mean storage cost held by each list, K_j = sum_i sigma_i pi_ij
listcost = NaN(1,h);
if ~isempty(sigma) && size(pij,1)==n && size(pij,2)==h+1
    listcost = cache_cost(gamma, m, sigma, costcap, pij);
end
runtime=toc(T0);
end
