function [QN,UN,RN,TN,CN,XN,lGN,runtime,iter,method,hitproblist,itemprob] = solver_mva_cache_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER,METHOD,HITPROBLIST] = SOLVER_MVA_CACHE_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0=tic;
QN = []; UN = [];
RN = []; TN = [];
CN = [];
XN = zeros(1,sn.nclasses);
lGN = NaN;
iter = NaN;

line_debug('MVA cache analyzer starting: method=%s, nclasses=%d', options.method, sn.nclasses);

source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
sourceRate = sn.rates(source_ist,:);
sourceRate(isnan(sourceRate)) = 0;
TN(source_ist,:) = sourceRate;

ch = sn.nodeparam{sn.nodetype == NodeType.Cache};

m = ch.itemcap;
n = ch.nitems;
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

Rcost = ch.accost;
if isempty(Rcost)
    % Default linear cache routing: items flow from list l to list l+1
    Rcost = cell(u, n);
    for v = 1:u
        for k = 1:n
            Rmat = diag(ones(1, h), 1);
            Rmat(h+1, h+1) = 1;
            Rcost{v, k} = Rmat;
        end
    end
end

gamma = cache_gamma_lp(lambda,Rcost);

% per-list hit probabilities are genuine only on the exact branch; see
% _kb/09-ldes-and-cache.md on cache-analyzer per-list/per-item reporting
pijlist = []; % genuine per-list occupancy (n x h), set in the exact branch
switch options.method
    case 'exact'
        line_debug('Using exact cache method');
        switch sn.nodeparam{sn.nodetype == NodeType.Cache}.replacestrat
            case {ReplacementStrategy.RR, ReplacementStrategy.FIFO}
                line_debug('Replacement strategy: RR/FIFO, calling cache_mva');
                [~,~,pij] = cache_mva(gamma, m);
                pij = [abs(1-sum(pij,2)),pij];
                pijlist = pij(:, 2:end);
            otherwise
                line_error(mfilename,'MVA does not support exact solution of the specified cache replacement policy.')
        end
    otherwise
        line_debug('Default method: using approximate cache method\n');
        line_debug('Using approximate cache method');
        switch sn.nodeparam{sn.nodetype == NodeType.Cache}.replacestrat
            case {ReplacementStrategy.RR, ReplacementStrategy.FIFO}
                line_debug('Replacement strategy: RR/FIFO, calling cache_prob_fpi');
                pij = cache_prob_fpi(gamma,m); % FPI method
            case ReplacementStrategy.LRU
                % Marked (MMAP) source is not IRM -> LRU(m)-MAP TTL approximation;
                % see _kb/09-ldes-and-cache.md on cache-analyzer reporting
                markedreaders = isfield(sn,'markidx') && ~isempty(sn.markidx) ...
                    && any(sn.markidx(source_ist,:) > 0);
                if markedreaders
                    Dcell = sn.proc{source_ist}{find(sn.markidx(source_ist,:)>0,1)};
                    D0 = Dcell{1}; D1agg = Dcell{2};
                    D0c = cell(1,n); D1c = cell(1,n);
                    allmarked = true;
                    for k=1:n
                        D1c{k} = zeros(size(D0));
                        for v=1:u
                            if ~isnan(ch.pread{v})
                                if sn.markidx(source_ist,v) > 0
                                    D1c{k} = D1c{k} + Dcell{2+sn.markidx(source_ist,v)} * ch.pread{v}(k);
                                elseif sourceRate(v) > 0 && any(ch.pread{v} > 0)
                                    allmarked = false; % unmarked reader mixed in
                                end
                            end
                        end
                        D0c{k} = D0 + D1agg - D1c{k};
                    end
                    if allmarked
                        line_debug('Replacement strategy: LRU with marked MAP source, calling cache_ttl_lrum_map');
                        pij = cache_ttl_lrum_map(D0c, D1c, m);
                    else
                        line_debug('Replacement strategy: LRU, calling cache_ttl_lrua');
                        pij = cache_ttl_lrua(lambda, Rcost, m);  % allows trees and access costs
                    end
                else
                    line_debug('Replacement strategy: LRU, calling cache_ttl_lrua');
                    pij = cache_ttl_lrua(lambda, Rcost, m);  % allows trees and access costs
                end
            case ReplacementStrategy.HLRU
                % h-LRU / LRU(m) characteristic-time approximation (linear
                % list topology; access-cost graphs are not supported)
                line_debug('Replacement strategy: HLRU, calling cache_ttl_hlru');
                pij = cache_ttl_hlru(lambda, m);
            otherwise
                line_error(mfilename,'MVA does not support approximate solution of the specified cache replacement policy.')
        end
end
missRate = zeros(1,u);
for v=1:u
    missRate(v) = lambda(v,:,1)*pij(:,1);
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

% per-item occupancy [nitems x (lists+1)] (col 1 = miss); derived from the exact
% cache_mva recursion (skipped, NaN, for >10 items); see
% _kb/09-ldes-and-cache.md on cache-analyzer per-list/per-item reporting
itemprob = [];
if ~isempty(pijlist)
    itemprob = pij; % exact branch: already [n x (h+1)] with col 1 = miss
elseif size(pij,2) == h+1
    switch sn.nodeparam{sn.nodetype == NodeType.Cache}.replacestrat
        case {ReplacementStrategy.RR, ReplacementStrategy.FIFO}
            if n > 10
                line_warning(mfilename, 'Per-item cache occupancy (getAvgItemTable) requires the exact algorithm for RR/FIFO and is skipped for caches with more than 10 items (%d items); reporting NaN.', n);
                itemprob = NaN(n, h+1);
            else
                [~,~,pij_ex] = cache_mva(gamma, m);
                itemprob = [abs(1-sum(pij_ex,2)), pij_ex];
            end
        otherwise
            itemprob = pij; % LRU-TTL: genuine per-list distribution
    end
end

for r = 1:sn.nclasses
    if length(ch.hitclass)>=r && ch.missclass(r)>0 && ch.hitclass(r)>0
        XN(ch.missclass(r)) = XN(ch.missclass(r)) + missRate(r);
        XN(ch.hitclass(r)) = XN(ch.hitclass(r)) + (sourceRate(r) - missRate(r));
    end
end

% Set the actual method used
if strcmp(options.method, 'exact')
    method = 'exact';
else
    switch sn.nodeparam{sn.nodetype == NodeType.Cache}.replacestrat
        case {ReplacementStrategy.RR, ReplacementStrategy.FIFO}
            method = 'fpi';
        case ReplacementStrategy.LRU
            method = 'ttl';
        otherwise
            method = options.method;
    end
end

runtime=toc(T0);
end
