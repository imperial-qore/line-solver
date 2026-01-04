function [QN,UN,RN,TN,CN,XN,lGN,runtime,iter,method] = solver_mva_cache_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER] = SOLVER_MVA_CACHE_ANALYZER(QN, OPTIONS)

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

gamma = cache_gamma_lp(lambda,Rcost);

switch options.method
    case 'exact'
        line_debug('Using exact cache method');
        switch sn.nodeparam{sn.nodetype == NodeType.Cache}.replacestrat
            case {ReplacementStrategy.RR, ReplacementStrategy.FIFO}
                line_debug('Replacement strategy: RR/FIFO, calling cache_mva');
                [~,~,pij] = cache_mva(gamma, m);
                pij = [abs(1-sum(pij,2)),pij];
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
                line_debug('Replacement strategy: LRU, calling cache_ttl_lrua');
                %pij = cache_ttl_lrum(lambda, m);  % linear topology only
                pij = cache_ttl_lrua(lambda, Rcost, m);  % allows trees and access costs
            %case ReplacementStrategy.HLRU
                % this is integrated but the cache_ttl_hlru script is not
                % working correctly at present
                %pij = cache_ttl_hlru(lambda, m);  % without considering different graph of different items linear
            otherwise
                line_error(mfilename,'MVA does not support approximate solution of the specified cache replacement policy.')
        end
end
missRate = zeros(1,u);
for v=1:u
    missRate(v) = lambda(v,:,1)*pij(:,1);
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
