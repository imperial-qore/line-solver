function [Q,U,R,T,C,X,lG,runtime,iter] = solver_mva_cache_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER] = SOLVER_MVA_CACHE_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

T0=tic;
Q = []; U = [];
R = []; T = [];
C = [];
X = zeros(1,sn.nclasses);
lG = NaN;
iter = NaN;

source_ist = sn.nodeToStation(sn.nodetype == NodeType.Source);
sourceRate = sn.rates(source_ist,:);
sourceRate(isnan(sourceRate)) = 0;
T(source_ist,:) = sourceRate;

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
        [~,~,pij] = cache_mva(gamma, m);
        pij = [abs(1-sum(pij,2)),pij];
        missRate = zeros(1,u);
        for v=1:u
            missRate(v) = lambda(v,:,1)*pij(:,1);
        end
    case 'ttl.lrum'
        pij = cache_ttl_lrum(lambda, m);  % without considering different graph of different items  linear
        missRate = zeros(1,u);
        for v=1:u
            missRate(v) = lambda(v,:,1)*pij(:,1);
        end        
    case 'ttl.hlru'
        pij = cache_ttl_hlru(lambda, m);  % without considering different graph of different items linear
        missRate = zeros(1,u);
        for v=1:u
            missRate(v) = lambda(v,:,1)*pij(:,1);
        end          
    case 'ttl.tree'
        pij = cache_ttl_tree(lambda, Rcost, m);  % considering different graphs of different items
        missRate = zeros(1,u);
        for v=1:u
            missRate(v) = lambda(v,:,1)*pij(:,1);
        end         
    otherwise
        pij = cache_prob_asy(gamma,m); % FPI method
        missRate = zeros(1,u);
        for v=1:u
            missRate(v) = lambda(v,:,1)*pij(:,1);
        end
end

for r = 1:sn.nclasses
    if length(ch.hitclass)>=r && ch.missclass(r)>0 && ch.hitclass(r)>0
        X(ch.missclass(r)) = X(ch.missclass(r)) + missRate(r);
        X(ch.hitclass(r)) = X(ch.hitclass(r)) + (sourceRate(r) - missRate(r));
    end
end
runtime=toc(T0);
end
