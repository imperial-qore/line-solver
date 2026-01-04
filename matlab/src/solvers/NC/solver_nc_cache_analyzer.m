function [QN,UN,RN,TN,CN,XN,lG,pij,runtime,method] = solver_nc_cache_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,PIJ,RUNTIME] = SOLVER_NC_CACHE_ANALYZER(QN, OPTIONS)

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
gamma = cache_gamma_lp(lambda,R);
switch options.method
    case 'exact'
        line_debug('Using exact method, calling cache_prob_erec');
        [pij] = cache_prob_erec(gamma, m);
        missRate = zeros(1,u);
        for v=1:u
            missRate(v) = lambda(v,:,1)*pij(:,1);
        end
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
runtime=toc(T0);
end
