function [QN,UN,RN,TN,CN,XN,lG,hitprob,missprob,runtime,it,method] = solver_nc_cacheqn_analyzer(self, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER] = SOLVER_NC_CACHEQN_ANALYZER(SELF, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

snorig = self.model.getStruct;
sn = snorig;
I = sn.nnodes;
K = sn.nclasses;

line_debug('NC cacheqn analyzer starting: method=%s, nclasses=%d', options.method, K);

statefulNodes = find(sn.isstateful)';
statefulNodesClasses = [];
for ind=statefulNodes %#ok<FXSET>
    statefulNodesClasses(end+1:end+K)= ((ind-1)*K+1):(ind*K);
end
lambda = zeros(1,K);
lambda_1 = zeros(1,K);
caches = find(sn.nodetype == NodeType.Cache);

hitprob = zeros(length(caches),K);
missprob = zeros(length(caches),K);

for it=1:options.iter_max
    for ind=caches
        ch = sn.nodeparam{ind};
        hitClass = ch.hitclass;
        missClass = ch.missclass;
        inputClass = find(hitClass);
        m = ch.itemcap;
        n = ch.nitems;
        if it == 1
            if n<sum(m)+2
                line_error(mfilename,'NC requires the number of items to exceed the cache capacity at least by 2.');
            end
            % initial random value of arrival rates to the cache
            lambda_1(inputClass) = rand(1,length(inputClass));
            lambda = lambda_1;
            sn.nodetype(ind) = NodeType.ClassSwitch;
        end

        % solution of isolated cache
        h = length(m);
        u = length(lambda);
        lambda_cache = zeros(u,n,h);

        for v=1:u
            for k=1:n
                for l=1:(h+1)
                    if ~isnan(ch.pread{v})
                        lambda_cache(v,k,l) = lambda(v) * ch.pread{v}(k);
                    end
                end
            end
        end

        Rcost = ch.accost;
        gamma = cache_gamma_lp(lambda_cache,Rcost);
        switch options.method
            case 'exact'
                [pij] = cache_prob_erec(gamma, m);
                missrate(ind,:) = zeros(1,u);
                for v=1:u
                    missrate(ind,v) = lambda_cache(v,:,1)*pij(:,1);
                end
                method = 'exact';
            otherwise
                line_debug('Default method: using SPM approximation for cache\n');
                line_debug('Using SPM approximation, calling cache_miss_spm');
                [~,missrate(ind,:)] = cache_miss_spm(gamma, m, lambda_cache);
                method = 'spm';
        end

        missprob(ind,:) = missrate(ind,:) ./ lambda; %  we set to NaN if no arrivals
        hitprob(ind,:) = 1 - missprob(ind,:);
        hitprob(isnan(hitprob)) = 0;
        missprob(isnan(missprob)) = 0;

        % bring back the isolated model results into the queueing model
        for r=inputClass
            sn.rtnodes((ind-1)*K+r,:) = 0;
            for jnd=1:I
                if sn.connmatrix(ind,jnd)
                    sn.rtnodes((ind-1)*K+r,(jnd-1)*K+hitClass(r)) = hitprob(ind,r);
                    sn.rtnodes((ind-1)*K+r,(jnd-1)*K+missClass(r)) = missprob(ind,r);
                end
            end
        end
        sn.rt = dtmc_stochcomp(sn.rtnodes,statefulNodesClasses);
    end
    [visits, nodevisits, sn] = sn_refresh_visits(sn, sn.chains, sn.rt, sn.rtnodes);
    sn.visits = visits;
    sn.nodevisits = nodevisits;

    if ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling)
        [QN,UN,RN,TN,CN,XN,lG,runtime] = solver_ncld_analyzer(sn, options);
    else
        [QN,UN,RN,TN,CN,XN,lG,runtime] = solver_nc_analyzer(sn, options);
    end

    nodevisits = cellsum(nodevisits);
    for ind=caches
        for r=inputClass
            c = find(sn.chains(:,r));
            inchain = find(sn.chains(c,:));
            if sn.refclass(c)>0
                lambda(r) = sum(XN(inchain)) * nodevisits(ind,r) / nodevisits(sn.stationToNode(sn.refstat(r)),sn.refclass(c));
            else
                lambda(r) = sum(XN(inchain)) * nodevisits(ind,r) / nodevisits(sn.stationToNode(sn.refstat(r)),r);
            end
        end
    end
    if norm(lambda-lambda_1,1) < options.iter_tol
        break
    end
    lambda_1 = lambda;
end
end
