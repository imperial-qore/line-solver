function [chains, visits, rt, nodevisits, rtnodes] = refreshChains(self, propagate)
% [CHAINS, VISITS, RT, NODEVISITS, RTNODES] = REFRESHCHAINS(PROPAGATE)
%
% PROPAGATE : true if the change needs to trigger recomputation of visits and station capacities (default: true)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<2 || nargin<3
    propagate = true;
end

%% obtain routing matrix
sn = self.sn;
rates = sn.rates;
[rt,~,rtnodes] = self.refreshRoutingMatrix(rates);
sn = self.sn;

%% determine class switching mask
stateful = find(sn.isstateful); stateful = stateful(:)';
K = sn.nclasses;
if isempty(self.csMatrix) % only for models created without a call to link()
    csmask = false(K,K);
    for r=1:K
        for s=1:K
            for isf=1:length(stateful) % source
                for jsf=1:length(stateful) % source
                    if rt((isf-1)*K+r, (jsf-1)*K+s) > 0
                        % this is to ensure that we use rt, which is
                        % the stochastic complement taken over the stateful
                        % nodes, otherwise sequences of cs can produce a wrong
                        % csmask
                        csmask(r,s) = true;
                    end
                end
            end
        end
    end
    
    for isf=1:length(stateful) % source
        % this is to ensure that also stateful cs like caches
        % are accounted
        ind = sn.statefulToNode(isf);
        isCS = sn.nodetype(ind) == NodeType.Cache | sn.nodetype(ind) == NodeType.ClassSwitch;
        for r=1:K
            csmask(r,r) = true;
            for s=1:K
                if r~=s
                    if isCS
                        if self.nodes{ind}.server.csFun(r,s,[],[])>0
                            csmask(r,s) = true;
                        end
                    end
                end
            end
        end
    end
    sn.csmask = csmask;
else
    % Start with csMatrix from link(), but also include any transitions
    % that may have been created by stochastic complement through
    % non-stateful nodes (e.g., auto-added ClassSwitch nodes)
    csmask = self.csMatrix;
    for r=1:K
        for s=1:K
            for isf=1:length(stateful)
                for jsf=1:length(stateful)
                    if rt((isf-1)*K+r, (jsf-1)*K+s) > 0
                        csmask(r,s) = true;
                    end
                end
            end
        end
    end
    sn.csmask = csmask;
end

if isfield(sn,'refclass') && length(sn.refclass)<sn.nchains
    % if the number of chains changed dynamically, extend refclass
    sn.refclass(end+1:sn.nchains) = 0; 
end

self.sn = sn;

%% compute visits
if propagate
    [visits, nodevisits, sn] = sn_refresh_visits(sn, sn.chains, rt, rtnodes);
end

self.sn = sn;

%% call dependent capacity refresh
if propagate
    refreshCapacity(self); % capacity depends on chains and rates
end

end
