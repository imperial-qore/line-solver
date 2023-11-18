function [ni, nir] = toMarginalAggr(sn, ind, state_i, K, Ks, space_buf, space_srv, space_var) %#ok<INUSD>
% [NI, NIR] = TOMARGINALAGGR(QN, IND, STATE_I, K, KS, SPACE_BUF, SPACE_SRV, SPACE_VAR) %#OK<INUSD>

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

if ~isstruct(sn) % the input can be a Network object too
    sn=sn.getStruct();
end
% ind: node index
ist = sn.nodeToStation(ind);
%isf = sn.nodeToStateful(ind);
R = sn.nclasses;
if ~sn.isstation(ind) && sn.isstateful(ind) % if stateful node
    ni = sum(state_i(1:(end-sum(sn.nvars(ind,:)))));
    nir = state_i(1:(end-sum(sn.nvars(ind,:))));
    return
end

if nargin < 5
    K = sn.phasessz(ist,:);
    Ks = sn.phaseshift(ist,:);
end

if nargin < 8
    space_var = state_i(:,(end-sum(sn.nvars(ind,:))+1):end); % server state
    space_srv = state_i(:,(end-sum(K)+1):(end-sum(sn.nvars(ind,:))));
    space_buf = state_i(:,1:(end-sum(K)));
end

nir = zeros(size(state_i,1),R); % class-r jobs in service
for r=1:R
    for k=1:K(r)
        nir(:,r) = nir(:,r) + space_srv(:,Ks(r)+k);
    end
end
switch sn.schedid(ist)
    case SchedStrategy.ID_EXT
        for r=1:R
            nir(:,r) = Inf;
        end
    case SchedStrategy.ID_FCFS
        for r=1:R
            nir(:,r) = nir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case SchedStrategy.ID_HOL
        for r=1:R
            nir(:,r) = nir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case SchedStrategy.ID_LCFS
        for r=1:R
            nir(:,r) = nir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case SchedStrategy.ID_SIRO
        for r=1:R
            nir(:,r) = nir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    case SchedStrategy.ID_SEPT
        for r=1:R
            nir(:,r) = nir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    case SchedStrategy.ID_LEPT
        for r=1:R
            nir(:,r) = nir(:,r) + space_buf(:,r); % class-r jobs in station
        end
        %otherwise % possibly other stateful nodes
        % no-op
end

for r=1:R
    if isnan(sn.rates(ist,r)) % if disabled
        nir(:,r) = 0;
    end
end

ni = sum(nir,2); % total jobs in station
end
