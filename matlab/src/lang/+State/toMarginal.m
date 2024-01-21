function [ni, nir, sir, kir] = toMarginal(sn, ind, state_i, phasesz, phaseshift, space_buf, space_srv, space_var) %#ok<INUSD>
% [NI, NIR, SIR, KIR] = TOMARGINAL(QN, IND, STATE_I, K, KS, SPACE_BUF, SPACE_SRV, SPACE_VAR) %#OK<INUSD>
% NI: total jobs in node IND
% NIR: total jobs per class in node IND
% SIR: total jobs in service per class in node IND
% KIR: totak jobs in service per class and per phase in node IND

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

if ~isstruct(sn) % the input can be a Network object too
    sn=sn.getStruct();
end

% ind: node index
if ~sn.isstation(ind) && sn.isstateful(ind) % if stateful node
    ni = sum(state_i(1:(end-sum(sn.nvars(ind,:)))));
    nir = state_i(1:(end-sum(sn.nvars(ind,:))));
    sir = nir; % jobs in service
    kir = sir; % jobs per phase
    return
end

R = sn.nclasses;
ist = sn.nodeToStation(ind);
%isf = sn.nodeToStateful(ind);

if nargin < 5
    phasesz = sn.phasessz(ist,:);
    phaseshift = sn.phaseshift(ist,:);
end

isExponential = false;
if max(phasesz)==1
    isExponential = true;
end

if nargin < 8
    space_var = state_i(:,(end-sum(sn.nvars(ind,:))+1):end); % server state
    space_srv = state_i(:,(end-sum(phasesz)-sum(sn.nvars(ind,:))+1):(end-sum(sn.nvars(ind,:))));
    space_buf = state_i(:,1:(end-sum(phasesz)-sum(sn.nvars(ind,:))));
end

if isExponential
    sir = space_srv;
    kir = space_srv;
else
    nir = zeros(size(state_i,1),R);
    sir = zeros(size(state_i,1),R); % class-r jobs in service
    kir = zeros(size(state_i,1),R,max(phasesz)); % class-r jobs in service in phase k
    for r=1:R
        for k=1:phasesz(r)
            kir(:,r,k) = space_srv(:,phaseshift(r)+k);
            sir(:,r) = sir(:,r) + kir(:,r,k);
        end
    end
end
switch sn.schedid(ist)
    case SchedStrategy.ID_INF
        for r=1:R
            nir(:,r) = sir(:,r); % class-r jobs in station
        end
    case SchedStrategy.ID_PS
        for r=1:R
            nir(:,r) = sir(:,r) ; % class-r jobs in station
        end
    case SchedStrategy.ID_EXT
        for r=1:R
            nir(:,r) = Inf;
        end
    case SchedStrategy.ID_FCFS
        for r=1:R
            nir(:,r) = sir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case SchedStrategy.ID_DPS
        for r=1:R
            nir(:,r) = sir(:,r) ; % class-r jobs in station
        end
    case SchedStrategy.ID_GPS
        for r=1:R
            nir(:,r) = sir(:,r) ; % class-r jobs in station
        end
    case SchedStrategy.ID_HOL
        for r=1:R
            nir(:,r) = sir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case SchedStrategy.ID_LCFS
        for r=1:R
            nir(:,r) = sir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case SchedStrategy.ID_LCFSPR
        if length(space_buf)>1
            space_buf = space_buf(1:2:end);
            %space_bufphase = space_buf(2:2:end);            
            for r=1:R
                nir(:,r) = sir(:,r) + sum(space_buf==r,2); % class-r jobs in station
            end
        else
            nir = sir;
        end
    case SchedStrategy.ID_SIRO
        for r=1:R
            nir(:,r) = sir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    case SchedStrategy.ID_SEPT
        for r=1:R
            nir(:,r) = sir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    case SchedStrategy.ID_LEPT
        for r=1:R
            nir(:,r) = sir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    otherwise % possibly other stateful nodes
        for r=1:R
            nir(:,r) = sir(:,r) ; % class-r jobs in station
        end
end

if sn.nodetype(ind) ~= NodeType.Place
    for r=1:R
        if isnan(sn.rates(ist,r)) % if disabled station
            nir(:,r) = 0;
            for k=1:phasesz(r)
                kir(:,r,k) = 0;
            end
            sir(:,r)=0;
        end
    end
end

ni = sum(nir,2); % total jobs in station
end
