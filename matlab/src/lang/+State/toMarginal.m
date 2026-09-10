function [ni, nir, sir, kir] = toMarginal(sn, ind, state_i, phasesz, phaseshift, space_buf, space_srv, space_var) %#ok<INUSD>
% TOMARGINAL Compute marginal distributions for a specific node
%
% [NI, NIR, SIR, KIR] = TOMARGINAL(SN, IND, STATE_I, PHASESZ, PHASESHIFT, SPACE_BUF, SPACE_SRV, SPACE_VAR)
%
% @brief Extracts marginal job statistics from global state information for a specific node
%
% This function processes the global network state to compute marginal 
% distributions and job counts for a specific node, considering different
% scheduling strategies, service phases, and node types.
%
% @param sn Network structure or Network object
% @param ind Node index to extract marginal information for
% @param state_i Global state matrix or vector
% @param phasesz Vector of phase sizes for each class
% @param phaseshift Phase shift parameters for state extraction
% @param space_buf Buffer space configuration
% @param space_srv Service space configuration  
% @param space_var Variable space configuration
%
% @return ni Total jobs in node IND
% @return nir Total jobs per class in node IND [matrix: states x classes]
% @return sir Total jobs in service per class in node IND [matrix: states x classes]
% @return kir Total jobs in service per class and per phase in node IND [3D array: states x classes x phases]

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isstruct(sn) % the input can be a Network object too
    sn=sn.getStruct();
end

% Join stations on FJ-augmented structs: plain per-class count state,
% no buffer/phase split and no service (jobs wait for synchronization)
if isfield(sn,'isfjaugmented') && sn.isfjaugmented && sn.nodetype(ind) == NodeType.Join
    R = sn.nclasses;
    if nargin < 3
        state_i = sn.state{sn.nodeToStateful(ind)};
    end
    nir = state_i(:,(end-R+1):end);
    sir = zeros(size(nir));
    kir = zeros(size(state_i,1),R,1);
    ni = sum(nir,2);
    return
end

% ind: node index
if ~sn.isstation(ind) && sn.isstateful(ind) % if stateful node
    if sn.nodetype(ind) == NodeType.Transition
        R = sn.nodeparam{ind}.nmodes;
        if nargin < 5
            % Initialize phasesz/phaseshift for Transition nodes from firingphases
            phasesz = sn.nodeparam{ind}.firingphases;
            if any(isnan(phasesz))
                phasesz = zeros(1, R);
                for m = 1:R
                    if iscell(sn.nodeparam{ind}.firingproc) && ~isempty(sn.nodeparam{ind}.firingproc{m})
                        phasesz(m) = size(sn.nodeparam{ind}.firingproc{m}{1}, 1);
                    else
                        phasesz(m) = 1;
                    end
                end
            end
            phaseshift = [0,cumsum(phasesz(1:end-1))];
        end
        if nargin < 8
            nmodes = R;
            space_buf = state_i(:,1:nmodes);
            space_srv = state_i(:,(nmodes+1):(nmodes+sum(phasesz)));
            space_var = state_i(:,(nmodes+sum(phasesz)+1):end);
        end
        sir = zeros(size(state_i,1),R); % class-r jobs in service
        kir = zeros(size(state_i,1),R,max(phasesz)); % class-r jobs in service in phase k
        for r=1:R
            for k=1:phasesz(r)
                kir(:,r,k) = space_srv(:,phaseshift(r)+k);
                sir(:,r) = sir(:,r) + kir(:,r,k);
            end
        end
        nir = sir;
        ni = sum(nir);
    else
        ni = sum(state_i(1:(end-sum(sn.nvars(ind,:)))));
        nir = state_i(1:(end-sum(sn.nvars(ind,:))));
        sir = nir; % jobs in service
        kir = sir; % jobs per phase
    end
    return
end

R = sn.nclasses;

ist = sn.nodeToStation(ind);
isf = sn.nodeToStateful(ind);

if nargin < 3
    state_i = sn.state{isf};
end

if nargin < 5
    phasesz = sn.phasessz(ist,:);
    phaseshift = sn.phaseshift(ist,:);
end

if sn.sched(ist) == SchedStrategy.PAS
    % PAS/OI: the state is the ordered list of class indices (no server split).
    Vp = sum(sn.nvars(ind,:));
    listcols = state_i(:,1:(end-Vp));
    nir = zeros(size(state_i,1),R);
    for r=1:R
        nir(:,r) = sum(listcols==r,2);
    end
    % Jobs in service = positions receiving a nonzero rate increment
    % Delta_mu(c1..cp) > 0 (the rest wait); utilization derives from this.
    muFun = sn.nodeparam{ind}.svcRateFun;
    sir = zeros(size(state_i,1),R);
    if ~isempty(muFun)
        for row=1:size(listcols,1)
            c = listcols(row, listcols(row,:)>0);
            muPrev = 0;
            for p=1:numel(c)
                muCur = muFun(c(1:p));
                if muCur - muPrev > 0
                    sir(row,c(p)) = sir(row,c(p)) + 1;
                end
                muPrev = muCur;
            end
        end
    else
        sir = nir;
    end
    kir = zeros(size(state_i,1),R,max(phasesz));
    kir(:,:,1) = sir;
    ni = sum(nir,2);
    return
end

isExponential = false;
if max(phasesz)==1
    isExponential = true;
end

if nargin < 8
    space_var = state_i(:,(end-sum(sn.nvars(ind,:))+1):end); % local variables
    space_srv = state_i(:,(end-sum(phasesz)-sum(sn.nvars(ind,:))+1):(end-sum(sn.nvars(ind,:)))); % server state
    space_buf = state_i(:,1:(end-sum(phasesz)-sum(sn.nvars(ind,:)))); % buffer state
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

switch sn.sched(ist)
    case SchedStrategy.INF
        for r=1:R
            nir(:,r) = sir(:,r); % class-r jobs in station
        end
    case {SchedStrategy.PS, SchedStrategy.PSPRIO}
        for r=1:R
            nir(:,r) = sir(:,r) ; % class-r jobs in station
        end
    case SchedStrategy.EXT
        for r=1:R
            nir(:,r) = Inf;
        end
    case {SchedStrategy.FCFS, SchedStrategy.FCFSPRIO, SchedStrategy.HOL}
        for r=1:R
            nir(:,r) = sir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case {SchedStrategy.DPS, SchedStrategy.DPSPRIO}
        for r=1:R
            nir(:,r) = sir(:,r) ; % class-r jobs in station
        end
    case {SchedStrategy.GPS, SchedStrategy.GPSPRIO}
        for r=1:R
            nir(:,r) = sir(:,r) ; % class-r jobs in station
        end
    case {SchedStrategy.LCFS, SchedStrategy.LCFSPRIO}
        for r=1:R
            nir(:,r) = sir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case {SchedStrategy.FCFSPI, SchedStrategy.FCFSPIPRIO, SchedStrategy.FCFSPR, SchedStrategy.FCFSPRPRIO, SchedStrategy.LCFSPI, SchedStrategy.LCFSPIPRIO, SchedStrategy.LCFSPR, SchedStrategy.LCFSPRPRIO}
        % buffer holds (class,phase) pairs: slice the class columns per row.
        % length()+linear indexing would flatten a multi-row state matrix.
        if size(space_buf,2)>1
            space_bufcls = space_buf(:,1:2:end);
            %space_bufphase = space_buf(:,2:2:end);
            for r=1:R
                nir(:,r) = sir(:,r) + sum(space_bufcls==r,2); % class-r jobs in station
            end
        else
            nir = sir;
        end
    case {SchedStrategy.SIRO, SchedStrategy.POLLING}
        for r=1:R
            nir(:,r) = sir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    case SchedStrategy.SEPT
        for r=1:R
            nir(:,r) = sir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    case SchedStrategy.LEPT
        for r=1:R
            nir(:,r) = sir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    case SchedStrategy.SRPT
        for r=1:R
            nir(:,r) = sir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    case SchedStrategy.SRPTPRIO
        for r=1:R
            nir(:,r) = sir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    otherwise % possibly other stateful nodes
        for r=1:R
            nir(:,r) = sir(:,r) ; % class-r jobs in station
        end
end

if sn.nodetype(ind) == NodeType.Place
    % A Place is a token container: the [buffer, server] split in its state
    % encoding is an artifact (a transition FIRE relocates the surviving
    % tokens from the server slot into the buffer slot). For the INF-family
    % scheduling used by Places the switch above set nir = sir (server slot
    % only), so fold the buffer slot back in. This makes the marginal token
    % count agree with State.toMarginalAggr (buf + srv), which drives the
    % transition enabling degree; without it the QLen measurement collapses
    % to the oscillating server slot while the dynamics stay correct.
    switch sn.sched(ist)
        case {SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.PSPRIO, ...
              SchedStrategy.DPS, SchedStrategy.DPSPRIO, ...
              SchedStrategy.GPS, SchedStrategy.GPSPRIO, SchedStrategy.LPS}
            if ~isempty(space_buf) && size(space_buf,2) >= R
                for r=1:R
                    nir(:,r) = sir(:,r) + space_buf(:,r);
                end
            end
    end
    % set all service and phase data to 0
    for r=1:R
        if isnan(sn.rates(ist,r)) % if disabled station
            for k=1:phasesz(r)
                kir(:,r,k) = 0;
            end
            sir(:,r)=0;
        end
    end
else
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
