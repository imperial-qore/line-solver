function [ni, nir] = toMarginalAggr(sn, ind, state_i, K, Ks, space_buf, space_srv, space_var) %#ok<INUSD>
% TOMARGINALAGGR Compute aggregate marginal distributions for a specific node
%
% [NI, NIR] = TOMARGINALAGGR(SN, IND, STATE_I, K, KS, SPACE_BUF, SPACE_SRV, SPACE_VAR)
%
% @brief Computes aggregate marginal job counts from state information without phase details
%
% This function provides a simplified version of toMarginal that computes
% aggregate job counts per node and per class, without considering individual
% service phases. It is more efficient when phase-level detail is not required.
%
% @param sn Network structure or Network object
% @param ind Node index to extract marginal information for
% @param state_i Global state matrix or vector
% @param K Vector of population for each class
% @param Ks Matrix of populations per chain and class
% @param space_buf Buffer space configuration
% @param space_srv Service space configuration  
% @param space_var Variable space configuration
%
% @return ni Total jobs in node IND (aggregate across all classes)
% @return nir Jobs per class in node IND [vector: classes]

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isstruct(sn) % the input can be a Network object too
    sn=sn.getStruct();
end
% ind: node index
ist = sn.nodeToStation(ind);
%isf = sn.nodeToStateful(ind);
R = sn.nclasses;

% Join stations on FJ-augmented structs: plain per-class count state
if isfield(sn,'isfjaugmented') && sn.isfjaugmented && sn.nodetype(ind) == NodeType.Join
    nir = state_i(:,(end-R+1):end);
    ni = sum(nir,2);
    return
end

if ~sn.isstation(ind) && sn.isstateful(ind) % if stateful node
    % Strip trailing nvars columns (e.g. RROBIN pointer on a Router) and
    % keep the buffer portion. If the buffer is empty (the node has no
    % per-class job buffer at all, e.g. a Router only carrying nvars
    % bookkeeping), return zeros of the expected (1, R) shape rather than
    % an empty slice, so callers can index nir(r) for r = 1..R.
    nvarsSum = sum(sn.nvars(ind,:));
    bufferLen = size(state_i, 2) - nvarsSum;
    if bufferLen <= 0
        ni  = 0;
        nir = zeros(1, R);
    else
        ni  = sum(state_i(1:bufferLen));
        nir = state_i(1:bufferLen);
        if numel(nir) < R
            nir(end+1:R) = 0;
        end
    end
    return
end

% PAS/OI stations: the state is the ordered list of class indices (no server split).
if sn.isstation(ind) && sn.sched(ist) == SchedStrategy.PAS
    Vp = sum(sn.nvars(ind,:));
    listcols = state_i(:,1:(end-Vp));
    nir = zeros(size(state_i,1),R);
    for r=1:R
        nir(:,r) = sum(listcols==r,2);
    end
    ni = sum(nir,2);
    return
end

% Place nodes: compute total jobs per class from state
if sn.nodetype(ind) == NodeType.Place
    state_len = size(state_i, 2);
    % Check if this is queue-based state (FCFS/LCFS/HOL) vs count-based
    % Queue-based: state contains class indices [1,1,2,1] meaning job arrivals
    % Count-based: state contains counts per class [n1, n2, ...] or [buf(R), srv]
    % Detect queue-based: length != R and length != 2*R and all values are valid class indices
    is_queue_based = state_len ~= R && state_len ~= 2*R;
    if is_queue_based && state_len > 0
        % Check all values are valid class indices (1 to R)
        all_vals = state_i(:);
        is_queue_based = all(all_vals >= 1 & all_vals <= R);
    end

    if is_queue_based
        % Queue-based format: count occurrences of each class
        nir = zeros(size(state_i,1), R);
        for r = 1:R
            nir(:,r) = sum(state_i == r, 2);
        end
    elseif nargin >= 4 && ~isempty(K)
        expected_len = R + sum(K);
        if state_len == expected_len
            % State has [buffer, server] format
            buf_part = state_i(:, 1:R);
            srv_part = state_i(:, (R+1):end);
            nir = buf_part;
            for r = 1:R
                for k = 1:K(r)
                    nir(:,r) = nir(:,r) + srv_part(:, Ks(r)+k);
                end
            end
        elseif state_len == R
            % State is just counts per class (simple SIRO format)
            nir = state_i;
        else
            % Unknown format - just take first R columns
            nir = state_i(:, 1:min(R, state_len));
            if size(nir, 2) < R
                nir = [nir, zeros(size(nir,1), R - size(nir,2))];
            end
        end
    else
        % No K provided - assume state is counts per class
        nir = state_i(:, 1:min(R, size(state_i, 2)));
        if size(nir, 2) < R
            nir = [nir, zeros(size(nir,1), R - size(nir,2))];
        end
    end
    ni = sum(nir,2);
    return
end

% Source nodes have infinite population but report queue length 0 (jobs are external)
if sn.nodetype(ind) == NodeType.Source
    nir = zeros(size(state_i,1), R);
    ni = zeros(size(state_i,1), 1);
    return
end

if nargin < 4
    K = sn.phasessz(ist,:);
end
    if nargin < 5
    Ks = sn.phaseshift(ist,:);
end

if nargin < 8
    % The local variables trail the server block, so both the server and the
    % buffer slices must be taken clear of them (mirrors State.toMarginal and
    % the slicing State.afterEvent hands to the explicit-argument callers).
    Vagg = sum(sn.nvars(ind,:));
    space_var = state_i(:,(end-Vagg+1):end); % local variables
    space_srv = state_i(:,(end-sum(K)-Vagg+1):(end-Vagg)); % server state
    space_buf = state_i(:,1:(end-sum(K)-Vagg)); % buffer state
end

nir = zeros(size(state_i,1),R); % class-r jobs in service
for r=1:R
    for k=1:K(r)
        nir(:,r) = nir(:,r) + space_srv(:,Ks(r)+k);
    end
end
switch sn.sched(ist)
    case SchedStrategy.EXT
        for r=1:R
            nir(:,r) = Inf;
        end
    case SchedStrategy.FCFS
        for r=1:R
            nir(:,r) = nir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case SchedStrategy.HOL
        for r=1:R
            nir(:,r) = nir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case SchedStrategy.LCFS
        for r=1:R
            nir(:,r) = nir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case {SchedStrategy.SIRO, SchedStrategy.POLLING}
        for r=1:R
            nir(:,r) = nir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    case SchedStrategy.SEPT
        for r=1:R
            nir(:,r) = nir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    case SchedStrategy.LEPT
        for r=1:R
            nir(:,r) = nir(:,r) + space_buf(:,r); % class-r jobs in station
        end
        %otherwise % possibly other stateful nodes
        % no-op
end

for r=1:R
    if isnan(sn.rates(ist,r)) && sn.nodetype(ind) ~= NodeType.Place % if disabled
        nir(:,r) = 0;
    end
end

ni = sum(nir,2); % total jobs in station
end
