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
if ~sn.isstation(ind) && sn.isstateful(ind) % if stateful node
    ni = sum(state_i(1:(end-sum(sn.nvars(ind,:)))));
    nir = state_i(1:(end-sum(sn.nvars(ind,:))));
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
    case SchedStrategy.SIRO
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
