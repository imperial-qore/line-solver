function [supported, reason, blocking] = solver_nc_mem_supports(sn)
% [SUPPORTED, REASON, BLOCKING] = SOLVER_NC_MEM_SUPPORTS(SN)
%
% Checks whether the Maximum Entropy Method (Kouvatsos 1994) supports the
% model described by the network structure SN. Returns SUPPORTED=true and
% an empty REASON when the model is a plain open queueing network
% (Section 3.2: Source, Queue, Delay and Sink nodes; GE/GE/1, GE/GE/c and
% GE/GE/inf building blocks), a plain closed queueing network
% (Section 3.3: Queue and Delay nodes; G/G/1 and G/G/inf building blocks
% only, so finite multiserver stations are rejected), or a mixed
% open/closed network (composition of the two algorithms; single-server
% and IS stations only), in all cases without class switching and with
% non-priority scheduling disciplines only; otherwise SUPPORTED=false and
% REASON explains the first unsupported feature found.
%
% An open model with finite station buffers is also supported, provided it
% is single class: it is then solved by the censored GE/GE/c/0;N building
% block of Section 4.1, under loss (drop rule DROP) or transfer blocking
% (drop rule BAS, holding-node expansion of Tahilramani, Manjunath and
% Bose 1999). BLOCKING returns true for such a model, so the caller can
% route it to ME_OQN_BLK instead of ME_OQN. The GE distribution is only
% defined for scv >= 1, so a finite-buffer model with a hypo-exponential
% service or arrival process is rejected rather than approximated.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

supported = false;
reason = '';
blocking = false;

isopen = sn_is_open_model(sn);
isclosed = sn_is_closed_model(sn);
ismixed = ~isopen && ~isclosed;

R = sn.nclasses;

% Node types: open models are Source/Queue/Delay/Sink networks, closed
% models are Queue/Delay networks
for ind = 1:sn.nnodes
    switch sn.nodetype(ind)
        case {NodeType.Queue, NodeType.Delay}
            % supported
        case {NodeType.Source, NodeType.Sink}
            if isclosed
                reason = 'MEM supports only Queue and Delay nodes in closed models.';
                return
            end
        otherwise
            reason = 'MEM supports only Source, Queue, Delay and Sink nodes.';
            return
    end
end

% Class switching is not part of the Kouvatsos (1994) network model
if any(any(sn.csmask & ~eye(R)))
    reason = 'MEM does not support class switching.';
    return
end

% Absorbing self-loops (p_ii=1) make the routing reducible and the
% geometric feedback transform 1/(1-p_ii) degenerate
for ist = 1:sn.nstations
    ind = sn.stationToNode(ist);
    for r = 1:R
        if sn.rtnodes((ind-1)*R + r, (ind-1)*R + r) >= 1 - 1e-9
            reason = 'MEM does not support absorbing self-loop routing (reducible network).';
            return
        end
    end
end

% Only non-priority disciplines are supported; the PR/HOL constraint
% formulae are not given in Kouvatsos (1994)
for ist = 1:sn.nstations
    switch sn.sched(ist)
        case {SchedStrategy.EXT, SchedStrategy.INF, SchedStrategy.FCFS, ...
                SchedStrategy.PS, SchedStrategy.SIRO, SchedStrategy.LCFS, ...
                SchedStrategy.LCFSPR}
            % supported
        otherwise
            reason = sprintf('MEM does not support the %s scheduling strategy.', SchedStrategy.toText(sn.sched(ist)));
            return
    end
end

if isopen || ismixed
    % A Source node must be present for the external arrival extraction
    hasSource = false;
    for ind = 1:sn.nnodes
        if sn.nodetype(ind) == NodeType.Source
            hasSource = true;
            break;
        end
    end
    if ~hasSource
        reason = 'MEM requires a Source node when open classes are present.';
        return
    end
end
if isclosed || ismixed
    % Closed classes build on G/G/1 and G/G/inf queues only (Section 3.3)
    for ist = 1:sn.nstations
        if isfinite(sn.nservers(ist)) && sn.nservers(ist) > 1
            reason = 'MEM does not support multiserver stations in closed or mixed models.';
            return
        end
    end
end

% Finite buffers admissible only in a single-class open model (the GE/GE/c/0;N
% building block is single class); see _kb/06-solver-catalog.md (mem.blocking note)
capped = false(sn.nstations, 1);
for ist = 1:sn.nstations
    if sn.sched(ist) == SchedStrategy.EXT
        continue
    end
    if isfinite(sn_get_buffer_size(sn, ist))
        capped(ist) = true;
    end
end
if any(capped)
    if ~isopen
        reason = 'MEM supports finite station buffers only in open models.';
        return
    end
    if R > 1
        reason = 'MEM supports finite station buffers only in single-class models: the censored GE/GE/c/0;N building block is single class.';
        return
    end
    for ist = 1:sn.nstations
        if ~capped(ist)
            continue
        end
        if ~isfinite(sn.nservers(ist)) || sn.nservers(ist) < 1
            reason = sprintf('MEM cannot apply a finite buffer to the infinite-server station %d.', ist);
            return
        end
        if sn.sched(ist) ~= SchedStrategy.FCFS
            reason = sprintf('MEM supports finite station buffers only under FCFS scheduling; station %d uses %s.', ist, SchedStrategy.toText(sn.sched(ist)));
            return
        end
        dr = DropStrategy.DROP;
        if ~isempty(sn.droprule) && size(sn.droprule, 1) >= ist
            dr = sn.droprule(ist, 1);
        end
        if dr ~= DropStrategy.DROP && dr ~= DropStrategy.BAS
            reason = sprintf('MEM supports the DROP and BAS drop rules at a finite buffer; station %d uses %s.', ist, DropStrategy.toText(dr));
            return
        end
        if isfinite(sn.scv(ist, 1)) && sn.scv(ist, 1) < 1 - 1e-12
            reason = sprintf('MEM with finite buffers needs a service scv of at least 1 at station %d: the GE distribution is not defined below 1.', ist);
            return
        end
    end
    for ist = 1:sn.nstations
        if sn.sched(ist) == SchedStrategy.EXT && isfinite(sn.scv(ist, 1)) && sn.scv(ist, 1) < 1 - 1e-12
            reason = 'MEM with finite buffers needs an external interarrival scv of at least 1: the GE distribution is not defined below 1.';
            return
        end
    end
    blocking = true;
end

supported = true;
end
