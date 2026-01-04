function [isSingleQueue, fluidInfo] = fluid_is_single_queue(sn)
%FLUID_IS_SINGLE_QUEUE Check if network is a valid single-queue topology for MFQ
%
% [isSingleQueue, fluidInfo] = FLUID_IS_SINGLE_QUEUE(sn)
%
% Validates that the network structure is suitable for the MFQ solver.
% MFQ requires: Source -> single Queue -> Sink topology with single server.
%
% Parameters:
%   sn (struct): Network structure from Network.getStruct()
%
% Returns:
%   isSingleQueue (logical): True if topology is valid for MFQ
%   fluidInfo (struct): Contains topology information and error message if invalid

    fluidInfo = struct();
    fluidInfo.errorMsg = '';
    fluidInfo.sourceIdx = [];
    fluidInfo.queueIdx = [];
    fluidInfo.sinkIdx = [];
    fluidInfo.sourceStation = [];
    fluidInfo.queueStation = [];

    % Check 1: Must have at least one open class
    if ~any(isinf(sn.njobs))
        fluidInfo.errorMsg = 'Not an open model - all classes are closed';
        isSingleQueue = false;
        return;
    end

    % Work at the node level (like fj_isfj.m does)
    % Find Source, Queue, and Sink nodes
    sourceNodeIdx = find(sn.nodetype == NodeType.Source);
    queueNodeIdx = find(sn.nodetype == NodeType.Queue);
    sinkNodeIdx = find(sn.nodetype == NodeType.Sink);

    % Check: exactly one source
    if length(sourceNodeIdx) == 0
        fluidInfo.errorMsg = 'No source node found';
        isSingleQueue = false;
        return;
    elseif length(sourceNodeIdx) > 1
        fluidInfo.errorMsg = sprintf('Multiple source nodes found (%d)', length(sourceNodeIdx));
        isSingleQueue = false;
        return;
    end

    % Check: exactly one sink
    if length(sinkNodeIdx) == 0
        fluidInfo.errorMsg = 'No sink node found';
        isSingleQueue = false;
        return;
    elseif length(sinkNodeIdx) > 1
        fluidInfo.errorMsg = sprintf('Multiple sink nodes found (%d)', length(sinkNodeIdx));
        isSingleQueue = false;
        return;
    end

    % Check: exactly one queue
    if length(queueNodeIdx) == 0
        fluidInfo.errorMsg = 'No queue node found';
        isSingleQueue = false;
        return;
    elseif length(queueNodeIdx) > 1
        fluidInfo.errorMsg = sprintf('Multiple queue nodes found (%d) - MFQ supports single queue only', length(queueNodeIdx));
        isSingleQueue = false;
        return;
    end

    % Get station indices (for stations that exist)
    sourceStation = sn.nodeToStation(sourceNodeIdx);
    queueStation = sn.nodeToStation(queueNodeIdx);

    % Check: Single-server or infinite-server constraint
    c = sn.nservers(queueStation);
    if c > 1 && c < Inf
        fluidInfo.errorMsg = sprintf('Multi-server queue (c=%d) not supported - MFQ requires c=1 or c=Inf', c);
        isSingleQueue = false;
        return;
    end

    % Check: No self-loops (independence assumption)
    % Check routing matrix for self-loops in the queue station
    K = sn.nclasses;
    for k = 1:K
        % Check if queue routes back to itself
        rtIdx = (queueStation - 1) * K + k;
        if sn.rt(rtIdx, rtIdx) > 0
            fluidInfo.errorMsg = 'Self-loop detected at queue - violates independence assumption';
            isSingleQueue = false;
            return;
        end
    end

    % All checks passed
    fluidInfo.sourceIdx = sourceNodeIdx;
    fluidInfo.queueIdx = queueNodeIdx;
    fluidInfo.sinkIdx = sinkNodeIdx;
    fluidInfo.sourceStation = sourceStation;
    fluidInfo.queueStation = queueStation;
    isSingleQueue = true;
end
