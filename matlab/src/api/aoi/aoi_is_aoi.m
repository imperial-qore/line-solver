function [isAoI, aoiInfo] = aoi_is_aoi(sn)
%AOI_IS_AOI Check if network is a valid AoI topology for aoi-fluid analysis
%
% [isAoI, aoiInfo] = AOI_IS_AOI(sn)
%
% Validates that the network structure is suitable for Age of Information
% analysis using the aoi-fluid MFQ solvers.
%
% Requirements:
% - Single open class
% - Single-queue open system: Source -> Queue -> Sink
% - Queue capacity = 1 (bufferless) or 2 (single-buffer)
% - Single server (nservers = 1)
% - Scheduling: FCFS, LCFS, or LCFSPR
% - For capacity=2: arrivals must be exponential (Poisson)
%
% Parameters:
%   sn (struct): Network structure from Network.getStruct()
%
% Returns:
%   isAoI (logical): True if topology is valid for AoI analysis
%   aoiInfo (struct): Contains topology information:
%       .sourceIdx, .queueIdx, .sinkIdx: Node indices
%       .sourceStation, .queueStation: Station indices
%       .capacity: Queue capacity (1 or 2)
%       .schedStrategy: Scheduling strategy
%       .systemType: 'bufferless' or 'singlebuffer'
%       .errorMsg: Error message if not valid

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

aoiInfo = struct();
aoiInfo.errorMsg = '';
aoiInfo.sourceIdx = [];
aoiInfo.queueIdx = [];
aoiInfo.sinkIdx = [];
aoiInfo.sourceStation = [];
aoiInfo.queueStation = [];
aoiInfo.capacity = [];
aoiInfo.schedStrategy = [];
aoiInfo.systemType = '';

% Check 1: Must have at least one open class
if ~any(isinf(sn.njobs))
    aoiInfo.errorMsg = 'Not an open model - all classes are closed';
    isAoI = false;
    return;
end

% Check 2: Must have exactly one open class
openClasses = find(isinf(sn.njobs));
if length(openClasses) > 1
    aoiInfo.errorMsg = sprintf('Multiple open classes found (%d) - AoI analysis requires single class', length(openClasses));
    isAoI = false;
    return;
end

% Find Source, Queue, and Sink nodes
sourceNodeIdx = find(sn.nodetype == NodeType.Source);
queueNodeIdx = find(sn.nodetype == NodeType.Queue);
sinkNodeIdx = find(sn.nodetype == NodeType.Sink);

% Check: exactly one source
if isempty(sourceNodeIdx)
    aoiInfo.errorMsg = 'No source node found';
    isAoI = false;
    return;
elseif length(sourceNodeIdx) > 1
    aoiInfo.errorMsg = sprintf('Multiple source nodes found (%d)', length(sourceNodeIdx));
    isAoI = false;
    return;
end

% Check: exactly one sink
if isempty(sinkNodeIdx)
    aoiInfo.errorMsg = 'No sink node found';
    isAoI = false;
    return;
elseif length(sinkNodeIdx) > 1
    aoiInfo.errorMsg = sprintf('Multiple sink nodes found (%d)', length(sinkNodeIdx));
    isAoI = false;
    return;
end

% Check: exactly one queue
if isempty(queueNodeIdx)
    aoiInfo.errorMsg = 'No queue node found';
    isAoI = false;
    return;
elseif length(queueNodeIdx) > 1
    aoiInfo.errorMsg = sprintf('Multiple queue nodes found (%d) - AoI analysis supports single queue only', length(queueNodeIdx));
    isAoI = false;
    return;
end

% Get station indices
sourceStation = sn.nodeToStation(sourceNodeIdx);
queueStation = sn.nodeToStation(queueNodeIdx);

aoiInfo.sourceIdx = sourceNodeIdx;
aoiInfo.queueIdx = queueNodeIdx;
aoiInfo.sinkIdx = sinkNodeIdx;
aoiInfo.sourceStation = sourceStation;
aoiInfo.queueStation = queueStation;

% Check: Single server
c = sn.nservers(queueStation);
if c ~= 1
    aoiInfo.errorMsg = sprintf('Queue has %d servers - AoI analysis requires single server (c=1)', c);
    isAoI = false;
    return;
end

% Check: Queue capacity (1 for bufferless, 2 for single-buffer)
% sn.cap is a vector indexed by station
cap = sn.cap(queueStation);

% Handle infinite or invalid capacity
if isinf(cap) || cap > 2 || cap < 1
    aoiInfo.errorMsg = sprintf('Queue capacity is %g - AoI analysis requires capacity 1 (bufferless) or 2 (single-buffer)', cap);
    isAoI = false;
    return;
end

aoiInfo.capacity = cap;
if cap == 1
    aoiInfo.systemType = 'bufferless';
else
    aoiInfo.systemType = 'singlebuffer';
end

% Check: Scheduling strategy (FCFS, LCFS, or LCFSPR)
schedStrategy = sn.sched(queueStation);
aoiInfo.schedStrategy = schedStrategy;

if schedStrategy ~= SchedStrategy.FCFS && ...
   schedStrategy ~= SchedStrategy.LCFS && ...
   schedStrategy ~= SchedStrategy.LCFSPR
    aoiInfo.errorMsg = sprintf('Unsupported scheduling strategy - AoI analysis supports FCFS, LCFS, or LCFSPR only');
    isAoI = false;
    return;
end

% For single-buffer (capacity=2): arrivals must be exponential (Poisson)
if cap == 2
    % Check the arrival process at the source
    arrivalProc = sn.proc{sourceStation}{openClasses(1)};
    if ~isempty(arrivalProc) && ~isnan(arrivalProc{1}(1))
        % Check if it's exponential (single phase)
        D0 = arrivalProc{1};  % MAP representation
        if size(D0, 1) > 1
            aoiInfo.errorMsg = 'Single-buffer (capacity=2) requires exponential arrivals (Poisson process)';
            isAoI = false;
            return;
        end
    end
end

% Check: No self-loops
K = sn.nclasses;
rtIdx = (queueStation - 1) * K + openClasses(1);
if sn.rt(rtIdx, rtIdx) > 0
    aoiInfo.errorMsg = 'Self-loop detected at queue - violates AoI model assumptions';
    isAoI = false;
    return;
end

% All checks passed
isAoI = true;
end
