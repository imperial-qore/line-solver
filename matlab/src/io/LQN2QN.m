function model = LQN2QN(lqn)
% LQN2QN Convert a LayeredNetwork (LQN) to a Network (QN) using REPLY signals
%
% model = LQN2QN(lqn) converts a LayeredNetwork model to an equivalent
% queueing network that uses REPLY signals to model synchronous call blocking.
%
% REPLY Signal Semantics:
% - Each synchCall creates a Request class and a Reply signal class
% - The caller queue blocks after completing service until Reply arrives
% - The callee processes the request and class-switches to Reply on completion
% - Reply signal unblocks the caller and continues downstream
%
% Network Topology:
%   Think -> CallerQueue -> CalleeQueue -> CallerQueue (reply) -> Think
%                 | blocks      | class switch to Reply    | unblocks
%
% This approach:
% - Correctly models BLOCKING semantics (server waits for reply)
% - Provides per-task queue metrics (queue length, utilization)
% - Supports multi-tier call chains (A -> B -> C)
% - Uses DES solver with REPLY signal support (via JAR)
%
% Example:
%   lqn = LayeredNetwork('MyLQN');
%   % ... define LQN model ...
%   model = LQN2QN(lqn);
%   SolverJMT(model).getAvgTable()  % or SolverDES via JAR
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% Get LQN structure
lsn = lqn.getStruct();

%% Create QN model
model = Network([lqn.getName(), '-QN']);

%% Identify reference tasks
refTaskIndices = find(lsn.isref);
nRefTasks = length(refTaskIndices);

if nRefTasks == 0
    line_error(mfilename, 'LQN must have at least one reference task.');
end

%% Detect phase-2 activities
hasPhase2 = isfield(lsn, 'actphase') && ~isempty(lsn.actphase) && any(lsn.actphase > 1);

%% Build task service demands (split by phase)
taskServiceByPhase = containers.Map('KeyType', 'double', 'ValueType', 'any');
for t = 1:lsn.ntasks
    tidx = lsn.tshift + t;
    [ph1Demand, ph2Demand] = getTaskServiceDemandByPhase(lsn, tidx, hasPhase2);
    taskServiceByPhase(tidx) = [ph1Demand, ph2Demand];
end

%% Build call graph: for each task, find all tasks it calls synchronously
synchCallsFrom = containers.Map('KeyType', 'double', 'ValueType', 'any');
for t = 1:lsn.ntasks
    tidx = lsn.tshift + t;
    calls = [];  % Array of structs: targetTidx, targetEidx, callMean

    entries = lsn.entriesof{tidx};
    for eidx = entries
        if eidx > length(lsn.actsof) || isempty(lsn.actsof{eidx})
            continue;
        end
        acts = lsn.actsof{eidx};
        for aidx = acts
            if lsn.type(aidx) ~= LayeredNetworkElement.ACTIVITY
                continue;
            end
            if aidx > length(lsn.callsof) || isempty(lsn.callsof{aidx})
                continue;
            end
            callList = lsn.callsof{aidx};
            for cidx = callList
                if full(lsn.calltype(cidx)) == CallType.SYNC
                    targetEidx = lsn.callpair(cidx, 2);
                    targetTidx = lsn.parent(targetEidx);
                    % Get call mean (number of calls)
                    callMean = 1.0;
                    if isfield(lsn, 'callproc') && ~isempty(lsn.callproc) && ...
                            cidx <= length(lsn.callproc) && isa(lsn.callproc{cidx}, 'Distribution')
                        callMean = lsn.callproc{cidx}.getMean();
                    end
                    callInfo = struct('targetTidx', targetTidx, 'targetEidx', targetEidx, 'callMean', callMean);
                    calls = [calls, callInfo]; %#ok<AGROW>
                end
            end
        end
    end
    synchCallsFrom(tidx) = calls;
end

%% Find all tasks in call chains starting from reference tasks
tasksInChain = [];
for rt = 1:nRefTasks
    refTidx = refTaskIndices(rt);
    tasksInChain = collectTasksInChain(refTidx, synchCallsFrom, tasksInChain);
end
tasksInChain = unique(tasksInChain);

%% Create nodes for all tasks
taskQueues = containers.Map('KeyType', 'double', 'ValueType', 'any');  % tidx -> Queue/Delay
thinkNodes = containers.Map('KeyType', 'double', 'ValueType', 'any');  % refTidx -> Think Delay

% Create think nodes for reference tasks
for rt = 1:nRefTasks
    refTidx = refTaskIndices(rt);
    taskName = lsn.names{refTidx};
    thinkNode = Delay(model, [taskName, '_Think']);
    thinkNodes(refTidx) = thinkNode;
end

% Create queues for all tasks in call chains
for i = 1:length(tasksInChain)
    tidx = tasksInChain(i);
    taskName = lsn.names{tidx};
    procIdx = lsn.parent(tidx);
    nServers = lsn.mult(procIdx);
    sched = lsn.sched(procIdx);

    if isinf(nServers) || sched == SchedStrategy.INF
        queue = Delay(model, taskName);
    else
        queue = Queue(model, taskName, sched);
        queue.setNumberOfServers(nServers);
    end
    taskQueues(tidx) = queue;
end

%% Create job classes for each reference task
% Each ref task gets: Request class + Reply signal class
requestClasses = containers.Map('KeyType', 'double', 'ValueType', 'any');
replySignals = containers.Map('KeyType', 'double', 'ValueType', 'any');

for rt = 1:nRefTasks
    refTidx = refTaskIndices(rt);
    taskName = lsn.names{refTidx};
    thinkNode = thinkNodes(refTidx);
    population = lsn.mult(refTidx);

    % Create request class (closed)
    requestClass = ClosedClass(model, [taskName, '_Req'], population, thinkNode);
    requestClasses(refTidx) = requestClass;

    % Create reply signal class
    replySignal = Signal(model, [taskName, '_Reply'], SignalType.REPLY);
    replySignal.forJobClass(requestClass);
    replySignals(refTidx) = replySignal;

    % IMPORTANT: Override the default RAND routing that Signal.m sets
    % We need DISABLED routing so link() can set the proper probabilistic routes
    for n = 1:length(model.nodes)
        if ~isa(model.nodes{n}, 'Sink')  % Don't change sink routing
            model.nodes{n}.setRouting(replySignal, RoutingStrategy.DISABLED);
        end
    end
end

%% Set service times for all nodes and classes
for rt = 1:nRefTasks
    refTidx = refTaskIndices(rt);
    requestClass = requestClasses(refTidx);
    replySignal = replySignals(refTidx);
    thinkNode = thinkNodes(refTidx);

    % Think time at think node
    thinkDist = lsn.think{refTidx};
    if isa(thinkDist, 'Immediate') || thinkDist.getMean() < GlobalConstants.FineTol
        thinkNode.setService(requestClass, Exp(1e8));
    else
        thinkNode.setService(requestClass, thinkDist);
    end
    % Reply passes through think instantly (shouldn't visit, but set for safety)
    thinkNode.setService(replySignal, Exp(1e9));

    % Set service times at all task queues
    for i = 1:length(tasksInChain)
        tidx = tasksInChain(i);
        queue = taskQueues(tidx);
        phaseDemands = taskServiceByPhase(tidx);
        serviceMean = phaseDemands(1) + phaseDemands(2);  % Total service

        if serviceMean > GlobalConstants.FineTol
            queue.setService(requestClass, Exp(1/serviceMean));
        else
            queue.setService(requestClass, Exp(1e8));
        end
        % Reply signal passes through instantly (unblocks at caller)
        queue.setService(replySignal, Exp(1e9));
    end
end

%% Build routing matrix with class switching for REPLY signals
P = model.initRoutingMatrix();

for rt = 1:nRefTasks
    refTidx = refTaskIndices(rt);
    requestClass = requestClasses(refTidx);
    replySignal = replySignals(refTidx);
    thinkNode = thinkNodes(refTidx);

    % Build the complete call chain from reference task to leaf
    % callChain = [refTidx, callee1, callee2, ..., leafTidx]
    callChain = buildCallChain(refTidx, synchCallsFrom);

    if length(callChain) == 1
        % No synch calls - just loop Think -> RefQueue -> Think
        if isKey(taskQueues, refTidx)
            refQueue = taskQueues(refTidx);
            P{requestClass, requestClass}(thinkNode, refQueue) = 1.0;
            P{requestClass, requestClass}(refQueue, thinkNode) = 1.0;
        else
            P{requestClass, requestClass}(thinkNode, thinkNode) = 1.0;
        end
        continue;
    end

    % Build Request path: Think -> task1 -> task2 -> ... -> leafTask
    % First: Think -> first task in chain
    firstTidx = callChain(1);
    if isKey(taskQueues, firstTidx)
        firstQueue = taskQueues(firstTidx);
        P{requestClass, requestClass}(thinkNode, firstQueue) = 1.0;
    else
        % Reference task has no queue (e.g., think-only task)
        % Skip directly to the first callee
        if length(callChain) > 1
            firstQueue = taskQueues(callChain(2));
            P{requestClass, requestClass}(thinkNode, firstQueue) = 1.0;
            callChain = callChain(2:end);  % Remove ref task from chain
        end
    end

    % Request path through the call chain
    for c = 1:(length(callChain) - 1)
        fromTidx = callChain(c);
        toTidx = callChain(c + 1);
        if isKey(taskQueues, fromTidx) && isKey(taskQueues, toTidx)
            fromQueue = taskQueues(fromTidx);
            toQueue = taskQueues(toTidx);
            P{requestClass, requestClass}(fromQueue, toQueue) = 1.0;
        end
    end

    % Class switch at the LEAF node (last in chain)
    % Reply path: leafTask -> ... -> task1 -> Think
    leafTidx = callChain(end);
    leafQueue = taskQueues(leafTidx);

    if length(callChain) == 1
        % Only one task in chain (the ref task itself)
        P{requestClass, replySignal}(leafQueue, thinkNode) = 1.0;
    else
        % Class switch at leaf, reply flows back through chain
        prevTidx = callChain(end - 1);
        prevQueue = taskQueues(prevTidx);
        P{requestClass, replySignal}(leafQueue, prevQueue) = 1.0;

        % Reply flows back through intermediate nodes
        for c = (length(callChain) - 1):-1:2
            fromTidx = callChain(c);
            toTidx = callChain(c - 1);
            fromQueue = taskQueues(fromTidx);
            toQueue = taskQueues(toTidx);
            P{replySignal, replySignal}(fromQueue, toQueue) = 1.0;
        end

        % Final hop: first task in chain -> Think (class switch back to Request)
        % Reply arriving at Think triggers class-switch back to Request to complete the cycle
        firstQueue = taskQueues(callChain(1));
        P{replySignal, requestClass}(firstQueue, thinkNode) = 1.0;
    end
end

model.link(P);

end

%% Helper: Build the call chain from a starting task to the leaf
function callChain = buildCallChain(startTidx, synchCallsFrom)
% Returns array of task indices from start to leaf (following first synch call at each level)

callChain = startTidx;
currentTidx = startTidx;

while true
    calls = synchCallsFrom(currentTidx);
    if isempty(calls)
        break;  % Reached leaf
    end
    % Follow the first synch call
    nextTidx = calls(1).targetTidx;
    if any(callChain == nextTidx)
        break;  % Avoid cycles
    end
    callChain = [callChain, nextTidx]; %#ok<AGROW>
    currentTidx = nextTidx;
end

end

%% Helper: Collects all tasks reachable via synchronous calls
function tasksInChain = collectTasksInChain(startTidx, synchCallsFrom, tasksInChain)

if any(tasksInChain == startTidx)
    return;
end
tasksInChain = [tasksInChain, startTidx];

calls = synchCallsFrom(startTidx);
if ~isempty(calls)
    for i = 1:length(calls)
        call = calls(i);
        tasksInChain = collectTasksInChain(call.targetTidx, synchCallsFrom, tasksInChain);
    end
end

end

%% Helper: Get service demand for a task split by phase
function [ph1Demand, ph2Demand] = getTaskServiceDemandByPhase(lsn, tidx, hasPhase2)
ph1Demand = 0;
ph2Demand = 0;

entries = lsn.entriesof{tidx};
for eidx = entries
    if eidx > length(lsn.actsof) || isempty(lsn.actsof{eidx})
        continue;
    end
    acts = lsn.actsof{eidx};
    for aidx = acts
        if lsn.type(aidx) == LayeredNetworkElement.ACTIVITY
            hostDem = lsn.hostdem{aidx};
            demand = 0;
            if isa(hostDem, 'Distribution')
                demand = hostDem.getMean();
            end

            if hasPhase2
                a = aidx - lsn.ashift;
                if a >= 1 && a <= length(lsn.actphase) && lsn.actphase(a) > 1
                    ph2Demand = ph2Demand + demand;
                else
                    ph1Demand = ph1Demand + demand;
                end
            else
                ph1Demand = ph1Demand + demand;
            end
        end
    end
end

end
