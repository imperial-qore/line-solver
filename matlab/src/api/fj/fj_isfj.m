%{ @file fj_isfj.m
 %  @brief Checks if network is a valid Fork-Join topology for FJ_codes analysis
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if network is a valid Fork-Join topology for FJ_codes analysis
 %
 % @details
 % Validates that the network structure matches the requirements for FJ_codes:
 % - Single Fork-Join pair
 % - K parallel queues between Fork and Join
 % - Homogeneous service distributions across parallel queues
 % - Supported distributions (Exp, HyperExp(2), Erlang(2), MAP(2))
 % - Open classes only
 %
 % @par Syntax:
 % @code
 % [isFJ, fjInfo] = fj_isfj(sn)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>isFJ<td>True if network is valid FJ topology for FJ_codes
 % <tr><td>fjInfo<td>Struct with fields: forkIdx, joinIdx, queueIdx, K, errorMsg
 % </table>
 %
 % @par Reference:
 % Z. Qiu, J.F. Pérez, and P. Harrison, "Beyond the Mean in Fork-Join Queues:
 % Efficient Approximation for Response-Time Tails", IFIP Performance 2015.
 % Copyright 2015 Imperial College London
%}
function [isFJ, fjInfo] = fj_isfj(sn)

% Initialize output
isFJ = false;
fjInfo = struct();
fjInfo.forkIdx = [];
fjInfo.joinIdx = [];
fjInfo.queueIdx = [];
fjInfo.K = 0;
fjInfo.errorMsg = '';

% Check if model has open classes only
if ~sn_is_open_model(sn)
    fjInfo.errorMsg = 'FJ_codes only supports open queueing models.';
    return;
end

% Check if network has fork-join
if ~sn_has_fork_join(sn)
    fjInfo.errorMsg = 'Network does not contain Fork-Join structure.';
    return;
end

% Find Fork and Join nodes
forkIndices = find(sn.nodetype == NodeType.Fork);
joinIndices = find(sn.nodetype == NodeType.Join);

if isempty(forkIndices) || isempty(joinIndices)
    fjInfo.errorMsg = 'Network must contain both Fork and Join nodes.';
    return;
end

% FJ_codes supports single Fork-Join pair
if length(forkIndices) > 1
    fjInfo.errorMsg = 'FJ_codes only supports a single Fork-Join pair. Found multiple Fork nodes.';
    return;
end

if length(joinIndices) > 1
    fjInfo.errorMsg = 'FJ_codes only supports a single Fork-Join pair. Found multiple Join nodes.';
    return;
end

forkIdx = forkIndices(1);
joinIdx = joinIndices(1);

% Check if Fork and Join are paired using sn.fj matrix
if sn.fj(forkIdx, joinIdx) == 0
    fjInfo.errorMsg = sprintf('Fork node %d and Join node %d are not paired.', forkIdx, joinIdx);
    return;
end

fjInfo.forkIdx = forkIdx;
fjInfo.joinIdx = joinIdx;

% Find queues between Fork and Join
% These are nodes that receive routing from Fork and route to Join
queueIdx = [];
for i = 1:sn.nnodes
    if sn.nodetype(i) == NodeType.Queue
        % Check if this queue is in a path from Fork to Join
        % Use rtnodes which is indexed by node indices (not station indices)
        hasForkInput = false;
        hasJoinOutput = false;

        % Check if Fork routes to this queue (using rtnodes for node-based routing)
        if sn.rtnodes(forkIdx, i) > 0
            hasForkInput = true;
        end
        % Check if this queue routes to Join
        if sn.rtnodes(i, joinIdx) > 0
            hasJoinOutput = true;
        end

        if hasForkInput && hasJoinOutput
            queueIdx = [queueIdx, i];
        end
    end
end

if isempty(queueIdx)
    fjInfo.errorMsg = 'No Queue nodes found between Fork and Join.';
    return;
end

K = length(queueIdx);
fjInfo.queueIdx = queueIdx;
fjInfo.K = K;

% Validate homogeneous service distributions across parallel queues
% For each class, all K queues must have the same service distribution
for r = 1:sn.nclasses
    % Get PH representation of first queue's service distribution
    firstQueueIdx = queueIdx(1);
    firstPH = sn.proc{sn.nodeToStation(firstQueueIdx)}{r};

    if isempty(firstPH) || isnan(firstPH{1}(1))
        fjInfo.errorMsg = sprintf('Queue %d has no valid service distribution for class %d.', ...
            firstQueueIdx, r);
        return;
    end

    % Check all other queues have the same distribution
    for k = 2:K
        queueIdx_k = queueIdx(k);
        ph_k = sn.proc{sn.nodeToStation(queueIdx_k)}{r};

        if isempty(ph_k) || isnan(ph_k{1}(1))
            fjInfo.errorMsg = sprintf('Queue %d has no valid service distribution for class %d.', ...
                queueIdx_k, r);
            return;
        end

        % Compare PH representations (must be identical)
        % Compare number of phases
        if length(ph_k{1}) ~= length(firstPH{1})
            fjInfo.errorMsg = sprintf('Queues have heterogeneous service distributions for class %d. FJ_codes requires homogeneous servers.', r);
            return;
        end

        % Compare initial probability vector and rate matrix
        if ~isequal(size(ph_k{1}), size(firstPH{1})) || ...
           ~isequal(size(ph_k{2}), size(firstPH{2})) || ...
           max(abs(ph_k{1} - firstPH{1})) > GlobalConstants.FineTol || ...
           max(max(abs(ph_k{2} - firstPH{2}))) > GlobalConstants.FineTol
            fjInfo.errorMsg = sprintf('Queues have heterogeneous service distributions for class %d. FJ_codes requires homogeneous servers.', r);
            return;
        end
    end
end

% Validate supported scheduling strategies (FCFS or PS)
for k = 1:K
    queueSt = sn.nodeToStation(queueIdx(k));
    if sn.sched(queueSt) ~= SchedStrategy.FCFS && sn.sched(queueSt) ~= SchedStrategy.PS
        fjInfo.errorMsg = sprintf('Queue %d has unsupported scheduling strategy. FJ_codes supports FCFS or PS only.', queueIdx(k));
        return;
    end
end

% All validations passed
isFJ = true;

end
