%{ @file fj_extract_params.m
 %  @brief Extract FJ_codes parameters from LINE network structure
 %
 %  @author LINE Development Team
%}

%{
 % @brief Extract FJ_codes parameters from LINE network structure
 %
 % @details
 % Extracts arrival process, service process, and K value from a validated
 % Fork-Join network structure for use with FJ_codes analysis.
 %
 % @par Syntax:
 % @code
 % [arrival, service, K, fjInfo] = fj_extract_params(sn, fjInfo)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure (after sn_nonmarkov_toph conversion)
 % <tr><td>fjInfo<td>FJ topology info from fj_isfj (forkIdx, joinIdx, queueIdx, K)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>arrival<td>Cell array of arrival structs (one per class) with fields: lambda, lambda0, lambda1, ma, Ia
 % <tr><td>service<td>Cell array of service structs (one per class) with fields: mu, ST, St, tau_st, SerChoice
 % <tr><td>K<td>Number of parallel queues
 % <tr><td>fjInfo<td>Updated fjInfo with distribution info
 % </table>
 %
 % @par Reference:
 % Z. Qiu, J.F. Pérez, and P. Harrison, "Beyond the Mean in Fork-Join Queues:
 % Efficient Approximation for Response-Time Tails", IFIP Performance 2015.
 % Copyright 2015 Imperial College London
%}
function [arrival, service, K, fjInfo] = fj_extract_params(sn, fjInfo)

% Get Fork-Join structure info
forkIdx = fjInfo.forkIdx;
joinIdx = fjInfo.joinIdx;
queueIdx = fjInfo.queueIdx;
K = fjInfo.K;

% Find Source node to extract arrival process
sourceIdx = find(sn.nodetype == NodeType.Source);
if isempty(sourceIdx)
    line_error(mfilename, 'No Source node found in network.');
end
sourceIdx = sourceIdx(1); % Take first source
sourceStat = sn.nodeToStation(sourceIdx);

% Initialize cell arrays for multi-class support
nClasses = sn.nclasses;
arrival = cell(1, nClasses);
service = cell(1, nClasses);

% Extract parameters for each class
for r = 1:nClasses
    % ===== ARRIVAL PROCESS =====
    % Get arrival MAP from Source node
    arrivalMAP = sn.proc{sourceStat}{r};

    if isempty(arrivalMAP) || isnan(arrivalMAP{1}(1))
        line_error(mfilename, 'Source has no valid arrival process for class %d.', r);
    end

    % Convert to FJ_codes format
    arrival{r} = fj_dist2fj(arrivalMAP, 'arrival', sn, sourceStat, r);

    % ===== SERVICE PROCESS =====
    % Get service MAP from first queue (all queues are homogeneous)
    firstQueue = queueIdx(1);
    firstQueueStat = sn.nodeToStation(firstQueue);
    serviceMAP = sn.proc{firstQueueStat}{r};

    if isempty(serviceMAP) || isnan(serviceMAP{1}(1))
        line_error(mfilename, 'Queue %d has no valid service process for class %d.', ...
            firstQueue, r);
    end

    % Convert to FJ_codes format
    service{r} = fj_dist2fj(serviceMAP, 'service', sn, firstQueueStat, r);

    % ===== STABILITY CHECK =====
    % System must be stable: lambda < mu for each class
    % Each queue sees arrival rate lambda and has utilization rho = lambda/mu
    if arrival{r}.lambda >= service{r}.mu
        line_warning(mfilename, ...
            'Class %d may be unstable: arrival rate (%.4f) >= service rate (%.4f). FJ_codes requires stable systems (rho = lambda/mu < 1).', ...
            r, arrival{r}.lambda, service{r}.mu);
    end
end

% Store distribution info in fjInfo for reference
fjInfo.arrival = arrival;
fjInfo.service = service;
fjInfo.nClasses = nClasses;

end
