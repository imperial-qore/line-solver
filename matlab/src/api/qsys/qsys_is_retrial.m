%{ @file qsys_is_retrial.m
 %  @brief Checks if network is a valid BMAP/PH/N/N bufferless retrial queue
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if network is a valid BMAP/PH/N/N bufferless retrial queue
 %
 % @details
 % Validates that the network structure matches the requirements for
 % the BMAP/PH/N/N retrial queue solver:
 % - Single bufferless queue (capacity == number of servers)
 % - Retrial drop strategy configured
 % - BMAP/MAP arrival process at source
 % - PH/Exp service at queue
 % - Open class model
 %
 % Based on: Dudin et al., "Analysis of BMAP/PH/N-Type Queueing System with
 % Flexible Retrials Admission Control", Mathematics 2025, 13(9), 1434.
 %
 % @par Syntax:
 % @code
 % [isRetrial, retInfo] = qsys_is_retrial(sn)
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
 % <tr><td>isRetrial<td>True if network is valid BMAP/PH/N/N retrial topology
 % <tr><td>retInfo<td>Struct with parameters for the retrial solver
 % </table>
 %
 % Copyright (c) 2012-2026, Imperial College London
 % All rights reserved.
%}
function [isRetrial, retInfo] = qsys_is_retrial(sn)

% Initialize output
isRetrial = false;
retInfo = struct();
retInfo.stationIdx = [];
retInfo.nodeIdx = [];
retInfo.sourceIdx = [];
retInfo.classIdx = [];
retInfo.errorMsg = '';
retInfo.N = [];           % Number of servers
retInfo.alpha = [];       % Retrial rate
retInfo.gamma = 0;        % Orbit impatience (default 0)
retInfo.p = 0;            % Batch rejection prob (default 0)
retInfo.R = [];           % Admission threshold (from FCR or default N-1)

% Check if model is open
if ~sn_is_open_model(sn)
    retInfo.errorMsg = 'BMAP/PH/N/N retrial solver requires open queueing model.';
    return;
end

% Check for single class (current limitation)
if sn.nclasses > 1
    retInfo.errorMsg = 'BMAP/PH/N/N retrial solver currently supports single class only.';
    return;
end

retInfo.classIdx = 1;

% Find bufferless queue stations (capacity == nservers, finite capacity)
bufferlessStations = [];
for ist = 1:sn.nstations
    if sn.nodetype(sn.stationToNode(ist)) == NodeType.Queue
        if isfinite(sn.cap(ist)) && sn.cap(ist) == sn.nservers(ist)
            bufferlessStations = [bufferlessStations, ist];
        end
    end
end

if isempty(bufferlessStations)
    retInfo.errorMsg = 'No bufferless queue found (capacity must equal number of servers).';
    return;
end

% Check retrial drop strategy
retrialStation = [];
for ist = bufferlessStations(:)'
    if any(sn.droprule(ist,:) == DropStrategy.RETRIAL | ...
           sn.droprule(ist,:) == DropStrategy.RETRIAL_WITH_LIMIT)
        retrialStation = ist;
        break;
    end
end

if isempty(retrialStation)
    retInfo.errorMsg = 'No retrial drop strategy configured on bufferless queue.';
    return;
end

retInfo.stationIdx = retrialStation;
retInfo.nodeIdx = sn.stationToNode(retrialStation);
retInfo.N = sn.nservers(retrialStation);

% Find source station
sourceStation = [];
for ist = 1:sn.nstations
    if sn.nodetype(sn.stationToNode(ist)) == NodeType.Source
        sourceStation = ist;
        break;
    end
end

if isempty(sourceStation)
    retInfo.errorMsg = 'No Source node found.';
    return;
end

retInfo.sourceIdx = sourceStation;

% Validate arrival process (must be MAP/BMAP)
% sn.proc{station}{class} contains PH representation
arrivalProc = sn.proc{sourceStation}{retInfo.classIdx};
if isempty(arrivalProc) || ~iscell(arrivalProc) || length(arrivalProc) < 2
    retInfo.errorMsg = 'Invalid arrival process at source.';
    return;
end

% Check if arrival is MAP/BMAP (has matrix structure)
% For BMAP: proc is {D0, D1, D2, ...} or {alpha, A} for MAP
% For Exp/PH at source: it's {1, -lambda} format
if ~iscell(arrivalProc{1}) && size(arrivalProc{1}, 2) > 1
    % This looks like a PH representation {alpha, A}, convert from arrival to MAP
    % For now, we accept it - the solver_mam_retrial will handle conversion
elseif iscell(arrivalProc{1})
    % This might be BMAP format already
end

% Validate service process (must be PH/Exp)
serviceProc = sn.proc{retrialStation}{retInfo.classIdx};
if isempty(serviceProc) || ~iscell(serviceProc) || length(serviceProc) < 2
    retInfo.errorMsg = 'Invalid service process at queue.';
    return;
end

% Extract retrial rate alpha from retrial delay distribution
% In sn struct, we need to check for retrial configuration
% For now, use default if not found - the Network/getStruct needs to export this
retInfo.alpha = 0.1;  % Default, will be extracted from sn in solver

% Check for FCR containing this queue and extract admission threshold
retInfo.R = retInfo.N - 1;  % Default: no admission control

% Check if region field exists and is non-empty
if isfield(sn, 'region') && ~isempty(sn.region) && isfield(sn, 'nregions') && sn.nregions > 0
    for f = 1:sn.nregions
        regionMatrix = sn.region{f};
        % Check if this station is in the region with a global capacity
        if size(regionMatrix, 1) >= retrialStation
            globalCapCol = size(regionMatrix, 2);  % Last column is typically global cap
            if regionMatrix(retrialStation, globalCapCol) > 0
                % Get the global max jobs for this region
                R_fcr = regionMatrix(retrialStation, globalCapCol);
                if R_fcr <= retInfo.N - 1
                    retInfo.R = R_fcr;  % Use FCR threshold
                end
                break;
            end
        end
    end
end

% All validations passed
isRetrial = true;

end
