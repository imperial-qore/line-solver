%{ @file fes_validate.m
 %  @brief Validates inputs for Flow-Equivalent Server (FES) aggregation
 %
 %  @author LINE Development Team
%}

%{
 % @brief Validates inputs for Flow-Equivalent Server (FES) aggregation
 %
 % @details
 % Checks that the network structure and station subset are valid for FES aggregation:
 % - Model is a closed product-form queueing network
 % - Station subset is non-empty and not the entire network
 % - Subset contains only Queue or Delay stations
 % - All station indices are valid
 %
 % @par Syntax:
 % @code
 % [isValid, errorMsg] = fes_validate(sn, subsetIndices)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure (from model.getStruct())
 % <tr><td>subsetIndices<td>Array of station indices to aggregate (1-based)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>isValid<td>True if inputs are valid for FES aggregation
 % <tr><td>errorMsg<td>Error message if validation fails, empty otherwise
 % </table>
%}
function [isValid, errorMsg] = fes_validate(sn, subsetIndices)

isValid = false;
errorMsg = '';

% Check sn is a struct
if ~isstruct(sn)
    errorMsg = 'First argument must be a network structure (sn).';
    return;
end

% Check required fields exist in sn
requiredFields = {'nstations', 'nclasses', 'njobs', 'nodetype', 'stationToNode'};
for i = 1:length(requiredFields)
    if ~isfield(sn, requiredFields{i})
        errorMsg = sprintf('Network structure missing required field: %s', requiredFields{i});
        return;
    end
end

% Check subsetIndices is a numeric array
if ~isnumeric(subsetIndices)
    errorMsg = 'Station subset must be a numeric array of station indices.';
    return;
end

% Check subset is non-empty
if isempty(subsetIndices)
    errorMsg = 'Station subset cannot be empty.';
    return;
end

% Check model has only closed classes (FES applies to closed networks)
if sn_is_open_model(sn)
    errorMsg = 'FES aggregation only applies to closed queueing networks. Model contains open classes.';
    return;
end

M = sn.nstations;

% Check subset is not the entire network
if length(subsetIndices) >= M
    errorMsg = 'Cannot aggregate all stations. The subset must be a proper subset of the network.';
    return;
end

% Validate each station index in the subset
for i = 1:length(subsetIndices)
    stationIdx = subsetIndices(i);

    % Check station index is valid
    if stationIdx < 1 || stationIdx > M || stationIdx ~= floor(stationIdx)
        errorMsg = sprintf('Station index %d is invalid. Must be an integer between 1 and %d.', stationIdx, M);
        return;
    end

    % Check station is a Queue or Delay (not Fork, Join, Source, Sink, Cache, etc.)
    nodeIdx = sn.stationToNode(stationIdx);
    nodeType = sn.nodetype(nodeIdx);

    if nodeType ~= NodeType.Queue && nodeType ~= NodeType.Delay
        errorMsg = sprintf('Station %d has type %s. FES aggregation only supports Queue and Delay stations.', ...
            stationIdx, NodeType.toText(nodeType));
        return;
    end
end

% Check for duplicate stations in subset
if length(unique(subsetIndices)) ~= length(subsetIndices)
    errorMsg = 'Station subset contains duplicate stations.';
    return;
end

% All validations passed
isValid = true;

end
