%{ @file sn_print.m
 %  @brief Prints comprehensive information about a NetworkStruct object
 %
 %  @author LINE Development Team
%}

%{
 % @brief Prints comprehensive information about a NetworkStruct object
 %
 % @details
 % This function displays all fields, matrices, lists, and maps in a formatted
 % manner useful for debugging and inspection of network structures.
 %
 % @par Syntax:
 % @code
 % sn_print(sn)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure to inspect
 % </table>
%}
function sn_print(sn)

% Suppress warnings about struct on objects
warning('off', 'MATLAB:structOnObject');

% Basic integer fields
fprintf('nstations: %d\n', sn.nstations);
fprintf('nstateful: %d\n', sn.nstateful);
fprintf('nnodes: %d\n', sn.nnodes);
fprintf('nclasses: %d\n', sn.nclasses);
fprintf('nclosedjobs: %d\n', sn.nclosedjobs);
fprintf('nchains: %d\n', sn.nchains);

% All matrix fields
printMatrix('refstat', sn.refstat);
printMatrix('njobs', sn.njobs);
printMatrix('nservers', sn.nservers);
printMatrix('connmatrix', sn.connmatrix);
printMatrix('scv', sn.scv);
printMatrix('isstation', sn.isstation);
printMatrix('isstateful', sn.isstateful);
printMatrix('isstatedep', sn.isstatedep);
printMatrix('nodeToStateful', sn.nodeToStateful);
printMatrix('nodeToStation', sn.nodeToStation);
printMatrix('stationToNode', sn.stationToNode);
printMatrix('stationToStateful', sn.stationToStateful);
printMatrix('statefulToStation', sn.statefulToStation);
printMatrix('statefulToNode', sn.statefulToNode);
printMatrix('rates', sn.rates);
printMatrix('classprio', sn.classprio);
printMatrix('phases', sn.phases);
printMatrix('phasessz', sn.phasessz);
printMatrix('phaseshift', sn.phaseshift);
printMatrix('schedparam', sn.schedparam);
printMatrix('chains', sn.chains);
printMatrix('rt', sn.rt);
printMatrix('nvars', sn.nvars);
printMatrix('rtnodes', sn.rtnodes);
printMatrix('csmask', sn.csmask);
printMatrix('isslc', sn.isslc);
printMatrix('cap', sn.cap);
printMatrix('classcap', sn.classcap);
printMatrix('refclass', sn.refclass);
printMatrix('lldscaling', sn.lldscaling);
printMatrix('fj', sn.fj);

% List fields with content
printNodeTypeList('nodetype', sn.nodetype);
if isfield(sn, 'classnames')
    if iscell(sn.classnames)
        printList('classnames', sn.classnames);
    else
        % Handle single string classname
        fprintf('classnames: ["%s"]\n', sn.classnames);
    end
end
if isfield(sn, 'nodenames')
    printList('nodenames', sn.nodenames);
end

% Map fields with detailed contents
printMapIfExists(sn, 'rtorig');
%printMapIfExists(sn, 'lst');
printMapIfExists(sn, 'state');
printMapIfExists(sn, 'stateprior');
printMapIfExists(sn, 'space');
printMapIfExists(sn, 'routing');
printMapIfExists(sn, 'procid');
printMapIfExists(sn, 'mu');
printMapIfExists(sn, 'phi');
printMapIfExists(sn, 'proc');
printMapIfExists(sn, 'pie');
printMapIfExists(sn, 'sched');
printMapIfExists(sn, 'inchain');
printMapIfExists(sn, 'visits');
printMapIfExists(sn, 'nodevisits');
printMapIfExists(sn, 'droprule');
printMapIfExists(sn, 'nodeparam');
% Special handling for sync field to avoid duplicate printing
if isfield(sn, 'sync')
    fprintf('sync: ');
    value = sn.sync;
    if iscell(value)
        % sync is a cell array of sync objects
        fprintf('{');
        for i = 1:length(value)
            if i > 1
                fprintf(', ');
            end
            fprintf('%d: ', i-1); % 0-based index
            
            syncObj = value{i};
            if isstruct(syncObj) && isfield(syncObj, 'active') && isfield(syncObj, 'passive')
                fprintf('{');
                % Print active events with indices
                fprintf('"active": {');
                if iscell(syncObj.active)
                    for j = 1:length(syncObj.active)
                        if j > 1
                            fprintf(', ');
                        end
                        fprintf('%d: ', j-1); % 0-based index
                        printValue(syncObj.active{j});
                    end
                else
                    fprintf('0: ');
                    printValue(syncObj.active);
                end
                fprintf('}');
                
                % Print passive events with indices
                fprintf(', "passive": {');
                if iscell(syncObj.passive)
                    for j = 1:length(syncObj.passive)
                        if j > 1
                            fprintf(', ');
                        end
                        fprintf('%d: ', j-1); % 0-based index
                        printValue(syncObj.passive{j});
                    end
                else
                    fprintf('0: ');
                    printValue(syncObj.passive);
                end
                fprintf('}');
                fprintf('}');
            else
                printValue(syncObj);
            end
        end
        fprintf('}\n');
    elseif isstruct(value)
        % Check if it's a map of sync objects
        fields = fieldnames(value);
        if ~isempty(fields)
            fprintf('{');
            for i = 1:length(fields)
                if i > 1
                    fprintf(', ');
                end
                key = fields{i};
                % Handle numeric keys like x0, x1
                if regexp(key, '^x\d+$')
                    fprintf('%s', key(2:end));
                else
                    fprintf('"%s"', key);
                end
                fprintf(': ');
                
                syncObj = value.(key);
                if isstruct(syncObj) && isfield(syncObj, 'active') && isfield(syncObj, 'passive')
                    fprintf('{');
                    % Print active events with indices
                    fprintf('"active": {');
                    if iscell(syncObj.active)
                        for j = 1:length(syncObj.active)
                            if j > 1
                                fprintf(', ');
                            end
                            fprintf('%d: ', j-1); % 0-based index
                            printValue(syncObj.active{j});
                        end
                    else
                        fprintf('0: ');
                        printValue(syncObj.active);
                    end
                    fprintf('}');
                    
                    % Print passive events with indices
                    fprintf(', "passive": {');
                    if iscell(syncObj.passive)
                        for j = 1:length(syncObj.passive)
                            if j > 1
                                fprintf(', ');
                            end
                            fprintf('%d: ', j-1); % 0-based index
                            printValue(syncObj.passive{j});
                        end
                    else
                        fprintf('0: ');
                        printValue(syncObj.passive);
                    end
                    fprintf('}');
                    fprintf('}');
                else
                    printValue(syncObj);
                end
            end
            fprintf('}\n');
        else
            fprintf('{}\n');
        end
    else
        printValue(value);
        fprintf('\n');
    end
else
    % Field doesn't exist, skip
end
printMapIfExists(sn, 'gsync');
printMapIfExists(sn, 'cdscaling');

% Restore warnings
warning('on', 'MATLAB:structOnObject');

end

function matStr = printMatrixCompact(matrix)
% Helper function to print matrix in compact format for maps
if isempty(matrix)
    matStr = '[]';
    return;
end

sb = '[';
[numRows, numCols] = size(matrix);
for i = 1:numRows
    if i > 1
        sb = [sb '; '];
    end
    for j = 1:numCols
        if j > 1
            sb = [sb ' '];
        end
        value = matrix(i, j);
        % Check if matrix contains only integers or if individual value is integer
        if (isnumeric(matrix) && all(matrix(:) == floor(matrix(:))) && all(~isinf(matrix(:)))) || ...
           (value == floor(value) && ~isinf(value))
            sb = [sb num2str(round(value))];
        else
            sb = [sb num2str(value)];
        end
    end
end
sb = [sb ']'];
matStr = sb;
end

function printNodeTypeList(name, nodeTypes)
% Helper function to print nodetype array with proper formatting
fprintf('%s: ', name);
if isempty(nodeTypes)
    fprintf('[]\n');
else
    fprintf('[');
    for i = 1:length(nodeTypes)
        if i > 1
            fprintf(', ');
        end
        fprintf('%s', NodeType.toText(nodeTypes(i)));
    end
    fprintf(']\n');
end
end

function printList(name, list)
% Helper function to print lists with proper formatting
fprintf('%s: ', name);
if isempty(list)
    fprintf('[]\n');
elseif iscell(list)
    fprintf('[');
    for i = 1:length(list)
        if i > 1
            fprintf(', ');
        end
        if ischar(list{i}) || isstring(list{i})
            fprintf('"%s"', char(list{i}));
        else
            fprintf('%s', mat2str(list{i}));
        end
    end
    fprintf(']\n');
elseif isnumeric(list)
    % Handle numeric arrays
    fprintf('%s\n', printMatrixCompact(list));
else
    % Handle other types
    fprintf('%s\n', mat2str(list));
end
end

function printMatrix(name, matrix)
% Helper function to print matrix with all values
if isempty(matrix)
    fprintf('%s: []\n', name);
elseif isnumeric(matrix) && all(isnan(matrix(:)))
    fprintf('%s: null\n', name);
else
    fprintf('%s: %s\n', name, printMatrixCompact(matrix));
end
end

function printMapIfExists(sn, fieldName)
% Helper function to print fields only if they exist
if isfield(sn, fieldName)
    value = sn.(fieldName);
    if isstruct(value)
        printMapContents(fieldName, value);
    elseif iscell(value)
        printCellContents(fieldName, value);
    else
        % For other types, use printMatrix
        printMatrix(fieldName, value);
    end
end
end

function printCellContents(name, cellArray)
% Helper function to print cell array contents
fprintf('%s: ', name);
if isempty(cellArray)
    fprintf('{}\n');
else
    % Check if this is a single-element cell array containing a matrix
    if numel(cellArray) == 1 && isnumeric(cellArray{1})
        % Just print the matrix directly without extra brackets
        fprintf('%s\n', printMatrixCompact(cellArray{1}));
        return;
    end
    
    % Check if this is a map-like cell array
    if numel(cellArray) > 0 && isstruct(cellArray{1})
        % Print as a map
        printMapContents(name, cellArray{1});
        return;
    end
    
    % Special handling for certain fields that are map-like
    if strcmp(name, 'state') || strcmp(name, 'stateprior') || strcmp(name, 'space') || ...
       strcmp(name, 'mu') || strcmp(name, 'phi') || strcmp(name, 'pie')
        % These fields use station names as keys
        fprintf('{');
        % Assume we have access to station names through evalin
        try
            nodenames = evalin('caller', 'sn.nodenames');
            for i = 1:numel(cellArray)
                if i > 1
                    fprintf(', ');
                end
                if i <= length(nodenames)
                    fprintf('"%s": ', nodenames{i});
                end
                printValue(cellArray{i});
            end
        catch
            % Fall back to simple array printing
            for i = 1:numel(cellArray)
                if i > 1
                    fprintf(', ');
                end
                printValue(cellArray{i});
            end
        end
        fprintf('}\n');
    else
        % Regular cell array - use [] for arrays
        fprintf('[');
        for i = 1:numel(cellArray)
            if i > 1
                fprintf(', ');
            end
            printValue(cellArray{i});
        end
        fprintf(']\n');
    end
end
end

function printMapContents(name, mapStruct)
% Helper function to print map contents with detailed structure
if isempty(mapStruct) || ~exist('mapStruct', 'var')
    if ~isempty(name)
        fprintf('%s: null\n', name);
    else
        fprintf('null\n');
    end
    return;
end

if ~isstruct(mapStruct)
    if ~isempty(name)
        fprintf('%s: null\n', name);
    else
        fprintf('null\n');
    end
    return;
end

fields = fieldnames(mapStruct);
if isempty(fields)
    if ~isempty(name)
        fprintf('%s: {}\n', name);
    else
        fprintf('{}\n');
    end
    return;
end

if ~isempty(name)
    fprintf('%s: {', name);
else
    fprintf('{');
end
first = true;
for i = 1:length(fields)
    if ~first
        fprintf(', ');
    end
    first = false;
    
    % Print key - handle different key types like the Kotlin version
    key = fields{i};
    % Check if key is numeric (like x0, x1, etc.)
    if regexp(key, '^x\d+$')
        % Extract the numeric part
        fprintf('%s', key(2:end));
    elseif ischar(key)
        fprintf('"%s"', key);
    else
        fprintf('"%s"', char(string(key)));
    end
    
    fprintf(': ');
    
    % Print value
    value = mapStruct.(key);
    printValue(value);
end
fprintf('}\n');
end

function printValue(value)
% Helper function to print individual values with appropriate formatting
if isempty(value)
    fprintf('null');
elseif isa(value, 'jline.lang.Event') || (isjava(value) && contains(class(value), 'Event'))
    % Handle Java Event objects
    fprintf('{');
    printed = false;
    
    % Try to print each field
    try
        nodeVal = value.getNode();
        fprintf('"node": %d', nodeVal);
        printed = true;
    catch
        % Can't get node
    end
    
    try
        eventVal = value.getEvent();
        if printed
            fprintf(', ');
        end
        fprintf('"event": ');
        try
            fprintf('%d', eventVal.getId());
        catch
            fprintf('"%s"', char(eventVal.toString()));
        end
        printed = true;
    catch
        % Can't get event
    end
    
    try
        classVal = value.getJobClass();
        if printed
            fprintf(', ');
        end
        fprintf('"class": %d', classVal);
        printed = true;
    end
    
    try
        probVal = value.getProb();
        if printed
            fprintf(', ');
        end
        fprintf('"prob": ');
        if isnan(probVal)
            fprintf('NaN');
        else
            fprintf('%g', probVal);
        end
        printed = true;
    catch
        % Skip prob
    end
    
    try
        stateVal = value.getState();
        if printed
            fprintf(', ');
        end
        fprintf('"state": ');
        if isempty(stateVal) || (isa(stateVal, 'jline.util.matrix.Matrix') && stateVal.isEmpty())
            fprintf('[]');
        else
            printValue(stateVal);
        end
        printed = true;
    catch
        % Skip state
    end
    
    try
        tVal = value.getT();
        if printed
            fprintf(', ');
        end
        fprintf('"t": ');
        if isnan(tVal)
            fprintf('NaN');
        else
            fprintf('%g', tVal);
        end
        printed = true;
    catch
        % Skip t
    end
    
    try
        jobVal = value.getJob();
        if printed
            fprintf(', ');
        end
        fprintf('"job": ');
        if isnan(jobVal)
            fprintf('NaN');
        else
            fprintf('%g', jobVal);
        end
        printed = true;
    catch
        % Skip job
    end
    
    % If nothing was printed, show it's an Event
    if ~printed
        fprintf('<Event>');
    end
    
    fprintf('}');
elseif isstruct(value)
    fields = fieldnames(value);
    if isempty(fields)
        % Check if this might be a Java object that appears as empty struct
        className = class(value);
        if contains(className, 'Event')
            fprintf('{<Event object - fields not accessible>}');
        else
            fprintf('{}');
        end
    else
        fprintf('{');
        for i = 1:length(fields)
            if i > 1
                fprintf(', ');
            end
            fprintf('"%s": ', fields{i});
            printValue(value.(fields{i}));
        end
        fprintf('}');
    end
elseif isnumeric(value)
    if numel(value) > 1
        fprintf('%s', printMatrixCompact(value));
    else
        if value == floor(value) && ~isinf(value)
            fprintf('%d', round(value));
        else
            fprintf('%g', value);
        end
    end
elseif ischar(value) || isstring(value)
    fprintf('"%s"', char(value));
elseif islogical(value)
    if value
        fprintf('true');
    else
        fprintf('false');
    end
elseif iscell(value)
    % Check if this is a nested cell array (like proc field)
    if numel(value) > 0 && iscell(value{1})
        % Handle nested cell array as indexed map
        fprintf('{');
        for i = 1:numel(value)
            if i > 1
                fprintf(', ');
            end
            fprintf('[%d]: ', i-1);  % Use 0-based indexing for consistency with Java
            printValue(value{i});
        end
        fprintf('}');
    else
        % Regular cell array
        fprintf('[');
        for i = 1:length(value)
            if i > 1
                fprintf(', ');
            end
            printValue(value{i});
        end
        fprintf(']');
    end
elseif isa(value, 'function_handle')
    % Display function handles similar to Java lambda notation
    funcStr = func2str(value);
    if contains(funcStr, '@')
        funcStr = strrep(funcStr, '@', '');
    end
    fprintf('%s', funcStr);
else
    % Handle objects and other types
    try
        % Check if it's a Sync object
        if isa(value, 'Sync')
            fprintf('{');
            fprintf('"active": ');
            printValue(value.active);
            fprintf(', "passive": ');
            printValue(value.passive);
            fprintf('}');
        else
            str = char(string(value));
            % Check if this is a sync string that needs special handling
            if contains(str, 'sync:')
                % This is likely the sync field being printed as a string
                % Extract the content after "sync:"
                syncContent = strtrim(extractAfter(str, 'sync:'));
                fprintf('%s', syncContent);
            elseif contains(str, 'x') && contains(str, ' ')
                % Clean up MATLAB object notation (e.g., 1x1 ClosedClass)
                parts = split(str, ' ');
                if length(parts) >= 2
                    fprintf('%s', parts{end});
                else
                    fprintf('%s', str);
                end
            else
                fprintf('%s', str);
            end
        end
    catch
        % If conversion fails, just display the class name
        fprintf('<%s>', class(value));
    end
end
end