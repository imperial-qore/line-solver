classdef WfCommonsLoader < handle
    % WfCommonsLoader - Load WfCommons JSON workflows into LINE Workflow objects.
    %
    % WfCommonsLoader provides static methods to load workflow traces from
    % the WfCommons format (https://github.com/wfcommons/workflow-schema) into
    % LINE Workflow objects for queueing analysis.
    %
    % Supported schema versions: 1.4, 1.5
    %
    % Example:
    %   wf = WfCommonsLoader.load('workflow.json');
    %   ph = wf.toPH();
    %
    % Example with options:
    %   options.distributionType = 'exp';
    %   options.defaultRuntime = 1.0;
    %   wf = WfCommonsLoader.load('workflow.json', options);
    %
    % options.resolveChildrenByName (default true): a child reference is looked up
    % among the task ids first and, failing that, among the task names. The
    % Pegasus traces of schema 1.4 (e.g. Montage) give each task an id such as
    % 'ID0000001' but list its children by name ('mBackground_ID0000013'), so an
    % id-only lookup drops every edge. A name shared by several tasks is
    % ambiguous and does not resolve. A reference that resolves to no task is
    % reported by line_warning. Set it false for the id-only lookup.
    %
    % Edges are the union of the `children` and `parents` lists, without duplicates, so a trace that
    % gives either form (the Montage dss trace lists parents only) yields the same graph. Parent
    % references resolve by the same id-then-unique-name rule.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        SUPPORTED_SCHEMA_VERSIONS = {'1.3', '1.4', '1.5'};
    end

    methods (Static)
        function wf = load(jsonFile, options)
            % LOAD Load a WfCommons JSON file into a Workflow object.
            %
            % WF = WFCOMMONSLOADER.LOAD(JSONFILE)
            % WF = WFCOMMONSLOADER.LOAD(JSONFILE, OPTIONS)
            %
            % Parameters:
            %   jsonFile - Path to WfCommons JSON file
            %   options  - Optional struct with:
            %     .distributionType - 'exp' (default), 'det', 'aph', 'hyperexp'
            %     .defaultSCV       - Default SCV for APH/HyperExp (default: 1.0)
            %     .defaultRuntime   - Default runtime when missing (default: 1.0)
            %     .useExecutionData - Use execution data if available (default: true)
            %     .storeMetadata    - Store WfCommons metadata (default: true)
            %
            % Returns:
            %   wf - Workflow object

            if nargin < 2
                options = struct();
            end
            options = WfCommonsLoader.parseOptions(options);

            % Read and parse JSON
            jsonStr = fileread(jsonFile);
            data = jsondecode(jsonStr);

            % Validate schema
            WfCommonsLoader.validateSchema(data);

            % Extract workflow name
            workflowName = WfCommonsLoader.extractName(data, jsonFile);

            % Build workflow
            wf = WfCommonsLoader.buildWorkflow(data, workflowName, options);
        end

        function wf = loadFromStruct(data, options)
            % LOADFROMSTRUCT Load from pre-parsed struct.
            %
            % WF = WFCOMMONSLOADER.LOADFROMSTRUCT(DATA)
            % WF = WFCOMMONSLOADER.LOADFROMSTRUCT(DATA, OPTIONS)
            %
            % Parameters:
            %   data    - Struct from jsondecode
            %   options - Options struct (see load method)
            %
            % Returns:
            %   wf - Workflow object

            if nargin < 2
                options = struct();
            end
            options = WfCommonsLoader.parseOptions(options);

            WfCommonsLoader.validateSchema(data);
            workflowName = WfCommonsLoader.extractName(data, 'struct_input');
            wf = WfCommonsLoader.buildWorkflow(data, workflowName, options);
        end

        function wf = loadFromUrl(urlString, options)
            % LOADFROMURL Load WfCommons JSON from a URL.
            %
            % WF = WFCOMMONSLOADER.LOADFROMURL(URLSTRING)
            % WF = WFCOMMONSLOADER.LOADFROMURL(URLSTRING, OPTIONS)
            %
            % Useful for loading workflows directly from repositories like
            % wfcommons/pegasus-instances.
            %
            % Parameters:
            %   urlString - URL pointing to WfCommons JSON file
            %   options   - Options struct (see load method)
            %
            % Returns:
            %   wf - Workflow object

            if nargin < 2
                options = struct();
            end
            options = WfCommonsLoader.parseOptions(options);

            % Fetch content from URL
            jsonStr = webread(urlString);

            % Parse JSON
            data = jsondecode(jsonStr);

            % Validate schema
            WfCommonsLoader.validateSchema(data);

            % Extract name from URL if not in JSON
            [~, defaultName, ~] = fileparts(urlString);

            % Build workflow
            workflowName = WfCommonsLoader.extractName(data, defaultName);
            wf = WfCommonsLoader.buildWorkflow(data, workflowName, options);
        end

        function isValid = validateFile(jsonFile)
            % VALIDATEFILE Check if file is valid WfCommons schema.
            %
            % ISVALID = WFCOMMONSLOADER.VALIDATEFILE(JSONFILE)
            %
            % Returns:
            %   isValid - true if file is valid WfCommons format

            try
                jsonStr = fileread(jsonFile);
                data = jsondecode(jsonStr);
                WfCommonsLoader.validateSchema(data);
                isValid = true;
            catch
                isValid = false;
            end
        end
    end

    methods (Static, Access = private)
        function options = parseOptions(options)
            % PARSEOPTIONS Set default option values.

            if ~isfield(options, 'distributionType')
                options.distributionType = 'exp';
            end
            if ~isfield(options, 'defaultSCV')
                options.defaultSCV = 1.0;
            end
            if ~isfield(options, 'defaultRuntime')
                options.defaultRuntime = 1.0;
            end
            if ~isfield(options, 'useExecutionData')
                options.useExecutionData = true;
            end
            if ~isfield(options, 'storeMetadata')
                options.storeMetadata = true;
            end
            if ~isfield(options, 'resolveChildrenByName')
                options.resolveChildrenByName = true;
            end
        end

        function validateSchema(data)
            % VALIDATESCHEMA Validate required fields in WfCommons schema.

            % Check schemaVersion
            if ~isfield(data, 'schemaVersion')
                line_error(mfilename, 'Missing schemaVersion field in WfCommons JSON.');
            end

            version = data.schemaVersion;
            if ~any(strcmp(version, WfCommonsLoader.SUPPORTED_SCHEMA_VERSIONS))
                line_warning(mfilename, 'Schema version %s may not be fully supported.', version);
            end

            % Check workflow structure
            if ~isfield(data, 'workflow')
                line_error(mfilename, 'Missing workflow field in WfCommons JSON.');
            end

            % Schema 1.4 and earlier use workflow.tasks directly
            % Schema 1.5+ uses workflow.specification.tasks
            hasTasks = false;
            if isfield(data.workflow, 'specification')
                if isfield(data.workflow.specification, 'tasks') && ~isempty(data.workflow.specification.tasks)
                    hasTasks = true;
                end
            elseif isfield(data.workflow, 'tasks') && ~isempty(data.workflow.tasks)
                % Schema 1.4 format: tasks directly under workflow
                hasTasks = true;
            end

            if ~hasTasks
                line_error(mfilename, 'Workflow must have at least one task.');
            end
        end

        function name = extractName(data, defaultName)
            % EXTRACTNAME Extract workflow name from data or use default.

            if isfield(data, 'name') && ~isempty(data.name)
                name = data.name;
            else
                [~, name, ~] = fileparts(defaultName);
            end
            % Sanitize name for MATLAB
            name = regexprep(name, '[^a-zA-Z0-9_]', '_');
            if isempty(name)
                name = 'Workflow';
            end
        end

        function wf = buildWorkflow(data, workflowName, options)
            % BUILDWORKFLOW Build Workflow object from parsed data.

            wf = Workflow(workflowName);

            % Detect schema format: 1.5+ uses specification.tasks, 1.4 uses tasks directly
            isLegacySchema = ~isfield(data.workflow, 'specification');
            if isLegacySchema
                % Schema 1.4 format: tasks directly under workflow with embedded runtime
                tasks = data.workflow.tasks;
            else
                % Schema 1.5+ format: tasks under specification
                spec = data.workflow.specification;
                tasks = spec.tasks;
            end

            % Build execution data map if available
            execMap = configureDictionary('string','cell');
            if options.useExecutionData
                if isLegacySchema
                    % Schema 1.4: runtime data is embedded in task objects
                    if isstruct(tasks)
                        for i = 1:length(tasks)
                            task = tasks(i);
                            taskId = WfCommonsLoader.getTaskId(task, i);
                            execMap{taskId} = task;
                        end
                    elseif iscell(tasks)
                        for i = 1:length(tasks)
                            task = tasks{i};
                            taskId = WfCommonsLoader.getTaskId(task, i);
                            execMap{taskId} = task;
                        end
                    end
                elseif isfield(data.workflow, 'execution')
                    % Schema 1.5+: separate execution section
                    exec = data.workflow.execution;
                    if isfield(exec, 'tasks')
                        execTasks = exec.tasks;
                        if isstruct(execTasks)
                            for i = 1:length(execTasks)
                                execMap{execTasks(i).id} = execTasks(i);
                            end
                        elseif iscell(execTasks)
                            for i = 1:length(execTasks)
                                execMap{execTasks{i}.id} = execTasks{i};
                            end
                        end
                    end
                end
            end

            % Phase 1: Create all activities
            taskMap = configureDictionary('string','cell');
            taskIdxMap = configureDictionary('string','double');
            n = length(tasks);

            for i = 1:n
                if isstruct(tasks)
                    task = tasks(i);
                else
                    task = tasks{i};
                end

                % Get task ID (schema 1.4 may use 'name' instead of 'id')
                taskId = WfCommonsLoader.getTaskId(task, i);

                % Get runtime from execution data
                runtime = options.defaultRuntime;
                if isKey(execMap, taskId)
                    execData = execMap{taskId};
                    if isfield(execData, 'runtimeInSeconds')
                        runtime = execData.runtimeInSeconds;
                    end
                end

                % Fit distribution
                dist = WfCommonsLoader.fitDistribution(runtime, options);

                % Create activity
                act = wf.addActivity(taskId, dist);
                taskMap{taskId} = act;
                taskIdxMap(taskId) = i;

                % Store metadata if requested
                if options.storeMetadata
                    act.metadata = WfCommonsLoader.extractMetadata(task, taskId, execMap);
                end
            end

            % Phase 2: Build adjacency structures
            [adjList, inDeg, outDeg, predList] = WfCommonsLoader.buildAdjacency(tasks, taskIdxMap, options);

            % Phase 3: Detect and add precedences
            WfCommonsLoader.addPrecedences(wf, tasks, taskMap, adjList, inDeg, outDeg, predList);
        end

        function dist = fitDistribution(runtime, options)
            % FITDISTRIBUTION Fit a distribution from runtime.

            if runtime <= GlobalConstants.FineTol
                dist = Immediate.getInstance();
            else
                switch lower(options.distributionType)
                    case 'exp'
                        dist = Exp.fitMean(runtime);
                    case 'det'
                        dist = Det(runtime);
                    case 'aph'
                        dist = APH.fitMeanAndSCV(runtime, options.defaultSCV);
                    case 'hyperexp'
                        if options.defaultSCV > 1.0
                            dist = HyperExp.fitMeanAndSCV(runtime, options.defaultSCV);
                        else
                            dist = Exp.fitMean(runtime);
                        end
                    otherwise
                        dist = Exp.fitMean(runtime);
                end
            end
        end

        function taskId = getTaskId(task, idx)
            % GETTASKID Extract task ID from task object.
            % Schema 1.4 may use 'name' as identifier, Schema 1.5+ uses 'id'
            if isfield(task, 'id')
                taskId = task.id;
            elseif isfield(task, 'name')
                taskId = task.name;
            else
                taskId = sprintf('task_%d', idx);
            end
        end

        function metadata = extractMetadata(task, taskId, execMap)
            % EXTRACTMETADATA Extract WfCommons metadata from task.

            metadata = struct();
            metadata.taskId = taskId;

            if isfield(task, 'name')
                metadata.name = task.name;
            end
            if isfield(task, 'inputFiles')
                metadata.inputFiles = task.inputFiles;
            end
            if isfield(task, 'outputFiles')
                metadata.outputFiles = task.outputFiles;
            end

            if isKey(execMap, taskId)
                exec = execMap{taskId};
                fields = {'executedAt', 'command', 'coreCount', 'avgCPU', ...
                          'readBytes', 'writtenBytes', 'memoryInBytes', ...
                          'energyInKWh', 'avgPowerInW', 'priority', 'machines'};
                for j = 1:length(fields)
                    f = fields{j};
                    if isfield(exec, f)
                        metadata.(f) = exec.(f);
                    end
                end
            end
        end

        function [adjList, inDeg, outDeg, predList] = buildAdjacency(tasks, taskIdxMap, options)
            % BUILDADJACENCY Build adjacency structures from tasks.
            % A child or parent resolves by task id, then (resolveChildrenByName) by a task name held by one task only.

            if nargin < 3
                options = struct('resolveChildrenByName', true);
            end
            n = length(tasks);
            nameIdxMap = configureDictionary('string','double');
            nameCount = configureDictionary('string','double');
            if options.resolveChildrenByName
                for i = 1:n
                    if isstruct(tasks)
                        task = tasks(i);
                    else
                        task = tasks{i};
                    end
                    if isfield(task, 'name') && (ischar(task.name) || isstring(task.name))
                        nm = char(task.name);
                        nameIdxMap(nm) = i;
                        if isKey(nameCount, nm)
                            nameCount(nm) = nameCount(nm) + 1;
                        else
                            nameCount(nm) = 1;
                        end
                    end
                end
            end
            unresolved = {};
            adjList = cell(n, 1);
            predList = cell(n, 1);
            inDeg = zeros(n, 1);
            outDeg = zeros(n, 1);

            for i = 1:n
                adjList{i} = [];
                predList{i} = [];
            end

            % Edges come from both `children` and `parents`, as a union without duplicates, so a trace that
            % lists either form (the Montage dss trace lists parents only) gives the same graph
            seen = false(n, n);
            for pass = 1:2
                if pass == 1, field = 'children'; else, field = 'parents'; end
                for i = 1:n
                    if isstruct(tasks)
                        task = tasks(i);
                    else
                        task = tasks{i};
                    end
                    if ~isfield(task, field) || isempty(task.(field))
                        continue
                    end
                    refs = task.(field);
                    if ischar(refs) || isstring(refs)
                        refs = cellstr(refs);
                    end
                    for j = 1:length(refs)
                        refId = char(refs{j});
                        refIdx = 0;
                        if isKey(taskIdxMap, refId)
                            refIdx = taskIdxMap(refId);
                        elseif options.resolveChildrenByName && isKey(nameIdxMap, refId) && nameCount(refId) == 1
                            refIdx = nameIdxMap(refId);
                        end
                        if refIdx == 0
                            unresolved{end+1} = refId; %#ok<AGROW>
                            continue
                        end
                        if pass == 1, from = i; to = refIdx; else, from = refIdx; to = i; end
                        if ~seen(from, to)
                            seen(from, to) = true;
                            adjList{from} = [adjList{from}, to];
                            predList{to} = [predList{to}, from];
                            outDeg(from) = outDeg(from) + 1;
                            inDeg(to) = inDeg(to) + 1;
                        end
                    end
                end
            end
            if ~isempty(unresolved)
                shown = unique(unresolved, 'stable');
                line_warning(mfilename, '%d child/parent reference(s) match no task id%s and were dropped: %s%s', ...
                    numel(unresolved), ternaryText(options.resolveChildrenByName, ' or unique task name', ''), ...
                    strjoin(shown(1:min(5, numel(shown))), ', '), ternaryText(numel(shown) > 5, ', ...', ''));
            end
        end

        function addPrecedences(wf, tasks, taskMap, adjList, inDeg, outDeg, predList)
            % ADDPRECEDENCES Detect patterns and add precedences to workflow.

            n = length(tasks);
            processed = false(n, n);  % Track processed edges

            % Step 1: Identify and process fork-join pairs
            for i = 1:n
                if outDeg(i) > 1
                    % Potential AND-fork
                    children = adjList{i};

                    % Find common join point
                    joinPoint = WfCommonsLoader.findCommonJoin(children, adjList, inDeg, n);

                    if ~isempty(joinPoint)
                        % Valid fork-join structure
                        if isstruct(tasks)
                            preTask = tasks(i);
                            postTask = tasks(joinPoint);
                        else
                            preTask = tasks{i};
                            postTask = tasks{joinPoint};
                        end

                        preTaskId = WfCommonsLoader.getTaskId(preTask, i);
                        preAct = taskMap{preTaskId};
                        postActs = cell(1, length(children));
                        for j = 1:length(children)
                            if isstruct(tasks)
                                childTask = tasks(children(j));
                            else
                                childTask = tasks{children(j)};
                            end
                            childTaskId = WfCommonsLoader.getTaskId(childTask, children(j));
                            postActs{j} = taskMap{childTaskId};
                        end
                        postTaskId = WfCommonsLoader.getTaskId(postTask, joinPoint);
                        postAct = taskMap{postTaskId};

                        % Add fork and join
                        wf.addPrecedence(Workflow.AndFork(preAct, postActs));
                        wf.addPrecedence(Workflow.AndJoin(postActs, postAct));

                        % Mark edges as processed
                        for j = 1:length(children)
                            processed(i, children(j)) = true;
                            processed(children(j), joinPoint) = true;
                        end
                    end
                end
            end

            % Step 2: Add remaining edges as serial connections
            for i = 1:n
                for j = adjList{i}
                    if ~processed(i, j)
                        if isstruct(tasks)
                            preTask = tasks(i);
                            postTask = tasks(j);
                        else
                            preTask = tasks{i};
                            postTask = tasks{j};
                        end

                        preTaskId = WfCommonsLoader.getTaskId(preTask, i);
                        postTaskId = WfCommonsLoader.getTaskId(postTask, j);
                        preAct = taskMap{preTaskId};
                        postAct = taskMap{postTaskId};
                        wf.addPrecedence(Workflow.Serial(preAct, postAct));
                    end
                end
            end
        end

        function joinPoint = findCommonJoin(children, adjList, inDeg, n)
            % FINDCOMMONJOIN Find a common join point for all children.

            if isempty(children)
                joinPoint = [];
                return;
            end

            % Use BFS from each child to find reachable nodes
            reachable = cell(1, length(children));
            for i = 1:length(children)
                reachable{i} = WfCommonsLoader.getReachableNodes(children(i), adjList, n);
            end

            % Find intersection
            common = reachable{1};
            for i = 2:length(children)
                common = intersect(common, reachable{i});
            end

            if isempty(common)
                joinPoint = [];
                return;
            end

            % Find the first common node that has all children as predecessors
            % (i.e., in-degree matches number of children from this fork)
            for node = common(:)'
                if inDeg(node) >= length(children)
                    % Check if all children can reach this node directly
                    allDirectChild = true;
                    for i = 1:length(children)
                        if ~ismember(node, adjList{children(i)})
                            allDirectChild = false;
                            break;
                        end
                    end
                    if allDirectChild
                        joinPoint = node;
                        return;
                    end
                end
            end

            joinPoint = [];
        end

        function reachable = getReachableNodes(start, adjList, n)
            % GETREACHABLENODES BFS to find all reachable nodes.

            visited = false(1, n);
            queue = start;
            visited(start) = true;
            reachable = [];

            while ~isempty(queue)
                current = queue(1);
                queue(1) = [];

                for next = adjList{current}
                    if ~visited(next)
                        visited(next) = true;
                        queue = [queue, next];
                        reachable = [reachable, next];
                    end
                end
            end
        end
    end
end

function s = ternaryText(cond, a, b)
% S = TERNARYTEXT(COND, A, B) A when COND holds, else B
if cond
    s = a;
else
    s = b;
end
end
