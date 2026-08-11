classdef Ensemble < Model
    % A model defined by a collection of sub-models.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        ensemble;
    end
    
    methods
        function self = Ensemble(models)
            % SELF = ENSEMBLE(MODELS)
            
            self@Model('Ensemble');
            self.ensemble = reshape(models,1,numel(models)); % flatten
        end
        
        function self = setEnsemble(self,ensemble)
            % SELF = SETENSEMBLE(SELF,ENSEMBLE)
            
            self.ensemble = ensemble;
        end
        
        function ensemble = getEnsemble(self)
            % ENSEMBLE = GETENSEMBLE()
            
            ensemble = self.ensemble;
        end
        
        function model = getModel(self, modelIdx)
            model = self.ensemble{modelIdx};
        end
    end
    
    methods(Access = protected)
        % Override copyElement method:
        function clone = copyElement(self)
            % CLONE = COPYELEMENT()

            % Make a shallow copy of all properties
            clone = copyElement@Copyable(self);
            % Make a deep copy of each ensemble object
            for e=1:length(self.ensemble)
                clone.ensemble{e} = copy(self.ensemble{e});
            end
        end
    end

    methods(Static)
        function unionNetwork = merge(ensemble)
            % MERGE Create a union Network from all Networks in an ensemble
            %
            % UNIONNETWORK = MERGE(ENSEMBLE) returns a single Network containing
            % all nodes and classes from each Network in the ensemble as
            % disconnected subnetworks. Node and class names are prefixed with
            % their originating model name to avoid collisions.
            %
            % Input:
            %   ensemble - Ensemble object or cell array of Network objects
            %
            % Output:
            %   unionNetwork - Network object containing merged subnetworks

            % Handle input - accept Ensemble object or cell array
            if isa(ensemble, 'Ensemble')
                models = ensemble.getEnsemble();
            elseif iscell(ensemble)
                models = ensemble;
            else
                line_error(mfilename, 'Input must be an Ensemble or cell array of Networks');
            end

            % Edge case: empty ensemble
            if isempty(models)
                line_warning(mfilename, 'Empty ensemble provided, returning empty Network');
                unionNetwork = Network('EmptyMergedNetwork');
                return;
            end

            % Create union network with combined name
            modelNames = cellfun(@(m) m.getName(), models, 'UniformOutput', false);
            unionName = strjoin(modelNames, '_');
            if length(unionName) > 50
                unionName = sprintf('MergedNetwork_%d', length(models));
            end
            unionNetwork = Network(unionName);

            % Check if any model has open classes (needs Source/Sink)
            hasOpenClasses = false;
            for m = 1:length(models)
                if models{m}.hasOpenClasses()
                    hasOpenClasses = true;
                    break;
                end
            end

            % Create single Source and Sink if needed
            unionSource = [];
            unionSink = [];
            if hasOpenClasses
                unionSource = Source(unionNetwork, 'MergedSource');
                unionSink = Sink(unionNetwork, 'MergedSink');
            end

            % Storage for mapping old nodes/classes to new ones
            nodeMap = cell(1, length(models));   % containers.Map for each model
            classMap = cell(1, length(models));  % containers.Map for each model
            joinPendingList = {};  % Join nodes need Forks to exist first

            % Phase 1: Create nodes for each model
            for m = 1:length(models)
                model = models{m};
                modelName = model.getName();
                prefix = [modelName, '_'];

                nodes = model.getNodes();
                nodeMap{m} = containers.Map('KeyType', 'char', 'ValueType', 'any');

                for n = 1:length(nodes)
                    oldNode = nodes{n};
                    newName = [prefix, oldNode.getName()];

                    if isa(oldNode, 'Source')
                        % Map to merged source
                        nodeMap{m}(oldNode.getName()) = unionSource;
                        continue;
                    elseif isa(oldNode, 'Sink')
                        % Map to merged sink
                        nodeMap{m}(oldNode.getName()) = unionSink;
                        continue;
                    elseif isa(oldNode, 'Delay')
                        newNode = Delay(unionNetwork, newName);
                    elseif isa(oldNode, 'Queue')
                        newNode = Queue(unionNetwork, newName, oldNode.schedStrategy);
                        if ~isinf(oldNode.numberOfServers) && oldNode.numberOfServers > 1
                            newNode.setNumberOfServers(oldNode.numberOfServers);
                        end
                        if ~isempty(oldNode.lldScaling)
                            newNode.setLoadDependence(oldNode.lldScaling);
                        end
                        if ~isinf(oldNode.cap)
                            newNode.setCapacity(oldNode.cap);
                        end
                    elseif isa(oldNode, 'Router')
                        newNode = Router(unionNetwork, newName);
                    elseif isa(oldNode, 'ClassSwitch')
                        % Skip auto-added ClassSwitch nodes (recreated by link())
                        if isprop(oldNode, 'autoAdded') && oldNode.autoAdded
                            continue;
                        end
                        newNode = ClassSwitch(unionNetwork, newName);
                    elseif isa(oldNode, 'Fork')
                        newNode = Fork(unionNetwork, newName);
                        if ~isempty(oldNode.output) && isprop(oldNode.output, 'tasksPerLink') && ~isempty(oldNode.output.tasksPerLink)
                            newNode.setTasksPerLink(oldNode.output.tasksPerLink);
                        end
                    elseif isa(oldNode, 'Join')
                        % Queue for later processing after Forks exist
                        joinPendingList{end+1} = struct('modelIdx', m, ...
                            'oldNode', oldNode, 'newName', newName);
                        continue;
                    elseif isa(oldNode, 'Logger')
                        newNode = Logger(unionNetwork, newName);
                    elseif isa(oldNode, 'Cache')
                        newNode = Cache(unionNetwork, newName, ...
                            oldNode.items.nitems, oldNode.itemLevelCap, ...
                            oldNode.replacementStrategy);
                    elseif isa(oldNode, 'Place')
                        newNode = Place(unionNetwork, newName);
                    elseif isa(oldNode, 'Transition')
                        newNode = Transition(unionNetwork, newName);
                    else
                        line_warning(mfilename, ...
                            sprintf('Unsupported node type: %s, skipping', class(oldNode)));
                        continue;
                    end

                    nodeMap{m}(oldNode.getName()) = newNode;
                end
            end

            % Phase 2: Process pending Join nodes (after Forks exist)
            for j = 1:length(joinPendingList)
                item = joinPendingList{j};
                m = item.modelIdx;
                oldNode = item.oldNode;
                newName = item.newName;

                if ~isempty(oldNode.joinOf)
                    forkName = oldNode.joinOf.getName();
                    if isKey(nodeMap{m}, forkName)
                        forkNode = nodeMap{m}(forkName);
                        newNode = Join(unionNetwork, newName, forkNode);
                    else
                        newNode = Join(unionNetwork, newName);
                    end
                else
                    newNode = Join(unionNetwork, newName);
                end
                nodeMap{m}(oldNode.getName()) = newNode;
            end

            % Phase 3: Create job classes for each model
            for m = 1:length(models)
                model = models{m};
                modelName = model.getName();
                prefix = [modelName, '_'];

                classes = model.getClasses();
                classMap{m} = containers.Map('KeyType', 'char', 'ValueType', 'any');

                for c = 1:length(classes)
                    oldClass = classes{c};
                    newName = [prefix, oldClass.getName()];

                    if isa(oldClass, 'OpenClass')
                        newClass = OpenClass(unionNetwork, newName, oldClass.priority);
                    elseif isa(oldClass, 'SelfLoopingClass')
                        % Map reference station to new node
                        oldRefStatName = oldClass.refstat.getName();
                        if isKey(nodeMap{m}, oldRefStatName)
                            newRefStat = nodeMap{m}(oldRefStatName);
                        else
                            line_error(mfilename, ...
                                sprintf('Reference station %s not found for class %s', oldRefStatName, oldClass.getName()));
                        end
                        newClass = SelfLoopingClass(unionNetwork, newName, ...
                            oldClass.population, newRefStat, oldClass.priority);
                    elseif isa(oldClass, 'ClosedClass')
                        % Map reference station to new node
                        oldRefStatName = oldClass.refstat.getName();
                        if isKey(nodeMap{m}, oldRefStatName)
                            newRefStat = nodeMap{m}(oldRefStatName);
                        else
                            line_error(mfilename, ...
                                sprintf('Reference station %s not found for class %s', oldRefStatName, oldClass.getName()));
                        end
                        newClass = ClosedClass(unionNetwork, newName, ...
                            oldClass.population, newRefStat, oldClass.priority);
                    else
                        line_warning(mfilename, ...
                            sprintf('Unknown class type: %s, skipping', class(oldClass)));
                        continue;
                    end

                    classMap{m}(oldClass.getName()) = newClass;
                end
            end

            % Phase 4: Set service and arrival distributions
            for m = 1:length(models)
                model = models{m};
                nodes = model.getNodes();
                classes = model.getClasses();

                for n = 1:length(nodes)
                    oldNode = nodes{n};

                    % Handle Source arrivals
                    if isa(oldNode, 'Source')
                        for c = 1:length(classes)
                            oldClass = classes{c};
                            if isa(oldClass, 'OpenClass')
                                if isKey(classMap{m}, oldClass.getName())
                                    newClass = classMap{m}(oldClass.getName());
                                    if ~isempty(oldNode.arrivalProcess) && ...
                                       length(oldNode.arrivalProcess) >= oldClass.index && ...
                                       ~isempty(oldNode.arrivalProcess{oldClass.index})
                                        dist = oldNode.arrivalProcess{oldClass.index};
                                        if ~isa(dist, 'Disabled')
                                            unionSource.setArrival(newClass, copy(dist));
                                        end
                                    end
                                end
                            end
                        end
                        continue;
                    end

                    % Skip Sink
                    if isa(oldNode, 'Sink')
                        continue;
                    end

                    % Skip auto-added ClassSwitch
                    if isa(oldNode, 'ClassSwitch') && isprop(oldNode, 'autoAdded') && oldNode.autoAdded
                        continue;
                    end

                    % Get mapped node
                    if ~isKey(nodeMap{m}, oldNode.getName())
                        continue;
                    end
                    newNode = nodeMap{m}(oldNode.getName());

                    % Set service distributions for Queue/Delay
                    if isa(oldNode, 'Queue') || isa(oldNode, 'Delay')
                        for c = 1:length(classes)
                            oldClass = classes{c};
                            if ~isKey(classMap{m}, oldClass.getName())
                                continue;
                            end
                            newClass = classMap{m}(oldClass.getName());

                            try
                                dist = oldNode.getService(oldClass);
                                if ~isempty(dist)
                                    % Copy all services including Disabled to ensure proper configuration
                                    newNode.setService(newClass, copy(dist));
                                end
                            catch
                                % Service not defined for this class, skip
                            end
                        end
                    end
                end
            end

            % Phase 5: Build routing matrix from linked routing matrices
            P = unionNetwork.initRoutingMatrix();
            nUnionClasses = unionNetwork.getNumberOfClasses();
            nUnionNodes = unionNetwork.getNumberOfNodes();

            % Initialize P cells if needed
            for r = 1:nUnionClasses
                for s = 1:nUnionClasses
                    if isempty(P{r,s})
                        P{r,s} = zeros(nUnionNodes);
                    end
                end
            end

            for m = 1:length(models)
                thisModel = models{m};
                nodes = thisModel.getNodes();
                classes = thisModel.getClasses();

                % Get the linked routing matrix from this model
                % This contains the full inter-class routing that was passed to link()
                Pm = [];
                try
                    Pm = thisModel.getLinkedRoutingMatrix();
                catch
                    % If model was never linked, try to get struct and use rtorig
                    try
                        sn = thisModel.getStruct(false);
                        Pm = sn.rtorig;
                    catch
                        Pm = [];
                    end
                end

                if isempty(Pm)
                    % Fall back to outputStrategy-based extraction if no linked matrix
                    for n = 1:length(nodes)
                        oldNode = nodes{n};

                        % Skip nodes without routing output
                        if ~isprop(oldNode, 'output') || isempty(oldNode.output) || ...
                           ~isprop(oldNode.output, 'outputStrategy')
                            continue;
                        end

                        % Skip auto-added ClassSwitch
                        if isa(oldNode, 'ClassSwitch') && isprop(oldNode, 'autoAdded') && oldNode.autoAdded
                            continue;
                        end

                        % Get mapped source node
                        if ~isKey(nodeMap{m}, oldNode.getName())
                            continue;
                        end
                        newSrcNode = nodeMap{m}(oldNode.getName());

                        for c = 1:length(classes)
                            oldClass = classes{c};

                            if ~isKey(classMap{m}, oldClass.getName())
                                continue;
                            end

                            if length(oldNode.output.outputStrategy) < oldClass.index || ...
                               isempty(oldNode.output.outputStrategy{oldClass.index})
                                continue;
                            end

                            newClass = classMap{m}(oldClass.getName());
                            strategy = oldNode.output.outputStrategy{oldClass.index};

                            % Check for probabilistic routing with destinations
                            if length(strategy) >= 3 && ~isempty(strategy{3})
                                probs = strategy{3};

                                for p = 1:length(probs)
                                    destNode = probs{p}{1};
                                    prob = probs{p}{2};

                                    if isKey(nodeMap{m}, destNode.getName())
                                        newDestNode = nodeMap{m}(destNode.getName());

                                        % Set probability in routing matrix
                                        % P{r,s}(i,j) - from class r at node i to class s at node j
                                        P{newClass.index, newClass.index}(...
                                            newSrcNode.index, newDestNode.index) = prob;
                                    end
                                end
                            end
                        end
                    end
                else
                    % Use the linked routing matrix (includes inter-class routing)
                    nModelClasses = length(classes);
                    nModelNodes = length(nodes);

                    for r = 1:nModelClasses
                        oldClassR = classes{r};
                        if ~isKey(classMap{m}, oldClassR.getName())
                            continue;
                        end
                        newClassR = classMap{m}(oldClassR.getName());

                        for s = 1:nModelClasses
                            oldClassS = classes{s};
                            if ~isKey(classMap{m}, oldClassS.getName())
                                continue;
                            end
                            newClassS = classMap{m}(oldClassS.getName());

                            % Check if this routing block exists in original matrix
                            if size(Pm,1) < r || size(Pm,2) < s || isempty(Pm{r,s})
                                continue;
                            end

                            Prs = Pm{r,s};

                            % Copy routing probabilities, mapping old node indices to new
                            for i = 1:min(nModelNodes, size(Prs,1))
                                oldNodeI = nodes{i};

                                % Skip auto-added ClassSwitch nodes
                                if isa(oldNodeI, 'ClassSwitch') && isprop(oldNodeI, 'autoAdded') && oldNodeI.autoAdded
                                    continue;
                                end

                                if ~isKey(nodeMap{m}, oldNodeI.getName())
                                    continue;
                                end
                                newNodeI = nodeMap{m}(oldNodeI.getName());

                                for j = 1:min(nModelNodes, size(Prs,2))
                                    if Prs(i,j) == 0
                                        continue;
                                    end

                                    oldNodeJ = nodes{j};

                                    % Skip auto-added ClassSwitch nodes
                                    if isa(oldNodeJ, 'ClassSwitch') && isprop(oldNodeJ, 'autoAdded') && oldNodeJ.autoAdded
                                        continue;
                                    end

                                    if ~isKey(nodeMap{m}, oldNodeJ.getName())
                                        continue;
                                    end
                                    newNodeJ = nodeMap{m}(oldNodeJ.getName());

                                    P{newClassR.index, newClassS.index}(newNodeI.index, newNodeJ.index) = Prs(i,j);
                                end
                            end
                        end
                    end
                end
            end

            % Disable checks before linking - merged networks from LQN ensembles may have
            % Delay nodes with all-Disabled services (e.g., async-only layers) which would
            % otherwise fail sanitize validation. This mirrors buildLayersRecursive behavior.
            unionNetwork.setChecks(false);

            % Link the network with the routing matrix
            unionNetwork.link(P);
        end

        function unionNetwork = toNetwork(ensemble)
            % TONETWORK Alias for merge
            %
            % UNIONNETWORK = TONETWORK(ENSEMBLE) - see MERGE

            unionNetwork = Ensemble.merge(ensemble);
        end
    end
end
