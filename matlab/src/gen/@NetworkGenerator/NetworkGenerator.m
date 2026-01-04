classdef NetworkGenerator < handle
    % A generator object that generates queueing network models
    % based on user specification. Characteristics of generated
    % models can be configured via the generator's properties. See
    % user guide in report for detailed usage instructions.

    properties
        schedStrat char
        routingStrat char
        distribution char
        cclassJobLoad char
        hasVaryingServiceRates (1, 1) logical
        hasMultiServerQueues (1, 1) logical
        hasRandomCSNodes (1, 1) logical
        hasMultiChainCS (1, 1) logical
        topologyFcn function_handle
    end

    properties (SetAccess = private)
        maxServers = 40
        highJobLoadRange = [31 40]
        medJobLoadRange = [11 20]
        lowJobLoadRange = [1 5]
    end

    methods
        function obj = NetworkGenerator(varargin)
            p = inputParser;
            addParameter(p, 'schedStrat', 'randomize'); % alternatives: 'fcfs', 'ps'
            addParameter(p, 'routingStrat', 'randomize'); % alternatives: 'Random', 'Probabilities'
            addParameter(p, 'distribution', 'randomize'); % alternatives: 'exp', 'erlang', 'hyperexp'
            addParameter(p, 'cclassJobLoad', 'randomize'); % alternatives: 'low', 'medium', 'high'
            addParameter(p, 'hasVaryingServiceRates', true);
            addParameter(p, 'hasMultiServerQueues', true);
            addParameter(p, 'hasRandomCSNodes', true);
            addParameter(p, 'hasMultiChainCS', true);
            addParameter(p, 'topologyFcn', @NetworkGenerator.randGraph);
            parse(p, varargin{:});

            obj.schedStrat = p.Results.schedStrat;
            obj.routingStrat = p.Results.routingStrat;
            obj.distribution = p.Results.distribution;
            obj.cclassJobLoad = p.Results.cclassJobLoad;
            obj.hasVaryingServiceRates = p.Results.hasVaryingServiceRates;
            obj.hasMultiServerQueues = p.Results.hasMultiServerQueues;
            obj.hasRandomCSNodes = p.Results.hasRandomCSNodes;
            obj.hasMultiChainCS = p.Results.hasMultiChainCS;
            obj.topologyFcn = p.Results.topologyFcn;
        end

        %% Main function to call.
        % Returns a generated QN model according to
        % specified properties of the NetworkGenerator object
        function model = generate(obj, numQueues, numDelays, numOClass, numCClass)
            if nargin < 2 || isempty(numQueues)
                numQueues = randi([1, 8]);
            end
            if nargin < 3 || isempty(numDelays)
                if numQueues > 1
                    numDelays = randi([0,1]);
                else
                    numDelays = 1;
                end
            end
            if nargin < 4 || isempty(numOClass)
                numOClass = 0;
            end
            if nargin < 5 || isempty(numCClass)
                numCClass = randi([1, 4]);
            end
            obj.validateArgs(numQueues, numDelays, numOClass, numCClass);
            model = Network('nw');
            stations = obj.createStations(model, numQueues, numDelays, numOClass);
            classes = obj.createClasses(model, numOClass, numCClass);
            obj.setServiceProcesses(stations, classes);
            obj.defineTopology(model, obj.topologyFcn(length(stations)));
            %model = NetworkGenerator.initDefaultCustom(model);
        end

        function set.schedStrat(obj, strat)
            if (strcmp(strat, 'fcfs') ... % All stations have FCFS
                    || strcmp(strat, 'ps') ... % All stations have PS
                    || strcmp(strat, 'inf') ...
                    || strcmp(strat, 'lcfs') ...
                    || strcmp(strat, 'lcfspr') ...
                    || strcmp(strat, 'siro') ...
                    || strcmp(strat, 'sjf') ...
                    || strcmp(strat, 'ljf') ...
                    || strcmp(strat, 'sept') ...
                    || strcmp(strat, 'lept') ...
                    || strcmp(strat, 'randomize')) % Randomize across stations
                obj.schedStrat = strat;
            else
                error('NG:schedStrat', ...
                    'Scheduling strategy does not exist or is not supported');
            end
        end

        function set.routingStrat(obj, strat)
            if (strcmp(strat, 'Probabilities') ... % All probabilistic routing
                    || strcmp(strat, 'Random') ... % All random routing
                    || strcmp(strat, 'randomize')) % Randomize across all classes and stations
                obj.routingStrat = strat;
            else
                error('NG:routingStrat', ...
                    'Routing strategy does not exist or is not supported');
            end
        end

        function set.distribution(obj, distrib)
            if (strcmpi(distrib, 'Exp') ... % All exponentially distributed service distributions
                    || strcmpi(distrib, 'HyperExp') ... % All hyperexponentially distributed
                    || strcmpi(distrib, 'Erlang') ... % All Erlang distributed
                    || strcmpi(distrib, 'randomize')) % Randomize across all classes and stations
                obj.distribution = distrib;
            else
                error('NG:distribution', ...
                    'Distribution does not exist or is not supported');
            end
        end

        function set.cclassJobLoad(obj, load)
            if (strcmpi(load, 'high') ... % All stations have high job load
                    || strcmpi(load, 'medium') ... % All stations have medium job load
                    || strcmpi(load, 'low') ... % All stations have low job load
                    || strcmpi(load, 'randomize')) % Randomize loads across stations
                obj.cclassJobLoad = load;
            else
                error('NG:cclassJobLoad', ...
                    'Model load can only be high, medium, or low');
            end
        end
        % Can be any function handle that takes an integer and returns a
        % digraph object
        function set.topologyFcn(obj, fcn)
            try
                graph = fcn(2);
            catch
                line_error(mfilename,'topologyFcn should take a positive integer and return a digraph');
            end

            if ~isa(graph, 'digraph') || numnodes(graph) ~= 2
                line_error(mfilename,'topologyFcn should take a positive integer and return a digraph');
            end

            obj.topologyFcn = fcn;
        end

    end

    methods (Access = private)
        % Validates that parameter values for the network are sound
        function validateArgs(~, numQueues, numDelays, numOClass, numCClass)
            if numQueues < 0 || numDelays < 0 || numOClass < 0 || numCClass < 0
                error('NG:negativeArgs', 'Arguments must be non-negative');
            elseif numQueues + numDelays <= 0 || numOClass + numCClass <= 0
                error('NG:noStationsOrJobs', 'At least one station and one job class required');
            end
        end

        % Creates the stations in the network, including source and sink
        function stations = createStations(obj, model, numQueues, numDelays, hasOClass)
            queues = cell(numQueues, 1);
            for i = 1 : numQueues
                queues{i} = Queue(model, obj.name('queue', i), obj.chooseSchedStrat);
                queues{i}.setNumberOfServers(chooseNumServers);
            end

            delays = cell(numDelays, 1);
            for i = 1 : numDelays
                delays{i} = Delay(model, obj.name('delay', i));
            end

            if hasOClass
                Source(model, 'source');
                Sink(model, 'sink');
            end
            % Source and sink intentionally excluded from station list
            stations = [queues; delays];

            function n = chooseNumServers
                if obj.hasMultiServerQueues
                    n = randi(obj.maxServers);
                else
                    n = 1;
                end
            end
        end

        % Creates the job classes that are serviced within the network
        % Note: Arrival processes for open classes are set here for simplicity
        function classes = createClasses(obj, model, numOClass, numCClass)
            openClasses = cell(numOClass, 1);
            for i = 1 : numOClass
                openClasses{i} = OpenClass(model, obj.name('OClass', i));
                model.getSource.setArrival(openClasses{i}, obj.chooseDistribution);
            end

            closedClasses = cell(numCClass, 1);
            refStat = model.stations{randi(model.getNumberOfStations - (numOClass > 0))};
            for i = 1 : numCClass
                closedClasses{i} = ClosedClass(model, ...
                    obj.name('CClass', i), chooseNumJobs, refStat);
            end

            classes = [openClasses; closedClasses];

            function numJobs = chooseNumJobs
                switch lower(obj.cclassJobLoad)
                    case 'high'
                        numJobs = randi(obj.highJobLoadRange);
                    case 'medium'
                        numJobs = randi(obj.medJobLoadRange);
                    case 'low'
                        numJobs = randi(obj.lowJobLoadRange);
                    case 'randomize'
                        numJobs = randi(obj.highJobLoadRange(2));
                end
            end
        end

        % Sets the service processes for all job classes at all stations
        % Note: Stations here excludes the Source node for simplicity
        function setServiceProcesses(obj, stations, classes)
            for i = 1 : length(classes)
                for j = 1 : length(stations)
                    stations{j}.setService(classes{i}, obj.chooseDistribution);
                end
            end
        end

        % Defines the topology of the network, adding random class switching nodes
        function defineTopology(obj, model, topology)
            numStations = model.getNumberOfStations;
            csMask = obj.genCSMask(model);
            % Add source and sink to topology graph
            if model.hasOpenClasses
                topology = addnode(topology, 2);
                topology = addedge(topology, model.getIndexSourceNode, randi(numStations - 1));
                topology = addedge(topology, randi(numStations - 1), model.getIndexSinkNode);
            end
            for i = 1 : numStations
                [~, destIDs] = outedges(topology, i);
                outgoingNodes = obj.addOutgoingLinks(model, i, sort(destIDs), csMask);
                obj.setRoutingStrategies(model, model.stations{i}, outgoingNodes);
            end
        end

        % Creates a mask that indicates which classes can switch with each other
        function mask = genCSMask(obj, model)
            numOClasses = length(model.getIndexOpenClasses);
            numCClasses = length(model.getIndexClosedClasses);
            mask = zeros(numOClasses + numCClasses);

            if ~obj.hasMultiChainCS
                mask(1 : numOClasses, 1 : numOClasses) = 1;
                mask(numOClasses + 1 : end, numOClasses + 1 : end) = 1;
                mask = logical(mask);
                return;
            end

            openChains = [];
            closedChains = [];
            if model.hasOpenClasses
                openChains = assignChains(numOClasses, randi(numOClasses));
            end
            if model.hasClosedClasses
                closedChains = assignChains(numCClasses, randi(numCClasses));
            end
            allChains = [openChains closedChains];

            startIdx = 1;
            for i = 1 : length(allChains)
                endIdx = startIdx + allChains(i) - 1;
                mask(startIdx : endIdx, startIdx : endIdx) = 1;
                startIdx = endIdx + 1;
            end

            mask = logical(mask);

            function chains = assignChains(numClasses, numChains)
                chains = randintfixedsum(numClasses, numChains);
                chains = chains(randperm(numChains));

                function res = randintfixedsum(s, n)
                    if n == 1
                        res = s;
                        return;
                    elseif s == n
                        res = ones(1, n);
                        return;
                    end
                    first = randi(s - n);
                    res = [first randintfixedsum(s - first, n - 1)];
                end
            end
        end

        % For a single node, add all outgoing links and random class switch nodes
        function outgoingNodes = addOutgoingLinks(obj, model, sourceID, allDestIDs, mask)
            sourceNode = model.nodes{sourceID};
            outgoingNodes = cell(length(allDestIDs), 1);

            for i = 1 : length(allDestIDs)
                destNode = model.nodes{allDestIDs(i)};
                if obj.hasRandomCSNodes ...
                        && randi([0 1]) == 1 ...
                        && ~(isa(sourceNode, 'Source') || isa(destNode, 'Sink'))
                    cs = obj.randClassSwitchNode(model, mask, sourceNode, destNode);
                    outgoingNodes{i} = cs;
                    model.addLink(sourceNode, cs);
                    model.addLink(cs, destNode);
                    for j = 1 : model.getNumberOfClasses
                        cs.setProbRouting(model.classes{j}, destNode, 1);
                    end
                else
                    outgoingNodes{i} = destNode;
                    model.addLink(sourceNode, destNode);
                end
            end
        end

        % Instantiates a class switch node with a random switching matrix
        function cs = randClassSwitchNode(obj, model, mask, sourceNode, destNode)
            name = obj.name('cs', sourceNode.getName, destNode.getName);
            cs = ClassSwitch(model, name, randClassSwitchMatrix);

            function matrix = randClassSwitchMatrix
                matrix = zeros(length(mask));
                for i = 1 : length(mask)
                    matrix(i, mask(i, :)) = obj.randfixedsumone(nnz(mask(i, :)));
                end
            end
        end

        % Set routing strategies for each job class at a specified station
        function setRoutingStrategies(obj, model, station, outgoingNodes)
            if isa(station, 'Source')
                openClasses = model.classes(model.getIndexOpenClasses);
                for i = 1 : length(openClasses)
                    for j = 1 : length(outgoingNodes)
                        station.setProbRouting(openClasses{i}, outgoingNodes{j}, 1);
                    end
                end
                return;
            end

            classes = model.classes;
            for i = 1 : length(classes)
                strat = obj.chooseRoutingStrat;
                station.setRouting(classes{i}, RoutingStrategy.fromText(strat));
                if strcmp(strat, 'Probabilities')
                    if isa(classes{i}, 'ClosedClass') && isa(outgoingNodes{end}, 'Sink')
                        station.setProbRouting(classes{i}, outgoingNodes{end}, 0);
                        probs = obj.randfixedsumone(length(outgoingNodes) - 1);
                    else
                        probs = obj.randfixedsumone(length(outgoingNodes));
                    end
                    for j = 1 : length(probs)
                        station.setProbRouting(classes{i}, outgoingNodes{j}, probs(j));
                    end
                end
            end
        end

        % Generate random probabilities that sum up to one
        function probs = randfixedsumone(obj, numElems)
            probs = randfixedsum(numElems, 1, 1, 0, 1);
            probs = ceil(probs * 1000) / 1000;
            [~, maxIdx] = max(probs);
            probs(maxIdx) = probs(maxIdx) - (sum(probs) - 1);
        end

        % Creates a name string for network elements via concatenation
        function name = name(~, str, varargin)
            if strcmpi(str, 'cs')
                name = strcat(str, '_', varargin{1}, '_', varargin{2});
            else
                name = strcat(str, int2str(varargin{1}));
            end
        end

        % Returns a random scheduling strategy for a queueing station
        function strat = chooseSchedStrat(obj)
            if strcmpi(obj.schedStrat, 'randomize')
                id = randi(2);
                switch id
                    case 1
                        strat = SchedStrategy.FCFS;
                    case 2
                        strat = SchedStrategy.PS;
                end
            else
                strat = categorical({obj.schedStrat});
            end
        end

        % Returns a random routing strategy for a network node
        function strat = chooseRoutingStrat(obj)
            if strcmpi(obj.routingStrat, 'randomize')
                id = randi(2);
                switch id
                    case 1
                        strat = 'Random';
                    case 2
                        strat = 'Probabilities';
                end
            else
                strat = obj.routingStrat;
            end
        end

        % Returns a random distribution object for a service/arrival process
        function dist = chooseDistribution(obj)
            switch lower(obj.distribution)
                case 'randomize'
                    id = randi(3);
                case 'exp'
                    id = 1;
                case 'erlang'
                    id = 2;
                case 'hyperexp'
                    id = 3;
            end

            if obj.hasVaryingServiceRates
                mean = 2 ^ randi([-6 6]);
            else
                mean = 1;
            end

            switch id
                case 1
                    dist = Exp.fitMeanAndSCV(mean, 1);
                case 2
                    dist = Erlang.fitMeanAndSCV(mean, 1/(2^randi([0 6])));
                case 3
                    dist = HyperExp.fitMeanAndSCV(mean, 2^randi([0 6]));
            end
        end
    end

    methods (Static)
        graph = randGraph(numVertices);
        model = initDefaultCustom(model, nodes);
    end
end