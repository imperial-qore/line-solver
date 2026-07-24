classdef LayeredNetworkGenerator < handle
    % A generator object that generates layered queueing network models
    % based on user specification. Characteristics of generated
    % models can be configured via the generator's properties. See
    % user guide in report for detailed usage instructions.
    
    properties
        populationRange (1, 2) double
        thinkTimeRange (1, 2) double
        taskInfProbability (1, 1) double
        procInfProbability (1, 1) double
        taskMultiRange (1, 2) double
        procMultiRange (1, 2) double
        hostDemandRange (1, 2) double
        synchCallRange (1, 2) double
    end
    
    properties (SetAccess = private)
        cActivities;
        cEntries;
        cTasks;
        cProcessors;
        activities;
        entries;
        tasks;
        processors;
        numTasksPerLevel;
        numTasksPerProcessor;
    end
    
    methods
        function obj = LayeredNetworkGenerator(varargin)
            p = inputParser();
            addParameter(p, 'populationRange', [1, 1]);
            addParameter(p, 'thinkTimeRange', [1, 1]);
            addParameter(p, 'taskInfProbability', 0);
            addParameter(p, 'procInfProbability', 0);
            addParameter(p, 'taskMultiRange', [1, 1]);
            addParameter(p, 'procMultiRange', [1, 1]);
            addParameter(p, 'hostDemandRange', [1, 1]);
            addParameter(p, 'synchCallRange', [1, 1]);
            parse(p, varargin{:});
            
            obj.populationRange = p.Results.populationRange;
            obj.thinkTimeRange = p.Results.thinkTimeRange;
            obj.taskInfProbability = p.Results.taskInfProbability;
            obj.procInfProbability = p.Results.procInfProbability;
            obj.taskMultiRange = p.Results.taskMultiRange;
            obj.procMultiRange = p.Results.procMultiRange;
            obj.hostDemandRange = p.Results.hostDemandRange;
            obj.synchCallRange = p.Results.synchCallRange;
        end
        
        function set.populationRange(obj, range)
            if range(1) > 0 && range(1) <= range(2)
                obj.populationRange(1) = range(1);
                obj.populationRange(2) = range(2);
            else
                error('NG:populationRange', 'Population range is not valid');
            end
        end
        
        function set.thinkTimeRange(obj, range)
            if range(1) >= 0 && range(1) <= range(2)
                obj.thinkTimeRange(1) = range(1);
                obj.thinkTimeRange(2) = range(2);
            else
                error('NG:thinkTimeRange', 'Think time range is not valid');
            end
        end
        
        function set.taskInfProbability(obj, probability)
            if probability >= 0 && probability <= 1
                obj.taskInfProbability = probability;
            else
                error('NG:taskInfProbability', 'Task infinite probability is not valid');
            end
        end
        
        function set.procInfProbability(obj, probability)
            if probability >= 0 && probability <= 1
                obj.procInfProbability = probability;
            else
                error('NG:procInfProbability', 'Processor infinite probability is not valid');
            end
        end
        
        function set.taskMultiRange(obj, range)
            if range(1) > 0 && range(1) <= range(2)
                obj.taskMultiRange(1) = range(1);
                obj.taskMultiRange(2) = range(2);
            else
                error('NG:taskMultiRange', 'Task multiplicity range is not valid');
            end
        end
        
        function set.procMultiRange(obj, range)
            if range(1) > 0 && range(1) <= range(2)
                obj.procMultiRange(1) = range(1);
                obj.procMultiRange(2) = range(2);
            else
                error('NG:procMultiRange', 'Processor multiplicity range is not valid');
            end
        end
        
        function set.hostDemandRange(obj, range)
            if range(1) >= 0 && range(1) <= range(2)
                obj.hostDemandRange(1) = range(1);
                obj.hostDemandRange(2) = range(2);
            else
                error('NG:hostDemandRange', 'Host demand range is not valid');
            end
        end
        
        function set.synchCallRange(obj, range)
            if range(1) > 0 && range(1) <= range(2)
                obj.synchCallRange(1) = range(1);
                obj.synchCallRange(2) = range(2);
            else
                error('NG:synchCallRange', 'Synchronous call range is not valid');
            end
        end
        
        % Main function to call. Returns a generated LQN model according to
        % specified properties of the LayeredNetworkGenerator object
        function model = generate(obj, numClients, numLevels, numTasks, numProcessors)
            obj.validateArgs(numClients, numLevels, numTasks, numProcessors);
            model = LayeredNetwork('lnw');
            obj.createClients(model, numClients);
            obj.createTasks(model, numTasks);
            obj.createProcessors(model, numProcessors);
            obj.assignTasks(numLevels, numTasks, numProcessors);
            obj.connectClientsToTasks(numClients);
            obj.connectTasksToTasks(numLevels);
            obj.connectTasksToProcessors(numProcessors);
        end
    end
    
    methods (Access = private)
        % Validates that parameter values for the network are sound
        function validateArgs(~, numClients, numLevels, numTasks, numProcessors)
            if numClients < 1
                error('NG:invalidArgs', 'Number of clients is less than one');
            elseif numLevels < 1
                error('NG:invalidArgs', 'Number of levels is less than one');
            elseif numTasks < 1
                error('NG:invalidArgs', 'Number of tasks is less than one');
            elseif numProcessors < 1
                error('NG:invalidArgs', 'Number of processors is less than one');
            elseif numLevels > numTasks
                error('NG:invalidArgs', 'Number of levels is greater than that of tasks');
            elseif numProcessors > numTasks
                error('NG:invalidArgs', 'Number of processors is greater than that of tasks');
            end
        end
        
        % Creates the clients in the layered network
        function createClients(obj, model, numClients)
            obj.cActivities = cell(numClients, 1);
            obj.cEntries = cell(numClients, 1);
            obj.cTasks = cell(numClients, 1);
            obj.cProcessors = cell(numClients, 1);
            for c = 1 : numClients
                population = LayeredNetworkGenerator.sampleIntegerValue(obj.populationRange);
                thinkTime = LayeredNetworkGenerator.sampleRealValue(obj.thinkTimeRange);
                
                obj.cActivities{c} = Activity(model, ['c_activity_', int2str(c)], thinkTime);
                obj.cEntries{c} = Entry(model, ['c_entry_', int2str(c)]);
                obj.cTasks{c} = Task(model, ['c_task_', int2str(c)], population, SchedStrategy.REF);
                obj.cProcessors{c} = Processor(model, ['c_processor_', int2str(c)], Inf, SchedStrategy.INF);
                
                obj.cActivities{c}.on(obj.cTasks{c}).boundTo(obj.cEntries{c});
                obj.cEntries{c}.on(obj.cTasks{c});
                obj.cTasks{c}.on(obj.cProcessors{c});
            end
        end
        
        % Creates the tasks in the layered network
        function createTasks(obj, model, numTasks)
            obj.activities = cell(numTasks, 1);
            obj.entries = cell(numTasks, 1);
            obj.tasks = cell(numTasks, 1);
            for t = 1 : numTasks
                isInfinite = LayeredNetworkGenerator.chooseBooleanValue(obj.taskInfProbability);
                if isInfinite
                    multiplicity = Inf;
                    scheduling = SchedStrategy.INF;
                else
                    multiplicity = LayeredNetworkGenerator.sampleIntegerValue(obj.taskMultiRange);
                    scheduling = SchedStrategy.FCFS;
                end
                hostDemand = LayeredNetworkGenerator.sampleRealValue(obj.hostDemandRange);
                
                obj.activities{t} = Activity(model, ['activity_', int2str(t)], hostDemand);
                obj.entries{t} = Entry(model, ['entry_', int2str(t)]);
                obj.tasks{t} = Task(model, ['task_', int2str(t)], multiplicity, scheduling);
                
                obj.activities{t}.on(obj.tasks{t}).boundTo(obj.entries{t}).repliesTo(obj.entries{t});
                obj.entries{t}.on(obj.tasks{t});
            end
        end
        
        % Creates the processors in the layered network
        function createProcessors(obj, model, numProcessors)
            obj.processors = cell(numProcessors, 1);
            for i = 1 : numProcessors
                isInfinite = LayeredNetworkGenerator.chooseBooleanValue(obj.procInfProbability);
                if isInfinite
                    multiplicity = Inf;
                    scheduling = SchedStrategy.INF;
                else
                    multiplicity = LayeredNetworkGenerator.sampleIntegerValue(obj.procMultiRange);
                    scheduling = SchedStrategy.PS;
                end
                
                obj.processors{i} = Processor(model, ['processor_', int2str(i)], multiplicity, scheduling);
            end
        end
        
        % Assigns the tasks to different levels and processors
        function assignTasks(obj, numLevels, numTasks, numProcessors)
            obj.numTasksPerLevel = LayeredNetworkGenerator.makeIntegerVector(numLevels, numTasks);
            obj.numTasksPerProcessor = LayeredNetworkGenerator.makeIntegerVector(numProcessors, numTasks);
        end
        
        % Connects the clients to the tasks
        function connectClientsToTasks(obj, numClients)
            clientConnected = false(1, numClients);
            for t = 1 : obj.numTasksPerLevel(1)
                synchCall = LayeredNetworkGenerator.sampleRealValue(obj.synchCallRange);
                
                c = LayeredNetworkGenerator.sampleIntegerValue([1, numClients]);
                obj.cActivities{c}.synchCall(obj.entries{t}, synchCall);
                clientConnected(c) = true;
            end
            
            for c = 1 : numClients
                if clientConnected(c)
                    continue
                end
                
                synchCall = LayeredNetworkGenerator.sampleRealValue(obj.synchCallRange);
                
                t = LayeredNetworkGenerator.sampleIntegerValue([1, obj.numTasksPerLevel(1)]);
                obj.cActivities{c}.synchCall(obj.entries{t}, synchCall);
                clientConnected(c) = true;
            end
        end
        
        % Connects the tasks between adjacent levels
        function connectTasksToTasks(obj, numLevels)
            numTasks = obj.numTasksPerLevel(1);
            for l = 2 : numLevels
                for t2 = numTasks + (1 : obj.numTasksPerLevel(l))
                    synchCall = LayeredNetworkGenerator.sampleRealValue(obj.synchCallRange);
                    
                    t1 = LayeredNetworkGenerator.sampleIntegerValue(numTasks - [obj.numTasksPerLevel(l - 1) - 1, 0]);
                    obj.activities{t1}.synchCall(obj.entries{t2}, synchCall);
                end
                numTasks = numTasks + obj.numTasksPerLevel(l);
            end
        end
        
        % Connects the tasks to the processors
        function connectTasksToProcessors(obj, numProcessors)
            numTasks = 0;
            for p = 1 : numProcessors
                for t = numTasks + (1 : obj.numTasksPerProcessor(p))
                    obj.tasks{t}.on(obj.processors{p});
                end
                numTasks = numTasks + obj.numTasksPerProcessor(p);
            end
        end
    end
    
    methods (Static)
        % Samples an integer value from a given range
        function value = sampleIntegerValue(range)
            value = randi([ceil(range(1)), floor(range(2))]);
        end
        
        % Samples a real value from a given range
        function value = sampleRealValue(range)
            value = range(1) + (range(2) - range(1)) * rand();
        end
        
        % Chooses a boolean value for a given probability
        function value = chooseBooleanValue(probability)
            value = rand() < probability;
        end
        
        % Makes an integer vector with given length and sum
        function vector = makeIntegerVector(length, sum)
            vector = ones(length, 1);
            for s = 1 : sum - length
                i = randi([1, length]);
                vector(i) = vector(i) + 1;
            end
        end
    end
    
end
