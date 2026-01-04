classdef Signal < JobClass
    % Signal Job class representing a signal (e.g., negative customer in G-networks)
    %
    % Signal is a placeholder class that automatically resolves to OpenSignal
    % or ClosedSignal based on the network structure. Users can simply use
    % Signal in both open and closed networks - the resolution happens when
    % the model is finalized (during getStruct/refreshStruct).
    %
    % @brief Job class for modeling signals in G-networks and related models
    %
    % Key characteristics:
    % - Automatically resolves to OpenSignal or ClosedSignal
    % - Supports different signal types (NEGATIVE, REPLY, CATASTROPHE)
    % - NEGATIVE signals remove jobs from destination queues
    % - CATASTROPHE signals reset the state of queues
    % - Used in G-networks (Gelenbe networks)
    %
    % Signal types:
    % - SignalType.NEGATIVE: Removes a job from the destination queue
    % - SignalType.REPLY: Triggers a reply action
    % - SignalType.CATASTROPHE: Resets destination queue to empty state
    %
    % Example (Open Network):
    % @code
    % model = Network('GNetwork');
    % source = Source(model, 'Source');
    % sink = Sink(model, 'Sink');
    % queue = Queue(model, 'Queue', SchedStrategy.FCFS);
    % posClass = OpenClass(model, 'Positive');       % Normal customers
    % negClass = Signal(model, 'Negative', SignalType.NEGATIVE);  % Resolves to OpenSignal
    % source.setArrival(posClass, Exp(1.0));
    % source.setArrival(negClass, Exp(0.3));
    % @endcode
    %
    % Example (Closed Network):
    % @code
    % model = Network('ClosedGNetwork');
    % delay = Delay(model, 'Think');
    % queue = Queue(model, 'Queue', SchedStrategy.FCFS);
    % jobClass = ClosedClass(model, 'Job', 5, delay);
    % replySignal = Signal(model, 'Reply', SignalType.REPLY).forJobClass(jobClass);  % Resolves to ClosedSignal
    % @endcode
    %
    % Reference: Gelenbe, E. (1991). "Product-form queueing networks with
    %            negative and positive customers", Journal of Applied Probability
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        signalType           % SignalType constant (NEGATIVE, REPLY, CATASTROPHE)
        targetJobClass       % JobClass that this signal is associated with (for REPLY: the class to unblock)
        removalDistribution  % DiscreteDistribution for number of removals (empty = remove exactly 1)
        removalPolicy        % RemovalPolicy constant (RANDOM, FCFS, LCFS)
        model                % Reference to the Network model
    end

    methods

        %Constructor
        function self = Signal(model, name, signalType, prio, removalDistribution, removalPolicy)
            % SIGNAL Create a signal class instance
            %
            % @brief Creates a Signal class for G-network modeling
            % @param model Network model to add the signal class to
            % @param name String identifier for the signal class
            % @param signalType SignalType constant (REQUIRED: NEGATIVE, REPLY, or CATASTROPHE)
            % @param prio Optional priority level (default: 0)
            % @param removalDistribution Optional discrete distribution for batch removals (default: [])
            % @param removalPolicy Optional RemovalPolicy constant (default: RemovalPolicy.RANDOM)
            % @return self Signal instance ready for arrival specification

            if nargin < 3 || isempty(signalType)
                line_error(mfilename, 'signalType is required. Use SignalType.NEGATIVE, SignalType.REPLY, or SignalType.CATASTROPHE.');
            end
            if nargin < 6 || isempty(removalPolicy)
                removalPolicy = RemovalPolicy.RANDOM;
            end
            if nargin < 5
                removalDistribution = [];
            end
            if nargin < 4 || isempty(prio)
                prio = 0;
            end

            % Initialize as JobClass (type will be determined during resolution)
            self@JobClass(JobClassType.OPEN, name);  % Default to OPEN, resolved later
            self.priority = prio;
            self.signalType = signalType;
            self.targetJobClass = [];
            self.removalDistribution = removalDistribution;
            self.removalPolicy = removalPolicy;
            self.model = model;

            % Register with the model
            model.addJobClass(self);

            % Set default routing for this class at all nodes
            for i = 1:length(model.nodes)
                if isa(model.nodes{i}, 'Join')
                    model.nodes{i}.setStrategy(self, JoinStrategy.STD);
                    model.nodes{i}.setRequired(self, -1);
                end
                if ~isempty(model.nodes{i})
                    model.nodes{i}.setRouting(self, RoutingStrategy.RAND);
                end
            end

            % Java interop is handled during resolution
        end

        function concrete = resolve(self, isOpen, refstat)
            % RESOLVE Resolve this Signal placeholder to OpenSignal or ClosedSignal
            %
            % @param isOpen true if the network is open (has Source node)
            % @param refstat Reference station for closed networks (ignored for open)
            % @return concrete OpenSignal or ClosedSignal instance

            % Save properties before removing placeholder
            savedIndex = self.index;
            savedName = self.name;
            savedTargetJobClass = self.targetJobClass;

            % Save outputStrategy entries for this signal at all nodes
            % These need to be restored at the new index after resolution
            savedOutputStrategies = cell(1, length(self.model.nodes));
            for n = 1:length(self.model.nodes)
                node = self.model.nodes{n};
                if ~isempty(node.output) && isprop(node.output, 'outputStrategy')
                    outStrat = node.output.outputStrategy;
                    % Check both 1D and 2D indexing patterns
                    if iscell(outStrat) && length(outStrat) >= savedIndex
                        savedOutputStrategies{n} = outStrat{savedIndex};
                    elseif iscell(outStrat) && size(outStrat, 2) >= savedIndex
                        savedOutputStrategies{n} = outStrat{1, savedIndex};
                    end
                end
            end

            % Remove placeholder from model's classes list to avoid duplicate
            % when concrete signal is created (its constructor calls addJobClass)
            placeholderIdx = find(cellfun(@(c) c == self, self.model.classes));
            if ~isempty(placeholderIdx)
                self.model.classes(placeholderIdx) = [];
                % Update indices of classes that come after
                for c = placeholderIdx:length(self.model.classes)
                    self.model.classes{c}.index = c;
                end
            end

            % Create concrete signal (constructor will add it to model.classes)
            if isOpen
                concrete = OpenSignal(self.model, savedName, self.signalType, self.priority);
            else
                concrete = ClosedSignal(self.model, savedName, self.signalType, refstat, ...
                    self.priority, self.removalDistribution, self.removalPolicy);
            end

            % Restore saved outputStrategy entries at the new index
            newIndex = concrete.index;
            for n = 1:length(self.model.nodes)
                if ~isempty(savedOutputStrategies{n})
                    node = self.model.nodes{n};
                    if ~isempty(node.output) && isprop(node.output, 'outputStrategy')
                        % Ensure outputStrategy has enough elements
                        outStrat = node.output.outputStrategy;
                        while length(outStrat) < newIndex
                            outStrat{end+1} = {}; %#ok<AGROW>
                        end
                        outStrat{newIndex} = savedOutputStrategies{n};
                        node.output.outputStrategy = outStrat;
                    end
                end
            end

            % Copy over targetJobClass association
            if ~isempty(savedTargetJobClass)
                concrete.forJobClass(savedTargetJobClass);
            end

            % Copy removal configuration (both OpenSignal and ClosedSignal have these now)
            concrete.removalDistribution = self.removalDistribution;
            concrete.removalPolicy = self.removalPolicy;
        end

        function type = getSignalType(self)
            % GETSIGNALTYPE Get the signal type
            %
            % @return type The SignalType of this signal class
            type = self.signalType;
        end

        function self = forJobClass(self, jobClass)
            % FORJOBCLASS Associate this signal with a job class
            %
            % self = FORJOBCLASS(self, jobClass) associates this signal with
            % the specified job class. For REPLY signals, this specifies which
            % job class's servers will be unblocked when this signal arrives.
            %
            % @param jobClass The JobClass to associate with this signal
            % @return self The modified Signal instance (for chaining)
            %
            % Example:
            %   replySignal = Signal(model, 'Reply', SignalType.REPLY).forJobClass(reqClass);

            self.targetJobClass = jobClass;

            % Also update the JobClass to know about this signal (bidirectional link)
            if ~isempty(jobClass)
                jobClass.replySignalClass = self;
            end
        end

        function jobClass = getTargetJobClass(self)
            % GETTARGETJOBCLASS Get the associated job class
            %
            % @return jobClass The JobClass associated with this signal
            jobClass = self.targetJobClass;
        end

        function idx = getTargetJobClassIndex(self)
            % GETTARGETJOBCLASSINDEX Get the index of the associated job class
            %
            % @return idx Index of the associated JobClass, or -1 if none
            if isempty(self.targetJobClass)
                idx = -1;
            else
                idx = self.targetJobClass.index;
            end
        end

        function dist = getRemovalDistribution(self)
            % GETREMOVALDISTRIBUTION Get the removal distribution
            %
            % @return dist The discrete distribution for batch removals
            dist = self.removalDistribution;
        end

        function self = setRemovalDistribution(self, dist)
            % SETREMOVALDISTRIBUTION Set the removal distribution
            %
            % @param dist DiscreteDistribution for number of removals
            self.removalDistribution = dist;
        end

        function policy = getRemovalPolicy(self)
            % GETREMOVALPOLICY Get the removal policy
            %
            % @return policy The RemovalPolicy constant
            policy = self.removalPolicy;
        end

        function self = setRemovalPolicy(self, policy)
            % SETREMOVALPOLICY Set the removal policy
            %
            % @param policy RemovalPolicy constant (RANDOM, FCFS, LCFS)
            self.removalPolicy = policy;
        end

        function b = isCatastrophe(self)
            % ISCATASTROPHE Check if this is a catastrophe signal
            %
            % @return b true if signalType is SignalType.CATASTROPHE
            b = (self.signalType == SignalType.CATASTROPHE);
        end

    end

end
