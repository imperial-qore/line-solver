classdef OpenSignal < OpenClass
    % OpenSignal Signal class for open queueing networks
    %
    % OpenSignal is a specialized OpenClass for modeling signals in open
    % queueing networks. Unlike regular customers, signals can have special
    % effects on queues they visit, such as removing jobs (negative signals)
    % or unblocking servers (reply signals).
    %
    % For closed networks, use ClosedSignal instead.
    %
    % Signal types:
    % - SignalType.NEGATIVE: Removes a job from the destination queue
    % - SignalType.REPLY: Unblocks servers waiting for a reply
    %
    % Example:
    % @code
    % model = Network('OpenModel');
    % source = Source(model, 'Source');
    % sink = Sink(model, 'Sink');
    % queue = Queue(model, 'Queue', SchedStrategy.FCFS);
    % reqClass = OpenClass(model, 'Request');
    % replySignal = OpenSignal(model, 'Reply', SignalType.REPLY).forJobClass(reqClass);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        signalType           % SignalType constant (NEGATIVE, REPLY)
        targetJobClass       % JobClass that this signal is associated with
        removalDistribution  % DiscreteDistribution for number of removals (empty = remove exactly 1)
        removalPolicy        % RemovalPolicy constant (RANDOM, FCFS, LCFS)
    end

    methods

        function self = OpenSignal(model, name, signalType, prio)
            % OPENSIGNAL Create an open signal class instance
            %
            % @param model Network model to add the signal class to
            % @param name String identifier for the signal class
            % @param signalType SignalType constant (default: SignalType.NEGATIVE)
            % @param prio Optional priority level (default: 0)
            % @return self OpenSignal instance

            if nargin < 4
                prio = 0;
            end
            if nargin < 3 || isempty(signalType)
                signalType = SignalType.NEGATIVE;
            end

            self@OpenClass(model, name, prio);
            self.signalType = signalType;
            self.targetJobClass = [];
            self.removalDistribution = [];
            self.removalPolicy = RemovalPolicy.RANDOM;

            if model.isJavaNative()
                % Replace the OpenClass Java object with a Signal object
                self.obj = jline.lang.Signal(model.obj, name, signalType, prio);
            end
        end

        function type = getSignalType(self)
            % GETSIGNALTYPE Get the signal type
            type = self.signalType;
        end

        function self = forJobClass(self, jobClass)
            % FORJOBCLASS Associate this signal with a job class
            %
            % For REPLY signals, this specifies which job class's servers
            % will be unblocked when this signal arrives.
            %
            % @param jobClass The JobClass to associate with this signal
            % @return self The modified Signal instance (for chaining)

            self.targetJobClass = jobClass;
            if ~isempty(jobClass)
                jobClass.replySignalClass = self;
            end
        end

        function jobClass = getTargetJobClass(self)
            % GETTARGETJOBCLASS Get the associated job class
            jobClass = self.targetJobClass;
        end

        function idx = getTargetJobClassIndex(self)
            % GETTARGETJOBCLASSINDEX Get the index of the associated job class
            if isempty(self.targetJobClass)
                idx = -1;
            else
                idx = self.targetJobClass.index;
            end
        end

        function b = isCatastrophe(self)
            % ISCATASTROPHE Check if this is a catastrophe signal
            %
            % @return b true if signalType is SignalType.CATASTROPHE
            b = (self.signalType == SignalType.CATASTROPHE);
        end

        function dist = getRemovalDistribution(self)
            % GETREMOVALDISTRIBUTION Get the removal distribution
            dist = self.removalDistribution;
        end

        function self = setRemovalDistribution(self, dist)
            % SETREMOVALDISTRIBUTION Set the removal distribution
            self.removalDistribution = dist;
        end

        function policy = getRemovalPolicy(self)
            % GETREMOVALPOLICY Get the removal policy
            policy = self.removalPolicy;
        end

        function self = setRemovalPolicy(self, policy)
            % SETREMOVALPOLICY Set the removal policy
            self.removalPolicy = policy;
        end

    end

end
