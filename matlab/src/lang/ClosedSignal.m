classdef ClosedSignal < ClosedClass
    % ClosedSignal Signal class for closed queueing networks
    %
    % ClosedSignal is a specialized ClosedClass for modeling signals in closed
    % queueing networks. Unlike regular customers, signals can have special
    % effects on queues they visit, such as removing jobs (negative signals)
    % or unblocking servers (reply signals).
    %
    % For open networks, use OpenSignal instead.
    %
    % ClosedSignal has zero population - signals are created dynamically
    % through class switching from the target job class.
    %
    % Signal types:
    % - SignalType.NEGATIVE: Removes a job from the destination queue
    % - SignalType.REPLY: Unblocks servers waiting for a reply
    %
    % Example:
    % @code
    % model = Network('ClosedModel');
    % delay = Delay(model, 'Think');
    % queue1 = Queue(model, 'Client', SchedStrategy.FCFS);
    % queue2 = Queue(model, 'Server', SchedStrategy.FCFS);
    % jobClass = ClosedClass(model, 'Job', 5, delay);
    % replySignal = ClosedSignal(model, 'Reply', SignalType.REPLY, delay).forJobClass(jobClass);
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

        function self = ClosedSignal(model, name, signalType, refstat, prio, removalDistribution, removalPolicy)
            % CLOSEDSIGNAL Create a closed signal class instance
            %
            % @param model Network model to add the signal class to
            % @param name String identifier for the signal class
            % @param signalType SignalType constant (default: SignalType.NEGATIVE)
            % @param refstat Reference station (should match target job class)
            % @param prio Optional priority level (default: 0)
            % @param removalDistribution Optional discrete distribution for batch removals (default: [])
            % @param removalPolicy Optional RemovalPolicy constant (default: RemovalPolicy.RANDOM)
            % @return self ClosedSignal instance

            if nargin < 7 || isempty(removalPolicy)
                removalPolicy = RemovalPolicy.RANDOM;
            end
            if nargin < 6
                removalDistribution = [];
            end
            if nargin < 5 || isempty(prio)
                prio = 0;
            end
            if nargin < 4
                line_error(mfilename, 'ClosedSignal requires a reference station.');
            end
            if nargin < 3 || isempty(signalType)
                signalType = SignalType.NEGATIVE;
            end

            % ClosedSignal has 0 population - signals are created by class switching
            self@ClosedClass(model, name, 0, refstat, prio);
            self.signalType = signalType;
            self.targetJobClass = [];
            self.removalDistribution = removalDistribution;
            self.removalPolicy = removalPolicy;
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

        function b = isCatastrophe(self)
            % ISCATASTROPHE Check if this is a catastrophe signal
            %
            % @return b true if signalType is SignalType.CATASTROPHE
            b = (self.signalType == SignalType.CATASTROPHE);
        end

        function summary(self)
            % SUMMARY()
            line_printf('Signal (%s): <strong>%s</strong> [%s]', self.type, self.getName, self.signalType);
        end

    end

end
