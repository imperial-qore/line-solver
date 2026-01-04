classdef JobClass < NetworkElement
    % An abstract class for a collection of indistinguishable jobs
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        priority;
        deadline; % relative deadline from arrival (Inf = no deadline)
        refstat; % reference station
        isrefclass; % is this a reference class within a chain?
        index; % node index
        type;
        completes; % true if passage through reference station is a completion
        replySignalClass; % Signal class that will unblock servers waiting for reply (for synchronous call semantics)
        patience; % Global patience distribution for this class (customers abandon queues after patience time)
        impatienceType; % Global impatience type for this class (ImpatienceType.RENEGING or ImpatienceType.BALKING)
        immediateFeedback; % Boolean: true if this class uses immediate feedback on self-loops globally
    end
    
    methods (Hidden)
        %Constructor
        function self = JobClass(type, name)
            % SELF = JOBCLASS(TYPE, NAME)

            self@NetworkElement(name);
            self.priority = 0;
            self.deadline = Inf;
            self.refstat = Node('Unallocated');
            self.isrefclass = false;
            self.index = 1;
            self.type=type;
            self.completes = true;
            self.replySignalClass = [];
            self.patience = [];
            self.impatienceType = [];
            self.immediateFeedback = false;
        end
        
        function self = setReferenceStation(self, source)
            % SELF = SETREFERENCESTATION(SOURCE)
            
            self.refstat = source;
        end
        
        function self = setReferenceClass(self, bool)
            % SELF = SETREFERENCECLASS(BOOL)            
            self.isrefclass = bool;
        end
        
        function boolIsa = isReferenceStation(self, node)
            % BOOLISA = ISREFERENCESTATION(NODE)
            
            boolIsa = strcmp(self.refstat.name,node.name);
        end
        
        function boolIsa = isReferenceClass(self)
            % BOOLISA = ISREFERENCECLASS()
            
            boolIsa = self.isrefclass;
        end
        
        %         function self = set.priority(self, priority)
        % SELF = SET.PRIORITY(PRIORITY)
        
        %             if ~(rem(priority,1) == 0 && priority >= 0)
        %                 line_error(mfilename,'Priority must be an integer.\n');
        %             end
        %             self.priority = priority;
        %         end
    end
    
    methods (Access=public)
        function ind = subsindex(self)
            % IND = SUBSINDEX()

            ind = double(self.index)-1; % 0 based
        end

        function summary(self)
            % SUMMARY()
            line_printf('Class (%s): <strong>%s</strong>',self.type,self.getName);
        end

        function self = setReplySignalClass(self, replyClass)
            % SETREPLYSIGNALCLASS Set the Signal class for synchronous call reply
            %
            % self = SETREPLYSIGNALCLASS(self, replyClass) configures this job
            % class to expect a reply signal from the specified Signal class.
            % When a job of this class completes service, the server will block
            % until receiving a REPLY signal from the specified class.
            %
            % This implements LQN-style synchronous call semantics where a
            % client sends a request, blocks waiting for a reply, and then
            % continues processing after the reply arrives.
            %
            % @param replyClass Signal object with SignalType.REPLY that will unblock the server
            % @return self The modified JobClass instance

            if ~isa(replyClass, 'Signal')
                line_error(mfilename, 'Reply class must be a Signal.');
            end
            if replyClass.signalType ~= SignalType.REPLY
                line_error(mfilename, 'Reply class must have SignalType.REPLY.');
            end
            self.replySignalClass = replyClass;
        end

        function idx = getReplySignalClassIndex(self)
            % GETREPLYSIGNALCLASSINDEX Get the index of the reply signal class
            %
            % idx = GETREPLYSIGNALCLASSINDEX(self) returns the index of the
            % Signal class that will unblock servers waiting for this class,
            % or -1 if no reply is expected.
            %
            % @return idx Index of reply signal class, or -1 if none

            if isempty(self.replySignalClass)
                idx = -1;
            else
                idx = self.replySignalClass.index;
            end
        end

        function tf = expectsReply(self)
            % EXPECTSREPLY Check if this class expects a reply signal
            %
            % tf = EXPECTSREPLY(self) returns true if this class has been
            % configured to expect a reply signal (via setReplySignalClass).
            %
            % @return tf True if a reply signal is expected

            tf = ~isempty(self.replySignalClass);
        end

        function self = setPatience(self, varargin)
            % SELF = SETPATIENCE(DISTRIBUTION) - Backwards compatible
            % SELF = SETPATIENCE(IMPATIENCETYPE, DISTRIBUTION) - Explicit type
            %
            % Sets the global impatience type and distribution for this job class.
            % This applies to all queues unless overridden by queue-specific settings.
            %
            % Parameters:
            %   impatienceType - (Optional) ImpatienceType constant (RENEGING or BALKING)
            %                  If omitted, defaults to ImpatienceType.RENEGING
            %   distribution - Any LINE distribution (Exp, Erlang, HyperExp, etc.)
            %                  excluding modulated processes (BMAP, MAP, MMPP2)
            %
            % Examples:
            %   jobclass.setPatience(Exp(0.1))  % Defaults to RENEGING
            %   jobclass.setPatience(ImpatienceType.RENEGING, Exp(0.1))
            %   jobclass.setPatience(ImpatienceType.BALKING, Det(5.0))

            % Handle backwards compatibility: 1 or 2 arguments
            if length(varargin) == 1
                % Old signature: setPatience(distribution)
                distribution = varargin{1};
                impatienceType = ImpatienceType.RENEGING;  % Default to RENEGING
            elseif length(varargin) == 2
                % New signature: setPatience(impatienceType, distribution)
                impatienceType = varargin{1};
                distribution = varargin{2};
            else
                line_error(mfilename, 'Invalid number of arguments. Use setPatience(distribution) or setPatience(impatienceType, distribution)');
            end

            if isa(distribution, 'BMAP') || isa(distribution, 'MAP') || isa(distribution, 'MMPP2')
                line_error(mfilename, 'Modulated processes (BMAP, MAP, MMPP2) are not supported for patience distributions.');
            end

            % Validate impatience type
            if impatienceType ~= ImpatienceType.RENEGING && impatienceType ~= ImpatienceType.BALKING
                line_error(mfilename, 'Invalid impatience type. Use ImpatienceType.RENEGING or ImpatienceType.BALKING.');
            end

            % Only RENEGING is currently supported
            if impatienceType == ImpatienceType.BALKING
                line_error(mfilename, 'BALKING impatience type is not yet supported. Use ImpatienceType.RENEGING.');
            end

            if isempty(self.obj)
                self.patience = distribution;
                self.impatienceType = impatienceType;
            else
                self.obj.setPatience(impatienceType, distribution.obj);
            end
        end

        function distribution = getPatience(self)
            % DISTRIBUTION = GETPATIENCE()
            %
            % Returns the global patience distribution for this job class.
            %
            % Returns:
            %   distribution - The patience distribution, or [] if not set

            if isempty(self.obj)
                if ~isempty(self.patience)
                    distribution = self.patience;
                else
                    distribution = [];
                end
            else
                distObj = self.obj.getPatience();
                if isempty(distObj)
                    distribution = [];
                else
                    % Convert Java object to MATLAB Distribution
                    distribution = Distribution.fromJavaObject(distObj);
                end
            end
        end

        function impatienceType = getImpatienceType(self)
            % IMPATIENCETYPE = GETIMPATIENCETYPE()
            %
            % Returns the global impatience type for this job class.
            %
            % Returns:
            %   impatienceType - The impatience type (ImpatienceType constant), or [] if not set

            if isempty(self.obj)
                if ~isempty(self.impatienceType)
                    impatienceType = self.impatienceType;
                else
                    impatienceType = [];
                end
            else
                impatienceTypeObj = self.obj.getImpatienceType();
                if isempty(impatienceTypeObj)
                    impatienceType = [];
                else
                    impatienceType = ImpatienceType.fromId(impatienceTypeObj.getID());
                end
            end
        end

        function tf = hasPatience(self)
            % TF = HASPATIENCE()
            %
            % Returns true if this class has a patience distribution set.

            dist = self.getPatience();
            tf = ~isempty(dist) && ~isa(dist, 'Disabled');
        end

        function self = setImmediateFeedback(self, value)
            % SETIMMEDIATEFEEDBACK Set immediate feedback for self-loops
            %
            % SETIMMEDIATEFEEDBACK(true) enables immediate feedback for this class globally
            % SETIMMEDIATEFEEDBACK(false) disables immediate feedback for this class
            %
            % When enabled, a job of this class that self-loops at any station stays
            % in service instead of going back to the queue.

            if isempty(self.obj)
                self.immediateFeedback = logical(value);
            else
                self.obj.setImmediateFeedback(logical(value));
            end
        end

        function tf = hasImmediateFeedback(self)
            % HASIMMEDIATEFEEDBACK Check if immediate feedback is enabled
            %
            % TF = HASIMMEDIATEFEEDBACK() returns true if immediate feedback is enabled

            if isempty(self.obj)
                tf = self.immediateFeedback;
            else
                tf = self.obj.hasImmediateFeedback();
            end
        end
    end

end
