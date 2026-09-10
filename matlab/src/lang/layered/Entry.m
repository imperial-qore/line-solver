classdef Entry < LayeredNetworkElement
    % An entry point of service for a Task.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        parent;
        replyActivity = {};
        arrival;  % Open arrival distribution
        scheduling = [];
        forwardingDests = {};
        forwardingProbs = [];
        type = 'PH1PH2';  % string, entry type (default: PH1PH2, can be NONE, etc.)
    end
    
    methods
        %public methods, including constructor
        
        %constructor
        function self = Entry(model, name)
            % SELF = ENTRY(MODEL, NAME)
            name = char(name);
            if nargin<2 %~exist('name','var')
                line_error(mfilename,'Constructor requires to specify at least a name.');
            end
            self@LayeredNetworkElement(name);

            if isa(model,'LayeredNetwork')
                self.arrival = [];
                model.entries{end+1} = self;
                self.model = model;
            elseif isa(model,'JLayeredNetwork')
                % JLayeredNetwork support would go here if it exists
                self.obj = jline.lang.layered.Entry(model.obj, name);
                model.addEntry(self);
            end
        end
        
        function self = on(self, parent)
            % SELF = ON(SELF, PARENT)

            parent.addEntry(self);
            self.parent = parent;
        end

        function self = setArrival(self, arvDist)
            % SETARRIVAL Sets the open arrival distribution for this entry
            %
            % Parameters:
            %   arvDist - Arrival distribution (e.g., Exp, Erlang, HyperExp)
            %
            % Examples:
            %   entry.setArrival(Exp(2.5));  % Exponential with rate 2.5
            %   entry.setArrival(Erlang(2, 0.4));  % Erlang-2 with rate 0.4

            if ~isa(arvDist, 'Distribution')
                line_error(mfilename, 'Arrival must be a Distribution object');
            end
            self.arrival = arvDist;
        end

        function self = setType(self, type)
            % SELF = SETTYPE(self, TYPE)
            % Set the entry type attribute (e.g., PH1PH2, NONE)

            self.type = type;
        end

        function self = forward(self, dest, prob)
            % SELF = FORWARD(SELF, DEST, PROB)
            %
            % Add a forwarding call to another entry with a specified probability.
            % Forwarding allows this entry to redirect the reply to another entry
            % instead of replying directly to the original caller.
            %
            % Parameters:
            %   dest - Destination entry object or entry name (string)
            %   prob - Probability of forwarding (0.0 to 1.0), default is 1.0
            %
            % Returns:
            %   self - This entry for method chaining

            if nargin < 3
                prob = 1.0;
            end

            % Validate probability
            if prob < 0.0 || prob > 1.0
                line_error(mfilename, sprintf('Forwarding probability must be between 0.0 and 1.0, got: %f', prob));
            end

            % Check sum of probabilities
            currentSum = sum(self.forwardingProbs);
            if currentSum + prob > 1.0 + 1e-6  % Small tolerance for floating point errors
                line_error(mfilename, sprintf('Sum of forwarding probabilities would exceed 1.0 (current: %.6f, adding: %.6f)', currentSum, prob));
            end

            % Get destination name
            if isa(dest, 'Entry')
                destName = dest.name;
            elseif ischar(dest) || isstring(dest)
                destName = char(dest);
            else
                line_error(mfilename, 'Destination must be an Entry object or a string');
            end

            % Add to forwarding lists
            self.forwardingDests{end+1} = destName;
            self.forwardingProbs(end+1) = prob;
        end

        % Getter methods for API consistency with Java/Python
        function val = getParent(obj)
            % GETPARENT Get the parent task
            val = obj.parent;
        end

        function val = getReplyActivity(obj)
            % GETREPLYACTIVITY Get the reply activities for this entry
            val = obj.replyActivity;
        end

        function val = getArrival(obj)
            % GETARRIVAL Get the arrival distribution
            val = obj.arrival;
        end

        function val = getForwardingDests(obj)
            % GETFORWARDINGDESTS Get the forwarding destinations
            val = obj.forwardingDests;
        end

        function val = getForwardingProbs(obj)
            % GETFORWARDINGPROBS Get the forwarding probabilities
            val = obj.forwardingProbs;
        end

    end

end
