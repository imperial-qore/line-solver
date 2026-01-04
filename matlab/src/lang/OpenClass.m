classdef OpenClass < JobClass
    % OpenClass Job class for external arrivals with infinite population
    %
    % OpenClass represents a job class where jobs arrive from an external source
    % with potentially infinite population. Jobs enter the network through a
    % Source node, traverse the network according to routing probabilities,
    % and exit through a Sink node. Open classes are essential for modeling
    % systems with external arrival streams.
    %
    % @brief Job class for modeling external arrivals with infinite population
    %
    % Key characteristics:
    % - Infinite external population
    % - Jobs arrive from Source nodes
    % - Jobs exit through Sink nodes
    % - Variable network population over time
    % - Arrival rate determines load intensity
    % - Priority-based service differentiation
    %
    % Open class features:
    % - External arrival process modeling
    % - Unlimited population size
    % - Dynamic network population
    % - Priority assignment for service
    % - Integration with routing strategies
    % - Performance metrics per class
    %
    % OpenClass is used for:
    % - Web server request modeling
    % - Call center customer arrivals
    % - Manufacturing job arrivals
    % - Network packet flows
    % - Service request streams
    %
    % Example:
    % @code
    % model = Network('OpenSystem');
    % source = Source(model, 'Arrivals');
    % sink = Sink(model, 'Departures');
    % job_class = OpenClass(model, 'WebRequests', 1); % Priority 1
    % source.setArrival(job_class, Exp(2.0)); % Poisson arrivals, rate 2
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods

        %Constructor
        function self = OpenClass(model, name, prio, deadline)
            % OPENCLASS Create an open job class instance
            %
            % @brief Creates an OpenClass for external arrival modeling
            % @param model Network model to add the open class to
            % @param name String identifier for the job class
            % @param prio Optional priority level (default: 0, higher = more priority)
            % @param deadline Optional relative deadline from arrival (default: Inf, no deadline)
            % @return self OpenClass instance ready for arrival specification

            self@JobClass(JobClassType.OPEN, name);
            if nargin<3
                prio = 0;
            end
            if nargin<4
                deadline = Inf;
            end
            if model.isMatlabNative()
                if isempty(model.getSource)
                    line_error(mfilename,'The model requires a Source prior to instantiating open classes.');
                end
                if isempty(model.getSink)
                    line_error(mfilename,'The model requires a Sink prior to instantiating open classes.');
                end
                if nargin >= 3 && isa(prio,'Source')
                    % user error, source passed as ref station
                    prio = 0;
                end
                self.priority = 0;
                self.deadline = Inf;
                if nargin>=3 %exist('prio','var')
                    self.priority = prio;
                end
                if nargin>=4
                    self.deadline = deadline;
                end
                model.addJobClass(self);
                setReferenceStation(self, model.getSource());

                % set default scheduling for this class at all nodes
                for i=1:length(model.nodes)
                    if isa(model.nodes{i},'Join')
                        model.nodes{i}.setStrategy(self,JoinStrategy.STD);
                        model.nodes{i}.setRequired(self,-1);
                    end
                    if ~isempty(model.nodes{i})
                        model.nodes{i}.setRouting(self, RoutingStrategy.RAND);
                    end
                end
            elseif model.isJavaNative()
                self.index = model.getNumberOfClasses + 1;
                self.obj = jline.lang.OpenClass(model.obj, name, prio, deadline);
            end
        end

        function setReferenceStation(class, source)
            % SETREFERENCESTATION(CLASS, SOURCE)

            if ~isa(source,'Source')
                line_error(mfilename,'The reference station for an open class must be a Source.');
            end
            setReferenceStation@JobClass(class, source);
        end
    end

end

