classdef Source < Station
    % Source External job arrival node for open queueing networks
    %
    % Source represents an external arrival node that generates jobs for open
    % classes according to specified arrival processes. It serves as the entry
    % point for jobs entering the network from the external environment, with
    % configurable arrival rates and distributions for each job class.
    %
    % @brief External arrival node generating jobs for open queueing networks
    %
    % Key characteristics:
    % - External job generation for open classes
    % - Class-dependent arrival processes
    % - Infinite capacity job source
    % - Configurable inter-arrival time distributions
    % - Integration with network routing
    %
    % Source node features:
    % - Multiple job class support
    % - Flexible arrival process specification
    % - Poisson, MAP, and general arrival processes
    % - Arrival rate configuration per class
    % - Disabled arrival capability for specific classes
    %
    % Source is used for:
    % - Web server client arrivals
    % - Manufacturing job arrivals
    % - Call center customer generation
    % - Network packet injection
    % - Open system workload modeling
    %
    % Example:
    % @code
    % model = Network('WebServer');
    % source = Source(model, 'ClientArrivals');
    % webClass = OpenClass(model, 'WebRequests', 1);
    % source.setArrival(webClass, Exp(2.0)); % Poisson arrivals, rate 2
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        schedStrategy;
        arrivalProcess;
    end

    methods
        %Constructor
        function self = Source(model, name)
            % SOURCE Create an external arrival source node
            %
            % @brief Creates a Source node for external job generation
            % @param model Network model to add the source node to
            % @param name String identifier for the source node
            % @return self Source instance ready for arrival process configuration
            self@Station(name);
            if model.isMatlabNative()
                self.numberOfServers = 1;
                if(model ~= 0)
                    classes = model.getClasses();
                    self.classCap = Inf*ones(1,length(classes));
                    self.output = Dispatcher(classes);
                    self.server = ServiceTunnel();
                    self.input = RandomSource(classes);
                    self.schedStrategy = SchedStrategy.EXT;
                    self.setModel(model);
                    model.addNode(self);                    
                end
            elseif model.isJavaNative()
                self.setModel(model);
                self.obj=jline.lang.nodes.Source(model.obj, name);
            end
        end

        function setArrival(self, class, distribution)
            % SETARRIVAL(CLASS, DISTRIBUTION)
            % distribution can be a Distribution object or a Workflow object

            % If Workflow, convert to PH distribution
            if isa(distribution, 'Workflow')
                distribution = distribution.toPH();
            end

            if isempty(self.obj)
                % Check if arrival was already configured
                if length(self.input.sourceClasses) >= class.index && ~isempty(self.input.sourceClasses{1, class.index})
                    % Note: We no longer invalidate hasStruct here as it causes severe performance
                    % issues in iterative solvers like LN. The refreshRates/refreshProcesses methods
                    % called during solver post-iteration phase handle updating procid appropriately.
                    self.model.setInitialized(false);
                end
                self.input.sourceClasses{1, class.index}{2} = ServiceStrategy.LI;
                self.input.sourceClasses{1, class.index}{3} = distribution;
                self.arrivalProcess{1,class.index} = distribution;
                if distribution.isDisabled()
                    self.classCap(class.index) = 0;
                else
                    self.classCap(class.index) = Inf;
                end
                % Update cached procid if struct exists to avoid stale values
                % This is needed because we don't invalidate hasStruct for performance
                if self.model.hasStruct && ~isempty(self.model.sn)
                    ist = self.model.getStationIndex(self);
                    c = class.index;
                    procTypeId = ProcessType.toId(ProcessType.fromText(builtin('class', distribution)));
                    self.model.sn.procid(ist, c) = procTypeId;
                end
            else
                self.obj.setArrival(class.obj, distribution.obj);
                % Also update MATLAB-side storage to keep in sync with Java object
                % This ensures getArrivalProcess returns the correct distribution
                % Check if arrival was already configured
                if length(self.input.sourceClasses) >= class.index && ~isempty(self.input.sourceClasses{1, class.index})
                    % Note: We no longer invalidate hasStruct here as it causes severe performance
                    % issues in iterative solvers like LN. The refreshRates/refreshProcesses methods
                    % called during solver post-iteration phase handle updating procid appropriately.
                    self.model.setInitialized(false);
                end
                self.input.sourceClasses{1, class.index}{2} = ServiceStrategy.LI;
                self.input.sourceClasses{1, class.index}{3} = distribution;
                self.arrivalProcess{1,class.index} = distribution;
                if distribution.isDisabled()
                    self.classCap(class.index) = 0;
                else
                    self.classCap(class.index) = Inf;
                end
                % Update cached procid if struct exists to avoid stale values
                if self.model.hasStruct && ~isempty(self.model.sn)
                    ist = self.model.getStationIndex(self);
                    c = class.index;
                    procTypeId = ProcessType.toId(ProcessType.fromText(builtin('class', distribution)));
                    self.model.sn.procid(ist, c) = procTypeId;
                end
            end
        end

        function distrib = getArrivalProcess(self, oclass)
            distrib = self.arrivalProcess{oclass};
        end

    end

end
