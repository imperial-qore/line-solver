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
        markedProcess;   % shared MarkedMAP driving the marked classes (empty if none)
        markedClasses;   % (1,K) class indexes bound to marks 1..K of markedProcess
        arrivalBatch;    % (1,K) cell of batch-size DiscreteDistribution per class (empty = single arrivals)
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

        function setArrivalBatch(self, class, batchSize)
            % SETARRIVALBATCH(CLASS, BATCHSIZE)
            %
            % Turns each arrival epoch of CLASS into the simultaneous release of
            % a batch of jobs. The interarrival distribution set by SETARRIVAL
            % keeps spacing the epochs; this decides how many jobs each epoch
            % releases. Geometric interarrivals with a Geometric batch size is
            % the Geo^X arrival stream, whose analytical counterpart is
            % QSYS_GEOXGEO1.
            %
            % The batch size must be supported on {1,2,...}: an epoch that
            % releases no job is not an arrival epoch, so a law that can return
            % zero is rejected rather than clamped. Pass [] to restore single
            % arrivals.
            if isempty(batchSize)
                if numel(self.arrivalBatch) >= class.index
                    self.arrivalBatch{1, class.index} = [];
                end
                return
            end
            if ~isa(batchSize, 'DiscreteDistribution')
                line_error(mfilename, 'arrival batch size must be a DiscreteDistribution');
            end
            if batchSize.getMean() < 1
                line_error(mfilename, sprintf(['arrival batch size for class ''%s'' has mean %g; ' ...
                    'a batch must carry at least one job, so its support must be {1,2,...}'], ...
                    class.name, batchSize.getMean()));
            end
            if isempty(self.obj)
                self.arrivalBatch{1, class.index} = batchSize;
            else
                self.obj.setArrivalBatch(class.obj, batchSize.obj);
            end
        end

        function self = removeJobClass(self, jobclass)
            % SELF = REMOVEJOBCLASS(JOBCLASS)
            %
            % Drop the arrival process of JOBCLASS, at the source and inside
            % its input section, on top of what Station and Node remove.

            remaining = self.remainingClassIndexes(jobclass);
            removeJobClass@Station(self, jobclass);
            K = numel(remaining) + 1;
            if numel(self.arrivalProcess) == K
                self.arrivalProcess = self.arrivalProcess(remaining);
            end
            if numel(self.arrivalBatch) == K
                self.arrivalBatch = self.arrivalBatch(remaining);
            end
            if numel(self.markedClasses) == K
                self.markedClasses = self.markedClasses(remaining);
            end
            if ~isempty(self.input) && isprop(self.input, 'sourceClasses') ...
                    && numel(self.input.sourceClasses) == K
                self.input.sourceClasses = self.input.sourceClasses(remaining);
            end
        end

        function batchSize = getArrivalBatch(self, class)
            % BATCHSIZE = GETARRIVALBATCH(CLASS)
            % Returns the batch-size law bound to CLASS, or [] for single arrivals.
            batchSize = [];
            idx = class;
            if isobject(class) && isprop(class, 'index')
                idx = class.index;
            end
            if numel(self.arrivalBatch) >= idx
                batchSize = self.arrivalBatch{1, idx};
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

        function setMarkedArrival(self, mmap, classes)
            % SETMARKEDARRIVAL(MMAP, CLASSES)
            %
            % Bind a MarkedMAP with K marks to K open classes: mark k emits
            % jobs of class CLASSES{k}, with all marks driven by one shared
            % modulating chain. CLASSES is a cell array or vector of K
            % distinct OpenClass handles, ordered by mark index.
            if ~isa(mmap, 'MarkedMAP')
                line_error(mfilename, 'setMarkedArrival requires a MarkedMAP/MMAP arrival process.');
            end
            if ~isempty(self.obj)
                line_error(mfilename, 'setMarkedArrival is not yet supported with the Java backend (lang=''java'').');
            end
            K = mmap.getNumberOfTypes;
            if iscell(classes)
                classList = classes;
            else
                classList = num2cell(classes);
            end
            if numel(classList) ~= K
                line_error(mfilename, sprintf('The MarkedMAP has %d types but %d classes were supplied.', K, numel(classList)));
            end
            markClasses = zeros(1, K);
            for k = 1:K
                cls = classList{k};
                if ~isa(cls, 'OpenClass')
                    line_error(mfilename, 'setMarkedArrival requires open classes.');
                end
                markClasses(k) = cls.index;
            end
            if numel(unique(markClasses)) ~= K
                line_error(mfilename, 'setMarkedArrival requires distinct classes for the marks.');
            end
            % Per-class binding: every marked class stores the shared MarkedMAP
            % object, so procid resolves to MMAP and rates/SCV are derived from
            % the per-mark marginal in the refresh layer.
            for k = 1:K
                self.setArrival(classList{k}, mmap);
            end
            self.markedProcess = mmap;
            self.markedClasses = markClasses;
        end

    end

end
