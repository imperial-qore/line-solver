classdef Queue < ServiceStation
    % A service station with queueing
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        setupTime;
        delayoffTime;
        switchoverTime;
        pollingType;
        pollingPar;
        impatienceTypes;
        % Balking properties
        balkingStrategies;        % Cell array: per-class BalkingStrategy constant
        balkingThresholds;        % Cell array: per-class balking thresholds (list of {minJobs, maxJobs, probability})
        % Retrial properties
        retrialDelays;            % Cell array: per-class retrial delay distributions
        retrialMaxAttempts;       % Array: per-class max retrial attempts (-1 = unlimited)
        orbitImpatienceDistributions;  % Cell array: per-class orbit abandonment distributions
        batchRejectProb;               % Array: per-class batch rejection probability [0,1]
        % Heterogeneous server properties
        serverTypes;              % Cell array of ServerType objects
        heteroSchedPolicy;        % HeteroSchedPolicy for server assignment
        heteroServiceDistributions;  % containers.Map: ServerType -> (containers.Map: JobClass -> Distribution)
        % Immediate feedback property
        immediateFeedback;        % Cell array of class indices, or 'all' for all classes
    end

    methods
        %Constructor
        function self = Queue(model, name, schedStrategy)
            % SELF = QUEUE(MODEL, NAME, SCHEDSTRATEGY)

            self@ServiceStation(name);

            if model.isMatlabNative()
                classes = model.getClasses();
                self.input = Buffer(classes);
                self.output = Dispatcher(classes);
                self.schedPolicy = SchedStrategyType.PR;
                self.schedStrategy = SchedStrategy.PS;
                self.serviceProcess = {};
                self.server = Server(classes);
                self.numberOfServers = 1;
                self.schedStrategyPar = zeros(1,length(model.getClasses()));
                self.setModel(model);
                self.model.addNode(self);
                self.dropRule = [];
                self.obj = [];
                self.setupTime = {};
                self.delayoffTime = {};
                self.pollingType = {};
                self.switchoverTime  = {};
                self.patienceDistributions = {};
                self.impatienceTypes = {};
                self.balkingStrategies = {};
                self.balkingThresholds = {};
                self.retrialDelays = {};
                self.retrialMaxAttempts = [];
                self.orbitImpatienceDistributions = {};
                self.batchRejectProb = [];
                self.serverTypes = {};
                self.heteroSchedPolicy = HeteroSchedPolicy.ORDER;
                self.heteroServiceDistributions = containers.Map();
                self.immediateFeedback = {};

                if nargin>=3 %exist('schedStrategy','var')
                    self.schedStrategy = schedStrategy;
                    switch SchedStrategy.toId(self.schedStrategy)
                        case {SchedStrategy.PS, SchedStrategy.DPS,SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO,SchedStrategy.GPSPRIO, SchedStrategy.LPS}
                            self.schedPolicy = SchedStrategyType.PR;
                            self.server = SharedServer(classes);
                        case {SchedStrategy.LCFSPR, SchedStrategy.LCFSPRPRIO, SchedStrategy.FCFSPR, SchedStrategy.FCFSPRPRIO, SchedStrategy.LCFSPI, SchedStrategy.LCFSPIPRIO, SchedStrategy.FCFSPI, SchedStrategy.FCFSPIPRIO, SchedStrategy.EDF}
                            self.schedPolicy = SchedStrategyType.PR;
                            self.server = PreemptiveServer(classes);
                        case {SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, SchedStrategy.SEPT, SchedStrategy.LEPT, SchedStrategy.SJF, SchedStrategy.LJF, SchedStrategy.EDD, SchedStrategy.SRPT, SchedStrategy.SRPTPRIO, SchedStrategy.PSJF, SchedStrategy.FB, SchedStrategy.LRPT}
                            self.schedPolicy = SchedStrategyType.NP;
                            self.server = Server(classes);
                        case SchedStrategy.INF
                            self.schedPolicy = SchedStrategyType.NP;
                            self.server = InfiniteServer(classes);
                            self.numberOfServers = Inf;
                        case {SchedStrategy.HOL, SchedStrategy.FCFSPRIO, SchedStrategy.LCFSPRIO}
                            self.schedPolicy = SchedStrategyType.NP;
                            self.server = Server(classes);
                        case SchedStrategy.POLLING
                            self.schedPolicy = SchedStrategyType.NP;
                            self.server = PollingServer(classes);
                        otherwise
                            line_error(mfilename,sprintf('The specified scheduling strategy (%s) is unsupported.',schedStrategy));
                    end
                end
            elseif model.isJavaNative()
                self.setModel(model);
                switch SchedStrategy.toId(schedStrategy)
                    case SchedStrategy.INF
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.INF);
                    case SchedStrategy.FCFS
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.FCFS);
                    case SchedStrategy.LCFS
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.LCFS);
                    case SchedStrategy.SIRO
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.SIRO);
                    case SchedStrategy.SJF
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.SJF);
                    case SchedStrategy.LJF
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.LJF);
                    case SchedStrategy.PS
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.PS);
                    case SchedStrategy.PSPRIO
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.PSPRIO);
                    case SchedStrategy.DPS
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.DPS);
                    case SchedStrategy.DPSPRIO
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.DPSPRIO);
                    case SchedStrategy.GPS
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.GPS);
                    case SchedStrategy.GPSPRIO
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.GPSPRIO);
                    case SchedStrategy.SEPT
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.SEPT);
                    case SchedStrategy.LEPT
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.LEPT);
                    case SchedStrategy.SRPT
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.SRPT);
                    case SchedStrategy.SRPTPRIO
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.SRPTPRIO);
                    case SchedStrategy.PSJF
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.PSJF);
                    case SchedStrategy.FB
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.FB);
                    case SchedStrategy.LRPT
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.LRPT);
                    case SchedStrategy.HOL
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.HOL);
                    case SchedStrategy.FORK
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.FORK);
                    case SchedStrategy.EXT
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.EXT);
                    case SchedStrategy.REF
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.REF);
                    case SchedStrategy.LCFSPR
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.LCFSPR);
                    case SchedStrategy.FCFSPR
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.FCFSPR);
                    case SchedStrategy.FCFSPRPRIO
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.FCFSPRPRIO);
                    case SchedStrategy.EDD
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.EDD);
                    case SchedStrategy.EDF
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.EDF);
                    case SchedStrategy.LPS
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.LPS);
                end
                self.obj.setNumberOfServers(1);
                self.index = model.obj.getNodeIndex(self.obj);
            end
        end

        function setLoadDependence(self, alpha)
            switch SchedStrategy.toId(self.schedStrategy)
                case {SchedStrategy.PS, SchedStrategy.FCFS}
                    setLimitedLoadDependence(self, alpha);
                otherwise
                    line_error(mfilename,'Load-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) stations.');
            end
        end

        function setClassDependence(self, beta)
            switch SchedStrategy.toId(self.schedStrategy)
                case {SchedStrategy.PS, SchedStrategy.FCFS}
                    setLimitedClassDependence(self, beta);
                otherwise
                    line_error(mfilename,'Class-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) stations.');
            end
        end

        function setJointDependence(self, scalingTable, cutoffs)
            % SETJOINTDEPENDENCE(SCALINGTABLE, CUTOFFS)
            %
            % Sets joint class-dependent scaling for service rates using a
            % lookup table indexed by per-class population vector.
            %
            % scalingTable: R-dimensional array or linearized 1D vector
            %   If R-dimensional: dimensions are [N1+1, N2+1, ..., NR+1]
            %   If linearized: index = 1 + n1 + n2*(N1+1) + n3*(N1+1)*(N2+1) + ...
            % cutoffs: (optional) per-class cutoffs [N1, N2, ..., NR]
            %   If omitted, uses default: ceil(6000^(1/(M*K)))

            switch SchedStrategy.toId(self.schedStrategy)
                case {SchedStrategy.PS, SchedStrategy.FCFS}
                    if nargin < 3
                        % Auto-compute cutoffs
                        M = length(self.model.stations);
                        K = length(self.model.classes);
                        defaultCutoff = ceil(6000^(1/(M*K)));
                        cutoffs = defaultCutoff * ones(1, K);
                    end

                    % Convert multi-dimensional to linearized if needed
                    if ~isvector(scalingTable)
                        scalingTable = scalingTable(:)'; % linearize
                    end

                    % Validate table size matches cutoffs
                    expectedSize = prod(cutoffs + 1);
                    if length(scalingTable) ~= expectedSize
                        line_error(mfilename, sprintf('Scaling table size (%d) does not match expected size from cutoffs (%d).', length(scalingTable), expectedSize));
                    end

                    setLimitedJointDependence(self, scalingTable, cutoffs);
                otherwise
                    line_error(mfilename,'Joint-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) stations.');
            end
        end

        function setJointClassDependence(self, scalingTables, cutoffs)
            % SETJOINTCLASSDEPENDENCE(SCALINGTABLES, CUTOFFS)
            %
            % Sets per-class joint class-dependent scaling for service rates.
            % Unlike setJointDependence which uses a single scaling factor for
            % all classes, this method allows each class to have its own
            % scaling table indexed by per-class population vector.
            %
            % This is essential for Flow-Equivalent Server (FES) aggregation
            % where the service rate for class c in state (n1,...,nK) equals
            % the throughput of class c in an isolated subnetwork.
            %
            % scalingTables: cell array {1,K} of linearized scaling vectors
            %   Each scalingTables{c} has length prod(cutoffs + 1)
            %   Index: idx = 1 + n1 + n2*(N1+1) + n3*(N1+1)*(N2+1) + ...
            % cutoffs: per-class cutoffs [N1, N2, ..., NR]

            switch SchedStrategy.toId(self.schedStrategy)
                case {SchedStrategy.PS, SchedStrategy.FCFS}
                    K = length(cutoffs);
                    expectedSize = prod(cutoffs + 1);

                    if ~iscell(scalingTables) || length(scalingTables) ~= K
                        line_error(mfilename, 'scalingTables must be cell array of length K (number of classes).');
                    end

                    for c = 1:K
                        % Convert multi-dimensional to linearized if needed
                        if ~isvector(scalingTables{c})
                            scalingTables{c} = scalingTables{c}(:)';
                        end

                        if length(scalingTables{c}) ~= expectedSize
                            line_error(mfilename, sprintf('scalingTables{%d} size (%d) does not match expected size from cutoffs (%d).', c, length(scalingTables{c}), expectedSize));
                        end
                    end

                    setLimitedJointClassDependence(self, scalingTables, cutoffs);
                otherwise
                    line_error(mfilename,'Joint-class-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) stations.');
            end
        end

        function setNumberOfServers(self, value)
            % SETNUMBEROFSERVERS(VALUE)
            if isempty(self.obj)
                switch SchedStrategy.toId(self.schedStrategy)
                    case SchedStrategy.INF
                        %line_warning(mfilename,'A request to change the number of servers in an infinite server node has been ignored.');
                        %ignore
                    otherwise
                        self.setNumServers(value);
                end
            else
                self.obj.setNumberOfServers(value);
            end
        end

        function setNumServers(self, value)
            % SETNUMSERVERS(VALUE)
            if isempty(self.obj)
                switch SchedStrategy.toId(self.schedStrategy)
                    case {SchedStrategy.DPS, SchedStrategy.GPS}
                        if value ~= 1
                            line_error(mfilename,sprintf('Cannot use multi-server stations with %s scheduling.', self.schedStrategy));
                        end
                    otherwise
                        self.numberOfServers = value;
                end
            else
                self.obj.setNumberOfServers(value);
            end
        end

        function self = setStrategyParam(self, class, weight)
            % SELF = SETSTRATEGYPARAM(CLASS, WEIGHT)

            % For LPS scheduling, schedStrategyPar(1) stores the limit set via setLimit()
            % Don't overwrite it with the default weight
            if SchedStrategy.toId(self.schedStrategy) == SchedStrategy.LPS
                % For LPS, only set weight if explicitly provided (not default 1.0)
                % or if this is not the first class (which would overwrite the limit)
                if class.index == 1 && weight == 1.0 && ~isempty(self.schedStrategyPar) && self.schedStrategyPar(1) > 1
                    % Preserve the LPS limit, don't overwrite with default weight
                    return;
                end
            end
            self.schedStrategyPar(class.index) = weight;
        end

        function distribution = getService(self, class)
            % DISTRIBUTION = GETSERVICE(CLASS)

            % return the service distribution assigned to the given class
            if nargin<2 %~exist('class','var')
                for s = 1:length(self.model.getClasses())
                    classes = self.model.getClasses();
                distribution{s} = self.server.serviceProcess{1, classes{s}}{3};
                end
            else
                try
                    distribution = self.server.serviceProcess{1, class.index}{3};
                catch ME                    
                    distribution = [];
                    line_warning(mfilename,'No distribution is available for the specified class.\n');
                end
            end
        end

        function setService(self, class, distribution, weight)
            % SETSERVICE(CLASS, DISTRIBUTION, WEIGHT)
            % distribution can be a Distribution object or a Workflow object

            if nargin<4 %~exist('weight','var')
                weight=1.0;
            end

            % If Workflow, convert to PH distribution
            if isa(distribution, 'Workflow')
                distribution = distribution.toPH();
            end

            if distribution.isImmediate()
                distribution = Immediate.getInstance();
            end
            if isa(class,'SelfLoopingClass') && class.refstat.index ~= self.index && ~isa(distribution,'Disabled')
                line_error(mfilename, 'For a self-looping class, service cannot be set on stations other than the reference station of the class.');
            end
            if isempty(self.obj)
                server = self.server; % by reference
                c = class.index;
                if length(server.serviceProcess) >= c && ~isempty(server.serviceProcess{1,c}) % if the distribution was already configured
                    % this is a forced state reset in case for example the number of phases changes
                    % appears to run faster without checks, probably due to
                    % isa being slow
                    %oldDistribution = server.serviceProcess{1, c}{3};
                    %isOldMarkovian = isa(oldDistribution,'Markovian');
                    %isNewMarkovian = isa(distribution,'Markovian');
                    %if distribution.getNumParams ~= oldDistribution.getNumParams
                    % %|| (isOldMarkovian && ~isNewMarkovian) || (~isOldMarkovian && isNewMarkovian) || (isOldMarkovian && isNewMarkovian && distribution.getNumberOfPhases ~= oldDistribution.getNumberOfPhases)
                    self.model.setInitialized(false); % this is a better way to invalidate to avoid that sequential calls to setService all trigger an initDefault
                    self.model.hasStruct = false; % invalidate struct so procid gets recomputed
                    self.state=[]; % reset the state vector
                    %end
                else % if first configuration
                    if length(self.classCap) < c
                        self.classCap((length(self.classCap)+1):c) = Inf;
                    end
                    self.setStrategyParam(class, weight);
                    % Default to Drop for finite capacity, WaitingQueue otherwise
                    if self.cap < intmax && ~isinf(self.cap)
                        self.dropRule(c) = DropStrategy.DROP;
                    else
                        self.dropRule(c) = DropStrategy.WAITQ;
                    end
                    server.serviceProcess{1, c}{2} = ServiceStrategy.LI;
                end
                server.serviceProcess{1, c}{3} = distribution;
                self.serviceProcess{c} = distribution;
            else
                self.obj.setService(class.obj, distribution.obj, weight);
                % Also update MATLAB-side storage to keep in sync with Java object
                % This ensures getService returns the correct distribution
                c = class.index;
                self.serviceProcess{c} = distribution;
            end
        end

        function setDelayOff(self, jobclass, setupTime, delayoffTime)
            c = jobclass.index;
            self.setupTime{1, c} = setupTime;
            self.delayoffTime{1, c} = delayoffTime;
        end

        function dist = getSetupTime(self, jobclass)
            c = jobclass.index;
            if c <= length(self.setupTime) && ~isempty(self.setupTime{1, c})
                dist = self.setupTime{1, c};
            else
                dist = [];
            end
        end

        function dist = getDelayOffTime(self, jobclass)
            c = jobclass.index;
            if c <= length(self.delayoffTime) && ~isempty(self.delayoffTime{1, c})
                dist = self.delayoffTime{1, c};
            else
                dist = [];
            end
        end

        function setSwitchover(self, varargin)
            if isempty(self.switchoverTime)
                if SchedStrategy.toId(self.schedStrategy) == SchedStrategy.POLLING
                    self.switchoverTime = cell(1,length(self.model.getClasses()));
                else
                    K = length(self.model.getClasses());
                    self.switchoverTime = cell(K,K);
                    for r=1:K
                        for s=1:K
                            self.switchoverTime{r,s} = Immediate();
                        end
                    end
                end
            end
            if length(varargin)==2
                jobclass = varargin{1};
                soTime = varargin{2};
                % time to switch from queue i to the next one
                if SchedStrategy.toId(self.schedStrategy) ~= SchedStrategy.POLLING
                    line_error(mfilename,'setSwitchover(jobclass, distrib) can only be invoked on queues with SchedStrategy.POLLING.\n');
                end
                c = jobclass.index;
                self.switchoverTime{1,c} = soTime;
            elseif length(varargin)==3
                jobclass_from = varargin{1};
                jobclass_to = varargin{2};
                soTime = varargin{3};
                f = jobclass_from.index;
                t = jobclass_to.index;
                self.switchoverTime{f,t} = soTime;
            end
        end

        function setPollingType(self, rule, par)
            if PollingType.toId(rule) ~= PollingType.KLIMITED
                par = [];
            elseif PollingType.toId(rule) == PollingType.KLIMITED && nargin<3
                line_error(mfilename,'K-Limited polling requires to specify the parameter K, e.g., setPollingType(PollingType.KLIMITED, 2).\n');
            end
            % support only identical polling type at each class buffer
            if SchedStrategy.toId(self.schedStrategy) ~= SchedStrategy.POLLING
                line_error(mfilename,'setPollingType can only be invoked on queues with SchedStrategy.POLLING.\n');
            end
            for r=1:length(self.model.getClasses())
                self.pollingType{1,r} = rule;
                self.pollingPar = par;
                classes = self.model.getClasses();
                setSwitchover(self, classes{r}, Immediate());
            end
        end

        function setPatience(self, class, varargin)
            % SETPATIENCE(CLASS, DISTRIBUTION) - Backwards compatible
            % SETPATIENCE(CLASS, PATIENCETYPE, DISTRIBUTION) - Explicit type
            %
            % Sets the patience type and distribution for a specific job class at this queue.
            % Jobs that wait longer than their patience time will abandon the queue.
            %
            % Parameters:
            %   class        - JobClass object
            %   impatienceType - (Optional) ImpatienceType constant (RENEGING or BALKING)
            %                  If omitted, defaults to ImpatienceType.RENEGING
            %   distribution - Any LINE distribution (Exp, Erlang, HyperExp, etc.)
            %                  excluding modulated processes (BMAP, MAP, MMPP2)
            %
            % Note: This setting takes precedence over the global class patience.
            %
            % Examples:
            %   queue.setPatience(jobclass, Exp(0.2))  % Defaults to RENEGING
            %   queue.setPatience(jobclass, ImpatienceType.RENEGING, Exp(0.2))
            %   queue.setPatience(jobclass, ImpatienceType.BALKING, Exp(0.5))

            % Handle backwards compatibility: 2 or 3 arguments
            if length(varargin) == 1
                % Old signature: setPatience(class, distribution)
                distribution = varargin{1};
                impatienceType = ImpatienceType.RENEGING;  % Default to RENEGING
            elseif length(varargin) == 2
                % New signature: setPatience(class, impatienceType, distribution)
                impatienceType = varargin{1};
                distribution = varargin{2};
            else
                line_error(mfilename, 'Invalid number of arguments. Use setPatience(class, distribution) or setPatience(class, impatienceType, distribution)');
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

            if distribution.isImmediate()
                distribution = Immediate.getInstance();
            end

            if isempty(self.obj)
                c = class.index;
                self.patienceDistributions{1, c} = distribution;
                self.impatienceTypes{1, c} = impatienceType;
            else
                self.obj.setPatience(class.obj, impatienceType, distribution.obj);
            end
        end

        function distribution = getPatience(self, class)
            % DISTRIBUTION = GETPATIENCE(CLASS)
            %
            % Returns the patience distribution for a specific job class.
            % Returns the queue-specific setting if available, otherwise
            % falls back to the global class patience.
            %
            % Parameters:
            %   class - JobClass object
            %
            % Returns:
            %   distribution - The patience distribution, or [] if not set

            if isempty(self.obj)
                c = class.index;
                % Check queue-specific patience first
                if c <= length(self.patienceDistributions) && ~isempty(self.patienceDistributions{1, c})
                    distribution = self.patienceDistributions{1, c};
                else
                    % Fall back to global class patience
                    distribution = class.getPatience();
                end
            else
                distObj = self.obj.getPatience(class.obj);
                if isempty(distObj)
                    distribution = [];
                else
                    distribution = Distribution.fromJavaObject(distObj);
                end
            end
        end

        function setLimit(self, limit)
            % SETLIMIT(LIMIT) Sets the maximum number of jobs for LPS scheduling
            %
            % Parameters:
            %   limit - Maximum number of jobs in PS (processor sharing) mode for LPS

            if isempty(self.obj)
                % MATLAB native implementation - store as a scheduling parameter
                if SchedStrategy.toId(self.schedStrategy) ~= SchedStrategy.LPS
                    line_warning(mfilename, 'setLimit is only applicable to LPS (Least Progress Scheduling) queues.');
                    return;
                end
                % Store limit in schedStrategyPar (use index 0 for queue-level parameter)
                if length(self.schedStrategyPar) < 1
                    self.schedStrategyPar = zeros(1, 1);
                end
                self.schedStrategyPar(1) = limit;
            else
                % JavaNative implementation
                if SchedStrategy.toId(self.schedStrategy) ~= SchedStrategy.LPS
                    line_warning(mfilename, 'setLimit is only applicable to LPS (Least Progress Scheduling) queues.');
                    return;
                end
                self.obj.setLimit(limit);
            end
        end

        function limit = getLimit(self)
            % LIMIT = GETLIMIT() Returns the maximum number of jobs for LPS scheduling
            %
            % Returns:
            %   limit - Maximum number of jobs in PS mode for LPS

            if isempty(self.obj)
                % MATLAB native implementation
                if length(self.schedStrategyPar) >= 1
                    limit = self.schedStrategyPar(1);
                else
                    limit = [];
                end
            else
                % JavaNative implementation
                limit = self.obj.getLimit();
            end
        end

        function impatienceType = getImpatienceType(self, class)
            % IMPATIENCETYPE = GETIMPATIENCETYPE(CLASS)
            %
            % Returns the impatience type for a specific job class.
            % Returns the queue-specific setting if available, otherwise
            % falls back to the global class impatience type.
            %
            % Parameters:
            %   class - JobClass object
            %
            % Returns:
            %   impatienceType - The impatience type (ImpatienceType constant), or [] if not set

            if isempty(self.obj)
                c = class.index;
                % Check queue-specific impatience type first
                if c <= length(self.impatienceTypes) && ~isempty(self.impatienceTypes{1, c})
                    impatienceType = self.impatienceTypes{1, c};
                else
                    % Fall back to global class impatience type
                    impatienceType = class.getImpatienceType();
                end
            else
                impatienceTypeId = self.obj.getImpatienceType(class.obj);
                if isempty(impatienceTypeId)
                    impatienceType = [];
                else
                    impatienceType = ImpatienceType.fromId(impatienceTypeId.getID());
                end
            end
        end

        function tf = hasPatience(self, class)
            % TF = HASPATIENCE(CLASS)
            %
            % Returns true if this class has patience configured at this queue
            % (either locally or globally).

            dist = self.getPatience(class);
            tf = ~isempty(dist) && ~isa(dist, 'Disabled');
        end

        function setBalking(self, class, strategy, thresholds)
            % SETBALKING(CLASS, STRATEGY, THRESHOLDS)
            %
            % Configures balking behavior for a specific job class at this queue.
            % When a customer arrives, they may refuse to join based on queue length.
            %
            % Parameters:
            %   class      - JobClass object
            %   strategy   - BalkingStrategy constant:
            %                QUEUE_LENGTH - Balk based on current queue length
            %                EXPECTED_WAIT - Balk based on expected waiting time
            %                COMBINED - Both conditions (OR logic)
            %   thresholds - Cell array of balking thresholds, each element is:
            %                {minJobs, maxJobs, probability}
            %                where probability is the chance to balk when queue
            %                length is in [minJobs, maxJobs] range.
            %
            % Example:
            %   % Balk with 30% probability when 5-10 jobs in queue,
            %   % 80% when 11-20 jobs, 100% when >20 jobs
            %   queue.setBalking(jobclass, BalkingStrategy.QUEUE_LENGTH, ...
            %       {{5, 10, 0.3}, {11, 20, 0.8}, {21, Inf, 1.0}});

            if isempty(self.obj)
                c = class.index;
                self.balkingStrategies{1, c} = strategy;
                self.balkingThresholds{1, c} = thresholds;
            else
                % Java native - convert thresholds to Java format
                jThresholds = jline.util.BalkingThresholdList();
                for i = 1:length(thresholds)
                    th = thresholds{i};
                    minJobs = th{1};
                    maxJobs = th{2};
                    if isinf(maxJobs)
                        maxJobs = java.lang.Integer.MAX_VALUE;
                    end
                    probability = th{3};
                    jThresholds.add(jline.lang.BalkingThreshold(minJobs, maxJobs, probability));
                end
                % Convert strategy to Java enum
                switch strategy
                    case BalkingStrategy.QUEUE_LENGTH
                        jStrategy = jline.lang.constant.BalkingStrategy.QUEUE_LENGTH;
                    case BalkingStrategy.EXPECTED_WAIT
                        jStrategy = jline.lang.constant.BalkingStrategy.EXPECTED_WAIT;
                    case BalkingStrategy.COMBINED
                        jStrategy = jline.lang.constant.BalkingStrategy.COMBINED;
                end
                self.obj.setBalking(class.obj, jStrategy, jThresholds);
            end
        end

        function [strategy, thresholds] = getBalking(self, class)
            % [STRATEGY, THRESHOLDS] = GETBALKING(CLASS)
            %
            % Returns the balking configuration for a specific job class.
            %
            % Parameters:
            %   class - JobClass object
            %
            % Returns:
            %   strategy   - BalkingStrategy constant, or [] if not configured
            %   thresholds - Cell array of {minJobs, maxJobs, probability} tuples

            if isempty(self.obj)
                c = class.index;
                if c <= length(self.balkingStrategies) && ~isempty(self.balkingStrategies{1, c})
                    strategy = self.balkingStrategies{1, c};
                    thresholds = self.balkingThresholds{1, c};
                else
                    strategy = [];
                    thresholds = {};
                end
            else
                jStrategy = self.obj.getBalkingStrategy(class.obj);
                if isempty(jStrategy)
                    strategy = [];
                    thresholds = {};
                else
                    strategy = BalkingStrategy.fromId(jStrategy.getId());
                    jThresholds = self.obj.getBalkingThresholds(class.obj);
                    thresholds = {};
                    if ~isempty(jThresholds)
                        for i = 0:(jThresholds.size()-1)
                            jTh = jThresholds.get(i);
                            maxJobs = jTh.getMaxJobs();
                            if maxJobs == java.lang.Integer.MAX_VALUE
                                maxJobs = Inf;
                            end
                            thresholds{end+1} = {jTh.getMinJobs(), maxJobs, jTh.getProbability()};
                        end
                    end
                end
            end
        end

        function tf = hasBalking(self, class)
            % TF = HASBALKING(CLASS)
            %
            % Returns true if this class has balking configured at this queue.

            [strategy, ~] = self.getBalking(class);
            tf = ~isempty(strategy);
        end

        function setRetrial(self, class, delayDistribution, maxAttempts)
            % SETRETRIAL(CLASS, DELAYDISTRIBUTION, MAXATTEMPTS)
            %
            % Configures retrial behavior for a specific job class at this queue.
            % When a customer is rejected (queue full), they move to an orbit
            % and retry after a random delay.
            %
            % Parameters:
            %   class             - JobClass object
            %   delayDistribution - Distribution for retrial delay (e.g., Exp(0.5))
            %   maxAttempts       - Maximum number of retrial attempts:
            %                       -1 = unlimited retries (default)
            %                       N  = drop after N failed attempts
            %
            % Example:
            %   % Retry with exponential delay, unlimited attempts
            %   queue.setRetrial(jobclass, Exp(0.5), -1);
            %
            %   % Retry up to 3 times with Erlang delay
            %   queue.setRetrial(jobclass, Erlang(2, 0.3), 3);

            if nargin < 4
                maxAttempts = -1;  % Unlimited by default
            end

            if isa(delayDistribution, 'BMAP') || isa(delayDistribution, 'MAP') || isa(delayDistribution, 'MMPP2')
                line_error(mfilename, 'Modulated processes (BMAP, MAP, MMPP2) are not supported for retrial delay distributions.');
            end

            if isempty(self.obj)
                c = class.index;
                self.retrialDelays{1, c} = delayDistribution;
                % Ensure array is large enough
                if length(self.retrialMaxAttempts) < c
                    self.retrialMaxAttempts(end+1:c) = -1;
                end
                self.retrialMaxAttempts(c) = maxAttempts;
                % Also set drop rule to RETRIAL or RETRIAL_WITH_LIMIT
                if maxAttempts < 0
                    self.dropRule(c) = DropStrategy.RETRIAL;
                else
                    self.dropRule(c) = DropStrategy.RETRIAL_WITH_LIMIT;
                end
            else
                self.obj.setRetrial(class.obj, delayDistribution.obj, maxAttempts);
            end
        end

        function [delayDistribution, maxAttempts] = getRetrial(self, class)
            % [DELAYDISTRIBUTION, MAXATTEMPTS] = GETRETRIAL(CLASS)
            %
            % Returns the retrial configuration for a specific job class.
            %
            % Parameters:
            %   class - JobClass object
            %
            % Returns:
            %   delayDistribution - Retrial delay distribution, or [] if not configured
            %   maxAttempts       - Maximum retrial attempts (-1 = unlimited)

            if isempty(self.obj)
                c = class.index;
                if c <= length(self.retrialDelays) && ~isempty(self.retrialDelays{1, c})
                    delayDistribution = self.retrialDelays{1, c};
                    if c <= length(self.retrialMaxAttempts)
                        maxAttempts = self.retrialMaxAttempts(c);
                    else
                        maxAttempts = -1;
                    end
                else
                    delayDistribution = [];
                    maxAttempts = -1;
                end
            else
                distObj = self.obj.getRetrialDelayDistribution(class.obj);
                if isempty(distObj)
                    delayDistribution = [];
                    maxAttempts = -1;
                else
                    delayDistribution = Distribution.fromJavaObject(distObj);
                    maxAttempts = self.obj.getMaxRetrialAttempts(class.obj);
                end
            end
        end

        function tf = hasRetrial(self, class)
            % TF = HASRETRIAL(CLASS)
            %
            % Returns true if this class has retrial configured at this queue.

            [dist, ~] = self.getRetrial(class);
            tf = ~isempty(dist) && ~isa(dist, 'Disabled');
        end

        function setOrbitImpatience(self, class, distribution)
            % SETORBITIMPATIENCE(CLASS, DISTRIBUTION)
            %
            % Sets the impatience (abandonment) rate for customers in the orbit.
            % This is separate from queue patience (reneging from waiting queue).
            % Used in BMAP/PH/N/N retrial queues where customers in the orbit
            % may abandon before successfully retrying.
            %
            % Parameters:
            %   class        - JobClass object
            %   distribution - Distribution for orbit abandonment time (e.g., Exp(gamma))
            %
            % Example:
            %   queue.setOrbitImpatience(jobclass, Exp(0.008));  % gamma = 0.008

            if isa(distribution, 'BMAP') || isa(distribution, 'MAP') || isa(distribution, 'MMPP2')
                line_error(mfilename, 'Modulated processes (BMAP, MAP, MMPP2) are not supported for orbit impatience distributions.');
            end

            if isempty(self.obj)
                c = class.index;
                self.orbitImpatienceDistributions{1, c} = distribution;
            else
                self.obj.setOrbitImpatience(class.obj, distribution.obj);
            end
        end

        function distribution = getOrbitImpatience(self, class)
            % DISTRIBUTION = GETORBITIMPATIENCE(CLASS)
            %
            % Returns the orbit impatience distribution for a specific job class.
            %
            % Parameters:
            %   class - JobClass object
            %
            % Returns:
            %   distribution - The orbit impatience distribution, or [] if not set

            if isempty(self.obj)
                c = class.index;
                if c <= length(self.orbitImpatienceDistributions) && ~isempty(self.orbitImpatienceDistributions{1, c})
                    distribution = self.orbitImpatienceDistributions{1, c};
                else
                    distribution = [];
                end
            else
                distObj = self.obj.getOrbitImpatience(class.obj);
                if isempty(distObj)
                    distribution = [];
                else
                    distribution = Distribution.fromJavaObject(distObj);
                end
            end
        end

        function tf = hasOrbitImpatience(self, class)
            % TF = HASORBITORBITIMPATIENCE(CLASS)
            %
            % Returns true if this class has orbit impatience configured at this queue.

            dist = self.getOrbitImpatience(class);
            tf = ~isempty(dist) && ~isa(dist, 'Disabled');
        end

        function setBatchRejectProbability(self, class, p)
            % SETBATCHREJECTPROBABILITY(CLASS, P)
            %
            % Sets the probability that an entire batch is rejected when it
            % cannot be fully admitted. Used in BMAP/PH/N/N retrial queues
            % with batch arrivals.
            %
            % When a batch of size k arrives and only m < k servers are free:
            %   - With probability p: entire batch is rejected to orbit
            %   - With probability (1-p): m customers are admitted, k-m go to orbit
            %
            % Parameters:
            %   class - JobClass object
            %   p     - Probability [0,1] that batch is rejected vs partially admitted
            %           Default is 0 (partial admission allowed)
            %
            % Example:
            %   queue.setBatchRejectProbability(jobclass, 0.4);

            if p < 0 || p > 1
                line_error(mfilename, 'Batch reject probability must be in [0, 1].');
            end

            if isempty(self.obj)
                c = class.index;
                % Ensure array is large enough
                if length(self.batchRejectProb) < c
                    self.batchRejectProb(end+1:c) = 0;
                end
                self.batchRejectProb(c) = p;
            else
                self.obj.setBatchRejectProbability(class.obj, p);
            end
        end

        function p = getBatchRejectProbability(self, class)
            % P = GETBATCHREJECTPROBABILITY(CLASS)
            %
            % Returns the batch reject probability for a specific job class.
            %
            % Parameters:
            %   class - JobClass object
            %
            % Returns:
            %   p - Batch reject probability [0,1], or 0 if not set

            if isempty(self.obj)
                c = class.index;
                if c <= length(self.batchRejectProb) && self.batchRejectProb(c) > 0
                    p = self.batchRejectProb(c);
                else
                    p = 0;  % Default: partial admission allowed
                end
            else
                p = self.obj.getBatchRejectProbability(class.obj);
            end
        end

        %        function distrib = getServiceProcess(self, oclass)
        %            distrib = self.serviceProcess{oclass};
        %        end

        % ==================== Heterogeneous Server Methods ====================

        function self = addServerType(self, serverType)
            % ADDSERVERTYPE Add a server type to this queue
            %
            % self = ADDSERVERTYPE(serverType) adds a ServerType to this queue
            % for heterogeneous multiserver configuration.
            %
            % When server types are added, the queue becomes a heterogeneous
            % multiserver queue where different server types can have different
            % service rates and serve different subsets of job classes.
            %
            % @param serverType The ServerType object to add

            if isempty(serverType)
                line_error(mfilename, 'Server type cannot be empty');
            end

            % Check if already added
            for i = 1:length(self.serverTypes)
                if self.serverTypes{i} == serverType
                    line_error(mfilename, 'Server type ''%s'' is already added to this queue', serverType.getName());
                end
            end

            if isempty(self.obj)
                % MATLAB native implementation
                serverType.setId(length(self.serverTypes));
                serverType.setParentQueue(self);
                self.serverTypes{end+1} = serverType;

                % Initialize service distribution map for this server type
                self.heteroServiceDistributions(serverType.getName()) = containers.Map();

                % Update total number of servers
                self.updateTotalServerCount();
            else
                % Java native - delegate to Java object
                self.obj.addServerType(serverType.obj);
                % Also store locally
                self.serverTypes{end+1} = serverType;
            end
        end

        function updateTotalServerCount(self)
            % UPDATETOTALSERVERCOUNT Update total server count from all types
            %
            % Internal method to recalculate numberOfServers.

            if isempty(self.serverTypes)
                return;
            end
            total = 0;
            for i = 1:length(self.serverTypes)
                total = total + self.serverTypes{i}.getNumOfServers();
            end
            self.numberOfServers = total;
        end

        function types = getServerTypes(self)
            % GETSERVERTYPES Get the list of server types
            %
            % types = GETSERVERTYPES() returns a cell array of ServerType objects.

            types = self.serverTypes;
        end

        function n = getNumServerTypes(self)
            % GETNUMSERVERTYPES Get the number of server types
            %
            % n = GETNUMSERVERTYPES() returns the number of server types,
            % or 0 if this is a homogeneous queue.

            n = length(self.serverTypes);
        end

        function result = isHeterogeneous(self)
            % ISHETEROGENEOUS Check if this is a heterogeneous multiserver queue
            %
            % result = ISHETEROGENEOUS() returns true if server types are defined.

            result = ~isempty(self.serverTypes);
        end

        function self = setHeteroSchedPolicy(self, policy)
            % SETHETEROSCHEDPOLICY Set the heterogeneous server scheduling policy
            %
            % self = SETHETEROSCHEDPOLICY(policy) sets the policy that determines
            % how jobs are assigned to server types when a job's class is
            % compatible with multiple server types.
            %
            % @param policy HeteroSchedPolicy constant (ORDER, ALIS, ALFS, FAIRNESS, FSF, RAIS)

            if isempty(self.obj)
                self.heteroSchedPolicy = policy;
            else
                % Convert to Java enum
                switch policy
                    case HeteroSchedPolicy.ORDER
                        jPolicy = jline.lang.constant.HeteroSchedPolicy.ORDER;
                    case HeteroSchedPolicy.ALIS
                        jPolicy = jline.lang.constant.HeteroSchedPolicy.ALIS;
                    case HeteroSchedPolicy.ALFS
                        jPolicy = jline.lang.constant.HeteroSchedPolicy.ALFS;
                    case HeteroSchedPolicy.FAIRNESS
                        jPolicy = jline.lang.constant.HeteroSchedPolicy.FAIRNESS;
                    case HeteroSchedPolicy.FSF
                        jPolicy = jline.lang.constant.HeteroSchedPolicy.FSF;
                    case HeteroSchedPolicy.RAIS
                        jPolicy = jline.lang.constant.HeteroSchedPolicy.RAIS;
                end
                self.obj.setHeteroSchedPolicy(jPolicy);
                self.heteroSchedPolicy = policy;
            end
        end

        function policy = getHeteroSchedPolicy(self)
            % GETHETEROSCHEDPOLICY Get the heterogeneous server scheduling policy
            %
            % policy = GETHETEROSCHEDPOLICY() returns the HeteroSchedPolicy.

            policy = self.heteroSchedPolicy;
        end

        function setHeteroService(self, jobClass, serverType, distribution)
            % SETHETEROSERVICE Set service distribution for a job class and server type
            %
            % SETHETEROSERVICE(jobClass, serverType, distribution) sets the
            % service time distribution for a specific job class when served
            % by a specific server type.
            %
            % @param jobClass The JobClass object
            % @param serverType The ServerType object
            % @param distribution The service time Distribution

            if isempty(jobClass)
                line_error(mfilename, 'Job class cannot be empty');
            end
            if isempty(serverType)
                line_error(mfilename, 'Server type cannot be empty');
            end
            if isempty(distribution)
                line_error(mfilename, 'Distribution cannot be empty');
            end

            % Check if server type is in this queue
            found = false;
            for i = 1:length(self.serverTypes)
                if self.serverTypes{i} == serverType
                    found = true;
                    break;
                end
            end
            if ~found
                line_error(mfilename, 'Server type ''%s'' is not added to this queue. Call addServerType() first.', serverType.getName());
            end

            if isempty(self.obj)
                % MATLAB native implementation
                if ~isKey(self.heteroServiceDistributions, serverType.getName())
                    self.heteroServiceDistributions(serverType.getName()) = containers.Map();
                end
                classMap = self.heteroServiceDistributions(serverType.getName());
                classMap(jobClass.getName()) = distribution;
                self.heteroServiceDistributions(serverType.getName()) = classMap;

                % Ensure compatibility
                if ~serverType.isCompatible(jobClass)
                    serverType.addCompatible(jobClass);
                end
            else
                % Java native - delegate to Java object
                self.obj.setService(jobClass.obj, serverType.obj, distribution.obj);
            end
        end

        function distribution = getHeteroService(self, jobClass, serverType)
            % GETHETEROSERVICE Get service distribution for a job class and server type
            %
            % distribution = GETHETEROSERVICE(jobClass, serverType) returns the
            % service time distribution for a specific job class and server type.
            %
            % @param jobClass The JobClass object
            % @param serverType The ServerType object
            % @return distribution The service time Distribution, or [] if not set

            if isempty(self.obj)
                if isKey(self.heteroServiceDistributions, serverType.getName())
                    classMap = self.heteroServiceDistributions(serverType.getName());
                    if isKey(classMap, jobClass.getName())
                        distribution = classMap(jobClass.getName());
                    else
                        distribution = [];
                    end
                else
                    distribution = [];
                end
            else
                distObj = self.obj.getService(jobClass.obj, serverType.obj);
                if isempty(distObj)
                    distribution = [];
                else
                    distribution = Distribution.fromJavaObject(distObj);
                end
            end
        end

        function st = getServerTypeById(self, id)
            % GETSERVERTYPEBYID Get a server type by its ID
            %
            % st = GETSERVERTYPEBYID(id) returns the ServerType with the given ID,
            % or [] if not found.

            if id >= 0 && id < length(self.serverTypes)
                st = self.serverTypes{id + 1};  % MATLAB 1-indexed
            else
                st = [];
            end
        end

        function st = getServerTypeByName(self, name)
            % GETSERVERTYPEBYNAME Get a server type by its name
            %
            % st = GETSERVERTYPEBYNAME(name) returns the ServerType with the given
            % name, or [] if not found.

            st = [];
            for i = 1:length(self.serverTypes)
                if strcmp(self.serverTypes{i}.getName(), name)
                    st = self.serverTypes{i};
                    return;
                end
            end
        end

        function result = validateCompatibility(self)
            % VALIDATECOMPATIBILITY Check all job classes have compatible server types
            %
            % result = VALIDATECOMPATIBILITY() returns true if all job classes
            % in the model have at least one compatible server type at this queue.

            if ~self.isHeterogeneous()
                result = true;
                return;
            end

            classes = self.model.getClasses();
            for c = 1:length(classes)
                jobClass = classes{c};
                hasCompatible = false;
                for s = 1:length(self.serverTypes)
                    if self.serverTypes{s}.isCompatible(jobClass)
                        hasCompatible = true;
                        break;
                    end
                end
                if ~hasCompatible
                    result = false;
                    return;
                end
            end
            result = true;
        end

        % ==================== Immediate Feedback Methods ====================

        function setImmediateFeedback(self, varargin)
            % SETIMMEDIATEFEEDBACK Set immediate feedback for self-loops
            %
            % SETIMMEDIATEFEEDBACK(true) enables immediate feedback for all classes
            % SETIMMEDIATEFEEDBACK(false) disables immediate feedback for all classes
            % SETIMMEDIATEFEEDBACK(jobClass) enables for a specific class
            % SETIMMEDIATEFEEDBACK({class1, class2}) enables for multiple classes
            %
            % When enabled, a job that self-loops at this station stays in service
            % instead of going back to the queue.

            if isempty(self.obj)
                % MATLAB native implementation
                if nargin == 2
                    arg = varargin{1};
                    if islogical(arg) || isnumeric(arg)
                        if arg
                            % Enable for all classes
                            self.immediateFeedback = 'all';
                        else
                            % Disable for all classes
                            self.immediateFeedback = {};
                        end
                    elseif isa(arg, 'JobClass')
                        % Single class
                        if isempty(self.immediateFeedback) || ischar(self.immediateFeedback)
                            self.immediateFeedback = {};
                        end
                        if ~any(cellfun(@(x) x == arg.index, self.immediateFeedback))
                            self.immediateFeedback{end+1} = arg.index;
                        end
                    elseif iscell(arg)
                        % Cell array of classes
                        self.immediateFeedback = {};
                        for i = 1:length(arg)
                            if isa(arg{i}, 'JobClass')
                                self.immediateFeedback{end+1} = arg{i}.index;
                            end
                        end
                    end
                end
            else
                % Java native implementation
                if nargin == 2
                    arg = varargin{1};
                    if islogical(arg) || isnumeric(arg)
                        self.obj.setImmediateFeedback(logical(arg));
                    elseif isa(arg, 'JobClass')
                        self.obj.setImmediateFeedback(arg.obj);
                    elseif iscell(arg)
                        classList = java.util.ArrayList();
                        for i = 1:length(arg)
                            if isa(arg{i}, 'JobClass')
                                classList.add(arg{i}.obj);
                            end
                        end
                        self.obj.setImmediateFeedbackForClasses(classList);
                    end
                end
            end
        end

        function tf = hasImmediateFeedback(self, varargin)
            % HASIMMEDIATEFEEDBACK Check if immediate feedback is enabled
            %
            % TF = HASIMMEDIATEFEEDBACK() returns true if enabled for any class
            % TF = HASIMMEDIATEFEEDBACK(jobClass) returns true if enabled for specific class

            if isempty(self.obj)
                % MATLAB native implementation
                if isempty(self.immediateFeedback)
                    tf = false;
                elseif ischar(self.immediateFeedback) && strcmp(self.immediateFeedback, 'all')
                    tf = true;
                elseif nargin == 1
                    % No class specified - check if any class has it enabled
                    tf = ~isempty(self.immediateFeedback);
                else
                    % Check specific class
                    jobClass = varargin{1};
                    if ischar(self.immediateFeedback) && strcmp(self.immediateFeedback, 'all')
                        tf = true;
                    else
                        tf = any(cellfun(@(x) x == jobClass.index, self.immediateFeedback));
                    end
                end
            else
                % Java native implementation
                if nargin == 1
                    tf = self.obj.hasImmediateFeedback();
                else
                    jobClass = varargin{1};
                    tf = self.obj.hasImmediateFeedback(jobClass.index - 1);  % Java 0-indexed
                end
            end
        end

        function classes = getImmediateFeedbackClasses(self)
            % GETIMMEDIATEFEEDBACKCLASSES Get list of class indices with immediate feedback
            %
            % CLASSES = GETIMMEDIATEFEEDBACKCLASSES() returns cell array of class indices

            if isempty(self.obj)
                if isempty(self.immediateFeedback)
                    classes = {};
                elseif ischar(self.immediateFeedback) && strcmp(self.immediateFeedback, 'all')
                    classes = 'all';
                else
                    classes = self.immediateFeedback;
                end
            else
                jClasses = self.obj.getImmediateFeedbackClasses();
                if isempty(jClasses)
                    classes = {};
                elseif jClasses.equals("all")
                    classes = 'all';
                else
                    classes = {};
                    for i = 0:(jClasses.size()-1)
                        classes{end+1} = jClasses.get(i) + 1;  % Convert to MATLAB 1-indexed
                    end
                end
            end
        end

    end
end
