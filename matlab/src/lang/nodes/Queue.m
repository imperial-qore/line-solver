classdef Queue < Station
    % A service station with queueing
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.

    properties
        schedPolicy;
        schedStrategy;
        schedStrategyPar;
        serviceProcess;
        setupTime;
        delayoffTime;
        switchoverTime;
        pollingType;
        pollingPar;
    end

    methods
        %Constructor
        function self = Queue(model, name, schedStrategy)
            % SELF = QUEUE(MODEL, NAME, SCHEDSTRATEGY)

            self@Station(name);

            if isa(model,'Network')
                classes = model.classes;
                self.input = Buffer(classes);
                self.output = Dispatcher(classes);
                self.schedPolicy = SchedStrategyType.PR;
                self.schedStrategy = SchedStrategy.PS;
                self.serviceProcess = {};
                self.server = Server(classes);
                self.numberOfServers = 1;
                self.schedStrategyPar = zeros(1,length(model.classes));
                self.setModel(model);
                self.model.addNode(self);
                self.dropRule = [];
                self.obj = [];
                self.setupTime = {};
                self.delayoffTime = {};
                self.pollingType = {};
                self.switchoverTime  = {};

                if nargin>=3 %exist('schedStrategy','var')
                    self.schedStrategy = schedStrategy;
                    switch SchedStrategy.toId(self.schedStrategy)
                        case {SchedStrategy.ID_PS, SchedStrategy.ID_DPS,SchedStrategy.ID_GPS}
                            self.schedPolicy = SchedStrategyType.PR;
                            self.server = SharedServer(classes);
                        case SchedStrategy.ID_LCFSPR
                            self.schedPolicy = SchedStrategyType.PR;
                            self.server = PreemptiveServer(classes);
                        case {SchedStrategy.ID_FCFS, SchedStrategy.ID_LCFS, SchedStrategy.ID_SIRO, SchedStrategy.ID_SEPT, SchedStrategy.ID_LEPT, SchedStrategy.ID_SJF, SchedStrategy.ID_LJF}
                            self.schedPolicy = SchedStrategyType.NP;
                            self.server = Server(classes);
                        case SchedStrategy.ID_INF
                            self.schedPolicy = SchedStrategyType.NP;
                            self.server = InfiniteServer(classes);
                            self.numberOfServers = Inf;
                        case SchedStrategy.ID_HOL
                            self.schedPolicy = SchedStrategyType.NP;
                            self.server = Server(classes);
                        case SchedStrategy.ID_POLLING
                            self.schedPolicy = SchedStrategyType.NP;
                            self.server = PollingServer(classes);
                        otherwise
                            line_error(mfilename,sprintf('The specified scheduling strategy (%s) is unsupported.',schedStrategy));
                    end
                end
            elseif isa(model,'JNetwork')
                self.setModel(model);
                switch SchedStrategy.toId(schedStrategy)
                    case SchedStrategy.ID_INF
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.INF);
                    case SchedStrategy.ID_FCFS
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.FCFS);
                    case SchedStrategy.ID_LCFS
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.LCFS);
                    case SchedStrategy.ID_SIRO
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.SIRO);
                    case SchedStrategy.ID_SJF
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.SJF);
                    case SchedStrategy.ID_LJF
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.LJF);
                    case SchedStrategy.ID_PS
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.PS);
                    case SchedStrategy.ID_DPS
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.DPS);
                    case SchedStrategy.ID_GPS
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.GPS);
                    case SchedStrategy.ID_SEPT
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.SEPT);
                    case SchedStrategy.ID_LEPT
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.LEPT);
                    case SchedStrategy.ID_HOL
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.HOL);
                    case SchedStrategy.ID_FORK
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.FORK);
                    case SchedStrategy.ID_EXT
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.EXT);
                    case SchedStrategy.ID_REF
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.REF);
                    case SchedStrategy.ID_LCFSPR
                        self.obj = jline.lang.nodes.Queue(model.obj, name, jline.lang.constant.SchedStrategy.LCFSPR);
                end
                self.obj.setNumberOfServers(1);
            end
        end

        function setLoadDependence(self, alpha)
            switch SchedStrategy.toId(self.schedStrategy)
                case {SchedStrategy.ID_PS, SchedStrategy.ID_FCFS}
                    setLimitedLoadDependence(self, alpha);
                otherwise
                    line_error(mfilename,'Load-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) stations.');
            end
        end

        function setClassDependence(self, beta)
            switch SchedStrategy.toId(self.schedStrategy)
                case {SchedStrategy.ID_PS, SchedStrategy.ID_FCFS}
                    setLimitedClassDependence(self, beta);
                otherwise
                    line_error(mfilename,'Class-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) stations.');
            end
        end

        function setNumberOfServers(self, value)
            % SETNUMBEROFSERVERS(VALUE)
            if isempty(self.obj)
                switch SchedStrategy.toId(self.schedStrategy)
                    case SchedStrategy.ID_INF
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
                    case {SchedStrategy.ID_DPS, SchedStrategy.ID_GPS}
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

            self.schedStrategyPar(class.index) = weight;
        end

        function distribution = getService(self, class)
            % DISTRIBUTION = GETSERVICE(CLASS)

            % return the service distribution assigned to the given class
            if nargin<2 %~exist('class','var')
                for s = 1:length(self.model.classes)
                    distribution{s} = self.server.serviceProcess{1, self.model.classes{s}}{3};
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

        function distrib = getServiceProcess(self, oclass)
            distrib = self.getService{oclass};
        end

        function setService(self, class, distribution, weight)
            % SETSERVICE(CLASS, DISTRIBUTION, WEIGHT)
            if nargin<4 %~exist('weight','var')
                weight=1.0;
            end
            if distribution.isImmediate()
                distribution = Immediate.getInstance();
            end
            if isempty(self.obj)
                server = self.server; % by reference
                c = class.index;
                if length(server.serviceProcess) >= c && ~isempty(server.serviceProcess{1,c}) % if the distribution was already configured
                    % this is a forced state reset in case for example the number of phases changes
                    % appears to run faster without checks, probably due to
                    % isa being slow
                    %oldDistribution = server.serviceProcess{1, c}{3};
                    %isOldMarkovian = isa(oldDistribution,'MarkovianDistribution');
                    %isNewMarkovian = isa(distribution,'MarkovianDistribution');
                    %if distribution.getNumParams ~= oldDistribution.getNumParams
                    % %|| (isOldMarkovian && ~isNewMarkovian) || (~isOldMarkovian && isNewMarkovian) || (isOldMarkovian && isNewMarkovian && distribution.getNumberOfPhases ~= oldDistribution.getNumberOfPhases)
                    self.model.setInitialized(false); % this is a better way to invalidate to avoid that sequential calls to setService all trigger an initDefault
                    self.state=[]; % reset the state vector
                    %end
                else % if first configuration
                    if length(self.classCap) < c
                        self.classCap((length(self.classCap)+1):c) = Inf;
                    end
                    self.input.inputJobClasses{c} = {class, self.schedPolicy, DropStrategy.WaitingQueue};
                    self.setStrategyParam(class, weight);
                    self.dropRule(c) = DropStrategy.ID_WAITQ;
                    server.serviceProcess{1, c}{2} = ServiceStrategy.LI;
                end
                server.serviceProcess{1, c}{3} = distribution;
                self.serviceProcess{c} = distribution;
            else
                self.obj.setService(class.obj, distribution.obj, weight);
            end
        end

        function setDelayOff(self, jobclass, setupTime, delayoffTime)
            c = jobclass.index;
            self.setupTime{1, c} = setupTime;
            self.delayoffTime{1, c} = delayoffTime;
        end

        function sections = getSections(self)
            % SECTIONS = GETSECTIONS()

            sections = {self.input, self.server, self.output};
        end

        function setSwitchover(self, varargin)
            if isempty(self.switchoverTime)
                if SchedStrategy.toId(self.schedStrategy) == SchedStrategy.ID_POLLING
                    self.switchoverTime = cell(1,length(self.model.classes));
                else
                    K = length(self.model.classes);
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
                if SchedStrategy.toId(self.schedStrategy) ~= SchedStrategy.ID_POLLING
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
            if PollingType.toId(rule) ~= PollingType.ID_KLIMITED
                par = [];
            elseif PollingType.toId(rule) == PollingType.ID_KLIMITED && nargin<3
                    line_error(mfilename,'K-Limited polling requires to specify the parameter K, e.g., setPollingType(PollingType.KLIMITED, 2).\n');
            end
            % support only identical polling type at each class buffer
            if SchedStrategy.toId(self.schedStrategy) ~= SchedStrategy.ID_POLLING
                line_error(mfilename,'setPollingType can only be invoked on queues with SchedStrategy.POLLING.\n');
            end
            for r=1:length(self.model.classes)
                self.pollingType{1,r} = rule;
                self.pollingPar = par;
                setSwitchover(self, self.model.classes{r}, Immediate());
            end
        end

        %        function distrib = getServiceProcess(self, oclass)
        %            distrib = self.serviceProcess{oclass};
        %        end

    end
end
