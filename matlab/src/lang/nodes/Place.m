classdef Place < Station
    % 
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        schedStrategies;
        schedStrategy;
        schedStrategyPar;
        queueing;             % true if this is a queueing place (has an embedded queue)
        serviceProcess;       % per-class service (queue) processes; empty for ordinary places
        departureDiscipline;  % per-class depository departure discipline
    end

    methods
        function self = Place(model,name,schedStrategy)
            % PLACE(MODEL, NAME)
            % PLACE(MODEL, NAME, SCHEDSTRATEGY)
            % SCHEDSTRATEGY (optional) turns the place into a queueing place with
            % an embedded queue served under the given scheduling strategy once a
            % service process is assigned via setService.

            self@Station(name);
            if model.isMatlabNative()
                classes = model.getClasses();
                self.input = Storage(classes);
                self.output = Linkage(classes);
                self.setModel(model);
                self.model.addNode(self);
                self.server = ServiceTunnel();

                %numOfClasses = [];
                % Places have infinite capacity (like delay nodes)
                self.numberOfServers = Inf;
                self.schedStrategy = SchedStrategy.INF;
                self.schedStrategyPar = [];

                self.classCap = [];
                self.cap = [];
                self.schedStrategies = [];
                self.dropRule = [];
                self.queueing = false;
                self.serviceProcess = {};
                self.departureDiscipline = [];
            elseif model.isJavaNative()
                self.setModel(model);
                if nargin>=3 && ~isempty(schedStrategy)
                    self.obj = jline.lang.nodes.Place(model.obj, name, jline.lang.constant.SchedStrategy.fromText(SchedStrategy.toText(schedStrategy)));
                else
                    self.obj = jline.lang.nodes.Place(model.obj, name);
                end
                self.index = model.obj.getNodeIndex(self.obj);
            end
            if nargin>=3 && ~isempty(schedStrategy)
                self.schedStrategy = schedStrategy;
            end
        end

        function init(self)
            numOfClasses = length(self.model.getClasses());
            % Preserve the embedded-queue scheduling strategy of a queueing
            % place (set via the constructor); only ordinary places, whose
            % strategy is nominal, default to FCFS here. Clobbering it would
            % turn e.g. an infinite-server think place into a single-server
            % FCFS station on the next refresh.
            if ~self.queueing
                self.schedStrategy = SchedStrategy.FCFS;
            end
            self.schedStrategyPar = zeros(1,numOfClasses);

            % Preserve a user-configured finite place capacity (setClassCapacity
            % is called before model.link, which invokes init): a bounded place is
            % what makes an open Source->Place->Transition net a finite ergodic
            % CTMC. Only (re)initialize entries the user has not set.
            if isempty(self.classCap)
                self.classCap = Inf(1,numOfClasses);
            elseif length(self.classCap) < numOfClasses
                tmp = Inf(1,numOfClasses);
                tmp(1:length(self.classCap)) = self.classCap;
                self.classCap = tmp;
            end
            self.cap = Inf;
            self.schedStrategies = ones(1, numOfClasses);
            for r=1:numOfClasses
                classes = self.model.getClasses();
                self.dropRule(classes{r}.index) = DropStrategy.WAITQ;
            end
        end

        function self = setClassCapacity(self, class, capacity)
            % SELF = SETCLASSCAPACITY(CLASS, CAPACITY)

            % Resolve a JobClass object to its column index (mirrors Source.m);
            % indexing classCap by the object itself silently fails to store the
            % capacity, so the Place stays unbounded and setClassCapacity is a
            % no-op for the CTMC marking bound.
            if isa(class, 'JobClass')
                r = class.index;
            else
                r = class;
            end
            self.classCap(r) = capacity;
        end

        function self = setSchedStrategies(self, class, strategy)
            % SELF = SETSCHEDSTRATEGIES(CLASS, STRATEGY)

            self.schedStrategies(class) = strategy;
        end

        function setService(self, class, distribution)
            % SETSERVICE(CLASS, DISTRIBUTION)
            % Assigns a service process to a token color, turning this ordinary
            % place into a queueing place. The embedded queue serves tokens under
            % the place's scheduling strategy; the first such call installs the
            % concrete server section for that strategy.

            if isa(distribution,'Workflow')
                distribution = distribution.toPH();
            end
            if distribution.isImmediate()
                distribution = Immediate.getInstance();
            end
            if isempty(self.obj)
                if ~self.queueing
                    self.installQueueServer();
                    self.queueing = true;
                end
                c = class.index;
                self.server.serviceProcess{1, c}{1} = class;
                self.server.serviceProcess{1, c}{2} = ServiceStrategy.LI;
                self.server.serviceProcess{1, c}{3} = distribution;
                self.serviceProcess{c} = distribution;
                if length(self.departureDiscipline) < c || self.departureDiscipline(c)==0
                    self.departureDiscipline(c) = DepartureDiscipline.NORMAL;
                end
                self.model.setInitialized(false);
            else
                self.obj.setService(class.obj, distribution.obj);
                self.serviceProcess{class.index} = distribution;
                self.queueing = true;
            end
        end

        function installQueueServer(self)
            % INSTALLQUEUESERVER()
            % Installs the concrete server section matching the place's scheduling
            % strategy, replacing the default ServiceTunnel used by ordinary places.

            classes = self.model.getClasses();
            switch SchedStrategy.toId(self.schedStrategy)
                case SchedStrategy.INF
                    self.server = InfiniteServer(classes);
                    self.numberOfServers = Inf;
                case {SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.LPS}
                    self.server = SharedServer(classes);
                    if isinf(self.numberOfServers); self.numberOfServers = 1; end
                otherwise
                    self.server = Server(classes);
                    if isinf(self.numberOfServers); self.numberOfServers = 1; end
            end
            % A queueing place is a station with finite/infinite capacity like a Queue;
            % initialize cap/classCap (left empty by the ordinary-place constructor).
            nclasses = length(classes);
            if isempty(self.cap); self.cap = Inf; end
            if isempty(self.classCap); self.classCap = Inf(1,nclasses); end
            if isempty(self.dropRule)
                for r=1:nclasses
                    self.dropRule(classes{r}.index) = DropStrategy.WAITQ;
                end
            end
        end

        function bool = isQueueing(self)
            % BOOL = ISQUEUEING()
            bool = self.queueing;
        end

        function distribution = getService(self, class)
            % DISTRIBUTION = GETSERVICE(CLASS)
            distribution = self.serviceProcess{class.index};
        end

        function setDepartureDiscipline(self, class, discipline)
            % SETDEPARTUREDISCIPLINE(CLASS, DISCIPLINE)
            self.departureDiscipline(class.index) = discipline;
            if ~isempty(self.obj)
                self.obj.setDepartureDiscipline(class.obj, jline.lang.constant.DepartureDiscipline.fromID(discipline));
            end
        end

        function self = setMarking(self, state)
            % SELF = SETMARKING(STATE)
            % Alias for setState using Petri-net terminology: sets the initial
            % token marking (number of tokens per class) of this place.

            self.setState(state);
        end

    end
end

