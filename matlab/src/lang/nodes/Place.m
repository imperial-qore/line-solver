classdef Place < Station
    % 
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties               
        schedStrategies;
        schedStrategy;
        schedStrategyPar;
    end

    methods
        function self = Place(model,name)
            % PLACE(MODEL, NAME)
            
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
            elseif model.isJavaNative()
                self.setModel(model);
                self.obj = jline.lang.nodes.Place(model.obj, name);
                self.index = model.obj.getNodeIndex(self.obj);
            end
        end

        function init(self)            
            numOfClasses = length(self.model.getClasses());
            self.schedStrategy = SchedStrategy.FCFS;
            self.schedStrategyPar = zeros(1,numOfClasses);

            self.classCap = Inf(1,numOfClasses);
            self.cap = Inf;
            self.schedStrategies = ones(1, numOfClasses);
            for r=1:numOfClasses
                classes = self.model.getClasses();
                self.dropRule(classes{r}.index) = DropStrategy.WAITQ;
            end
        end

        function self = setClassCapacity(self, class, capacity)
            % SELF = SETCLASSCAPACITY(CLASS, CAPACITY)
            
            self.classCap(class) = capacity;
        end

        function self = setSchedStrategies(self, class, strategy)
            % SELF = SETSCHEDSTRATEGIES(CLASS, STRATEGY)
            
            self.schedStrategies(class) = strategy;
        end

    end
end

