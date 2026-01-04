classdef Task < LayeredNetworkElement
    % A software server in a LayeredNetwork.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        parent;
        multiplicity;       %int
        replication;       %int
        scheduling;         %string
        priority = 0;       %int, priority level (0 = default/no priority)
        fanInSource = '';   %string, source task for fan-in
        fanInValue = 0;     %int, fan-in value (load distribution count)
        thinkTime;
        thinkTimeMean;      %double
        thinkTimeSCV;       %double
        setupTime;
        setupTimeMean;      %double
        setupTimeSCV;       %double
        delayOffTime;
        delayOffTimeMean;   %double
        delayOffTimeSCV;    %double
        entries = [];
        activities = [];
        precedences = [];
        replyEntry;
    end
    
    
    methods
        %public methods, including constructor
        
        %constructor
        function self = Task(model, name, multiplicity, scheduling, thinkTime)
            % self = TASK(MODEL, NAME, MULTIPLICITY, SCHEDULING, THINKTIME)
            name = char(name);
            if nargin<2%~exist('name','var')
                line_error(mfilename,'Constructor requires to specify at least a name.');
            end
            self@LayeredNetworkElement(name);
            
            if nargin<3%~exist('multiplicity','var')
                multiplicity = 1;
            end
            if nargin<4%~exist('scheduling','var')
                scheduling = SchedStrategy.INF;
            end
            if nargin<5%~exist('thinkTime','var')
                thinkTime = GlobalConstants.FineTol;
            end
            self.replication = 1;
            self.multiplicity = multiplicity;
            switch scheduling
                case SchedStrategy.INF
                if isfinite(multiplicity)
                    line_warning(mfilename,'Finite multiplicity is not allowed with INF scheduling. Setting it to INF.\n');
                    self.multiplicity = Inf;
                end
            end
            self.scheduling = SchedStrategy.toText(scheduling);
            self.setThinkTime(thinkTime);
            self.setSetupTime(Immediate());
            self.setDelayOffTime(Immediate());
            self.parent = [];

            if isa(model,'LayeredNetwork')
                model.tasks{end+1} = self;
                switch scheduling
                    case 'ref'
                        model.reftasks{end+1} = self;
                end            
                self.model = model;
            elseif isa(model,'JLayeredNetwork')
                % JLayeredNetwork support would go here if it exists
                self.obj = jline.lang.layered.Task(model.obj, name, multiplicity, scheduling, thinkTime);
                model.addTask(self);
            end
        end
        
        function self = setReplication(self, replication)
            self.replication = replication;
        end
                
        function self = on(self, parent)
            % self = ON(self, PARENT)
            if ~isa(parent,'Host') && ~isa(parent,'Processor')
                line_error(mfilename,'Invalid .on() argument: expected a host processor.')
            end            
            if isempty(self.parent)
                self.parent = parent;
                parent.addTask(self);
            else
                line_error(mfilename,'Parent processor already defined.')
            end
        end
        
        function self = setAsReferenceTask(self)
            % self = SETASREFERENCETASK(self)
            
            self.scheduling = SchedStrategy.REF;
        end
        
        function self = setThinkTime(self, thinkTime)
            % self = SETTHINKTIME(self, THINKTIME)

            if isnumeric(thinkTime)
                if thinkTime <= GlobalConstants.FineTol
                    self.thinkTime = Immediate.getInstance();
                    self.thinkTimeMean = GlobalConstants.FineTol;
                    self.thinkTimeSCV = GlobalConstants.FineTol;
                else
                    self.thinkTime = Exp(1/thinkTime);
                    self.thinkTimeMean = thinkTime;
                    self.thinkTimeSCV = 1.0;
                end
            elseif isa(thinkTime,'Distribution')
                self.thinkTime = thinkTime;
                self.thinkTimeMean = thinkTime.getMean();
                self.thinkTimeSCV = thinkTime.getSCV();
            end
        end

        function self = setSetupTime(self, setupTime)
            % self = SETSETUPTIME(self, SETUPTIME)
            % Set the setup time (cold start time) for the task.

            if isnumeric(setupTime)
                if setupTime <= GlobalConstants.FineTol
                    self.setupTime = Immediate.getInstance();
                    self.setupTimeMean = GlobalConstants.FineTol;
                    self.setupTimeSCV = GlobalConstants.FineTol;
                else
                    self.setupTime = Exp(1/setupTime);
                    self.setupTimeMean = setupTime;
                    self.setupTimeSCV = 1.0;
                end
            elseif isa(setupTime,'Distribution')
                self.setupTime = setupTime;
                self.setupTimeMean = setupTime.getMean();
                self.setupTimeSCV = setupTime.getSCV();
            end
        end

        function self = setDelayOffTime(self, delayOffTime)
            % self = SETDELAYOFFTIME(self, DELAYOFFTIME)
            % Set the delay-off time (teardown time) for the task.

            if isnumeric(delayOffTime)
                if delayOffTime <= GlobalConstants.FineTol
                    self.delayOffTime = Immediate.getInstance();
                    self.delayOffTimeMean = GlobalConstants.FineTol;
                    self.delayOffTimeSCV = GlobalConstants.FineTol;
                else
                    self.delayOffTime = Exp(1/delayOffTime);
                    self.delayOffTimeMean = delayOffTime;
                    self.delayOffTimeSCV = 1.0;
                end
            elseif isa(delayOffTime,'Distribution')
                self.delayOffTime = delayOffTime;
                self.delayOffTimeMean = delayOffTime.getMean();
                self.delayOffTimeSCV = delayOffTime.getSCV();
            end
        end

        function result = hasSetupDelayoff(self)
            % result = HASSETUPDELAYOFF(self)
            % Check if this task has setup/delayoff configured (i.e., non-trivial values).

            hasSetup = ~isempty(self.setupTime) && ~isa(self.setupTime, 'Immediate') ...
                       && self.setupTimeMean > GlobalConstants.FineTol;
            hasDelayoff = ~isempty(self.delayOffTime) && ~isa(self.delayOffTime, 'Immediate') ...
                          && self.delayOffTimeMean > GlobalConstants.FineTol;
            result = hasSetup || hasDelayoff;
        end

        %addEntry
        function self = addEntry(self, newEntry)
            % self = ADDENTRY(self, NEWENTRY)
            
            self.entries = [self.entries; newEntry];
        end
        
        %addActivity
        function self = addActivity(self, newAct)
            % self = ADDACTIVITY(self, NEWACT)
            
            newAct.setParent(self.name);
            self.activities = [self.activities; newAct];
        end
        
        %setActivity
        function self = setActivity(self, newAct, index)
            % self = SETACTIVITY(self, NEWACT, INDEX)
            
            self.activities(index,1) = newAct;
        end
        
        %removeActivity
        function self = removeActivity(self, index)
            % self = REMOVEACTIVITY(self, INDEX)
            
            idxToKeep = [1:index-1,index+1:length(self.activities)];
            self.activities = self.activities(idxToKeep);
            self.actNames = self.actNames(idxToKeep);
        end
        
        %addPrecedence
        function self = addPrecedence(self, newPrec)
            % self = ADDPRECEDENCE(self, NEWPREC)
            
            if iscell(newPrec)
                for m=1:length(newPrec)
                    self.precedences = [self.precedences; newPrec{m}];
                end
            else
                self.precedences = [self.precedences; newPrec];
            end
        end
        
        %setReplyEntry
        function self = setReplyEntry(self, newReplyEntry)
            % self = SETREPLYENTRY(self, NEWREPLYENTRY)
            
            self.replyEntry = newReplyEntry;
        end
        
        function meanHostDemand = getMeanHostDemand(self, entryName)
            % MEANHOSTDEMAND = GETMEANHOSTDEMAND(self, ENTRYNAME)

            % determines the demand posed by the entry entryName
            % the demand is located in the activity of the corresponding entry

            meanHostDemand = -1;
            for j = 1:length(self.entries)
                if strcmp(self.entries(j).name, entryName)
                    meanHostDemand = self.entries(j).activities(1).hostDemandMean;
                    break;
                end
            end
        end

        % Getter methods for API consistency with Java/Python
        function val = getMultiplicity(obj)
            % GETMULTIPLICITY Get the multiplicity (number of task instances)
            val = obj.multiplicity;
        end

        function val = getReplication(obj)
            % GETREPLICATION Get the replication factor
            val = obj.replication;
        end

        function val = getScheduling(obj)
            % GETSCHEDULING Get the scheduling strategy
            val = obj.scheduling;
        end

        function val = getThinkTimeMean(obj)
            % GETTHINKTIMEMEAN Get the mean think time
            val = obj.thinkTimeMean;
        end

        function val = getThinkTimeSCV(obj)
            % GETTHINKTIMESCV Get the SCV of think time
            val = obj.thinkTimeSCV;
        end

        function val = getParent(obj)
            % GETPARENT Get the parent host/processor
            val = obj.parent;
        end

        function val = getPrecedences(obj)
            % GETPRECEDENCES Get the list of activity precedences
            val = obj.precedences;
        end

        function val = getSetupTimeMean(obj)
            % GETSETUPTIMEMEAN Get the mean setup time
            val = obj.setupTimeMean;
        end

        function val = getDelayOffTimeMean(obj)
            % GETDELAYOFFTIMEMEAN Get the mean delay-off time
            val = obj.delayOffTimeMean;
        end

        % setPriority
        function self = setPriority(self, priority)
            % self = SETPRIORITY(self, PRIORITY)
            self.priority = priority;
        end

        % setFanIn
        function self = setFanIn(self, source, value)
            % self = SETFANIN(self, SOURCE, VALUE)
            self.fanInSource = source;
            self.fanInValue = value;
        end

    end

end