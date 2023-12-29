classdef Task < LayeredNetworkElement
    % A software server in a LayeredNetwork.
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    properties
        parent;
        multiplicity;       %int
        replication;       %int
        scheduling;         %string
        thinkTime;
        thinkTimeMean;      %double
        thinkTimeSCV;       %double
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
            self.scheduling = scheduling;
            self.setThinkTime(thinkTime);
            self.parent = [];
            model.tasks{end+1} = self;
            switch scheduling
                case 'ref'
                    model.reftasks{end+1} = self;
            end            
            self.model = model;
        end
        
        function self = setReplication(self, replication)
            self.replication = replication;
        end
                
        function self = on(self, parent)
            % self = ON(self, PARENT)
            self.parent = parent;
            parent.addTask(self);
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
        
    end
    
end