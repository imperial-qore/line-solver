classdef Host  < LayeredNetworkElement
    % A hardware server in a LayeredNetwork.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        multiplicity;       %int
        replication;        %int
        scheduling;         %char: ps, fcfs, inf, ref
        quantum;            %double
        speedFactor;        %double
        tasks = [];         %list of tasks
        ID;                 %int
    end
    
    methods
        %public methods, including constructor
        
        %constructor
        function self = Host(model, name, multiplicity, scheduling, quantum, speedFactor)
            % self = HOST(MODEL, NAME, MULTIPLICITY, SCHEDULING, QUANTUM, SPEEDFACTOR)
            
            if nargin<2 %~exist('name','var')
                line_error(mfilename,'Constructor requires to specify at least a name.');
            end
            self@LayeredNetworkElement(name);
            
            if nargin<3 %~exist('multiplicity','var')
                multiplicity = 1;
            end
            if nargin<4 %~exist('scheduling','var')
                scheduling = SchedStrategy.PS;
            end
            if nargin<5 %~exist('quantum','var')
                quantum = 0.001;
            end
            if nargin<6 %~exist('speedFactor','var')
                speedFactor = 1;
            end
            self.replication = 1;            
            self.multiplicity = multiplicity;
            self.scheduling = SchedStrategy.toText(scheduling);
            self.quantum = quantum;
            self.speedFactor = speedFactor;
            
            if isa(model,'LayeredNetwork')
                model.hosts{end+1} = self;
                self.model = model;
            elseif isa(model,'JLayeredNetwork')
                % JLayeredNetwork support would go here if it exists
                self.obj = jline.lang.layered.Host(model.obj, name, multiplicity, scheduling, quantum, speedFactor);
                model.addHost(self);
            end
        end
        
        function self=setReplication(self, replication)
            self.replication = replication;
        end
        
        
        %addTask
        function self = addTask(self, newTask)
            % self = ADDTASK(self, NEWTASK)
            self.tasks = [self.tasks; newTask];
            newTask.parent = self;
        end

        % Getter methods for API consistency with Java/Python
        function val = getMultiplicity(obj)
            % GETMULTIPLICITY Get the multiplicity (number of processor instances)
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

        function val = getQuantum(obj)
            % GETQUANTUM Get the time quantum for scheduling
            val = obj.quantum;
        end

        function val = getSpeedFactor(obj)
            % GETSPEEDFACTOR Get the speed factor
            val = obj.speedFactor;
        end

    end

end
