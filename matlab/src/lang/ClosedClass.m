classdef ClosedClass < JobClass
    % ClosedClass Job class with fixed population circulating in the network
    %
    % ClosedClass represents a job class with a fixed number of jobs that
    % perpetually circulate within the network. Jobs never leave the system
    % and the total population remains constant. This is essential for modeling
    % systems with limited resources or finite user populations.
    %
    % @brief Job class with fixed population circulating within the network
    %
    % Key characteristics:
    % - Fixed finite population of jobs
    % - Jobs never enter or leave the network
    % - Constant total network population
    % - Reference station for performance metrics
    % - Priority-based service differentiation
    % - Think time modeling for user behavior
    %
    % Closed class features:
    % - Fixed population constraint enforcement
    % - Reference station designation for metrics
    % - Think time specification for user delays
    % - Priority assignment for service
    % - Circulation-based performance analysis
    % - Saturation and throughput modeling
    %
    % ClosedClass is used for:
    % - Terminal-based computer systems
    % - Time-sharing system modeling
    % - Manufacturing systems with fixed workpieces
    % - Batch processing systems
    % - Systems with resource constraints
    %
    % Example:
    % @code
    % model = Network('ClosedSystem');
    % cpu = Queue(model, 'CPU', SchedStrategy.PS);
    % disk = Queue(model, 'Disk', SchedStrategy.FCFS);
    % think = Delay(model, 'ThinkTime');
    % users = ClosedClass(model, 'Users', 20, think, 1); % 20 users
    % think.setService(users, Exp(0.1)); % Think time
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        population;
    end

    methods

        %Constructor
        function self = ClosedClass(model, name, njobs, refstat, prio, deadline)
            % CLOSEDCLASS Create a closed job class instance
            %
            % @brief Creates a ClosedClass with fixed population circulating in network
            % @param model Network model to add the closed class to
            % @param name String identifier for the job class
            % @param njobs Number of jobs in the class (fixed population)
            % @param refstat Reference station for performance measurement
            % @param prio Optional priority level (default: 0)
            % @param deadline Optional relative deadline from arrival (default: Inf, no deadline)
            % @return self ClosedClass instance with specified population

            self@JobClass(JobClassType.CLOSED, name);

            %global GlobalConstants.CoarseTol

            if nargin<5
                prio = 0;
            end
            if nargin<6
                deadline = Inf;
            end

            if model.isMatlabNative()
                self.type = JobClassType.CLOSED;
                self.population = njobs;
                if abs(njobs-round(njobs))>GlobalConstants.CoarseTol
                    line_warning(mfilename,sprintf('The number of jobs in class %s should be an integer, some solvers might fail.\n', name));
                end
                self.priority = 0;
                self.deadline = Inf;
                if nargin>=5 %exist('prio','var')
                    self.priority = prio;
                end
                if nargin>=6
                    self.deadline = deadline;
                end
                model.addJobClass(self);
                if ~isa(refstat, 'Station')
                    if isa(refstat, 'Node')
                        line_error(mfilename,sprintf('The reference station of class %s needs to be a station, not a node.\n', name));
                    else
                        line_error(mfilename,sprintf('The parameter for the reference station of class %s is not a valid object.\n', name));
                    end
                end
                setReferenceStation(self, refstat);

                % set default scheduling for this class at all nodes
                for i=1:length(model.nodes)
                    model.nodes{i}.setRouting(self, RoutingStrategy.RAND);
                    if isa(model.nodes{i},'Join')
                        model.nodes{i}.setStrategy(self, JoinStrategy.STD);
                        model.nodes{i}.setRequired(self, -1);
                    end
                end
            elseif model.isJavaNative()
                self.index = model.obj.getNumberOfClasses + 1;
                if ~isa(self,'SelfLoopingClass')
                    self.obj = jline.lang.ClosedClass(model.obj, name, njobs, refstat.obj, prio, deadline);
                end
                setReferenceStation(self, refstat);
            end
        end

        function setReferenceStation(class, source)
            % SETREFERENCESTATION(CLASS, SOURCE)
            setReferenceStation@JobClass(class, source);
        end

        function summary(self)
            % SUMMARY()
            line_printf('Class (%s): <strong>%s</strong> [%d jobs]',self.type,self.getName,self.population);
        end
    end

end

