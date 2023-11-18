classdef ClosedClass < JobClass
    % A class of jobs that perpetually cycle inside the model.
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    properties
        population;
    end

    methods

        %Constructor
        function self = ClosedClass(model, name, njobs, refstat, prio)
            % SELF = CLOSEDCLASS(MODEL, NAME, NJOBS, REFSTAT, PRIO)

            self@JobClass('closed', name);

            %global GlobalConstants.CoarseTol
            
            if nargin<5
                prio = 0;
            end

            if isa(model,'Network')
                self.type = JobClassType.CLOSED;
                self.population = njobs;
                if abs(njobs-round(njobs))>GlobalConstants.CoarseTol
                    line_warning(mfilename,sprintf('The number of jobs in class %s should be an integer, some solvers might fail.\n', name));
                end
                self.priority = 0;
                if nargin>=5 %exist('prio','var')
                    self.priority = prio;
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
                    %if ~isempty(model.nodes{i})  %&& ~isa(model.nodes{i},'Join')&& ~isa(model.nodes{i},'Cache') %&& ~isa(model.nodes{i},'Source') && ~isa(model.nodes{i},'Sink')
                    %    model.nodes{i}.setScheduling(self,SchedStrategy.FCFS);
                    %end
                    %if ~isempty(model.nodes{i})
                    %                    && (isa(model.nodes{i},'Queue') || isa(model.nodes{i},'Router'))
                    %model.nodes{i}.setRouting(self, RoutingStrategy.DISABLED);
                    model.nodes{i}.setRouting(self, RoutingStrategy.RAND);
                    if isa(model.nodes{i},'Join')
                        model.nodes{i}.setStrategy(self, JoinStrategy.STD);
                        model.nodes{i}.setRequired(self, -1);
                    end
                    %end
                end
            elseif isa(model,'JNetwork')
                self.index = model.getNumberOfClasses + 1;
                self.obj = jline.lang.ClosedClass(model.obj, name, njobs, refstat.obj, prio);
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

