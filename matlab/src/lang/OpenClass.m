classdef OpenClass < JobClass
    % A class of jobs that arrive from the external world
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    methods

        %Constructor
        function self = OpenClass(model, name, prio)
            % SELF = OPENCLASS(MODEL, NAME, PRIO)

            self@JobClass('open', name);
            if nargin<3
                prio = 0;
            end
            if isa(model,'Network')
                if isempty(model.getSource)
                    line_error(mfilename,'The model requires a Source prior to instantiating open classes.');
                end
                if isempty(model.getSink)
                    line_error(mfilename,'The model requires a Sink prior to instantiating open classes.');
                end
                if nargin >= 3 && isa(prio,'Source')
                    % user error, source passed as ref station
                    prio = 0;
                end
                self.type = JobClassType.OPEN;
                self.priority = 0;
                if nargin>=3 %exist('prio','var')
                    self.priority = prio;
                end
                addJobClass(model, self);
                setReferenceStation(self, model.getSource());

                % set default scheduling for this class at all nodes
                for i=1:length(model.nodes)
                    if ~isempty(model.nodes{i}) && ~isa(model.nodes{i},'Source') && ~isa(model.nodes{i},'Sink') && ~isa(model.nodes{i},'Join')
                        %&& ~isa(model.nodes{i},'Fork')
                        model.nodes{i}.setScheduling(self, SchedStrategy.FCFS);
                    end
                    if isa(model.nodes{i},'Join')
                        model.nodes{i}.setStrategy(self,JoinStrategy.STD);
                        model.nodes{i}.setRequired(self,-1);
                    end
                    if ~isempty(model.nodes{i})
                        % this should be disabled, but it causes a problem with JMT
                        %if ~isa(model.nodes{i},'Source')
                        %    model.nodes{i}.setRouting(self, RoutingStrategy.DISABLED);
                        %else
                            model.nodes{i}.setRouting(self, RoutingStrategy.RAND);
                        %end
                    end
                end
            elseif isa(model,'JNetwork')
                self.index = model.getNumberOfClasses + 1;
                self.obj = jline.lang.OpenClass(model.obj, name, prio);
            end
        end

        function setReferenceStation(class, source)
            % SETREFERENCESTATION(CLASS, SOURCE)

            if ~isa(source,'Source')
                line_error(mfilename,'The reference station for an open class must be a Source.');
            end
            setReferenceStation@JobClass(class, source);
        end
    end

end

