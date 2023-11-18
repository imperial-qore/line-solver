classdef Sink < Node
    % A node to let jobs in open classes depart the model
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    properties
        schedStrategy;
    end

    methods
        %Constructor
        function self = Sink(model, name)
            % SELF = SINK(MODEL, NAME)

            self@Node(name);            
            if isa(model,'Network')
                if model ~= 0
                    self.input = '';
                    self.output = '';
                    self.server = Section('JobSink');
                    self.setModel(model);
                    self.model.addNode(self);
                    self.schedStrategy = SchedStrategy.EXT;
                    if length(model.classes)>1 % if Sink created after some closed classes
                        for r=1:model.getNumberOfClasses
                            self.setRouting(model.classes{r},RoutingStrategy.DISABLED);
                        end                        
                    end
                end
            elseif isa(model,'JNetwork')
                self.setModel(model);
                self.obj=jline.lang.nodes.Sink(model.obj, name);
            end

        end

        function sections = getSections(self)
            % SECTIONS = GETSECTIONS()

            sections = {'', self.server, ''};
        end

        function self = setScheduling(self, class, strategy)
            %noop
        end

    end

    methods(Access = protected)
        % Override copyElement method:
        function clone = copyElement(self)
            % CLONE = COPYELEMENT()

            % Make a shallow copy of all properties
            clone = copyElement@Copyable(self);
            % Make a deep copy of each object
            clone.input = self.input;
            clone.server = self.server.copy;
            clone.output = self.output;
        end

    end

end
