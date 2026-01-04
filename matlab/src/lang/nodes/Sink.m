classdef Sink < Node
    % External job departure node for open queueing networks
    %
    % Represents network exit point where jobs are removed from open classes.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        schedStrategy;
    end

    methods
        %Constructor
        function self = Sink(model, name)
            % SINK Create an external departure sink node
            %
            % @brief Creates a Sink node for external job removal
            % @param model Network model to add the sink node to
            % @param name String identifier for the sink node
            % @return self Sink instance ready for job absorption

            self@Node(name);            
            if model.isMatlabNative()
                if model ~= 0
                    self.input = '';
                    self.output = '';
                    self.server = Section('JobSink');
                    self.setModel(model);
                    self.model.addNode(self);
                    self.schedStrategy = SchedStrategy.EXT;
                    if length(model.getClasses())>1 % if Sink created after some closed classes
                        for r=1:model.getNumberOfClasses
                            classes = model.getClasses();
                    self.setRouting(classes{r},RoutingStrategy.DISABLED);
                        end                        
                    end
                end
            elseif model.isJavaNative()
                self.setModel(model);
                self.obj=jline.lang.nodes.Sink(model.obj, name);
            end

        end

        function sections = getSections(self)
            % SECTIONS = GETSECTIONS()

            sections = {'', self.server, ''};
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
