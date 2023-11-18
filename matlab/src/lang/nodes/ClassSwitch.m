classdef ClassSwitch < Node
    % A node to change the class of visiting jobs
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    properties (Hidden)
        autoAdded;
    end

    properties
        cap;
        schedPolicy;
        schedStrategy;
    end

    methods
        %Constructor
        function self = ClassSwitch(model, name, csMatrix)
            % SELF = CLASSSWITCH(MODEL, NAME, CSMATRIX)
            classes = model.classes;
            if nargin < 3
                csMatrix = eye(length(classes));
            end
            self@Node(name);

            self.autoAdded = false;
            self.input = Buffer(classes);
            self.output = Dispatcher(classes);
            self.cap = Inf;
            self.schedPolicy = SchedStrategyType.NP;
            self.schedStrategy = SchedStrategy.FCFS;
            self.server = StatelessClassSwitcher(classes, csMatrix);
            self.setModel(model);
            self.model.addNode(self);
        end

        function C = initClassSwitchMatrix(self)
            % C = INITCLASSSWITCHMATRIX()

            K = self.model.getNumberOfClasses;
            C = zeros(K,K);
        end

        function setClassSwitchingMatrix(self, csMatrix)
            self.server.updateClasses(self.model.classes);
            self.server.updateClassSwitch(csMatrix);
        end

        function setProbRouting(self, class, destination, probability)
            % SETPROBROUTING(CLASS, DESTINATION, PROBABILITY)

            setRouting(self, class, RoutingStrategy.PROB, destination, probability);
        end

        function sections = getSections(self)
            % SECTIONS = GETSECTIONS()

            sections = {self.input, self.server, self.output};
        end

        function summary(self)
            % SUMMARY()

            line_printf('\nNode: <strong>%s</strong>',self.getName);
            for r=1:length(self.output.outputStrategy)
                %line_printf('Routing %s: %s',self.model.classes{r}.name,self.output.outputStrategy{r}{2});
                for s=1:length(self.output.outputStrategy)
                    if self.server.csMatrix(r,s)>0
                        line_printf('Routing %s->%s: %g',self.model.classes{r}.name,self.model.classes{s}.name,self.server.csMatrix(r,s));
                    end
                end
            end
            %            self.input.summary;
            %            self.server.summary;
            %            self.output.summary;
        end
    end

end
