classdef ClassSwitch < Node
    % A node to change the class of visiting jobs
    %
    % Copyright (c) 2012-2026, Imperial College London
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
            
            self@Node(name);
            if isa(model,'Network')
                % Handle the new delegation pattern
                if model.isMatlabNative()
                    classes = model.getClasses();
                    if nargin < 3
                        csMatrix = eye(length(classes));
                    end
                    
                    self.autoAdded = false;
                    self.input = Buffer(classes);
                    self.output = Dispatcher(classes);
                    self.cap = Inf;
                    self.schedPolicy = SchedStrategyType.NP;
                    self.schedStrategy = SchedStrategy.FCFS;
                    self.server = StatelessClassSwitcher(classes, csMatrix);
                    self.setModel(model);
                    model.addNode(self);
                else
                    % Java implementation through delegation
                    self.setModel(model);
                    if nargin < 3
                        self.obj = jline.lang.nodes.ClassSwitch(model.implementation.obj, name);
                    else
                        self.obj = jline.lang.nodes.ClassSwitch(model.implementation.obj, name, csMatrix);
                    end
                    self.index = model.implementation.obj.getNodeIndex(self.obj);
                end
            elseif isa(model,'MNetwork') || (isa(model,'Network') && model.isMatlabNative())
                classes = model.getClasses();
                if nargin < 3
                    csMatrix = eye(length(classes));
                end
                
                self.autoAdded = false;
                self.input = Buffer(classes);
                self.output = Dispatcher(classes);
                self.cap = Inf;
                self.schedPolicy = SchedStrategyType.NP;
                self.schedStrategy = SchedStrategy.FCFS;
                self.server = StatelessClassSwitcher(classes, csMatrix);
                self.setModel(model);
                self.model.addNode(self);
            elseif isa(model,'JNetwork') || (isa(model,'Network') && model.isJavaNative())
                self.setModel(model);
                if nargin < 3
                    self.obj = jline.lang.nodes.ClassSwitch(model.obj, name);
                else
                    self.obj = jline.lang.nodes.ClassSwitch(model.obj, name, csMatrix);
                end
                self.index = model.obj.getNodeIndex(self.obj);
            end
        end

        function C = initClassSwitchMatrix(self)
            % C = INITCLASSSWITCHMATRIX()

            K = self.model.getNumberOfClasses;
            C = zeros(K,K);
        end

        function setClassSwitchingMatrix(self, csMatrix)
            self.server.updateClasses(self.model.getClasses());
            self.server.updateClassSwitch(csMatrix);
        end

        function setProbRouting(self, class, destination, probability)
            % SETPROBROUTING(CLASS, DESTINATION, PROBABILITY)

            setRouting(self, class, RoutingStrategy.PROB, destination, probability);
        end

        function summary(self)
            % SUMMARY()

            line_printf('\nNode: <strong>%s</strong>',self.getName);
            for r=1:length(self.output.outputStrategy)
                %classes = self.model.getClasses();
            %line_printf('Routing %s: %s',classes{r}.name,self.output.outputStrategy{r}{2});
                for s=1:length(self.output.outputStrategy)
                    if self.server.csMatrix(r,s)>0
                        classes = self.model.getClasses();
            line_printf('Routing %s->%s: %g',classes{r}.name,classes{s}.name,self.server.csMatrix(r,s));
                    end
                end
            end
            %            self.input.summary;
            %            self.server.summary;
            %            self.output.summary;
        end
    end

end
