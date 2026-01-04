classdef Join < Station
    % Task synchronization node for fork-join models
    %
    % Combines parallel sibling tasks created by Fork nodes, waiting until all arrive.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        joinStrategy;
        joinOf;
    end
    
    methods
        %Constructor
        function self = Join(model, name, fork)
            % JOIN Create a Join node instance
            %
            % @brief Creates a Join node for synchronizing parallel tasks
            % @param model Network model to add the join node to
            % @param name String identifier for the join node  
            % @param fork Optional Fork node this join synchronizes with
            % @return self Join instance configured for task synchronization
            %
            % The constructor creates a Join node with appropriate joiners,
            % dispatchers, and service tunnels. If a fork parameter is provided,
            % the join is associated with that specific fork node.
            
            self@Station(name);
            if model.isMatlabNative()
                if(model ~= 0)
                    classes = model.getClasses();
                    self.input = Joiner(classes);
                    self.output = Dispatcher(classes);
                    self.server = ServiceTunnel();
                    self.dropRule = [];
                    self.numberOfServers = Inf;
                    self.setModel(model);
                    self.joinOf = fork;
                    model.addNode(self);
                end
            elseif model.isJavaNative()
                self.setModel(model);
                if nargin >= 3 && ~isempty(fork)
                    self.obj = jline.lang.nodes.Join(model.obj, name, fork.obj);
                else
                    self.obj = jline.lang.nodes.Join(model.obj, name);
                end
                self.index = model.obj.getNodeIndex(self.obj);
            end
            %             if ~exist('joinstrategy','var')
            %                 joinstrategy = JoinStrategy.STD;
            %             end
            %             setStrategy(joinstrategy);
        end
    end
    
    methods
        function self = setStrategy(self, class, strategy)
            % SELF = SETSTRATEGY(CLASS, STRATEGY)
            
            self.input.setStrategy(class,strategy);
        end
        
        function self = setRequired(self, class, njobs)
            % SELF = SETREQUIRED(CLASS, NJOBS)
            
            self.input.setRequired(class,njobs);
        end
        
        function self = setProbRouting(self, class, destination, probability)
            % SELF = SETPROBROUTING(CLASS, DESTINATION, PROBABILITY)
            
            setRouting(self, class, RoutingStrategy.PROB, destination, probability);
        end

        function summary(self)
            % SUMMARY()

            line_printf('\nNode: <strong>%s</strong>',self.getName);            
            for r=1:length(self.output.outputStrategy)
                classes = self.model.getClasses();
                line_printf('Routing %s: %s\n', classes{r}.name, self.output.outputStrategy{r}{2});
            end
        end
        
    end
    
end
