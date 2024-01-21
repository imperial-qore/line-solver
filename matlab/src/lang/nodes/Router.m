classdef Router < StatefulNode
    % A node to route jobs towards other nodes
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties
        cap;
        numberOfServers;
        schedPolicy;
        schedStrategy;
    end
    
    methods
        %Constructor
        % This is a node and not a Station because jobs cannot station
        % inside it
        function self = Router(model, name)
            % SELF = ROUTER(MODEL, NAME)
            
            self@StatefulNode(name);
            
            classes = model.classes;
            self.schedPolicy = SchedStrategyType.NP;
            self.schedStrategy = SchedStrategy.FCFS;
            self.input = Buffer(classes);
            self.output = Dispatcher(classes);
            self.cap = Inf;
            self.schedPolicy = SchedStrategyType.NP;
            self.server = ServiceTunnel();
            self.numberOfServers = 1;
            self.setModel(model);
            if ~self.model.addNode(self) % if not a replacement
                for r=1:length(classes)
                    self.setRouting(RoutingStrategy.RAND,classes{r});
                end
            end
        end
        
        function setProbRouting(self, class, destination, probability)
            % SETPROBROUTING(CLASS, DESTINATION, PROBABILITY)
            
            setRouting(self, class, RoutingStrategy.PROB, destination, probability);
        end
        
        function setScheduling(self, class, strategy)
            % SETSCHEDULING(CLASS, STRATEGY)
            
            self.input.inputJobClasses{1, class.index}{2} = strategy;
        end
        
        function setService(self, class, distribution)
            % SETSERVICE(CLASS, DISTRIBUTION)
            
            self.server.serviceProcess{1, class.index}{2} = ServiceStrategy.LI;
            self.server.serviceProcess{1, class.index}{3} = distribution;
        end
        
        function sections = getSections(self)
            % SECTIONS = GETSECTIONS()
            
            sections = {self.input, self.server, self.output};
        end
    end
    
end
