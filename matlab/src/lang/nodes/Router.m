classdef Router < StatefulNode
    % Router Probabilistic job routing node for queueing networks
    %
    % Router is a stateful node that routes jobs towards other nodes based on
    % probabilistic routing strategies. Unlike service stations, jobs do not
    % queue inside a Router node but are immediately routed to downstream nodes
    % according to configured routing probabilities or strategies.
    %
    % @brief Probabilistic routing node that directs jobs to downstream destinations
    %
    % Key characteristics:
    % - Jobs flow through without queueing or service delays
    % - Supports probabilistic routing matrices
    % - Maintains state for routing decisions
    % - Can implement complex routing strategies (random, round-robin, etc.)
    % - Essential for modeling load balancing and traffic distribution
    %
    % Router nodes are commonly used in:
    % - Load balancing scenarios
    % - Multi-path routing
    % - Traffic distribution models
    % - Complex network topologies with branching
    %
    % Example:
    % @code
    % router = Router(model, 'LoadBalancer');
    % model.addNode(router);
    % % Set routing probabilities in routing matrix
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        cap;
        numberOfServers;
        schedPolicy;
        schedStrategy;
    end

    methods
        % Constructor
        function self = Router(model, name)
            % ROUTER Create a Router node instance
            %
            % @brief Creates a Router node for probabilistic job routing
            % @param model Network model to add the router to
            % @param name String identifier for the router node
            % @return self Router instance configured for the given model
            %
            % Note: This is a node and not a Station because jobs cannot queue
            % inside it - they flow through immediately to downstream nodes.

            self@StatefulNode(name);
            self.setModel(model);
            if isa(model,'Network')
                % Handle the new delegation pattern
                if model.isMatlabNative()
                    classes = model.getClasses();
                    self.schedPolicy = SchedStrategyType.NP;
                    self.schedStrategy = SchedStrategy.FCFS;
                    self.input = Buffer(classes);
                    self.output = Dispatcher(classes);
                    self.cap = Inf;
                    self.schedPolicy = SchedStrategyType.NP;
                    self.server = ServiceTunnel();
                    self.numberOfServers = 1;
                    if ~model.addNode(self) % if not a replacement
                        for r=1:length(classes)
                            self.setRouting(RoutingStrategy.RAND,classes{r});
                        end
                    end
                else
                    % Java implementation through delegation
                    self.obj = jline.lang.nodes.Router(model.implementation.obj, name);
                    self.index = model.implementation.obj.getNodeIndex(self.obj);
                end
            elseif model.isMatlabNative()
                classes = model.getClasses();
                self.schedPolicy = SchedStrategyType.NP;
                self.schedStrategy = SchedStrategy.FCFS;
                self.input = Buffer(classes);
                self.output = Dispatcher(classes);
                self.cap = Inf;
                self.schedPolicy = SchedStrategyType.NP;
                self.server = ServiceTunnel();
                self.numberOfServers = 1;
                if ~self.model.addNode(self) % if not a replacement
                    for r=1:length(classes)
                        self.setRouting(RoutingStrategy.RAND,classes{r});
                    end
                end
            elseif model.isJavaNative()
                self.obj=jline.lang.nodes.Router(model.obj, name);
            end
        end

        function setProbRouting(self, class, destination, probability)
            % SETPROBROUTING(CLASS, DESTINATION, PROBABILITY)

            setRouting(self, class, RoutingStrategy.PROB, destination, probability);
        end

        function setService(self, class, distribution)
            % SETSERVICE(CLASS, DISTRIBUTION)

            self.server.serviceProcess{1, class.index}{2} = ServiceStrategy.LI;
            self.server.serviceProcess{1, class.index}{3} = distribution;
        end

    end

end
