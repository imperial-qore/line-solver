classdef Dispatcher < OutputSection
    % Dispatcher Output section for routing jobs to downstream destinations
    %
    % Dispatcher implements an output section that routes completed jobs to
    % their next destinations according to specified routing strategies. It
    % manages class-dependent routing decisions and coordinates job flow
    % through the queueing network topology.
    %
    % @brief Output section for class-dependent job routing and dispatching
    %
    % Key characteristics:
    % - Class-dependent routing strategies
    % - Destination selection and job dispatching
    % - Support for various routing policies
    % - Integration with network topology
    % - Probabilistic and deterministic routing
    %
    % Dispatcher features:
    % - Routing strategy management per class
    % - Destination probability specification
    % - Random, round-robin, and custom routing
    % - Class-based routing differentiation
    % - Network flow coordination
    %
    % Dispatcher is used in:
    % - Queue output sections for job routing
    % - Router nodes for traffic distribution
    % - Load balancing implementations
    % - Multi-path network routing
    % - Class-based service differentiation
    %
    % Example:
    % @code
    % classes = {OpenClass(model, 'WebReq'), OpenClass(model, 'DBReq')};
    % dispatcher = Dispatcher(classes);
    % dispatcher.setStrategy(classes{1}, RoutingStrategy.PROB);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = Dispatcher(customerClasses)
            % SELF = DISPATCHER(CUSTOMERCLASSES)
            
            self@OutputSection('Dispatcher');
            self.outputStrategy = {};
            initDispatcherJobClasses(self, customerClasses);
        end
    end
    
    methods
        function initDispatcherJobClasses(self, customerClasses)
            % INITDISPATCHERJOBCLASSES(CUSTOMERCLASSES)
            
            for r = 1 : length(customerClasses)
                self.outputStrategy{r} = {customerClasses{r}.name, RoutingStrategy.toText(RoutingStrategy.DISABLED)};
            end
        end
        
        function setStrategy(self, customerClasses, strategy)
            % SETSTRATEGY(CUSTOMERCLASSES, STRATEGY)
            
            if length(strategy)==1 %#ok<ISCL>
                self.outputStrategy{customerClasses{r}.index} = {customerClasses{r}.name, strategy};
            else
                for r = 1 : length(customerClasses)
                    self.outputStrategy{customerClasses{r}.index} = {customerClasses{r}.name, strategy{r}};
                end
            end
        end
    end
    
end
