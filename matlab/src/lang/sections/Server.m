classdef Server < ServiceSection
    % Server Basic non-preemptive service section for job processing
    %
    % Server implements a basic service section that processes jobs according
    % to non-preemptive disciplines. It manages service processes for different
    % job classes and coordinates with scheduling sections to provide complete
    % queueing station functionality.
    %
    % @brief Non-preemptive service section with class-dependent service processes
    %
    % Key characteristics:
    % - Non-preemptive job service
    % - Class-dependent service time distributions
    % - Single or multiple server support
    % - Load-independent service strategies
    % - Integration with scheduling disciplines
    %
    % Server features:
    % - Service process management per class
    % - Service time distribution configuration
    % - Server capacity specification
    % - Load-independent operation
    % - Deep copy support for model cloning
    %
    % Server is used in:
    % - Basic FCFS and LCFS queues
    % - Priority queue service sections
    % - Multi-server queue implementations
    % - Load-independent service modeling
    % - Class-differentiated service systems
    %
    % Example:
    % @code
    % classes = {OpenClass(model, 'Class1'), ClosedClass(model, 'Class2', 5)};
    % server = Server(classes);
    % server.setServiceProcess(1, 1, Exp(2.0)); % Class1 service
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = Server(classes)
            % SERVER Create a basic server section instance
            %
            % @brief Creates a Server section for non-preemptive job processing
            % @param classes Cell array of JobClass objects to be served
            % @return self Server instance configured for the specified classes
            
            self@ServiceSection('Server');
            self.numberOfServers = 1;
            self.serviceProcess = {};
            initServers(self, classes);
        end
    end
    
    methods (Access = 'private')
        function initServers(self, classes)
            % INITSERVERS(CLASSES)
            
            for i = 1 : length(classes)
                self.serviceProcess{1, i} = {classes{i}.name, ServiceStrategy.LI, Exp(0.0)};
            end
        end
    end
    
    methods(Access = protected)
        % Override copyElement method:
        function clone = copyElement(self)
            % CLONE = COPYELEMENT()
            
            % Make a shallow copy of all properties
            clone = copyElement@Copyable(self);
            % Make a deep copy of each object
            for i = 1 : length(self.serviceProcess)
                if ishandle(clone.serviceProcess{1,i}{3})
                    % this is a problem if one modifies the classes in the
                    % model because the one below is not an handle so it
                    % will not be modified
                    clone.serviceProcess{1,i}{3} = self.serviceProcess{1,i}{3}.copy;
                end
            end
        end
    end
    
end

