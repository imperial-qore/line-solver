classdef ServiceSection < Section
    % An abstract class for the service section of a node.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        numberOfServers;
        serviceProcess;
    end
    
    methods(Hidden)
        %Constructor
        function self = ServiceSection(className)
            % SELF = SERVICESECTION(CLASSNAME)
            
            self@Section(className);
        end
    end
end
