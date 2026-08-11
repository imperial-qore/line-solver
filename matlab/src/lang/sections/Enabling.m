classdef Enabling < InputSection
    % An input section for jobs to wait for service.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        size;
    end
    
    methods
        %Constructor
        function self = Enabling(classes)
            % SELF = ENABLING(CLASSES)
            
            self@InputSection('Enabling');
            self.size = -1;
            self.schedPolicy = SchedStrategyType.PR;
        end
    end
       
end

