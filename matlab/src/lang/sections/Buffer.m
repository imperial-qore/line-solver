classdef Buffer < InputSection
    % An input section for jobs to wait for service.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        size;
    end
    
    methods
        %Constructor
        function self = Buffer(classes)
            % SELF = BUFFER(CLASSES)
            
            self@InputSection('Buffer');
            self.size = -1;
            self.schedPolicy = SchedStrategyType.NP;
        end
        
        function self = setSize(self, size)
            self.size = size;
        end
    end
    
end

