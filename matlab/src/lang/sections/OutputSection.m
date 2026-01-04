classdef OutputSection < Section
    % An abstract class for the output section of a node.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        outputStrategy;
    end
    
    methods(Hidden)
        %Constructor
        function self = OutputSection(className)
            % SELF = OUTPUTSECTION(CLASSNAME)
            
            self@Section(className);
        end
    end
end
