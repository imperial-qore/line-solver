classdef InputSection < Section
    % An abstract class for the input section of a node.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        schedPolicy;        
    end
    
    methods(Hidden)
        %Constructor
        function self = InputSection(className)
            % SELF = INPUTSECTION(CLASSNAME)
            
            self@Section(className);
        end
    end
end
