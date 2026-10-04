classdef Section < NetworkElement
    % An abstract class for sections that define the behavior of a Node
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        className;
    end
    
    methods(Hidden)
        %Constructor
        function self = Section(className)
            % SELF = SECTION(CLASSNAME)
            
            self@NetworkElement('Section');
            self.className = className;
        end
    end
    
    methods
        function summary(self)
            % SUMMARY()
            
            %            line_printf('\n%s',class(self));
        end
    end
end
