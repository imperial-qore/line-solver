classdef Firing < OutputSection
    % An output section for Transitions
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = Firing(customerClasses)
            % SELF = FORKER(CUSTOMERCLASSES)
            
            self@OutputSection('Firing');
            initDispatcherJobClasses(self, customerClasses);
        end
    end
    
    methods (Access = 'private')
        function initDispatcherJobClasses(self, customerClasses)
            % INITDISPATCHERJOBCLASSES(CUSTOMERCLASSES)
            
            for i = 1 : length(customerClasses)
                self.outputStrategy{i} = {customerClasses{i}.name, RoutingStrategy.RAND};
            end
        end
    end

end
