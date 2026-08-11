classdef Linkage < OutputSection
    % An output section dispatching jobs to a destination
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = Linkage(customerClasses)
            % SELF = DISPATCHER(CUSTOMERCLASSES)
            
            self@OutputSection('Linkage');
            self.outputStrategy = {};
            initDispatcherJobClasses(self, customerClasses);
        end
    end
    
    methods
        function initDispatcherJobClasses(self, customerClasses)
            % INITDISPATCHERJOBCLASSES(CUSTOMERCLASSES)
            
            for r = 1 : length(customerClasses)
                self.outputStrategy{r} = {customerClasses{r}.name, RoutingStrategy.RAND};
            end
        end
        
        function setStrategy(self, customerClasses, strategy)
            % SETSTRATEGY(CUSTOMERCLASSES, STRATEGY)
            
            if isscalar(strategy)
                self.outputStrategy{customerClasses{r}.index} = {customerClasses{r}.name, strategy};
            else
                for r = 1 : length(customerClasses)
                    self.outputStrategy{customerClasses{r}.index} = {customerClasses{r}.name, strategy{r}};
                end
            end
        end
    end
    
end
