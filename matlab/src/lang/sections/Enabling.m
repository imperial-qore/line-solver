classdef Enabling < InputSection
    % An input section for jobs to wait for service.
    %
    % Copyright (c) 2012-2024, Imperial College London
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
            self.inputJobClasses = {};
            initQueueJobClasses(self, classes);
        end
    end
    
    methods (Access = 'private')
        function initQueueJobClasses(self, customerClasses)
            % INITQUEUEJOBCLASSES(CUSTOMERCLASSES)
            
            for i = 1 : length(customerClasses)
                self.inputJobClasses{i} = {customerClasses{i}, SchedStrategy.FCFS, DropStrategy.WaitingQueue};
            end
        end
    end
    
end

