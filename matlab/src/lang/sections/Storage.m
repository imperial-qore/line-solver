classdef Storage < InputSection
    % An input section for jobs to wait for service.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = Storage(classes)
            % SELF = STORAGE(CLASSES)
            
            self@InputSection('Storage');
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

