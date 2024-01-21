classdef Buffer < InputSection
    % An input section for jobs to wait for service.
    %
    % Copyright (c) 2012-2024, Imperial College London
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
            self.inputJobClasses = {};
            initQueueJobClasses(self, classes);
        end
        
        function self = setSize(self, size)
            self.size = size;
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
    
    methods(Access = protected)
        % Override copyElement method:
        function clone = copyElement(self)
            % CLONE = COPYELEMENT()
            
            % Make a shallow copy of all properties
            clone = copyElement@Copyable(self);
            % Make a deep copy of each object
            for i = 1 : length(self.inputJobClasses)
                if ishandle(clone.inputJobClasses{i}{1})
                    % this is a problem if one modifies the classes in the
                    % model because the one below is not an handle so it
                    % will not be modified
                    clone.inputJobClasses{i}{1} = self.inputJobClasses{i}{1}.copy;
                end
            end
        end
    end
    
end

