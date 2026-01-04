classdef Forker < OutputSection
    % An output section forking jobs into sibling tasks
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        tasksPerLink;
    end
    
    methods
        %Constructor
        function self = Forker(customerClasses)
            % SELF = FORKER(CUSTOMERCLASSES)
            
            self@OutputSection('Forker');
            self.tasksPerLink = 1.0;
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
    
    methods(Access = protected)
        % Override copyElement method:
        function clone = copyElement(self)
            % CLONE = COPYELEMENT()
            
            % Make a shallow copy of all properties
            clone = copyElement@Copyable(self);
            % Make a deep copy of each object
            for i = 1 : length(self.outputStrategy)
                if ishandle(clone.outputStrategy{i}{1})
                    % this is a problem if one modifies the classes in the
                    % model because the one below is not an handle so it
                    % will not be modified
                    clone.outputStrategy{i}{1} = self.outputStrategy{i}{1}.copy;
                end
            end
        end
    end
end
