classdef RandomSource < InputSection
    % An input section generating arrivals from the external world at random
    % times
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        sourceClasses;
    end
    
    methods
        %Constructor
        function self = RandomSource(classes)
            % SELF = RANDOMSOURCE(CLASSES)
            
            self@InputSection('RandomSource');
            for i = 1 : length(classes)
                self.sourceClasses{1,i} = {[],ServiceStrategy.LI,Disabled.getInstance()};
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
            for i = 1 : length(self.sourceClasses)
                clone.sourceClasses{1,i}{3} = self.sourceClasses{1,i}{3}.copy;
            end
        end
    end
    
end

