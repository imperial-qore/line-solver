classdef StatefulClassSwitcher < ClassSwitcher
    % An abstract class for a state-dependent class switch
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = StatefulClassSwitcher(classes, name)
            % SELF = STATEFULCLASSSWITCHER(CLASSES, NAME)
            
            self@ClassSwitcher(classes, name);
            self.csFun = @(r, s, state, statep) StatefulClassSwitcher.classHolderFun(r, s, state, statep); % do nothing by default
        end
    end
    
    methods (Static)
        function prob = classHolderFun(r, s, state, statep)
            % PROB = CLASSHOLDERFUN(R, S, STATE, STATEP)
            
            % Class holder: keep the class (identity switch), independent of state
            if r == s
                prob = 1;
            else
                prob = 0;
            end
        end
    end
end
