classdef (Sealed) TimingStrategy
    % Enumeration of timing polcies in petri nets transitions.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        TIMED = 'timed';
        IMMEDIATE = 'immediate';
        
        ID_TIMED = 0;
        ID_IMMEDIATE = 1;
    end
    
    methods (Static)
        
        function id = toId(type)
            % ID = TOOD(TYPE)
            
            switch type
                case TimingStrategy.TIMED
                    id = TimingStrategy.ID_TIMED;
                case TimingStrategy.IMMEDIATE
                    id = TimingStrategy.ID_IMMEDIATE;
            end
            
        end
        
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            
            switch type
                case TimingStrategy.Timed
                    text = 'Timed Transition';
                case TimingStrategy.Immediate
                    text = 'Immediate Transition';
            end
        end
    end
end

