classdef (Sealed) TimingStrategy
    % Enumeration of timing policies in Petri nets transitions.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        TIMED = 0;
        IMMEDIATE = 1;
    end
    
    methods (Static)
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            switch type
                case TimingStrategy.TIMED
                    text = 'timed';
                case TimingStrategy.IMMEDIATE
                    text = 'immediate';
                otherwise
                    line_error(mfilename, 'Unrecognized timing strategy type.');
            end
        end
    end
end
