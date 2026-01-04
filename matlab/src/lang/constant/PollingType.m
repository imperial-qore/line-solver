classdef (Sealed) PollingType
    % Enumeration of polling service types
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        GATED = 0;
        EXHAUSTIVE = 1;
        KLIMITED = 2;
    end
    
    methods (Static)        
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            switch type
                case PollingType.GATED
                    text = 'Gated';
                case PollingType.EXHAUSTIVE
                    text = 'Exhaustive';
                case PollingType.KLIMITED
                    text = 'K-Limited';
            end
        end

        function id = toId(type)
            % ID = TOID(TYPE)
            id = type;
        end
    end
end
