classdef (Sealed) PollingType
    % Enumeration of polling service types
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        ID_GATED = 0;
        ID_EXHAUSTIVE = 1;
        ID_KLIMITED = 2;
        
        GATED = 'gated'; % gated service type
        EXHAUSTIVE = 'exhaustive'; % exhaustive service type
        KLIMITED = 'klimited'; % Ki-limited service type at buffer i
    end
    
    methods (Access = private)
        %private so that it cannot be instatiated.
        function out = PollingType
            % OUT = POLLINGTYPE
            
        end
    end
    
    methods (Static)
        function typeId = getTypeId(strategyType)
            % TYPEID = GETTYPEID(STRATEGYTYPE)
            % Classifies the polling service type
            switch strategyType
                case GATED
                    typeId = PollingType.ID_GATED;
                case EXHAUSTIVE
                    typeId = PollingType.ID_EXHAUSTIVE;
                case KLIMITED
                    typeId = PollingType.ID_KLIMITED;
                otherwise
                    line_error(mfilename,'Unrecognized polling type.');
            end
        end
        
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
            switch type
                case PollingType.GATED
                    id = PollingType.ID_GATED;
                case PollingType.EXHAUSTIVE
                    id = PollingType.ID_EXHAUSTIVE;
                case PollingType.KLIMITED
                    id = PollingType.ID_KLIMITED;
            end
        end
    end
end
