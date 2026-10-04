classdef (Sealed) PollingType
    % Enumeration of polling service types
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        GATED = 0;
        EXHAUSTIVE = 1;
        KLIMITED = 2;
        DECREMENTING = 3;
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
                case PollingType.DECREMENTING
                    text = 'Decrementing';
            end
        end

        function id = toId(type)
            % ID = TOID(TYPE)
            id = type;
        end

        function name = toName(type)
            % NAME = TONAME(TYPE) - Canonical name used on the JSON wire.
            % The numeric ids agree with the Java enum but not with the Python
            % one, which assigns them via auto(), so names are used to
            % interchange polling types across codebases.
            switch PollingType.toId(type)
                case PollingType.GATED
                    name = 'GATED';
                case PollingType.EXHAUSTIVE
                    name = 'EXHAUSTIVE';
                case PollingType.KLIMITED
                    name = 'KLIMITED';
                case PollingType.DECREMENTING
                    name = 'DECREMENTING';
                otherwise
                    line_error(mfilename, sprintf('Unable to return a PollingType name for value: %d.\n', PollingType.toId(type)));
            end
        end

        function id = fromName(name)
            % ID = FROMNAME(NAME) - Inverse of toName.
            switch upper(name)
                case 'GATED'
                    id = PollingType.GATED;
                case 'EXHAUSTIVE'
                    id = PollingType.EXHAUSTIVE;
                case 'KLIMITED'
                    id = PollingType.KLIMITED;
                case 'DECREMENTING'
                    id = PollingType.DECREMENTING;
                otherwise
                    line_error(mfilename, sprintf('Unable to return a PollingType for name: %s.\n', name));
            end
        end
    end
end
