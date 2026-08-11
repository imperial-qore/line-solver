classdef (Sealed) SignalType
    % Enumeration of signal types for signal classes.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        NEGATIVE = 1;     % Negative customer signal (removes jobs)
        REPLY = 0;        % Reply signal (triggers unblocking)
        CATASTROPHE = 2;  % Catastrophe signal (removes ALL jobs)
    end

    methods (Static)
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            switch type
                case SignalType.NEGATIVE
                    text = 'negative';
                case SignalType.REPLY
                    text = 'reply';
                case SignalType.CATASTROPHE
                    text = 'catastrophe';
                otherwise
                    line_error(mfilename, 'Unrecognized signal type.');
            end
        end

        function type = fromText(text)
            % TYPE = FROMTEXT(TEXT)
            if iscell(text)
                text = text{:};
            end
            switch lower(text)
                case 'negative'
                    type = SignalType.NEGATIVE;
                case 'reply'
                    type = SignalType.REPLY;
                case 'catastrophe'
                    type = SignalType.CATASTROPHE;
                otherwise
                    line_error(mfilename, 'Unrecognized signal type text.');
            end
        end
    end
end
