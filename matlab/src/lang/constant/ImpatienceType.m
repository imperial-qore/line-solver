classdef (Sealed) ImpatienceType
    % Enumeration of customer impatience types.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        RENEGING = 1;  % Customer abandons after joining the queue (timer-based)
        BALKING = 2;   % Customer refuses to join based on queue state
        RETRIAL = 3;   % Customer moves to orbit and retries after delay
    end

    methods (Static)
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            %
            % Convert impatience type constant to text description
            switch type
                case ImpatienceType.RENEGING
                    text = 'reneging';
                case ImpatienceType.BALKING
                    text = 'balking';
                case ImpatienceType.RETRIAL
                    text = 'retrial';
                otherwise
                    line_error(mfilename, 'Unrecognized impatience type.');
            end
        end

        function id = toId(type)
            % ID = TOID(TYPE)
            %
            % Convert impatience type constant to numeric ID
            id = type;
        end

        function type = fromId(id)
            % TYPE = FROMID(ID)
            %
            % Convert numeric ID to impatience type constant
            switch id
                case 1
                    type = ImpatienceType.RENEGING;
                case 2
                    type = ImpatienceType.BALKING;
                case 3
                    type = ImpatienceType.RETRIAL;
                otherwise
                    line_error(mfilename, sprintf('Unrecognized impatience type ID: %d', id));
            end
        end
    end
end
