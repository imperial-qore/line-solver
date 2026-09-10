classdef (Sealed) DepartureDiscipline
    % Departure disciplines for the depository of a queueing place (QPN semantics).
    %
    % A queueing place serves tokens in its embedded queue and, on service
    % completion, moves them to a depository from which they become available to
    % the output transitions. The departure discipline governs the order in which
    % depository tokens become available.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        NORMAL = 0;   % Tokens available immediately upon service completion (standard QPN)
        FIFO = 1;     % Tokens available in their order of arrival to the depository
    end

    methods (Static)
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            switch type
                case DepartureDiscipline.NORMAL
                    text = 'normal';
                case DepartureDiscipline.FIFO
                    text = 'fifo';
                otherwise
                    line_error(mfilename, 'Unrecognized departure discipline type.');
            end
        end
    end
end
