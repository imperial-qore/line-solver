classdef (Sealed) BalkingStrategy
    % Enumeration of balking strategies.
    %
    % BalkingStrategy defines how customers decide whether to balk (refuse
    % to join a queue):
    %   QUEUE_LENGTH - Balk based on current queue length ranges
    %   EXPECTED_WAIT - Balk based on expected waiting time
    %   COMBINED - Balk if either condition is met (OR logic)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        QUEUE_LENGTH = 1;   % Balk based on queue length ranges with probability
        EXPECTED_WAIT = 2;  % Balk based on expected waiting time
        COMBINED = 3;       % Both conditions (OR logic)
    end

    methods (Static)
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            %
            % Convert balking strategy constant to text description
            switch type
                case BalkingStrategy.QUEUE_LENGTH
                    text = 'queue length';
                case BalkingStrategy.EXPECTED_WAIT
                    text = 'expected wait';
                case BalkingStrategy.COMBINED
                    text = 'combined';
                otherwise
                    line_error(mfilename, 'Unrecognized balking strategy type.');
            end
        end

        function id = toId(type)
            % ID = TOID(TYPE)
            %
            % Convert balking strategy constant to numeric ID
            id = type;
        end

        function type = fromId(id)
            % TYPE = FROMID(ID)
            %
            % Convert numeric ID to balking strategy constant
            switch id
                case 1
                    type = BalkingStrategy.QUEUE_LENGTH;
                case 2
                    type = BalkingStrategy.EXPECTED_WAIT;
                case 3
                    type = BalkingStrategy.COMBINED;
                otherwise
                    line_error(mfilename, sprintf('Unrecognized balking strategy ID: %d', id));
            end
        end
    end
end
