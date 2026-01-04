classdef (Sealed) DropStrategy
    % Enumeration of drop policies in stations.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        WAITQ = -1;
        DROP = 1;
        BAS = 2;
        BBS = 3;
        RSRD = 4;
        RETRIAL = 5;            % Job moves to orbit and retries after delay (unlimited attempts)
        RETRIAL_WITH_LIMIT = 6; % Job retries up to max attempts, then drops
    end

    methods (Static)
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            switch type
                case DropStrategy.WAITQ
                    text = 'waiting queue';
                case DropStrategy.DROP
                    text = 'drop';
                case DropStrategy.BAS
                    text = 'BAS blocking';
                case DropStrategy.BBS
                    text = 'BBS blocking';
                case DropStrategy.RSRD
                    text = 'RSRD blocking';
                case DropStrategy.RETRIAL
                    text = 'retrial';
                case DropStrategy.RETRIAL_WITH_LIMIT
                    text = 'retrial with limit';
                otherwise
                    line_error(mfilename, 'Unrecognized drop strategy type.');
            end
        end
    end
end


