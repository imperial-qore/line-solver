classdef ActivityPrecedenceType
    % Activity Precedence types.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        PRE_SEQ    = 1;
        PRE_AND    = 2;
        PRE_OR     = 3;
        POST_SEQ   = 11;
        POST_AND   = 12;
        POST_OR    = 13;
        POST_LOOP  = 14;
        POST_CACHE = 15;
    end

    methods (Static)
        function txt = toText(precedence)
            % TXT = TOTEXT(PRECEDENCE)
            % Convert numeric ID to string label
            switch precedence
                case ActivityPrecedenceType.PRE_SEQ
                    txt = 'pre';
                case ActivityPrecedenceType.PRE_AND
                    txt = 'pre-AND';
                case ActivityPrecedenceType.PRE_OR
                    txt = 'pre-OR';
                case ActivityPrecedenceType.POST_SEQ
                    txt = 'post';
                case ActivityPrecedenceType.POST_AND
                    txt = 'post-AND';
                case ActivityPrecedenceType.POST_OR
                    txt = 'post-OR';
                case ActivityPrecedenceType.POST_LOOP
                    txt = 'post-LOOP';
                case ActivityPrecedenceType.POST_CACHE
                    txt = 'post-CACHE';
                otherwise
                    line_error(mfilename, 'Unrecognized precedence type.');
            end
        end
    end
end
