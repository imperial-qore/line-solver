classdef ActivityPrecedenceType
    % Activity Precedence types.
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    properties (Constant)
        PRE_SEQ  = 'pre';
        PRE_AND = 'pre-AND';
        PRE_OR = 'pre-OR';
        POST_SEQ = 'post';
        POST_AND = 'post-AND';
        POST_OR = 'post-OR';
        POST_LOOP = 'post-LOOP';
        POST_CACHE = 'post-CACHE';

        ID_PRE_SEQ  = 1;
        ID_PRE_AND = 2;
        ID_PRE_OR = 3;
        ID_POST_SEQ = 11;
        ID_POST_AND = 12;
        ID_POST_OR = 13;
        ID_POST_LOOP = 14;
        ID_POST_CACHE = 15;
    end

    methods (Static)
        function typeId = toId(precedence)
            % TYPEID = TOID(PRECEDENCE)
            % Classifies the activity precedence
            switch precedence
                case ActivityPrecedenceType.PRE_SEQ
                    typeId = ActivityPrecedenceType.ID_PRE_SEQ;
                case ActivityPrecedenceType.PRE_AND
                    typeId = ActivityPrecedenceType.ID_PRE_AND;
                case ActivityPrecedenceType.PRE_OR
                    typeId = ActivityPrecedenceType.ID_PRE_OR;
                case ActivityPrecedenceType.POST_SEQ
                    typeId = ActivityPrecedenceType.ID_POST_SEQ;
                case ActivityPrecedenceType.POST_AND
                    typeId = ActivityPrecedenceType.ID_POST_AND;
                case ActivityPrecedenceType.POST_OR
                    typeId = ActivityPrecedenceType.ID_POST_OR;
                case ActivityPrecedenceType.POST_LOOP
                    typeId = ActivityPrecedenceType.ID_POST_LOOP;
                case ActivityPrecedenceType.POST_CACHE
                    typeId = ActivityPrecedenceType.ID_POST_CACHE;
                otherwise
                    line_error(mfilename,'Unrecognized precedence type.');
            end
        end
        function feat = toFeature(precedence)
            % TYPEID = TOFEATURE(PRECEDENCE)
            % Return the precedence label
            switch precedence
                case ActivityPrecedenceType.ID_PRE_SEQ
                    feat = 'ActivityPrecedenceType_PRE_SEQ';
                case ActivityPrecedenceType.ID_PRE_AND
                    feat = 'ActivityPrecedenceType_PRE_AND';
                case ActivityPrecedenceType.ID_PRE_OR
                    feat = 'ActivityPrecedenceType_PRE_OR';
                case ActivityPrecedenceType.ID_POST_SEQ
                    feat = 'ActivityPrecedenceType_POST_SEQ';
                case ActivityPrecedenceType.ID_POST_AND
                    feat = 'ActivityPrecedenceType_POST_AND';
                case ActivityPrecedenceType.ID_POST_OR
                    feat = 'ActivityPrecedenceType_POST_OR';
                case ActivityPrecedenceType.ID_POST_LOOP
                    feat = 'ActivityPrecedenceType_POST_LOOP';
                case ActivityPrecedenceType.ID_POST_CACHE
                    feat = 'ActivityPrecedenceType_POST_CACHE';
                otherwise
                    line_error(mfilename,'Unrecognized precedence type.');
            end
        end
    end
end
