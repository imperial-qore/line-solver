classdef (Sealed) ReplacementStrategy
    % Enumeration of cache replacement strategies
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        RR = 0;
        FIFO = 1;
        SFIFO = 2; % strict fifo
        LRU = 3;
        %HLRU = 4; % also called k-LRU (Martina et al, Gast & van Houdt; not supported at present)
    end

    methods (Static)

        function text = toString(type)
            % TEXT = TOSTRING(ID)
            text = ReplacementStrategy.toText(type);
        end

        function text = toText(type)
            % TEXT = TOTEXT(ID)
            switch type
                case ReplacementStrategy.RR
                    text = 'rr';
                case ReplacementStrategy.FIFO
                    text = 'fifo';
                case ReplacementStrategy.SFIFO
                    text = 'strict-fifo';
                case ReplacementStrategy.LRU
                    text = 'lru';
            end
        end

        function text = toFeature(type)
            % TEXT = TOFEATURE(TYPE)

            switch type
                case ReplacementStrategy.RR
                    text = 'ReplacementStrategy_RR';
                case ReplacementStrategy.FIFO
                    text = 'ReplacementStrategy_FIFO';
                case ReplacementStrategy.SFIFO
                    text = 'ReplacementStrategy_SFIFO';
                case ReplacementStrategy.LRU
                    text = 'ReplacementStrategy_LRU';
                    %case ReplacementStrategy.HLRU
                    %    text = 'ReplacementStrategy_HLRU';
            end
        end
    end
end
