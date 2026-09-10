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
        HLRU = 4;  % hierarchical/k-LRU: h lists, LRU discipline, promote i->i+1 on hit
        CLIMB = 5; % move-up-one-position on hit (transposition rule)
        QLRU = 6;  % q-LRU: LRU discipline with probabilistic admission q on a miss
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
                case ReplacementStrategy.HLRU
                    text = 'hlru';
                case ReplacementStrategy.CLIMB
                    text = 'climb';
                case ReplacementStrategy.QLRU
                    text = 'qlru';
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
                case ReplacementStrategy.HLRU
                    text = 'ReplacementStrategy_HLRU';
                case ReplacementStrategy.CLIMB
                    text = 'ReplacementStrategy_CLIMB';
                case ReplacementStrategy.QLRU
                    text = 'ReplacementStrategy_QLRU';
            end
        end
    end
end
