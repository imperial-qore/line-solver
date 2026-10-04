classdef (Sealed) RemovalPolicy
    % Enumeration of removal policies for negative signals in G-networks.
    %
    % When a negative signal (or catastrophe) arrives at a queue and removes
    % positive customers, the removal policy determines which customers are
    % selected for removal.
    %
    % Policies:
    % - RANDOM: Uniform random selection from all jobs at the station
    % - FCFS: Remove oldest jobs first (first-come-first-served order)
    % - LCFS: Remove newest jobs first (last-come-first-served order)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        RANDOM = 0;  % Uniform random selection
        FCFS = 1;    % Remove oldest jobs first
        LCFS = 2;    % Remove newest jobs first
    end

    methods (Static)
        function text = toText(policy)
            % TEXT = TOTEXT(POLICY) Convert policy to text representation
            switch policy
                case RemovalPolicy.RANDOM
                    text = 'random';
                case RemovalPolicy.FCFS
                    text = 'fcfs';
                case RemovalPolicy.LCFS
                    text = 'lcfs';
                otherwise
                    line_error(mfilename, 'Unrecognized removal policy.');
            end
        end

        function policy = fromText(text)
            % POLICY = FROMTEXT(TEXT) Convert text to policy constant
            if iscell(text)
                text = text{:};
            end
            switch lower(text)
                case 'random'
                    policy = RemovalPolicy.RANDOM;
                case 'fcfs'
                    policy = RemovalPolicy.FCFS;
                case 'lcfs'
                    policy = RemovalPolicy.LCFS;
                otherwise
                    line_error(mfilename, 'Unrecognized removal policy text.');
            end
        end
    end
end
