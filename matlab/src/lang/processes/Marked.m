classdef Marked < MarkovModulated
    % An abstract class for Markov-modulated marked processes
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods (Hidden)
        %Constructor
        function self = MarkedMarkovModulated(name, numParam)
            % SELF = MARKOVMODULATED(NAME, NUMPARAM)

            self@PointProcess(name, numParam);
        end
    end

    methods
        function X = sample(self, n)
            % X = SAMPLE(N)
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end

        function P = getEmbedded(self)
            % P = GETEMBEDDED()
            %
            % Get DTMC embedded at event arrival times

            P = map_embedded(self.getProcess);
        end

        function pie = getEmbeddedProb(self)
            % PIE = GETEMBEDDEDRPOB()
            %
            % Solve DTMC embedded embedded at event arrival times

            pie = map_pie(self.getProcess);
        end
    end

    methods (Abstract)
        phases = getNumberOfPhases(self);
    end

end

