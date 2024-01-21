classdef MarkovModulated < PointProcess
    % An abstract class for Markov-modulated processes
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    methods (Hidden)
        %Constructor
        function self = MarkovModulated(name, numParam)
            % SELF = MARKOVMODULATED(NAME, NUMPARAM)
            
            self@PointProcess(name, numParam);
        end
    end
    
    methods
        function X = sample(self, n)
            % X = SAMPLE(N)
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end
        
        function P = getEmbDTMC(self)
            % P = GETEMBDTMC()
            %
            % Get DTMC embedded at event arrival times
            
            P = map_embedded(self.getRepres);
        end

        function pie = getEmbProb(self)
            % PIE = GETEMBPROB()
            %
            % Solve DTMC embedded embedded at event arrival times
            
            pie = map_pie(self.getRepres);
        end
    end
    
    methods %(Abstract) % implemented with errors for Octave compatibility
        function phases = getNumberOfPhases(self)
            % PHASES = GETNUMBEROFPHASES()
            
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end
        function MAP = getRepres(self)
            % MAP = GETREPRESENTATION()
            
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end
    end
    
    methods (Static)
        function cx = fit(MEAN, SCV)
            % CX = FIT(MEAN, SCV)
            
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end
    end
    
end

