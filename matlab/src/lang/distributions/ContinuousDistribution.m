classdef ContinuousDistribution < Distribution
    % An abstract class for continuous distributions
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    methods (Hidden)
        function self = ContinuousDistribution(name, numParam, support)
            % SELF = CONTINUOUSDISTRIB(NAME, NUMPARAM, SUPPORT)
            
            % Construct a continuous distribution from name, number of
            % parameters, and range
            self@Distribution(name,numParam,support);
        end
    end
    
end
