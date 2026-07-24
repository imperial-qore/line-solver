classdef DiscreteDistribution < Distribution
    % An abstract class for continuous distributions
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods (Hidden)
        function self = DiscreteDistribution(name, numParam, support)
            % SELF = DISCRETEDISTRIB(NAME, NUMPARAM, SUPPORT)
            
            % Construct a continuous distribution from name, number of
            % parameters, and range
            self@Distribution(name,numParam,support);
        end
    end
    
end
