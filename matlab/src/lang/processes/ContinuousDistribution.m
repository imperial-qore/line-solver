classdef ContinuousDistribution < Distribution
    % ContinuousDistribution Abstract base class for continuous distributions
    %
    % ContinuousDistribution provides the common interface and functionality
    % for all continuous-valued statistical distributions. It extends the base
    % Distribution class with methods specific to continuous random variables
    % such as Laplace-Stieltjes transform evaluation.
    %
    % @brief Abstract base class for continuous-valued distributions
    %
    % Key characteristics:
    % - Abstract interface for continuous distributions
    % - Support over continuous intervals (often (0,∞) or ℝ)
    % - Provides Laplace-Stieltjes transform interface
    % - Foundation for exponential, gamma, uniform, etc.
    % - Integrates with queueing theory analysis methods
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods (Hidden)
        function self = ContinuousDistribution(name, numParam, support)
            % SELF = CONTINUOUSDISTRIB(NAME, NUMPARAM, SUPPORT)
            
            % Construct a continuous distribution from name, number of
            % parameters, and range
            self@Distribution(name,numParam,support);
        end
    end

    methods
        function L = evalLST(self, s)
           line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end
    end
    
end
