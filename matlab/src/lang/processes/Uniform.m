classdef Uniform < ContinuousDistribution
    % Uniform Continuous uniform distribution over a bounded interval
    %
    % Uniform represents the continuous uniform distribution with equal probability
    % density over a specified interval [min, max]. This distribution models
    % scenarios where all values within a range are equally likely, making it
    % useful for modeling bounded random processes with no preferred values.
    %
    % @brief Uniform distribution with equal probability over [min, max] interval
    %
    % Key characteristics:
    % - Two parameters: minimum (min) and maximum (max) values
    % - Mean = (min + max) / 2
    % - Variance = (max - min)² / 12
    % - Constant probability density over [min, max]
    % - Support: [min, max]
    %
    % The uniform distribution is used for:
    % - Modeling bounded processes with no bias
    % - Random number generation in simulations
    % - Representing equally likely outcomes
    % - Initial parameter estimation
    % - Load testing with uniform traffic patterns
    %
    % Example:
    % @code
    % service_dist = Uniform(1.0, 3.0);  % Service time between 1 and 3 units
    % % Mean = 2.0, all values in [1.0, 3.0] equally likely
    % samples = service_dist.sample(1000);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = Uniform(minVal, maxVal)
            % UNIFORM Create a Uniform distribution instance
            %
            % @brief Creates a Uniform distribution over the interval [minVal, maxVal]
            % @param minVal Minimum value of the distribution support
            % @param maxVal Maximum value of the distribution support (must be > minVal)
            % @return self Uniform distribution instance
            self@ContinuousDistribution('Uniform',2,[minVal,maxVal]);
            setParam(self, 1, 'min', minVal);
            setParam(self, 2, 'max', maxVal);
        end
        
        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Get distribution mean
            ex = (self.getParam(2).paramValue+self.getParam(1).paramValue) / 2;
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            var = (self.getParam(2).paramValue-self.getParam(1).paramValue)^2 / 12;
            SCV = var/getMean(self)^2;
        end
        
        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            
            % Evaluate the cumulative distribution function at t
            % AT T
            
            minVal = self.getParam(1).paramValue;
            maxVal = self.getParam(2).paramValue;
            if t < minVal
                Ft = 0;
            elseif t > maxVal
                Ft = 0;
            else
                Ft = 1/(maxVal-minVal);
            end
        end
        
        function L = evalLST(self, s)
            % L = EVALST(S)
            
            % Evaluate the Laplace-Stieltjes transform of the distribution function at t
            
            minVal = self.getParam(1).paramValue;
            maxVal = self.getParam(2).paramValue;
            
            L = ( exp(-s*minVal) - exp(-s*maxVal) )/(s*(maxVal - minVal));
        end
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            
            % Get n samples from the distribution
            if nargin<2 %~exist('n','var'), 
                n = 1; 
            end
            minVal = self.getParam(1).paramValue;
            maxVal = self.getParam(2).paramValue;
            X = minVal + (maxVal-minVal)*rand(n,1);
        end
        
        function proc = getProcess(self)
            % PROC = GETPROCESS()

            % Get process representation with actual distribution parameters
            % Returns {min, max} for uniform distribution
            proc = {self.getParam(1).paramValue, self.getParam(2).paramValue};
        end
    end
    
end

