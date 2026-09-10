classdef Gamma < ContinuousDistribution
    % Gamma General-purpose continuous distribution with shape and scale
    %
    % Gamma represents the gamma distribution with shape (alpha) and scale (beta)
    % parameters. This versatile distribution can model various phenomena and
    % includes the exponential and Erlang distributions as special cases. It is
    % commonly used for modeling service times and inter-arrival times.
    %
    % @brief Gamma distribution with flexible shape and scale parameters
    %
    % Key characteristics:
    % - Two parameters: shape (alpha) and scale (beta)
    % - Mean = alpha * beta, Variance = alpha * beta²
    % - SCV = 1/alpha (decreases with shape parameter)
    % - Support: (0, ∞)
    % - Includes exponential (alpha=1) and Erlang as special cases
    %
    % The gamma distribution is used for:
    % - Service time modeling with flexible variability
    % - Waiting time distributions
    % - Inter-arrival time modeling
    % - Reliability and survival analysis
    % - Approximating other continuous distributions
    %
    % Example:
    % @code
    % service_dist = Gamma(2.0, 1.5);  % Shape=2, Scale=1.5, Mean=3.0, SCV=0.5
    % samples = service_dist.sample(1000);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = Gamma(shape, scale)
            % GAMMA Create a Gamma distribution instance
            %
            % @brief Creates a Gamma distribution with specified shape and scale
            % @param shape Shape parameter alpha (must be positive)
            % @param scale Scale parameter beta (must be positive)
            % @return self Gamma distribution instance
            self@ContinuousDistribution('Gamma',2,[0,Inf]);
            setParam(self, 1, 'alpha', shape);
            setParam(self, 2, 'beta', scale);
        end
    end
    
    methods
        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Get distribution mean
            shape = self.getParam(1).paramValue;
            scale = self.getParam(2).paramValue;
            ex = shape*scale;
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            
            shape = self.getParam(1).paramValue;
            SCV = 1 / shape;
        end
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            
            % Get n samples from the distribution
            if nargin<2 %~exist('n','var'), 
                n = 1; 
            end
            shape = self.getParam(1).paramValue;
            scale = self.getParam(2).paramValue;
            X = gamrnd(shape, scale, n, 1);
        end
        
        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            
            % Evaluate the cumulative distribution function at t
            % AT T
            
            shape = self.getParam(1).paramValue;
            scale = self.getParam(2).paramValue;
            Ft = gamcdf(t,shape,scale);
        end
        
        function L = evalLST(self, s)
            % L = EVALLAPLACETRANSFORM(S)
            
            % Evaluate the Laplace transform of the distribution function at t
            
            shape = self.getParam(1).paramValue; %alpha
            scale = self.getParam(2).paramValue; %beta
            L = ((1/scale) / (s+1/scale))^shape;
        end
        
        function proc = getProcess(self)
            % PROC = GETPROCESS()

            % Get process representation with actual distribution parameters
            % Returns {shape, scale} for gamma distribution
            proc = {self.getParam(1).paramValue, self.getParam(2).paramValue};
        end
    end
    
    methods(Static)
        function gm = fitMeanAndSCV(MEAN, SCV)
            % GM = FITMEANANDSCV(MEAN, SCV)
            
            % Fit distribution from mean and squared coefficient of
            % variation
            shape = 1 / SCV;
            scale = MEAN / shape;
            gm = Gamma(shape, scale);
        end
    end
    
end

