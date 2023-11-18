classdef Gamma < ContinuousDistribution
    % The gamma statistical distribution
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    methods
        function self = Gamma(shape, scale)
            % SELF = GAMMA(SHAPE, SCALE)
            
            % Constructs a gamma distribution from shape and scale
            % parameters
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

