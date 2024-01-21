classdef Uniform < ContinuousDistribution
    % The uniform statistical distribution
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    methods
        function self = Uniform(minVal, maxVal)
            % SELF = UNIFORM(MINVAL, MAXVAL)
            
            % Constructs an uniform distribution with specified minimum and
            % maximum values
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
    end
    
end

