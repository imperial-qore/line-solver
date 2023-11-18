classdef Det < ContinuousDistribution & DiscreteDistribution
    % Deterministic distribution
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    properties
    end
    
    methods
        function self = Det(t)
            % SELF = DET(T)
            
            % Construct a deterministic distribution with value t
            self@ContinuousDistribution('Det',1,[t,t]);
            self@DiscreteDistribution('Det',1,[t,t]);
            setParam(self, 1, 't', t);
        end
        
        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Get distribution mean
            ex = self.getParam(1).paramValue;
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            SCV = 0;
        end
        
        function L = evalLST(self, s)
            % L = EVALST(S)
            
            % Evaluate the Laplace-Stieltjes transform of the distribution function at t            
            
            t = self.getParam(1).paramValue;
            
            L = exp(-s*t);
        end
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            
            % Get n samples from the distribution
            X = self.getParam(1).paramValue * ones(n,1);
        end
        
        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            
            % Evaluate the cumulative distribution function at t
            % AT T
            
            if t < self.getParam(1).paramValue
                Ft = 0;
            else
                Ft = 1;
            end
        end
    end
    
end

