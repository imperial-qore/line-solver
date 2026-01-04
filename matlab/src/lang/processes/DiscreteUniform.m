classdef DiscreteUniform < DiscreteDistribution
    % The uniform statistical distribution
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = DiscreteUniform(minVal, maxVal)
            % SELF = UNIFORM(MINVAL, MAXVAL)
            
            % Constructs an uniform distribution with specified minimum and
            % maximum values
            self@DiscreteDistribution('DiscreteUniform',2,[minVal,maxVal]);
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
            var = ((self.getParam(2).paramValue-self.getParam(1).paramValue+1)^2-1) / 12;
            SCV = var/getMean(self)^2;
        end
        
        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            
            % Evaluate the cumulative distribution function at T
            
            minVal = self.getParam(1).paramValue;
            maxVal = self.getParam(2).paramValue;
            if t < minVal
                Ft = 0;
            elseif t > maxVal
                Ft = 0;
            else
                Ft = (floor(t)-minVal+1)/(maxVal-minVal+1);
            end
        end

        function L = evalLST(self, s)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform of the distribution function at s
            % For DiscreteUniform(a, b), LST(s) = (e^(-as) - e^(-(b+1)s)) / ((b-a+1)(1 - e^(-s)))

            a = self.getParam(1).paramValue;
            b = self.getParam(2).paramValue;
            
            if abs(s) < 1e-10
                % Handle s ≈ 0 case to avoid division by zero
                L = 1.0;
            else
                e_neg_s = exp(-s);
                numerator = exp(-a * s) - exp(-(b + 1) * s);
                denominator = (b - a + 1) * (1 - e_neg_s);
                L = numerator / denominator;
            end
        end
        
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            
            % Get n samples from the distribution
            if nargin<2 %~exist('n','var'), 
                n = 1; 
            end
            minVal = self.getParam(1).paramValue;
            maxVal = self.getParam(2).paramValue;
            X = round(minVal + (maxVal-minVal)*rand(n,1));
        end
        
        function proc = getProcess(self)
            % PROC = GETPROCESS()
            
            % Get process representation for non-Markovian distribution
            % Returns [mean, SCV] pair for use in network analysis
            proc = [self.getMean(), self.getSCV()];
        end
    end
    
end

