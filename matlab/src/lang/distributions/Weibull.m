classdef Weibull < ContinuousDistribution
    % The Weibull statistical distribution
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    methods
        function self = Weibull(shape, scale)
            % SELF = WEIBULL(SHAPE, SCALE)
            
            % Constructs a Weibull distribution with given shape and scale
            % parameters
            self@ContinuousDistribution('Weibull',2,[0,Inf]);
            if shape < 0
                line_error(mfilename,'shape parameter must be >= 0.0');
            end
            setParam(self, 1, 'alpha', scale);
            setParam(self, 2, 'r', shape);
        end
        
        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Get distribution mean
            alpha = self.getParam(1).paramValue;
            r = self.getParam(2).paramValue;
            ex = alpha * gamma(1+1/r);
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            alpha = self.getParam(1).paramValue;
            r = self.getParam(2).paramValue;
            VAR = alpha^2*(gamma(1+2/r)-(gamma(1+1/r))^2);
            ex = alpha * gamma(1+1/r);
            SCV = VAR / ex^2;
        end
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            
            % Get n samples from the distribution
            if nargin<2 %~exist('n','var'), 
                n = 1; 
            end
            alpha = self.getParam(1).paramValue;
            r = self.getParam(2).paramValue;
            X = wblrnd(alpha,r,n,1);
        end
        
        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            
            % Evaluate the cumulative distribution function at t
            % AT T            
            alpha = self.getParam(1).paramValue;
            r = self.getParam(2).paramValue;
            Ft = wblcdf(t,alpha,r);
        end
        
        function L = evalLST(self, sl)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform of the distribution function at t
            error('Laplace-Stieltjes transform of the Weibull distribution not available yet.');
        end
    end
    
    methods (Static)
        function pa = fitMeanAndSCV(MEAN, SCV)
            % PA = FITMEANANDSCV(MEAN, SCV)
            
            % Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
            c = sqrt(SCV);
            r = c^(-1.086); % Justus approximation (1976)
            alpha = MEAN / gamma(1+1/r);
            pa = Weibull(r,alpha);
        end
    end
    
end
