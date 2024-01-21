classdef Lognormal < ContinuousDistribution
    % The Lognormal statistical distribution
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    methods
        function self = Lognormal(mu, sigma)
            % SELF = LOGNORMAL(MU, SIGMA)
            
            % Constructs a Lognormal distribution with given mu and sigma
            % parameters
            self@ContinuousDistribution('Lognormal',2,[0,Inf]);
            if sigma < 0
                line_error(mfilename,'sigma parameter must be >= 0.0');
            end
            setParam(self, 1, 'mu', mu);
            setParam(self, 2, 'sigma', sigma);
        end
        
        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Get distribution mean
            mu = self.getParam(1).paramValue;
            sigma = self.getParam(2).paramValue;
            ex =  exp(mu+sigma*sigma/2);
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            mu = self.getParam(1).paramValue;
            sigma = self.getParam(2).paramValue;
            ex =  exp(mu+sigma*sigma/2);
            VAR = (exp(sigma*sigma)-1)*exp(2*mu+sigma*sigma);
            SCV = VAR / ex^2;
        end
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            
            % Get n samples from the distribution
            if nargin<2 %~exist('n','var'), 
                n = 1; 
            end
            mu = self.getParam(1).paramValue;
            sigma = self.getParam(2).paramValue;
            X = lognrnd(mu,sigma,n,1);
        end
        
        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            
            % Evaluate the cumulative distribution function at t
            % AT T            
            mu = self.getParam(1).paramValue;
            sigma = self.getParam(2).paramValue;
            Ft = logncdf(t,mu,sigma);
        end
        
        function L = evalLST(self, sl)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform of the distribution function at t
            error('Laplace-Stieltjes transform of the Lognormal distribution not available yet.');
        end
    end
    
    methods (Static)
        function pa = fitMeanAndSCV(MEAN, SCV)
            % PA = FITMEANANDSCV(MEAN, SCV)
            
            % Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
            c = sqrt(SCV);
            mu = log(MEAN  / sqrt(c*c + 1));
            sigma = sqrt(log(c*c + 1));            
            pa = Lognormal(mu,sigma);
        end
    end
    
end
