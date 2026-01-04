classdef Lognormal < ContinuousDistribution
    % The Lognormal statistical distribution
    %
    % Copyright (c) 2012-2026, Imperial College London
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
        
        function L = evalLST(self, s)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform of the distribution function at s
            % The LST of Lognormal distribution doesn't have a simple closed form.
            % We compute it numerically using the definition: E[e^(-sX)] = ∫₀^∞ e^(-sx) f(x) dx
            % where f(x) is the lognormal PDF
            
            mu = self.getParam(1).paramValue;
            sigma = self.getParam(2).paramValue;
            
            % Use numerical integration
            % For practical purposes, we integrate up to a reasonable upper bound
            upperBound = exp(mu + 5 * sigma); % covers most of the distribution
            n = 1000;
            dx = upperBound / n;
            L = 0.0;
            
            for i = 1:n
                x = i * dx;
                if x > 0
                    logx = log(x);
                    % Lognormal PDF: exp(-(log(x)-μ)²/(2σ²)) / (x*σ*√(2π))
                    pdf = exp(-(logx - mu)^2 / (2 * sigma^2)) / (x * sigma * sqrt(2 * pi));
                    L = L + exp(-s * x) * pdf;
                end
            end
            
            L = L * dx;
        end
        
        function proc = getProcess(self)
            % PROC = GETPROCESS()

            % Get process representation with actual distribution parameters
            % Returns {mu, sigma} for lognormal distribution
            proc = {self.getParam(1).paramValue, self.getParam(2).paramValue};
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
