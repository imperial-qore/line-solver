classdef Weibull < ContinuousDistribution
    % Weibull Distribution for reliability and failure time modeling
    %
    % Weibull represents the Weibull distribution with shape and scale parameters.
    % This distribution is widely used in reliability engineering and survival
    % analysis for modeling failure times, life distributions, and extreme value
    % phenomena. It can model increasing, decreasing, or constant hazard rates.
    %
    % @brief Weibull distribution for reliability and failure time analysis
    %
    % Key characteristics:
    % - Two parameters: scale (alpha) and shape (r)
    % - Mean = alpha * Γ(1 + 1/r) where Γ is the gamma function
    % - Flexible hazard rate behavior based on shape parameter
    % - Support: (0, ∞)
    % - Includes exponential (r=1) as special case
    %
    % Hazard rate behavior:
    % - r < 1: Decreasing failure rate (infant mortality)
    % - r = 1: Constant failure rate (exponential, random failures)  
    % - r > 1: Increasing failure rate (wear-out failures)
    %
    % The Weibull distribution is used for:
    % - Reliability and survival analysis
    % - Failure time modeling
    % - Lifetime data analysis
    % - Wind speed modeling
    % - Service time distributions with varying hazard rates
    %
    % Example:
    % @code
    % failure_dist = Weibull(2.0, 1000);  % Scale=1000, Shape=2 (wear-out)
    % reliability_time = failure_dist.sample(100);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = Weibull(shape, scale)
            % WEIBULL Create a Weibull distribution instance
            %
            % @brief Creates a Weibull distribution with specified shape and scale
            % @param shape Shape parameter r (must be positive)
            % @param scale Scale parameter alpha (must be positive)  
            % @return self Weibull distribution instance
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
        
        function L = evalLST(self, s)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform of the distribution function at s
            % The LST of Weibull distribution doesn't have a simple closed form.
            % We compute it numerically using the definition: E[e^(-sX)] = ∫₀^∞ e^(-sx) f(x) dx
            % where f(x) is the Weibull PDF
            
            alpha = self.getParam(1).paramValue; % scale parameter
            r = self.getParam(2).paramValue; % shape parameter
            
            % Use numerical integration
            % For practical purposes, we integrate up to a reasonable upper bound
            upperBound = alpha * (-log(1e-10))^(1/r); % covers most of the distribution
            n = 1000;
            dx = upperBound / n;
            L = 0.0;
            
            for i = 1:n
                x = i * dx;
                if x > 0
                    % Weibull PDF: (r/α)(x/α)^(r-1) * exp(-(x/α)^r)
                    pdf = (r / alpha) * (x / alpha)^(r - 1) * exp(-(x / alpha)^r);
                    L = L + exp(-s * x) * pdf;
                end
            end
            
            L = L * dx;
        end
        
        function proc = getProcess(self)
            % PROC = GETPROCESS()

            % Get process representation with actual distribution parameters
            % Returns {shape (r), scale (alpha)} for Weibull distribution
            proc = {self.getParam(2).paramValue, self.getParam(1).paramValue};
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
