classdef Pareto < ContinuousDistribution
    % The Pareto statistical distribution
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = Pareto(shape, scale)
            % SELF = PARETO(SHAPE, SCALE)
            
            % Constructs a Pareto distribution with given shape and scale
            % parameters
            self@ContinuousDistribution('Pareto',2,[0,Inf]);
            if shape < 2
                line_error(mfilename,'shape parameter must be >= 2.0');
            end
            setParam(self, 1, 'alpha', shape);
            setParam(self, 2, 'k', scale);
        end
        
        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Get distribution mean
            shape = self.getParam(1).paramValue;
            scale = self.getParam(2).paramValue;
            ex = shape * scale / (shape - 1);
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            shape = self.getParam(1).paramValue;
            scale = self.getParam(2).paramValue;
            VAR = scale^2 * shape / (shape - 1)^2 / (shape - 2);
            ex = shape * scale / (shape - 1);
            SCV = VAR / ex^2;
        end
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            
            % Get n samples from the distribution
            if nargin<2 %~exist('n','var'), 
                n = 1; 
            end
            shape = self.getParam(1).paramValue;
            scale = self.getParam(2).paramValue;
            k = 1/shape;
            sigma = scale * k;
            X = gprnd(k, sigma, sigma/k, n, 1);
        end
        
        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            
            % Evaluate the cumulative distribution function at t
            % AT T
            
            shape = self.getParam(1).paramValue;
            scale = self.getParam(2).paramValue;
            k = 1/shape;
            sigma = scale * k;
            Ft = gpcdf(t, k, sigma, sigma/k);
        end
        
        function L = evalLST(self, s)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform of the distribution function at s
            % The LST of Pareto distribution doesn't have a simple closed form.
            % We compute it numerically using the definition: E[e^(-sX)] = ∫_{k}^∞ e^(-sx) f(x) dx
            % where f(x) is the Pareto PDF: α*k^α / x^(α+1) for x >= k
            
            alpha = self.getParam(1).paramValue; % shape parameter
            k = self.getParam(2).paramValue; % scale parameter
            
            % Use numerical integration
            % For practical purposes, we integrate up to a reasonable upper bound
            upperBound = k * (1000)^(1/alpha); % covers most of the distribution
            n = 1000;
            dx = (upperBound - k) / n;
            L = 0.0;
            
            for i = 1:n
                x = k + i * dx;
                if x >= k
                    % Pareto PDF: α*k^α / x^(α+1)
                    pdf = alpha * k^alpha / x^(alpha + 1);
                    L = L + exp(-s * x) * pdf;
                end
            end
            
            L = L * dx;
            
            % OLD IMPLEMENTATION (using integral function):
            % Saralees Nadarajah & Samuel Kotz. 
            % On the Laplace transform of the Pareto distribution
            % Queueing Syst (2006) 54:243-244
            % DOI 10.1007/s11134-006-0299-1                                   
            % L = [];
            % for s=sl(:)'
            %     alpha = self.getParam(1).paramValue;
            %     k = self.getParam(2).paramValue;            
            %     IL = integral(@(t)t.^(-alpha-1).*exp(-t),s*k,Inf); 
            %     L(end+1) = alpha*k^alpha*IL*s^(1+alpha)/s;
            %     %a = shape; b = scale; L = integral(@(x) a*b^a*exp(-s*x)./(x).^(a+1),b,1e6);
            % end
        end
        
        function proc = getProcess(self)
            % PROC = GETPROCESS()

            % Get process representation with actual distribution parameters
            % Returns {shape (alpha), scale (k)} for Pareto distribution
            proc = {self.getParam(1).paramValue, self.getParam(2).paramValue};
        end
    end
    
    methods (Static)
        function pa = fitMeanAndSCV(MEAN, SCV)
            % PA = FITMEANANDSCV(MEAN, SCV)

            % Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
            % For Pareto distribution with shape alpha and scale k:
            %   Mean = alpha*k / (alpha-1)
            %   SCV = 1 / (alpha*(alpha-2))
            %
            % Solving for alpha from SCV:
            %   alpha*(alpha-2) = 1/SCV
            %   alpha^2 - 2*alpha - 1/SCV = 0
            %   alpha = 1 + sqrt(1 + 1/SCV)  (taking positive root, need alpha > 2)
            %
            % Then scale k from mean:
            %   k = MEAN * (alpha-1) / alpha

            shape = 1 + sqrt(1 + 1/SCV);
            scale = MEAN * (shape - 1) / shape;
            pa = Pareto(shape, scale);
        end
    end
    
end
