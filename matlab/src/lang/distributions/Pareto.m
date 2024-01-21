classdef Pareto < ContinuousDistribution
    % The Pareto statistical distribution
    %
    % Copyright (c) 2012-2024, Imperial College London
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
        
        function L = evalLST(self, sl)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform of the distribution function at t
            
            % Saralees Nadarajah & Samuel Kotz. 
            % On the Laplace transform of the Pareto distribution
            % Queueing Syst (2006) 54:243-244
            % DOI 10.1007/s11134-006-0299-1                                   
            L = [];
            for s=sl(:)'
                alpha = self.getParam(1).paramValue;
                k = self.getParam(2).paramValue;            
                IL = integral(@(t)t.^(-alpha-1).*exp(-t),s*k,Inf); 
                L(end+1) = alpha*k^alpha*IL*s^(1+alpha)/s;
                %a = shape; b = scale; L = integral(@(x) a*b^a*exp(-s*x)./(x).^(a+1),b,1e6);
            end
        end
    end
    
    methods (Static)
        function pa = fitMeanAndSCV(MEAN, SCV)
            % PA = FITMEANANDSCV(MEAN, SCV)
            
            % Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
            shape = (SCV*MEAN + MEAN*(SCV*(SCV + 1))^(1/2))/(SCV*MEAN);
            scale = MEAN + SCV*MEAN - MEAN*(SCV*(SCV + 1))^(1/2);
            pa = Pareto(shape,scale);
        end
    end
    
end
