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
            %
            % A*(s) = E[e^{-sX}] = int_k^Inf e^{-sx} alpha k^alpha x^{-(alpha+1)} dx.
            % Substituting x = k/u maps the infinite tail onto a unit interval and
            % cancels the scale exactly:
            %
            %     A*(s) = alpha * int_0^1 u^(alpha-1) exp(-s*k/u) du
            %
            % This is the same transform as the closed form of Nadarajah & Kotz,
            % A*(s) = alpha*(s*k)^alpha*Gamma(-alpha, s*k) = alpha*E_{alpha+1}(s*k)
            % (Queueing Syst (2006) 54:243-244, DOI 10.1007/s11134-006-0299-1),
            % but in a form that stays accurate as s -> 0, where the incomplete-gamma
            % product underflows to 0/Inf. Here s = 0 gives alpha*int_0^1 u^(alpha-1)
            % du = 1 exactly, and the integrand is bounded and C^Inf on a FINITE
            % interval for alpha >= 2 (the shape floor the constructor enforces).
            %
            % Accuracy: adaptive Gauss-Kronrod at RelTol 1e-12, i.e. ~1e-12 relative,
            % verified against mpmath to 1e-15. The previous implementation was a
            % 1000-point right-endpoint rectangle sum truncated at k*1000^(1/alpha);
            % it lost the mass beyond the truncation point and biased the transform
            % low by ~3.1% at alpha=2.0078 (it returned A*(0)=0.96914, not 1).

            alpha = self.getParam(1).paramValue; % shape parameter
            k = self.getParam(2).paramValue; % scale parameter

            L = zeros(size(s));
            for i = 1:numel(s)
                si = s(i);
                if si == 0
                    L(i) = 1; % A*(0) = 1 exactly; skip the quadrature
                else
                    % u=0 is an essential zero of the integrand (exp(-s*k/u) and all
                    % its derivatives vanish there), so the guard only avoids 0/0.
                    % AbsTol is set at the underflow floor rather than a "small"
                    % value on purpose: A*(s) spans hundreds of decades (9.3e-220
                    % at alpha=50, k=0.5, s=1000), so any meaningful absolute floor
                    % lets the quadrature stop before doing real work and return a
                    % value with no correct digits.
                    g = @(u) (u > 0) .* u.^(alpha-1) .* exp(-si .* k ./ max(u, realmin));
                    L(i) = alpha * integral(g, 0, 1, 'RelTol', 1e-12, 'AbsTol', realmin);
                end
            end
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
