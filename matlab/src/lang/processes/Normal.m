classdef Normal < ContinuousDistribution
    % Normal Normal (Gaussian) distribution
    %
    % Normal represents a normal (Gaussian) distribution with mean mu and
    % standard deviation sigma. This distribution is widely used in modeling
    % continuous phenomena due to the central limit theorem.
    %
    % @brief Normal distribution with specified mean and standard deviation
    %
    % Key characteristics:
    % - Two parameters: mean (mu) and standard deviation (sigma)
    % - Variance = sigma^2
    % - SCV = sigma^2 / mu^2
    % - Support: (-Inf, Inf)
    % - Symmetric with zero skewness
    %
    % Note: For queueing applications, a truncated or shifted version may
    % be needed since normal distributions can take negative values.
    %
    % Example:
    % @code
    % dist = Normal(5.0, 1.0);  % Mean=5, StdDev=1
    % samples = dist.sample(1000);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = Normal(mu, sigma)
            % NORMAL Create a Normal distribution instance
            %
            % @brief Creates a Normal distribution with specified mean and std dev
            % @param mu Mean of the distribution
            % @param sigma Standard deviation of the distribution (must be positive)
            % @return self Normal distribution instance
            self@ContinuousDistribution('Normal', 2, [-Inf, Inf]);
            if sigma < 0
                line_error(mfilename, 'sigma parameter must be >= 0.0');
            end
            setParam(self, 1, 'mu', mu);
            setParam(self, 2, 'sigma', sigma);
        end

        function ex = getMean(self)
            % EX = GETMEAN()

            % Get distribution mean
            ex = self.getParam(1).paramValue;
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()

            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            mu = self.getParam(1).paramValue;
            sigma = self.getParam(2).paramValue;
            if abs(mu) < GlobalConstants.FineTol
                SCV = Inf;
            else
                SCV = (sigma * sigma) / (mu * mu);
            end
        end

        function VAR = getVar(self)
            % VAR = GETVAR()

            % Get distribution variance
            sigma = self.getParam(2).paramValue;
            VAR = sigma * sigma;
        end

        function STD = getStd(self)
            % STD = GETSTD()

            % Get distribution standard deviation
            STD = self.getParam(2).paramValue;
        end

        function SKEW = getSkewness(self)
            % SKEW = GETSKEWNESS()

            % Get distribution skewness (always 0 for normal)
            SKEW = 0.0;
        end

        function Ft = evalCDF(self, t)
            % FT = EVALCDF(SELF,T)

            % Evaluate the cumulative distribution function at t
            mu = self.getParam(1).paramValue;
            sigma = self.getParam(2).paramValue;
            Ft = 0.5 * (1 + erf((t - mu) / (sigma * sqrt(2))));
        end

        function L = evalLST(self, s)
            % L = EVALST(S)

            % Evaluate the Laplace-Stieltjes transform of the distribution
            % function at s: E[e^(-sX)] = exp(mu*s + 0.5*sigma^2*s^2)
            % Note: this is the moment-generating function evaluated at -s
            mu = self.getParam(1).paramValue;
            sigma = self.getParam(2).paramValue;
            L = exp(mu * s + 0.5 * sigma * sigma * s * s);
        end

        function X = sample(self, n)
            % X = SAMPLE(N)

            % Get n samples from the distribution
            if nargin < 2
                n = 1;
            end
            mu = self.getParam(1).paramValue;
            sigma = self.getParam(2).paramValue;
            X = mu + sigma * randn(n, 1);
        end

        function proc = getProcess(self)
            % PROC = GETPROCESS()

            % Get process representation with actual distribution parameters
            % Returns {mu, sigma} for normal distribution
            proc = {self.getParam(1).paramValue, self.getParam(2).paramValue};
        end
    end

    methods (Static)
        function nm = fitMean(MEAN)
            % NM = FITMEAN(MEAN)

            % Fit distribution with given mean and unit standard deviation
            nm = Normal(MEAN, 1.0);
        end

        function nm = fitMeanAndStd(MEAN, STD)
            % NM = FITMEANANDSTD(MEAN, STD)

            % Fit distribution with given mean and standard deviation
            nm = Normal(MEAN, max(GlobalConstants.FineTol, STD));
        end

        function nm = fitMeanAndVar(MEAN, VAR)
            % NM = FITMEANANDVAR(MEAN, VAR)

            % Fit distribution with given mean and variance
            nm = Normal(MEAN, max(GlobalConstants.FineTol, sqrt(VAR)));
        end

        function nm = fitMeanAndSCV(MEAN, SCV)
            % NM = FITMEANANDSCV(MEAN, SCV)

            % Fit distribution with given mean and squared coefficient of
            % variation (SCV = variance / mean^2)
            VAR = SCV * MEAN^2;
            nm = Normal(MEAN, max(GlobalConstants.FineTol, sqrt(VAR)));
        end
    end

end
