classdef Poisson < DiscreteDistribution
    % Poisson Discrete distribution for counting rare events
    %
    % Poisson represents the Poisson distribution with rate parameter lambda.
    % This distribution models the number of events occurring in a fixed interval
    % when events occur independently at a constant average rate. It is commonly
    % used for modeling arrival counts, defects, and other rare event phenomena.
    %
    % @brief Poisson distribution for modeling discrete event counts
    %
    % Key characteristics:
    % - Single parameter: rate lambda (average number of events)
    % - Mean = Variance = lambda
    % - SCV = 1/lambda (decreases as lambda increases)
    % - Support: {0, 1, 2, 3, ...}
    % - Limit of binomial distribution as n→∞, p→0, np→lambda
    %
    % The Poisson distribution is used for:
    % - Counting arrivals in fixed time intervals
    % - Modeling defects or failures
    % - Call center arrival modeling
    % - Network packet arrival counts
    % - Population dynamics and birth processes
    %
    % Example:
    % @code
    % arrival_count = Poisson(3.5);  % Average 3.5 arrivals per interval
    % samples = arrival_count.sample(1000);  % Generate counts
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = Poisson(lambda)
            % POISSON Create a Poisson distribution instance
            %
            % @brief Creates a Poisson distribution with specified rate parameter
            % @param lambda Rate parameter (average number of events, must be positive)
            % @return self Poisson distribution instance
            self@DiscreteDistribution('Poisson', 1, [0,Inf]);
            setParam(self, 1, 'lambda', lambda);
            self.immediate = (1/lambda) < GlobalConstants.FineTol;
            self.obj = jline.lang.processes.Exp(lambda);
        end

        function X = sample(self, n)
            % X = SAMPLE(N)
            % Get n samples from the distribution
            lambda = self.getParam(1).paramValue;

            % Generate uniform random variables
            u = rand(n, 1);

            % Initialize output array
            X = zeros(n, 1);

            % Compute Poisson samples using inverse transform sampling
            for i = 1:n
                p = exp(-lambda);
                F = p;
                k = 0;
                while u(i) >= F
                    k = k + 1;
                    p = p * lambda / k;
                    F = F + p;
                end
                X(i) = k;
            end
        end

        function Ft = evalCDF(self,k)
            % FT = EVALCDF(SELF,K)
            % Evaluate the cumulative distribution function at K

            lambda = self.getParam(1).paramValue;
            Ft = poisscdf(k, lambda);
        end

        function L = evalLST(self, s)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform of the distribution function at s
            % For Poisson(λ), LST(s) = exp(λ(e^(-s) - 1))

            lambda = self.getParam(1).paramValue;
            L = exp(lambda * (exp(-s) - 1));
        end

        function mean = getMean(self)
            lambda = self.getParam(1).paramValue;
            mean = lambda;
        end

        function rate = getRate(self)
            lambda = self.getParam(1).paramValue;
            rate = 1 / lambda;
        end

        function scv = getSCV(self)
            lambda = self.getParam(1).paramValue;
            scv = 1.0 / lambda;
        end

        function scv = getVar(self)
            lambda = self.getParam(1).paramValue;
            scv = lambda;
        end        

        function skew = getSkewness(self)
            lambda = self.getParam(1).paramValue;
            skew = sqrt(lambda);
        end        

        function proc = getProcess(self)
            % PROC = GETPROCESS()
            
            % Get process representation for non-Markovian distribution
            % Returns [mean, SCV] pair for use in network analysis
            proc = [self.getMean(), self.getSCV()];
        end

    end

end