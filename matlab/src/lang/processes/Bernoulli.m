classdef Bernoulli < DiscreteDistribution
    % Bernoulli Binary outcome distribution for success/failure modeling
    %
    % Bernoulli represents a Bernoulli distribution with probability parameter p.
    % This is the simplest discrete distribution, modeling a single trial with
    % two possible outcomes: success (1) with probability p, and failure (0)
    % with probability 1-p. It forms the foundation for binomial distributions.
    %
    % @brief Bernoulli distribution for binary success/failure outcomes
    %
    % Key characteristics:
    % - Single parameter: success probability p ∈ [0,1]
    % - Mean = p, Variance = p(1-p)
    % - SCV = (1-p)/(p) 
    % - Support: {0, 1}
    % - Building block for binomial and geometric distributions
    %
    % The Bernoulli distribution is used for:
    % - Modeling binary outcomes (success/failure, yes/no)
    % - Cache hit/miss modeling
    % - Component reliability analysis
    % - Binary decision processes
    % - Foundation for more complex discrete distributions
    %
    % Example:
    % @code
    % coin_flip = Bernoulli(0.5);      % Fair coin (50% success)
    % biased_coin = Bernoulli(0.7);    % Biased coin (70% success)
    % samples = coin_flip.sample(1000);
    % @endcode
    %
    % Copyright (c) 2018-2022, Imperial College London
    % All rights reserved.

    methods
        function self = Bernoulli(p)
            % BERNOULLI Create a Bernoulli distribution instance
            %
            % @brief Creates a Bernoulli distribution with success probability p
            % @param p Success probability (must be in [0,1])
            % @return self Bernoulli distribution instance

            self@DiscreteDistribution('Bernoulli',1,[0,Inf]);
            setParam(self, 1, 'p', p);
        end

        function ex = getMean(self)
            % EX = GETMEAN()

            % Get distribution mean
            n = 1;
            p = self.getParam(1).paramValue;
            ex = p*n;
        end

        function V = getVar(self)
            % V = GETVAR()
            n = 1;
            p = self.getParam(1).paramValue;
            V = n*p*(1-p);
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()

            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            p = self.getParam(1).paramValue;
            n = 1;
            SCV = (1-p) / (n*p);
        end

        function X = sample(self, nsamples)
            % X = SAMPLE(N)

            X = zeros(1, nsamples);
            p = self.getParam(1).paramValue;
            n = 1;

            for k = 1:nsamples
                acc = 0;
                for i = 1:n
                    if rand <= p
                        acc = acc + 1;
                    end
                end
                X(k) = acc;
            end
        end

        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)

            % Evaluate the cumulative distribution function at T

            p = self.getParam(1).paramValue;
            n = 1;
            Ft = 0.0;
            if t >= n
                Ft = 1.0;
            elseif t >= 0
                Ft = 1.0 - betainc(p, t + 1.0, n - t);
            end
        end

        function L = evalLST(self, s)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform of the distribution function at s
            % For Bernoulli(p), LST(s) = (1 - p + p*e^(-s))

            p = self.getParam(1).paramValue;
            L = 1 - p + p * exp(-s);
        end

        function P = evalPMF(self, k)
            % P = EVALPMF(K)

            % Evaluate the probability mass function at k            

            p = self.getParam(1).paramValue;
            n = 1;
            P = binopdf(k, n, p);
        end
        
        function proc = getProcess(self)
            % PROC = GETPROCESS()
            
            % Get process representation for non-Markovian distribution
            % Returns [mean, SCV] pair for use in network analysis
            proc = [self.getMean(), self.getSCV()];
        end
    end

end

