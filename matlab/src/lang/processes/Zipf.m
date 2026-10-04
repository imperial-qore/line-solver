classdef Zipf < DiscreteDistribution
    % Zipf Power-law distribution for modeling popularity and ranking
    %
    % Zipf represents the Zipf distribution commonly used for modeling popularity
    % distributions, where items are ranked by frequency. The probability of an
    % item is inversely proportional to its rank raised to the shape parameter.
    % This distribution is essential for cache modeling and content popularity.
    %
    % @brief Zipf distribution for popularity and ranking-based phenomena
    %
    % Key characteristics:
    % - Two parameters: shape (s) and number of items (n)
    % - Probability ∝ 1/rank^s for each rank
    % - Heavy-tailed distribution (popular items dominate)
    % - Mean = H(s-1,n)/H(s,n) where H is generalized harmonic number
    % - Support: {1, 2, ..., n}
    %
    % Shape parameter behavior:
    % - s > 1: Heavy tail, few items very popular
    % - s = 1: Classic Zipf law
    % - s < 1: Less skewed, more uniform popularity
    %
    % The Zipf distribution is used for:
    % - Web cache popularity modeling
    % - Content access patterns
    % - Word frequency in natural language
    % - City population distributions
    % - Internet traffic modeling
    %
    % Example:
    % @code
    % web_popularity = Zipf(1.2, 1000);  % 1000 items, shape=1.2
    % access_pattern = web_popularity.sample(10000);
    % @endcode
    %
    % Copyright (c) 2018-2022, Imperial College London
    % All rights reserved.

    methods
        function self = Zipf(s, n)
            % ZIPF Create a Zipf distribution instance
            %
            % @brief Creates a Zipf distribution with shape parameter and item count
            % @param s Shape parameter (s > 0, controls skewness)
            % @param n Number of items (positive integer)
            % @return self Zipf distribution instance with specified parameters
            
            self@DiscreteDistribution('Zipf',4,[1,n]);
            p = 1./((1:n).^s)/Zipf.genHarmonic(s,n);
            x = 1:n;
            setParam(self, 1, 'p', p(:)');
            setParam(self, 2, 'x', x(:)');
            setParam(self, 3, 's', s);
            setParam(self, 4, 'n', n);
        end

        function ex = getMean(self)
            % EX = GETMEAN()

            % Get distribution mean
            s = self.getParam(3).paramValue;
            n = self.getParam(4).paramValue;
            ex = self.genHarmonic(s-1,n) / self.genHarmonic(s,n);
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()

            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            s = self.getParam(3).paramValue;
            n = self.getParam(4).paramValue;
            ex = getMean(self);
            var = self.genHarmonic(s-2,n) / self.genHarmonic(s,n) - ex^2;
            SCV = var / ex^2;
        end

        function X = sample(self, n)
            % X = SAMPLE(N)
            % to be checked

            alpha = self.getParam(3).paramValue;
            ranks = self.getParam(2).paramValue;

            pmf = (1./ranks.^alpha) / sum(1./ranks.^alpha);
            cdf = cumsum(pmf);
            uniformRandomNumbers = rand(n, 1);
            X = arrayfun(@(x) find(cdf >= x, 1, 'first'), uniformRandomNumbers);
        end

        function Ft = evalCDF(self,k)
            % FT = EVALCDF(SELF,K)

            % Evaluate the cumulative distribution function at t
            % AT T

            s = self.getParam(3).paramValue;
            n = self.getParam(4).paramValue;
            Ft = self.genHarmonic(s,k) / self.genHarmonic(s,n);
        end

        function L = evalLST(self, s)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform of the distribution function at s
            % For Zipf distribution, LST(s) = Σ(k=1 to n) P(X=k) * e^(-s*k)

            shape = self.getParam(3).paramValue;
            n = self.getParam(4).paramValue;
            Hns = self.genHarmonic(shape, n);
            
            L = 0.0;
            for k = 1:n
                prob = 1.0 / k^shape / Hns;
                L = L + prob * exp(-s * k);
            end
        end

        function p = evalPMF(self, k)
            % P = EVALPMF(K)

            % Evaluate the probability mass function at k
            % AT K

            s = self.getParam(3).paramValue;
            n = self.getParam(4).paramValue;
            if nargin<2 %~exist('k','var')
                k = 1:n;
            end
            Hns = Zipf.genHarmonic(s,n);
            p = 1./(k.^s)/Hns;
        end
        
        function proc = getProcess(self)
            % PROC = GETPROCESS()
            
            % Get process representation for non-Markovian distribution
            % Returns [mean, SCV] pair for use in network analysis
            proc = [self.getMean(), self.getSCV()];
        end
    end

    methods (Static)
        function Hnm = genHarmonic(s,n)
            % HNM = GENHARMONIC(S,N)

            % Generate harmonic numbers to normalize a Zipf-like distribution
            % on n items with shape parameter s
            Hnm = 0;
            for k=1:n
                Hnm = Hnm + 1/k^s;
            end
        end
    end

end

