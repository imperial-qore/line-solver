classdef Geometric < DiscreteDistribution
    % A Geometric probability distribution
    %
    % The distribution of the number of Bernoulli trials needed to get 
    % one success.
    %
    % Copyright (c) 2018-2022, Imperial College London
    % All rights reserved.
    
    methods
        function self = Geometric(p)
            % SELF = GEOMETRIC(P)
            self@DiscreteDistribution('Geometric',1,[1,Inf]);
            % Construct a geometric distribution with probability p
            
            setParam(self, 1, 'p', p);
        end
        
        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Get distribution mean
            p = self.getParam(1).paramValue;

            ex = 1 / p;
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            p = self.getParam(1).paramValue;
            
            SCV = 1 - p;
        end
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            if nargin < 2
                n = 1;
            end
            % Get n samples from the distribution
            p = self.getParam(1).paramValue;
            r = rand(n,1);
            if p >= 1
                % Degenerate case: the first trial always succeeds. The
                % inversion below cannot express it, because log(1-p) is -Inf
                % and the quotient rounds to 0, outside the declared support
                % {1,2,...}. The draws above are still consumed so a parameter
                % sweep stays stream-synchronized.
                X = ones(n,1);
                return
            end
            X = ceil(log(1-r) ./ log(1-p));
        end
        
        function Ft = evalCDF(self,k)
            % FT = EVALCDF(SELF,K)
            
            % Evaluate the cumulative distribution function at t
            % AT T
            
            p = self.getParam(1).paramValue;
            Ft = 1 - (1-p)^k;
        end

        function L = evalLST(self, s)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform of the distribution function at s
            % For Geometric(p), LST(s) = p*e^(-s) / (1 - (1-p)*e^(-s))

            p = self.getParam(1).paramValue;
            e_neg_s = exp(-s);
            L = (p * e_neg_s) / (1 - (1 - p) * e_neg_s);
        end
        
        function pr = evalPMF(self, k)
            % PR = EVALPMF(K)

            % Evaluate the probability mass function at k
            % AT K

            p = self.getParam(1).paramValue;
            % The support is {1,2,...} (the trial index of the first success), so
            % the mass vanishes below it. Without this the formula continues
            % analytically to k=0 and returns p/(1-p), which is not a probability.
            if k < 1
                pr = 0;
                return
            end
            pr = (1-p)^(k-1)*p;
        end
        
        function proc = getProcess(self)
            % PROC = GETPROCESS()
            
            % Get process representation for non-Markovian distribution
            % Returns [mean, SCV] pair for use in network analysis
            proc = [self.getMean(), self.getSCV()];
        end
    end
    
end

