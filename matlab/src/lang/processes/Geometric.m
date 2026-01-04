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

