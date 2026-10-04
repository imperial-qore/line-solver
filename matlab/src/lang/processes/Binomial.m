classdef Binomial < DiscreteDistribution
    % A Binomial distribution
    %
    % Copyright (c) 2018-2022, Imperial College London
    % All rights reserved.

    methods
        function self = Binomial(n, p)
            % SELF = BINOMIAL(N, P)

            self@DiscreteDistribution('Binomial',2,[0,Inf]);
            setParam(self, 1, 'n', n);
            setParam(self, 2, 'p', p);
        end

        function ex = getMean(self)
            % EX = GETMEAN()

            % Get distribution mean
            n = self.getParam(1).paramValue;
            p = self.getParam(2).paramValue;
            ex = p*n;
        end

        function V = getVar(self)
            % V = GETVAR()
            n = self.getParam(1).paramValue;
            p = self.getParam(2).paramValue;
            V = n*p*(1-p);
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()

            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            p = self.getParam(2).paramValue;
            n = self.getParam(1).paramValue;
            SCV = (1-p) / (n*p);
        end

        function X = sample(self, nsamples)
            % X = SAMPLE(N)

            X = zeros(1, nsamples);
            p = self.getParam(2).paramValue;
            n = self.getParam(1).paramValue;

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

            p = self.getParam(2).paramValue;
            n = self.getParam(1).paramValue;
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
            % For Binomial(n, p), LST(s) = (1 - p + p*e^(-s))^n

            p = self.getParam(2).paramValue;
            n = self.getParam(1).paramValue;
            L = (1 - p + p * exp(-s))^n;
        end

        function P = evalPMF(self, k)
            % P = EVALPMF(K)

            % Evaluate the probability mass function at k            

            p = self.getParam(2).paramValue;
            n = self.getParam(1).paramValue;
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

