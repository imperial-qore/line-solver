classdef Zipf < DiscreteDistribution
    % A Zipf-like popularity distribution
    %
    % Copyright (c) 2018-2022, Imperial College London
    % All rights reserved.
    
    methods
        function self = Zipf(s, n)
            % SELF = ZIPF(S, N)
            
            % Construct a Zipf-like distribution on n items with given
            % shape parameter s
            if nargin<2 %~exist('n','var')
                n = GlobalConstants.Immediate;
            end
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
            
            % Get n samples from the distribution
            X = self.getParam(3).paramValue * ones(n,1);
        end
        
        function Ft = evalCDF(self,k)
            % FT = EVALCDF(SELF,K)
            
            % Evaluate the cumulative distribution function at t
            % AT T
            
            s = self.getParam(3).paramValue;
            n = self.getParam(4).paramValue;
            Ft = self.genHarmonic(s,k) / self.genHarmonic(s,n);
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

