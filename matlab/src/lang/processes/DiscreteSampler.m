classdef DiscreteSampler < DiscreteDistribution
    % A class for discrete distributions specified from the probability mass function
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = DiscreteSampler(p,x)
            % SELF = DISCRETESAMPLER(P,X)
            
            % Construct a discrete distribution from a finite probability
            % vector p at the points specificied in vector x
            
            % p(i) is the probability of item i
            % x(i) is the value of item i
            n = length(p);
            if nargin<2 %~exist('x','var')
                x=1:n;
            end
            self@DiscreteDistribution('DiscreteSampler',3,[min(x),max(x)]);
            setParam(self, 1, 'p', p(:)');
            setParam(self, 2, 'x', x(:)');
            setParam(self, 3, 'f', cumsum(p(:)')/sum(p));
        end
        
        function ex = getMean(self)
            % EX = GETMEAN()

            % Get distribution mean.
            % Weighted by the SUPPORT x, not by the rank 1:n. The two agree on
            % the default support x=1:n, which is why the rank form went
            % unnoticed; on any other support it returned the mean of a
            % distribution the caller never asked for.
            p = self.getParam(1).paramValue;
            x = self.getParam(2).paramValue;
            ex = sum(p(:)'.*x(:)') / sum(p);
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()

            % Get distribution squared coefficient of variation (SCV = variance / mean^2)

            p = self.getParam(1).paramValue;
            x = self.getParam(2).paramValue;
            % `ex` was never assigned here, so every call errored out.
            ex = getMean(self);
            e2 = sum(p(:)'.*(x(:)'.^2)) / sum(p);
            var = e2 - ex^2;
            SCV = var / ex^2;
        end
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            
            x = self.getParam(2).paramValue;
            f = self.getParam(3).paramValue;
            r = rand(n,1);
            X = x(maxpos(min(repmat(r,1,size(f,2))<=f,2)')');
        end
        
        function Ft = evalCDF(self,k)
            % FT = EVALCDF(SELF,K)
            
            f = self.getParam(3).paramValue;
            if k>=1 && k<length(f)
                Ft = f(k);
            else
                Ft = 0;
            end
        end
        
        function pk = evalPMF(self, v)
            % PK = EVALPMF(V)
            
            p = self.getParam(1).paramValue;
            x = self.getParam(2).paramValue;
            if nargin<2 %~exist('v','var')
                pk = p;
                return
            end
            pk = zeros(size(v)); % a point outside the support carries no mass
            for i=1:length(v)
                idx = find(x==v(i),1);
                if ~isempty(idx)
                    pk(i) = p(idx);
                end
            end
        end
        
        function bool = isDisabled(self)
            % BOOL = ISDISABLED()
            
            p = self.getParam(1).paramValue;
            bool = any(isnan(p));
        end
    end
    
end

