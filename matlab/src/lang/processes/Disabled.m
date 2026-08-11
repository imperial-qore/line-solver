classdef Disabled < ContinuousDistribution & DiscreteDistribution
    % A distribution that is not configured
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = Disabled()
            % SELF = DISABLED()
            
            % Constructs a disabled distribution
            self@ContinuousDistribution('Disabled',1,[NaN,NaN]);
            self@DiscreteDistribution('Disabled',1,[NaN,NaN]);
            setParam(self, 1, 'value', NaN);
            % not implemented yet
            self.obj = jline.lang.processes.Disabled();
        end
    end
    
    methods
        function bool = isContinuous(self)
            % BOOL = ISCONTINUOUS()
            
            % Returns true is the distribution is continuous
            bool = true;
        end
        
        function bool = isDisabled(self)
            % BOOL = ISDISABLED()
            
            % Returns true is the distribution is disabled
            bool = true;
        end
        
        function bool = isDiscrete(self)
            % BOOL = ISDISCRETE()
            
            % Returns true is the distribution is discrete
            bool = true;
        end
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            
            % Get n samples from the distribution
            if nargin<2 %~exist('n','var'), 
                n = 1; 
            end
            X = nan(n,1);
        end
        
        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Get distribution mean
            ex = NaN;
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            SCV = NaN;
        end
        
        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            
            % Evaluate the cumulative distribution function at t
            % AT T
            
            Ft = NaN;
        end
        
        function p = evalPMF(self, k)
            % P = EVALPMF(K)
            
            % Evaluate the probability mass function at k
            % AT K
            
            p = 0*k;
        end
    end
    
    methods (Static)
        function singleton = getInstance
            persistent staticDisabled
            if isempty(staticDisabled)
                staticDisabled = Disabled();
            end
            singleton = staticDisabled;
        end
    end
    
end

