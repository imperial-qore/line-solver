classdef Immediate < Distribution
    % A distribution with probability mass entirely at zero
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods (Hidden)
        %Constructor
        function self = Immediate()
            % SELF = IMMEDIATE()

            self@Distribution('Immediate', 0,[0,0]);
            self.immediate = true;
            self.obj = jline.lang.processes.Immediate();            
        end
    end

    methods
        function bool = isDisabled(self)
            % BOOL = ISDISABLED()

            bool = false;
        end

        function X = sample(self, n)
            % X = SAMPLE(N)

            if nargin<2 %~exist('n','var'),
                n = 1;
            end
            X = zeros(n,1);
        end

        function rate = getRate(self)
            %global GlobalConstants.Immediate
            rate = GlobalConstants.Immediate;
        end

        function ex = getMean(self)
            % EX = GETMEAN()

            % Get distribution mean
            ex = 0;
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()

            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            % Return 1 to match native Python behavior which produces correct LQNS results

            SCV = 1;
        end

        function mu = getMu(self)
            % MU = GETMU()
            %global GlobalConstants.Immediate
            % Return total outgoing rate from each state
            mu = GlobalConstants.Immediate;
        end

        function phi = getPhi(self)
            % PHI = GETPHI()

            % Return the probability that a transition out of a state is
            % absorbing
            phi = 1.0;
        end

        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)

            % Evaluate the cumulative distribution function at t
            % AT T

            Ft = 1;
        end

        function L = evalLST(self, s)
            % L = EVALST(S)

            % Evaluate the Laplace-Stieltjes transform of the distribution function at t

            L = 1; % as in Det(0)
        end

        function PH = getProcess(self)
            PH = {[-GlobalConstants.Immediate] ,[GlobalConstants.Immediate]};
        end

        function bool = isImmediate(self)
            % BOOL = ISIMMEDIATE()
            % Check if the distribution is equivalent to an Immediate
            % distribution
            % Overrides Distribution.isImmediate(self)
            bool = true;
        end

    end

    methods (Static)
        function [singleton, javasingleton] = getInstance
            persistent staticImmediate
            if isempty(staticImmediate)
                staticImmediate = Immediate();
            end
            singleton = staticImmediate;
            javasingleton = staticImmediate.obj;
        end
    end
end