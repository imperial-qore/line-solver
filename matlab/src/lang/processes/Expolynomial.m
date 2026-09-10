classdef Expolynomial < ContinuousDistribution
    % Expolynomial distribution with density f(x) = sum ci * x^ai * exp(-li*x)
    %
    % Represents an expolynomial density over a bounded domain [eft, lft],
    % matching the GEN expolynomial format of external stochastic Petri net tools.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = Expolynomial(density, eft, lft)
            % EXPOLYNOMIAL Create an Expolynomial distribution instance
            %
            % @param density Density expression string in expolynomial (GEN) format
            % @param eft Earliest firing time (lower bound of support)
            % @param lft Latest firing time (upper bound of support, can be Inf)
            % @return self Expolynomial distribution instance
            self@ContinuousDistribution('Expolynomial', 3, [NaN, NaN, NaN]);
            setParam(self, 1, 'density', density);
            setParam(self, 2, 'eft', eft);
            setParam(self, 3, 'lft', lft);
        end

        function ex = getMean(self)
            % EX = GETMEAN()

            % Get distribution mean (returns NaN - numerical integration not supported)
            ex = NaN;
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()

            % Get distribution SCV (returns NaN - numerical integration not supported)
            SCV = NaN;
        end

        function Ft = evalCDF(self, t)
            % FT = EVALCDF(SELF,T)

            % Evaluate the CDF at t (returns NaN - not supported)
            Ft = NaN;
        end

        function X = sample(self, n)
            % X = SAMPLE(N)

            % Get n samples from the distribution (returns NaN - not supported)
            if nargin < 2
                n = 1;
            end
            X = NaN(n, 1);
        end

        function proc = getProcess(self)
            % PROC = GETPROCESS()

            % Get process representation
            proc = {self.getParam(1).paramValue, self.getParam(2).paramValue, self.getParam(3).paramValue};
        end
    end

end
