classdef Erlang < MarkovianDistribution
    % The Erlang statistical distribution
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    methods
        function self = Erlang(phaseRate, nphases)
            % SELF = ERLANG(PHASERATE, NPHASES)

            % Constructs an erlang distribution from the rate in each state
            % and the number of phases
            self@MarkovianDistribution('Erlang',2);
            setParam(self, 1, 'alpha', phaseRate); % rate in each state
            setParam(self, 2, 'r', round(nphases)); % number of phases
            self.obj = jline.lang.distributions.Erlang(phaseRate, nphases);
        end

        function phases = getNumberOfPhases(self)
            % PHASES = GETNUMBEROFPHASES()

            % Get number of phases in the underpinnning phase-type
            % representation
            phases  = self.getParam(2).paramValue; %r
        end

        function ex = getMean(self)
            % EX = GETMEAN()

            % Get distribution mean
            alpha = self.getParam(1).paramValue;
            r = self.getParam(2).paramValue;
            ex = r/alpha;
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()
            % Get the squared coefficient of variation of the distribution (SCV = variance / mean^2)
            r = self.getParam(2).paramValue;
            SCV = 1/r;
        end

        function updateMean(self,MEAN)
            % UPDATEMEAN(SELF,MEAN)
            % Update parameters to match the given mean
            updateMean@MarkovianDistribution(self,MEAN);
            r = self.getParam(2).paramValue;
            self.setParam(1, 'alpha', r/MEAN);
        end


        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)

            % Evaluate the cumulative distribution function at t
            % AT T

            alpha = self.getParam(1).paramValue; % rate
            r = self.getParam(2).paramValue; % stages
            Ft = 1;
            for j=0:(r-1)
                Ft = Ft - exp(-alpha*t).*(alpha*t).^j/factorial(j);
            end
        end

        function PH = getPH(self)
            % PH = GETREPRESENTATION()

            % Return the renewal process associated to the distribution
            r = self.getParam(2).paramValue;
            PH = map_erlang(getMean(self),r);
        end

        function L = evalLST(self, s)
            % L = EVALLAPLACETRANSFORM(S)

            % Evaluate the Laplace transform of the distribution function at t
            % AT T

            alpha = self.getParam(1).paramValue; % rate
            r = self.getParam(2).paramValue; % stages
            L = (alpha / (alpha + s))^r;
        end
    end

    methods(Static)
        function er = fit(MEAN, SCV, SKEW)
            % ER = FITCENTRAL(MEAN, SCV, SKEW)

            % Fit distribution from first three central moments (mean,
            % SCV, skewness)
            er = Erlang.fitMeanAndSCV(MEAN,SCV);
            er.immediate = MEAN < GlobalConstants.CoarseTol;
        end

        function er = fitRate(RATE)
            % ER = FITRATE(RATE)

            % Fit distribution with given rate
            line_warning(mfilename,'The Erlang distribution is underspecified by the rate, setting the number of phases to 2.\n');
            er = Erlang.fitMeanAndOrder(1/RATE, 2);
            er.immediate = 1/RATE < GlobalConstants.CoarseTol;
        end

        function er = fitMean(MEAN)
            % ER = FITMEAN(MEAN)

            % Fit distribution with given mean
            line_warning(mfilename,'The Erlang distribution is underspecified by the mean, setting the number of phases to 2.\n');
            er = Erlang.fitMeanAndOrder(MEAN, 2);
            er.immediate = MEAN < GlobalConstants.CoarseTol;
        end

        function er = fitMeanAndSCV(MEAN, SCV)
            % ER = FITMEANANDSCV(MEAN, SCV)
            if SCV>1
                line_error(mfilename,'The Erlang distribution requires squared coefficient of vairation <= 1.');
            end
            % Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
            r = ceil(1/SCV);
            alpha = r/MEAN;
            er = Erlang(alpha, r);
            er.immediate = MEAN < GlobalConstants.CoarseTol;
        end

        function er = fitMeanAndOrder(MEAN, n)
            % ER = FITMEANANDORDER(MEAN, N)

            % Fit distribution with given mean and number of phases
            SCV = 1/n;
            r = ceil(1/SCV);
            alpha = r/MEAN;
            er = Erlang(alpha, r);
            er.immediate = MEAN < GlobalConstants.CoarseTol;
        end
    end

end