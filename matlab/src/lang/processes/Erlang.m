classdef Erlang < Markovian
    % Erlang Multi-phase exponential distribution with reduced variability
    %
    % Erlang represents the distribution of the sum of r independent exponential
    % random variables with the same rate. It models processes with reduced
    % variability compared to exponential distributions and is commonly used
    % for service times with more predictable durations.
    %
    % @brief Erlang distribution with phase rate and number of phases
    %
    % Key characteristics:
    % - Sum of r independent exponential phases
    % - Two parameters: phase rate (alpha) and number of phases (r)
    % - Mean = r/alpha, Variance = r/alpha²
    % - SCV = 1/r (always ≤ 1, decreasing with r)
    % - Multi-phase Markovian representation
    %
    % The Erlang distribution is used for:
    % - Service times with reduced variability
    % - Multi-stage service processes
    % - Modeling deterministic processes approximately
    % - Building blocks for phase-type distributions
    % - Systems where variability decreases with stages
    %
    % Example:
    % @code
    % service_dist = Erlang(3.0, 4);    % 4 phases, rate 3.0 each
    % % Mean = 4/3.0 = 1.33, SCV = 1/4 = 0.25
    % samples = service_dist.sample(1000);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = Erlang(phaseRate, nphases)
            % ERLANG Create an Erlang distribution instance
            %
            % @brief Creates an Erlang distribution with specified phase rate and phases
            % @param phaseRate Rate parameter for each exponential phase  
            % @param nphases Number of sequential exponential phases (positive integer)
            % @return self Erlang distribution instance
            self@Markovian('Erlang',2);
            setParam(self, 1, 'alpha', phaseRate); % rate in each state
            setParam(self, 2, 'r', round(nphases)); % number of phases
            self.obj = jline.lang.processes.Erlang(phaseRate, nphases);
            r = self.getParam(2).paramValue;
            self.process = map_erlang(getMean(self),r);
            self.nPhases = length(self.process{1});
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

        function setMean(self,MEAN)
            % UPDATEMEAN(SELF,MEAN)
            % Update parameters to match the given mean
            setMean@Markovian(self,MEAN);
            r = self.getParam(2).paramValue;
            self.setParam(1, 'alpha', r/MEAN);
            self.process = map_erlang(getMean(self),r);
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
            er.immediate = MEAN < GlobalConstants.FineTol;
        end

        % function er = fitRate(RATE)
        %     % ER = FITRATE(RATE)
        % 
        %     % Fit distribution with given rate
        %     line_warning(mfilename,'The Erlang distribution is underspecified by the rate, setting the number of phases to 2.\n');
        %     er = Erlang.fitMeanAndOrder(1/RATE, 2);
        %     er.immediate = 1/RATE < GlobalConstants.FineTol;
        % end

        % function er = fitMean(MEAN)
        %     % ER = FITMEAN(MEAN)
        % 
        %     % Fit distribution with given mean
        %     line_warning(mfilename,'The Erlang distribution is underspecified by the mean, setting the number of phases to 2.\n');
        %     er = Erlang.fitMeanAndOrder(MEAN, 2);
        %     er.immediate = MEAN < GlobalConstants.FineTol;
        % end

        function er = fitMeanAndSCV(MEAN, SCV)
            % ER = FITMEANANDSCV(MEAN, SCV)
            if SCV>1
                line_error(mfilename,'The Erlang distribution requires squared coefficient of vairation <= 1.');
            end
            % Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
            r = ceil(1/SCV);
            alpha = r/MEAN;
            er = Erlang(alpha, r);
            er.immediate = MEAN < GlobalConstants.FineTol;
        end

        function er = fitMeanAndOrder(MEAN, n)
            % ER = FITMEANANDORDER(MEAN, N)

            % Fit distribution with given mean and number of phases
            SCV = 1/n;
            r = ceil(1/SCV);
            alpha = r/MEAN;
            er = Erlang(alpha, r);
            er.immediate = MEAN < GlobalConstants.FineTol;
        end
    end

end