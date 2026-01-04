classdef Det < ContinuousDistribution & DiscreteDistribution
    % Det Deterministic distribution with constant value
    %
    % Det represents a deterministic distribution that always produces the same
    % constant value. This distribution has zero variance and is commonly used
    % for modeling constant service times, fixed delays, or deterministic
    % processing requirements in queueing systems.
    %
    % @brief Deterministic distribution with constant value and zero variance
    %
    % Key characteristics:
    % - Single parameter: constant value t
    % - Zero variance (SCV = 0)
    % - Mean = Variance = t
    % - All samples equal to t
    % - Degenerate case of all distributions
    %
    % The deterministic distribution is used for:
    % - Constant service times
    % - Fixed processing delays
    % - Deterministic timeouts
    % - Modeling perfectly predictable processes
    % - Lower bound analysis in queueing systems
    %
    % Example:
    % @code
    % constant_service = Det(2.5);     % Always takes exactly 2.5 time units
    % samples = constant_service.sample(100); % All samples = 2.5
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
    end

    methods
        function self = Det(t)
            % DET Create a deterministic distribution instance
            %
            % @brief Creates a deterministic distribution with constant value t
            % @param t Constant value (must be non-negative for time-based processes)
            % @return self Det distribution instance with constant value t
            self@ContinuousDistribution('Det',1,[t,t]);
            self@DiscreteDistribution('Det',1,[t,t]);
            setParam(self, 1, 't', t);
        end

        function ex = getMean(self)
            % EX = GETMEAN()

            % Get distribution mean
            ex = self.getParam(1).paramValue;
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()

            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            SCV = 0;
        end

        function L = evalLST(self, s)
            % L = EVALST(S)

            % Evaluate the Laplace-Stieltjes transform of the distribution function at t

            t = self.getParam(1).paramValue;

            L = exp(-s*t);
        end

        function X = sample(self, n)
            % X = SAMPLE(N)

            % Get n samples from the distribution
            X = self.getParam(1).paramValue * ones(n,1);
        end

        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)

            % Evaluate the cumulative distribution function at t
            % AT T

            if t < self.getParam(1).paramValue
                Ft = 0;
            else
                Ft = 1;
            end
        end

        function proc = getProcess(self)
            % PROC = GETPROCESS()

            % Get process representation with actual distribution parameters
            % Returns {t} for deterministic value
            proc = {self.getParam(1).paramValue};
        end
    end

    methods (Static)
        function d = fitMean(t)
            d=Det(t);
        end
    end

end

