classdef Prior < Distribution
    % Prior Discrete prior distribution over alternative distributions
    %
    % Prior represents parameter uncertainty by specifying a discrete set of
    % alternative distributions with associated probabilities. When used with
    % setService or setArrival, it causes the Posterior solver to expand the
    % model into a family of networks, one for each alternative.
    %
    % This is NOT a mixture distribution - each alternative represents a
    % separate model realization with its associated prior probability.
    %
    % @brief Discrete prior for Bayesian-style parameter uncertainty modeling
    %
    % Key characteristics:
    % - Discrete set of alternative distributions
    % - Probability-weighted alternatives (must sum to 1)
    % - Used with Posterior solver for Bayesian analysis
    % - Supports single Prior per model (current limitation)
    %
    % Example:
    % @code
    % % Service time with uncertain rate
    % prior = Prior({Exp(1.0), Exp(2.0), Erlang(2,1.5)}, [0.4, 0.35, 0.25]);
    % queue.setService(class, prior);
    %
    % % Solve with Posterior wrapper
    % post = Posterior(model, @SolverMVA);
    % avgTable = post.getAvgTable();       % Prior-weighted expectations
    % postTable = post.getPosteriorTable(); % Per-alternative breakdown
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        distributions;  % Cell array of alternative distributions
        probabilities;  % Vector of prior probabilities (sum to 1)
    end

    methods
        function self = Prior(distributions, probabilities)
            % PRIOR Create a prior distribution instance
            %
            % @brief Creates a Prior with specified alternatives and probabilities
            % @param distributions Cell array of Distribution objects
            % @param probabilities Vector of probabilities (must sum to 1)
            % @return self Prior instance

            self@Distribution('Prior', 2, [0, Inf]);

            % Validate distributions input
            if ~iscell(distributions)
                line_error(mfilename, 'distributions must be a cell array');
            end
            if isempty(distributions)
                line_error(mfilename, 'distributions cannot be empty');
            end
            for i = 1:length(distributions)
                if ~isa(distributions{i}, 'Distribution')
                    line_error(mfilename, sprintf('Element %d is not a Distribution object', i));
                end
            end

            % Validate probabilities input
            if length(distributions) ~= length(probabilities)
                line_error(mfilename, 'Number of distributions must match number of probabilities');
            end
            if abs(sum(probabilities) - 1) > GlobalConstants.CoarseTol
                line_error(mfilename, sprintf('Probabilities must sum to 1 (current sum: %f)', sum(probabilities)));
            end
            if any(probabilities < 0)
                line_error(mfilename, 'Probabilities must be non-negative');
            end

            self.distributions = distributions(:)';  % row cell array
            self.probabilities = probabilities(:)';  % row vector

            setParam(self, 1, 'distributions', distributions);
            setParam(self, 2, 'probabilities', probabilities);
        end

        function n = getNumAlternatives(self)
            % N = GETNUMALTERNATIVES()
            % Return number of alternative distributions
            n = length(self.distributions);
        end

        function dist = getAlternative(self, idx)
            % DIST = GETALTERNATIVE(IDX)
            % Return the distribution at index idx
            if idx < 1 || idx > self.getNumAlternatives()
                line_error(mfilename, 'Index out of bounds');
            end
            dist = self.distributions{idx};
        end

        function p = getProbability(self, idx)
            % P = GETPROBABILITY(IDX)
            % Return the probability of alternative idx
            if idx < 1 || idx > self.getNumAlternatives()
                line_error(mfilename, 'Index out of bounds');
            end
            p = self.probabilities(idx);
        end

        function MEAN = getMean(self)
            % MEAN = GETMEAN()
            % Get prior-weighted mean (expected mean over alternatives)
            %
            % E[X] = sum_i p_i * E[X_i]
            MEAN = 0;
            for i = 1:self.getNumAlternatives()
                MEAN = MEAN + self.probabilities(i) * self.distributions{i}.getMean();
            end
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()
            % Get prior-weighted SCV using law of total variance
            %
            % Var(X) = E[Var(X|D)] + Var(E[X|D])
            % SCV = Var(X) / E[X]^2

            E_mean = 0;       % E[E[X|D]]
            E_var = 0;        % E[Var(X|D)]
            E_mean_sq = 0;    % E[E[X|D]^2]

            for i = 1:self.getNumAlternatives()
                m = self.distributions{i}.getMean();
                v = self.distributions{i}.getSCV() * m^2;  % Var(X|D=i)
                E_mean = E_mean + self.probabilities(i) * m;
                E_var = E_var + self.probabilities(i) * v;
                E_mean_sq = E_mean_sq + self.probabilities(i) * m^2;
            end

            % Total variance = E[Var(X|D)] + Var(E[X|D])
            % Var(E[X|D]) = E[E[X|D]^2] - E[E[X|D]]^2
            total_var = E_var + (E_mean_sq - E_mean^2);
            SCV = total_var / E_mean^2;
        end

        function SKEW = getSkewness(self)
            % SKEW = GETSKEWNESS()
            % Get prior-weighted skewness (approximation using mixture formula)

            % For mixture: use law of total cumulance (simplified)
            % This is an approximation - exact formula is more complex
            mu = self.getMean();
            sigma2 = self.getVar();
            sigma = sqrt(sigma2);

            if sigma < GlobalConstants.FineTol
                SKEW = 0;
                return;
            end

            % E[(X - mu)^3] via mixture
            third_central = 0;
            for i = 1:self.getNumAlternatives()
                mi = self.distributions{i}.getMean();
                vi = self.distributions{i}.getVar();
                si = sqrt(vi);
                skewi = self.distributions{i}.getSkewness();

                % E[(Xi - mu)^3] = E[(Xi - mi + mi - mu)^3]
                % Using binomial expansion
                delta = mi - mu;
                % Third central moment of Xi around its own mean
                m3i = skewi * si^3;
                % Third central moment of Xi around global mu
                m3_shifted = m3i + 3*vi*delta + delta^3;

                third_central = third_central + self.probabilities(i) * m3_shifted;
            end

            SKEW = third_central / sigma^3;
        end

        function X = sample(self, n)
            % X = SAMPLE(N)
            % Sample from prior (mixture sampling)
            %
            % Samples are drawn from the mixture distribution where each
            % sample comes from one of the alternatives selected according
            % to the prior probabilities.

            if nargin < 2
                n = 1;
            end
            X = zeros(n, 1);

            % Generate samples
            cumprob = cumsum(self.probabilities);
            for i = 1:n
                % Select alternative based on probabilities
                r = rand();
                idx = find(r <= cumprob, 1, 'first');
                X(i) = self.distributions{idx}.sample(1);
            end
        end

        function Ft = evalCDF(self, t)
            % FT = EVALCDF(T)
            % Evaluate mixture CDF at t
            %
            % F(t) = sum_i p_i * F_i(t)

            Ft = 0;
            for i = 1:self.getNumAlternatives()
                Ft = Ft + self.probabilities(i) * self.distributions{i}.evalCDF(t);
            end
        end

        function L = evalLST(self, s)
            % L = EVALLST(S)
            % Evaluate mixture Laplace-Stieltjes transform
            %
            % L(s) = sum_i p_i * L_i(s)

            L = 0;
            for i = 1:self.getNumAlternatives()
                L = L + self.probabilities(i) * self.distributions{i}.evalLST(s);
            end
        end

        function bool = isPrior(self)
            % BOOL = ISPRIOR()
            % Return true (used for detection by Posterior solver)
            bool = true;
        end
    end

    methods (Static)
        function bool = isPriorDistribution(dist)
            % BOOL = ISPRIORDISTRIBUTION(DIST)
            % Check if a distribution is a Prior
            %
            % @param dist Distribution object to check
            % @return bool True if dist is a Prior
            bool = isa(dist, 'Prior');
        end
    end
end
