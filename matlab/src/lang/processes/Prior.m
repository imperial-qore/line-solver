classdef Prior < Distribution
    % Prior Discrete prior distribution over alternative distributions
    %
    % Prior represents parameter uncertainty by specifying a discrete set of
    % alternative distributions with associated probabilities. When used with
    % setService or setArrival, it causes the UQ solver to expand the
    % model into a family of networks, one for each alternative.
    %
    % This is NOT a mixture distribution - each alternative represents a
    % separate model realization with its associated prior probability.
    %
    % Two forms are supported:
    % - Discrete: an explicit set of alternative distributions with weights.
    % - Continuous: a density f(theta) over a scalar parameter theta plus a
    %   factory mapping theta to a Distribution. This is the form required by
    %   the epistemic uncertainty propagation of Trivedi and Bobbio (2017),
    %   Sec. 3.4, where the unconditional measure is the integral of the
    %   conditional measure against f(theta). The continuous form is reduced to
    %   a weighted alternative set by discretize(), so both forms are consumed
    %   identically downstream.
    %
    % @brief Discrete or continuous prior for parameter uncertainty modeling
    %
    % Key characteristics:
    % - Discrete set of alternative distributions, or a continuous parameter density
    % - Probability-weighted alternatives (must sum to 1)
    % - Used with UQ solver for Bayesian analysis
    %
    % Example:
    % @code
    % % Discrete form: service time with uncertain rate
    % prior = Prior({Exp(1.0), Exp(2.0), Erlang(2,1.5)}, [0.4, 0.35, 0.25]);
    % queue.setService(class, prior);
    %
    % % Continuous form: rate is itself Erlang-distributed
    % prior = Prior(Erlang(10, 3), @(lambda) Exp(lambda));
    %
    % % Continuous form from k lifetime observations summing to s
    % % (Jeffreys posterior of Trivedi-Bobbio Eq. 3.71)
    % prior = Prior.fromSample(10, 5.0);
    %
    % % Solve with UQ wrapper
    % post = UQ(model, @SolverMVA);
    % avgTable = post.getAvgTable();       % Prior-weighted expectations
    % postTable = post.getPosteriorTable(); % Per-alternative breakdown
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        distributions;  % Cell array of alternative distributions (discrete form)
        probabilities;  % Vector of prior probabilities (sum to 1) (discrete form)
        paramDist;      % Distribution over the scalar parameter (continuous form)
        distFactory;    % Handle theta -> Distribution (continuous form)
        kind;           % 'discrete' or 'continuous'
    end

    methods
        function self = Prior(varargin)
            % PRIOR Create a prior distribution instance
            %
            % @brief Creates a Prior in either discrete or continuous form
            %
            % PRIOR(DISTRIBUTIONS, PROBABILITIES) discrete form.
            % @param distributions Cell array of Distribution objects
            % @param probabilities Vector of probabilities (must sum to 1)
            %
            % PRIOR(PARAMDIST, DISTFACTORY) continuous form.
            % @param paramDist Distribution of the scalar parameter theta
            % @param distFactory Handle theta -> Distribution
            %
            % @return self Prior instance

            self@Distribution('Prior', 2, [0, Inf]);

            if nargin ~= 2
                line_error(mfilename, 'Prior requires two arguments: (distributions, probabilities) or (paramDist, distFactory)');
            end

            if isa(varargin{1}, 'Distribution') && isa(varargin{2}, 'function_handle')
                % Continuous form
                self.kind = 'continuous';
                self.paramDist = varargin{1};
                self.distFactory = varargin{2};
                self.distributions = {};
                self.probabilities = [];
                setParam(self, 1, 'paramDist', self.paramDist);
                setParam(self, 2, 'distFactory', self.distFactory);
                return;
            end

            % Discrete form
            self.kind = 'discrete';
            distributions = varargin{1};
            probabilities = varargin{2};

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

        function bool = isContinuousPrior(self)
            % BOOL = ISCONTINUOUSPRIOR()
            % Return true if the prior is specified by a parameter density.
            % Note: the base class isContinuous() refers to the support of the
            % distribution itself, not to the form of the prior.
            bool = strcmp(self.kind, 'continuous');
        end

        function [dists, weights] = discretize(self, n, method)
            % [DISTS, WEIGHTS] = DISCRETIZE(N, METHOD)
            % Reduce the prior to N weighted alternatives.
            %
            % The method is honoured for both forms of prior:
            %   'quadrature'  For a discrete prior, the alternatives and their
            %                 probabilities unchanged, and N is ignored: the
            %                 set is already exact. For a continuous prior,
            %                 stratified quantile midpoints with weights 1/N.
            %                 Each node is the conditional median of an
            %                 equal-mass stratum, so the rule integrates the
            %                 parameter density in probability space and needs
            %                 only evalCDF, which every Distribution provides.
            %   'montecarlo'  N i.i.d. draws, weights 1/N. For a discrete prior
            %                 the draws are of the alternative index against
            %                 its probabilities; returning the alternatives
            %                 unweighted here would silently drop the prior.
            %
            % @param n Number of alternatives (ignored by a discrete quadrature)
            % @param method 'quadrature' (default) or 'montecarlo'
            % @return dists Cell array of Distribution objects
            % @return weights Row vector of weights summing to 1

            if nargin < 2 || isempty(n)
                n = 11;
            end
            if nargin < 3 || isempty(method)
                method = 'quadrature';
            end
            if ~any(strcmp(method, {'quadrature', 'montecarlo'}))
                line_error(mfilename, sprintf('Unknown discretization method: %s', method));
            end

            if strcmp(self.kind, 'discrete')
                switch method
                    case 'quadrature'
                        dists = self.distributions;
                        weights = self.probabilities;
                    case 'montecarlo'
                        cumprob = cumsum(self.probabilities);
                        dists = cell(1, n);
                        for i = 1:n
                            idx = find(rand() <= cumprob, 1, 'first');
                            dists{i} = self.distributions{idx};
                        end
                        weights = ones(1, n) / n;
                end
                return;
            end

            switch method
                case 'quadrature'
                    % Midpoint of each equal-probability stratum
                    p = ((1:n) - 0.5) / n;
                    theta = zeros(1, n);
                    for i = 1:n
                        theta(i) = Prior.quantile(self.paramDist, p(i));
                    end
                case 'montecarlo'
                    s = self.paramDist.sample(n);
                    theta = s(:)';
            end

            weights = ones(1, n) / n;
            dists = cell(1, n);
            for i = 1:n
                dists{i} = self.distFactory(theta(i));
                if ~isa(dists{i}, 'Distribution')
                    line_error(mfilename, 'distFactory must return a Distribution object');
                end
            end
        end

        function n = getNumAlternatives(self)
            % N = GETNUMALTERNATIVES()
            % Return number of alternative distributions.
            % A continuous prior has no alternatives until discretize() is
            % called, so this returns NaN to force callers to discretize.
            if strcmp(self.kind, 'continuous')
                n = NaN;
                return;
            end
            n = length(self.distributions);
        end

        function dist = getAlternative(self, idx)
            % DIST = GETALTERNATIVE(IDX)
            % Return the distribution at index idx
            self.assertDiscrete('getAlternative');
            if idx < 1 || idx > self.getNumAlternatives()
                line_error(mfilename, 'Index out of bounds');
            end
            dist = self.distributions{idx};
        end

        function p = getProbability(self, idx)
            % P = GETPROBABILITY(IDX)
            % Return the probability of alternative idx
            self.assertDiscrete('getProbability');
            if idx < 1 || idx > self.getNumAlternatives()
                line_error(mfilename, 'Index out of bounds');
            end
            p = self.probabilities(idx);
        end

        function assertDiscrete(self, caller)
            % ASSERTDISCRETE(CALLER)
            % Reject enumeration of a continuous prior.
            %
            % getNumAlternatives returns NaN for a continuous prior, and every
            % comparison against NaN is false, so a bounds check alone would
            % pass and the caller would fault on an empty array instead.
            if strcmp(self.kind, 'continuous')
                line_error(mfilename, sprintf(['%s applies to a discrete prior only. ', ...
                    'A continuous prior has no alternatives until discretize() is called.'], caller));
            end
        end

        function MEAN = getMean(self)
            % MEAN = GETMEAN()
            % Get prior-weighted mean (expected mean over alternatives)
            %
            % E[X] = sum_i p_i * E[X_i]
            [dists, probs] = self.discretize();
            MEAN = 0;
            for i = 1:length(dists)
                MEAN = MEAN + probs(i) * dists{i}.getMean();
            end
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()
            % Get prior-weighted SCV using law of total variance
            %
            % Var(X) = E[Var(X|D)] + Var(E[X|D])
            % SCV = Var(X) / E[X]^2

            [dists, probs] = self.discretize();
            E_mean = 0;       % E[E[X|D]]
            E_var = 0;        % E[Var(X|D)]
            E_mean_sq = 0;    % E[E[X|D]^2]

            for i = 1:length(dists)
                m = dists{i}.getMean();
                v = dists{i}.getSCV() * m^2;  % Var(X|D=i)
                E_mean = E_mean + probs(i) * m;
                E_var = E_var + probs(i) * v;
                E_mean_sq = E_mean_sq + probs(i) * m^2;
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
            [dists, probs] = self.discretize();
            third_central = 0;
            for i = 1:length(dists)
                mi = dists{i}.getMean();
                vi = dists{i}.getVar();
                si = sqrt(vi);
                skewi = dists{i}.getSkewness();

                % E[(Xi - mu)^3] = E[(Xi - mi + mi - mu)^3]
                % Using binomial expansion
                delta = mi - mu;
                % Third central moment of Xi around its own mean
                m3i = skewi * si^3;
                % Third central moment of Xi around global mu
                m3_shifted = m3i + 3*vi*delta + delta^3;

                third_central = third_central + probs(i) * m3_shifted;
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

            % A continuous prior is sampled exactly, by drawing the parameter
            % and then the variate: discretizing first would return the law of
            % a quadrature approximation rather than of the prior itself.
            if strcmp(self.kind, 'continuous')
                theta = self.paramDist.sample(n);
                theta = theta(:)';
                testSample = self.distFactory(theta(1)).sample(1);
                X = zeros(n, numel(testSample));
                for i = 1:n
                    s = self.distFactory(theta(i)).sample(1);
                    X(i,:) = s(:)';
                end
                return;
            end

            [dists, probs] = self.discretize();

            % Determine dimensionality from first distribution
            testSample = dists{1}.sample(1);
            d = numel(testSample);

            X = zeros(n, d);

            % Generate samples
            cumprob = cumsum(probs);
            for i = 1:n
                % Select alternative based on probabilities
                r = rand();
                idx = find(r <= cumprob, 1, 'first');
                s = dists{idx}.sample(1);
                X(i,:) = s(:)';
            end
        end

        function Ft = evalCDF(self, t)
            % FT = EVALCDF(T)
            % Evaluate mixture CDF at t
            %
            % F(t) = sum_i p_i * F_i(t)

            [dists, probs] = self.discretize();
            Ft = 0;
            for i = 1:length(dists)
                Ft = Ft + probs(i) * dists{i}.evalCDF(t);
            end
        end

        function L = evalLST(self, s)
            % L = EVALLST(S)
            % Evaluate mixture Laplace-Stieltjes transform
            %
            % L(s) = sum_i p_i * L_i(s)

            [dists, probs] = self.discretize();
            L = 0;
            for i = 1:length(dists)
                L = L + probs(i) * dists{i}.evalLST(s);
            end
        end

        function bool = isPrior(self)
            % BOOL = ISPRIOR()
            % Return true (used for detection by UQ solver)
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

        function self = fromSample(k, s, distFactory)
            % SELF = FROMSAMPLE(K, S, DISTFACTORY)
            % Continuous prior for a rate estimated from lifetime data.
            %
            % Given K i.i.d. observations of an exponential random variable
            % summing to S, the Jeffreys improper prior f(lambda) = s/lambda
            % yields the posterior density of the rate
            %
            %   f(lambda|s) = lambda^(k-1) s^k exp(-lambda s) / (k-1)!
            %
            % which is an Erlang density with K phases and phase rate S.
            % See Trivedi and Bobbio (2017), Eq. (3.71). The posterior has mean
            % K/S, i.e. the maximum-likelihood rate estimate, and variance
            % K/S^2, so it concentrates on the estimate as K grows.
            %
            % @param k Number of observations (positive integer)
            % @param s Sum of the observed lifetimes (positive)
            % @param distFactory Handle theta -> Distribution, default @(lambda) Exp(lambda)
            % @return self Prior instance in continuous form

            if nargin < 3 || isempty(distFactory)
                distFactory = @(lambda) Exp(lambda);
            end
            if ~(isscalar(k) && k >= 1 && k == round(k))
                line_error(mfilename, 'k must be a positive integer number of observations');
            end
            if ~(isscalar(s) && s > 0)
                line_error(mfilename, 's must be a positive sum of observed lifetimes');
            end
            self = Prior(Erlang(s, k), distFactory);
        end

        function x = quantile(dist, p)
            % X = QUANTILE(DIST, P)
            % Numerical inverse CDF by bisection.
            %
            % Uses only evalCDF, so it applies to any Distribution. Bracketing
            % starts from the mean and doubles outward, which terminates for
            % any distribution with finite mean.
            %
            % @param dist Distribution object
            % @param p Probability level in (0,1)
            % @return x Value with F(x) = p

            if p <= 0 || p >= 1
                line_error(mfilename, 'p must lie strictly between 0 and 1');
            end

            lo = 0;
            hi = max(dist.getMean(), GlobalConstants.FineTol);
            maxExpand = 200;
            for i = 1:maxExpand
                if dist.evalCDF(hi) >= p
                    break;
                end
                hi = hi * 2;
                if i == maxExpand
                    line_error(mfilename, 'Failed to bracket the requested quantile');
                end
            end

            for i = 1:200
                mid = (lo + hi) / 2;
                if dist.evalCDF(mid) < p
                    lo = mid;
                else
                    hi = mid;
                end
                if (hi - lo) <= GlobalConstants.FineTol * max(1, hi)
                    break;
                end
            end
            x = (lo + hi) / 2;
        end
    end
end
