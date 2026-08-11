classdef HyperExp < Markovian
    % HyperExp Hyper-exponential distribution for high-variability processes
    %
    % HyperExp represents a mixture of exponential distributions, where jobs
    % select one of multiple exponential phases with given probabilities.
    % This distribution has high variability (SCV > 1) and is commonly used
    % for modeling service times with large variation or mixed workload types.
    %
    % @brief Hyper-exponential mixture distribution with multiple phases
    %
    % Key characteristics:
    % - Mixture of multiple exponential distributions
    % - Probabilistic selection of exponential phases
    % - High variability (SCV ≥ 1)
    % - Supports both 2-phase and n-phase variants
    % - Parallel phase structure (choose one of n phases)
    %
    % The hyper-exponential distribution is used for:
    % - Service times with high variability
    % - Mixed workload modeling (fast/slow jobs)
    % - Modeling systems with multiple service classes
    % - Approximating heavy-tailed distributions
    % - Building phase-type distributions with SCV > 1
    %
    % Example:
    % @code
    % % Two-phase: 70% fast (rate=5), 30% slow (rate=0.5)
    % mixed_service = HyperExp(0.7, 5.0, 0.5);
    % % n-phase: equal probability, different rates
    % multi_service = HyperExp([0.3, 0.4, 0.3], [3.0, 1.0, 0.2]);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = HyperExp(varargin)
            % HYPEREXP Create a hyper-exponential distribution instance
            %
            % @brief Creates a hyper-exponential distribution with specified phases
            % @param varargin Variable arguments: (p, lambda1, lambda2) or (prob_vec, rate_vec)
            % @return self HyperExp distribution instance
            %
            % Usage: HyperExp(p1, lambda1, lambda2) for 2-phase
            %        HyperExp(prob_vector, rate_vector) for n-phase
            self@Markovian('HyperExp',nargin);
            if length(varargin)==2
                p = varargin{1};
                lambda = varargin{2};
                if isscalar(p)
                    % Two phases sharing one rate: HyperExp(p, lambda).
                    setParam(self, 1, 'p', p);
                    setParam(self, 2, 'lambda1', lambda);
                    setParam(self, 3, 'lambda2', lambda);
                    self.obj = jline.lang.processes.HyperExp(p, lambda, lambda);
                else
                    % n phases: p and lambda are equal-length vectors of the
                    % phase probabilities and phase rates.
                    %
                    % param2 holds the whole rate vector and param3 is left
                    % empty. Assigning the full vector to BOTH, as this did,
                    % made linemodel_save emit a lambda of length 2n, and let
                    % the moment methods below read a rate vector where they
                    % expect a scalar. No jline object is built: the Java
                    % HyperExp constructors take scalars, and a native APH
                    % likewise leaves obj unset.
                    p = p(:).';
                    lambda = lambda(:).';
                    if numel(p) ~= numel(lambda)
                        line_error(mfilename, sprintf(['HyperExp: p and lambda must have ', ...
                            'the same length (got %d and %d).'], numel(p), numel(lambda)));
                    end
                    setParam(self, 1, 'p', p);
                    setParam(self, 2, 'lambda1', lambda);
                    setParam(self, 3, 'lambda2', []);
                end
            elseif length(varargin)==3
                p1 = varargin{1};
                lambda1 = varargin{2};
                lambda2 = varargin{3};
                setParam(self, 1, 'p', p1);
                setParam(self, 2, 'lambda1', lambda1);
                setParam(self, 3, 'lambda2', lambda2);
                self.obj = jline.lang.processes.HyperExp(p1, lambda1, lambda2);
            end
            p = self.getParam(1).paramValue;
            n = length(p);
            if n == 1
                mu1 = self.getParam(2).paramValue;
                mu2 = self.getParam(3).paramValue;
                PH={[-mu1,0;0,-mu2],[mu1*p,mu1*(1-p);mu2*p,mu2*(1-p)]};
            else
                mu = self.getParam(2).paramValue;
                D0 = -diag(mu);
                % After completing phase i (at rate mu(i)) the process restarts
                % in phase j with probability p(j), so D1(i,j) = mu(i)*p(j).
                %
                % This is `(-D0*ones(n,1)) * p(:).'`, an outer product of the
                % exit rates with the entry vector. Writing it as
                % `-D0*p(:)*ones(1,n)` instead associates to the LEFT, giving
                % ((-D0)*p(:))*ones(1,n), i.e. D1(i,j) = mu(i)*p(i) replicated
                % across every column. That is not a generator: the rows of
                % D0+D1 then sum to mu(i)*(n*p(i)-1) rather than to 0, and any
                % consumer recovering the entry vector from D1 (JMT export via
                % sn.pie/sn.impatiencePie, MAM) read back p(i)*ones(1,n)
                % instead of p. It matched only for uniform p, where n*p(i)=1.
                D1 = (-D0*ones(n,1))*p(:).';
                PH = {D0,D1};
            end
            self.process = PH;
            self.nPhases = length(self.process{1});
        end

        function ex = getMean(self)
            % EX = GETMEAN()
            % Get distribution mean
            %
            % The branch tests the SCALAR form, not getNumberOfPhases==2: a
            % two-phase vector form (HyperExp([p1 p2],[l1 l2])) also has two
            % phases, and reading its vector parameters as scalars silently
            % produced a vector "mean". The n-phase closed form is exact, so
            % neither branch needs the base-class fallback, which was called
            % without assigning the output and therefore always errored.
            p = self.getParam(1).paramValue;
            if isscalar(p)
                mu1 = self.getParam(2).paramValue;
                mu2 = self.getParam(3).paramValue;
                ex = p/mu1 + (1-p)/mu2;
            else
                mu = self.getParam(2).paramValue;
                ex = sum(p(:)./mu(:));
            end
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()
            % Get the squared coefficient of variation of the distribution (SCV = variance / mean^2)
            p = self.getParam(1).paramValue;
            if isscalar(p)
                mu1 = self.getParam(2).paramValue;
                mu2 = self.getParam(3).paramValue;
                SCV = (2*(p/mu1^2 + (1-p)/mu2^2) - (p/mu1 + (1-p)/mu2)^2)/(p/mu1 + (1-p)/mu2)^2;
            else
                % E[X]=sum(p_i/mu_i), E[X^2]=sum(2 p_i/mu_i^2), SCV=E[X^2]/E[X]^2-1
                mu = self.getParam(2).paramValue;
                m1 = sum(p(:)./mu(:));
                m2 = sum(2*p(:)./(mu(:).^2));
                SCV = m2/m1^2 - 1;
            end
        end

        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            % Evaluate the cumulative distribution function at t
            % AT T

            p = self.getParam(1).paramValue;
            if isscalar(p)
                mu1 = self.getParam(2).paramValue;
                mu2 = self.getParam(3).paramValue;
                Ft = p*(1-exp(-mu1*t))+(1-p)*(1-exp(-mu2*t));
            else
                % Mixture of exponentials: F(t) = sum_i p_i (1 - exp(-mu_i t))
                mu = self.getParam(2).paramValue;
                Ft = zeros(size(t));
                for i = 1:numel(p)
                    Ft = Ft + p(i)*(1-exp(-mu(i)*t));
                end
            end
        end
    end

    methods(Static)
        function he = fit(MEAN, SCV, SKEW)
            % HE = FIT(MEAN, SCV, SKEW)
            % Fit distribution from first three standard moments
            e1 = MEAN;
            e2 = (1+SCV)*e1^2;
            e3 = -(2*e1^3-3*e1*e2-SKEW*(e2-e1^2)^(3/2));
            for n=[2] % larger ph requires more moments
                PH = kpcfit_ph_prony([e1,e2,e3],n);
                if map_isfeasible(PH)
                    p = map_pie(PH);
                    lambda = -diag(PH{1});
                    he = HyperExp(p,lambda);
                    return
                end
            end
            he = HyperExp.fitMeanAndSCV(MEAN,SCV);
            he.immediate = MEAN < GlobalConstants.CoarseTol;
        end

        function he = fitRate(RATE)
            % HE = FITRATE(RATE)
            % Fit distribution with given rate
            % both phases share the rate, so the mixing probability is
            % immaterial; p was never defined here and every call errored
            p = 0.5;
            he = HyperExp(p, RATE, RATE);
            he.immediate = 1/RATE < GlobalConstants.CoarseTol;
        end

        function he = fitMean(MEAN)
            % HE = FITMEAN(MEAN)
            % Fit distribution with given mean
            % as in fitRate: equal rates, so p is immaterial but must be set
            p = 0.5;
            he = HyperExp(p, 1/MEAN, 1/MEAN);
            he.immediate = MEAN < GlobalConstants.CoarseTol;
        end

        function he = fitMeanAndSCV(MEAN, SCV)
            % HE = FITMEANANDSCV(MEAN, SCV)
            % Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
            [~,mu1,mu2,p] = map_hyperexp(MEAN,SCV);
            he = HyperExp(p, mu1, mu2);
            he.immediate = MEAN < GlobalConstants.CoarseTol;
        end

        function he = fitMeanAndSCVBalanced(MEAN, SCV)
            % HE = FITMEANANDSCV(MEAN, SCV)
            % Fit distribution with given mean and squared coefficient of
            % variation (SCV=variance/mean^2) and balanced means, i.e.,
            % p/mu1 = (1-p)/mu2
            mu1 =  -(2*(((SCV - 1)/(SCV + 1))^(1/2)/2 - 1/2))/MEAN;
            p= 1/2 - ((SCV - 1)/(SCV + 1))^(1/2)/2;
            if mu1<0 || p<0 || p>1
                p = ((SCV - 1)/(SCV + 1))^(1/2)/2 + 1/2;
                mu1 = (2*(((SCV - 1)/(SCV + 1))^(1/2)/2 + 1/2))/MEAN;
            end
            mu2=(1-p)/p*mu1;
            p = real(p);
            mu1 = real(mu1);
            mu2 = real(mu2);
            he = HyperExp(p, mu1, mu2);
            he.immediate = MEAN < GlobalConstants.CoarseTol;
        end
    end

end

