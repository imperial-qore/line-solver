classdef ME < Markovian
    % Matrix Exponential (ME) distribution
    %
    % ME distributions are characterized by an initial vector alpha and
    % a matrix parameter A. They generalize Phase-Type (PH) distributions
    % by allowing alpha to have entries outside [0,1] and A to have
    % arbitrary structure (not necessarily a valid sub-generator).
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        alpha;  % Initial vector
        A;      % Matrix parameter
    end

    methods
        function self = ME(alpha, A, checkDensity)
            % ME Create a Matrix Exponential distribution instance
            %
            % @brief Creates an ME distribution with the given initial vector and matrix parameter
            % @param alpha Initial vector (may have negative entries or sum != 1)
            % @param A Matrix parameter (must have all eigenvalues with negative real parts)
            % @param checkDensity Scan the density for a negative value (default true).
            %        Set to false only by subclasses whose representation is a density by
            %        construction, such as CME, where the scan cannot fire and costs O(1e5)
            %        propagations of a large matrix.
            % @return self ME distribution instance

            if nargin < 3 || isempty(checkDensity)
                checkDensity = true;
            end

            % Call superclass constructor
            self@Markovian('ME', 2);

            % Validate using BuTools
            if ~CheckMERepresentation(alpha, A)
                error('Invalid ME representation: Check that A is square, alpha and A have compatible dimensions, all eigenvalues of A have negative real parts, and the dominant eigenvalue is real.');
            end

            % A valid ME representation may still have a density that goes
            % negative, which is not a distribution, so the density is scanned
            % for an actual negative value. Only a witness is reported: the scan
            % warns when it has found a point where f(t) < 0, and stays silent
            % otherwise, since no cheap test establishes the converse. This is a
            % warning and not an error, so that a representation whose density
            % only dips below zero at the level of round-off remains
            % constructible. See ME.scanNegativeDensity for why
            % CheckMEPositiveDensity is not used here.
            isNegDensity = false; fmin = 0; tmin = 0;
            if checkDensity
                [isNegDensity, fmin, tmin] = ME.scanNegativeDensity(alpha, A);
            end
            if isNegDensity
                line_warning(mfilename, 'The ME representation has a negative density: f(t) = -alpha*expm(A*t)*A*e reaches %g at t = %g. Moments and transforms remain well defined, but evalPDF returns negative values and sample() will not reproduce a proper distribution.\n', fmin, tmin);
            end

            % Store parameters
            self.alpha = alpha;
            self.A = A;
            self.nPhases = length(alpha);

            % Set parameters
            setParam(self, 1, 'alpha', alpha);
            setParam(self, 2, 'A', A);

            % Create Java object
            alphaMatrix = jline.util.matrix.Matrix(alpha);
            AMatrix = jline.util.matrix.Matrix(A);
            self.obj = jline.lang.processes.ME(alphaMatrix, AMatrix);

            % Build process representation: {D0=A, D1=-A*e*alpha'}
            % where e is column vector of ones
            e = ones(self.nPhases, 1);
            self.process = {A, -A * e * alpha};

            self.immediate = false;
        end

        function X = sample(self, n)
            % X = SAMPLE(N)
            % Get n samples from the distribution using inverse CDF interpolation

            if nargin < 2
                n = 1;
            end

            % Use me_sample for accurate sampling
            X = me_sample(self.process, n);
        end

        function phases = getNumberOfPhases(self)
            % PHASES = GETNUMBEROFPHASES()
            % Get number of phases in the ME representation
            phases = self.nPhases;
        end

        function Ft = evalCDF(self, t)
            % FT = EVALCDF(SELF, T)
            % Evaluate the cumulative distribution function at t
            %
            % For ME distribution: CDF(t) = 1 - alpha * exp(A*t) * e

            Ft = map_cdf(self.process, t);
        end

        function ft = evalPDF(self, t)
            % FT = EVALPDF(SELF, T)
            % Evaluate the probability density function at t
            %
            % For ME distribution: PDF(t) = -alpha * exp(A*t) * A * e

            ft = map_pdf(self.process, t);
        end

        function L = evalLST(self, s)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform at s

            % LST(s) = alpha * (sI - A)^(-1) * (-A) * e
            e = ones(self.nPhases, 1);
            sI = s * eye(self.nPhases);
            L = self.alpha * ((sI - self.A) \ (-self.A * e));
        end

        function mean_val = getMean(self)
            % MEAN_VAL = GETMEAN()
            % Get mean of the ME distribution, m1 = -alpha*inv(A)*e
            %
            % The definition is used rather than map_mean, which obtains the
            % rate from the stationary vector of D0+D1. That vector is a
            % probabilistic object of a Markovian process, and solving for it on
            % an ME degrades with the oscillation of A: a CME of order 101 came
            % out with a relative error of 5.6e-8 in the mean and 2.6e-4 in the
            % SCV, where the definition is exact to 1e-13. Native Python
            % ME.getMean and jline.lang.processes.ME do the same.

            e = ones(self.nPhases, 1);
            mean_val = -self.alpha * (self.A \ e);
        end

        function var_val = getVar(self)
            % VAR_VAL = GETVAR()
            % Get variance of the ME distribution, m2 - m1^2 with m2 = 2*alpha*inv(A)^2*e

            e = ones(self.nPhases, 1);
            Ainve = self.A \ e;
            m1 = -self.alpha * Ainve;
            m2 = 2 * self.alpha * (self.A \ Ainve);
            var_val = m2 - m1^2;
        end

        function scv = getSCV(self)
            % SCV = GETSCV()
            % Get squared coefficient of variation, var/mean^2

            scv = self.getVar() / self.getMean()^2;
        end

        function proc = getProcess(self)
            % PROC = GETPROCESS()
            % Get process representation {D0, D1}

            proc = self.process;
        end
    end

    methods(Static)
        function me = fitMoments(moms)
            % ME = FITMOMENTS(MOMS)
            % Create ME distribution by fitting the given moments
            % Uses BuTools MEFromMoments algorithm
            %
            % @param moms Array of moments (requires 2*M-1 moments for order M)
            % @return me ME distribution matching the given moments

            [alpha, A] = MEFromMoments(moms);
            me = ME(alpha, A);
        end

        function me = fromExp(rate)
            % ME = FROMEXP(RATE)
            % Create ME distribution from exponential distribution
            % Convenience method showing that Exp is a special case of ME
            %
            % @param rate Rate parameter (lambda)
            % @return me ME distribution equivalent to Exp(rate)

            alpha = 1.0;
            A = -rate;
            me = ME(alpha, A);
        end

        function me = fromErlang(k, rate)
            % ME = FROMERLANG(K, RATE)
            % Create ME distribution from Erlang distribution
            % Convenience method showing that Erlang is a special case of ME
            %
            % @param k Number of phases
            % @param rate Rate parameter for each phase
            % @return me ME distribution equivalent to Erlang(k, rate)

            alpha = zeros(1, k);
            alpha(1) = 1.0;  % alpha = [1, 0, 0, ..., 0]

            A = zeros(k, k);
            for i = 1:k
                A(i, i) = -rate;  % diagonal
                if i < k
                    A(i, i+1) = rate;  % super-diagonal
                end
            end

            me = ME(alpha, A);
        end

        function me = fromHyperExp(p, rates)
            % ME = FROMHYPEREXP(P, RATES)
            % Create ME distribution from HyperExponential distribution
            % Convenience method showing that HyperExp is a special case of ME
            %
            % @param p Array of probabilities for each branch
            % @param rates Array of rates for each branch
            % @return me ME distribution equivalent to HyperExp(p, rates)

            if length(p) ~= length(rates)
                error('p and rates must have the same length');
            end

            k = length(p);
            alpha = p;  % alpha = p

            A = zeros(k, k);
            for i = 1:k
                A(i, i) = -rates(i);  % diagonal matrix of rates
            end

            me = ME(alpha, A);
        end

        function [isNeg, fmin, tmin] = scanNegativeDensity(alpha, A)
            % [ISNEG, FMIN, TMIN] = SCANNEGATIVEDENSITY(ALPHA, A)
            % Search the density f(t) = -alpha*expm(A*t)*A*e for a negative
            % value.
            %
            % A negative value found here is a witness: it proves that the
            % representation is not a distribution. Finding none proves
            % nothing, so the caller must not report the converse.
            %
            % This replaces CheckMEPositiveDensity as the trigger for the
            % construction-time warning. That routine searches for a Markovian
            % monocyclic equivalent, which is a sufficient condition only, and
            % its verdict depends on the representation rather than on the
            % distribution: for alpha=[1,0,0], A=[-0.5 0 0; 0 -1 w; 0 -w -1]
            % the distribution is Exp(0.5) for every w, yet the search fails
            % once w >= 2*pi. It also costs of the order of a second per call
            % at search order 1000, which is far too slow for a constructor.
            %
            % The horizon covers all but scanTail of the mass, using the
            % dominant (least negative) eigenvalue of A; the sampling rate
            % resolves the fastest oscillation present, taken from the largest
            % imaginary part. The constants match native Python
            % line_solver.distributions.markovian and jline.lang.processes.ME.

            scanTail = 1e-12;        % residual mass left beyond the horizon
            scanHorizonCap = 1e4;    % cap on the horizon for near-degenerate A
            scanPerPeriod = 20;      % samples per period of fastest oscillation
            scanMinPts = 2001;
            scanMaxPts = 200001;
            scanRelTol = 1e-10;      % negative only if below -reltol*max|f|

            isNeg = false; fmin = 0; tmin = 0;

            alpha = reshape(alpha, 1, numel(alpha));
            n = size(A, 1);

            lambda = eig(A);
            decay = max(real(lambda));
            if ~isfinite(decay) || decay >= 0
                % Not a valid ME; CheckMERepresentation has already rejected it.
                return;
            end
            horizon = min(-log(scanTail) / abs(decay), scanHorizonCap);

            npts = scanMinPts;
            wmax = max(abs(imag(lambda)));
            if wmax > 0
                npts = max(npts, ceil(horizon * wmax * scanPerPeriod / (2*pi)) + 1);
            end
            npts = min(npts, scanMaxPts);

            step = horizon / (npts - 1);
            % One matrix exponential, then propagate: v_k = alpha*expm(A*k*step)
            E = expm(A * step);
            ve = -(A * ones(n, 1));

            v = alpha;
            fmin = Inf;
            tmin = 0;
            fabsmax = 0;
            for k = 0:npts-1
                f = v * ve;
                if f < fmin
                    fmin = f;
                    tmin = k * step;
                end
                if abs(f) > fabsmax
                    fabsmax = abs(f);
                end
                v = v * E;
            end

            isNeg = fmin < -scanRelTol * fabsmax;
        end
    end
end
