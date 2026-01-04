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
        function self = ME(alpha, A)
            % ME Create a Matrix Exponential distribution instance
            %
            % @brief Creates an ME distribution with the given initial vector and matrix parameter
            % @param alpha Initial vector (may have negative entries or sum != 1)
            % @param A Matrix parameter (must have all eigenvalues with negative real parts)
            % @return self ME distribution instance

            % Call superclass constructor
            self@Markovian('ME', 2);

            % Validate using BuTools
            if ~CheckMERepresentation(alpha, A)
                error('Invalid ME representation: Check that A is square, alpha and A have compatible dimensions, all eigenvalues of A have negative real parts, and the dominant eigenvalue is real.');
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
            % Get mean of the ME distribution

            mean_val = map_mean(self.process);
        end

        function var_val = getVar(self)
            % VAR_VAL = GETVAR()
            % Get variance of the ME distribution

            var_val = map_var(self.process);
        end

        function scv = getSCV(self)
            % SCV = GETSCV()
            % Get squared coefficient of variation

            scv = map_scv(self.process);
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
    end
end
