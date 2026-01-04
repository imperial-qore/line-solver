classdef RAP < Markovian
    % Rational Arrival Process (RAP) distribution
    %
    % RAP is a generalization of the Markovian Arrival Process (MAP) where
    % the matrices H0 and H1 represent hidden and visible transitions respectively,
    % but with relaxed constraints compared to MAP.
    %
    % Representation:
    % - H0: matrix of hidden transition rates (transitions without arrivals)
    % - H1: matrix of visible transition rates (transitions with arrivals)
    % - H0 + H1 must form a valid infinitesimal generator (row sums = 0)
    % - All eigenvalues of H0 must have negative real parts
    % - Dominant eigenvalue of H0 must be negative and real
    %
    % The marginal distribution of inter-arrival times is a Matrix Exponential (ME).
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        H0;  % Hidden transition matrix
        H1;  % Visible transition matrix
    end

    methods
        function self = RAP(H0, H1)
            % RAP Create a Rational Arrival Process instance
            %
            % @brief Creates a RAP with the given H0 and H1 matrices
            % @param H0 Hidden transition matrix (square matrix)
            % @param H1 Visible transition matrix (square matrix, same size as H0)
            % @return self RAP distribution instance

            % Call superclass constructor
            self@Markovian('RAP', 2);

            % Validate using BuTools
            if ~CheckRAPRepresentation(H0, H1)
                error('Invalid RAP representation: Check that H0 and H1 are square matrices of the same size, H0 + H1 forms a valid infinitesimal generator (row sums = 0), all eigenvalues of H0 have negative real parts, and the dominant eigenvalue of H0 is real.');
            end

            % Store parameters
            self.H0 = H0;
            self.H1 = H1;
            self.nPhases = size(H0, 1);

            % Set parameters
            setParam(self, 1, 'H0', H0);
            setParam(self, 2, 'H1', H1);

            % Create Java object
            H0Matrix = jline.util.matrix.Matrix(H0);
            H1Matrix = jline.util.matrix.Matrix(H1);
            self.obj = jline.lang.processes.RAP(H0Matrix, H1Matrix);

            % Build process representation: {D0=H0, D1=H1}
            % RAP uses same format as MAP
            self.process = {H0, H1};

            self.immediate = false;
        end

        function X = sample(self, n)
            % X = SAMPLE(N)
            % Get n samples from the distribution using rap_sample

            if nargin < 2
                n = 1;
            end

            % Use rap_sample for accurate sampling
            X = rap_sample(self.process, n);
        end

        function phases = getNumberOfPhases(self)
            % PHASES = GETNUMBEROFPHASES()
            % Get number of phases in the RAP representation
            phases = self.nPhases;
        end

        function Ft = evalCDF(self, t)
            % FT = EVALCDF(SELF, T)
            % Evaluate the cumulative distribution function at t
            %
            % For RAP, the marginal CDF is same as MAP

            Ft = map_cdf(self.process, t);
        end

        function L = evalLST(self, s)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform at s

            % For RAP marginal: LST(s) = pie * (sI - H0)^(-1) * (-H0) * e
            % where pie is the stationary distribution of H0 + H1
            pie = map_pie(self.process);
            n = self.nPhases;
            e = ones(n, 1);
            sI = s * eye(n);
            L = pie * ((sI - self.H0) \ (-self.H0 * e));
        end

        function mean_val = getMean(self)
            % MEAN_VAL = GETMEAN()
            % Get mean of the RAP distribution

            mean_val = map_mean(self.process);
        end

        function var_val = getVar(self)
            % VAR_VAL = GETVAR()
            % Get variance of the RAP distribution

            var_val = map_var(self.process);
        end

        function scv = getSCV(self)
            % SCV = GETSCV()
            % Get squared coefficient of variation

            scv = map_scv(self.process);
        end

        function lam = getRate(self)
            % LAM = GETRATE()
            % Get arrival rate (lambda) of the RAP

            lam = map_lambda(self.process);
        end

        function acf = getACF(self, lags)
            % ACF = GETACF(SELF, LAGS)
            % Get autocorrelation function at specified lags

            acf = map_acf(self.process, lags);
        end

        function idc = getIDC(self, t)
            % IDC = GETIDC(SELF, T)
            % Get index of dispersion for counts
            % If t is not provided, returns asymptotic IDC

            if nargin < 2
                idc = map_idc(self.process);
            else
                idc = map_count_var(self.process, t) / map_count_mean(self.process, t);
            end
        end

        function proc = getProcess(self)
            % PROC = GETPROCESS()
            % Get process representation {H0, H1}

            proc = self.process;
        end

        function H0_out = getH0(self)
            % H0_OUT = GETH0()
            % Get H0 matrix (hidden transitions)

            H0_out = self.H0;
        end

        function H1_out = getH1(self)
            % H1_OUT = GETH1()
            % Get H1 matrix (visible transitions)

            H1_out = self.H1;
        end
    end

    methods(Static)
        function rap = fromPoisson(rate)
            % RAP = FROMPOISSON(RATE)
            % Create RAP from exponential renewal process (Poisson)
            % Convenience method showing that Poisson is a special case of RAP
            %
            % @param rate Arrival rate (lambda)
            % @return rap RAP distribution equivalent to Poisson(rate)

            H0 = -rate;
            H1 = rate;
            rap = RAP(H0, H1);
        end

        function rap = fromErlang(k, rate)
            % RAP = FROMERLANG(K, RATE)
            % Create RAP from Erlang renewal process
            % Convenience method showing that Erlang is a special case of RAP
            %
            % @param k Number of phases
            % @param rate Rate parameter for each phase
            % @return rap RAP distribution equivalent to Erlang(k, rate)

            H0 = zeros(k, k);
            H1 = zeros(k, k);

            for i = 1:k
                H0(i, i) = -rate;  % diagonal
                if i < k
                    H0(i, i+1) = rate;  % transitions to next phase (no arrival)
                end
            end
            % Last phase transition with arrival back to first phase
            H1(k, 1) = rate;

            rap = RAP(H0, H1);
        end

        function rap = fromMAP(map)
            % RAP = FROMMAP(MAP)
            % Create RAP from Markovian Arrival Process
            % Convenience method showing that MAP is a special case of RAP
            %
            % @param map MAP distribution instance
            % @return rap RAP distribution equivalent to the given MAP

            if isa(map, 'MAP')
                proc = map.getProcess();
                D0 = proc{1};
                D1 = proc{2};
            elseif iscell(map)
                D0 = map{1};
                D1 = map{2};
            else
                error('Input must be a MAP object or a cell array {D0, D1}');
            end

            rap = RAP(D0, D1);
        end

        function rap = fitMoments(moms)
            % RAP = FITMOMENTS(MOMS)
            % Create RAP by fitting the given moments
            % Uses BuTools RAPFromMoments algorithm
            %
            % @param moms Array of moments
            % @return rap RAP distribution matching the given moments

            error('RAP.fitMoments() requires RAPFromMoments from BuTools. Use the RAP(H0, H1) constructor directly for now.');
        end

        function rap = fitMomentsAndCorrelations(moms, corrs)
            % RAP = FITMOMENTSANDCORRELATIONS(MOMS, CORRS)
            % Create RAP by fitting the given moments and correlations
            % Uses BuTools RAPFromMomentsAndCorrelations algorithm
            %
            % @param moms Array of moments
            % @param corrs Array of lag-k correlations
            % @return rap RAP distribution matching the given moments and correlations

            error('RAP.fitMomentsAndCorrelations() requires RAPFromMomentsAndCorrelations from BuTools. Use the RAP(H0, H1) constructor directly for now.');
        end
    end
end
