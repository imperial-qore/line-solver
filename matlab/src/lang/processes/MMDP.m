classdef MMDP < Markovian
    % Markov-Modulated Deterministic Process for fluid queue modeling
    %
    % Models fluid flow with deterministic rates modulated by a background
    % Markov chain. Suitable for arrival and service processes in
    % Markovian fluid queues analyzed by the mfq method of SolverFLD.
    %
    % The (Q, R) parameterization follows BUTools conventions:
    % - Q: Generator matrix of the modulating CTMC
    % - R: Diagonal matrix of deterministic rates per state
    %
    % MMDP is the deterministic analogue of MMPP:
    % - MMPP: Poisson arrival rates modulated by a Markov chain
    % - MMDP: Deterministic rates modulated by a Markov chain
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = MMDP(Q, R)
            % MMDP Create a Markov-Modulated Deterministic Process
            %
            % @brief Creates an MMDP with specified generator Q and rate matrix R
            % @param Q n×n generator matrix (row sums = 0)
            % @param R n×n diagonal matrix of rates, OR n-vector
            % @return self MMDP instance with specified matrices

            self@Markovian('MMDP', 2);

            % Convert vector to diagonal matrix if needed
            if isvector(R)
                R = diag(R(:));
            end

            setParam(self, 1, 'Q', Q);
            setParam(self, 2, 'R', R);

            % Store in process format: {Q, R}
            self.process = {Q, R};

            % Validate
            if ~mmdp_isfeasible(Q, R)
                line_warning(mfilename, 'MMDP is infeasible.\n');
            end
        end

        function Q_mat = Q(self)
            % Q Return the generator matrix
            %
            % @return Q_mat n×n generator matrix of the modulating CTMC

            Q_mat = self.getParam(1).paramValue;
        end

        function R_mat = R(self)
            % R Return the rate matrix (diagonal)
            %
            % @return R_mat n×n diagonal matrix of deterministic rates

            R_mat = self.getParam(2).paramValue;
        end

        function r_vec = r(self)
            % r Return the rate vector (diagonal of R)
            %
            % @return r_vec n-vector of deterministic rates per phase

            r_vec = diag(self.R());
        end

        function n = getNumberOfPhases(self)
            % GETNUMBEROFPHASES Return the number of phases
            %
            % @return n Number of phases in the modulating CTMC

            n = size(self.Q(), 1);
        end

        function rate = getMeanRate(self)
            % GETMEANRATE Compute stationary mean rate
            %
            % Computes E[r] = π * r, where π is the stationary distribution
            % of the modulating CTMC with generator Q.
            %
            % @return rate Stationary mean deterministic rate

            Q_mat = self.Q();
            r_vec = self.r();
            pi = ctmc_solve(Q_mat);
            rate = pi * r_vec;
        end

        function rate = getRate(self)
            % GETRATE Alias for getMeanRate
            %
            % @return rate Stationary mean deterministic rate

            rate = self.getMeanRate();
        end

        function mean = getMean(self)
            % GETMEAN Return mean inter-arrival time (inverse of mean rate)
            %
            % @return mean Mean inter-arrival time (1/rate), or Inf if rate=0

            rate = self.getMeanRate();
            if rate > 0
                mean = 1 / rate;
            else
                mean = Inf;
            end
        end

        function scv = getSCV(self)
            % GETSCV Return squared coefficient of variation of rates
            %
            % Computes Var[r]/E[r]^2 where expectation is over the
            % stationary distribution of the modulating CTMC.
            %
            % @return scv Squared coefficient of variation

            Q_mat = self.Q();
            r_vec = self.r();
            pi = ctmc_solve(Q_mat);
            mean_rate = pi * r_vec;
            var_rate = pi * (r_vec.^2) - mean_rate^2;
            if mean_rate > 0
                scv = var_rate / (mean_rate^2);
            else
                scv = Inf;
            end
        end

        function bool = isImmediate(self)
            % ISIMMEDIATE Check if the process has infinite rate
            %
            % @return bool True if mean rate is effectively infinite

            bool = self.getMean < GlobalConstants.FineTol;
        end
    end

    methods (Static)
        function mmdp = fromMAP(map_obj)
            % FROMMAP Convert a MAP to MMDP (deterministic representation)
            %
            % Converts a Markovian Arrival Process to a Markov-Modulated
            % Deterministic Process by extracting the full generator and
            % using row sums of D1 as the deterministic rates.
            %
            % @param map_obj MAP object to convert
            % @return mmdp MMDP representation of the MAP

            D0 = map_obj.D(0);
            D1 = map_obj.D(1);
            Q = D0 + D1;
            R = diag(sum(D1, 2));  % Row sums of D1 as diagonal
            mmdp = MMDP(Q, R);
        end

        function mmdp = fromMMPP2(lambda0, lambda1, sigma0, sigma1)
            % FROMMMPP2 Create MMDP from MMPP2 parameters
            %
            % Creates a 2-state MMDP using the same parameterization as MMPP2.
            %
            % @param lambda0 Rate in state 0
            % @param lambda1 Rate in state 1
            % @param sigma0 Transition rate from state 0 to state 1
            % @param sigma1 Transition rate from state 1 to state 0
            % @return mmdp 2-state MMDP

            Q = [-sigma0, sigma0; sigma1, -sigma1];
            R = diag([lambda0; lambda1]);
            mmdp = MMDP(Q, R);
        end
    end
end
