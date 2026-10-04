classdef MMDP2 < MMDP
    % 2-state Markov-Modulated Deterministic Process
    %
    % A specialized MMDP with exactly 2 phases, using a convenient
    % parameterization analogous to MMPP2.
    %
    % Parameterization:
    %   r0, r1: Deterministic rates in states 0 and 1
    %   sigma0: Transition rate from state 0 to state 1
    %   sigma1: Transition rate from state 1 to state 0
    %
    % The generator matrix is:
    %   Q = [-sigma0, sigma0; sigma1, -sigma1]
    %
    % The rate matrix is:
    %   R = diag([r0, r1])
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = MMDP2(r0, r1, sigma0, sigma1)
            % MMDP2 Create a 2-state Markov-Modulated Deterministic Process
            %
            % @brief Creates a 2-state MMDP with specified rates and transitions
            % @param r0 Deterministic rate in state 0
            % @param r1 Deterministic rate in state 1
            % @param sigma0 Transition rate from state 0 to state 1
            % @param sigma1 Transition rate from state 1 to state 0
            % @return self MMDP2 instance

            Q = [-sigma0, sigma0; sigma1, -sigma1];
            R = diag([r0; r1]);
            self@MMDP(Q, R);
            self.name = 'MMDP2';

            % Override parameters with original scalar values for convenience
            self.params = cell(1, 4);
            setParam(self, 1, 'r0', r0);
            setParam(self, 2, 'r1', r1);
            setParam(self, 3, 'sigma0', sigma0);
            setParam(self, 4, 'sigma1', sigma1);
        end

        function rate = getMeanRate(self)
            % GETMEANRATE Compute stationary mean rate (closed-form)
            %
            % For a 2-state MMDP, the mean rate has the closed form:
            %   E[r] = (r0*σ1 + r1*σ0) / (σ0 + σ1)
            %
            % @return rate Stationary mean deterministic rate

            r0 = self.getParam(1).paramValue;
            r1 = self.getParam(2).paramValue;
            sigma0 = self.getParam(3).paramValue;
            sigma1 = self.getParam(4).paramValue;
            rate = (r0*sigma1 + r1*sigma0) / (sigma0 + sigma1);
        end

        function scv = getSCV(self)
            % GETSCV Compute squared coefficient of variation (closed-form)
            %
            % For a 2-state MMDP, the SCV has a closed form based on
            % the variance of rates over the stationary distribution.
            %
            % @return scv Squared coefficient of variation

            r0 = self.getParam(1).paramValue;
            r1 = self.getParam(2).paramValue;
            sigma0 = self.getParam(3).paramValue;
            sigma1 = self.getParam(4).paramValue;

            % Stationary probabilities
            pi0 = sigma1 / (sigma0 + sigma1);
            pi1 = sigma0 / (sigma0 + sigma1);

            % Mean and variance
            mean_rate = pi0 * r0 + pi1 * r1;
            var_rate = pi0 * r0^2 + pi1 * r1^2 - mean_rate^2;

            if mean_rate > 0
                scv = var_rate / (mean_rate^2);
            else
                scv = Inf;
            end
        end

        function n = getNumberOfPhases(~)
            % GETNUMBEROFPHASES Return the number of phases
            %
            % @return n Always 2 for MMDP2

            n = 2;
        end

        function Q_mat = Q(self)
            % Q Return the generator matrix
            %
            % @return Q_mat 2×2 generator matrix

            sigma0 = self.getParam(3).paramValue;
            sigma1 = self.getParam(4).paramValue;
            Q_mat = [-sigma0, sigma0; sigma1, -sigma1];
        end

        function R_mat = R(self)
            % R Return the rate matrix (diagonal)
            %
            % @return R_mat 2×2 diagonal rate matrix

            r0 = self.getParam(1).paramValue;
            r1 = self.getParam(2).paramValue;
            R_mat = diag([r0; r1]);
        end

        function r_vec = r(self)
            % r Return the rate vector
            %
            % @return r_vec 2-vector of rates [r0; r1]

            r0 = self.getParam(1).paramValue;
            r1 = self.getParam(2).paramValue;
            r_vec = [r0; r1];
        end
    end
end
