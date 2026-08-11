function prOt = overtake_prob(self, eidx)
% OVERTAKE_PROB Compute overtaking probability using transient Markov chain
%
% PROT = OVERTAKE_PROB(SELF, EIDX) computes the probability that a new
% arrival to entry EIDX finds the server in phase-2 (post-reply processing).
%
% This uses a 3-state Continuous Time Markov Chain (CTMC):
%   State 0: Server idle
%   State 1: Server in phase-1 (caller is blocked)
%   State 2: Server in phase-2 (caller has been released)
%
% Transitions:
%   0 -> 1: arrival (rate lambda)
%   1 -> 2: phase-1 completion (rate mu1 = 1/S1)
%   2 -> 0: phase-2 completion (rate mu2 = 1/S2)
%
% By PASTA (Poisson Arrivals See Time Averages), the overtaking probability
% equals the steady-state probability of being in phase-2.
%
% For multi-server systems, an approximation is used.
%
% References:
%   Franks, G. and Woodside, M. "Effectiveness of early replies in
%   client-server systems", Performance Evaluation, 36(1):165-184, 1999.

    lqn = self.lqn;
    e = eidx - lqn.eshift;
    tidx = lqn.parent(eidx);

    S1 = self.servt_ph1(eidx);
    S2 = self.servt_ph2(eidx);

    % Get throughput - use entry if available, otherwise use task
    if self.tput(eidx) > GlobalConstants.FineTol
        lambda = self.tput(eidx);
    elseif self.tput(tidx) > GlobalConstants.FineTol
        lambda = self.tput(tidx);
    else
        lambda = 0;
    end

    c = lqn.mult(tidx);  % number of servers

    % Handle degenerate cases
    if S2 < GlobalConstants.FineTol || lambda < GlobalConstants.FineTol || S1 < GlobalConstants.FineTol
        prOt = 0;
        return;
    end

    mu1 = 1/S1;
    mu2 = 1/S2;

    if c == 1
        % Single server: exact CTMC solution
        % Generator matrix Q for states {idle, phase-1, phase-2}
        %
        % Q = | -lambda   lambda      0     |
        %     |    0      -mu1       mu1    |
        %     |   mu2       0       -mu2    |
        %
        % Steady-state: pi * Q = 0, sum(pi) = 1

        % Solve using the augmented system [Q'; 1 1 1] * pi = [0; 0; 0; 1]
        Q = [-lambda, lambda, 0;
             0, -mu1, mu1;
             mu2, 0, -mu2];

        A = [Q'; ones(1,3)];
        b = [0; 0; 0; 1];

        % Use least squares for numerical stability
        pi = A \ b;

        % PASTA: P(find in phase-2) = pi(3)
        prOt = max(0, min(1, pi(3)));
    else
        % Multi-server approximation
        % Based on utilization decomposition
        rho = lambda * (S1 + S2) / c;

        if rho >= 1
            % Saturated system: probability proportional to phase-2 fraction
            prOt = S2 / (S1 + S2);
        else
            % Probability a random server is in phase-2
            % Approximation: fraction of time in phase-2 weighted by utilization
            prOt = (S2 / (S1 + S2)) * rho;
        end

        % Bound the result
        prOt = max(0, min(1, prOt));
    end
end
