function [zlb, x, info] = npfqn_bnd_bpt(lambda0, mu, P, stationOf, c)
% [ZLB, X, INFO] = NPFQN_BND_BPT(LAMBDA0, MU, P, STATIONOF, C)
%
% First-order linear-programming relaxation of the achievable region of a
% multiclass open Markovian queueing network. Returns a LOWER bound ZLB on
% sum_r C(r)*x_r, where x_r is the mean sojourn time of class r, valid for
% EVERY non-idling scheduling policy.
%
% A "class" here is a buffer with its own exponential service rate and its own
% Markovian routing, so a station that serves several customer types owns one
% class per type, and a customer type that visits a station twice owns two
% classes. The network is open: class r receives external Poisson arrivals at
% rate LAMBDA0(r) and, on completing service, becomes class r' with probability
% P(r,r') or leaves with the deficit probability 1 - sum_r' P(r,r').
%
% METHOD. Uniformize the chain and let R(t) = sum_r f(r) n_r(t) for an
% arbitrary vector f. Writing the steady-state balance of E[R^2] gives an
% identity that is quadratic in f; since it holds for every f, the coefficient
% matrices of the two sides agree entrywise. The diagonal entries give one
% equation per class and the off-diagonal entries one per unordered pair, in
% the variables
%
%   x_r    = E[T_r],  the mean sojourn time of class r,
%   I(r,l) = E[1{server sigma(r) busy with class r} * n_l],
%   N(i,l) = E[1{server i idle} * n_l].
%
% A third block states that the events "station i serves class r", r in C_i,
% and "station i idle" are mutually exclusive and exhaustive, so their I and N
% terms sum to E[n_l] = lambda_l x_l. Minimizing C'x over this polyhedron is a
% relaxation of the true achievable region, hence a lower bound. This is the
% nonparametric derivation of the reference (obtained independently by Kumar
% and Kumar 1994); it dominates the parametric potential-function bound and
% needs O(K^2) variables and constraints rather than O(2^K) constraints.
%
% EXACT ON M/M/1. With one class the LP reads mu*I11 - lambda^2*x = lambda and
% I11 + N11 = lambda*x with N11 >= 0, whence x >= 1/(mu-lambda) with equality.
%
% NOT INCLUDED, DELIBERATELY. The inequality I(r,r) >= rho_r (a server busy
% with class r holds at least one class-r customer) is valid and would tighten
% the relaxation, but it is not part of the reference's characterization, and
% reproducing the reference's published bounds is the acceptance test here.
%
% Inputs:
%   LAMBDA0   K x 1, external Poisson arrival rate into each class (0 if none)
%   MU        K x 1, exponential service rate of each class
%   P         K x K, P(r,r') = probability class r becomes class r' after
%             service; row sums must not exceed 1
%   STATIONOF K x 1, index in 1..M of the station serving each class
%   C         K x 1, objective weights (default: all ones)
%
% Outputs:
%   ZLB       lower bound on sum_r C(r)*x_r
%   X         K x 1, the x block of the LP optimizer. Only the objective value
%             is a bound; an individual X(r) is a vertex coordinate, not a
%             bound on class r, unless C is the r-th unit vector
%   INFO      struct with fields lambda (effective rates), rho (per-class
%             utilizations), rhoStation, exitflag, nvars, nrows
%
% Reference: D. Bertsimas, I. Paschalidis, J. Tsitsiklis (1994). Optimization
% of multiclass queueing networks: polyhedral and nonlinear characterizations
% of achievable performance. Annals of Applied Probability 4(1), 43-75. See
% also D. Bertsimas (1995), Queueing Systems 21, 337-389, Theorem 9, which
% restates the same characterization.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

lambda0 = lambda0(:);
mu = mu(:);
stationOf = stationOf(:);
K = numel(lambda0);
if nargin < 5 || isempty(c)
    c = ones(K,1);
end
c = c(:);
if numel(mu) ~= K || numel(stationOf) ~= K || numel(c) ~= K
    line_error(mfilename, 'lambda0, mu, stationOf and c must all have %d entries.', K);
end
if size(P,1) ~= K || size(P,2) ~= K
    line_error(mfilename, 'P must be %dx%d.', K, K);
end
if any(mu <= 0)
    line_error(mfilename, 'Every class needs a strictly positive service rate.');
end
if any(sum(P,2) > 1 + 1e-9)
    line_error(mfilename, 'The routing matrix has a row summing above one.');
end
M = max(stationOf);

% ----- traffic equations, lambda = lambda0 + P'*lambda -----
lambda = (eye(K) - P') \ lambda0;
if any(lambda < -1e-9)
    line_error(mfilename, 'The traffic equations have no nonnegative solution.');
end
lambda = max(lambda, 0);
rho = lambda ./ mu;
rhoStation = zeros(M,1);
for i = 1:M
    rhoStation(i) = sum(rho(stationOf == i));
end
if any(rhoStation >= 1 - 1e-12)
    line_error(mfilename, ...
        'Station %d is saturated (rho=%.6g): no policy stabilizes the network.', ...
        find(rhoStation >= 1 - 1e-12, 1), max(rhoStation));
end

% ----- variable layout -----
%   x(r)    -> r
%   I(r,l)  -> K + (r-1)*K + l
%   N(i,l)  -> K + K*K + (i-1)*K + l
ix = @(r) r;
iI = @(r,l) K + (r-1)*K + l;
iN = @(i,l) K + K*K + (i-1)*K + l;
nv = K + K*K + M*K;

ri = []; cj = []; vv = []; beq = [];
nrow = 0;

% ----- (a) diagonal equations, one per class (test function n_r^2) -----
% 2*mu_r*I(r,r) - 2*sum_w mu_w*P(w,r)*I(w,r) - 2*lambda0_r*lambda_r*x_r
%   = 2*lambda_r*(1 - P(r,r))
for r = 1:K
    ci = []; cv = [];
    ci(end+1) = iI(r,r); cv(end+1) = 2*mu(r);
    for w = 1:K
        if P(w,r) ~= 0
            ci(end+1) = iI(w,r); cv(end+1) = -2*mu(w)*P(w,r); %#ok<AGROW>
        end
    end
    ci(end+1) = ix(r); cv(end+1) = -2*lambda0(r)*lambda(r); %#ok<AGROW>
    nrow = nrow + 1;
    ri = [ri, nrow*ones(1,numel(ci))]; cj = [cj, ci]; vv = [vv, cv];
    beq(nrow,1) = 2*lambda(r)*(1 - P(r,r));
end

% ----- (b) off-diagonal equations, one per unordered pair (n_r*n_s) -----
% mu_r*I(r,s) + mu_s*I(s,r) - sum_w mu_w*P(w,r)*I(w,s) - sum_w mu_w*P(w,s)*I(w,r)
%   - lambda0_r*lambda_s*x_s - lambda0_s*lambda_r*x_r
%   = -lambda_r*P(r,s) - lambda_s*P(s,r)
for r = 2:K
    for s = 1:(r-1)
        ci = []; cv = [];
        ci(end+1) = iI(r,s); cv(end+1) = mu(r); %#ok<AGROW>
        ci(end+1) = iI(s,r); cv(end+1) = mu(s); %#ok<AGROW>
        for w = 1:K
            if P(w,r) ~= 0
                ci(end+1) = iI(w,s); cv(end+1) = -mu(w)*P(w,r); %#ok<AGROW>
            end
            if P(w,s) ~= 0
                ci(end+1) = iI(w,r); cv(end+1) = -mu(w)*P(w,s); %#ok<AGROW>
            end
        end
        ci(end+1) = ix(s); cv(end+1) = -lambda0(r)*lambda(s); %#ok<AGROW>
        ci(end+1) = ix(r); cv(end+1) = -lambda0(s)*lambda(r); %#ok<AGROW>
        nrow = nrow + 1;
        ri = [ri, nrow*ones(1,numel(ci))]; cj = [cj, ci]; vv = [vv, cv];
        beq(nrow,1) = -lambda(r)*P(r,s) - lambda(s)*P(s,r);
    end
end

% ----- (c) exhaustiveness at each station: sum_{r in C_i} I(r,l) + N(i,l) = E[n_l] -----
for i = 1:M
    Ci = find(stationOf == i);
    for l = 1:K
        ci = []; cv = [];
        for r = Ci(:)'
            ci(end+1) = iI(r,l); cv(end+1) = 1; %#ok<AGROW>
        end
        ci(end+1) = iN(i,l); cv(end+1) = 1; %#ok<AGROW>
        ci(end+1) = ix(l); cv(end+1) = -lambda(l); %#ok<AGROW>
        nrow = nrow + 1;
        ri = [ri, nrow*ones(1,numel(ci))]; cj = [cj, ci]; vv = [vv, cv];
        beq(nrow,1) = 0;
    end
end

% ----- assemble and solve -----
Aeq = sparse(ri, cj, vv, nrow, nv);

f = zeros(nv,1);
f(1:K) = c;
lb = zeros(nv,1);

lpopt = optimoptions('linprog', 'Display', 'off');
[sol, fval, exitflag] = linprog(f, [], [], Aeq, beq, lb, [], lpopt);

if exitflag ~= 1 || isempty(sol)
    line_error(mfilename, ...
        'The achievable-region LP did not solve to optimality (exitflag %d).', exitflag);
end

zlb = fval;
x = sol(1:K);
info = struct('lambda', lambda, 'rho', rho, 'rhoStation', rhoStation, ...
    'exitflag', exitflag, 'nvars', nv, 'nrows', numel(beq));
end
