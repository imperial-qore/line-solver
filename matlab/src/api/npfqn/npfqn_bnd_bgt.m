function [Qub, info] = npfqn_bnd_bgt(lambda, mu, sigma, J)
% [QUB, INFO] = NPFQN_BND_BGT(LAMBDA, MU, SIGMA, J)
%
% Piecewise-linear Lyapunov UPPER bound on the steady-state queue lengths of a
% multitype (deterministic-routing) multiclass Markovian queueing network,
% valid for EVERY work-conserving Markovian policy.
%
% MODEL. J single-server stations; I customer types; type i arrives as a
% Poisson stream of rate LAMBDA(i) and passes through stages
% k = 1..numel(MU{i}), stage k being served at station SIGMA{i}(k) at
% exponential rate MU{i}(k). Class (i,k) is the buffer of type i at stage k;
% N = sum_i numel(MU{i}) is the number of classes.
%
% METHOD. Solve the Down-Meyn global-stability linear program GLP[dm], eq.
% (25)-(28) of the reference, in the piecewise-linear Lyapunov function
% Phi(x) = max_j L^j'x:
%
%   L^j(i,1) lambda_i + mu(i,k) (L^j(i,k+1) - L^j(i,k)) + V_j <= -gamma
%                                                    for (i,k) in station j
%   mu(i,k) (L^j(i,k+1) - L^j(i,k)) <= V_j           for (i,k) not in j
%   (1/(J-1)) sum_{j' ~= j} L^j'(i,k) >= L^j(i,k)    for (i,k) not in j
%   L, V, gamma >= 0
%
% with the convention L^j(i,Ji+1) = 0. A feasible solution with gamma > 0
% certifies that EVERY work-conserving policy is stable, and a smoothed
% Phi is then a Lyapunov function with drift gamma/4 and an explicit exception
% parameter, which turns into the bound of the reference's Theorem 4:
%
%   E[L^j'Q] <= 16 N J^2 (J-1) (Lmax+gamma)^3 / gamma^2
%               + 8 (Lmax + gamma/2)^2 / gamma        =: U,
%
% Lmax = max L^j(i,k), for every j. Since L >= 0 and Q >= 0, this gives
% E[Q(i,k)] <= U / max_j L^j(i,k) per class, and the geometric tail
%
%   P( L^j'Q >= B + 2(Lmax+gamma/2) m ) <= ((Lmax+gamma/2)/(Lmax+3gamma/4))^m.
%
% THE RATES ARE RESCALED so that sum_i lambda_i + sum_{i,k} mu(i,k) = 1, the
% uniformization the reference imposes before Theorem 4. Queue lengths are
% counts and are unaffected by the time scale, so QUB is in the original units.
%
% NORMALIZATION, WHICH THE REFERENCE LEAVES OPEN. GLP[dm] is homogeneous, and
% so is the bound: scaling (L,V,gamma) by t > 0 scales U by t and every
% denominator L^j(i,k) by t. This routine therefore fixes L^j(i,k) <= 1 and
% MAXIMIZES gamma, then breaks ties among gamma-optimal solutions by
% maximizing sum L -- a degenerate optimum can otherwise leave some
% L^j(i,k) = 0 and report an infinite bound for a class for no reason.
%
% THE BOUND IS LOOSE, and knowingly so: the exception parameter carries
% (Lmax+gamma)^3/gamma^2 and dominates as soon as J > 1. On M/M/1 (J=1, where
% that term vanishes) it is about 23x the exact mean queue length; on a
% two-station tandem it is three orders of magnitude above it. What is sharp is
% the STABILITY CERTIFICATE gamma > 0 and the geometric tail RATE; the constant
% in front of the tail is not.
%
% Inputs:
%   LAMBDA  I x 1, Poisson arrival rate of each type
%   MU      1 x I cell, MU{i} is the vector of stage service rates of type i
%   SIGMA   1 x I cell, SIGMA{i}(k) in 1..J is the station of stage k
%   J       number of stations (default: max over SIGMA)
%
% Outputs:
%   QUB     1 x I cell, QUB{i}(k) upper bound on E[Q(i,k)]; Inf where the LP
%           optimum leaves max_j L^j(i,k) = 0
%   INFO    struct with fields gamma, Lmax, L (J x N), V (J x 1), B (exception
%           parameter), U (the Theorem 4 bound on E[L^j'Q]), tailRatio (the
%           geometric decay ratio), tailStep (2(Lmax+gamma/2)), rho (per-class
%           load), rhoStation (J x 1), scale (the uniformization divisor),
%           classType, classStage, classStation (N x 1 index maps)
%
% Reference: D. Bertsimas, D. Gamarnik, J. N. Tsitsiklis (2001). Performance of
% multiclass Markovian queueing networks via piecewise linear Lyapunov
% functions. Annals of Applied Probability 11(4), 1384-1428, Section 5.1
% (GLP[dm] of Down and Meyn 1997, and Theorem 4).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

lambda = lambda(:);
I = numel(lambda);
if numel(mu) ~= I || numel(sigma) ~= I
    line_error(mfilename, 'MU and SIGMA must both have %d entries.', I);
end
if nargin < 4 || isempty(J)
    J = 0;
    for i = 1:I
        J = max(J, max(sigma{i}));
    end
end

% ----- flatten (i,k) into a class index -----
classType = []; classStage = []; classStation = []; muc = [];
firstOf = zeros(I,1);
for i = 1:I
    Ji = numel(mu{i});
    if numel(sigma{i}) ~= Ji
        line_error(mfilename, 'MU{%d} and SIGMA{%d} have different lengths.', i, i);
    end
    if Ji == 0
        line_error(mfilename, 'Type %d has no stage.', i);
    end
    firstOf(i) = numel(classType) + 1;
    for k = 1:Ji
        classType(end+1,1) = i; %#ok<AGROW>
        classStage(end+1,1) = k; %#ok<AGROW>
        classStation(end+1,1) = sigma{i}(k); %#ok<AGROW>
        muc(end+1,1) = mu{i}(k); %#ok<AGROW>
    end
end
N = numel(muc);
if any(muc <= 0)
    line_error(mfilename, 'Every stage needs a strictly positive service rate.');
end
if any(lambda <= 0)
    line_error(mfilename, 'Every type needs a strictly positive arrival rate.');
end
nextOf = zeros(N,1);
for c = 1:N
    if c < N && classType(c+1) == classType(c)
        nextOf(c) = c + 1;
    end
end

% ----- loads -----
rho = zeros(N,1);
for c = 1:N
    rho(c) = lambda(classType(c)) / muc(c);
end
rhoStation = zeros(J,1);
for j = 1:J
    rhoStation(j) = sum(rho(classStation == j));
end
if any(rhoStation >= 1)
    line_error(mfilename, ...
        'Station %d is saturated (rho=%.6g): the load condition of the reference fails.', ...
        find(rhoStation >= 1, 1), max(rhoStation));
end

% ----- uniformization: rescale so sum(lambda) + sum(mu) = 1 -----
scale = sum(lambda) + sum(muc);
lam = lambda / scale;
mus = muc / scale;

% ----- LP layout: L(j,c) -> (j-1)*N + c ; V(j) -> J*N + j ; gamma -> J*N+J+1
iL = @(j,c) (j-1)*N + c;
iV = @(j) J*N + j;
ig = J*N + J + 1;
nv = J*N + J + 1;

A = []; b = [];
    rowsI = []; rowsJ = []; rowsV = []; nrow = 0;
for j = 1:J
    for c = 1:N
        ci = []; cv = [];
        if classStation(c) == j
            % (25): the station's own service term, the arrival term of the
            % type this class belongs to, V_j and gamma
            ci(end+1) = iL(j, firstOf(classType(c))); cv(end+1) = lam(classType(c)); %#ok<AGROW>
            ci(end+1) = iL(j,c); cv(end+1) = -mus(c); %#ok<AGROW>
            if nextOf(c) > 0
                ci(end+1) = iL(j,nextOf(c)); cv(end+1) = mus(c); %#ok<AGROW>
            end
            ci(end+1) = iV(j); cv(end+1) = 1; %#ok<AGROW>
            ci(end+1) = ig; cv(end+1) = 1; %#ok<AGROW>
            nrow = nrow + 1;
            rowsI = [rowsI, nrow*ones(1,numel(ci))]; rowsJ = [rowsJ, ci]; rowsV = [rowsV, cv];
            b(nrow,1) = 0; %#ok<AGROW>
        else
            % (26): a class served elsewhere contributes at most V_j
            ci(end+1) = iL(j,c); cv(end+1) = -mus(c); %#ok<AGROW>
            if nextOf(c) > 0
                ci(end+1) = iL(j,nextOf(c)); cv(end+1) = mus(c); %#ok<AGROW>
            end
            ci(end+1) = iV(j); cv(end+1) = -1; %#ok<AGROW>
            nrow = nrow + 1;
            rowsI = [rowsI, nrow*ones(1,numel(ci))]; rowsJ = [rowsJ, ci]; rowsV = [rowsV, cv];
            b(nrow,1) = 0; %#ok<AGROW>
            % (27): the averaging condition across the other stations
            if J > 1
                ci = iL(j,c); cv = 1;
                for jp = 1:J
                    if jp ~= j
                        ci(end+1) = iL(jp,c); cv(end+1) = -1/(J-1); %#ok<AGROW>
                    end
                end
                nrow = nrow + 1;
                rowsI = [rowsI, nrow*ones(1,numel(ci))]; rowsJ = [rowsJ, ci]; rowsV = [rowsV, cv];
                b(nrow,1) = 0; %#ok<AGROW>
            end
        end
    end
end
A = sparse(rowsI, rowsJ, rowsV, nrow, nv);

lb = zeros(nv,1);
ub = inf(nv,1);
ub(1:J*N) = 1;   % the homogeneous normalization Lmax <= 1

lpopt = optimoptions('linprog', 'Display', 'off');
f = zeros(nv,1); f(ig) = -1;   % maximize gamma
[sol, fval, exitflag] = linprog(f, A, b, [], [], lb, ub, lpopt);
if exitflag ~= 1 || isempty(sol)
    line_error(mfilename, ...
        'GLP[dm] did not solve to optimality (exitflag %d).', exitflag);
end
gamma = -fval;
if gamma <= 0
    line_error(mfilename, ...
        ['GLP[dm] has no solution with gamma > 0: this network is not certified ' ...
         'globally stable, so no finite piecewise-linear Lyapunov bound exists.']);
end

% Tie-break among gamma-optimal solutions: maximize sum L, so a degenerate
% vertex does not report an infinite bound for a class it zeroed arbitrarily.
A2 = [A; sparse(1, ig, -1, 1, nv)];
b2 = [b; -gamma];
f2 = zeros(nv,1); f2(1:J*N) = -1;
[sol2, ~, ef2] = linprog(f2, A2, b2, [], [], lb, ub, lpopt);
if ef2 == 1 && ~isempty(sol2)
    sol = sol2;
    gamma = sol(ig);
end

L = reshape(sol(1:J*N), N, J)';
V = sol(J*N+(1:J));
Lmax = max(L(:));

% ----- Theorem 4 -----
B = 16*N*J^2*(J-1)*(Lmax + gamma)^3 / gamma^2;
U = B + 8*(Lmax + gamma/2)^2 / gamma;
tailStep = 2*(Lmax + gamma/2);
tailRatio = (Lmax + gamma/2) / (Lmax + 0.75*gamma);

% E[Q(i,k)] <= U / max_j L^j(i,k), since L >= 0 and Q >= 0
Lbest = max(L, [], 1)';
qub = inf(N,1);
pos = Lbest > 0;
qub(pos) = U ./ Lbest(pos);

Qub = cell(1,I);
for i = 1:I
    Qub{i} = qub(classType == i)';
end

info = struct('gamma', gamma, 'Lmax', Lmax, 'L', L, 'V', V, 'B', B, 'U', U, ...
    'tailRatio', tailRatio, 'tailStep', tailStep, 'rho', rho, ...
    'rhoStation', rhoStation, 'scale', scale, 'classType', classType, ...
    'classStage', classStage, 'classStation', classStation);
end
