function [x, pi, Q, stats] = solver_mam_autocat(R, AP, varargin)
% SOLVER_MAM_AUTOCAT Native MATLAB implementation of autocat
%
% [X, PI, Q, STATS] = SOLVER_MAM_AUTOCAT(R, AP, OPTIONS)
%
% Searches for RCAT product-form solutions using various relaxation methods.
% Native MATLAB implementation without YALMIP dependency.
%
% Input:
%   R   - Cell array of rate matrices
%   AP  - Action-Process mapping (A x 2)
%   OPTIONS (name-value pairs):
%     'relaxation' - Relaxation method (default: 'auto')
%                    'lpr'      - LP relaxation with McCormick envelopes
%                    'tlpr0inf' - Tightened LP with 0-level + infinity
%                    'tlpr1inf' - Tightened LP with 1-level + infinity
%                    'tlprUinf' - Tightened LP with U-level + infinity
%                    'ens'      - Exact nonlinear system (fmincon)
%                    'qcp'      - Quadratically constrained (fmincon)
%                    'ma'       - Mean Approximation (first moment matching)
%                    'va'       - Variance Approximation (slack variables)
%                    'minres'   - Minimum Residual (L1-norm of RC3)
%                    'inap'     - Iterative fixed-point approximation
%                    'tma1inf'  - Mean Approx + Tightened LP with 1-level + infinity
%                    'zpr'      - Zero potential relaxation with cutting planes
%                    'tzpr0inf' - Tightened zero potential (0-level)
%                    'tzpr1inf' - Tightened zero potential (1-level)
%                    'auto'     - Automatic selection
%     'policy'    - Bound selection policy (default: 'vol')
%                   'vol', 'lu', 'l', 'u'
%     'maxiter'   - Maximum iterations (default: 100)
%     'tol'       - Convergence tolerance (default: 1e-4)
%     'verbose'   - Verbosity level 0/1/2 (default: 1)
%
% Output:
%   x     - Action rate parameters
%   pi    - Equilibrium distributions
%   Q     - Generator matrices
%   stats - Statistics struct
%
% Copyright (c) 2012-2025, Imperial College London
% All rights reserved.

%% Parse options
p = inputParser;
addParameter(p, 'relaxation', 'auto');
addParameter(p, 'policy', 'vol');
addParameter(p, 'maxiter', 100);
addParameter(p, 'tol', 1e-4);
addParameter(p, 'verbose', 1);
parse(p, varargin{:});
opts = p.Results;

RELAXATION = opts.relaxation;
POLICY = opts.policy;
MAXITER = opts.maxiter;
PFTOL = opts.tol;
GAPTOL = 0.01;
verbose = opts.verbose;

%% Extract model parameters
A = size(AP, 1);
ACT = AP(:, 1);
PSV = AP(:, 2);
M = max(AP(:));

% Extract rate matrices
Aa = cell(1, A);
Pb = cell(1, A);
for a = 1:A
    Aa{a} = R{a, 1};
    Pb{a} = R{a, 2};
end

L = cell(1, M);
for k = 1:M
    L{k} = R{A+1, k};
end

% State space sizes
N = zeros(1, M);
for k = 1:M
    N(k) = size(L{k}, 1);
end

%% Scale rates for numerical stability
scale = 1.0;
if ismember(RELAXATION, {'auto', 'tlpr0inf', 'tlpr1inf', 'tlprUinf', 'tzpr0inf', 'tzpr1inf', 'zpr', 'tma1inf'})
    for a = 1:A
        scale = max([scale, max(sum(abs(Aa{a}), 2))]);
    end
    for k = 1:M
        scale = max([scale, max(sum(abs(L{k}), 2))]);
    end
    scale = 1.1 * scale;
end

% Apply scaling
Aa_scaled = cell(1, A);
L_scaled = cell(1, M);
for a = 1:A
    Aa_scaled{a} = Aa{a} / scale;
end
for k = 1:M
    L_scaled{k} = L{k} / scale;
end

%% Initialize bounds
xL = zeros(A, 1);
xU = zeros(A, 1);
for a = 1:A
    xU(a) = max(Aa_scaled{a}(:)) * 1.1;
    xL(a) = max(PFTOL, min(Aa_scaled{a}(Aa_scaled{a} > 0)) * 0.1);
end

piL = cell(M, 1);
piU = cell(M, 1);
for k = 1:M
    piL{k} = zeros(1, N(k));
    piU{k} = ones(1, N(k));
end

%% Statistics
stats = struct();
stats.iter = 0;
stats.action_seq = [];
stats.bound_seq = [];
stats.tot_time = 0;
stats.relaxation_used = {};

%% Main iteration loop
gap = Inf(A, 1);
currentRelaxation = RELAXATION;

for iter = 1:MAXITER
    stats.iter = iter;

    % Select relaxation based on iteration (for 'auto' mode)
    if strcmpi(RELAXATION, 'auto')
        if iter <= A
            currentRelaxation = 'lpr';
        elseif iter <= 2*A
            currentRelaxation = 'lpr';
        elseif iter <= 3*A
            currentRelaxation = 'tlpr1inf';
        else
            currentRelaxation = 'tlpr1inf';
        end
    end

    % Select action and bound type using policy
    [targetAction, boundType] = select_action_bound(stats, A, xL, xU, POLICY);
    stats.action_seq(end+1) = targetAction;
    stats.bound_seq(end+1) = boundType;
    stats.relaxation_used{end+1} = currentRelaxation;

    T0 = tic;

    % Solve using selected relaxation
    piOpt = {};
    switch currentRelaxation
        case 'lpr'
            [xOpt, piOpt, fval, exitflag] = solve_lpr(targetAction, boundType, ...
                Aa_scaled, Pb, L_scaled, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL);
        case {'tlpr0inf', 'tlpr1inf', 'tlprUinf'}
            U = 1;  % Level for tightening
            if strcmpi(currentRelaxation, 'tlpr0inf')
                U = 0;
            elseif strcmpi(currentRelaxation, 'tlprUinf')
                U = min(2, floor(iter/A));
            end
            [xOpt, piOpt, fval, exitflag] = solve_tlpr(targetAction, boundType, ...
                Aa_scaled, Pb, L_scaled, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL, U);
        case 'ens'
            [xOpt, piOpt, fval, exitflag] = solve_ens(targetAction, boundType, ...
                Aa_scaled, Pb, L_scaled, ACT, PSV, M, A, N, xL, xU, PFTOL);
        case 'qcp'
            [xOpt, piOpt, fval, exitflag] = solve_qcp(...
                Aa_scaled, Pb, L_scaled, ACT, PSV, M, A, N, xL, xU, PFTOL);
        case 'ma'
            [xOpt, piOpt, fval, exitflag] = solve_ma(targetAction, boundType, ...
                Aa_scaled, Pb, L_scaled, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL);
        case 'va'
            [xOpt, piOpt, fval, exitflag] = solve_va(...
                Aa_scaled, Pb, L_scaled, ACT, PSV, M, A, N, xL, xU, PFTOL);
        case 'minres'
            [xOpt, piOpt, fval, exitflag] = solve_minres(targetAction, boundType, ...
                Aa_scaled, Pb, L_scaled, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL);
        case 'inap'
            [xOpt, piOpt, fval, exitflag] = solve_inap(...
                Aa_scaled, Pb, L_scaled, ACT, PSV, M, A, N, xL, xU, PFTOL);
        case 'tma1inf'
            [xOpt, piOpt, fval, exitflag] = solve_tma1inf(targetAction, boundType, ...
                Aa_scaled, Pb, L_scaled, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL);
        case 'zpr'
            [xOpt, piOpt, fval, exitflag] = solve_zpr(targetAction, boundType, ...
                Aa_scaled, Pb, L_scaled, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL);
        case {'tzpr0inf', 'tzpr1inf'}
            U = 1;  % Level for tightening
            if strcmpi(currentRelaxation, 'tzpr0inf')
                U = 0;
            end
            [xOpt, piOpt, fval, exitflag] = solve_tzpr(targetAction, boundType, ...
                Aa_scaled, Pb, L_scaled, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL, U);
        otherwise
            error('Unknown relaxation: %s', currentRelaxation);
    end

    stats.tot_time = stats.tot_time + toc(T0);

    if exitflag <= 0
        if verbose >= 1
            fprintf('iter %3d: x(%d) %c - infeasible (%s)\n', ...
                iter, targetAction, char('L' + boundType*('U'-'L')), currentRelaxation);
        end
        continue;
    end

    % Update x bounds
    newBnd = abs(xOpt(targetAction));
    if boundType == 0  % Lower bound
        oldBnd = xL(targetAction);
        xL(targetAction) = max(xL(targetAction), min(newBnd, xU(targetAction)));
    else  % Upper bound
        oldBnd = xU(targetAction);
        xU(targetAction) = min(xU(targetAction), max(newBnd, xL(targetAction)));
    end

    % Update pi bounds from LP solution (tightens McCormick envelopes)
    if ~isempty(piOpt)
        piTol = 0.1;  % Relaxed tolerance for pi bounds
        for k = 1:M
            if ~isempty(piOpt{k})
                piL{k} = max(piL{k}, piOpt{k} - piTol);
                piU{k} = min(piU{k}, piOpt{k} + piTol);
            end
        end
    end

    gap(targetAction) = abs(1 - xU(targetAction)/xL(targetAction));

    if verbose >= 1
        fprintf('iter %3d: x(%d) %c = %8.5f  gap = %8.5f  [%s]\n', ...
            iter, targetAction, char('L' + boundType*('U'-'L')), ...
            newBnd * scale, gap(targetAction), currentRelaxation);
    end

    % Check convergence
    if max(gap) < GAPTOL
        if verbose >= 1
            fprintf('\nConverged: max gap = %.6f\n', max(gap));
        end
        break;
    end
end

%% Compute final solution
x = (xL + xU) / 2 * scale;
[pi, Q] = compute_equilibrium(x/scale, Aa_scaled, Pb, L_scaled, ACT, PSV, M, A, N);

% Rescale Q matrices
for k = 1:M
    Q{k} = Q{k} * scale;
end

%% Check product-form conditions
rc3 = check_rc3(x/scale, pi, Aa_scaled, ACT, A, N) * scale;
if verbose >= 1
    fprintf('\nFinal solution:\n');
    for a = 1:A
        fprintf('  x(%d) = %.6f  [%.6f, %.6f]\n', a, x(a), xL(a)*scale, xU(a)*scale);
    end
    fprintf('  RC3 residual = %.6f\n', rc3);
    fprintf('  Total time = %.2f sec\n', stats.tot_time);
end

stats.xL = xL * scale;
stats.xU = xU * scale;
stats.gap = gap;
stats.rc3 = rc3;
stats.scale = scale;

end

%% Action/Bound Selection Policy
function [action, boundType] = select_action_bound(stats, A, xL, xU, policy)
    if isempty(stats.action_seq)
        action = 1;
        boundType = 0;
        return;
    end

    switch policy
        case 'vol'
            % Volume-based: prioritize largest gap, alternate L/U
            if stats.bound_seq(end) == 0
                action = stats.action_seq(end);
                boundType = 1;
            else
                [~, pos] = sort(xU - xL, 'descend');
                for a = pos(:)'
                    recent = stats.action_seq(max(1, end-2*(A-1)+1):end);
                    if isempty(find(a == recent, 1))
                        action = a;
                        boundType = 0;
                        return;
                    end
                end
                % Fallback
                action = mod(stats.action_seq(end), A) + 1;
                boundType = 0;
            end
        case 'lu'
            % Alternate lower/upper for each action
            if stats.bound_seq(end) == 0
                action = stats.action_seq(end);
                boundType = 1;
            else
                action = mod(stats.action_seq(end), A) + 1;
                boundType = 0;
            end
        case 'l'
            % Lower bounds only
            action = mod(stats.action_seq(end), A) + 1;
            boundType = 0;
        case 'u'
            % Upper bounds only
            action = mod(stats.action_seq(end), A) + 1;
            boundType = 1;
        otherwise
            action = mod(stats.action_seq(end), A) + 1;
            boundType = mod(stats.bound_seq(end) + 1, 2);
    end
end

%% LP Relaxation (lpr)
function [xOpt, piOpt, fval, exitflag] = solve_lpr(targetAction, boundType, ...
    Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL)

% Variables: [x(1:A); pi{1}(:); pi{2}(:); ...; z{1,1}(:); z{1,2}(:); ...]
nX = A;
nPi = sum(N);
nZ = 0;
for a = 1:A
    nZ = nZ + N(ACT(a)) + N(PSV(a));
end
nVar = nX + nPi + nZ;

% Variable indexing
xIdx = 1:A;
piIdx = cell(M, 1);
offset = A;
for k = 1:M
    piIdx{k} = offset + (1:N(k));
    offset = offset + N(k);
end
zIdx = cell(A, 2);
for a = 1:A
    zIdx{a, 1} = offset + (1:N(ACT(a)));
    offset = offset + N(ACT(a));
    zIdx{a, 2} = offset + (1:N(PSV(a)));
    offset = offset + N(PSV(a));
end

% Objective
f = zeros(nVar, 1);
if boundType == 0
    f(targetAction) = 1;
else
    f(targetAction) = -1;
end

% Build constraints
Aeq = [];
beq = [];
Aineq = [];
bineq = [];

% pi{k} * 1 = 1
for k = 1:M
    row = zeros(1, nVar);
    row(piIdx{k}) = 1;
    Aeq = [Aeq; row];
    beq = [beq; 1];
end

% sum(z{a,type}) = x(a)
for a = 1:A
    row = zeros(1, nVar);
    row(zIdx{a, 1}) = 1;
    row(xIdx(a)) = -1;
    Aeq = [Aeq; row];
    beq = [beq; 0];

    row = zeros(1, nVar);
    row(zIdx{a, 2}) = 1;
    row(xIdx(a)) = -1;
    Aeq = [Aeq; row];
    beq = [beq; 0];
end

% pi{k} * Q{k} = 0 (balance equations)
for k = 1:M
    Qlocal = L{k} - diag(sum(L{k}, 2));

    for s = 1:N(k)
        row = zeros(1, nVar);
        row(piIdx{k}) = Qlocal(:, s)';

        for a = 1:A
            if PSV(a) == k
                Qp = Pb{a} - diag(sum(Pb{a}, 2));
                row(zIdx{a, 2}) = row(zIdx{a, 2}) + Qp(:, s)';
            elseif ACT(a) == k
                Qa = Aa{a} - diag(sum(Aa{a}, 2));
                row(piIdx{k}) = row(piIdx{k}) + Qa(:, s)';
            end
        end

        Aineq = [Aineq; row; -row];
        bineq = [bineq; PFTOL; PFTOL];
    end
end

% McCormick envelopes
for a = 1:A
    k_act = ACT(a);
    k_psv = PSV(a);

    % Active z
    for i = 1:N(k_act)
        xLa = xL(a); xUa = xU(a);
        pL = piL{k_act}(i); pU = piU{k_act}(i);

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = -1;
        row(xIdx(a)) = pL;
        row(piIdx{k_act}(i)) = xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; pL * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(xIdx(a)) = -pL;
        row(piIdx{k_act}(i)) = -xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pL * xUa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(xIdx(a)) = -pU;
        row(piIdx{k_act}(i)) = -xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pU * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = -1;
        row(xIdx(a)) = pU;
        row(piIdx{k_act}(i)) = xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; pU * xUa];
    end

    % Passive z
    for i = 1:N(k_psv)
        xLa = xL(a); xUa = xU(a);
        pL = piL{k_psv}(i); pU = piU{k_psv}(i);

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = -1;
        row(xIdx(a)) = pL;
        row(piIdx{k_psv}(i)) = xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; pL * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = 1;
        row(xIdx(a)) = -pL;
        row(piIdx{k_psv}(i)) = -xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pL * xUa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = 1;
        row(xIdx(a)) = -pU;
        row(piIdx{k_psv}(i)) = -xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pU * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = -1;
        row(xIdx(a)) = pU;
        row(piIdx{k_psv}(i)) = xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; pU * xUa];
    end
end

% Bounds
lb = zeros(nVar, 1);
ub = Inf(nVar, 1);
lb(xIdx) = xL;
ub(xIdx) = xU;
for k = 1:M
    lb(piIdx{k}) = piL{k};
    ub(piIdx{k}) = piU{k};
end
for a = 1:A
    lb(zIdx{a, 1}) = 0;
    ub(zIdx{a, 1}) = xU(a);
    lb(zIdx{a, 2}) = 0;
    ub(zIdx{a, 2}) = xU(a);
end

% Solve LP
options = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'off');
[sol, fval, exitflag] = linprog(f, Aineq, bineq, Aeq, beq, lb, ub, options);

if exitflag > 0
    xOpt = sol(xIdx);
    piOpt = cell(M, 1);
    for k = 1:M
        piOpt{k} = sol(piIdx{k})';
    end
else
    xOpt = [];
    piOpt = {};
end

end

%% Tightened LP Relaxation (tlpr)
function [xOpt, piOpt, fval, exitflag] = solve_tlpr(targetAction, boundType, ...
    Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL, U)

% Start with basic LPR
[xOpt, piOpt, fval, exitflag] = solve_lpr(targetAction, boundType, ...
    Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL);

if exitflag <= 0
    return;
end

% Add tightening constraints iteratively
% The tightening uses pi*Q*A^u = 0 and pi*Q*(I-A)^{-1} = 0
% For simplicity, we iterate the basic LPR with updated bounds

for u = 1:U
    % Update pi bounds based on current solution
    if ~isempty(piOpt)
        for k = 1:M
            piL{k} = max(piL{k}, piOpt{k} - PFTOL);
            piU{k} = min(piU{k}, piOpt{k} + PFTOL);
        end
    end

    % Re-solve with tightened bounds
    [xOpt, piOpt, fval, exitflag] = solve_lpr(targetAction, boundType, ...
        Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL);

    if exitflag <= 0
        break;
    end
end

end

%% Exact Nonlinear System (ens) using fmincon
function [xOpt, piOpt, fval, exitflag] = solve_ens(targetAction, boundType, ...
    Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, PFTOL)

% ENS finds product-form solutions by minimizing RC3 violations
% Uses penalty method: minimize |pi*A - x*pi|^2 subject to normalization and balance

% Variables: [x(1:A); pi{1}(:); pi{2}(:); ...]
nX = A;
nPi = sum(N);
nVar = nX + nPi;

% Variable indexing
xIdx = 1:A;
piIdx = cell(M, 1);
offset = A;
for k = 1:M
    piIdx{k} = offset + (1:N(k));
    offset = offset + N(k);
end

% Better initial guess: use geometric distribution approximation
% For M/M/1 queue, pi(n) = (1-rho) * rho^n where rho = lambda/mu
x0 = zeros(nVar, 1);
x0(xIdx) = (xL + xU) / 2;
for k = 1:M
    % Estimate utilization from rate matrices
    rho_est = 0.5;  % Default estimate
    % Try to extract arrival/service rates from L
    arrivals = sum(diag(L{k}, 1));
    departures = sum(diag(L{k}, -1));
    if arrivals > 0 && departures > 0
        rho_est = min(0.95, arrivals / departures);
    end
    % Geometric distribution
    pi_init = (1 - rho_est) * rho_est.^(0:(N(k)-1));
    pi_init = pi_init / sum(pi_init);  % Normalize
    x0(piIdx{k}) = pi_init;
end

% Bounds
lb = zeros(nVar, 1);
ub = Inf(nVar, 1);
lb(xIdx) = xL;
ub(xIdx) = xU;
for k = 1:M
    lb(piIdx{k}) = 1e-12;  % Small positive to avoid log issues
    ub(piIdx{k}) = 1;
end

% Objective: minimize RC3 violation (squared sum of residuals)
% This is better conditioned than using as hard constraints
    function obj = objfun_penalty(v)
        xv = v(xIdx);
        rc3_penalty = 0;

        for kk = 1:M
            piv = v(piIdx{kk})';

            % RC3: pi * Aa = x(a) * pi for active processes
            for aa = 1:A
                if ACT(aa) == kk
                    residual = piv * Aa{aa} - xv(aa) * piv;
                    rc3_penalty = rc3_penalty + sum(residual.^2);
                end
            end
        end

        % Add optimization direction as small perturbation
        if boundType == 0  % minimizing
            obj = rc3_penalty + 1e-6 * v(targetAction);
        else  % maximizing
            obj = rc3_penalty - 1e-6 * v(targetAction);
        end
    end

% Equality constraints: only normalization and balance (not RC3)
    function [c, ceq] = nlcon(v)
        xv = v(xIdx);
        c = [];
        ceq = [];

        for kk = 1:M
            piv = v(piIdx{kk})';

            % Normalization: sum(pi) = 1
            ceq = [ceq; sum(piv) - 1];

            % Balance: pi * Q = 0 (only N-1 are independent, but fmincon handles this)
            Qk = L{kk} - diag(sum(L{kk}, 2));
            for aa = 1:A
                if PSV(aa) == kk
                    Qk = Qk + xv(aa) * (Pb{aa} - diag(sum(Pb{aa}, 2)));
                elseif ACT(aa) == kk
                    Qk = Qk + Aa{aa} - diag(sum(Aa{aa}, 2));
                end
            end
            % Use only first N-1 balance equations (last is dependent on normalization)
            balance = piv * Qk;
            ceq = [ceq; balance(1:end-1)'];
        end
    end

% Solve with relaxed tolerances and interior-point algorithm
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off', ...
    'MaxIterations', 2000, 'MaxFunctionEvaluations', 50000, ...
    'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-6, ...
    'StepTolerance', 1e-12);

[sol, fval, exitflag] = fmincon(@objfun_penalty, x0, [], [], [], [], lb, ub, @nlcon, options);

% Check if RC3 is satisfied (exitflag > 0 only indicates convergence)
if exitflag > 0 && fval < 1e-4  % RC3 residual is small enough
    xOpt = sol(xIdx);
    piOpt = cell(M, 1);
    for k = 1:M
        piOpt{k} = sol(piIdx{k})';
    end
elseif exitflag > 0
    % Converged but RC3 not satisfied - still return result for analysis
    xOpt = sol(xIdx);
    piOpt = cell(M, 1);
    for k = 1:M
        piOpt{k} = sol(piIdx{k})';
    end
else
    xOpt = [];
    piOpt = {};
end

end

%% Quadratically Constrained Program (qcp) using fmincon
function [xOpt, piOpt, fval, exitflag] = solve_qcp(...
    Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, PFTOL)

% QCP minimizes RC3 violations directly (no slack variables needed)
% This is a simpler reformulation that works better with fmincon

% Variables: [x(1:A); pi{1}(:); pi{2}(:); ...]
nX = A;
nPi = sum(N);
nVar = nX + nPi;

% Variable indexing
xIdx = 1:A;
piIdx = cell(M, 1);
offset = A;
for k = 1:M
    piIdx{k} = offset + (1:N(k));
    offset = offset + N(k);
end

% Better initial guess: use geometric distribution approximation
x0 = zeros(nVar, 1);
x0(xIdx) = (xL + xU) / 2;
for k = 1:M
    % Estimate utilization from rate matrices
    rho_est = 0.5;
    arrivals = sum(diag(L{k}, 1));
    departures = sum(diag(L{k}, -1));
    if arrivals > 0 && departures > 0
        rho_est = min(0.95, arrivals / departures);
    end
    pi_init = (1 - rho_est) * rho_est.^(0:(N(k)-1));
    pi_init = pi_init / sum(pi_init);
    x0(piIdx{k}) = pi_init;
end

% Bounds
lb = zeros(nVar, 1);
ub = Inf(nVar, 1);
lb(xIdx) = xL;
ub(xIdx) = xU;
for k = 1:M
    lb(piIdx{k}) = 1e-12;
    ub(piIdx{k}) = 1;
end

% Objective: minimize sum of squared RC3 violations (quadratic form)
    function obj = objfun_qcp(v)
        xv = v(xIdx);
        rc3_sum = 0;

        for kk = 1:M
            piv = v(piIdx{kk})';

            for aa = 1:A
                if ACT(aa) == kk
                    residual = piv * Aa{aa} - xv(aa) * piv;
                    rc3_sum = rc3_sum + sum(residual.^2);
                end
            end
        end
        obj = rc3_sum;
    end

% Equality constraints: normalization and balance only
    function [c, ceq] = nlcon(v)
        xv = v(xIdx);
        c = [];
        ceq = [];

        for kk = 1:M
            piv = v(piIdx{kk})';

            % Normalization: sum(pi) = 1
            ceq = [ceq; sum(piv) - 1];

            % Balance: pi * Q = 0
            Qk = L{kk} - diag(sum(L{kk}, 2));
            for aa = 1:A
                if PSV(aa) == kk
                    Qk = Qk + xv(aa) * (Pb{aa} - diag(sum(Pb{aa}, 2)));
                elseif ACT(aa) == kk
                    Qk = Qk + Aa{aa} - diag(sum(Aa{aa}, 2));
                end
            end
            % Use only first N-1 balance equations
            balance = piv * Qk;
            ceq = [ceq; balance(1:end-1)'];
        end
    end

% Solve with interior-point (better for quadratic objectives)
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off', ...
    'MaxIterations', 2000, 'MaxFunctionEvaluations', 50000, ...
    'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-6, ...
    'StepTolerance', 1e-12);

[sol, fval, exitflag] = fmincon(@objfun_qcp, x0, [], [], [], [], lb, ub, @nlcon, options);

if exitflag > 0
    xOpt = sol(xIdx);
    piOpt = cell(M, 1);
    for k = 1:M
        piOpt{k} = sol(piIdx{k})';
    end
else
    xOpt = [];
    piOpt = {};
end

end

%% Compute equilibrium
function [pi, Q] = compute_equilibrium(x, Aa, Pb, L, ACT, PSV, M, A, N)
Q = cell(1, M);
pi = cell(1, M);

for k = 1:M
    Qk = L{k} - diag(sum(L{k}, 2));

    for a = 1:A
        if PSV(a) == k
            Qk = Qk + x(a) * (Pb{a} - diag(sum(Pb{a}, 2)));
        elseif ACT(a) == k
            Qk = Qk + Aa{a} - diag(sum(Aa{a}, 2));
        end
    end

    Q{k} = Qk;
    Qk_inf = ctmc_makeinfgen(Qk);
    pi{k} = ctmc_solve(Qk_inf);
end
end

%% Check RC3 condition
function rc3 = check_rc3(x, pi, Aa, ACT, A, N)
rc3 = 0;
for a = 1:A
    k = ACT(a);
    residual = abs(pi{k} * Aa{a} - x(a) * pi{k});
    rc3 = rc3 + mean(residual);
end
rc3 = rc3 / A;
end

%% Mean Approximation (ma)
% Relaxes RC3 to only enforce first moment matching: sum(z) = sum(pi*Aa)
% This is equivalent to x = E[eigenvalue] from the active rate matrix
function [xOpt, piOpt, fval, exitflag] = solve_ma(targetAction, boundType, ...
    Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL)

% Variables: [x(1:A); pi{1}(:); pi{2}(:); ...; z{1,1}(:); z{1,2}(:); ...]
nX = A;
nPi = sum(N);
nZ = 0;
for a = 1:A
    nZ = nZ + N(ACT(a)) + N(PSV(a));
end
nVar = nX + nPi + nZ;

% Variable indexing
xIdx = 1:A;
piIdx = cell(M, 1);
offset = A;
for k = 1:M
    piIdx{k} = offset + (1:N(k));
    offset = offset + N(k);
end
zIdx = cell(A, 2);
for a = 1:A
    zIdx{a, 1} = offset + (1:N(ACT(a)));  % z_active
    offset = offset + N(ACT(a));
    zIdx{a, 2} = offset + (1:N(PSV(a)));  % z_passive
    offset = offset + N(PSV(a));
end

% Objective
f = zeros(nVar, 1);
if boundType == 0
    f(targetAction) = 1;
else
    f(targetAction) = -1;
end

% Build constraints
Aeq = [];
beq = [];
Aineq = [];
bineq = [];

% pi{k} * 1 = 1 (normalization)
for k = 1:M
    row = zeros(1, nVar);
    row(piIdx{k}) = 1;
    Aeq = [Aeq; row];
    beq = [beq; 1];
end

% Mean approximation: sum(z{a,active}) = sum(pi*Aa) = x(a) * sum(pi) = x(a)
% and sum(z{a,passive}) = x(a)
for a = 1:A
    % z_active sums to x
    row = zeros(1, nVar);
    row(zIdx{a, 1}) = 1;
    row(xIdx(a)) = -1;
    Aeq = [Aeq; row];
    beq = [beq; 0];

    % z_passive sums to x
    row = zeros(1, nVar);
    row(zIdx{a, 2}) = 1;
    row(xIdx(a)) = -1;
    Aeq = [Aeq; row];
    beq = [beq; 0];
end

% Balance equations: pi{k} * Q{k} = 0 (relaxed to tolerance)
for k = 1:M
    Qlocal = L{k} - diag(sum(L{k}, 2));

    for s = 1:N(k)
        row = zeros(1, nVar);
        row(piIdx{k}) = Qlocal(:, s)';

        for a = 1:A
            if PSV(a) == k
                Qp = Pb{a} - diag(sum(Pb{a}, 2));
                row(zIdx{a, 2}) = row(zIdx{a, 2}) + Qp(:, s)';
            elseif ACT(a) == k
                Qa = Aa{a} - diag(sum(Aa{a}, 2));
                row(piIdx{k}) = row(piIdx{k}) + Qa(:, s)';
            end
        end

        Aineq = [Aineq; row; -row];
        bineq = [bineq; PFTOL; PFTOL];
    end
end

% McCormick envelopes for z = x * pi
for a = 1:A
    k_act = ACT(a);
    k_psv = PSV(a);

    % Active z
    for i = 1:N(k_act)
        xLa = xL(a); xUa = xU(a);
        pL = piL{k_act}(i); pU = piU{k_act}(i);

        % z >= pL*x + xL*pi - pL*xL
        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = -1;
        row(xIdx(a)) = pL;
        row(piIdx{k_act}(i)) = xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; pL * xLa];

        % z <= pL*x + xU*pi - pL*xU
        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(xIdx(a)) = -pL;
        row(piIdx{k_act}(i)) = -xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pL * xUa];

        % z <= pU*x + xL*pi - pU*xL
        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(xIdx(a)) = -pU;
        row(piIdx{k_act}(i)) = -xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pU * xLa];

        % z >= pU*x + xU*pi - pU*xU
        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = -1;
        row(xIdx(a)) = pU;
        row(piIdx{k_act}(i)) = xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; pU * xUa];
    end

    % Passive z
    for i = 1:N(k_psv)
        xLa = xL(a); xUa = xU(a);
        pL = piL{k_psv}(i); pU = piU{k_psv}(i);

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = -1;
        row(xIdx(a)) = pL;
        row(piIdx{k_psv}(i)) = xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; pL * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = 1;
        row(xIdx(a)) = -pL;
        row(piIdx{k_psv}(i)) = -xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pL * xUa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = 1;
        row(xIdx(a)) = -pU;
        row(piIdx{k_psv}(i)) = -xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pU * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = -1;
        row(xIdx(a)) = pU;
        row(piIdx{k_psv}(i)) = xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; pU * xUa];
    end
end

% Relaxed RC3 (mean approximation): sum(z*I - pi*Aa) = 0
% This only enforces first moment, not pointwise equality
for a = 1:A
    k = ACT(a);
    % sum(z{a,active}) - sum(pi * Aa) should be close to 0
    % but sum(z{a,active}) = x(a), so this is: x(a) - sum(pi*Aa) close to 0
    row = zeros(1, nVar);
    row(xIdx(a)) = 1;
    row(piIdx{k}) = -sum(Aa{a}, 2)';
    Aineq = [Aineq; row; -row];
    bineq = [bineq; PFTOL; PFTOL];
end

% Bounds
lb = zeros(nVar, 1);
ub = Inf(nVar, 1);
lb(xIdx) = xL;
ub(xIdx) = xU;
for k = 1:M
    lb(piIdx{k}) = piL{k};
    ub(piIdx{k}) = piU{k};
end
for a = 1:A
    lb(zIdx{a, 1}) = 0;
    ub(zIdx{a, 1}) = xU(a);
    lb(zIdx{a, 2}) = 0;
    ub(zIdx{a, 2}) = xU(a);
end

% Solve LP
options = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'off');
[sol, fval, exitflag] = linprog(f, Aineq, bineq, Aeq, beq, lb, ub, options);

if exitflag > 0
    xOpt = sol(xIdx);
    piOpt = cell(M, 1);
    for k = 1:M
        piOpt{k} = sol(piIdx{k})';
    end
else
    xOpt = [];
    piOpt = {};
end

end

%% Variance Approximation (va)
% Uses fixed x from bounds, solves for pi with slack variables on RC3
function [xOpt, piOpt, fval, exitflag] = solve_va(...
    Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, PFTOL)

% Use midpoint of bounds as fixed x
x_fixed = (xL + xU) / 2;

% Variables: [pi{1}(:); pi{2}(:); ...; slack{1}; slack{2}; ...]
nPi = sum(N);
nSlack = 0;
for a = 1:A
    nSlack = nSlack + N(ACT(a));
end
nVar = nPi + nSlack;

% Variable indexing
piIdx = cell(M, 1);
offset = 0;
for k = 1:M
    piIdx{k} = offset + (1:N(k));
    offset = offset + N(k);
end
slackIdx = cell(A, 1);
for a = 1:A
    slackIdx{a} = offset + (1:N(ACT(a)));
    offset = offset + N(ACT(a));
end

% Objective: minimize sum of absolute slacks (L1 norm)
f = zeros(nVar, 1);
for a = 1:A
    f(slackIdx{a}) = 1;  % Minimize slack
end

% Build constraints
Aeq = [];
beq = [];
Aineq = [];
bineq = [];

% pi{k} * 1 = 1 (normalization)
for k = 1:M
    row = zeros(1, nVar);
    row(piIdx{k}) = 1;
    Aeq = [Aeq; row];
    beq = [beq; 1];
end

% Balance equations: pi{k} * Q{k} = 0
for k = 1:M
    Qk = L{k} - diag(sum(L{k}, 2));
    for a = 1:A
        if PSV(a) == k
            Qk = Qk + x_fixed(a) * (Pb{a} - diag(sum(Pb{a}, 2)));
        elseif ACT(a) == k
            Qk = Qk + Aa{a} - diag(sum(Aa{a}, 2));
        end
    end

    for s = 1:(N(k)-1)  % Only N-1 are independent
        row = zeros(1, nVar);
        row(piIdx{k}) = Qk(:, s)';
        Aeq = [Aeq; row];
        beq = [beq; 0];
    end
end

% RC3 with slack: pi*Aa - x*pi + slack >= 0, pi*Aa - x*pi - slack <= 0
% This means |pi*Aa - x*pi| <= slack
for a = 1:A
    k = ACT(a);
    for i = 1:N(k)
        % (pi*Aa - x*pi)_i + slack_i >= 0
        row = zeros(1, nVar);
        row(piIdx{k}) = Aa{a}(:, i)' - x_fixed(a) * (1:N(k) == i);
        row(slackIdx{a}(i)) = 1;
        Aineq = [Aineq; -row];  % -(...) <= 0
        bineq = [bineq; 0];

        % (pi*Aa - x*pi)_i - slack_i <= 0
        row = zeros(1, nVar);
        row(piIdx{k}) = Aa{a}(:, i)' - x_fixed(a) * (1:N(k) == i);
        row(slackIdx{a}(i)) = -1;
        Aineq = [Aineq; row];
        bineq = [bineq; 0];
    end
end

% Bounds
lb = zeros(nVar, 1);
ub = Inf(nVar, 1);
for k = 1:M
    lb(piIdx{k}) = 0;
    ub(piIdx{k}) = 1;
end
for a = 1:A
    lb(slackIdx{a}) = 0;  % Slack must be non-negative
end

% Solve LP
options = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'off');
[sol, fval, exitflag] = linprog(f, Aineq, bineq, Aeq, beq, lb, ub, options);

if exitflag > 0
    xOpt = x_fixed;  % Return the fixed x
    piOpt = cell(M, 1);
    for k = 1:M
        piOpt{k} = sol(piIdx{k})';
    end
else
    xOpt = [];
    piOpt = {};
end

end

%% Minimum Residual (minres)
% Minimizes L1 norm of RC3 violations with McCormick envelopes
function [xOpt, piOpt, fval, exitflag] = solve_minres(targetAction, boundType, ...
    Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL)

% Variables: [x; pi{1}(:); ...; z{1,1}(:); z{1,2}(:); ...; slack_pos; slack_neg]
nX = A;
nPi = sum(N);
nZ = 0;
for a = 1:A
    nZ = nZ + N(ACT(a)) + N(PSV(a));
end
nSlack = 0;
for a = 1:A
    nSlack = nSlack + 2 * N(ACT(a));  % Positive and negative slack for each RC3 component
end
nVar = nX + nPi + nZ + nSlack;

% Variable indexing
xIdx = 1:A;
piIdx = cell(M, 1);
offset = A;
for k = 1:M
    piIdx{k} = offset + (1:N(k));
    offset = offset + N(k);
end
zIdx = cell(A, 2);
for a = 1:A
    zIdx{a, 1} = offset + (1:N(ACT(a)));
    offset = offset + N(ACT(a));
    zIdx{a, 2} = offset + (1:N(PSV(a)));
    offset = offset + N(PSV(a));
end
slackPosIdx = cell(A, 1);
slackNegIdx = cell(A, 1);
for a = 1:A
    slackPosIdx{a} = offset + (1:N(ACT(a)));
    offset = offset + N(ACT(a));
    slackNegIdx{a} = offset + (1:N(ACT(a)));
    offset = offset + N(ACT(a));
end

% Objective: minimize L1 norm of slack (minimize RC3 residuals)
f = zeros(nVar, 1);
for a = 1:A
    f(slackPosIdx{a}) = 1;
    f(slackNegIdx{a}) = 1;
end

% Build constraints
Aeq = [];
beq = [];
Aineq = [];
bineq = [];

% pi{k} * 1 = 1
for k = 1:M
    row = zeros(1, nVar);
    row(piIdx{k}) = 1;
    Aeq = [Aeq; row];
    beq = [beq; 1];
end

% sum(z{a,type}) = x(a)
for a = 1:A
    row = zeros(1, nVar);
    row(zIdx{a, 1}) = 1;
    row(xIdx(a)) = -1;
    Aeq = [Aeq; row];
    beq = [beq; 0];

    row = zeros(1, nVar);
    row(zIdx{a, 2}) = 1;
    row(xIdx(a)) = -1;
    Aeq = [Aeq; row];
    beq = [beq; 0];
end

% Balance equations: pi{k} * Q{k} = 0
for k = 1:M
    Qlocal = L{k} - diag(sum(L{k}, 2));

    for s = 1:N(k)
        row = zeros(1, nVar);
        row(piIdx{k}) = Qlocal(:, s)';

        for a = 1:A
            if PSV(a) == k
                Qp = Pb{a} - diag(sum(Pb{a}, 2));
                row(zIdx{a, 2}) = row(zIdx{a, 2}) + Qp(:, s)';
            elseif ACT(a) == k
                Qa = Aa{a} - diag(sum(Aa{a}, 2));
                row(piIdx{k}) = row(piIdx{k}) + Qa(:, s)';
            end
        end

        Aineq = [Aineq; row; -row];
        bineq = [bineq; PFTOL; PFTOL];
    end
end

% McCormick envelopes for z = x * pi
for a = 1:A
    k_act = ACT(a);
    k_psv = PSV(a);

    for i = 1:N(k_act)
        xLa = xL(a); xUa = xU(a);
        pL = piL{k_act}(i); pU = piU{k_act}(i);

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = -1;
        row(xIdx(a)) = pL;
        row(piIdx{k_act}(i)) = xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; pL * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(xIdx(a)) = -pL;
        row(piIdx{k_act}(i)) = -xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pL * xUa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(xIdx(a)) = -pU;
        row(piIdx{k_act}(i)) = -xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pU * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = -1;
        row(xIdx(a)) = pU;
        row(piIdx{k_act}(i)) = xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; pU * xUa];
    end

    for i = 1:N(k_psv)
        xLa = xL(a); xUa = xU(a);
        pL = piL{k_psv}(i); pU = piU{k_psv}(i);

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = -1;
        row(xIdx(a)) = pL;
        row(piIdx{k_psv}(i)) = xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; pL * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = 1;
        row(xIdx(a)) = -pL;
        row(piIdx{k_psv}(i)) = -xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pL * xUa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = 1;
        row(xIdx(a)) = -pU;
        row(piIdx{k_psv}(i)) = -xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pU * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = -1;
        row(xIdx(a)) = pU;
        row(piIdx{k_psv}(i)) = xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; pU * xUa];
    end
end

% RC3 with slack: z_i - (pi*Aa)_i = slack_pos_i - slack_neg_i
% We use z{a,active} as the bilinear term x*pi
for a = 1:A
    k = ACT(a);
    for i = 1:N(k)
        % z_i - (pi*Aa)_i - slack_pos + slack_neg = 0
        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(piIdx{k}) = -Aa{a}(:, i)';
        row(slackPosIdx{a}(i)) = -1;
        row(slackNegIdx{a}(i)) = 1;
        Aeq = [Aeq; row];
        beq = [beq; 0];
    end
end

% Bounds
lb = zeros(nVar, 1);
ub = Inf(nVar, 1);
lb(xIdx) = xL;
ub(xIdx) = xU;
for k = 1:M
    lb(piIdx{k}) = piL{k};
    ub(piIdx{k}) = piU{k};
end
for a = 1:A
    lb(zIdx{a, 1}) = 0;
    ub(zIdx{a, 1}) = xU(a);
    lb(zIdx{a, 2}) = 0;
    ub(zIdx{a, 2}) = xU(a);
    lb(slackPosIdx{a}) = 0;
    lb(slackNegIdx{a}) = 0;
end

% Solve LP
options = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'off');
[sol, fval, exitflag] = linprog(f, Aineq, bineq, Aeq, beq, lb, ub, options);

if exitflag > 0
    xOpt = sol(xIdx);
    piOpt = cell(M, 1);
    for k = 1:M
        piOpt{k} = sol(piIdx{k})';
    end
else
    xOpt = [];
    piOpt = {};
end

end

%% Iterative Fixed-Point Approximation (inap)
% Uses fixed-point iteration: x(a) = mean(pi*Aa ./ pi)
function [xOpt, piOpt, fval, exitflag] = solve_inap(...
    Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, PFTOL)

% Initialize x randomly within bounds
x = xL + (xU - xL) .* rand(A, 1);

% Compute initial equilibrium
[pi, Q] = compute_eq_internal(x, Aa, Pb, L, ACT, PSV, M, A, N);

% Fixed-point iteration
maxIter = 50;
tol = 1e-6;

for iter = 1:maxIter
    xprev = x;

    for a = 1:A
        k = ACT(a);
        % LAMBDA = (pi * Aa) ./ pi (element-wise)
        piAa = pi{k} * Aa{a};
        LAMBDA = piAa ./ pi{k};
        LAMBDA = LAMBDA(~isnan(LAMBDA) & ~isinf(LAMBDA) & LAMBDA > 0);

        if ~isempty(LAMBDA)
            x(a) = mean(LAMBDA);
            % Clamp to bounds
            x(a) = max(xL(a), min(xU(a), x(a)));
        end
    end

    % Recompute equilibrium
    [pi, Q] = compute_eq_internal(x, Aa, Pb, L, ACT, PSV, M, A, N);

    % Check convergence
    if norm(x - xprev, inf) < tol
        break;
    end
end

% Check if RC3 is satisfied
rc3_satisfied = true;
for a = 1:A
    k = ACT(a);
    piAa = pi{k} * Aa{a};
    LAMBDA = piAa ./ pi{k};
    LAMBDA = LAMBDA(~isnan(LAMBDA) & ~isinf(LAMBDA) & LAMBDA > 0);
    if ~isempty(LAMBDA) && (max(LAMBDA) - min(LAMBDA)) > 0.1
        rc3_satisfied = false;
    end
end

if rc3_satisfied
    xOpt = x;
    piOpt = pi;
    fval = 0;
    exitflag = 1;
else
    xOpt = x;
    piOpt = pi;
    fval = 1;  % Indicates approximate solution
    exitflag = 1;  % Still return solution for use
end

end

%% Internal equilibrium computation for inap
function [pi, Q] = compute_eq_internal(x, Aa, Pb, L, ACT, PSV, M, A, N)
Q = cell(1, M);
pi = cell(1, M);

for k = 1:M
    Qk = L{k} - diag(sum(L{k}, 2));

    for a = 1:A
        if PSV(a) == k
            Qk = Qk + x(a) * (Pb{a} - diag(sum(Pb{a}, 2)));
        elseif ACT(a) == k
            Qk = Qk + Aa{a} - diag(sum(Aa{a}, 2));
        end
    end

    Q{k} = Qk;
    Qk_inf = ctmc_makeinfgen(Qk);
    pi{k} = ctmc_solve(Qk_inf);
end
end

%% Tightened Mean Approximation with 1-level + infinity (tma1inf)
% Combines mean approximation (RC3 first moment only) with tightened LP constraints
function [xOpt, piOpt, fval, exitflag] = solve_tma1inf(targetAction, boundType, ...
    Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL)

% First run mean approximation to get initial solution
[xMA, piMA, ~, exitflagMA] = solve_ma(targetAction, boundType, ...
    Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL);

if exitflagMA <= 0
    xOpt = [];
    piOpt = {};
    fval = Inf;
    exitflag = -1;
    return;
end

% Use MA solution to warm-start tightened LP with additional constraints
% Variables: [x(1:A); pi{1}(:); pi{2}(:); ...; z{1,1}(:); z{1,2}(:); ...]
nX = A;
nPi = sum(N);
nZ = 0;
for a = 1:A
    nZ = nZ + N(ACT(a)) + N(PSV(a));
end
nVar = nX + nPi + nZ;

% Variable indexing
xIdx = 1:A;
piIdx = cell(M, 1);
offset = A;
for k = 1:M
    piIdx{k} = offset + (1:N(k));
    offset = offset + N(k);
end
zIdx = cell(A, 2);
for a = 1:A
    zIdx{a, 1} = offset + (1:N(ACT(a)));
    offset = offset + N(ACT(a));
    zIdx{a, 2} = offset + (1:N(PSV(a)));
    offset = offset + N(PSV(a));
end

% Objective
f = zeros(nVar, 1);
if boundType == 0
    f(targetAction) = 1;
else
    f(targetAction) = -1;
end

% Build constraints
Aeq = [];
beq = [];
Aineq = [];
bineq = [];

% pi{k} * 1 = 1
for k = 1:M
    row = zeros(1, nVar);
    row(piIdx{k}) = 1;
    Aeq = [Aeq; row];
    beq = [beq; 1];
end

% RC3 first moment constraint: sum(z{a,active}) = x(a) (mean approximation)
for a = 1:A
    row = zeros(1, nVar);
    row(zIdx{a, 1}) = 1;
    row(xIdx(a)) = -1;
    Aeq = [Aeq; row];
    beq = [beq; 0];

    row = zeros(1, nVar);
    row(zIdx{a, 2}) = 1;
    row(xIdx(a)) = -1;
    Aeq = [Aeq; row];
    beq = [beq; 0];
end

% pi{k} * Q{k} = 0 (balance equations)
for k = 1:M
    Qlocal = L{k} - diag(sum(L{k}, 2));

    for s = 1:N(k)
        row = zeros(1, nVar);
        row(piIdx{k}) = Qlocal(:, s)';

        for a = 1:A
            if PSV(a) == k
                Qp = Pb{a} - diag(sum(Pb{a}, 2));
                row(zIdx{a, 2}) = row(zIdx{a, 2}) + Qp(:, s)';
            elseif ACT(a) == k
                Qa = Aa{a} - diag(sum(Aa{a}, 2));
                row(piIdx{k}) = row(piIdx{k}) + Qa(:, s)';
            end
        end

        Aineq = [Aineq; row; -row];
        bineq = [bineq; PFTOL; PFTOL];
    end
end

% Tightened LP constraints: pi*Q*Aa^u = 0 for u=0,1,...,U and pi*Q*Aa^inf = 0
U = 1;  % Tightening level
for k = 1:M
    for b = 1:A
        if ACT(b) == k
            Qlocal = L{k} - diag(sum(L{k}, 2));

            % Add constraints for u = 1
            for u = 1:U
                Aa_pow = Aa{b}^u;

                for s = 1:N(k)
                    row = zeros(1, nVar);

                    % z{b,active} * Aa^(u-1) * Q_local
                    coeffZ = (Aa{b}^(u-1)) * Qlocal;
                    row(zIdx{b, 1}) = coeffZ(:, s)';

                    for c = 1:A
                        if PSV(c) == k
                            Qp = Pb{c} - diag(sum(Pb{c}, 2));
                            coeffZp = Aa_pow * Qp;
                            row(zIdx{c, 2}) = row(zIdx{c, 2}) + coeffZp(:, s)';
                        elseif ACT(c) == k
                            Qa = Aa{c} - diag(sum(Aa{c}, 2));
                            coeffZa = (Aa{b}^(u-1)) * Qa;
                            row(zIdx{b, 1}) = row(zIdx{b, 1}) + coeffZa(:, s)';
                        end
                    end

                    Aineq = [Aineq; row; -row];
                    bineq = [bineq; PFTOL; PFTOL];
                end
            end

            % Add infinity constraint: z * (I-Aa)^{-1} * Q = 0
            try
                Aa_inv = inv(eye(N(k)) - Aa{b});

                for s = 1:N(k)
                    row = zeros(1, nVar);

                    coeffZ = Aa_inv * Qlocal;
                    row(zIdx{b, 1}) = coeffZ(:, s)';

                    for c = 1:A
                        if PSV(c) == k
                            Qp = Pb{c} - diag(sum(Pb{c}, 2));
                            coeffZp = Aa{b} * Aa_inv * Qp;
                            row(zIdx{c, 2}) = row(zIdx{c, 2}) + coeffZp(:, s)';
                        elseif ACT(c) == k
                            Qa = Aa{c} - diag(sum(Aa{c}, 2));
                            coeffZa = Aa_inv * Qa;
                            row(zIdx{b, 1}) = row(zIdx{b, 1}) + coeffZa(:, s)';
                        end
                    end

                    Aineq = [Aineq; row; -row];
                    bineq = [bineq; PFTOL; PFTOL];
                end
            catch
                % Matrix not invertible, skip infinity constraints
            end
        end
    end
end

% McCormick envelopes
for a = 1:A
    k_act = ACT(a);
    k_psv = PSV(a);

    % Active z
    for i = 1:N(k_act)
        xLa = xL(a); xUa = xU(a);
        pL = piL{k_act}(i); pU = piU{k_act}(i);

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = -1;
        row(xIdx(a)) = pL;
        row(piIdx{k_act}(i)) = xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; pL * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(xIdx(a)) = -pL;
        row(piIdx{k_act}(i)) = -xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pL * xUa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(xIdx(a)) = -pU;
        row(piIdx{k_act}(i)) = -xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pU * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = -1;
        row(xIdx(a)) = pU;
        row(piIdx{k_act}(i)) = xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; pU * xUa];
    end

    % Passive z
    for i = 1:N(k_psv)
        xLa = xL(a); xUa = xU(a);
        pL = piL{k_psv}(i); pU = piU{k_psv}(i);

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = -1;
        row(xIdx(a)) = pL;
        row(piIdx{k_psv}(i)) = xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; pL * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = 1;
        row(xIdx(a)) = -pL;
        row(piIdx{k_psv}(i)) = -xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pL * xUa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = 1;
        row(xIdx(a)) = -pU;
        row(piIdx{k_psv}(i)) = -xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pU * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = -1;
        row(xIdx(a)) = pU;
        row(piIdx{k_psv}(i)) = xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; pU * xUa];
    end
end

% Bounds
lb = zeros(nVar, 1);
ub = Inf(nVar, 1);
lb(xIdx) = xL;
ub(xIdx) = xU;
for k = 1:M
    lb(piIdx{k}) = piL{k};
    ub(piIdx{k}) = piU{k};
end
for a = 1:A
    lb(zIdx{a, 1}) = 0;
    ub(zIdx{a, 1}) = xU(a);
    lb(zIdx{a, 2}) = 0;
    ub(zIdx{a, 2}) = xU(a);
end

% Solve LP
options = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'off');
[sol, fval, exitflag] = linprog(f, Aineq, bineq, Aeq, beq, lb, ub, options);

if exitflag > 0
    xOpt = sol(xIdx);
    piOpt = cell(M, 1);
    for k = 1:M
        piOpt{k} = sol(piIdx{k})';
    end
else
    xOpt = [];
    piOpt = {};
end

end

%% Zero Potential Relaxation (zpr)
% Uses potential matrix to add cutting plane constraints
function [xOpt, piOpt, fval, exitflag] = solve_zpr(targetAction, boundType, ...
    Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL)

% First get initial solution using inap or from midpoint
x0 = (xL + xU) / 2;
[pi0, Q0] = compute_eq_zpr(x0, Aa, Pb, L, ACT, PSV, M, A, N);

% Compute potential matrices using Jacobi iteration
g = cell(M, max(N));
for k = 1:M
    for n = 1:N(k)
        f_vec = zeros(N(k), 1);
        f_vec(n) = 1;
        g{k, n} = potential_jacobi(Q0{k}, pi0{k}, f_vec);
    end
end

% Variables: [x(1:A); pi{1}(:); pi{2}(:); ...; z{1,1}(:); z{1,2}(:); ...]
nX = A;
nPi = sum(N);
nZ = 0;
for a = 1:A
    nZ = nZ + N(ACT(a)) + N(PSV(a));
end
nVar = nX + nPi + nZ;

% Variable indexing
xIdx = 1:A;
piIdx = cell(M, 1);
offset = A;
for k = 1:M
    piIdx{k} = offset + (1:N(k));
    offset = offset + N(k);
end
zIdx = cell(A, 2);
for a = 1:A
    zIdx{a, 1} = offset + (1:N(ACT(a)));
    offset = offset + N(ACT(a));
    zIdx{a, 2} = offset + (1:N(PSV(a)));
    offset = offset + N(PSV(a));
end

% Objective
f = zeros(nVar, 1);
if boundType == 0
    f(targetAction) = 1;
else
    f(targetAction) = -1;
end

% Build constraints
Aeq = [];
beq = [];
Aineq = [];
bineq = [];

% pi{k} * 1 = 1
for k = 1:M
    row = zeros(1, nVar);
    row(piIdx{k}) = 1;
    Aeq = [Aeq; row];
    beq = [beq; 1];
end

% sum(z{a,type}) = x(a)
for a = 1:A
    row = zeros(1, nVar);
    row(zIdx{a, 1}) = 1;
    row(xIdx(a)) = -1;
    Aeq = [Aeq; row];
    beq = [beq; 0];

    row = zeros(1, nVar);
    row(zIdx{a, 2}) = 1;
    row(xIdx(a)) = -1;
    Aeq = [Aeq; row];
    beq = [beq; 0];
end

% pi{k} * Q{k} = 0 (balance equations)
for k = 1:M
    Qlocal = L{k} - diag(sum(L{k}, 2));

    for s = 1:N(k)
        row = zeros(1, nVar);
        row(piIdx{k}) = Qlocal(:, s)';

        for a = 1:A
            if PSV(a) == k
                Qp = Pb{a} - diag(sum(Pb{a}, 2));
                row(zIdx{a, 2}) = row(zIdx{a, 2}) + Qp(:, s)';
            elseif ACT(a) == k
                Qa = Aa{a} - diag(sum(Aa{a}, 2));
                row(piIdx{k}) = row(piIdx{k}) + Qa(:, s)';
            end
        end

        Aineq = [Aineq; row; -row];
        bineq = [bineq; PFTOL; PFTOL];
    end
end

% RC3 constraint: z{a,active}*I - pi*Aa = 0
for a = 1:A
    k = ACT(a);
    for i = 1:N(k)
        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(piIdx{k}) = -Aa{a}(:, i)';

        Aineq = [Aineq; row; -row];
        bineq = [bineq; PFTOL; PFTOL];
    end
end

% Potential constraints (cutting planes from zpr)
for k = 1:M
    for n = 1:N(k)
        % delta constraint: pi0 - pi + correction terms = 0
        row = zeros(1, nVar);
        row(piIdx{k}(n)) = -1;
        rhs = -pi0{k}(n);

        for b = 1:A
            if PSV(b) == k
                Qp = Pb{b} - diag(sum(Pb{b}, 2));
                % z{b,passive} * Qp * g{k,n}
                coeff = Qp * g{k, n};
                row(zIdx{b, 2}) = row(zIdx{b, 2}) + coeff';
                % -x0(b) * pi * Qp * g{k,n}
                coeff2 = Qp * g{k, n};
                row(piIdx{k}) = row(piIdx{k}) - x0(b) * coeff2';
            end
        end

        Aineq = [Aineq; row; -row];
        bineq = [bineq; rhs + PFTOL; -rhs + PFTOL];
    end
end

% McCormick envelopes
for a = 1:A
    k_act = ACT(a);
    k_psv = PSV(a);

    for i = 1:N(k_act)
        xLa = xL(a); xUa = xU(a);
        pL = piL{k_act}(i); pU = piU{k_act}(i);

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = -1;
        row(xIdx(a)) = pL;
        row(piIdx{k_act}(i)) = xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; pL * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(xIdx(a)) = -pL;
        row(piIdx{k_act}(i)) = -xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pL * xUa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(xIdx(a)) = -pU;
        row(piIdx{k_act}(i)) = -xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pU * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = -1;
        row(xIdx(a)) = pU;
        row(piIdx{k_act}(i)) = xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; pU * xUa];
    end

    for i = 1:N(k_psv)
        xLa = xL(a); xUa = xU(a);
        pL = piL{k_psv}(i); pU = piU{k_psv}(i);

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = -1;
        row(xIdx(a)) = pL;
        row(piIdx{k_psv}(i)) = xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; pL * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = 1;
        row(xIdx(a)) = -pL;
        row(piIdx{k_psv}(i)) = -xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pL * xUa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = 1;
        row(xIdx(a)) = -pU;
        row(piIdx{k_psv}(i)) = -xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pU * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = -1;
        row(xIdx(a)) = pU;
        row(piIdx{k_psv}(i)) = xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; pU * xUa];
    end
end

% Bounds
lb = zeros(nVar, 1);
ub = Inf(nVar, 1);
lb(xIdx) = xL;
ub(xIdx) = xU;
for k = 1:M
    lb(piIdx{k}) = piL{k};
    ub(piIdx{k}) = piU{k};
end
for a = 1:A
    lb(zIdx{a, 1}) = 0;
    ub(zIdx{a, 1}) = xU(a);
    lb(zIdx{a, 2}) = 0;
    ub(zIdx{a, 2}) = xU(a);
end

% Solve LP
options = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'off');
[sol, fval, exitflag] = linprog(f, Aineq, bineq, Aeq, beq, lb, ub, options);

if exitflag > 0
    xOpt = sol(xIdx);
    piOpt = cell(M, 1);
    for k = 1:M
        piOpt{k} = sol(piIdx{k})';
    end
else
    xOpt = [];
    piOpt = {};
end

end

%% Tightened Zero Potential Relaxation (tzpr0inf/tzpr1inf)
% Combines zero potential relaxation with tightened LP constraints
% U = tightening level (0 for tzpr0inf, 1 for tzpr1inf)
function [xOpt, piOpt, fval, exitflag] = solve_tzpr(targetAction, boundType, ...
    Aa, Pb, L, ACT, PSV, M, A, N, xL, xU, piL, piU, PFTOL, U)

% First get initial solution from midpoint
x0 = (xL + xU) / 2;
[pi0, Q0] = compute_eq_zpr(x0, Aa, Pb, L, ACT, PSV, M, A, N);

% Compute potential matrices
g = cell(M, max(N));
for k = 1:M
    for n = 1:N(k)
        f_vec = zeros(N(k), 1);
        f_vec(n) = 1;
        g{k, n} = potential_jacobi(Q0{k}, pi0{k}, f_vec);
    end
end

% Variables: [x(1:A); pi{1}(:); pi{2}(:); ...; z{1,1}(:); z{1,2}(:); ...]
nX = A;
nPi = sum(N);
nZ = 0;
for a = 1:A
    nZ = nZ + N(ACT(a)) + N(PSV(a));
end
nVar = nX + nPi + nZ;

% Variable indexing
xIdx = 1:A;
piIdx = cell(M, 1);
offset = A;
for k = 1:M
    piIdx{k} = offset + (1:N(k));
    offset = offset + N(k);
end
zIdx = cell(A, 2);
for a = 1:A
    zIdx{a, 1} = offset + (1:N(ACT(a)));
    offset = offset + N(ACT(a));
    zIdx{a, 2} = offset + (1:N(PSV(a)));
    offset = offset + N(PSV(a));
end

% Objective
f = zeros(nVar, 1);
if boundType == 0
    f(targetAction) = 1;
else
    f(targetAction) = -1;
end

% Build constraints
Aeq = [];
beq = [];
Aineq = [];
bineq = [];

% pi{k} * 1 = 1
for k = 1:M
    row = zeros(1, nVar);
    row(piIdx{k}) = 1;
    Aeq = [Aeq; row];
    beq = [beq; 1];
end

% sum(z{a,type}) = x(a)
for a = 1:A
    row = zeros(1, nVar);
    row(zIdx{a, 1}) = 1;
    row(xIdx(a)) = -1;
    Aeq = [Aeq; row];
    beq = [beq; 0];

    row = zeros(1, nVar);
    row(zIdx{a, 2}) = 1;
    row(xIdx(a)) = -1;
    Aeq = [Aeq; row];
    beq = [beq; 0];
end

% pi{k} * Q{k} = 0 (balance equations)
for k = 1:M
    Qlocal = L{k} - diag(sum(L{k}, 2));

    for s = 1:N(k)
        row = zeros(1, nVar);
        row(piIdx{k}) = Qlocal(:, s)';

        for a = 1:A
            if PSV(a) == k
                Qp = Pb{a} - diag(sum(Pb{a}, 2));
                row(zIdx{a, 2}) = row(zIdx{a, 2}) + Qp(:, s)';
            elseif ACT(a) == k
                Qa = Aa{a} - diag(sum(Aa{a}, 2));
                row(piIdx{k}) = row(piIdx{k}) + Qa(:, s)';
            end
        end

        Aineq = [Aineq; row; -row];
        bineq = [bineq; PFTOL; PFTOL];
    end
end

% RC3 constraint
for a = 1:A
    k = ACT(a);
    for i = 1:N(k)
        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(piIdx{k}) = -Aa{a}(:, i)';

        Aineq = [Aineq; row; -row];
        bineq = [bineq; PFTOL; PFTOL];
    end
end

% Potential constraints (from zpr)
for k = 1:M
    for n = 1:N(k)
        row = zeros(1, nVar);
        row(piIdx{k}(n)) = -1;
        rhs = -pi0{k}(n);

        for b = 1:A
            if PSV(b) == k
                Qp = Pb{b} - diag(sum(Pb{b}, 2));
                coeff = Qp * g{k, n};
                row(zIdx{b, 2}) = row(zIdx{b, 2}) + coeff';
                coeff2 = Qp * g{k, n};
                row(piIdx{k}) = row(piIdx{k}) - x0(b) * coeff2';
            end
        end

        Aineq = [Aineq; row; -row];
        bineq = [bineq; rhs + PFTOL; -rhs + PFTOL];
    end
end

% Tightened LP constraints (from tlpr)
% U is passed as parameter: 0 for tzpr0inf, 1 for tzpr1inf
for k = 1:M
    for b = 1:A
        if ACT(b) == k
            Qlocal = L{k} - diag(sum(L{k}, 2));

            for u = 1:U
                Aa_pow = Aa{b}^u;

                for s = 1:N(k)
                    row = zeros(1, nVar);
                    coeffZ = (Aa{b}^(u-1)) * Qlocal;
                    row(zIdx{b, 1}) = coeffZ(:, s)';

                    for c = 1:A
                        if PSV(c) == k
                            Qp = Pb{c} - diag(sum(Pb{c}, 2));
                            coeffZp = Aa_pow * Qp;
                            row(zIdx{c, 2}) = row(zIdx{c, 2}) + coeffZp(:, s)';
                        elseif ACT(c) == k
                            Qa = Aa{c} - diag(sum(Aa{c}, 2));
                            coeffZa = (Aa{b}^(u-1)) * Qa;
                            row(zIdx{b, 1}) = row(zIdx{b, 1}) + coeffZa(:, s)';
                        end
                    end

                    Aineq = [Aineq; row; -row];
                    bineq = [bineq; PFTOL; PFTOL];
                end
            end

            % Infinity constraint
            try
                Aa_inv = inv(eye(N(k)) - Aa{b});

                for s = 1:N(k)
                    row = zeros(1, nVar);
                    coeffZ = Aa_inv * Qlocal;
                    row(zIdx{b, 1}) = coeffZ(:, s)';

                    for c = 1:A
                        if PSV(c) == k
                            Qp = Pb{c} - diag(sum(Pb{c}, 2));
                            coeffZp = Aa{b} * Aa_inv * Qp;
                            row(zIdx{c, 2}) = row(zIdx{c, 2}) + coeffZp(:, s)';
                        elseif ACT(c) == k
                            Qa = Aa{c} - diag(sum(Aa{c}, 2));
                            coeffZa = Aa_inv * Qa;
                            row(zIdx{b, 1}) = row(zIdx{b, 1}) + coeffZa(:, s)';
                        end
                    end

                    Aineq = [Aineq; row; -row];
                    bineq = [bineq; PFTOL; PFTOL];
                end
            catch
                % Matrix not invertible
            end
        end
    end
end

% McCormick envelopes
for a = 1:A
    k_act = ACT(a);
    k_psv = PSV(a);

    for i = 1:N(k_act)
        xLa = xL(a); xUa = xU(a);
        pL = piL{k_act}(i); pU = piU{k_act}(i);

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = -1;
        row(xIdx(a)) = pL;
        row(piIdx{k_act}(i)) = xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; pL * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(xIdx(a)) = -pL;
        row(piIdx{k_act}(i)) = -xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pL * xUa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = 1;
        row(xIdx(a)) = -pU;
        row(piIdx{k_act}(i)) = -xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pU * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 1}(i)) = -1;
        row(xIdx(a)) = pU;
        row(piIdx{k_act}(i)) = xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; pU * xUa];
    end

    for i = 1:N(k_psv)
        xLa = xL(a); xUa = xU(a);
        pL = piL{k_psv}(i); pU = piU{k_psv}(i);

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = -1;
        row(xIdx(a)) = pL;
        row(piIdx{k_psv}(i)) = xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; pL * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = 1;
        row(xIdx(a)) = -pL;
        row(piIdx{k_psv}(i)) = -xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pL * xUa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = 1;
        row(xIdx(a)) = -pU;
        row(piIdx{k_psv}(i)) = -xLa;
        Aineq = [Aineq; row];
        bineq = [bineq; -pU * xLa];

        row = zeros(1, nVar);
        row(zIdx{a, 2}(i)) = -1;
        row(xIdx(a)) = pU;
        row(piIdx{k_psv}(i)) = xUa;
        Aineq = [Aineq; row];
        bineq = [bineq; pU * xUa];
    end
end

% Bounds
lb = zeros(nVar, 1);
ub = Inf(nVar, 1);
lb(xIdx) = xL;
ub(xIdx) = xU;
for k = 1:M
    lb(piIdx{k}) = piL{k};
    ub(piIdx{k}) = piU{k};
end
for a = 1:A
    lb(zIdx{a, 1}) = 0;
    ub(zIdx{a, 1}) = xU(a);
    lb(zIdx{a, 2}) = 0;
    ub(zIdx{a, 2}) = xU(a);
end

% Solve LP
options = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'off');
[sol, fval, exitflag] = linprog(f, Aineq, bineq, Aeq, beq, lb, ub, options);

if exitflag > 0
    xOpt = sol(xIdx);
    piOpt = cell(M, 1);
    for k = 1:M
        piOpt{k} = sol(piIdx{k})';
    end
else
    xOpt = [];
    piOpt = {};
end

end

%% Compute equilibrium for zpr
function [pi, Q] = compute_eq_zpr(x, Aa, Pb, L, ACT, PSV, M, A, N)
Q = cell(1, M);
pi = cell(1, M);

for k = 1:M
    Qk = L{k} - diag(sum(L{k}, 2));

    for a = 1:A
        if PSV(a) == k
            Qk = Qk + x(a) * (Pb{a} - diag(sum(Pb{a}, 2)));
        elseif ACT(a) == k
            Qk = Qk + Aa{a} - diag(sum(Aa{a}, 2));
        end
    end

    Q{k} = Qk;
    Qk_inf = ctmc_makeinfgen(Qk);
    pi{k} = ctmc_solve(Qk_inf);
end
end

%% Potential matrix computation using Jacobi iteration
function g = potential_jacobi(Q, pi, f)
% Compute potential g such that Q*g = f - pi*f*e
% Uses Jacobi iteration method

n = length(pi);
tol = 1e-10;
maxIter = 1000;

% Right-hand side
b = f - (pi * f) * ones(n, 1);

% Initialize g
g = zeros(n, 1);

% Extract diagonal
d = diag(Q);

% Jacobi iteration
for iter = 1:maxIter
    g_old = g;

    for i = 1:n
        if abs(d(i)) > 1e-12
            sum_term = Q(i, :) * g_old - d(i) * g_old(i);
            g(i) = (b(i) - sum_term) / d(i);
        end
    end

    if norm(g - g_old, inf) < tol
        break;
    end
end

% Normalize to satisfy pi*g = 0
g = g - (pi * g);

end
