function [QN, UN, RN, TN, CN, XN, iter] = solver_ag(sn, options)
% SOLVER_AG RCAT agent-based methods for SolverAG
%
% [QN, UN, RN, TN, CN, XN, ITER] = SOLVER_AG(SN, OPTIONS)
%
% Uses RCAT (Reversed Compound Agent Theorem) to find product-form
% solutions for queueing networks.
%
% Each (station, class) pair becomes an isolated component, and the components
% are coupled only through the reversed rates of the synchronizing actions. A
% component is a QBD whose LEVEL is the queue length and whose PHASE is the
% pair (arrival phase, service phase), laid out in the Kronecker order of
% qbd_mapmap1: an arrival moves the level up carrying kron(D1^a, I), a service
% completion moves it down carrying kron(I, D1^s), the busy levels evolve under
% krons(D0^a, D0^s) and level zero under kron(D0^a, I), because no server is
% running there. With exponential processes every block is 1x1 and the QBD
% collapses to the scalar birth-death chain this analyzer used before, entry
% for entry.
%
% Methods:
%   'inap'     - Iterative Numerical Approximation Procedure (default, fast)
%   'inapplus' - Improved INAP with weighted rates (no normalization)
%   'inapinf'  - INAP with matrix-geometric solution of the isolated open
%                components (no state-space truncation), per Marin, Rota Bulo,
%                Balsamo, "A Numerical Algorithm for the Decomposition of
%                Cooperating Structured Markov Processes", MASCOTS 2012.
%   'exact'    - Not available (see solver_ag_autocat, unreachable)
%
% Copyright (c) 2012-2025, Imperial College London
% All rights reserved.

M = sn.nstations;
K = sn.nclasses;

% Set default max states for truncation
if isfield(options, 'config') && isfield(options.config, 'maxStates')
    maxStates = options.config.maxStates;
else
    maxStates = 100;
end

% Set default tolerances
if isfield(options, 'iter_tol') && ~isempty(options.iter_tol)
    tol = options.iter_tol;
else
    tol = 1e-6;
end

if isfield(options, 'iter_max') && ~isempty(options.iter_max)
    maxiter = options.iter_max;
else
    maxiter = 1000;
end

% Infinite-server stations: the single-server component is exact only under
% the rule of AG_INF_SUPPORTS, which SolverAG.supportsModelMethod also asks.
[infOk, infWhy] = ag_inf_supports(sn);
if ~infOk
    line_error(mfilename, infWhy);
end

% Build RCAT model from network structure
[R, AP, processMap, actionMap, N, meta] = build_rcat(sn, maxStates);

% Execution backend of the fixed point. The agents are solved in ISOLATION
% given the reversed rates, so the sweep is order-free and every backend walks
% the same iterates; 'parallel' and 'cluster' change who evaluates an agent, not
% what the agent evaluates to. See ag_exec_resolve for the contract.
meta.exec = ag_exec_resolve(options);
if strcmp(meta.exec.mode, 'cluster') && strcmpi(options.method, 'inapinf')
    % The remote worker implements the FINITE agent solve. 'inapinf' replaces it
    % with the matrix-geometric treatment of an open agent -- Neuts' R matrix and
    % the scalar-tail detection that precedes it -- which the worker does not
    % carry, and answering with the finite solve instead would silently change
    % the method. Refuse by name rather than substitute.
    line_error(mfilename, ['The ''cluster'' execution backend does not carry the ' ...
        '''inapinf'' agent solve (the matrix-geometric tail of an open agent runs ' ...
        'on the coordinator only). Use exec ''serial'' or ''parallel'' with ' ...
        '''inapinf'', or method ''inap''/''inapplus'' with ''cluster''.']);
end

% Check if we have a valid model
numProcesses = max(processMap(:));
numActions = size(AP, 1);

% Return early only if no processes found
if numProcesses == 0
    line_warning(mfilename, 'Network could not be mapped to RCAT format (no processes found).\n');
    QN = zeros(M, K);
    UN = zeros(M, K);
    RN = zeros(M, K);
    TN = zeros(M, K);
    CN = zeros(1, K);
    XN = zeros(1, K);
    iter = 0;
    return;
end

% If no actions but we have processes, solve using local rates only
% This handles single-queue G-networks (Source -> Queue -> Sink)
if numActions == 0
    % No inter-station actions: solve equilibrium using only L matrices
    x = [];
    pi = cell(1, numProcesses);
    Q = cell(1, numProcesses);
    for p = 1:numProcesses
        L = R{1, p};  % Local rate matrix is in R{numActions+1, p} = R{1, p} when numActions=0
        % Convert to valid generator matrix
        Qp = L - diag(L * ones(size(L, 1), 1));
        Q{p} = ctmc_makeinfgen(Qp);
        % Solve for equilibrium through the same dispatcher the fixed point
        % uses: this branch carries a whole M/PH/1 on its own, whose marginal
        % spans tens of orders of magnitude over the truncation, and the level
        % recursions are stable there where a null-space solve is not.
        pi{p} = ag_solve_component(Q{p}, meta.mph(p), meta.nlev(p), meta.level{p});
    end
    iter = 0;
    [QN, UN, RN, TN, CN, XN] = rcat_metrics(sn, x, pi, Q, processMap, actionMap, N, meta);
    return;
end

% Open/closed flag per process (open classes have infinite population and
% are the ones the matrix-geometric 'inapinf' method solves without truncation).
isOpenProc = false(1, numProcesses);
for p = 1:numProcesses
    [ipst, ipr] = find(processMap == p);
    if ~isempty(ipst)
        isOpenProc(p) = isinf(sn.njobs(ipr(1)));
    end
end

% Choose solver method
method = options.method;
if strcmp(method, 'default')
    method = 'inap';
end

% Per-process geometric-tail decay (set only by 'inapinf'); empty => metrics
% are computed from the explicit stationary vectors pi.
rhoProc = [];
isGeomProc = [];
geomData = {};

switch method
    case 'inap'
        % Fast iterative heuristic
        [x, pi, Q, iter] = inap(R, AP, meta, tol, maxiter, 'inap');

    case 'inapplus'
        % Improved INAP with weighted rates (no normalization)
        [x, pi, Q, iter] = inap(R, AP, meta, tol, maxiter, 'inapplus');

    case 'inapinf'
        % Matrix-geometric INAP: solve isolated open components exactly on the
        % infinite state space (geometric tail), no truncation.
        [x, pi, Q, iter, rhoProc, isGeomProc, geomData, rcatRes] = ...
            inap_inf(R, AP, meta, tol, maxiter, isOpenProc);
        line_debug('inapinf: RCAT product-form residual = %.3e (iter=%d)', rcatRes, iter);

    case 'exact'
        % Optimization-based solver using autocat (not available in this version)
        line_warning(mfilename, '''exact'' method not available. Falling back to inap.\n');
        [x, pi, Q, iter] = inap(R, AP, meta, tol, maxiter, 'inap');

    otherwise
        line_error(mfilename, 'Unknown method: %s\n', method);
end

% Convert RCAT solution to LINE metrics
[QN, UN, RN, TN, CN, XN] = rcat_metrics(sn, x, pi, Q, processMap, actionMap, N, meta, ...
    rhoProc, isGeomProc, geomData);

end

%% Local Functions

function [x, pi, Q, iter] = inap(R, AP, meta, tol, maxiter, method)
% INAP Iterative Numerical Approximation Procedure for RCAT
%
% Methods:
%   'inap':     x(a) = mean over the support of (pi Aa)_j / pi_j
%   'inapplus': x(a) = sum_ij Aa(i,j) pi(i)

if nargin < 6 || isempty(method)
    method = 'inap';
end

% Parse R and AP
A = size(AP, 1);
ACT = AP(:, 1);
PSV = AP(:, 2);
numProcesses = max(AP(:));

% Extract rate matrices
Aa = cell(1, A);
Pb = cell(1, A);
for a = 1:A
    Aa{a} = R{a, 1};
    Pb{a} = R{a, 2};
end

% Extract local rates
L = cell(1, numProcesses);
for k = 1:numProcesses
    L{k} = R{A+1, k};
end

% Get state space sizes
N = zeros(1, numProcesses);
for k = 1:numProcesses
    N(k) = size(L{k}, 1);
end

% notBirthDeath selects INAP+ rate-conservation estimator (catastrophe/batch
% removal) vs INAP mean-of-ratios; see _kb/06-solver-catalog.md for rationale.
% The test is on the LEVEL distance, not the state distance: with a phase block
% per level the within-level phase transitions of a PH sit far off the
% diagonal and are not a departure from birth-death structure.
notBirthDeath = false(1, numProcesses);
for k = 1:numProcesses
    lvl = meta.level{k};
    for n = 1:N(k)
        for m = 1:N(k)
            if L{k}(n, m) > 0 && abs(lvl(n) - lvl(m)) > 1
                notBirthDeath(k) = true;
            end
        end
    end
    % A PHASE-EXPANDED component takes the same estimator, for the same reason.
    % On a birth-death chain every state-wise reversed rate equals lambda, so
    % their mean is exact; with a phase block per level they do not, the deep
    % truncation levels dominate the unweighted mean, and the mean-of-ratios
    % overestimates the departure rate exactly as it does on a catastrophe
    % (measured on a tandem with Erlang(2) service at Q1: the reversed rate came
    % out 1.27 against the exact 0.5, so flow was not conserved). Rate
    % conservation has no such failure mode.
    if meta.mph(k) > 1
        notBirthDeath(k) = true;
    end
end

% Columns of each active matrix that carry any rate. Aa does not depend on x,
% so this is fixed for the whole fixed point.
activeCols = cell(1, A);
for a = 1:A
    activeCols{a} = find(sum(Aa{a}, 1) > 0);
end

% Deterministic initial guess (reproducibility); see _kb/06-solver-catalog.md for rationale
x = (1:A)' / (A + 1);

% Compute initial equilibrium
[pi, Q] = compute_equilibrium(x, Aa, Pb, L, ACT, PSV, numProcesses, A, N, meta);

% reversed-rate fixed point on the isolated-component equilibria, driven
% by the generic DA driver
fpopts = struct('iter_max', maxiter, 'iter_tol', tol);
fpopts.config.da_norm = @pi_blocknorm;
[~, iter, cvg] = da_fpi(@inap_sweep, pi, fpopts);
if ~cvg
    iter = iter + 1; % legacy while-loop exited with the counter past the cap
end

    function [xnew, xref] = inap_sweep(picur, ~)
    xref = picur;

    % Update each action rate
    for a = 1:A
        k = ACT(a);
        v = pi{k}(:)' * Aa{a};

        if strcmp(method, 'inapplus') || notBirthDeath(k)
            % inapplus: x(a) = sum_ij Aa(i,j) pi(i), the departure rate of the
            % active component.
            LAMBDA_sum = sum(v);
            if LAMBDA_sum > 0
                x(a) = LAMBDA_sum;
            end
        else
            % inap: x(a) = mean over the support of the STATE-WISE reversed rate
            % (pi Aa)_j / pi_j, which RCAT requires to be independent of j. On a
            % birth-death component every column of Aa holds one entry, so this
            % is the reference's entrywise mean of Aa(i,j) pi(i) / pi(j) term
            % for term; with a phase block per level a column holds one entry
            % per phase, and only the column form is the reversed rate.
            cols = activeCols{a};
            num = v(cols);
            den = pi{k}(cols);
            ok = den(:)' > 0 & num > 0;
            if any(ok)
                ratio = num(ok) ./ den(ok);
                ratio = ratio(isfinite(ratio));
                if ~isempty(ratio)
                    x(a) = mean(ratio);
                end
            end
        end
    end

    % Recompute equilibrium with new x
    [pi, Q] = compute_equilibrium(x, Aa, Pb, L, ACT, PSV, numProcesses, A, N, meta);
    xnew = pi;
    end

    function e = pi_blocknorm(xn, xr)
    e = 0;
    for kk = 1:numProcesses
        e = max(e, norm(xn{kk} - xr{kk}, 1));
    end
    end

end

function [pi, Q] = compute_equilibrium(x, Aa, Pb, L, ACT, PSV, numProcesses, A, N, meta)
% Compute equilibrium distribution for each process given action rates x.
%
% Agent k is solved in ISOLATION: its generator reads the rest of the model
% only through the scalar reversed rates x, and it writes only its own slot of
% pi and Q. That is what makes the loop a fan-out rather than a recurrence, so
% the three backends below differ in WHO evaluates agent_solve and not in what
% it returns. Keep it that way: any cross-agent read added here would silently
% make 'parallel' and 'cluster' race, and the identity assertions in
% test_ag_exec_backends are the only thing that would catch it.

Q = cell(1, numProcesses);
pi = cell(1, numProcesses);

exec = meta.exec;
switch exec.mode
    case 'parallel'
        parfor k = 1:numProcesses
            [Q{k}, pi{k}] = agent_solve(k, x, Aa, Pb, L, ACT, PSV, A, N, meta);
        end
    case 'cluster'
        [pi, Q] = ag_cluster_sweep(x, Aa, Pb, L, ACT, PSV, numProcesses, A, N, meta);
    otherwise
        for k = 1:numProcesses
            [Q{k}, pi{k}] = agent_solve(k, x, Aa, Pb, L, ACT, PSV, A, N, meta);
        end
end

end

function [Qk, pik] = agent_solve(k, x, Aa, Pb, L, ACT, PSV, A, N, meta)
% One agent's generator and stationary vector, given the reversed rates.
%
% The single point where an agent is evaluated. Every execution backend routes
% through it -- serial, parfor and the remote worker alike -- so there is one
% definition of what an agent's answer is and the backends cannot drift.

Qk = ag_agent_generator(k, x, Aa, Pb, L, ACT, PSV, A, N);
pik = ag_solve_component(Qk, meta.mph(k), meta.nlev(k), meta.level{k});

end

function [x, pi, Q, iter, rhoProc, isGeomProc, geomData, rcatRes] = ...
    inap_inf(R, AP, meta, tol, maxiter, isOpenProc)
% INAP_QBD Matrix-geometric INAP for RCAT product forms (no truncation).
%
% Same fixed-point iteration over the reversed rates x_l as INAP, but each
% isolated OPEN component is solved directly on its infinite state space. With
% one phase per level that is a scalar matrix-geometric (QBD / catastrophe)
% decomposition: the marginal is geometric pi_n = (1-rho) rho^n with rho the
% sub-unit root of the QBD characteristic equation, and any catastrophe drain
% to the empty state is folded into the local outflow (it produces no interior
% inflow, so the geometric form is preserved). With a phase block per level it
% is Neuts' rate matrix R from qbd_R, pi_(n+1) = pi_n R, closed by the boundary
% equations of the level-0 and level-1 blocks. Closed components remain finite
% and are solved on the explicit state space. Reversed rates are updated by the
% weighted-mean formula Eq. (4) evaluated in closed form on the tail, and the
% RCAT product-form residual (Remark 2) is returned as a diagnostic.
%
% Reference: A. Marin, S. Rota Bulo, S. Balsamo, "A Numerical Algorithm for
% the Decomposition of Cooperating Structured Markov Processes", MASCOTS 2012.

A = size(AP, 1);
ACT = AP(:, 1);
PSV = AP(:, 2);
numProcesses = max(AP(:));

% Extract rate matrices
Aa = cell(1, A);
Pb = cell(1, A);
for a = 1:A
    Aa{a} = R{a, 1};
    Pb{a} = R{a, 2};
end

% Extract local rates
L = cell(1, numProcesses);
for k = 1:numProcesses
    L{k} = R{A+1, k};
end

% State space sizes and active-transition row sums (rate of the active label
% out of each state of the active component)
N = zeros(1, numProcesses);
for k = 1:numProcesses
    N(k) = size(L{k}, 1);
end
aRowSum = cell(1, A);
for a = 1:A
    aRowSum{a} = sum(Aa{a}, 2);
end

% Deterministic initial guess (see inap): reproducible across back-ends.
x = (1:A)' / (A + 1);

[pi, Q, rhoProc, isGeomProc, geomData] = ...
    compute_equilibrium_qbd(x, Aa, Pb, L, ACT, PSV, numProcesses, A, N, meta, isOpenProc);

% reversed-rate fixed point on the isolated-component equilibria (matrix-
% geometric variant), driven by the generic DA driver
fpopts = struct('iter_max', maxiter, 'iter_tol', tol);
fpopts.config.da_norm = @pi_blocknorm_trunc;
[~, iter, cvg] = da_fpi(@inapinf_sweep, pi, fpopts);
if ~cvg
    iter = iter + 1; % legacy while-loop exited with the counter past the cap
end
rcat_residual();

    function [xnew, xref] = inapinf_sweep(picur, ~)
    xref = picur;

    % Reversed-rate update, Eq. (4): x_l = pi^(alpha_l) T^(l) e.
    for a = 1:A
        k = ACT(a);
        if isGeomProc(k)
            if meta.mph(k) == 1
                % Geometric tail: the active label fires only in occupied
                % states, so x_l = (per-occupied-state active rate) * P(occupied).
                occ = aRowSum{a}(min(2, N(k)));
                x(a) = occ * rhoProc(k);
            else
                % Matrix-geometric tail: sum_{n>=1} pi_n = pi_1 (I - R)^-1, and
                % the active label has the same row sums at every busy level.
                occ = aRowSum{a}(ag_blk(1, meta.mph(k)));
                x(a) = geomData{k}.busy * occ(:);
            end
        else
            v = pi{k}(:);
            x(a) = v' * aRowSum{a};
        end
    end

    [pi, Q, rhoProc, isGeomProc, geomData] = ...
        compute_equilibrium_qbd(x, Aa, Pb, L, ACT, PSV, numProcesses, A, N, meta, isOpenProc);
    xnew = pi;
    end

    function e = pi_blocknorm_trunc(xn, xr)
    e = 0;
    for kk = 1:numProcesses
        m = min(length(xn{kk}), length(xr{kk}));
        e = max(e, norm(xn{kk}(1:m) - xr{kk}(1:m), 1));
    end
    end

    function rcat_residual()
    % RCAT product-form residual (Remark 2): max_l || pi^(alpha_l) (x_l I - T^(l)) ||,
    % where T^(l) is the active rate matrix Aa{a}. Zero iff the reversed rate is
    % state-independent, i.e. an exact product-form solution was found.
    rcatRes = 0;
    for a = 1:A
        k = ACT(a);
        v = pi{k}(:)';
        resVec = x(a) * v - v * Aa{a};
        rcatRes = max(rcatRes, norm(resVec, 2));
    end
    end

end

function [pi, Q, rhoProc, isGeomProc, geomData] = ...
    compute_equilibrium_qbd(x, Aa, Pb, L, ACT, PSV, numProcesses, A, N, meta, isOpenProc)
% Solve isolated components (open: matrix-geometric; closed: finite solve);
% see _kb/06-solver-catalog.md for rationale
Q = cell(1, numProcesses);
pi = cell(1, numProcesses);
rhoProc = zeros(1, numProcesses);
isGeomProc = false(1, numProcesses);
geomData = cell(1, numProcesses);

exec = meta.exec;
switch exec.mode
    case 'parallel'
        parfor k = 1:numProcesses
            [Q{k}, pi{k}, rhoProc(k), isGeomProc(k), geomData{k}] = ...
                agent_solve_qbd(k, x, Aa, Pb, L, ACT, PSV, A, N, meta, isOpenProc);
        end
    otherwise
        for k = 1:numProcesses
            [Q{k}, pi{k}, rhoProc(k), isGeomProc(k), geomData{k}] = ...
                agent_solve_qbd(k, x, Aa, Pb, L, ACT, PSV, A, N, meta, isOpenProc);
        end
end

end

function [Qk_out, pik, rho_out, isGeom, geom] = ...
    agent_solve_qbd(k, x, Aa, Pb, L, ACT, PSV, A, N, meta, isOpenProc)
% One agent under the matrix-geometric ('inapinf') treatment: its generator,
% its stationary vector and, when the agent is open and its tail is genuinely
% geometric, the decay that stands in for the truncated levels.
%
% The 'inapinf' twin of agent_solve, and the same rule holds: an agent reads
% the rest of the model only through x, so every execution backend routes
% through this one definition.

rho_out = 0;
isGeom = false;
geom = [];
pik = [];
Nk = N(k);
mph = meta.mph(k);
nlev = meta.nlev(k);

% Assemble strictly off-diagonal rate matrix for component k
Off = L{k} - diag(diag(L{k}));
for c = 1:A
    if PSV(c) == k
        Off = Off + x(c) * Pb{c};
    elseif ACT(c) == k
        Off = Off + Aa{c};
    end
end
Off = Off - diag(diag(Off));

Qk = Off - diag(sum(Off, 2));
Qk_out = ctmc_makeinfgen(Qk);

solvedGeom = false;
if isOpenProc(k) && nlev >= 5
    if mph == 1
        % Read the homogeneous interior rates one level below the truncation
        % boundary (avoids the reflecting boundary artefact of Off).
        s0 = Nk - 1;                 % interior state index (level s0-1)
        row = Off(s0, :);
        f  = row(s0 + 1);            % up-1  rate (arrival)
        b  = row(s0 - 1);            % down-1 rate (service + single removal)
        g0 = row(1);                 % drain to empty state (catastrophe)
        % Transitions to strictly-interior lower levels (batch removal to a
        % non-empty state) break the scalar-QBD structure; detect and defer.
        interDown = 0;
        if s0 - 2 >= 2
            interDown = sum(row(2:s0-2));
        end
        if interDown <= 1e-11 && f > 0
            rho = qbd_scalar_rho(f, b, g0);
            if isfinite(rho) && rho > 0 && rho < 1 - 1e-12
                rho_out = rho;
                isGeom = true;
                pik = (1 - rho) * rho .^ (0:Nk-1);
                solvedGeom = true;
            end
        end
    elseif ag_is_block_tridiagonal(Qk_out, meta.level{k})
        % Read the homogeneous interior blocks one level below the
        % truncation boundary, for the same reason the scalar branch reads
        % the interior row there.
        s0 = nlev - 2;
        A0 = Qk_out(ag_blk(s0, mph), ag_blk(s0 + 1, mph));   % up: arrival
        A1 = Qk_out(ag_blk(s0, mph), ag_blk(s0, mph));       % local, with diagonal
        A2 = Qk_out(ag_blk(s0, mph), ag_blk(s0 - 1, mph));   % down: departure
        if any(sum(A0, 2) > 0)
            g = qbd_matrix_tail(Qk_out, A0, A1, A2, mph);
            if ~isempty(g)
                isGeom = true;
                geom = g;
                pik = qbd_tail_expand(g, nlev, mph);
                solvedGeom = true;
            end
        end
    end
end

if ~solvedGeom
    pik = ag_solve_component(Qk_out, mph, nlev, meta.level{k});
end

end

function g = qbd_matrix_tail(Qk, A0, A1, A2, mph)
% Neuts' matrix-geometric solution of one open component with MPH phases per
% level: R from qbd_R, then the boundary equations of levels 0 and 1,
%   pi_0 B00 + pi_1 A2 = 0,   pi_0 A0 + pi_1 (A1 + R A2) = 0,
% normalized by pi_0 e + pi_1 (I - R)^-1 e = 1. Returns [] when R has no
% sub-unit spectral radius, i.e. when the isolated component is unstable and
% has no stationary tail to report.

g = [];
% Logarithmic reduction rather than successive substitutions: this runs once
% per component per fixed-point sweep, and the quadratic convergence is what
% keeps that affordable. The cap is small for the same reason -- logarithmic
% reduction needs a few dozen steps when the component is stable, and an
% unstable one must bail rather than grind to the library default of 1e5.
R = qbd_R_logred(A2, A1, A0, 500);
% The minimal solution of a QBD is NON-NEGATIVE; anything else is the iteration
% having failed rather than a rate matrix.
if any(~isfinite(R(:))) || any(R(:) < -1e-12)
    return;
end

IR = eye(mph) - R;
tailMass = IR \ ones(mph, 1);
% STABILITY WITHOUT AN EIGENSOLVER. (I-R)^-1 = I + R + R^2 + ... converges
% exactly when the spectral radius is below one, and every row of that series is
% e_i plus non-negative terms, so (I-R)^-1 e >= 1 entrywise. When the isolated
% component is unstable the series diverges and the inverse picks up negative
% entries, so this is the spectral condition without an eigensolve (the C++ twin
% has no LAPACK to call).
if any(~isfinite(tailMass)) || any(tailMass < 1 - 1e-9)
    return;
end

B00 = Qk(ag_blk(0, mph), ag_blk(0, mph));
Sys = [B00, A0; A2, A1 + R * A2];
% Replace one column by the normalization pi_0 e + pi_1 (I - R)^-1 e = 1.
Sys(:, 1) = [ones(mph, 1); tailMass];
rhs = zeros(1, 2 * mph);
rhs(1) = 1;
v = rhs / Sys;
if any(~isfinite(v))
    return;
end

g = struct();
g.R = R;
g.pi0 = v(1:mph);
g.pi1 = v(mph+1:end);
% sum_{n>=1} pi_n, the row vector every tail moment is read off.
g.busy = g.pi1 / IR;
% E[N] = sum_{n>=1} n pi_n e = pi_1 (I - R)^-2 e.
g.qlen = (g.pi1 / IR) / IR * ones(mph, 1);
end

function pi = qbd_tail_expand(g, nlev, mph)
% Materialize the matrix-geometric tail over NLEV levels, so the block norm of
% the fixed point and the RCAT residual read one vector shape for every
% component. The metrics use the closed forms in G instead.
pi = zeros(1, nlev * mph);
pi(ag_blk(0, mph)) = g.pi0;
v = g.pi1;
for n = 1:(nlev - 1)
    pi(ag_blk(n, mph)) = v;
    v = v * g.R;
end
end

function rho = qbd_scalar_rho(f, b, g)
% Sub-unit root rho of the scalar QBD characteristic equation
%   b*rho^2 - (f+b+g)*rho + f = 0,
% where f is the up-1 rate, b the down-1 rate, and g the extra local outflow
% (catastrophe drain to the empty state). This is the block-size-1 instance
% of Neuts' rate matrix R. For b == 0 the equation degenerates to the
% catastrophe-stabilised ratio rho = f/(f+g).

if b <= 1e-14
    if f + g <= 0
        rho = Inf;
    else
        rho = f / (f + g);
    end
    return;
end

c1 = -(f + b + g);
disc = c1^2 - 4 * b * f;
if disc < 0
    rho = Inf;
    return;
end
sq = sqrt(disc);
r1 = (-c1 - sq) / (2 * b);
r2 = (-c1 + sq) / (2 * b);
cands = sort([r1, r2]);
if cands(1) > 0
    rho = cands(1);
else
    rho = cands(2);
end

end

function [R, AP, processMap, actionMap, N, meta] = build_rcat(sn, maxStates)
% BUILD_RCAT Convert LINE network structure to RCAT format
%
% META carries the QBD shape of every component: NLEV levels of MPH phases
% each, LEVEL the level index of every state, and SVCRATE the service
% completion rate out of every state (zero on level 0, where no server runs).

if nargin < 2
    maxStates = 100;
end

M = sn.nstations;
K = sn.nclasses;
rt = sn.rt;  % (M*K) x (M*K) routing table

% Identify station types
sourceStations = [];
sinkStations = [];
queueStations = [];

for ist = 1:M
    nodeIdx = sn.stationToNode(ist);
    if sn.nodetype(nodeIdx) == NodeType.Source
        sourceStations(end+1) = ist;
    elseif sn.nodetype(nodeIdx) == NodeType.Sink
        sinkStations(end+1) = ist;
    else
        % Queue, Delay, or other service stations
        queueStations(end+1) = ist;
    end
end

% Create process mapping: each (station, class) pair at queue stations
% Note: Signal classes (negative customers) don't create separate processes
% as they only modify the state of positive customer processes
processIdx = 0;
processMap = zeros(M, K);
for ist = queueStations
    for r = 1:K
        % Skip Signal classes - they don't have their own queue state
        if sn.issignal(r)
            continue;
        end
        % Check if this station serves this class
        if ~isnan(sn.rates(ist, r)) && sn.rates(ist, r) > 0
            processIdx = processIdx + 1;
            processMap(ist, r) = processIdx;
        end
    end
end
numProcesses = processIdx;

if numProcesses == 0
    R = {};
    AP = [];
    actionMap = [];
    N = [];
    meta = struct('nlev', [], 'mph', [], 'level', {{}}, 'svcrate', {{}}, 'svcdown', {{}});
    return;
end

% G-network signals modify positive-customer processes, not their own;
% see _kb/06-solver-catalog.md for rationale

% Identify sink nodes (nodetype = -1 = NodeType.Sink)
% Use row vector to ensure for-loop doesn't execute when empty
sinkNodes = find(sn.nodetype == NodeType.Sink)';
if isempty(sinkNodes)
    sinkNodes = [];  % Ensure empty row vector, not column
end

% Markovian shape of each process: the levels it spans, the arrival MAP of the
% external streams reaching it and the service MAP of its station.
N = zeros(1, numProcesses);
meta = struct();
meta.nlev = zeros(1, numProcesses);
meta.mph = zeros(1, numProcesses);
meta.level = cell(1, numProcesses);
meta.svcrate = cell(1, numProcesses);
meta.svcdown = cell(1, numProcesses);
pinfo = cell(1, numProcesses);
for p = 1:numProcesses
    [ist, r] = find(processMap == p);
    ist = ist(1); r = r(1);
    pi_p = struct();
    pi_p.ist = ist;
    pi_p.r = r;
    [pi_p.Ds0, pi_p.Ds1] = rcat_proc_map(sn, ist, r);
    [pi_p.Da0, pi_p.Da1, pi_p.lamNeg, pi_p.lamCat, pi_p.batch] = ...
        rcat_arrival_map(sn, ist, r, rt, sourceStations, K);
    pi_p.ns = size(pi_p.Ds0, 1);
    pi_p.na = size(pi_p.Da0, 1);
    pi_p.mph = pi_p.na * pi_p.ns;
    if sn.njobs(r) < Inf  % Closed class
        pi_p.nlev = sn.njobs(r) + 1;  % Levels 0, 1, ..., njobs
    else  % Open class
        pi_p.nlev = maxStates;  % Truncate at maxStates
    end
    % Service completion: level down, arrival phase untouched (qbd_mapmap1).
    pi_p.Dsvc = kron(eye(pi_p.na), pi_p.Ds1);
    pi_p.N = pi_p.nlev * pi_p.mph;
    pinfo{p} = pi_p;

    N(p) = pi_p.N;
    meta.nlev(p) = pi_p.nlev;
    meta.mph(p) = pi_p.mph;
    meta.level{p} = reshape(repmat(0:(pi_p.nlev-1), pi_p.mph, 1), 1, pi_p.N);
    svcrow = sum(pi_p.Dsvc, 2)';
    meta.svcrate{p} = [zeros(1, pi_p.mph), repmat(svcrow, 1, pi_p.nlev - 1)];
    meta.svcdown{p} = svcrow;
end

% Count actions: each routing transition (i,r) -> (j,s) where P > 0
actionIdx = 0;
actionMap = struct('from_station', {}, 'from_class', {}, ...
                   'to_station', {}, 'to_class', {}, 'prob', {}, ...
                   'isNegative', {}, 'isCatastrophe', {}, 'removalDistribution', {});

for ist = queueStations
    for r = 1:K
        if processMap(ist, r) > 0
            % Check if class r is a removal signal class (NEGATIVE or
            % CATASTROPHE; the two are distinct SignalType values, so both
            % must be tested).
            isNegativeClass = false;
            isCatastropheClass = false;
            removalDist = [];
            if sn.issignal(r) && ~isnan(sn.signaltype{r}) && ...
                    (sn.signaltype{r} == SignalType.NEGATIVE || sn.signaltype{r} == SignalType.CATASTROPHE)
                isNegativeClass = true;
                % Check if class r is a catastrophe signal
                if (isfield(sn, 'iscatastrophe') && ~isempty(sn.iscatastrophe) && sn.iscatastrophe(r) > 0) ...
                        || sn.signaltype{r} == SignalType.CATASTROPHE
                    isCatastropheClass = true;
                end
                % Get removal distribution for this class
                if isfield(sn, 'signalremdist') && ~isempty(sn.signalremdist) && r <= length(sn.signalremdist)
                    removalDist = sn.signalremdist{r};
                end
            end

            for jst = queueStations
                for s = 1:K
                    if processMap(jst, s) > 0
                        % Get routing probability
                        prob_ij_rs = rt((ist-1)*K + r, (jst-1)*K + s);
                        if prob_ij_rs > 0 && (ist ~= jst || r ~= s)
                            % This is an action (departure from i,r triggers arrival at j,s)
                            actionIdx = actionIdx + 1;
                            actionMap(actionIdx).from_station = ist;
                            actionMap(actionIdx).from_class = r;
                            actionMap(actionIdx).to_station = jst;
                            actionMap(actionIdx).to_class = s;
                            actionMap(actionIdx).prob = prob_ij_rs;
                            actionMap(actionIdx).isNegative = isNegativeClass;
                            actionMap(actionIdx).isCatastrophe = isCatastropheClass;
                            actionMap(actionIdx).removalDistribution = removalDist;
                        end
                    end
                end
            end
        end
    end
end
numActions = actionIdx;

% Initialize R and AP
R = cell(numActions + 1, max(numProcesses, 2));
if numActions > 0
    AP = zeros(numActions, 2);
else
    AP = zeros(0, 2);  % Empty matrix when no actions
end

% Build local/hidden rate matrices L for each process (R{end,k})
for p = 1:numProcesses
    R{numActions + 1, p} = build_local_rates(sn, pinfo{p}, rt, sinkNodes, K);
end

% Build active and passive matrices for each action
for a = 1:numActions
    am = actionMap(a);

    % Active process (departure)
    ist = am.from_station;
    r = am.from_class;
    p_active = processMap(ist, r);
    AP(a, 1) = p_active;

    pa = pinfo{p_active};
    prob = am.prob;

    % Active matrix: level n -> n-1 carrying the service completion block
    % kron(I, D1^s), scaled by the routing probability of this action.
    Aa = zeros(pa.N);
    for n = 1:(pa.nlev - 1)
        Aa(ag_blk(n, pa.mph), ag_blk(n - 1, pa.mph)) = pa.Dsvc * prob;
    end
    % Boundary self-loop physical only for closed class (open-truncation bias);
    % see _kb/06-solver-catalog.md for rationale. It is written on the diagonal
    % so it stays inert in the generator while still contributing the
    % pi(i)/pi(i) = 1 ratio the INAP estimators read off the top level.
    if sn.njobs(r) < Inf
        top = ag_blk(pa.nlev - 1, pa.mph);
        Aa(top, top) = Aa(top, top) + diag(sum(pa.Dsvc, 2) * prob);
    end
    R{a, 1} = Aa;

    % Passive process (arrival or signal effect)
    jst = am.to_station;
    s = am.to_class;
    p_passive = processMap(jst, s);
    AP(a, 2) = p_passive;

    pp = pinfo{p_passive};
    Im = eye(pp.mph);
    Pb = zeros(pp.N);
    if am.isNegative
        % NEGATIVE: Job removal at destination (G-network negative customer)
        if am.isCatastrophe
            % CATASTROPHE: All jobs are removed - every level drops to level 0
            for n = 0:(pp.nlev - 1)
                Pb(ag_blk(n, pp.mph), ag_blk(0, pp.mph)) = Im;
            end
        elseif ~isempty(am.removalDistribution)
            % BATCH REMOVAL: Remove a random number of jobs based on distribution
            % P[n, m] = probability of transition from n to m jobs
            Pb = add_batch_removal(Pb, am.removalDistribution, 1, pp.nlev, pp.mph);
            Pb(ag_blk(0, pp.mph), ag_blk(0, pp.mph)) = Im;  % empty queue absorbs the signal
        else
            % DEFAULT: Remove exactly 1 job (original behavior)
            % Empty queue: no effect
            Pb(ag_blk(0, pp.mph), ag_blk(0, pp.mph)) = Im;
            % Non-empty queues: decrement
            for n = 1:(pp.nlev - 1)
                Pb(ag_blk(n, pp.mph), ag_blk(n - 1, pp.mph)) = Im;
            end
        end
    else
        % POSITIVE: Normal job arrival at destination. The phase is untouched:
        % a job joining does not restart the server, and the service phase
        % frozen at level 0 is the one the last completion left behind, which
        % for a phase-type is already its entry distribution.
        for n = 0:(pp.nlev - 2)
            Pb(ag_blk(n, pp.mph), ag_blk(n + 1, pp.mph)) = Im;
        end
        % Boundary: at max capacity
        top = ag_blk(pp.nlev - 1, pp.mph);
        Pb(top, top) = Im;
    end
    R{a, 2} = Pb;
end

end

function [D0, D1] = rcat_proc_map(sn, ist, r)
% (D0,D1) of the process at (IST,R), or the exponential pair built from
% sn.rates when the struct carries no usable matrix representation.
%
% A non-Markovian pair (a RAP or an ME) would assemble a rational generator
% rather than a CTMC, so it is refused here and answered as its mean rate; the
% solver-level gate rejects those models before they reach this point.
D0 = [];
D1 = [];
if isfield(sn, 'proc') && ~isempty(sn.proc) && ist <= numel(sn.proc) && ...
        ~isempty(sn.proc{ist}) && iscell(sn.proc{ist}) && r <= numel(sn.proc{ist})
    pr = sn.proc{ist}{r};
    if iscell(pr) && numel(pr) >= 2
        c0 = full(pr{1});
        c1 = full(pr{2});
        if mam_is_markovian_map(c0, c1)
            D0 = c0;
            D1 = c1;
        end
    end
end
if isempty(D0)
    rate = sn.rates(ist, r);
    if isnan(rate) || rate <= 0
        rate = 0;
    end
    D0 = -rate;
    D1 = rate;
end
end

function [Da0, Da1, lamNeg, lamCat, batch] = rcat_arrival_map(sn, ist, r, rt, sourceStations, K)
% External (Source) streams reaching (IST,R), as one arrival MAP for the
% positive customers plus the scalar rates of the removal signals.
%
% Each stream is thinned by its routing probability -- a MAP thinned with
% probability p is (D0 + (1-p) D1, p D1) -- and the streams are superposed by
% the Kronecker sum, so several Poisson sources still collapse to the single
% rate sum this analyzer used before. Removal signals stay scalar: a signal is
% a trigger with no service, and its arrival process is required exponential.

Da0 = 0;
Da1 = 0;
haveArrival = false;
lamNeg = 0;   % negative arrivals with single removal (default)
lamCat = 0;   % catastrophe arrivals (remove all)
batch = {};   % batch removal arrivals: {rate, distribution} pairs

for isrc = sourceStations
    for s_src = 1:K
        % Check if source class s_src is a signal
        isSignal = sn.issignal(s_src);

        if isSignal
            % For signals: they route to themselves (Signal -> Signal), but their effect
            % is on positive customers at the destination station. We check if the signal
            % routes to ANY class at this station (not just class r).
            prob_src = 0;
            for s_dst = 1:K
                prob_src = prob_src + rt((isrc-1)*K + s_src, (ist-1)*K + s_dst);
            end
        else
            % For regular classes: direct routing to (ist, r)
            prob_src = rt((isrc-1)*K + s_src, (ist-1)*K + r);
        end

        if prob_src > 0 && ~isnan(sn.rates(isrc, s_src))
            srcRate = sn.rates(isrc, s_src);
            % Check if source class s_src is a negative or catastrophe signal
            if isSignal && ~isnan(sn.signaltype{s_src}) && ...
                    (sn.signaltype{s_src} == SignalType.NEGATIVE || sn.signaltype{s_src} == SignalType.CATASTROPHE)
                % Check if it's a catastrophe (either via iscatastrophe flag or signaltype)
                isCat = (isfield(sn, 'iscatastrophe') && ~isempty(sn.iscatastrophe) && sn.iscatastrophe(s_src)) || ...
                        sn.signaltype{s_src} == SignalType.CATASTROPHE;
                if isCat
                    lamCat = lamCat + srcRate * prob_src;
                else
                    % Check if it has a removal distribution
                    removalDist = [];
                    if isfield(sn, 'signalremdist') && ~isempty(sn.signalremdist) && s_src <= length(sn.signalremdist)
                        removalDist = sn.signalremdist{s_src};
                    end
                    if ~isempty(removalDist)
                        batch{end+1} = {srcRate * prob_src, removalDist};
                    else
                        lamNeg = lamNeg + srcRate * prob_src;
                    end
                end
            elseif srcRate > 0
                [S0, S1] = rcat_proc_map(sn, isrc, s_src);
                if prob_src < 1
                    S0 = S0 + (1 - prob_src) * S1;
                    S1 = prob_src * S1;
                end
                if haveArrival
                    Da0 = krons(Da0, S0);
                    Da1 = krons(Da1, S1);
                else
                    Da0 = S0;
                    Da1 = S1;
                    haveArrival = true;
                end
            end
        end
    end
end
end

function B = add_batch_removal(B, dist, rate, nlev, mph)
% Accumulate the level n -> level m block of a batch removal into B, scaled by
% RATE. Landing on the empty level absorbs the whole upper tail of the pmf,
% which is what keeps the block stochastic once the batch exceeds the queue
% length. The phase is untouched: a removal takes a waiting job, not the one in
% service.
Im = eye(mph);
for n = 1:(nlev - 1)
    for m = 1:n
        p = dist.evalPMF(n - m);
        if p > 0
            B(ag_blk(n, mph), ag_blk(m, mph)) = B(ag_blk(n, mph), ag_blk(m, mph)) + rate * p * Im;
        end
    end
    cdf = 0;
    for j = 0:(n-1)
        cdf = cdf + dist.evalPMF(j);
    end
    tail = 1 - cdf;
    if tail > 0
        B(ag_blk(n, mph), ag_blk(0, mph)) = B(ag_blk(n, mph), ag_blk(0, mph)) + rate * tail * Im;
    end
end
end

function L = build_local_rates(sn, p, rt, sinkNodes, K)
% Build local/hidden transition matrix for the component P
% Note: sinkNodes contains node indices (not station indices) for Sink nodes

ist = p.ist;
r = p.r;
mph = p.mph;
nlev = p.nlev;
Im = eye(mph);
L = zeros(p.N);

% Level-local blocks: the arrival phase always runs, the service phase only
% while the server is busy (qbd_mapmap1's Lbar = kron(D0^a, I) at level 0 and
% L = krons(D0^a, D0^s) above it). With one phase each these are pure
% diagonals, which ctmc_makeinfgen discards and rebuilds from the row sums.
L(ag_blk(0, mph), ag_blk(0, mph)) = kron(p.Da0, eye(p.ns));
Lbusy = krons(p.Da0, p.Ds0);
for n = 1:(nlev - 1)
    L(ag_blk(n, mph), ag_blk(n, mph)) = Lbusy;
end

% Positive arrival transitions: level n -> n+1, carrying kron(D1^a, I).
Aup = kron(p.Da1, eye(p.ns));
for n = 0:(nlev - 2)
    L(ag_blk(n, mph), ag_blk(n + 1, mph)) = L(ag_blk(n, mph), ag_blk(n + 1, mph)) + Aup;
end
% At the truncation the job is lost but the arrival process still moves on, so
% the block stays on the top level. With a single arrival phase this is a pure
% diagonal and is discarded, exactly as before.
top = ag_blk(nlev - 1, mph);
L(top, top) = L(top, top) + Aup;

% Catastrophe arrival transitions: every busy level drops to level 0
if p.lamCat > 0
    for n = 1:(nlev - 1)
        L(ag_blk(n, mph), ag_blk(0, mph)) = L(ag_blk(n, mph), ag_blk(0, mph)) + p.lamCat * Im;
    end
end

% Batch removal arrival transitions: level n -> m at rate lambda*P(remove n-m)
for b = 1:length(p.batch)
    L = add_batch_removal(L, p.batch{b}{2}, p.batch{b}{1}, nlev, mph);
end

% Single removal negative arrival transitions: level n -> n-1 (busy levels only)
if p.lamNeg > 0
    for n = 1:(nlev - 1)
        L(ag_blk(n, mph), ag_blk(n - 1, mph)) = L(ag_blk(n, mph), ag_blk(n - 1, mph)) + p.lamNeg * Im;
    end
end

% Service completions that are not synchronizing actions: departures to a Sink
% (level down) and self-routing (level unchanged, service restarted).
rate = sn.rates(ist, r);
if ~isnan(rate) && rate > 0
    % Departures to sink (use rtnodes with node indices)
    nodeIdx = sn.stationToNode(ist);

    prob_sink = 0;
    if isfield(sn, 'rtnodes') && ~isempty(sn.rtnodes)
        for jsnk = sinkNodes
            for s = 1:K
                % rtnodes indices: (nodeIdx-1)*K + classIdx
                fromIdx = (nodeIdx - 1) * K + r;
                toIdx = (jsnk - 1) * K + s;
                if fromIdx <= size(sn.rtnodes, 1) && toIdx <= size(sn.rtnodes, 2)
                    prob_sink = prob_sink + sn.rtnodes(fromIdx, toIdx);
                end
            end
        end
    end

    % Self-routing (stays at same station, same class)
    prob_self = rt((ist-1)*K + r, (ist-1)*K + r);

    if prob_sink > 0
        for n = 1:(nlev - 1)
            L(ag_blk(n, mph), ag_blk(n - 1, mph)) = ...
                L(ag_blk(n, mph), ag_blk(n - 1, mph)) + p.Dsvc * prob_sink;
        end
    end
    if prob_self > 0
        for n = 1:(nlev - 1)
            L(ag_blk(n, mph), ag_blk(n, mph)) = L(ag_blk(n, mph), ag_blk(n, mph)) + p.Dsvc * prob_self;
        end
    end
end

end

function [QN, UN, RN, TN, CN, XN] = rcat_metrics(sn, x, pi, Q, processMap, actionMap, N, meta, ...
    rhoProc, isGeomProc, geomData)
% RCAT_METRICS Convert RCAT solution to LINE performance metrics
%
% When RHOPROC/ISGEOMPROC/GEOMDATA are supplied (matrix-geometric 'inapinf'
% method), processes flagged geometric use the exact closed-form moments of the
% infinite marginal instead of the truncated explicit vector pi.

if nargin < 9, rhoProc = []; end
if nargin < 10, isGeomProc = []; end
if nargin < 11, geomData = {}; end

M = sn.nstations;
K = sn.nclasses;

QN = zeros(M, K);
UN = zeros(M, K);
RN = zeros(M, K);
TN = zeros(M, K);

% Compute metrics for each (station, class) pair
for ist = 1:M
    for r = 1:K
        p = processMap(ist, r);
        if p > 0 && ~isempty(pi) && p <= length(pi) && ~isempty(pi{p})
            mph = meta.mph(p);

            if ~isempty(isGeomProc) && p <= numel(isGeomProc) && isGeomProc(p)
                if mph == 1
                    % Infinite geometric marginal pi_n = (1-rho) rho^n:
                    %   E[N] = rho/(1-rho),  P(N>0) = rho.
                    rho = rhoProc(p);
                    QN(ist, r) = rho / (1 - rho);
                    UN(ist, r) = rho;
                    mu_ir = sn.rates(ist, r);
                    if ~isnan(mu_ir) && mu_ir > 0
                        TN(ist, r) = mu_ir * rho;
                    end
                else
                    % Matrix-geometric tail pi_(n+1) = pi_n R.
                    g = geomData{p};
                    QN(ist, r) = g.qlen;
                    UN(ist, r) = sum(g.busy);
                    TN(ist, r) = g.busy * meta.svcdown{p}(:);
                end
            else
                v = pi{p}(:)';

                % Queue length: E[N] = sum over states of level(state)*pi(state)
                QN(ist, r) = meta.level{p} * v(:);

                % Utilization: P(N > 0) = 1 - P(level 0)
                UN(ist, r) = 1 - sum(v(1:mph));

                % Throughput: the rate of service completions, i.e. the
                % phase-dependent departure rate averaged over the marginal.
                % With one phase this is the mean rate times P(N>0).
                TN(ist, r) = meta.svcrate{p} * v(:);
            end
        end
    end
end

% Handle self-looping classes: they always stay at their reference station
% and share the server with other classes under PS scheduling.
% Only override if the method did not already compute SLC metrics (QN == 0).
if isfield(sn, 'isslc') && any(sn.isslc)
    for r = 1:K
        if sn.isslc(r)
            refst = sn.refstat(r);
            if refst > 0 && refst <= M && QN(refst, r) == 0
                % Self-looping class: all jobs stay at reference station
                QN(refst, r) = sn.njobs(r);

                % Service rate for this class
                mu_ir = sn.rates(refst, r);
                if ~isnan(mu_ir) && mu_ir > 0
                    nservers = sn.nservers(refst);
                    if isinf(nservers)
                        % Delay (infinite server): no capacity constraint,
                        % each job gets dedicated service
                        UN(refst, r) = QN(refst, r);
                        TN(refst, r) = mu_ir * QN(refst, r);
                    else
                        % Queue (finite server): capacity constraint applies
                        % Get utilization from other classes at this station
                        other_util = 0;
                        for s = 1:K
                            if s ~= r && ~sn.isslc(s)
                                other_util = other_util + UN(refst, s);
                            end
                        end

                        % Remaining capacity is shared with SLC
                        remaining_capacity = max(0, 1 - other_util);

                        % SLC utilization: min(demand, remaining capacity)
                        slc_demand = QN(refst, r) / mu_ir;
                        UN(refst, r) = min(slc_demand, remaining_capacity);
                        TN(refst, r) = mu_ir * UN(refst, r);
                    end
                end
            end
        end
    end
end

% Response times from Little's law: R = Q / T
for ist = 1:M
    for r = 1:K
        if TN(ist, r) > 0
            RN(ist, r) = QN(ist, r) / TN(ist, r);
        else
            RN(ist, r) = 0;
        end
    end
end

% System metrics
CN = zeros(1, K);
XN = zeros(1, K);

% For open classes: system throughput = arrival rate, system response time = sum of response times
for r = 1:K
    if sn.njobs(r) >= Inf  % Open class
        % System throughput equals arrival rate (from source)
        for ist = 1:M
            nodeIdx = sn.stationToNode(ist);
            if sn.nodetype(nodeIdx) == NodeType.Source
                XN(r) = sn.rates(ist, r);
                break;
            end
        end
        % System response time = sum over all stations
        CN(r) = sum(RN(:, r));
    else
        % Closed class: use reference station
        refst = sn.refstat(r);
        if refst > 0 && refst <= M
            XN(r) = TN(refst, r);
            if XN(r) > 0
                CN(r) = sn.njobs(r) / XN(r);
            end
        end
    end
end

end
