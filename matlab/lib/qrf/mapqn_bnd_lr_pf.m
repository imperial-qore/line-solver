function [result, x, fval, exitflag] = mapqn_bnd_lr_pf(params, objective_queue, sense)
% MAPQN_BND_LR_PF - Product-Form Linear Reduction Bounds
%
% MATLAB port of bnd_lr_pf.py
%
% Computes lower/upper bounds on utilization for closed queueing networks
% using a linear programming formulation that exploits product-form
% structure (no phases, scalar service rates per queue).
%
% Usage:
%   [result, x, fval, exitflag] = mapqn_bnd_lr_pf(params)
%   [result, x, fval, exitflag] = mapqn_bnd_lr_pf(params, objective_queue)
%   [result, x, fval, exitflag] = mapqn_bnd_lr_pf(params, objective_queue, sense)
%
% Inputs:
%   params          - Structure with model parameters:
%                     .M       - Number of queues
%                     .N       - Total population
%                     .mu      - [M x 1] Service rates (scalar per queue)
%                     .r       - [M x M] Routing probability matrix
%                     .verbose - (optional) boolean, default true
%
%   objective_queue - (optional) Queue index for objective (1-based), default 1
%   sense           - (optional) 'min' (default) or 'max'
%
% Outputs:
%   result   - Structure with results:
%              .objective - Objective function value
%              .exitflag  - Solver exit flag
%              .U         - [M x 1] Utilizations
%              .Q         - [M x 1] Mean queue lengths
%   x        - Raw solution vector
%   fval     - Objective function value
%   exitflag - Solver exit flag

    if nargin < 2 || isempty(objective_queue)
        objective_queue = 1;
    end
    if nargin < 3 || isempty(sense)
        sense = 'min';
    end

    % Extract parameters
    M = params.M;
    N = params.N;
    mu = params.mu(:);
    r = params.r;
    if isfield(params, 'verbose')
        verbose = params.verbose;
    else
        verbose = true;
    end

    % Compute transition rates q(i,j) = r(i,j) * mu(i)
    q = zeros(M, M);
    for i = 1:M
        for j = 1:M
            q(i,j) = r(i,j) * mu(i);
        end
    end

    %% Build variable indexing
    % Variables (all 1-based):
    %   U(i)          for i in 1:M
    %   Q(i)          for i in 1:M
    %   C(j,i)        for j in 1:M, i in 1:M
    %   p1(j,i,ni)    for j in 1:M, i in 1:M, ni in 0:N
    %   p1c(j,i,ni)   for j in 1:M, i in 1:M, ni in 0:N

    if verbose; fprintf('Building variable index map...\n'); end
    varCount = 0;

    % U variables: utilization
    Uidx = zeros(M, 1);
    for i = 1:M
        varCount = varCount + 1;
        Uidx(i) = varCount;
    end

    % Q variables: mean queue length
    Qidx = zeros(M, 1);
    for i = 1:M
        varCount = varCount + 1;
        Qidx(i) = varCount;
    end

    % C variables: conditional queue lengths C(j,i)
    Cidx = zeros(M, M);
    for j = 1:M
        for i = 1:M
            varCount = varCount + 1;
            Cidx(j,i) = varCount;
        end
    end

    % p1 variables: marginal probabilities p1(j,i,ni), ni in 0:N
    p1idx = zeros(M, M, N+1);
    for j = 1:M
        for i = 1:M
            for ni = 0:N
                varCount = varCount + 1;
                p1idx(j, i, ni+1) = varCount;
            end
        end
    end

    % p1c variables: complementary marginal p1c(j,i,ni), ni in 0:N
    p1cidx = zeros(M, M, N+1);
    for j = 1:M
        for i = 1:M
            for ni = 0:N
                varCount = varCount + 1;
                p1cidx(j, i, ni+1) = varCount;
            end
        end
    end

    nVars = varCount;
    if verbose; fprintf('Total variables: %d\n', nVars); end

    %% Initialize bounds
    lb = zeros(nVars, 1);
    ub = inf(nVars, 1);

    % U bounds: [0, 1]
    for i = 1:M
        ub(Uidx(i)) = 1;
    end

    % Q bounds: [0, N]
    for i = 1:M
        ub(Qidx(i)) = N;
    end

    % C bounds: [0, N]
    for j = 1:M
        for i = 1:M
            ub(Cidx(j,i)) = N;
        end
    end

    % p1 bounds: [0, 1]
    for j = 1:M
        for i = 1:M
            for ni = 0:N
                ub(p1idx(j, i, ni+1)) = 1;
            end
        end
    end

    % p1c bounds: [0, 1]
    for j = 1:M
        for i = 1:M
            for ni = 0:N
                ub(p1cidx(j, i, ni+1)) = 1;
            end
        end
    end

    %% Build constraints
    Aeq = [];
    beq = [];

    if verbose; fprintf('Building constraints...\n'); end

    %% ZER1: p1(j,j,0) = 0 for all j (via upper bound)
    if verbose; fprintf('  ZER1 constraints...\n'); end
    for j = 1:M
        ub(p1idx(j, j, 0+1)) = 0;
    end

    %% ZER3: p1(j,i,N) = 0 for j ~= i (via upper bound)
    if verbose; fprintf('  ZER3 constraints...\n'); end
    for j = 1:M
        for i = 1:M
            if j ~= i
                ub(p1idx(j, i, N+1)) = 0;
            end
        end
    end

    %% CEQU: C(j,j) = Q(j) for all j
    if verbose; fprintf('  CEQU constraints...\n'); end
    for j = 1:M
        row = zeros(1, nVars);
        row(Cidx(j,j)) = 1;
        row(Qidx(j)) = -1;
        Aeq = [Aeq; row];
        beq = [beq; 0];
    end

    %% ONE1: sum over ni of (p1(j,i,ni) + p1c(j,i,ni)) = 1 for each (j,i)
    if verbose; fprintf('  ONE1 constraints...\n'); end
    for j = 1:M
        for i = 1:M
            row = zeros(1, nVars);
            for ni = 0:N
                row(p1idx(j, i, ni+1)) = 1;
                row(p1cidx(j, i, ni+1)) = 1;
            end
            Aeq = [Aeq; row];
            beq = [beq; 1];
        end
    end

    %% UTIL: U(i) = sum over t,nt of p1(i,t,nt) for each (i,t)
    if verbose; fprintf('  UTIL constraints...\n'); end
    for i = 1:M
        for t = 1:M
            row = zeros(1, nVars);
            row(Uidx(i)) = 1;
            for nt = 0:N
                row(p1idx(i, t, nt+1)) = -1;
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% QLEN: Q(i) = sum over ni of ni*p1(i,i,ni) for each i
    if verbose; fprintf('  QLEN constraints...\n'); end
    for i = 1:M
        row = zeros(1, nVars);
        row(Qidx(i)) = 1;
        for ni = 0:N
            row(p1idx(i, i, ni+1)) = -ni;
        end
        Aeq = [Aeq; row];
        beq = [beq; 0];
    end

    %% CLEN: C(j,i) = sum over ni of ni*p1(j,i,ni) for each (j,i)
    if verbose; fprintf('  CLEN constraints...\n'); end
    for j = 1:M
        for i = 1:M
            row = zeros(1, nVars);
            row(Cidx(j,i)) = 1;
            for ni = 0:N
                row(p1idx(j, i, ni+1)) = -ni;
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% MPCB: sum over i of C(j,i) = N*U(j) for each j
    if verbose; fprintf('  MPCB constraints...\n'); end
    for j = 1:M
        row = zeros(1, nVars);
        for i = 1:M
            row(Cidx(j,i)) = 1;
        end
        row(Uidx(j)) = -N;
        Aeq = [Aeq; row];
        beq = [beq; 0];
    end

    %% POPC: sum over i of Q(i) = N
    if verbose; fprintf('  POPC constraint...\n'); end
    row = zeros(1, nVars);
    for i = 1:M
        row(Qidx(i)) = 1;
    end
    Aeq = [Aeq; row];
    beq = [beq; N];

    %% GFFL0: Global flow for ni=0
    % For each i: sum over j~=i of q(j,i)*p1(j,i,0) = sum over j~=i of q(i,j)*p1(i,i,1)
    if verbose; fprintf('  GFFL0 constraints...\n'); end
    for i = 1:M
        row = zeros(1, nVars);
        for j = 1:M
            if j ~= i
                row(p1idx(j, i, 0+1)) = row(p1idx(j, i, 0+1)) + q(j,i);
                row(p1idx(i, i, 1+1)) = row(p1idx(i, i, 1+1)) - q(i,j);
            end
        end
        Aeq = [Aeq; row];
        beq = [beq; 0];
    end

    %% GFFL: Global flow for ni in 1..N-1
    % For each (i,ni): sum over j~=i of q(j,i)*p1(j,i,ni) = sum over j~=i of q(i,j)*p1(i,i,ni+1)
    if verbose; fprintf('  GFFL constraints...\n'); end
    for i = 1:M
        for ni = 1:(N-1)
            row = zeros(1, nVars);
            for j = 1:M
                if j ~= i
                    row(p1idx(j, i, ni+1)) = row(p1idx(j, i, ni+1)) + q(j,i);
                    row(p1idx(i, i, ni+1+1)) = row(p1idx(i, i, ni+1+1)) - q(i,j);
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% UJNT: Joint probability symmetry
    % sum over ni=1..N of p1(j,i,ni) = sum over nj=1..N of p1(i,j,nj) for each (i,j)
    if verbose; fprintf('  UJNT constraints...\n'); end
    for i = 1:M
        for j = 1:M
            row = zeros(1, nVars);
            for ni = 1:N
                row(p1idx(j, i, ni+1)) = row(p1idx(j, i, ni+1)) + 1;
            end
            for nj = 1:N
                row(p1idx(i, j, nj+1)) = row(p1idx(i, j, nj+1)) - 1;
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% QBAL: Queue balance
    % For each i: sum over j~=i of (q(i,j)*U(i) - q(j,i)*sum_nj(p1(i,j,nj)) - q(j,i)*p1(j,i,0)) = 0
    if verbose; fprintf('  QBAL constraints...\n'); end
    for i = 1:M
        row = zeros(1, nVars);
        for j = 1:M
            if j ~= i
                row(Uidx(i)) = row(Uidx(i)) + q(i,j);
                for nj = 1:N
                    row(p1idx(i, j, nj+1)) = row(p1idx(i, j, nj+1)) - q(j,i);
                end
                row(p1idx(j, i, 0+1)) = row(p1idx(j, i, 0+1)) - q(j,i);
            end
        end
        Aeq = [Aeq; row];
        beq = [beq; 0];
    end

    %% Build objective function
    if verbose; fprintf('Building objective function...\n'); end
    c = zeros(nVars, 1);
    c(Uidx(objective_queue)) = 1;

    if strcmp(sense, 'max')
        c = -c;
    end

    %% Solve LP
    if verbose; fprintf('Solving LP with %d variables and %d equality constraints...\n', ...
        nVars, size(Aeq, 1)); end

    % LP algorithm. R2025a's default 'dual-simplex-highs' is broken in some
    % installs (errors "Unrecognized field name optimstatus"), so an
    % interior-point variant is required here as elsewhere in this family.
    %
    % THIS FILE DELIBERATELY DIFFERS FROM qrf_bas.m AND THE OTHER mapqn_bnd_*
    % FILES, WHICH DEFAULT TO 'interior-point-legacy'. DO NOT HARMONIZE THEM.
    % The sign of the tradeoff is OPPOSITE here. On the badly-scaled QRF/BAS
    % instances 'interior-point-legacy' is markedly the more accurate of the
    % two (see the comment in qrf_bas.m for the measured errors). On the LPs
    % this file builds it FAILS OUTRIGHT: for the M=2, N=2, mu=[1,1] tandem it
    % returns exitflag -2, "no feasible point found", on an LP that is
    % demonstrably feasible -- glpsol solves this function's own assembled Aeq,
    % beq, lb, ub to 2/3 in BOTH senses, matching the MVA reference and the
    % Python bnd_lr_pf, and the exact product-form point satisfies every
    % constraint block by hand. 'interior-point' returns 2/3 with exitflag 1.
    % It is also exact on this file's other test instances (the N=3 symmetric
    % tandem at 0.75, and the 3-queue asymmetric case, reproduced to ~1e-6), so
    % it is the correct primary here, not a fallback.
    %
    % Consequence worth knowing: a negative exitflag from this family is not
    % evidence that the model is infeasible. The [0.5, 1.0] that
    % 'interior-point-legacy' reported on the N=2 tandem was the residue of a
    % failed solve, not a loose bound. Cross-check against glpsol before
    % believing an exitflag -2.
    %
    % Override via params.lpAlgorithm.
    if isfield(params, 'lpAlgorithm') && ~isempty(params.lpAlgorithm)
        lpAlgorithm = params.lpAlgorithm;
    else
        lpAlgorithm = 'interior-point';
    end
    if verbose
        options = optimoptions('linprog', 'Display', 'final', 'Algorithm', lpAlgorithm);
    else
        options = optimoptions('linprog', 'Display', 'off', 'Algorithm', lpAlgorithm);
    end

    [x, fval, exitflag] = linprog(c, [], [], Aeq, beq, lb, ub, options);

    if strcmp(sense, 'max')
        fval = -fval;
    end

    %% Extract results
    result = struct();
    result.objective = fval;
    result.exitflag = exitflag;

    % Whether the metric fields below can be filled is a property of the
    % solution vector, not of exitflag. exitflag 0 (iteration limit) still
    % returns a usable interior point -- that is the normal outcome on badly
    % scaled instances such as the reference BAS network -- whereas a solve that
    % breaks down returns an empty or non-finite x, which must never be
    % consumed. Predicate on x accordingly, and say so rather than skipping
    % silently. (Observed here: exitflag -4 comes back with x empty.)
    hasSolution = ~isempty(x) && all(isfinite(x));
    if ~isempty(x) && ~hasSolution
        warning('mapqn_bnd_lr_pf:nonFiniteSolution', ...
            'linprog returned a non-finite solution (exitflag %d); U and the other metric fields are left unpopulated.', ...
            exitflag);
    end

    if hasSolution
        result.U = zeros(M, 1);
        result.Q = zeros(M, 1);

        for i = 1:M
            result.U(i) = x(Uidx(i));
            result.Q(i) = x(Qidx(i));
        end
    end

    if verbose; fprintf('\n=== Results ===\n'); end
    if verbose; fprintf('Objective value: %f\n', fval); end
    if verbose; fprintf('Exit flag: %d\n', exitflag); end
    if hasSolution && verbose
        fprintf('\nUtilizations:\n');
        for i = 1:M
            fprintf('  Queue %d: U = %.6f\n', i, result.U(i));
        end
        fprintf('\nMean queue lengths:\n');
        for i = 1:M
            fprintf('  Queue %d: Q = %.6f\n', i, result.Q(i));
        end
    end
end
