function result = mapqn_bnd_lr_mva(params, objective_queue, objective_level, sense, objective_var)
% MAPQN_BND_LR_MVA - MVA-based Linear Reduction Bounds for MAP Queueing Networks
%
% MATLAB port of bnd_lr_mva.py
%
% Usage:
%   result = mapqn_bnd_lr_mva(params)
%   result = mapqn_bnd_lr_mva(params, objective_queue)
%   result = mapqn_bnd_lr_mva(params, objective_queue, objective_level)
%   result = mapqn_bnd_lr_mva(params, objective_queue, objective_level, sense)
%
% Inputs:
%   params          - Structure with model parameters:
%                     .M      - Number of queues
%                     .N      - Total population
%                     .K      - Number of levels (scalar integer)
%                     .muM    - [M-1 x 1] service rates for non-MAP queues 1..M-1
%                     .muMAP  - [K x K] service rate matrix for MAP queue (queue M)
%                     .r      - [M x M] routing probability matrix
%                     .v      - [K x K] level change rate matrix
%                     .verbose - (optional) boolean
%
%   objective_queue - (optional) 1-based queue index, default 1
%   objective_level - (optional) 1-based level index, default 1. Pass 0 to
%                     optimize the AGGREGATE over levels, sum_k X(queue,k),
%                     which is the quantity the paper's bounds are stated on:
%                     U_i(N) = sum_k U_i^k(N) is the utilization of station i,
%                     while U_i^k alone is its utilization while the MAP sits
%                     in phase k. Maximizing the K terms separately and adding
%                     them is also an upper bound but a strictly looser one,
%                     since the phases cannot all peak at once.
%   sense           - (optional) 'min' or 'max', default 'max'
%   objective_var   - (optional) 'UN' (default) or 'QN', the variable family
%                     the objective is taken over
%
% Outputs:
%   result   - Structure with results:
%              .objective - Objective function value
%              .exitflag  - Solver exit flag
%              .UN        - [M x K] utilization matrix
%              .QN        - [M x K] queue length matrix
%
% Reference:
%   G. Casale, E. Smirni, "MAP-AMVA: Approximate Mean Value Analysis of Bursty
%   Systems", IEEE/IFIP DSN 2009, pp. 409-418. The LP assembled here is the
%   paper's MAP-AMVA optimization program: the population constraint (1), the
%   utilization bound (2), the MAP phase balance (3), the flow balance (4), the
%   generalized horizontal cut (12), the vertical-cut MVA relation (13) in the
%   linearized form (18)-(19) whose B(j,k,i) variables are the E_i^{j,k} of
%   Theorem 4, and the two auxiliary families QN <= N*UN and sum_w QN >= N*UN.

    if nargin < 2 || isempty(objective_queue)
        objective_queue = 1;
    end
    if nargin < 3 || isempty(objective_level)
        objective_level = 1;
    end
    if nargin < 4 || isempty(sense)
        sense = 'max';
    end
    if nargin < 5 || isempty(objective_var)
        objective_var = 'UN';
    end

    % Extract parameters
    M = params.M;
    N = params.N;
    K = params.K;
    muM = params.muM(:);
    muMAP = params.muMAP;
    r = params.r;
    v = params.v;
    if isfield(params, 'verbose')
        verbose = params.verbose;
    else
        verbose = true;
    end

    %% q function (all indices 1-based)
    % q(i,j,k,h) computes transition rate from queue i to queue j
    % at level k going to level h
    function val = q(i, j, k, h)
        if i < M
            % Non-MAP queue: scalar service rate
            if k == h
                val = r(i,j) * muM(i);
            else
                val = 0;
            end
        else
            % MAP queue (i == M)
            if j < M
                val = r(M,j) * muMAP(k,h);
            else
                % j == M
                if k ~= h
                    val = v(k,h) + r(M,M) * muMAP(k,h);
                else
                    val = 0;
                end
            end
        end
    end

    %% Build variable indexing
    % Variables are laid out as:
    %   UN(i,k)    for i=1..M, k=1..K   -> nVarsUN = M*K
    %   QN(i,k)    for i=1..M, k=1..K   -> nVarsQN = M*K
    %   B(j,k,i)   for j=1..M, k=1..K, i=1..M -> nVarsB = M*K*M

    nVarsUN = M * K;
    nVarsQN = M * K;
    nVarsB  = M * K * M;
    nVars   = nVarsUN + nVarsQN + nVarsB;

    % Index functions (1-based)
    getUNIdx = @(i, k) (i - 1) * K + k;
    getQNIdx = @(i, k) nVarsUN + (i - 1) * K + k;
    getBIdx  = @(j, k, i) nVarsUN + nVarsQN + (j - 1) * K * M + (k - 1) * M + i;

    if verbose; fprintf('Total variables: %d (UN=%d, QN=%d, B=%d)\n', nVars, nVarsUN, nVarsQN, nVarsB); end

    %% Initialize bounds
    lb = zeros(nVars, 1);
    ub = zeros(nVars, 1);

    % UN bounds: 0 <= UN(i,k) <= 1
    for i = 1:M
        for k = 1:K
            ub(getUNIdx(i, k)) = 1;
        end
    end

    % QN bounds: 0 <= QN(i,k) <= N
    for i = 1:M
        for k = 1:K
            ub(getQNIdx(i, k)) = N;
        end
    end

    % B bounds: 0 <= B(j,k,i) <= N
    for j = 1:M
        for k = 1:K
            for i = 1:M
                ub(getBIdx(j, k, i)) = N;
            end
        end
    end

    %% Build constraints using sparse triplets
    % We collect row, col, val triplets for equality and inequality matrices
    eq_rows = [];
    eq_cols = [];
    eq_vals = [];
    beq = [];
    nEq = 0;

    ineq_rows = [];
    ineq_cols = [];
    ineq_vals = [];
    bineq = [];
    nIneq = 0;

    if verbose; fprintf('Building constraints...\n'); end

    %% 1. QNB: QN(i,k) - B(j,k,i) >= 0 for all (i,k,j)
    %    Written as: -QN(i,k) + B(j,k,i) <= 0
    if verbose; fprintf('  QNB constraints...\n'); end
    for i = 1:M
        for k = 1:K
            for j = 1:M
                nIneq = nIneq + 1;
                ineq_rows = [ineq_rows; nIneq; nIneq];
                ineq_cols = [ineq_cols; getQNIdx(i, k); getBIdx(j, k, i)];
                ineq_vals = [ineq_vals; -1; 1];
                bineq = [bineq; 0];
            end
        end
    end

    %% 2. UMAX: sum_k UN(i,k) <= 1 for each i
    if verbose; fprintf('  UMAX constraints...\n'); end
    for i = 1:M
        nIneq = nIneq + 1;
        for k = 1:K
            ineq_rows = [ineq_rows; nIneq];
            ineq_cols = [ineq_cols; getUNIdx(i, k)];
            ineq_vals = [ineq_vals; 1];
        end
        bineq = [bineq; 1];
    end

    %% 3. POPCONSTR: sum over i,k of QN(i,k) = N
    if verbose; fprintf('  POPCONSTR constraint...\n'); end
    nEq = nEq + 1;
    for i = 1:M
        for k = 1:K
            eq_rows = [eq_rows; nEq];
            eq_cols = [eq_cols; getQNIdx(i, k)];
            eq_vals = [eq_vals; 1];
        end
    end
    beq = [beq; N];

    %% 4. FLOW: for each i: sum over k,m,w of (q(w,i,k,m)*UN(w,k) - q(i,w,m,k)*UN(i,m)) = 0
    if verbose; fprintf('  FLOW constraints...\n'); end
    for i = 1:M
        nEq = nEq + 1;
        row_coeffs = zeros(1, nVars);
        for k = 1:K
            for m = 1:K
                for w = 1:M
                    q_in = q(w, i, k, m);
                    q_out = q(i, w, m, k);
                    row_coeffs(getUNIdx(w, k)) = row_coeffs(getUNIdx(w, k)) + q_in;
                    row_coeffs(getUNIdx(i, m)) = row_coeffs(getUNIdx(i, m)) - q_out;
                end
            end
        end
        idx_nz = find(row_coeffs);
        eq_rows = [eq_rows; nEq * ones(length(idx_nz), 1)];
        eq_cols = [eq_cols; idx_nz(:)];
        eq_vals = [eq_vals; row_coeffs(idx_nz)'];
        beq = [beq; 0];
    end

    %% 5. UBAL: for each k: sum over h~=k,w of (q(M,w,k,h)*UN(M,k) - q(M,w,h,k)*UN(M,h)) = 0
    if verbose; fprintf('  UBAL constraints...\n'); end
    for k = 1:K
        nEq = nEq + 1;
        row_coeffs = zeros(1, nVars);
        for h = 1:K
            if h ~= k
                for w = 1:M
                    q_out = q(M, w, k, h);
                    q_in  = q(M, w, h, k);
                    row_coeffs(getUNIdx(M, k)) = row_coeffs(getUNIdx(M, k)) + q_out;
                    row_coeffs(getUNIdx(M, h)) = row_coeffs(getUNIdx(M, h)) - q_in;
                end
            end
        end
        idx_nz = find(row_coeffs);
        eq_rows = [eq_rows; nEq * ones(length(idx_nz), 1)];
        eq_cols = [eq_cols; idx_nz(:)];
        eq_vals = [eq_vals; row_coeffs(idx_nz)'];
        beq = [beq; 0];
    end

    %% 6. QBAL: for each k
    if verbose; fprintf('  QBAL constraints...\n'); end
    for k = 1:K
        nEq = nEq + 1;
        row_coeffs = zeros(1, nVars);

        % sum over h~=k,w of q(M,w,k,h)*QN(M,k)
        for h = 1:K
            if h ~= k
                for w = 1:M
                    qval = q(M, w, k, h);
                    row_coeffs(getQNIdx(M, k)) = row_coeffs(getQNIdx(M, k)) + qval;
                end
            end
        end

        % + sum over m,j<M of q(M,j,m,k)*UN(M,m)
        for m = 1:K
            for j = 1:(M-1)
                qval = q(M, j, m, k);
                row_coeffs(getUNIdx(M, m)) = row_coeffs(getUNIdx(M, m)) + qval;
            end
        end

        % - sum over j<M of q(j,M,k,k)*UN(j,k)
        for j = 1:(M-1)
            qval = q(j, M, k, k);
            row_coeffs(getUNIdx(j, k)) = row_coeffs(getUNIdx(j, k)) - qval;
        end

        % - sum over h~=k,w of q(M,w,h,k)*QN(M,h)
        for h = 1:K
            if h ~= k
                for w = 1:M
                    qval = q(M, w, h, k);
                    row_coeffs(getQNIdx(M, h)) = row_coeffs(getQNIdx(M, h)) - qval;
                end
            end
        end

        idx_nz = find(row_coeffs);
        eq_rows = [eq_rows; nEq * ones(length(idx_nz), 1)];
        eq_cols = [eq_cols; idx_nz(:)];
        eq_vals = [eq_vals; row_coeffs(idx_nz)'];
        beq = [beq; 0];
    end

    %% 7. MCC: for each i
    if verbose; fprintf('  MCC constraints...\n'); end
    for i = 1:M
        nEq = nEq + 1;
        row_coeffs = zeros(1, nVars);

        % sum over k,m,w~=i of q(i,w,k,m)*QN(i,k)
        for k = 1:K
            for m = 1:K
                for w = 1:M
                    if w ~= i
                        qval = q(i, w, k, m);
                        row_coeffs(getQNIdx(i, k)) = row_coeffs(getQNIdx(i, k)) + qval;
                    end
                end
            end
        end

        % + sum over k,m,j~=i of q(j,i,k,m)*QN(j,k)
        for k = 1:K
            for m = 1:K
                for j = 1:M
                    if j ~= i
                        qval = q(j, i, k, m);
                        row_coeffs(getQNIdx(j, k)) = row_coeffs(getQNIdx(j, k)) + qval;
                    end
                end
            end
        end

        % + sum over k,m,j~=i,wp~=i,wp~=j of q(j,i,k,m)*B(j,k,wp)
        for k = 1:K
            for m = 1:K
                for j = 1:M
                    if j ~= i
                        qval = q(j, i, k, m);
                        for wp = 1:M
                            if wp ~= i && wp ~= j
                                row_coeffs(getBIdx(j, k, wp)) = row_coeffs(getBIdx(j, k, wp)) + qval;
                            end
                        end
                    end
                end
            end
        end

        % - (N+1)*sum over k,m,j~=i of q(j,i,k,m)*UN(j,k)
        for k = 1:K
            for m = 1:K
                for j = 1:M
                    if j ~= i
                        qval = q(j, i, k, m);
                        row_coeffs(getUNIdx(j, k)) = row_coeffs(getUNIdx(j, k)) - (N + 1) * qval;
                    end
                end
            end
        end

        idx_nz = find(row_coeffs);
        eq_rows = [eq_rows; nEq * ones(length(idx_nz), 1)];
        eq_cols = [eq_cols; idx_nz(:)];
        eq_vals = [eq_vals; row_coeffs(idx_nz)'];
        beq = [beq; 0];
    end

    %% 8. MCC2: for each i
    if verbose; fprintf('  MCC2 constraints...\n'); end
    for i = 1:M
        nEq = nEq + 1;
        row_coeffs = zeros(1, nVars);

        % sum over k,m,w~=i of q(i,w,k,m)*QN(i,k)
        for k = 1:K
            for m = 1:K
                for w = 1:M
                    if w ~= i
                        qval = q(i, w, k, m);
                        row_coeffs(getQNIdx(i, k)) = row_coeffs(getQNIdx(i, k)) + qval;
                    end
                end
            end
        end

        % - sum over k,m,j~=i of q(j,i,k,m)*(B(j,k,i) + UN(j,k))
        for k = 1:K
            for m = 1:K
                for j = 1:M
                    if j ~= i
                        qval = q(j, i, k, m);
                        row_coeffs(getBIdx(j, k, i)) = row_coeffs(getBIdx(j, k, i)) - qval;
                        row_coeffs(getUNIdx(j, k)) = row_coeffs(getUNIdx(j, k)) - qval;
                    end
                end
            end
        end

        idx_nz = find(row_coeffs);
        eq_rows = [eq_rows; nEq * ones(length(idx_nz), 1)];
        eq_cols = [eq_cols; idx_nz(:)];
        eq_vals = [eq_vals; row_coeffs(idx_nz)'];
        beq = [beq; 0];
    end

    %% 9. QMAX: QN(w,k) - N*UN(w,k) <= 0 for all (w,k)
    if verbose; fprintf('  QMAX constraints...\n'); end
    for w = 1:M
        for k = 1:K
            nIneq = nIneq + 1;
            ineq_rows = [ineq_rows; nIneq; nIneq];
            ineq_cols = [ineq_cols; getQNIdx(w, k); getUNIdx(w, k)];
            ineq_vals = [ineq_vals; 1; -N];
            bineq = [bineq; 0];
        end
    end

    %% 10. QMIN: sum_w QN(w,k) - N*UN(j,k) >= 0 for each (k,j)
    %     Written as: -sum_w QN(w,k) + N*UN(j,k) <= 0
    if verbose; fprintf('  QMIN constraints...\n'); end
    for k = 1:K
        for j = 1:M
            nIneq = nIneq + 1;
            for w = 1:M
                ineq_rows = [ineq_rows; nIneq];
                ineq_cols = [ineq_cols; getQNIdx(w, k)];
                ineq_vals = [ineq_vals; -1];
            end
            ineq_rows = [ineq_rows; nIneq];
            ineq_cols = [ineq_cols; getUNIdx(j, k)];
            ineq_vals = [ineq_vals; N];
            bineq = [bineq; 0];
        end
    end

    %% Assemble sparse constraint matrices
    if verbose; fprintf('Assembling sparse matrices...\n'); end

    if nEq > 0
        Aeq = sparse(eq_rows, eq_cols, eq_vals, nEq, nVars);
    else
        Aeq = sparse(0, nVars);
        beq = [];
    end

    if nIneq > 0
        Aineq = sparse(ineq_rows, ineq_cols, ineq_vals, nIneq, nVars);
    else
        Aineq = sparse(0, nVars);
        bineq = [];
    end

    %% Build objective function
    if verbose; fprintf('Building objective function...\n'); end
    switch upper(objective_var)
        case 'UN'
            getObjIdx = getUNIdx;
        case 'QN'
            getObjIdx = getQNIdx;
        otherwise
            error('mapqn_bnd_lr_mva:objectiveVar', ...
                'objective_var must be ''UN'' or ''QN'', got ''%s''.', objective_var);
    end
    if objective_queue < 1 || objective_queue > M
        error('mapqn_bnd_lr_mva:objectiveQueue', ...
            'objective_queue must be in 1..%d, got %d.', M, objective_queue);
    end
    if objective_level < 0 || objective_level > K
        error('mapqn_bnd_lr_mva:objectiveLevel', ...
            'objective_level must be in 0..%d (0 aggregates over levels), got %d.', ...
            K, objective_level);
    end
    c = zeros(nVars, 1);
    if objective_level == 0
        for k = 1:K
            c(getObjIdx(objective_queue, k)) = 1;
        end
    else
        c(getObjIdx(objective_queue, objective_level)) = 1;
    end

    if strcmp(sense, 'max')
        c = -c;
    end

    %% Solve LP
    if verbose
        fprintf('Solving LP with %d variables and %d equality + %d inequality constraints...\n', ...
            nVars, nEq, nIneq);
    end

    % LP algorithm. R2025a's default 'dual-simplex-highs' is broken in some
    % installs (errors "Unrecognized field name optimstatus"), so an
    % interior-point variant is required. On the badly-scaled QRF instances
    % 'interior-point-legacy' is markedly more accurate than 'interior-point';
    % on the paper's BAS model it cuts the error against the published GLPK
    % optimum from 4.3e-05/8.0e-06 to 8.5e-07/9.2e-10. The residual is solver
    % tolerance, not formulation: the minimum lands above and the maximum below
    % the GLPK vertex, and both gaps shrink together as accuracy rises.
    % Override via params.lpAlgorithm.
    if isfield(params, 'lpAlgorithm') && ~isempty(params.lpAlgorithm)
        lpAlgorithm = params.lpAlgorithm;
    else
        lpAlgorithm = 'interior-point-legacy';
    end
    if verbose
        options = optimoptions('linprog', 'Display', 'final', 'Algorithm', lpAlgorithm);
    else
        options = optimoptions('linprog', 'Display', 'off', 'Algorithm', lpAlgorithm);
    end

    [x, fval, exitflag] = linprog(c, Aineq, bineq, Aeq, beq, lb, ub, options);

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
        warning('mapqn_bnd_lr_mva:nonFiniteSolution', ...
            'linprog returned a non-finite solution (exitflag %d); U and the other metric fields are left unpopulated.', ...
            exitflag);
    end

    if hasSolution
        % Extract UN and QN matrices
        result.UN = zeros(M, K);
        result.QN = zeros(M, K);

        for i = 1:M
            for k = 1:K
                result.UN(i, k) = x(getUNIdx(i, k));
                result.QN(i, k) = x(getQNIdx(i, k));
            end
        end
    else
        result.UN = [];
        result.QN = [];
    end

    if verbose; fprintf('\n=== Results ===\n'); end
    if verbose; fprintf('Objective value: %f\n', fval); end
    if verbose; fprintf('Exit flag: %d\n', exitflag); end
    if hasSolution && verbose
        fprintf('\nUtilizations (UN):\n');
        for i = 1:M
            fprintf('  Queue %d: [', i);
            for k = 1:K
                fprintf('%.6f ', result.UN(i, k));
            end
            fprintf(']\n');
        end
        fprintf('\nQueue lengths (QN):\n');
        for i = 1:M
            fprintf('  Queue %d: [', i);
            for k = 1:K
                fprintf('%.6f ', result.QN(i, k));
            end
            fprintf(']\n');
        end
    end
end
