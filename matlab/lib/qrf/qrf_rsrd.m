function [result, x, fval, exitflag] = qrf_rsrd(params, objective, sense)
% QRF_RSRD - Quadratic Reduction Framework for RS-RD blocking networks
%
% MATLAB port of the AMPL model qrboundsrsrd_skel.mod
%
% Usage:
%   [result, x, fval, exitflag] = qrf_rsrd(params)
%   [result, x, fval, exitflag] = qrf_rsrd(params, objective)
%   [result, x, fval, exitflag] = qrf_rsrd(params, objective, sense)
%
% Inputs:
%   params    - Structure with model parameters:
%               .M     - Number of queues
%               .N     - Total population
%               .F     - [M x 1] Capacity of each queue
%               .K     - [M x 1] Number of phases for each queue
%               .mu    - {M x 1} cell, each mu{i} is K(i) x K(i) completion rates
%               .v     - {M x 1} cell, each v{i} is K(i) x K(i) background rates
%               .r     - [M x M] Routing probabilities
%               .alpha - (optional) {M x 1} cell, each alpha{i} is [N x 1] load-dependent rates
%
%   objective - (optional) 'U1min' (default), 'U1max', or queue index 1..M
%   sense     - (optional) 'min' (default) or 'max'
%
% Outputs:
%   result   - Structure with results:
%              .U       - [M x 1] Utilization of each queue
%              .Ueff    - [M x 1] Effective utilization
%              .pb      - [M x 1] Blocking probability
%              .p2      - Decision variable tensor (marginal probabilities)
%   x        - Raw solution vector
%   fval     - Objective function value
%   exitflag - Solver exit flag

    if nargin < 2 || isempty(objective)
        objective = 'U1min';
    end
    if nargin < 3 || isempty(sense)
        sense = 'min';
    end

    % Extract parameters
    M = params.M;
    N = params.N;
    F = params.F(:);
    K = params.K(:);
    mu = params.mu;
    v = params.v;
    r = params.r;

    % Default alpha (load-independent)
    if isfield(params, 'alpha') && ~isempty(params.alpha)
        alpha = params.alpha;
    else
        alpha = cell(M, 1);
        for i = 1:M
            alpha{i} = ones(N, 1);
        end
    end

    % Compute transition rates q(i,j,k,h,n)
    % q{i,j} is a K(i) x K(i) x (N+1) array
    q = cell(M, M);
    for i = 1:M
        for j = 1:M
            q{i,j} = zeros(K(i), K(i), N+1);
            for ki = 1:K(i)
                for hi = 1:K(i)
                    for n = 1:N  % n=0 gives q=0
                        if j ~= i
                            q{i,j}(ki, hi, n+1) = r(i,j) * mu{i}(ki, hi) * alpha{i}(n);
                        else
                            q{i,j}(ki, hi, n+1) = (v{i}(ki, hi) + r(i,i) * mu{i}(ki, hi)) * alpha{i}(n);
                        end
                    end
                end
            end
        end
    end

    %% Build variable indexing
    % p2(j, nj, kj, i, ni, hi) for j in 1:M, nj in 0:F(j), kj in 1:K(j),
    %                              i in 1:M, ni in 0:F(i), hi in 1:K(i)

    % Count variables and build index map
    fprintf('Building variable index map...\n');
    varCount = 0;
    p2idx = cell(M, 1);
    for j = 1:M
        p2idx{j} = cell(F(j)+1, K(j), M);
        for nj = 0:F(j)
            for kj = 1:K(j)
                for i = 1:M
                    p2idx{j}{nj+1, kj, i} = zeros(F(i)+1, K(i));
                    for ni = 0:F(i)
                        for hi = 1:K(i)
                            varCount = varCount + 1;
                            p2idx{j}{nj+1, kj, i}(ni+1, hi) = varCount;
                        end
                    end
                end
            end
        end
    end
    nP2Vars = varCount;

    % U and Ueff variables: U(i,k,ni) and Ueff(i,k,ni) for ni >= 1
    Uidx = cell(M, 1);
    Ueffidx = cell(M, 1);
    for i = 1:M
        Uidx{i} = zeros(K(i), F(i));
        Ueffidx{i} = zeros(K(i), F(i));
        for ki = 1:K(i)
            for ni = 1:F(i)
                varCount = varCount + 1;
                Uidx{i}(ki, ni) = varCount;
                varCount = varCount + 1;
                Ueffidx{i}(ki, ni) = varCount;
            end
        end
    end

    nVars = varCount;
    fprintf('Total variables: %d (p2: %d, U/Ueff: %d)\n', nVars, nP2Vars, nVars - nP2Vars);

    %% Build constraints
    Aeq = [];
    beq = [];
    Aineq = [];
    bineq = [];

    % Helper to get variable index
    getIdx = @(j, nj, kj, i, ni, hi) p2idx{j}{nj+1, kj, i}(ni+1, hi);
    getUIdx = @(i, ki, ni) Uidx{i}(ki, ni);
    getUeffIdx = @(i, ki, ni) Ueffidx{i}(ki, ni);

    fprintf('Building constraints...\n');

    %% ONE: Normalization - sum over nj, kj equals 1 for each j
    fprintf('  ONE constraints...\n');
    for j = 1:M
        row = zeros(1, nVars);
        for nj = 0:F(j)
            for kj = 1:K(j)
                idx = getIdx(j, nj, kj, j, nj, kj);
                row(idx) = 1;
            end
        end
        Aeq = [Aeq; row];
        beq = [beq; 1];
    end

    %% ZERO constraints - fix infeasible states to zero
    fprintf('  ZERO constraints...\n');
    lb = zeros(nVars, 1);
    ub = ones(nVars, 1);
    % U and Ueff bounds are [0, 1]
    for i = 1:M
        for ki = 1:K(i)
            for ni = 1:F(i)
                ub(getUIdx(i, ki, ni)) = 1;
                ub(getUeffIdx(i, ki, ni)) = 1;
            end
        end
    end

    for j = 1:M
        for nj = 0:F(j)
            for kj = 1:K(j)
                for i = 1:M
                    for ni = 0:F(i)
                        for hi = 1:K(i)
                            idx = getIdx(j, nj, kj, i, ni, hi);

                            % ZERO1: i==j, nj==ni, h<>k
                            if i == j && nj == ni && hi ~= kj
                                ub(idx) = 0;
                            end

                            % ZERO2: i==j, nj<>ni
                            if i == j && nj ~= ni
                                ub(idx) = 0;
                            end

                            % ZERO3: i<>j, nj+ni > N
                            if i ~= j && nj + ni > N
                                ub(idx) = 0;
                            end

                            % ZERO6: i<>j, N-nj-ni > sum of other capacities
                            if i ~= j
                                sumOtherF = sum(F) - F(i) - F(j);
                                if N - nj - ni > sumOtherF
                                    ub(idx) = 0;
                                end
                            end
                        end
                    end
                end
            end
        end

        % ZERO7: N-nj > sum of other capacities (for diagonal)
        for nj = 0:F(j)
            for kj = 1:K(j)
                sumOtherF = sum(F) - F(j);
                if N - nj > sumOtherF
                    idx = getIdx(j, nj, kj, j, nj, kj);
                    ub(idx) = 0;
                end
            end
        end
    end

    %% SYMMETRY: p2(i,ni,hi,j,nj,kj) = p2(j,nj,kj,i,ni,hi)
    fprintf('  SYMMETRY constraints...\n');
    for j = 1:M
        for nj = 0:F(j)
            for kj = 1:K(j)
                for i = 1:M
                    if i <= j
                        continue; % avoid duplicate constraints
                    end
                    for ni = 0:F(i)
                        for hi = 1:K(i)
                            idx1 = getIdx(j, nj, kj, i, ni, hi);
                            idx2 = getIdx(i, ni, hi, j, nj, kj);
                            if idx1 ~= idx2
                                row = zeros(1, nVars);
                                row(idx1) = 1;
                                row(idx2) = -1;
                                Aeq = [Aeq; row];
                                beq = [beq; 0];
                            end
                        end
                    end
                end
            end
        end
    end

    %% MARGINALS: p2(j,nj,kj,j,nj,kj) = sum over ni,hi of p2(j,nj,kj,i,ni,hi)
    fprintf('  MARGINALS constraints...\n');
    for j = 1:M
        for kj = 1:K(j)
            for nj = 0:F(j)
                for i = 1:M
                    if i == j
                        continue;
                    end
                    row = zeros(1, nVars);
                    idx_diag = getIdx(j, nj, kj, j, nj, kj);
                    row(idx_diag) = 1;
                    for ni = 0:F(i)
                        for hi = 1:K(i)
                            idx = getIdx(j, nj, kj, i, ni, hi);
                            row(idx) = row(idx) - 1;
                        end
                    end
                    Aeq = [Aeq; row];
                    beq = [beq; 0];
                end
            end
        end
    end

    %% UCLASSIC: U(i,k,ni) = p2(i,ni,k,i,ni,k)
    fprintf('  UCLASSIC constraints...\n');
    for i = 1:M
        for ki = 1:K(i)
            for ni = 1:F(i)
                row = zeros(1, nVars);
                idx_U = getUIdx(i, ki, ni);
                idx_p2 = getIdx(i, ni, ki, i, ni, ki);
                row(idx_U) = 1;
                row(idx_p2) = -1;
                Aeq = [Aeq; row];
                beq = [beq; 0];
            end
        end
    end

    %% UEFFS: Ueff(i,k,ni) = p2(i,ni,k,i,ni,k) - sum_{j: r(i,j)>0} r(i,j)*p2(i,ni,k,j,F(j),h)
    fprintf('  UEFFS constraints...\n');
    for i = 1:M
        for ki = 1:K(i)
            for ni = 1:F(i)
                row = zeros(1, nVars);
                idx_Ueff = getUeffIdx(i, ki, ni);
                idx_p2_diag = getIdx(i, ni, ki, i, ni, ki);
                row(idx_Ueff) = 1;
                row(idx_p2_diag) = -1;
                % Add blocking terms
                for j = 1:M
                    if j ~= i && r(i,j) > 0
                        for hj = 1:K(j)
                            idx_block = getIdx(i, ni, ki, j, F(j), hj);
                            row(idx_block) = r(i,j);
                        end
                    end
                end
                Aeq = [Aeq; row];
                beq = [beq; 0];
            end
        end
    end

    %% THM2: Phase balance (Theorem 1 - THM:sdeffective)
    % sum_{ni} (sum_{j<>i, h<>k} q(i,j,k,h,ni)*Ueff(i,k,ni) + sum_{h<>k} q(i,i,k,h,ni)*U(i,k,ni))
    % = sum_{ni} (sum_{j<>i, h<>k} q(i,j,h,k,ni)*Ueff(i,h,ni) + sum_{h<>k} q(i,i,h,k,ni)*U(i,h,ni))
    fprintf('  THM2 (Phase balance) constraints...\n');
    for i = 1:M
        for ki = 1:K(i)
            row = zeros(1, nVars);
            % LHS terms
            for ni = 1:F(i)
                % Terms with Ueff (j<>i)
                for j = 1:M
                    if j == i
                        continue;
                    end
                    for hi = 1:K(i)
                        if hi == ki
                            continue;
                        end
                        coef = q{i,j}(ki, hi, ni+1);
                        idx_Ueff = getUeffIdx(i, ki, ni);
                        row(idx_Ueff) = row(idx_Ueff) + coef;
                    end
                end
                % Terms with U (j==i, h<>k)
                for hi = 1:K(i)
                    if hi == ki
                        continue;
                    end
                    coef = q{i,i}(ki, hi, ni+1);
                    idx_p2 = getIdx(i, ni, ki, i, ni, ki);
                    row(idx_p2) = row(idx_p2) + coef;
                end
            end
            % RHS terms (subtract)
            for ni = 1:F(i)
                % Terms with Ueff (j<>i)
                for j = 1:M
                    if j == i
                        continue;
                    end
                    for hi = 1:K(i)
                        if hi == ki
                            continue;
                        end
                        coef = q{i,j}(hi, ki, ni+1);
                        idx_Ueff = getUeffIdx(i, hi, ni);
                        row(idx_Ueff) = row(idx_Ueff) - coef;
                    end
                end
                % Terms with U (j==i, h<>k)
                for hi = 1:K(i)
                    if hi == ki
                        continue;
                    end
                    coef = q{i,i}(hi, ki, ni+1);
                    idx_p2 = getIdx(i, ni, hi, i, ni, hi);
                    row(idx_p2) = row(idx_p2) - coef;
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% THM1: Population constraint (Theorem 2)
    % sum_i sum_ni sum_hi ni * p2(j,nj,kj,i,ni,hi) = N * p2(j,nj,kj,j,nj,kj)
    fprintf('  THM1 (Population) constraints...\n');
    for j = 1:M
        for kj = 1:K(j)
            for nj = 0:F(j)
                row = zeros(1, nVars);
                % RHS: -N * p2(j,nj,kj,j,nj,kj)
                idx_diag = getIdx(j, nj, kj, j, nj, kj);
                row(idx_diag) = -N;
                % LHS: sum
                for i = 1:M
                    for ni = 1:F(i)  % ni >= 1
                        for hi = 1:K(i)
                            idx = getIdx(j, nj, kj, i, ni, hi);
                            row(idx) = row(idx) + ni;
                        end
                    end
                end
                Aeq = [Aeq; row];
                beq = [beq; 0];
            end
        end
    end

    %% THM3a: Marginal balance for ni in 1:F(i)-1
    fprintf('  THM3a (Marginal balance) constraints...\n');
    for i = 1:M
        for ni = 1:(F(i)-1)
            row = zeros(1, nVars);
            % LHS: arrivals to queue i
            for j = 1:M
                if j == i
                    continue;
                end
                for kj = 1:K(j)
                    for hj = 1:K(j)
                        for ui = 1:K(i)
                            for nj = 1:F(j)
                                idx = getIdx(j, nj, kj, i, ni, ui);
                                row(idx) = row(idx) + q{j,i}(kj, hj, nj+1);
                            end
                        end
                    end
                end
            end
            % RHS: departures from queue i at ni+1
            for j = 1:M
                if j == i
                    continue;
                end
                for ki = 1:K(i)
                    for hi = 1:K(i)
                        for uj = 1:K(j)
                            for nj = 0:(F(j)-1)
                                idx = getIdx(i, ni+1, ki, j, nj, uj);
                                row(idx) = row(idx) - q{i,j}(ki, hi, ni+2);
                            end
                        end
                    end
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% THM3b: Marginal balance for ni=0, per phase
    fprintf('  THM3b (Marginal balance ni=0) constraints...\n');
    for i = 1:M
        for ui = 1:K(i)
            row = zeros(1, nVars);
            % LHS: arrivals to queue i at ni=0
            for j = 1:M
                if j == i
                    continue;
                end
                for kj = 1:K(j)
                    for hj = 1:K(j)
                        for nj = 1:F(j)
                            idx = getIdx(j, nj, kj, i, 0, ui);
                            row(idx) = row(idx) + q{j,i}(kj, hj, nj+1);
                        end
                    end
                end
            end
            % RHS: departures from queue i at ni=1
            for j = 1:M
                if j == i
                    continue;
                end
                for ki = 1:K(i)
                    for nj = 0:(F(j)-1)
                        for hj = 1:K(j)
                            idx = getIdx(i, 1, ki, j, nj, hj);
                            row(idx) = row(idx) - q{i,j}(ki, ui, 2);
                        end
                    end
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% QBAL: Queue balance constraint
    % This is a complex balance equation that tightens the bounds
    fprintf('  QBAL (Queue balance) constraints...\n');
    for i = 1:M
        for ki = 1:K(i)
            row = zeros(1, nVars);

            % LHS Term 1: sum{h<>k} sum{j<>i} sum{ni} sum{u} sum{nj} q[i,j,k,h,ni]*ni*p2[i,ni,k,j,nj,u]
            for hi = 1:K(i)
                if hi == ki
                    continue;
                end
                for j = 1:M
                    if j == i
                        continue;
                    end
                    for ni = 1:F(i)
                        for uj = 1:K(j)
                            for nj = 0:(F(j)-1)
                                coef = q{i,j}(ki, hi, ni+1) * ni;
                                idx = getIdx(i, ni, ki, j, nj, uj);
                                row(idx) = row(idx) + coef;
                            end
                        end
                    end
                end
            end

            % LHS Term 2: sum{h<>k} sum{ni} q[i,i,k,h,ni]*ni*p2[i,ni,k,i,ni,k]
            for hi = 1:K(i)
                if hi == ki
                    continue;
                end
                for ni = 1:F(i)
                    coef = q{i,i}(ki, hi, ni+1) * ni;
                    idx = getIdx(i, ni, ki, i, ni, ki);
                    row(idx) = row(idx) + coef;
                end
            end

            % LHS Term 3: sum{j<>i} sum{h} sum{ni} sum{u} sum{nj<=min(F(j)-1,N-ni)} q[i,j,h,k,ni]*p2[i,ni,h,j,nj,u]
            for j = 1:M
                if j == i
                    continue;
                end
                for hi = 1:K(i)
                    for ni = 1:F(i)
                        for uj = 1:K(j)
                            maxNj = min(F(j)-1, N-ni);
                            for nj = 0:maxNj
                                coef = q{i,j}(hi, ki, ni+1);
                                idx = getIdx(i, ni, hi, j, nj, uj);
                                row(idx) = row(idx) + coef;
                            end
                        end
                    end
                end
            end

            % RHS Term 1: sum{j<>i} sum{h} sum{ni<=F(i)-1} sum{u} sum{nj>=1} q[j,i,h,u,nj]*p2[i,ni,k,j,nj,h]
            for j = 1:M
                if j == i
                    continue;
                end
                for hj = 1:K(j)
                    for ni = 0:(F(i)-1)
                        for uj = 1:K(j)
                            for nj = 1:F(j)
                                coef = q{j,i}(hj, uj, nj+1);
                                idx = getIdx(i, ni, ki, j, nj, hj);
                                row(idx) = row(idx) - coef;
                            end
                        end
                    end
                end
            end

            % RHS Term 2: sum{h<>k} sum{ni} q[i,i,h,k,ni]*ni*p2[i,ni,h,i,ni,h]
            for hi = 1:K(i)
                if hi == ki
                    continue;
                end
                for ni = 1:F(i)
                    coef = q{i,i}(hi, ki, ni+1) * ni;
                    idx = getIdx(i, ni, hi, i, ni, hi);
                    row(idx) = row(idx) - coef;
                end
            end

            % RHS Term 3: sum{h<>k} sum{j<>i} sum{ni} sum{u} sum{nj} q[i,j,h,k,ni]*ni*p2[i,ni,h,j,nj,u]
            for hi = 1:K(i)
                if hi == ki
                    continue;
                end
                for j = 1:M
                    if j == i
                        continue;
                    end
                    for ni = 1:F(i)
                        for uj = 1:K(j)
                            for nj = 0:(F(j)-1)
                                coef = q{i,j}(hi, ki, ni+1) * ni;
                                idx = getIdx(i, ni, hi, j, nj, uj);
                                row(idx) = row(idx) - coef;
                            end
                        end
                    end
                end
            end

            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% THM4: Queue-length bound inequality (Theorem 5)
    % sum_t sum_ht sum_nj sum_nt nt*p2(j,nj,kj,t,nt,ht) >= N * sum_hi sum_nj sum_ni p2(j,nj,kj,i,ni,hi)
    fprintf('  THM4 (Queue-length bound) constraints...\n');
    for j = 1:M
        for kj = 1:K(j)
            for i = 1:M
                row = zeros(1, nVars);
                % LHS: sum_t sum_ht sum_nj sum_nt nt * p2
                for t = 1:M
                    for ht = 1:K(t)
                        for nj = 0:F(j)
                            for nt = 1:F(t)  % nt >= 1
                                idx = getIdx(j, nj, kj, t, nt, ht);
                                row(idx) = row(idx) + nt;
                            end
                        end
                    end
                end
                % RHS: -N * sum
                for hi = 1:K(i)
                    for nj = 0:F(j)
                        for ni = 1:F(i)  % ni >= 1
                            idx = getIdx(j, nj, kj, i, ni, hi);
                            row(idx) = row(idx) - N;
                        end
                    end
                end
                Aineq = [Aineq; -row];  % >= becomes <= with negation
                bineq = [bineq; 0];
            end
        end
    end

    %% Build objective function
    fprintf('Building objective function...\n');
    c = zeros(nVars, 1);

    if ischar(objective)
        if strcmp(objective, 'U1min') || strcmp(objective, 'U1max')
            targetQueue = 1;
        else
            error('Unknown objective: %s', objective);
        end
    else
        targetQueue = objective;
    end

    % Utilization = sum over k, n of p2(i,n,k,i,n,k)
    for ki = 1:K(targetQueue)
        for ni = 1:F(targetQueue)
            idx = getIdx(targetQueue, ni, ki, targetQueue, ni, ki);
            c(idx) = 1;
        end
    end

    if strcmp(sense, 'max')
        c = -c;
    end

    %% Solve LP
    fprintf('Solving LP with %d variables and %d equality + %d inequality constraints...\n', ...
        nVars, size(Aeq, 1), size(Aineq, 1));

    options = optimoptions('linprog', 'Display', 'final', 'Algorithm', 'interior-point');

    [x, fval, exitflag] = linprog(c, Aineq, bineq, Aeq, beq, lb, ub, options);

    if strcmp(sense, 'max')
        fval = -fval;
    end

    %% Extract results
    result = struct();
    result.objective = fval;
    result.exitflag = exitflag;

    if exitflag > 0
        % Compute utilizations from U and Ueff variables
        result.U = zeros(M, 1);
        result.Ueff = zeros(M, 1);
        result.pb = zeros(M, 1);

        for i = 1:M
            % U(i) = sum over k, ni of U(i,k,ni)
            for ki = 1:K(i)
                for ni = 1:F(i)
                    idx_U = getUIdx(i, ki, ni);
                    idx_Ueff = getUeffIdx(i, ki, ni);
                    result.U(i) = result.U(i) + x(idx_U);
                    result.Ueff(i) = result.Ueff(i) + x(idx_Ueff);
                end
            end
            result.pb(i) = result.U(i) - result.Ueff(i);
        end

        % Store p2 as a function handle for easy access
        result.getP2 = @(j, nj, kj, i, ni, hi) x(getIdx(j, nj, kj, i, ni, hi));
    end

    fprintf('\n=== Results ===\n');
    fprintf('Objective value: %f\n', fval);
    fprintf('Exit flag: %d\n', exitflag);
    if exitflag > 0
        fprintf('\nUtilizations:\n');
        for i = 1:M
            fprintf('  Queue %d: U = %.6f, Ueff = %.6f, pb = %.6f\n', ...
                i, result.U(i), result.Ueff(i), result.pb(i));
        end
    end
end
