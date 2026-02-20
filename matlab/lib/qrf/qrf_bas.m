function [result, x, fval, exitflag] = qrf_bas(params, objective, sense)
% QRF_BAS - Quadratic Reduction Framework for BAS (Blocking-After-Service) networks
%
% MATLAB port of the AMPL model qrboundsbas_skel.mod
%
% Usage:
%   [result, x, fval, exitflag] = qrf_bas(params)
%   [result, x, fval, exitflag] = qrf_bas(params, objective)
%   [result, x, fval, exitflag] = qrf_bas(params, objective, sense)
%
% Inputs:
%   params    - Structure with model parameters:
%               .M     - Number of queues
%               .N     - Total population
%               .f     - Index of finite capacity queue (1-based)
%               .F     - [M x 1] Capacity of each queue
%               .K     - [M x 1] Number of phases for each queue
%               .mu    - {M x 1} cell, each mu{i} is K(i) x K(i) completion rates
%               .v     - {M x 1} cell, each v{i} is K(i) x K(i) background rates
%               .r     - [M x M] Routing probabilities
%               .MR    - Number of blocking configurations
%               .BB    - [MR x M] Blocking state (0/1)
%               .MM    - [MR x 2] Blocking order (queue indices)
%               .ZZ    - [MR x 1] Number of blocked queues in each config
%               .ZM    - Maximum blocking depth
%               .MM1   - [MR x M] Extended blocking order info
%
%   objective - (optional) 'U1min' (default), 'U1max', or queue index 1..M
%   sense     - (optional) 'min' (default) or 'max'
%
% Outputs:
%   result   - Structure with results:
%              .U       - [M x 1] Utilization of each queue
%              .e       - [M x max(K)] Effective utilization by phase
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
    f = params.f;  % finite capacity queue index
    F = params.F(:);
    K = params.K(:);
    mu = params.mu;
    v = params.v;
    r = params.r;
    MR = params.MR;
    BB = params.BB;
    MM = params.MM;
    ZZ = params.ZZ(:);
    ZM = params.ZM;
    MM1 = params.MM1;

    % Compute transition rates q(i,j,k,h)
    % q{i,j} is a K(i) x K(i) array (not load-dependent for BAS)
    q = cell(M, M);
    for i = 1:M
        for j = 1:M
            q{i,j} = zeros(K(i), K(i));
            for ki = 1:K(i)
                for hi = 1:K(i)
                    if j ~= i
                        q{i,j}(ki, hi) = r(i,j) * mu{i}(ki, hi);
                    else
                        q{i,j}(ki, hi) = v{i}(ki, hi) + r(i,i) * mu{i}(ki, hi);
                    end
                end
            end
        end
    end

    %% Build variable indexing
    % p2(j, nj, kj, i, ni, hi, m) for j in 1:M, nj in 0:N, kj in 1:K(j),
    %                                 i in 1:M, ni in 0:N, hi in 1:K(i), m in 1:MR
    % e(i, ki) for i in 1:M, ki in 1:K(i)

    fprintf('Building variable index map...\n');
    varCount = 0;

    % p2 variables
    p2idx = cell(M, 1);
    for j = 1:M
        p2idx{j} = cell(N+1, K(j), M, MR);
        for nj = 0:N
            for kj = 1:K(j)
                for i = 1:M
                    for m = 1:MR
                        p2idx{j}{nj+1, kj, i, m} = zeros(N+1, K(i));
                        for ni = 0:N
                            for hi = 1:K(i)
                                varCount = varCount + 1;
                                p2idx{j}{nj+1, kj, i, m}(ni+1, hi) = varCount;
                            end
                        end
                    end
                end
            end
        end
    end

    % e variables
    eidx = cell(M, 1);
    for i = 1:M
        eidx{i} = zeros(K(i), 1);
        for ki = 1:K(i)
            varCount = varCount + 1;
            eidx{i}(ki) = varCount;
        end
    end

    nVars = varCount;
    fprintf('Total variables: %d\n', nVars);

    %% Build constraints
    Aeq = [];
    beq = [];
    Aineq = [];
    bineq = [];

    % Helper functions
    getP2Idx = @(j, nj, kj, i, ni, hi, m) p2idx{j}{nj+1, kj, i, m}(ni+1, hi);
    getEIdx = @(i, ki) eidx{i}(ki);

    fprintf('Building constraints...\n');

    %% Initialize bounds
    lb = zeros(nVars, 1);
    ub = inf(nVars, 1);

    %% ZERO constraints - fix infeasible states
    fprintf('  ZERO constraints...\n');
    for j = 1:M
        for nj = 0:N
            for kj = 1:K(j)
                for i = 1:M
                    for ni = 0:N
                        for hi = 1:K(i)
                            for m = 1:MR
                                idx = getP2Idx(j, nj, kj, i, ni, hi, m);

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

                                % ZERO6: nj > F(j)
                                if nj > F(j)
                                    ub(idx) = 0;
                                end

                                % ZERO5: BB(m,j)==1 and nj==0
                                if m >= 2 && BB(m, j) == 1 && nj == 0
                                    ub(idx) = 0;
                                end

                                % ZERO7: BB(m,j)==1 and i<>j and i<>f and ni+nj+F(f)>N
                                if m >= 2 && BB(m, j) == 1 && i ~= j && i ~= f && ni + nj + F(f) > N
                                    ub(idx) = 0;
                                end

                                % ZERO8: finite queue not at capacity in blocking config
                                if j == f && nj >= 1 && nj <= F(f)-1 && m >= 2
                                    ub(idx) = 0;
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    % ZERO4: For m>=2 and j<>f, p2(j,nj,k,f,nf,h,m)=0 when nf < F(f)
    for j = 1:M
        if j == f
            continue;
        end
        for nj = 0:N
            for kj = 1:K(j)
                for m = 2:MR
                    for nf = 0:(F(f)-1)
                        for hf = 1:K(f)
                            idx = getP2Idx(j, nj, kj, f, nf, hf, m);
                            ub(idx) = 0;
                        end
                    end
                end
            end
        end
    end

    %% ONE: Normalization
    fprintf('  ONE constraints...\n');
    for j = 1:M
        row = zeros(1, nVars);
        for nj = 0:N
            for kj = 1:K(j)
                for m = 1:MR
                    idx = getP2Idx(j, nj, kj, j, nj, kj, m);
                    row(idx) = 1;
                end
            end
        end
        Aeq = [Aeq; row];
        beq = [beq; 1];
    end

    %% SYMMETRY
    fprintf('  SYMMETRY constraints...\n');
    for j = 1:M
        for nj = 0:min(N, F(j))
            for kj = 1:K(j)
                for i = 1:M
                    if i <= j
                        continue;
                    end
                    for ni = 0:min(N, F(i))
                        if i ~= j && nj + ni > N
                            continue;
                        end
                        for hi = 1:K(i)
                            for m = 1:MR
                                idx1 = getP2Idx(j, nj, kj, i, ni, hi, m);
                                idx2 = getP2Idx(i, ni, hi, j, nj, kj, m);
                                if ub(idx1) == 0 && ub(idx2) == 0
                                    continue;
                                end
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
    end

    %% MARGINALS
    fprintf('  MARGINALS constraints...\n');
    for j = 1:M
        for kj = 1:K(j)
            for nj = 0:min(N, F(j))
                for i = 1:M
                    if i == j
                        continue;
                    end
                    for m = 1:MR
                        row = zeros(1, nVars);
                        idx_diag = getP2Idx(j, nj, kj, j, nj, kj, m);
                        row(idx_diag) = 1;
                        for ni = 0:min(N-nj, F(i))
                            for hi = 1:K(i)
                                idx = getP2Idx(j, nj, kj, i, ni, hi, m);
                                row(idx) = row(idx) - 1;
                            end
                        end
                        Aeq = [Aeq; row];
                        beq = [beq; 0];
                    end
                end
            end
        end
    end

    %% UEFF: e(i,ki) = sum of p2 where queue i is not blocked
    fprintf('  UEFF constraints...\n');
    for i = 1:M
        for ki = 1:K(i)
            row = zeros(1, nVars);
            idx_e = getEIdx(i, ki);
            row(idx_e) = -1;
            for j = 1:M
                for nj = 0:min(N, F(j))
                    for kj = 1:K(j)
                        for m = 1:MR
                            if BB(m, i) == 0
                                for ni = 1:min(N, F(i))
                                    idx = getP2Idx(j, nj, kj, i, ni, ki, m);
                                    row(idx) = row(idx) + 1;
                                end
                            end
                        end
                    end
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% THM1: Phase balance (Theorem 1)
    % sum {j, h: j<>i or h<>k} q(i,j,k,h)*e(i,k) = sum {j, h: j<>i or h<>k} q(i,j,h,k)*e(i,h)
    fprintf('  THM1 (Phase balance) constraints...\n');
    for i = 1:M
        for ki = 1:K(i)
            row = zeros(1, nVars);
            % LHS
            for j = 1:M
                for hi = 1:K(i)
                    if j ~= i || hi ~= ki
                        idx_e = getEIdx(i, ki);
                        row(idx_e) = row(idx_e) + q{i,j}(ki, hi);
                    end
                end
            end
            % RHS (subtract)
            for j = 1:M
                for hi = 1:K(i)
                    if j ~= i || hi ~= ki
                        idx_e = getEIdx(i, hi);
                        row(idx_e) = row(idx_e) - q{i,j}(hi, ki);
                    end
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% THM2: Population constraint (Theorem 2)
    fprintf('  THM2 (Population) constraints...\n');
    for j = 1:M
        for kj = 1:K(j)
            for nj = 0:F(j)
                for m = 1:MR
                    row = zeros(1, nVars);
                    % RHS: -N * p2(j,nj,kj,j,nj,kj,m)
                    idx_diag = getP2Idx(j, nj, kj, j, nj, kj, m);
                    row(idx_diag) = -N;
                    % LHS: sum
                    for i = 1:M
                        for ni = 1:F(i)
                            for ki = 1:K(i)
                                idx = getP2Idx(j, nj, kj, i, ni, ki, m);
                                row(idx) = row(idx) + ni;
                            end
                        end
                    end
                    Aeq = [Aeq; row];
                    beq = [beq; 0];
                end
            end
        end
    end

    %% COR1: Second moment constraint (Corollary to Theorem 2)
    % sum_{m,i,j,nj,ni,ki,kj} ni*nj*p2(j,nj,kj,i,ni,ki,m) = N^2
    fprintf('  COR1 (Second moment) constraint...\n');
    row = zeros(1, nVars);
    for m = 1:MR
        for i = 1:M
            for j = 1:M
                for nj = 1:F(j)
                    for ni = 1:F(i)
                        for ki = 1:K(i)
                            for kj = 1:K(j)
                                idx = getP2Idx(j, nj, kj, i, ni, ki, m);
                                row(idx) = row(idx) + ni * nj;
                            end
                        end
                    end
                end
            end
        end
    end
    Aeq = [Aeq; row];
    beq = [beq; N^2];

    %% THM30: Marginal balance for ni=0 (per phase), i<>f
    fprintf('  THM30 (Marginal balance ni=0) constraints...\n');
    for i = 1:M
        if i == f
            continue;
        end
        for ui = 1:K(i)
            row = zeros(1, nVars);
            % LHS: arrivals from j<>i,j<>f with BB(m,j)==0
            for j = 1:M
                if j == i || j == f
                    continue;
                end
                for nj = 1:F(j)
                    for kj = 1:K(j)
                        for hj = 1:K(j)
                            for m = 1:MR
                                if BB(m, j) == 0
                                    idx = getP2Idx(j, nj, kj, i, 0, ui, m);
                                    row(idx) = row(idx) + q{j,i}(kj, hj);
                                end
                            end
                        end
                    end
                end
            end
            % LHS: arrivals from j==f with MM(m,1)<>i
            for nj = 1:F(f)
                for kj = 1:K(f)
                    for hj = 1:K(f)
                        for m = 1:MR
                            if MM(m, 1) ~= i
                                idx = getP2Idx(f, nj, kj, i, 0, ui, m);
                                row(idx) = row(idx) + q{f,i}(kj, hj);
                            end
                        end
                    end
                end
            end

            % RHS: departures from i at ni=1 to j<>i,j<>f with BB(m,i)==0
            for j = 1:M
                if j == i || j == f
                    continue;
                end
                for nj = 0:F(j)
                    for ki = 1:K(i)
                        for hj = 1:K(j)
                            for m = 1:MR
                                if BB(m, i) == 0
                                    idx = getP2Idx(j, nj, hj, i, 1, ki, m);
                                    row(idx) = row(idx) - q{i,j}(ki, ui);
                                end
                            end
                        end
                    end
                end
            end
            % RHS: departures to j==f with BB(m,i)==0 and nj<F(f)
            for nj = 0:(F(f)-1)
                for ki = 1:K(i)
                    for hj = 1:K(f)
                        for m = 1:MR
                            if BB(m, i) == 0
                                idx = getP2Idx(f, nj, hj, i, 1, ki, m);
                                row(idx) = row(idx) - q{i,f}(ki, ui);
                            end
                        end
                    end
                end
            end
            % RHS: unblocking when BB(m,i)==1 and MM(m,1)==i
            for m = 1:MR
                if BB(m, i) == 1 && MM(m, 1) == i
                    for kf = 1:K(f)
                        for pf = 1:K(f)
                            for w = 1:M
                                if w ~= f && w ~= i
                                    idx = getP2Idx(f, F(f), kf, i, 1, ui, m);
                                    row(idx) = row(idx) - q{f,w}(kf, pf);
                                end
                            end
                        end
                    end
                end
            end

            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% THM3: Marginal balance for ni in 1:F(i)-1, i<>f
    fprintf('  THM3 (Marginal balance) constraints...\n');
    for i = 1:M
        if i == f
            continue;
        end
        for ni = 1:(F(i)-1)
            row = zeros(1, nVars);
            % LHS: arrivals from j<>i,j<>f with BB(m,j)==0
            for j = 1:M
                if j == i || j == f
                    continue;
                end
                for nj = 1:F(j)
                    for kj = 1:K(j)
                        for hj = 1:K(j)
                            for ui = 1:K(i)
                                for m = 1:MR
                                    if BB(m, j) == 0
                                        idx = getP2Idx(j, nj, kj, i, ni, ui, m);
                                        row(idx) = row(idx) + q{j,i}(kj, hj);
                                    end
                                end
                            end
                        end
                    end
                end
            end
            % LHS: arrivals from j==f with MM(m,1)<>i
            for nj = 1:F(f)
                for kj = 1:K(f)
                    for hj = 1:K(f)
                        for ui = 1:K(i)
                            for m = 1:MR
                                if MM(m, 1) ~= i
                                    idx = getP2Idx(f, nj, kj, i, ni, ui, m);
                                    row(idx) = row(idx) + q{f,i}(kj, hj);
                                end
                            end
                        end
                    end
                end
            end

            % RHS: departures from i at ni+1 to j<>i,j<>f with BB(m,i)==0
            for j = 1:M
                if j == i || j == f
                    continue;
                end
                for nj = 0:F(j)
                    for ki = 1:K(i)
                        for hi = 1:K(i)
                            for uj = 1:K(j)
                                for m = 1:MR
                                    if BB(m, i) == 0
                                        idx = getP2Idx(j, nj, uj, i, ni+1, ki, m);
                                        row(idx) = row(idx) - q{i,j}(ki, hi);
                                    end
                                end
                            end
                        end
                    end
                end
            end
            % RHS: departures to j==f with BB(m,i)==0 and nj<F(f)
            for nj = 0:(F(f)-1)
                for ki = 1:K(i)
                    for hi = 1:K(i)
                        for uj = 1:K(f)
                            for m = 1:MR
                                if BB(m, i) == 0
                                    idx = getP2Idx(f, nj, uj, i, ni+1, ki, m);
                                    row(idx) = row(idx) - q{i,f}(ki, hi);
                                end
                            end
                        end
                    end
                end
            end
            % RHS: unblocking when BB(m,i)==1 and MM(m,1)==i
            for m = 1:MR
                if BB(m, i) == 1 && MM(m, 1) == i
                    for ki = 1:K(i)
                        for kf = 1:K(f)
                            for pf = 1:K(f)
                                for w = 1:M
                                    if w ~= f && w ~= i
                                        idx = getP2Idx(f, F(f), kf, i, ni+1, ki, m);
                                        row(idx) = row(idx) - q{f,w}(kf, pf);
                                    end
                                end
                            end
                        end
                    end
                end
            end

            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% THM3f: Marginal balance for i==f
    fprintf('  THM3f (Marginal balance for finite queue) constraints...\n');
    for ni = 0:(F(f)-1)
        row = zeros(1, nVars);
        % LHS: arrivals from j<>f with BB(m,j)==0 and ni<F(f)
        if ni < F(f)
            for j = 1:M
                if j == f
                    continue;
                end
                for nj = 1:F(j)
                    for kj = 1:K(j)
                        for hj = 1:K(j)
                            for uf = 1:K(f)
                                for m = 1:MR
                                    if BB(m, j) == 0
                                        idx = getP2Idx(j, nj, kj, f, ni, uf, m);
                                        row(idx) = row(idx) + q{j,f}(kj, hj);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        % RHS: departures from f at ni+1 to j<>f (only m=1, no blocking)
        for j = 1:M
            if j == f
                continue;
            end
            for nj = 0:F(j)
                for kf = 1:K(f)
                    for hf = 1:K(f)
                        for uj = 1:K(j)
                            if ni < F(f)
                                idx = getP2Idx(j, nj, uj, f, ni+1, kf, 1);
                                row(idx) = row(idx) - q{f,j}(kf, hf);
                            end
                        end
                    end
                end
            end
        end

        Aeq = [Aeq; row];
        beq = [beq; 0];
    end

    %% THM3I: Blocking depth balance (Theorem 4)
    fprintf('  THM3I (Blocking depth balance) constraints...\n');
    for z = 0:(ZM-1)
        row = zeros(1, nVars);
        % LHS: arrivals to f at F(f) from j<>f with BB(m,j)==0 and ZZ(m)==z
        for j = 1:M
            if j == f
                continue;
            end
            for nj = 1:F(j)
                for kj = 1:K(j)
                    for hj = 1:K(j)
                        for uf = 1:K(f)
                            for m = 1:MR
                                if BB(m, j) == 0 && ZZ(m) == z
                                    idx = getP2Idx(j, nj, kj, f, F(f), uf, m);
                                    row(idx) = row(idx) + q{j,f}(kj, hj);
                                end
                            end
                        end
                    end
                end
            end
        end
        % RHS: departures from f to j<>f with ZZ(m)==z+1
        for j = 1:M
            if j == f
                continue;
            end
            for nj = 0:F(j)
                for kf = 1:K(f)
                    for hf = 1:K(f)
                        for uj = 1:K(j)
                            for m = 1:MR
                                if ZZ(m) == z+1
                                    idx = getP2Idx(j, nj, uj, f, F(f), kf, m);
                                    row(idx) = row(idx) - q{f,j}(kf, hf);
                                end
                            end
                        end
                    end
                end
            end
        end

        Aeq = [Aeq; row];
        beq = [beq; 0];
    end

    %% THM3L: Maximum blocking depth constraint
    fprintf('  THM3L (Max blocking depth) constraints...\n');
    for m = 1:MR
        if ZZ(m) ~= ZM - 1
            continue;
        end
        row = zeros(1, nVars);
        % LHS: arrivals from j<>f with BB(m,j)==0 and MM1(m,j)>0
        for j = 1:M
            if j == f || BB(m, j) ~= 0 || MM1(m, j) <= 0
                continue;
            end
            for nj = 1:F(j)
                for kj = 1:K(j)
                    for hj = 1:K(j)
                        for uf = 1:K(f)
                            idx = getP2Idx(j, nj, kj, f, F(f), uf, m);
                            row(idx) = row(idx) + q{j,f}(kj, hj);
                        end
                    end
                end
            end
        end
        % RHS: uses MM1(m,j) to index into blocking configuration
        for j = 1:M
            if j == f || BB(m, j) ~= 0 || MM1(m, j) <= 0
                continue;
            end
            mp = MM1(m, j);  % blocking configuration index
            for kf = 1:K(f)
                for uf = 1:K(f)
                    for w = 1:M
                        if w ~= f
                            idx = getP2Idx(f, F(f), kf, f, F(f), kf, mp);
                            row(idx) = row(idx) - q{f,w}(kf, uf);
                        end
                    end
                end
            end
        end

        Aeq = [Aeq; row];
        beq = [beq; 0];
    end

    %% THM4: Queue-length bound inequality (Theorem 5)
    fprintf('  THM4 (Queue-length bound) constraints...\n');
    for j = 1:M
        for kj = 1:K(j)
            for i = 1:M
                for m = 1:MR
                    row = zeros(1, nVars);
                    % LHS: sum_t sum_ht sum_nj sum_nt nt * p2
                    for t = 1:M
                        for ht = 1:K(t)
                            for nj = 0:F(j)
                                for nt = 1:F(t)
                                    idx = getP2Idx(j, nj, kj, t, nt, ht, m);
                                    row(idx) = row(idx) + nt;
                                end
                            end
                        end
                    end
                    % RHS: -N * sum
                    for hi = 1:K(i)
                        for nj = 0:F(j)
                            for ni = 1:F(i)
                                idx = getP2Idx(j, nj, kj, i, ni, hi, m);
                                row(idx) = row(idx) - N;
                            end
                        end
                    end
                    Aineq = [Aineq; -row];
                    bineq = [bineq; 0];
                end
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

    % Utilization = sum over m, k, n of p2(i,n,k,i,n,k,m)
    for m = 1:MR
        for ki = 1:K(targetQueue)
            for ni = 1:F(targetQueue)
                idx = getP2Idx(targetQueue, ni, ki, targetQueue, ni, ki, m);
                c(idx) = 1;
            end
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
        % Compute utilizations
        result.U = zeros(M, 1);
        result.e = zeros(M, max(K));

        for i = 1:M
            for ki = 1:K(i)
                result.e(i, ki) = x(getEIdx(i, ki));
            end

            for m = 1:MR
                for ki = 1:K(i)
                    for ni = 1:F(i)
                        idx = getP2Idx(i, ni, ki, i, ni, ki, m);
                        result.U(i) = result.U(i) + x(idx);
                    end
                end
            end
        end
    end

    fprintf('\n=== Results ===\n');
    fprintf('Objective value: %f\n', fval);
    fprintf('Exit flag: %d\n', exitflag);
    if exitflag > 0
        fprintf('\nUtilizations:\n');
        for i = 1:M
            fprintf('  Queue %d: U = %.6f\n', i, result.U(i));
        end
        fprintf('\nEffective utilizations by phase:\n');
        for i = 1:M
            fprintf('  Queue %d: e = [', i);
            for ki = 1:K(i)
                fprintf('%.6f ', result.e(i, ki));
            end
            fprintf(']\n');
        end
    end
end
