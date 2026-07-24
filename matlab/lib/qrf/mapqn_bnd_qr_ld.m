function [result, x, fval, exitflag] = mapqn_bnd_qr_ld(params, objective_queue, objective_phase, objective_n, sense)
% MAPQN_BND_QR_LD - Quadratic Reduction Bounds for Load-Dependent Systems
%
% MATLAB port of bnd_qr_ld.py (Quadratic Reduction Bounds for Load-Dependent
% MAP queueing networks with load-dependent service rates).
%
% Usage:
%   [result, x, fval, exitflag] = mapqn_bnd_qr_ld(params, objective_queue, objective_phase, objective_n)
%   [result, x, fval, exitflag] = mapqn_bnd_qr_ld(params, objective_queue, objective_phase, objective_n, sense)
%
% Inputs:
%   params          - Structure with model parameters:
%                     .M     - Number of queues
%                     .N     - Total population
%                     .K     - [M x 1] Number of phases for each queue
%                     .mu    - {M x 1} cell, each mu{i} is K(i) x K(i) completion rates
%                     .v     - {M x 1} cell, each v{i} is K(i) x K(i) background rates
%                     .alpha - [M x N] load-dependent rates
%                     .r     - [M x M] Routing probabilities
%                     .verbose - (optional) boolean, default true
%
%   objective_queue - Queue index for objective (1-based, 1..M)
%   objective_phase - Phase index for objective (1-based, 1..K(objective_queue))
%   objective_n     - Population level for objective (0..N)
%   sense           - (optional) 'max' (default) or 'min'
%
% Outputs:
%   result   - Structure with results:
%              .objective    - Objective function value
%              .exitflag     - Solver exit flag
%              .p2marginals  - Marginal probabilities
%   x        - Raw solution vector
%   fval     - Objective function value
%   exitflag - Solver exit flag

    if nargin < 5 || isempty(sense)
        sense = 'max';
    end

    % Extract parameters
    M = params.M;
    N = params.N;
    K = params.K(:);
    mu = params.mu;
    v = params.v;
    alpha = params.alpha;
    r = params.r;
    if isfield(params, 'verbose')
        verbose = params.verbose;
    else
        verbose = true;
    end

    % Validate objective indices
    if objective_queue < 1 || objective_queue > M
        error('objective_queue must be in range 1..%d', M);
    end
    if objective_phase < 1 || objective_phase > K(objective_queue)
        error('objective_phase must be in range 1..%d', K(objective_queue));
    end
    if objective_n < 0 || objective_n > N
        error('objective_n must be in range 0..%d', N);
    end

    %% Define q function (1-based indices, load-dependent)
    % q_val(i,j,k,h,n): transition rate from queue i phase k to queue j phase h
    %                    when queue i has n customers
    q_val = @(i, j, k, h, n) q_func(i, j, k, h, n, mu, v, alpha, r);

    %% Build variable indexing
    % p2(j, nj, k, i, ni, h): joint probabilities
    % j in 1:M, nj in 0:N, k in 1:K(j), i in 1:M, ni in 0:N, h in 1:K(i)
    if verbose; fprintf('Building variable index map...\n'); end
    varCount = 0;

    % p2 variables
    p2idx = cell(M, 1);
    for j = 1:M
        p2idx{j} = cell(N+1, K(j), M);
        for nj = 0:N
            for kj = 1:K(j)
                for i = 1:M
                    p2idx{j}{nj+1, kj, i} = zeros(N+1, K(i));
                    for ni = 0:N
                        for hi = 1:K(i)
                            varCount = varCount + 1;
                            p2idx{j}{nj+1, kj, i}(ni+1, hi) = varCount;
                        end
                    end
                end
            end
        end
    end

    nVars = varCount;
    if verbose; fprintf('Total variables: %d\n', nVars); end

    %% Helper function
    getP2Idx = @(j, nj, kj, i, ni, hi) p2idx{j}{nj+1, kj, i}(ni+1, hi);

    %% Initialize bounds
    lb = zeros(nVars, 1);
    ub = ones(nVars, 1);

    %% Build constraints
    Aeq = [];
    beq = [];
    Aineq = [];
    bineq = [];

    if verbose; fprintf('Building constraints...\n'); end

    %% ZERO constraints - fix infeasible states via ub=0
    if verbose; fprintf('  ZERO constraints...\n'); end
    for j = 1:M
        for nj = 0:N
            for kj = 1:K(j)
                for i = 1:M
                    for ni = 0:N
                        for hi = 1:K(i)
                            idx = getP2Idx(j, nj, kj, i, ni, hi);

                            % ZERO1: i==j, nj==ni, h~=k
                            if i == j && nj == ni && hi ~= kj
                                ub(idx) = 0;
                            end

                            % ZERO2: i==j, nj~=ni
                            if i == j && nj ~= ni
                                ub(idx) = 0;
                            end

                            % ZERO3: i~=j, nj+ni > N
                            if i ~= j && nj + ni > N
                                ub(idx) = 0;
                            end
                        end
                    end
                end
            end
        end
    end

    %% ONE: Normalization
    % sum{nj,k} p2(j,nj,k,j,nj,k) = 1 for each j
    if verbose; fprintf('  ONE constraints...\n'); end
    for j = 1:M
        row = zeros(1, nVars);
        for nj = 0:N
            for kj = 1:K(j)
                idx = getP2Idx(j, nj, kj, j, nj, kj);
                row(idx) = 1;
            end
        end
        Aeq = [Aeq; row];
        beq = [beq; 1];
    end

    %% SYMMETRY: p2(i,ni,h,j,nj,k) = p2(j,nj,k,i,ni,h) (only one ordering)
    if verbose; fprintf('  SYMMETRY constraints...\n'); end
    for j = 1:M
        for nj = 0:N
            for kj = 1:K(j)
                for i = 1:M
                    if i <= j && ~(i == j)
                        % skip: we only do i > j for ordering
                    end
                    if i <= j
                        continue;
                    end
                    for ni = 0:N
                        if i ~= j && nj + ni > N
                            continue;
                        end
                        for hi = 1:K(i)
                            idx1 = getP2Idx(j, nj, kj, i, ni, hi);
                            idx2 = getP2Idx(i, ni, hi, j, nj, kj);
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

    %% MARGINALS: p2(j,nj,k,j,nj,k) = sum{ni=0..N-nj, h} p2(j,nj,k,i,ni,h) for i~=j
    if verbose; fprintf('  MARGINALS constraints...\n'); end
    for j = 1:M
        for kj = 1:K(j)
            for nj = 0:N
                for i = 1:M
                    if i == j
                        continue;
                    end
                    row = zeros(1, nVars);
                    idx_diag = getP2Idx(j, nj, kj, j, nj, kj);
                    row(idx_diag) = 1;
                    for ni = 0:(N - nj)
                        for hi = 1:K(i)
                            idx = getP2Idx(j, nj, kj, i, ni, hi);
                            row(idx) = row(idx) - 1;
                        end
                    end
                    Aeq = [Aeq; row];
                    beq = [beq; 0];
                end
            end
        end
    end

    %% THM1: Queue length theorem
    % for each (j,k): sum{i,nj>=1,ni>=1,h} ni*p2(j,nj,k,i,ni,h)
    %                  - N*sum{nj>=1} p2(j,nj,k,j,nj,k) = 0
    if verbose; fprintf('  THM1 (Queue length) constraints...\n'); end
    for j = 1:M
        for kj = 1:K(j)
            row = zeros(1, nVars);
            % LHS: sum over i, nj>=1, ni>=1, h of ni*p2
            for i = 1:M
                for nj = 1:N
                    for ni = 1:N
                        for hi = 1:K(i)
                            idx = getP2Idx(j, nj, kj, i, ni, hi);
                            row(idx) = row(idx) + ni;
                        end
                    end
                end
            end
            % RHS: -N * sum{nj>=1} p2(j,nj,k,j,nj,k)
            for nj = 1:N
                idx = getP2Idx(j, nj, kj, j, nj, kj);
                row(idx) = row(idx) - N;
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% THM1c: Empty queue theorem
    % for each (j,k): sum{i,ni>=1,h} ni*p2(j,0,k,i,ni,h) - N*p2(j,0,k,j,0,k) = 0
    if verbose; fprintf('  THM1c (Empty queue) constraints...\n'); end
    for j = 1:M
        for kj = 1:K(j)
            row = zeros(1, nVars);
            % LHS
            for i = 1:M
                for ni = 1:N
                    for hi = 1:K(i)
                        idx = getP2Idx(j, 0, kj, i, ni, hi);
                        row(idx) = row(idx) + ni;
                    end
                end
            end
            % RHS
            idx = getP2Idx(j, 0, kj, j, 0, kj);
            row(idx) = row(idx) - N;
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% PC2: Second moment constraint
    % sum{i,j,ni>=1,nj>=1,h,k} nj*ni*p2(j,nj,k,i,ni,h) = N^2
    if verbose; fprintf('  PC2 (Second moment) constraint...\n'); end
    row = zeros(1, nVars);
    for i = 1:M
        for j = 1:M
            for ni = 1:N
                for nj = 1:N
                    for hi = 1:K(i)
                        for kj = 1:K(j)
                            idx = getP2Idx(j, nj, kj, i, ni, hi);
                            row(idx) = row(idx) + nj * ni;
                        end
                    end
                end
            end
        end
    end
    Aeq = [Aeq; row];
    beq = [beq; N^2];

    %% THM2: Phase balance
    % for each (i,k): sum{j,h: j~=i or h~=k, ni>=1}
    %   (q(i,j,k,h,ni)*p2(i,ni,k,i,ni,k) - q(i,j,h,k,ni)*p2(i,ni,h,i,ni,h)) = 0
    if verbose; fprintf('  THM2 (Phase balance) constraints...\n'); end
    for i = 1:M
        for ki = 1:K(i)
            row = zeros(1, nVars);
            for j = 1:M
                for hi = 1:K(i)
                    if ~(hi == ki && j == i)
                        for ni = 1:N
                            q_out = q_val(i, j, ki, hi, ni);
                            q_in = q_val(i, j, hi, ki, ni);
                            idx_k = getP2Idx(i, ni, ki, i, ni, ki);
                            idx_h = getP2Idx(i, ni, hi, i, ni, hi);
                            row(idx_k) = row(idx_k) + q_out;
                            row(idx_h) = row(idx_h) - q_in;
                        end
                    end
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% THM3a: Flow balance for ni = 1..N-1
    % for each (i, ni in 1..N-1):
    %   sum{j~=i, k, h, u, nj=1..N-ni} q(j,i,k,h,nj)*p2(j,nj,k,i,ni,u)
    %   - sum{j~=i, k, h} q(i,j,k,h,ni+1)*p2(i,ni+1,k,i,ni+1,k) = 0
    if verbose; fprintf('  THM3a (Flow balance ni=1..N-1) constraints...\n'); end
    for i = 1:M
        for ni = 1:(N-1)
            row = zeros(1, nVars);
            % Incoming flow
            for j = 1:M
                if j ~= i
                    for kj = 1:K(j)
                        for hj = 1:K(j)
                            for u = 1:K(i)
                                for nj = 1:(N - ni)
                                    qv = q_val(j, i, kj, hj, nj);
                                    idx = getP2Idx(j, nj, kj, i, ni, u);
                                    row(idx) = row(idx) + qv;
                                end
                            end
                        end
                    end
                end
            end
            % Outgoing flow
            for j = 1:M
                if j ~= i
                    for ki = 1:K(i)
                        for hi = 1:K(i)
                            qv = q_val(i, j, ki, hi, ni + 1);
                            idx = getP2Idx(i, ni + 1, ki, i, ni + 1, ki);
                            row(idx) = row(idx) - qv;
                        end
                    end
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% THM3b: Flow balance for ni = 0
    % for each (i, u):
    %   sum{j~=i, k, h, nj=1..N} q(j,i,k,h,nj)*p2(j,nj,k,i,0,u)
    %   - sum{j~=i, k} q(i,j,k,u,1)*p2(i,1,k,i,1,k) = 0
    if verbose; fprintf('  THM3b (Flow balance ni=0) constraints...\n'); end
    for i = 1:M
        for u = 1:K(i)
            row = zeros(1, nVars);
            % Incoming to empty queue
            for j = 1:M
                if j ~= i
                    for kj = 1:K(j)
                        for hj = 1:K(j)
                            for nj = 1:N
                                qv = q_val(j, i, kj, hj, nj);
                                idx = getP2Idx(j, nj, kj, i, 0, u);
                                row(idx) = row(idx) + qv;
                            end
                        end
                    end
                end
            end
            % Outgoing from single customer
            for j = 1:M
                if j ~= i
                    for ki = 1:K(i)
                        qv = q_val(i, j, ki, u, 1);
                        idx = getP2Idx(i, 1, ki, i, 1, ki);
                        row(idx) = row(idx) - qv;
                    end
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% QBAL: Queue balance (from AMPL bnd_quadraticreduction_ld.mod)
    % for each (i, k):
    %   LHS1: sum{h~=k, j, ni>=1} q(i,j,k,h,ni)*ni*p2(i,ni,k,i,ni,k)
    %   + LHS2: sum{j~=i, h, ni>=1} q(i,j,h,k,ni)*p2(i,ni,h,i,ni,h)
    %   = RHS1: sum{j~=i, u in K(j), w in K(j), nj>=1} q(j,i,u,w,nj)*(p2(j,nj,u,i,0,k) + sum{ni>=1} p2(i,ni,k,j,nj,u))
    %   + RHS2: sum{h~=k, j, ni>=1} q(i,j,h,k,ni)*ni*p2(i,ni,h,i,ni,h)
    if verbose; fprintf('  QBAL (Queue balance) constraints...\n'); end
    for i = 1:M
        for ki = 1:K(i)
            row = zeros(1, nVars);
            % LHS1: sum{h~=k, j, ni>=1} q(i,j,k,h,ni)*ni*p2(i,ni,k,i,ni,k)
            for hi = 1:K(i)
                if hi ~= ki
                    for j = 1:M
                        for ni = 1:N
                            qv = q_val(i, j, ki, hi, ni);
                            idx = getP2Idx(i, ni, ki, i, ni, ki);
                            row(idx) = row(idx) + qv * ni;
                        end
                    end
                end
            end
            % LHS2: sum{j~=i, h, ni>=1} q(i,j,h,k,ni)*p2(i,ni,h,i,ni,h)
            for j = 1:M
                if j ~= i
                    for hi = 1:K(i)
                        for ni = 1:N
                            qv = q_val(i, j, hi, ki, ni);
                            idx = getP2Idx(i, ni, hi, i, ni, hi);
                            row(idx) = row(idx) + qv;
                        end
                    end
                end
            end
            % -RHS1 part a: sum{j~=i, u in K(j), w in K(j), nj>=1} q(j,i,u,w,nj)*p2(j,nj,u,i,0,k)
            for j = 1:M
                if j ~= i
                    for u = 1:K(j)
                        for w = 1:K(j)
                            for nj = 1:N
                                qv = q_val(j, i, u, w, nj);
                                idx = getP2Idx(j, nj, u, i, 0, ki);
                                row(idx) = row(idx) - qv;
                            end
                        end
                    end
                end
            end
            % -RHS1 part b: sum{j~=i, u in K(j), w in K(j), nj>=1, ni>=1} q(j,i,u,w,nj)*p2(i,ni,k,j,nj,u)
            for j = 1:M
                if j ~= i
                    for u = 1:K(j)
                        for w = 1:K(j)
                            for nj = 1:N
                                qv = q_val(j, i, u, w, nj);
                                for ni = 1:N
                                    idx = getP2Idx(i, ni, ki, j, nj, u);
                                    row(idx) = row(idx) - qv;
                                end
                            end
                        end
                    end
                end
            end
            % -RHS2: sum{h~=k, j, ni>=1} q(i,j,h,k,ni)*ni*p2(i,ni,h,i,ni,h)
            for hi = 1:K(i)
                if hi ~= ki
                    for j = 1:M
                        for ni = 1:N
                            qv = q_val(i, j, hi, ki, ni);
                            idx = getP2Idx(i, ni, hi, i, ni, hi);
                            row(idx) = row(idx) - qv * ni;
                        end
                    end
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% COR1a: Correlation constraint (from AMPL bnd_quadraticreduction_ld.mod)
    % for each (i, kstar, ni_c in 0..N-2)
    if verbose; fprintf('  COR1a constraints...\n'); end
    for i = 1:M
        for kstar = 1:K(i)
            for ni_c = 0:(N-2)
                row = zeros(1, nVars);
                % A: sum{j~=i,kj,hj,u~=kstar,nj=1..N-ni_c} q(j,i,kj,hj,nj)*p2(j,nj,kj,i,ni_c,u)
                for j = 1:M
                    if j ~= i
                        for kj = 1:K(j)
                            for hj = 1:K(j)
                                for u = 1:K(i)
                                    if u ~= kstar
                                        for nj = 1:(N-ni_c)
                                            qv = q_val(j, i, kj, hj, nj);
                                            idx = getP2Idx(j, nj, kj, i, ni_c, u);
                                            row(idx) = row(idx) + qv;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                % B: sum{j~=i,kj,hj,nj=1..N-ni_c} q(j,i,kj,hj,nj)*p2(j,nj,kj,i,ni_c+1,kstar)
                for j = 1:M
                    if j ~= i
                        for kj = 1:K(j)
                            for hj = 1:K(j)
                                for nj = 1:(N-ni_c)
                                    qv = q_val(j, i, kj, hj, nj);
                                    idx = getP2Idx(j, nj, kj, i, ni_c+1, kstar);
                                    row(idx) = row(idx) + qv;
                                end
                            end
                        end
                    end
                end
                % C: sum{k~=kstar} q(i,i,kstar,k,ni_c+1)*p2(i,ni_c+1,kstar,i,ni_c+1,kstar)
                for k2 = 1:K(i)
                    if k2 ~= kstar
                        qv = q_val(i, i, kstar, k2, ni_c+1);
                        idx = getP2Idx(i, ni_c+1, kstar, i, ni_c+1, kstar);
                        row(idx) = row(idx) + qv;
                    end
                end
                % -D: -sum{j~=i,k~=kstar} q(i,j,k,k,ni_c+1)*p2(i,ni_c+1,k,i,ni_c+1,k)
                for j = 1:M
                    if j ~= i
                        for k2 = 1:K(i)
                            if k2 ~= kstar
                                qv = q_val(i, j, k2, k2, ni_c+1);
                                idx = getP2Idx(i, ni_c+1, k2, i, ni_c+1, k2);
                                row(idx) = row(idx) - qv;
                            end
                        end
                    end
                end
                % -E: -sum{j~=i,k~=kstar,h~=k} q(i,j,k,h,ni_c+1)*p2(i,ni_c+1,k,i,ni_c+1,k)
                for j = 1:M
                    if j ~= i
                        for k2 = 1:K(i)
                            if k2 ~= kstar
                                for h2 = 1:K(i)
                                    if h2 ~= k2
                                        qv = q_val(i, j, k2, h2, ni_c+1);
                                        idx = getP2Idx(i, ni_c+1, k2, i, ni_c+1, k2);
                                        row(idx) = row(idx) - qv;
                                    end
                                end
                            end
                        end
                    end
                end
                % -F: -sum{j~=i,k~=kstar} q(i,j,k,kstar,ni_c+2)*p2(i,ni_c+2,k,i,ni_c+2,k)
                for j = 1:M
                    if j ~= i
                        for k2 = 1:K(i)
                            if k2 ~= kstar
                                qv = q_val(i, j, k2, kstar, ni_c+2);
                                idx = getP2Idx(i, ni_c+2, k2, i, ni_c+2, k2);
                                row(idx) = row(idx) - qv;
                            end
                        end
                    end
                end
                % -G: -sum{j~=i} q(i,j,kstar,kstar,ni_c+2)*p2(i,ni_c+2,kstar,i,ni_c+2,kstar)
                for j = 1:M
                    if j ~= i
                        qv = q_val(i, j, kstar, kstar, ni_c+2);
                        idx = getP2Idx(i, ni_c+2, kstar, i, ni_c+2, kstar);
                        row(idx) = row(idx) - qv;
                    end
                end
                % -H: -sum{k~=kstar} q(i,i,k,kstar,ni_c+1)*p2(i,ni_c+1,k,i,ni_c+1,k)
                for k2 = 1:K(i)
                    if k2 ~= kstar
                        qv = q_val(i, i, k2, kstar, ni_c+1);
                        idx = getP2Idx(i, ni_c+1, k2, i, ni_c+1, k2);
                        row(idx) = row(idx) - qv;
                    end
                end
                Aeq = [Aeq; row];
                beq = [beq; 0];
            end
        end
    end

    %% COR1b: Correlation constraint boundary (from AMPL bnd_quadraticreduction_ld.mod)
    % for each (i, kstar)
    if verbose; fprintf('  COR1b constraints...\n'); end
    for i = 1:M
        for kstar = 1:K(i)
            row = zeros(1, nVars);
            % A': sum{j~=i,kj,hj,u~=kstar} q(j,i,kj,hj,1)*p2(j,1,kj,i,N-1,u)
            for j = 1:M
                if j ~= i
                    for kj = 1:K(j)
                        for hj = 1:K(j)
                            for u = 1:K(i)
                                if u ~= kstar
                                    qv = q_val(j, i, kj, hj, 1);
                                    idx = getP2Idx(j, 1, kj, i, N-1, u);
                                    row(idx) = row(idx) + qv;
                                end
                            end
                        end
                    end
                end
            end
            % C': sum{k~=kstar} q(i,i,kstar,k,N)*p2(i,N,kstar,i,N,kstar)
            for k2 = 1:K(i)
                if k2 ~= kstar
                    qv = q_val(i, i, kstar, k2, N);
                    idx = getP2Idx(i, N, kstar, i, N, kstar);
                    row(idx) = row(idx) + qv;
                end
            end
            % -D': -sum{j~=i,k~=kstar} q(i,j,k,k,N)*p2(i,N,k,i,N,k)
            for j = 1:M
                if j ~= i
                    for k2 = 1:K(i)
                        if k2 ~= kstar
                            qv = q_val(i, j, k2, k2, N);
                            idx = getP2Idx(i, N, k2, i, N, k2);
                            row(idx) = row(idx) - qv;
                        end
                    end
                end
            end
            % -E': -sum{j~=i,k~=kstar,h~=k} q(i,j,k,h,N)*p2(i,N,k,i,N,k)
            for j = 1:M
                if j ~= i
                    for k2 = 1:K(i)
                        if k2 ~= kstar
                            for h2 = 1:K(i)
                                if h2 ~= k2
                                    qv = q_val(i, j, k2, h2, N);
                                    idx = getP2Idx(i, N, k2, i, N, k2);
                                    row(idx) = row(idx) - qv;
                                end
                            end
                        end
                    end
                end
            end
            % -H': -sum{k~=kstar} q(i,i,k,kstar,N)*p2(i,N,k,i,N,k)
            for k2 = 1:K(i)
                if k2 ~= kstar
                    qv = q_val(i, i, k2, kstar, N);
                    idx = getP2Idx(i, N, k2, i, N, k2);
                    row(idx) = row(idx) - qv;
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% THM4: QMIN inequality constraint
    % for each (j, k, i):
    %   sum{t, h, nj, nt} nt*p2(j,nj,k,t,nt,h) - N*sum{h,nj,ni} p2(j,nj,k,i,ni,h) >= 0
    if verbose; fprintf('  THM4 (QMIN) constraints...\n'); end
    for j = 1:M
        for kj = 1:K(j)
            for i = 1:M
                row = zeros(1, nVars);
                % LHS: sum{t,h,nj,nt} nt*p2(j,nj,k,t,nt,h)
                for t = 1:M
                    for ht = 1:K(t)
                        for nj = 0:N
                            for nt = 0:N
                                idx = getP2Idx(j, nj, kj, t, nt, ht);
                                row(idx) = row(idx) + nt;
                            end
                        end
                    end
                end
                % RHS: -N * sum{h,nj,ni} p2(j,nj,k,i,ni,h)
                for hi = 1:K(i)
                    for nj = 0:N
                        for ni = 0:N
                            idx = getP2Idx(j, nj, kj, i, ni, hi);
                            row(idx) = row(idx) - N;
                        end
                    end
                end
                % >= 0 means -row*x <= 0
                Aineq = [Aineq; -row];
                bineq = [bineq; 0];
            end
        end
    end

    %% Build objective function
    if verbose; fprintf('Building objective function...\n'); end
    c = zeros(nVars, 1);

    % Objective: p2(objective_queue, objective_n, objective_phase,
    %               objective_queue, objective_n, objective_phase)
    idx_obj = getP2Idx(objective_queue, objective_n, objective_phase, ...
                       objective_queue, objective_n, objective_phase);
    c(idx_obj) = 1;

    if strcmp(sense, 'max')
        c = -c;  % linprog minimizes, so negate for maximization
    end

    %% Solve LP
    if verbose
        fprintf('Solving LP with %d variables and %d equality + %d inequality constraints...\n', ...
            nVars, size(Aeq, 1), size(Aineq, 1));
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
        warning('mapqn_bnd_qr_ld:nonFiniteSolution', ...
            'linprog returned a non-finite solution (exitflag %d); U and the other metric fields are left unpopulated.', ...
            exitflag);
    end

    if hasSolution
        % Extract marginal probabilities p2(j,nj,k,j,nj,k)
        result.p2marginals = zeros(M, N+1, max(K));
        for j = 1:M
            for nj = 0:N
                for kj = 1:K(j)
                    idx = getP2Idx(j, nj, kj, j, nj, kj);
                    result.p2marginals(j, nj+1, kj) = x(idx);
                end
            end
        end
    end

    if verbose; fprintf('\n=== Results ===\n'); end
    if verbose; fprintf('Objective value: %f\n', fval); end
    if verbose; fprintf('Exit flag: %d\n', exitflag); end
    if hasSolution && verbose
        fprintf('\nMarginal probabilities p2(j,nj,k,j,nj,k):\n');
        for j = 1:M
            for nj = 0:N
                for kj = 1:K(j)
                    val = result.p2marginals(j, nj+1, kj);
                    if val > 1e-8
                        fprintf('  p2(%d,%d,%d) = %.6f\n', j, nj, kj, val);
                    end
                end
            end
        end
    end
end

function qv = q_func(i, j, k, h, n, mu, v, alpha, r)
% Q_FUNC - Compute load-dependent transition rate (1-based indices)
%
%   i: source queue (1-based)
%   j: destination queue (1-based)
%   k: source phase (1-based)
%   h: destination phase (1-based)
%   n: population at queue i (0 returns 0)
    if n == 0
        qv = 0;
        return;
    end

    % Load-dependent scaling factor
    if n <= size(alpha, 2)
        alpha_val = alpha(i, n);
    else
        alpha_val = 1;
    end

    if j ~= i
        qv = r(i, j) * mu{i}(k, h) * alpha_val;
    else
        qv = (v{i}(k, h) + r(i, i) * mu{i}(k, h)) * alpha_val;
    end
end
