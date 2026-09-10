function [result, x, fval, exitflag] = mapqn_bnd_qr(params, objective_queue, objective_phase, sense)
% MAPQN_BND_QR - General Quadratic Reduction Bounds for MAP Queueing Networks
%
% MATLAB port of the Python file bnd_qr.py
%
% Usage:
%   [result, x, fval, exitflag] = mapqn_bnd_qr(params)
%   [result, x, fval, exitflag] = mapqn_bnd_qr(params, objective_queue)
%   [result, x, fval, exitflag] = mapqn_bnd_qr(params, objective_queue, objective_phase)
%   [result, x, fval, exitflag] = mapqn_bnd_qr(params, objective_queue, objective_phase, sense)
%
% Inputs:
%   params          - Structure with model parameters:
%                     .M       - Number of queues
%                     .N       - Total population
%                     .K       - [M x 1] Number of phases for each queue
%                     .mu      - {M x 1} cell, each mu{i} is K(i) x K(i) completion rates
%                     .v       - {M x 1} cell, each v{i} is K(i) x K(i) background rates
%                     .r       - [M x M] Routing probabilities
%                     .verbose - (optional) boolean, default true
%
%   objective_queue - (optional) Queue index to optimize (1-based), default 1
%   objective_phase - (optional) Phase index to optimize (1-based), default 1
%   sense           - (optional) 'min' or 'max', default 'max'
%
% Outputs:
%   result   - Structure with results:
%              .objective - Objective function value
%              .exitflag  - Solver exit flag
%              .U         - [M x max(K)] Utilization matrix
%              .IT        - [M x max(K)] Idle time matrix
%              .Q         - [M x max(K)] Queue length matrix
%              .getP2     - Function handle to extract p2 values
%   x        - Raw solution vector
%   fval     - Objective function value
%   exitflag - Solver exit flag

    if nargin < 2 || isempty(objective_queue)
        objective_queue = 1;
    end
    if nargin < 3 || isempty(objective_phase)
        objective_phase = 1;
    end
    if nargin < 4 || isempty(sense)
        sense = 'max';
    end

    % Extract parameters
    M = params.M;
    N = params.N;
    K = params.K(:);
    mu = params.mu;
    v = params.v;
    r = params.r;
    if isfield(params, 'verbose')
        verbose = params.verbose;
    else
        verbose = true;
    end

    maxK = max(K);

    % Compute transition rates q{i,j}(k,h)
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
    if verbose; fprintf('Building variable index map...\n'); end
    varCount = 0;

    % U(i,k) variables: utilization at queue i, phase k
    Uidx = zeros(M, maxK);
    for i = 1:M
        for k = 1:K(i)
            varCount = varCount + 1;
            Uidx(i, k) = varCount;
        end
    end

    % IT(i,k) variables: idle time at queue i, phase k
    ITidx = zeros(M, maxK);
    for i = 1:M
        for k = 1:K(i)
            varCount = varCount + 1;
            ITidx(i, k) = varCount;
        end
    end

    % UP(j,k,i,h) variables: utilization products
    UPidx = zeros(M, maxK, M, maxK);
    for j = 1:M
        for kj = 1:K(j)
            for i = 1:M
                for hi = 1:K(i)
                    varCount = varCount + 1;
                    UPidx(j, kj, i, hi) = varCount;
                end
            end
        end
    end

    % QP(j,k,i,h) variables: queue-length products
    QPidx = zeros(M, maxK, M, maxK);
    for j = 1:M
        for kj = 1:K(j)
            for i = 1:M
                for hi = 1:K(i)
                    varCount = varCount + 1;
                    QPidx(j, kj, i, hi) = varCount;
                end
            end
        end
    end

    % Q(i,k) variables: mean queue length
    Qidx = zeros(M, maxK);
    for i = 1:M
        for k = 1:K(i)
            varCount = varCount + 1;
            Qidx(i, k) = varCount;
        end
    end

    % C(j,k,i) variables: conditional queue lengths
    Cidx = zeros(M, maxK, M);
    for j = 1:M
        for kj = 1:K(j)
            for i = 1:M
                varCount = varCount + 1;
                Cidx(j, kj, i) = varCount;
            end
        end
    end

    % I_var(j,k,i) variables: conditional idle lengths
    Iidx = zeros(M, maxK, M);
    for j = 1:M
        for kj = 1:K(j)
            for i = 1:M
                varCount = varCount + 1;
                Iidx(j, kj, i) = varCount;
            end
        end
    end

    % p1(j,k,i,ni,h) variables: marginal probabilities, ni from 0 to N
    % Stored as p1idx(j, k, i, ni+1, h)
    p1idx = zeros(M, maxK, M, N+1, maxK);
    for j = 1:M
        for kj = 1:K(j)
            for i = 1:M
                for ni = 0:N
                    for hi = 1:K(i)
                        varCount = varCount + 1;
                        p1idx(j, kj, i, ni+1, hi) = varCount;
                    end
                end
            end
        end
    end

    % p1c(j,k,i,ni,h) variables: complementary marginal probabilities, ni from 0 to N
    % Stored as p1cidx(j, k, i, ni+1, h)
    p1cidx = zeros(M, maxK, M, N+1, maxK);
    for j = 1:M
        for kj = 1:K(j)
            for i = 1:M
                for ni = 0:N
                    for hi = 1:K(i)
                        varCount = varCount + 1;
                        p1cidx(j, kj, i, ni+1, hi) = varCount;
                    end
                end
            end
        end
    end

    % p2(j,nj,k,i,ni,h) variables: joint probabilities, nj from 0 to N, ni from 0 to N
    % Stored as p2idx(j, nj+1, k, i, ni+1, h)
    p2idx = zeros(M, N+1, maxK, M, N+1, maxK);
    for j = 1:M
        for nj = 0:N
            for kj = 1:K(j)
                for i = 1:M
                    for ni = 0:N
                        for hi = 1:K(i)
                            varCount = varCount + 1;
                            p2idx(j, nj+1, kj, i, ni+1, hi) = varCount;
                        end
                    end
                end
            end
        end
    end

    nVars = varCount;
    if verbose; fprintf('Total variables: %d\n', nVars); end

    %% Initialize bounds
    lb = zeros(nVars, 1);
    ub = inf(nVars, 1);

    % Set upper bounds for each variable type
    for i = 1:M
        for k = 1:K(i)
            ub(Uidx(i, k)) = 1;
            ub(ITidx(i, k)) = 1;
            ub(Qidx(i, k)) = N;
        end
    end
    for j = 1:M
        for kj = 1:K(j)
            for i = 1:M
                for hi = 1:K(i)
                    ub(UPidx(j, kj, i, hi)) = 1;
                    ub(QPidx(j, kj, i, hi)) = N;
                end
            end
        end
    end
    for j = 1:M
        for kj = 1:K(j)
            for i = 1:M
                ub(Cidx(j, kj, i)) = N;
                ub(Iidx(j, kj, i)) = N;
            end
        end
    end
    for j = 1:M
        for kj = 1:K(j)
            for i = 1:M
                for ni = 0:N
                    for hi = 1:K(i)
                        ub(p1idx(j, kj, i, ni+1, hi)) = 1;
                        ub(p1cidx(j, kj, i, ni+1, hi)) = 1;
                    end
                end
            end
        end
    end
    for j = 1:M
        for nj = 0:N
            for kj = 1:K(j)
                for i = 1:M
                    for ni = 0:N
                        for hi = 1:K(i)
                            ub(p2idx(j, nj+1, kj, i, ni+1, hi)) = 1;
                        end
                    end
                end
            end
        end
    end

    %% Build constraints
    Aeq = [];
    beq = [];
    Aineq = [];
    bineq = [];

    if verbose; fprintf('Building constraints...\n'); end

    %% ZER1: p1(j,k,j,0,k) = 0 for all j,k (via upper bounds)
    if verbose; fprintf('  ZER1 constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            idx = p1idx(j, k, j, 0+1, k);
            ub(idx) = 0;
        end
    end

    %% ZER2: p1(j,k,j,nj,h) = 0 for h ~= k, all j,k,nj (via upper bounds)
    if verbose; fprintf('  ZER2 constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            for nj = 0:N
                for h = 1:K(j)
                    if h ~= k
                        idx = p1idx(j, k, j, nj+1, h);
                        ub(idx) = 0;
                    end
                end
            end
        end
    end

    %% ZER3: p1(j,k,i,N,h) = 0 for j ~= i, all j,k,i,h (via upper bounds)
    if verbose; fprintf('  ZER3 constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            for i = 1:M
                if j ~= i
                    for h = 1:K(i)
                        idx = p1idx(j, k, i, N+1, h);
                        ub(idx) = 0;
                    end
                end
            end
        end
    end

    %% ZER4: p1c(j,k,j,nj,h) = 0 for nj >= 1, all j,k,nj,h (via upper bounds)
    if verbose; fprintf('  ZER4 constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            for nj = 1:N
                for h = 1:K(j)
                    idx = p1cidx(j, k, j, nj+1, h);
                    ub(idx) = 0;
                end
            end
        end
    end

    %% ZER5: p2(j,nj,k,j,nj,h) = 0 for h ~= k (AMPL ZERO1) (via upper bounds)
    % The joint variables carry the same structural zeros as their p1
    % projections: a station cannot be in two phases at once, cannot hold two
    % different populations at once, and two stations cannot jointly hold more
    % than the closed population.
    if verbose; fprintf('  ZER5 (p2 phase consistency) constraints...\n'); end
    for j = 1:M
        for nj = 0:N
            for k = 1:K(j)
                for h = 1:K(j)
                    if h ~= k
                        ub(p2idx(j, nj+1, k, j, nj+1, h)) = 0;
                    end
                end
            end
        end
    end

    %% ZER6: p2(j,nj,k,j,ni,h) = 0 for ni ~= nj (AMPL ZERO2) (via upper bounds)
    if verbose; fprintf('  ZER6 (p2 population consistency) constraints...\n'); end
    for j = 1:M
        for nj = 0:N
            for k = 1:K(j)
                for ni = 0:N
                    if ni ~= nj
                        for h = 1:K(j)
                            ub(p2idx(j, nj+1, k, j, ni+1, h)) = 0;
                        end
                    end
                end
            end
        end
    end

    %% ZER7: p2(j,nj,k,i,ni,h) = 0 for i ~= j and nj+ni > N (AMPL ZERO3)
    if verbose; fprintf('  ZER7 (p2 population cap) constraints...\n'); end
    for j = 1:M
        for nj = 0:N
            for k = 1:K(j)
                for i = 1:M
                    if i ~= j
                        for ni = 0:N
                            if nj + ni > N
                                for h = 1:K(i)
                                    ub(p2idx(j, nj+1, k, i, ni+1, h)) = 0;
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    %% CEQU: C(j,k,j) = Q(j,k) for all j,k
    if verbose; fprintf('  CEQU constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            row = zeros(1, nVars);
            row(Cidx(j, k, j)) = 1;
            row(Qidx(j, k)) = -1;
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% ONE1: sum over kj,hi,ni of (p1 + p1c) = 1 for each (j,i)
    if verbose; fprintf('  ONE1 constraints...\n'); end
    for j = 1:M
        for i = 1:M
            row = zeros(1, nVars);
            for kj = 1:K(j)
                for hi = 1:K(i)
                    for ni = 0:N
                        row(p1idx(j, kj, i, ni+1, hi)) = 1;
                        row(p1cidx(j, kj, i, ni+1, hi)) = 1;
                    end
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 1];
        end
    end

    %% UTLB: U(i,k) = sum over t,nt,h of p1(i,k,t,nt,h) for each (i,k,t)
    if verbose; fprintf('  UTLB constraints...\n'); end
    for i = 1:M
        for k = 1:K(i)
            for t = 1:M
                row = zeros(1, nVars);
                row(Uidx(i, k)) = 1;
                for nt = 0:N
                    for h = 1:K(t)
                        row(p1idx(i, k, t, nt+1, h)) = -1;
                    end
                end
                Aeq = [Aeq; row];
                beq = [beq; 0];
            end
        end
    end

    %% UTLC: IT(i,k) = sum over t,nt,h of p1c(i,k,t,nt,h) for each (i,k,t)
    if verbose; fprintf('  UTLC constraints...\n'); end
    for i = 1:M
        for k = 1:K(i)
            for t = 1:M
                row = zeros(1, nVars);
                row(ITidx(i, k)) = 1;
                for nt = 0:N
                    for h = 1:K(t)
                        row(p1cidx(i, k, t, nt+1, h)) = -1;
                    end
                end
                Aeq = [Aeq; row];
                beq = [beq; 0];
            end
        end
    end

    %% QLEN: Q(i,k) = sum over ni of ni*p1(i,k,i,ni,k) for each (i,k)
    if verbose; fprintf('  QLEN constraints...\n'); end
    for i = 1:M
        for k = 1:K(i)
            row = zeros(1, nVars);
            row(Qidx(i, k)) = 1;
            for ni = 0:N
                row(p1idx(i, k, i, ni+1, k)) = row(p1idx(i, k, i, ni+1, k)) - ni;
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% SRVB: Service balance
    % sum{j,h} q{i,j}(k,h)*U(i,k) = sum{j,h} q{i,j}(h,k)*U(i,h) for each (i,k)
    % Stations with a single phase give an identically zero row, so skip them.
    if verbose; fprintf('  SRVB constraints...\n'); end
    for i = 1:M
        if K(i) < 2
            continue;
        end
        for k = 1:K(i)
            row = zeros(1, nVars);
            for j = 1:M
                for h = 1:K(i)
                    row(Uidx(i, k)) = row(Uidx(i, k)) + q{i,j}(k, h);
                    row(Uidx(i, h)) = row(Uidx(i, h)) - q{i,j}(h, k);
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% POPC: sum over i,k of Q(i,k) = N
    if verbose; fprintf('  POPC constraint...\n'); end
    row = zeros(1, nVars);
    for i = 1:M
        for k = 1:K(i)
            row(Qidx(i, k)) = 1;
        end
    end
    Aeq = [Aeq; row];
    beq = [beq; N];

    %% ONE: sum over k of (U(j,k) + IT(j,k)) = 1 for each j
    if verbose; fprintf('  ONE constraints...\n'); end
    for j = 1:M
        row = zeros(1, nVars);
        for k = 1:K(j)
            row(Uidx(j, k)) = 1;
            row(ITidx(j, k)) = 1;
        end
        Aeq = [Aeq; row];
        beq = [beq; 1];
    end

    %% PCL2: sum over i,j,ni>=1,nj>=1,h,k of ni*nj*p2(i,ni,h,j,nj,k) = N^2
    if verbose; fprintf('  PCL2 constraint...\n'); end
    row = zeros(1, nVars);
    for i = 1:M
        for j = 1:M
            for ni = 1:N
                for nj = 1:N
                    for h = 1:K(i)
                        for k = 1:K(j)
                            idx = p2idx(i, ni+1, h, j, nj+1, k);
                            row(idx) = row(idx) + ni * nj;
                        end
                    end
                end
            end
        end
    end
    Aeq = [Aeq; row];
    beq = [beq; N^2];

    %% PI21: p1(j,k,i,ni,h) = sum over nj=1..N of p2(j,nj,k,i,ni,h)
    if verbose; fprintf('  PI21 constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            for i = 1:M
                for ni = 0:N
                    for h = 1:K(i)
                        row = zeros(1, nVars);
                        row(p1idx(j, k, i, ni+1, h)) = 1;
                        for nj = 1:N
                            idx = p2idx(j, nj+1, k, i, ni+1, h);
                            row(idx) = row(idx) - 1;
                        end
                        Aeq = [Aeq; row];
                        beq = [beq; 0];
                    end
                end
            end
        end
    end

    %% PI22: p1c(j,k,i,ni,h) = p2(j,0,k,i,ni,h)
    if verbose; fprintf('  PI22 constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            for i = 1:M
                for ni = 0:N
                    for h = 1:K(i)
                        row = zeros(1, nVars);
                        row(p1cidx(j, k, i, ni+1, h)) = 1;
                        row(p2idx(j, 0+1, k, i, ni+1, h)) = -1;
                        Aeq = [Aeq; row];
                        beq = [beq; 0];
                    end
                end
            end
        end
    end

    %% PI23: p2(i,ni,h,j,nj,k) = p2(j,nj,k,i,ni,h) (symmetry)
    % Only generate for j < i, or (j == i and nj < ni), to avoid redundancy
    if verbose; fprintf('  PI23 (symmetry) constraints...\n'); end
    for j = 1:M
        for nj = 0:N
            for k = 1:K(j)
                for i = 1:M
                    for ni = 0:N
                        for h = 1:K(i)
                            % Only generate constraint if (j,nj,k) < (i,ni,h) in lex order
                            if j < i || (j == i && nj < ni) || (j == i && nj == ni && k < h)
                                idx1 = p2idx(i, ni+1, h, j, nj+1, k);
                                idx2 = p2idx(j, nj+1, k, i, ni+1, h);
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

    %% CLEN: C(j,k,i) = sum over ni,h of ni*p1(j,k,i,ni,h) for each (j,k,i)
    % Without this the C variables are defined only on the diagonal by CEQU,
    % which leaves CUB1/CUB2 vacuous for i ~= j.
    if verbose; fprintf('  CLEN constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            for i = 1:M
                row = zeros(1, nVars);
                row(Cidx(j, k, i)) = 1;
                for ni = 0:N
                    for h = 1:K(i)
                        idx = p1idx(j, k, i, ni+1, h);
                        row(idx) = row(idx) - ni;
                    end
                end
                Aeq = [Aeq; row];
                beq = [beq; 0];
            end
        end
    end

    %% MARG: p2(j,nj,k,j,nj,k) = sum over ni<=N-nj,h of p2(j,nj,k,i,ni,h)
    % AMPL MARGINALS. ONE1 only imposes the aggregate over nj, so the
    % per-population marginal consistency is a separate family.
    if verbose; fprintf('  MARG constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            for nj = 0:N
                for i = 1:M
                    if i ~= j
                        row = zeros(1, nVars);
                        row(p2idx(j, nj+1, k, j, nj+1, k)) = 1;
                        for ni = 0:(N-nj)
                            for h = 1:K(i)
                                idx = p2idx(j, nj+1, k, i, ni+1, h);
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

    %% THM2: sum over i,ni>=1,h of ni*p2(j,nj,k,i,ni,h) = N*p2(j,nj,k,j,nj,k)
    % AMPL THM2, the queue-length theorem conditioned on (j,nj,k). Summing it
    % over nj >= 1 recovers the aggregate form sum_i C(j,k,i) = N*U(j,k).
    if verbose; fprintf('  THM2 constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            for nj = 0:N
                row = zeros(1, nVars);
                for i = 1:M
                    for ni = 1:N
                        for h = 1:K(i)
                            idx = p2idx(j, nj+1, k, i, ni+1, h);
                            row(idx) = row(idx) + ni;
                        end
                    end
                end
                idx = p2idx(j, nj+1, k, j, nj+1, k);
                row(idx) = row(idx) - N;
                Aeq = [Aeq; row];
                beq = [beq; 0];
            end
        end
    end

    %% THM30: level-crossing balance at an empty station, per arrival phase
    % AMPL THM30. Rate into {n_i = 0, phase_i = u} from a busy neighbour equals
    % the rate out of {n_i = 1} through a completion at i.
    if verbose; fprintf('  THM30 constraints...\n'); end
    for i = 1:M
        for u = 1:K(i)
            row = zeros(1, nVars);
            for j = 1:M
                if j ~= i
                    for nj = 1:N
                        for k = 1:K(j)
                            for h = 1:K(j)
                                idx = p2idx(j, nj+1, k, i, 0+1, u);
                                row(idx) = row(idx) + q{j,i}(k, h);
                            end
                        end
                    end
                    for nj = 0:N
                        for k = 1:K(i)
                            for h = 1:K(j)
                                idx = p2idx(j, nj+1, h, i, 1+1, k);
                                row(idx) = row(idx) - q{i,j}(k, u);
                            end
                        end
                    end
                end
            end
            Aeq = [Aeq; row];
            beq = [beq; 0];
        end
    end

    %% THM3: level-crossing balance between n_i and n_i+1
    % AMPL THM3. This is the family that ties station i's arrival rate to its
    % departure rate; without it nothing prevents a station from being idle with
    % probability one, which is why the minimum utilization collapsed to zero.
    if verbose; fprintf('  THM3 constraints...\n'); end
    for i = 1:M
        for ni = 0:(N-1)
            row = zeros(1, nVars);
            for j = 1:M
                if j ~= i
                    for nj = 1:N
                        for k = 1:K(j)
                            for h = 1:K(j)
                                for u = 1:K(i)
                                    idx = p2idx(j, nj+1, k, i, ni+1, u);
                                    row(idx) = row(idx) + q{j,i}(k, h);
                                end
                            end
                        end
                    end
                    for nj = 0:N
                        for k = 1:K(i)
                            for u = 1:K(j)
                                for h = 1:K(i)
                                    idx = p2idx(j, nj+1, u, i, ni+1+1, k);
                                    row(idx) = row(idx) - q{i,j}(k, h);
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

    %% UUB1: sum over k of U(i,k) <= 1 for each i (inequality)
    if verbose; fprintf('  UUB1 constraints...\n'); end
    for i = 1:M
        row = zeros(1, nVars);
        for k = 1:K(i)
            row(Uidx(i, k)) = 1;
        end
        Aineq = [Aineq; row];
        bineq = [bineq; 1];
    end

    %% QUB1: Q(j,k) <= N*U(j,k) for each (j,k) (inequality)
    if verbose; fprintf('  QUB1 constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            row = zeros(1, nVars);
            row(Qidx(j, k)) = 1;
            row(Uidx(j, k)) = -N;
            Aineq = [Aineq; row];
            bineq = [bineq; 0];
        end
    end

    %% CUB1: C(j,k,i) <= sum over h of Q(i,h) for each (j,k,i) (inequality)
    if verbose; fprintf('  CUB1 constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            for i = 1:M
                row = zeros(1, nVars);
                row(Cidx(j, k, i)) = 1;
                for h = 1:K(i)
                    row(Qidx(i, h)) = -1;
                end
                Aineq = [Aineq; row];
                bineq = [bineq; 0];
            end
        end
    end

    %% CUB2: C(j,k,i) <= N*U(j,k) for each (j,k,i) (inequality)
    if verbose; fprintf('  CUB2 constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            for i = 1:M
                row = zeros(1, nVars);
                row(Cidx(j, k, i)) = 1;
                row(Uidx(j, k)) = -N;
                Aineq = [Aineq; row];
                bineq = [bineq; 0];
            end
        end
    end

    %% THM4: N*P(j busy in k, i nonempty) <= sum_t C(j,k,t) (inequality)
    % AMPL THM4.
    if verbose; fprintf('  THM4 constraints...\n'); end
    for j = 1:M
        for k = 1:K(j)
            for i = 1:M
                row = zeros(1, nVars);
                for h = 1:K(i)
                    for nj = 0:N
                        for ni = 1:N
                            idx = p2idx(j, nj+1, k, i, ni+1, h);
                            row(idx) = row(idx) + N;
                        end
                    end
                end
                for t = 1:M
                    for h = 1:K(t)
                        for nj = 0:N
                            for nt = 0:N
                                idx = p2idx(j, nj+1, k, t, nt+1, h);
                                row(idx) = row(idx) - nt;
                            end
                        end
                    end
                end
                Aineq = [Aineq; row];
                bineq = [bineq; 0];
            end
        end
    end

    %% Build objective function
    if verbose; fprintf('Building objective function...\n'); end
    c = zeros(nVars, 1);
    c(Uidx(objective_queue, objective_phase)) = 1;

    if strcmp(sense, 'max')
        c = -c;
    end

    %% Solve LP
    if verbose
        fprintf('Solving LP with %d variables and %d equality + %d inequality constraints...\n', ...
            nVars, size(Aeq, 1), size(Aineq, 1));
    end

    % LP algorithm. R2025a's default 'dual-simplex-highs' is broken in some
    % installs (errors "Unrecognized field name optimstatus"), so an
    % interior-point variant is required. 'interior-point-legacy' is more
    % accurate on loosely constrained instances, but once the balance families
    % are present it declares the (feasible) system infeasible with exitflag -2
    % and returns the bound box -- e.g. U1 in [0,1] instead of the exact
    % [0.75,0.75] on the K=1 symmetric tandem. 'interior-point' solves those to
    % ~1e-6 and agrees with the GLPK optimum of the reference AMPL model, so it
    % is the default. Override via params.lpAlgorithm.
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
        warning('mapqn_bnd_qr:nonFiniteSolution', ...
            'linprog returned a non-finite solution (exitflag %d); U and the other metric fields are left unpopulated.', ...
            exitflag);
    end

    if hasSolution
        % Compute utilizations, idle times, queue lengths
        result.U = zeros(M, maxK);
        result.IT = zeros(M, maxK);
        result.Q = zeros(M, maxK);

        for i = 1:M
            for k = 1:K(i)
                result.U(i, k) = x(Uidx(i, k));
                result.IT(i, k) = x(ITidx(i, k));
                result.Q(i, k) = x(Qidx(i, k));
            end
        end

        % Function handle to extract p2 values: p2(j,nj,k,i,ni,h)
        result.getP2 = @(j, nj, k, i, ni, h) x(p2idx(j, nj+1, k, i, ni+1, h));
    end

    if verbose; fprintf('\n=== Results ===\n'); end
    if verbose; fprintf('Objective value: %f\n', fval); end
    if verbose; fprintf('Exit flag: %d\n', exitflag); end
    if hasSolution && verbose
        fprintf('\nUtilizations:\n');
        for i = 1:M
            fprintf('  Queue %d: U = [', i);
            for k = 1:K(i)
                fprintf('%.6f ', result.U(i, k));
            end
            fprintf(']\n');
        end
        fprintf('\nIdle times:\n');
        for i = 1:M
            fprintf('  Queue %d: IT = [', i);
            for k = 1:K(i)
                fprintf('%.6f ', result.IT(i, k));
            end
            fprintf(']\n');
        end
        fprintf('\nQueue lengths:\n');
        for i = 1:M
            fprintf('  Queue %d: Q = [', i);
            for k = 1:K(i)
                fprintf('%.6f ', result.Q(i, k));
            end
            fprintf(']\n');
        end
    end
end
