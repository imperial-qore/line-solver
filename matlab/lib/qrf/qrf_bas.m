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
%               .lpAlgorithm - (optional) linprog Algorithm string;
%                       defaults to 'interior-point-legacy' (see note at the
%                       LP solve for why, and for the R2025a HiGHS caveat)
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
    if isfield(params, 'verbose')
        verbose = params.verbose;
    else
        verbose = true;
    end

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

    if verbose; fprintf('Building variable index map...\n'); end
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
    if verbose; fprintf('Total variables: %d\n', nVars); end

    %% Build constraints
    % Sparse triplet accumulation. The previous builder used
    %   row = zeros(1,nVars); row(idx) = row(idx) + v; Aeq = [Aeq; row];
    % which is O(nrows^2) in memory traffic: every append reallocated and copied
    % the whole matrix, and every row allocated nVars doubles. For the BAS
    % example (nVars ~6e4, ~1.1e5 rows) that does not terminate in practice.
    % We emit (row,col,val) triplets and build one sparse matrix at the end.
    % sparse() SUMS duplicate (i,j) entries, which reproduces the
    % row(idx)=row(idx)+v accumulation exactly. The plain assignments below
    % (SYMMETRY +-1, UEFF row(idx_e)=-1, THM2 row(idx_diag)=-N, ONE/MARGINALS
    % =1) are each the FIRST write to their index within their row, so emitting
    % them as triplets is equivalent; SYMMETRY additionally guards idx1~=idx2.
    eqNnzCap = 1048576;
    eqI = zeros(eqNnzCap, 1); eqJ = zeros(eqNnzCap, 1); eqV = zeros(eqNnzCap, 1);
    eqNnz = 0; nEq = 0;
    inNnzCap = 65536;
    inI = zeros(inNnzCap, 1); inJ = zeros(inNnzCap, 1); inV = zeros(inNnzCap, 1);
    inNnz = 0; nIn = 0;
    beq = zeros(65536, 1);
    bineq = zeros(4096, 1);
    rCap = 4096; rI = zeros(rCap, 1); rV = zeros(rCap, 1); rN = 0;

    % Helper functions
    getP2Idx = @(j, nj, kj, i, ni, hi, m) p2idx{j}{nj+1, kj, i, m}(ni+1, hi);
    getEIdx = @(i, ki) eidx{i}(ki);

    if verbose; fprintf('Building constraints...\n'); end

    %% Initialize bounds
    lb = zeros(nVars, 1);
    ub = inf(nVars, 1);

    %% ZERO constraints - fix infeasible states
    if verbose; fprintf('  ZERO constraints...\n'); end
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
    if verbose; fprintf('  ONE constraints...\n'); end
    for j = 1:M
        rN = 0;
        for nj = 0:N
            for kj = 1:K(j)
                for m = 1:MR
                    idx = getP2Idx(j, nj, kj, j, nj, kj, m);
                    rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=1;
                end
            end
        end
        nEq=nEq+1; if eqNnz+rN>numel(eqI), eqI(2*(eqNnz+rN))=0; eqJ(2*(eqNnz+rN))=0; eqV(2*(eqNnz+rN))=0; end; eqI(eqNnz+(1:rN))=nEq; eqJ(eqNnz+(1:rN))=rI(1:rN); eqV(eqNnz+(1:rN))=rV(1:rN); eqNnz=eqNnz+rN; if nEq>numel(beq), beq(2*nEq)=0; end
        beq(nEq)=1;
    end

    %% SYMMETRY
    if verbose; fprintf('  SYMMETRY constraints...\n'); end
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
                                    rN = 0;
                                    rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx1; rV(rN)=1;
                                    rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx2; rV(rN)=-1;
                                    nEq=nEq+1; if eqNnz+rN>numel(eqI), eqI(2*(eqNnz+rN))=0; eqJ(2*(eqNnz+rN))=0; eqV(2*(eqNnz+rN))=0; end; eqI(eqNnz+(1:rN))=nEq; eqJ(eqNnz+(1:rN))=rI(1:rN); eqV(eqNnz+(1:rN))=rV(1:rN); eqNnz=eqNnz+rN; if nEq>numel(beq), beq(2*nEq)=0; end
                                    beq(nEq)=0;
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    %% MARGINALS
    if verbose; fprintf('  MARGINALS constraints...\n'); end
    for j = 1:M
        for kj = 1:K(j)
            for nj = 0:min(N, F(j))
                for i = 1:M
                    if i == j
                        continue;
                    end
                    for m = 1:MR
                        rN = 0;
                        idx_diag = getP2Idx(j, nj, kj, j, nj, kj, m);
                        rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx_diag; rV(rN)=1;
                        for ni = 0:min(N-nj, F(i))
                            for hi = 1:K(i)
                                idx = getP2Idx(j, nj, kj, i, ni, hi, m);
                                rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=-(1);
                            end
                        end
                        nEq=nEq+1; if eqNnz+rN>numel(eqI), eqI(2*(eqNnz+rN))=0; eqJ(2*(eqNnz+rN))=0; eqV(2*(eqNnz+rN))=0; end; eqI(eqNnz+(1:rN))=nEq; eqJ(eqNnz+(1:rN))=rI(1:rN); eqV(eqNnz+(1:rN))=rV(1:rN); eqNnz=eqNnz+rN; if nEq>numel(beq), beq(2*nEq)=0; end
                        beq(nEq)=0;
                    end
                end
            end
        end
    end

    %% UEFF: e(i,ki) = sum of p2 where queue i is not blocked
    if verbose; fprintf('  UEFF constraints...\n'); end
    for i = 1:M
        for ki = 1:K(i)
            rN = 0;
            idx_e = getEIdx(i, ki);
            rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx_e; rV(rN)=-1;
            for j = 1:M
                for nj = 0:min(N, F(j))
                    for kj = 1:K(j)
                        for m = 1:MR
                            if BB(m, i) == 0
                                for ni = 1:min(N, F(i))
                                    idx = getP2Idx(j, nj, kj, i, ni, ki, m);
                                    rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=1;
                                end
                            end
                        end
                    end
                end
            end
            nEq=nEq+1; if eqNnz+rN>numel(eqI), eqI(2*(eqNnz+rN))=0; eqJ(2*(eqNnz+rN))=0; eqV(2*(eqNnz+rN))=0; end; eqI(eqNnz+(1:rN))=nEq; eqJ(eqNnz+(1:rN))=rI(1:rN); eqV(eqNnz+(1:rN))=rV(1:rN); eqNnz=eqNnz+rN; if nEq>numel(beq), beq(2*nEq)=0; end
            beq(nEq)=0;
        end
    end

    %% THM1: Phase balance (Theorem 1)
    % sum {j, h: j<>i or h<>k} q(i,j,k,h)*e(i,k) = sum {j, h: j<>i or h<>k} q(i,j,h,k)*e(i,h)
    if verbose; fprintf('  THM1 (Phase balance) constraints...\n'); end
    for i = 1:M
        for ki = 1:K(i)
            rN = 0;
            % LHS
            for j = 1:M
                for hi = 1:K(i)
                    if j ~= i || hi ~= ki
                        idx_e = getEIdx(i, ki);
                        rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx_e; rV(rN)=q{i,j}(ki, hi);
                    end
                end
            end
            % RHS (subtract)
            for j = 1:M
                for hi = 1:K(i)
                    if j ~= i || hi ~= ki
                        idx_e = getEIdx(i, hi);
                        rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx_e; rV(rN)=-(q{i,j}(hi, ki));
                    end
                end
            end
            nEq=nEq+1; if eqNnz+rN>numel(eqI), eqI(2*(eqNnz+rN))=0; eqJ(2*(eqNnz+rN))=0; eqV(2*(eqNnz+rN))=0; end; eqI(eqNnz+(1:rN))=nEq; eqJ(eqNnz+(1:rN))=rI(1:rN); eqV(eqNnz+(1:rN))=rV(1:rN); eqNnz=eqNnz+rN; if nEq>numel(beq), beq(2*nEq)=0; end
            beq(nEq)=0;
        end
    end

    %% THM2: Population constraint (Theorem 2)
    if verbose; fprintf('  THM2 (Population) constraints...\n'); end
    for j = 1:M
        for kj = 1:K(j)
            for nj = 0:F(j)
                for m = 1:MR
                    rN = 0;
                    % RHS: -N * p2(j,nj,kj,j,nj,kj,m)
                    idx_diag = getP2Idx(j, nj, kj, j, nj, kj, m);
                    rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx_diag; rV(rN)=-N;
                    % LHS: sum
                    for i = 1:M
                        for ni = 1:F(i)
                            for ki = 1:K(i)
                                idx = getP2Idx(j, nj, kj, i, ni, ki, m);
                                rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=ni;
                            end
                        end
                    end
                    nEq=nEq+1; if eqNnz+rN>numel(eqI), eqI(2*(eqNnz+rN))=0; eqJ(2*(eqNnz+rN))=0; eqV(2*(eqNnz+rN))=0; end; eqI(eqNnz+(1:rN))=nEq; eqJ(eqNnz+(1:rN))=rI(1:rN); eqV(eqNnz+(1:rN))=rV(1:rN); eqNnz=eqNnz+rN; if nEq>numel(beq), beq(2*nEq)=0; end
                    beq(nEq)=0;
                end
            end
        end
    end

    %% COR1: Second moment constraint (Corollary to Theorem 2)
    % sum_{m,i,j,nj,ni,ki,kj} ni*nj*p2(j,nj,kj,i,ni,ki,m) = N^2
    if verbose; fprintf('  COR1 (Second moment) constraint...\n'); end
    rN = 0;
    for m = 1:MR
        for i = 1:M
            for j = 1:M
                for nj = 1:F(j)
                    for ni = 1:F(i)
                        for ki = 1:K(i)
                            for kj = 1:K(j)
                                idx = getP2Idx(j, nj, kj, i, ni, ki, m);
                                rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=ni * nj;
                            end
                        end
                    end
                end
            end
        end
    end
    nEq=nEq+1; if eqNnz+rN>numel(eqI), eqI(2*(eqNnz+rN))=0; eqJ(2*(eqNnz+rN))=0; eqV(2*(eqNnz+rN))=0; end; eqI(eqNnz+(1:rN))=nEq; eqJ(eqNnz+(1:rN))=rI(1:rN); eqV(eqNnz+(1:rN))=rV(1:rN); eqNnz=eqNnz+rN; if nEq>numel(beq), beq(2*nEq)=0; end
    beq(nEq)=N^2;

    %% THM30: Marginal balance for ni=0 (per phase), i<>f
    if verbose; fprintf('  THM30 (Marginal balance ni=0) constraints...\n'); end
    for i = 1:M
        if i == f
            continue;
        end
        for ui = 1:K(i)
            rN = 0;
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
                                    rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=q{j,i}(kj, hj);
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
                                rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=q{f,i}(kj, hj);
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
                                    rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=-(q{i,j}(ki, ui));
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
                                rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=-(q{i,f}(ki, ui));
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
                                    rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=-(q{f,w}(kf, pf));
                                end
                            end
                        end
                    end
                end
            end

            nEq=nEq+1; if eqNnz+rN>numel(eqI), eqI(2*(eqNnz+rN))=0; eqJ(2*(eqNnz+rN))=0; eqV(2*(eqNnz+rN))=0; end; eqI(eqNnz+(1:rN))=nEq; eqJ(eqNnz+(1:rN))=rI(1:rN); eqV(eqNnz+(1:rN))=rV(1:rN); eqNnz=eqNnz+rN; if nEq>numel(beq), beq(2*nEq)=0; end
            beq(nEq)=0;
        end
    end

    %% THM3: Marginal balance for ni in 1:F(i)-1, i<>f
    if verbose; fprintf('  THM3 (Marginal balance) constraints...\n'); end
    for i = 1:M
        if i == f
            continue;
        end
        for ni = 1:(F(i)-1)
            rN = 0;
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
                                        rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=q{j,i}(kj, hj);
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
                                    rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=q{f,i}(kj, hj);
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
                                        rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=-(q{i,j}(ki, hi));
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
                                    rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=-(q{i,f}(ki, hi));
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
                                        rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=-(q{f,w}(kf, pf));
                                    end
                                end
                            end
                        end
                    end
                end
            end

            nEq=nEq+1; if eqNnz+rN>numel(eqI), eqI(2*(eqNnz+rN))=0; eqJ(2*(eqNnz+rN))=0; eqV(2*(eqNnz+rN))=0; end; eqI(eqNnz+(1:rN))=nEq; eqJ(eqNnz+(1:rN))=rI(1:rN); eqV(eqNnz+(1:rN))=rV(1:rN); eqNnz=eqNnz+rN; if nEq>numel(beq), beq(2*nEq)=0; end
            beq(nEq)=0;
        end
    end

    %% THM3f: Marginal balance for i==f
    if verbose; fprintf('  THM3f (Marginal balance for finite queue) constraints...\n'); end
    for ni = 0:(F(f)-1)
        rN = 0;
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
                                        rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=q{j,f}(kj, hj);
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
                                rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=-(q{f,j}(kf, hf));
                            end
                        end
                    end
                end
            end
        end

        nEq=nEq+1; if eqNnz+rN>numel(eqI), eqI(2*(eqNnz+rN))=0; eqJ(2*(eqNnz+rN))=0; eqV(2*(eqNnz+rN))=0; end; eqI(eqNnz+(1:rN))=nEq; eqJ(eqNnz+(1:rN))=rI(1:rN); eqV(eqNnz+(1:rN))=rV(1:rN); eqNnz=eqNnz+rN; if nEq>numel(beq), beq(2*nEq)=0; end
        beq(nEq)=0;
    end

    %% THM3I: Blocking depth balance (Theorem 4)
    if verbose; fprintf('  THM3I (Blocking depth balance) constraints...\n'); end
    for z = 0:(ZM-1)
        rN = 0;
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
                                    rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=q{j,f}(kj, hj);
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
                                    rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=-(q{f,j}(kf, hf));
                                end
                            end
                        end
                    end
                end
            end
        end

        nEq=nEq+1; if eqNnz+rN>numel(eqI), eqI(2*(eqNnz+rN))=0; eqJ(2*(eqNnz+rN))=0; eqV(2*(eqNnz+rN))=0; end; eqI(eqNnz+(1:rN))=nEq; eqJ(eqNnz+(1:rN))=rI(1:rN); eqV(eqNnz+(1:rN))=rV(1:rN); eqNnz=eqNnz+rN; if nEq>numel(beq), beq(2*nEq)=0; end
        beq(nEq)=0;
    end

    %% THM3L: Maximum blocking depth constraint
    if verbose; fprintf('  THM3L (Max blocking depth) constraints...\n'); end
    for m = 1:MR
        if ZZ(m) ~= ZM - 1
            continue;
        end
        rN = 0;
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
                            rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=q{j,f}(kj, hj);
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
                            rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=-(q{f,w}(kf, uf));
                        end
                    end
                end
            end
        end

        nEq=nEq+1; if eqNnz+rN>numel(eqI), eqI(2*(eqNnz+rN))=0; eqJ(2*(eqNnz+rN))=0; eqV(2*(eqNnz+rN))=0; end; eqI(eqNnz+(1:rN))=nEq; eqJ(eqNnz+(1:rN))=rI(1:rN); eqV(eqNnz+(1:rN))=rV(1:rN); eqNnz=eqNnz+rN; if nEq>numel(beq), beq(2*nEq)=0; end
        beq(nEq)=0;
    end

    %% THM4: Queue-length bound inequality (Theorem 5)
    if verbose; fprintf('  THM4 (Queue-length bound) constraints...\n'); end
    for j = 1:M
        for kj = 1:K(j)
            for i = 1:M
                for m = 1:MR
                    rN = 0;
                    % LHS: sum_t sum_ht sum_nj sum_nt nt * p2
                    for t = 1:M
                        for ht = 1:K(t)
                            for nj = 0:F(j)
                                for nt = 1:F(t)
                                    idx = getP2Idx(j, nj, kj, t, nt, ht, m);
                                    rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=nt;
                                end
                            end
                        end
                    end
                    % RHS: -N * sum
                    for hi = 1:K(i)
                        for nj = 0:F(j)
                            for ni = 1:F(i)
                                idx = getP2Idx(j, nj, kj, i, ni, hi, m);
                                rN=rN+1; if rN>rCap, rCap=2*rCap; rI(rCap)=0; rV(rCap)=0; end; rI(rN)=idx; rV(rN)=-(N);
                            end
                        end
                    end
                    nIn=nIn+1; if inNnz+rN>numel(inI), inI(2*(inNnz+rN))=0; inJ(2*(inNnz+rN))=0; inV(2*(inNnz+rN))=0; end; inI(inNnz+(1:rN))=nIn; inJ(inNnz+(1:rN))=rI(1:rN); inV(inNnz+(1:rN))=-rV(1:rN); inNnz=inNnz+rN; if nIn>numel(bineq), bineq(2*nIn)=0; end
                    bineq(nIn)=0;
                end
            end
        end
    end

    % Materialize the sparse constraint matrices from the triplets.
    Aeq = sparse(eqI(1:eqNnz), eqJ(1:eqNnz), eqV(1:eqNnz), nEq, nVars);
    beq = beq(1:nEq);
    Aineq = sparse(inI(1:inNnz), inJ(1:inNnz), inV(1:inNnz), nIn, nVars);
    bineq = bineq(1:nIn);

    %% Build objective function
    if verbose; fprintf('Building objective function...\n'); end
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
    if verbose; fprintf('Solving LP with %d variables and %d equality + %d inequality constraints...\n', ...
        nVars, size(Aeq, 1), size(Aineq, 1)); end

    % LP algorithm. R2025a's default 'dual-simplex-highs' is broken in some
    % installs (errors "Unrecognized field name optimstatus"), so an
    % interior-point variant is required here. On the paper's BAS instance,
    % which is badly scaled (mu spans 1.016186 down to 2.585708e-05),
    % 'interior-point-legacy' is markedly more accurate than 'interior-point':
    % against the published GLPK optimum it gives |dU1min|=8.5e-07 and
    % |dU1max|=9.2e-10, versus 4.3e-05 and 8.0e-06 for 'interior-point'. The
    % residual is pure solver tolerance, not formulation: the minimum lands
    % above and the maximum below the GLPK vertex, and both gaps shrink
    % together as the solver becomes more accurate. Override via
    % params.lpAlgorithm if a particular model needs a different method.
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
        warning('qrf_bas:nonFiniteSolution', ...
            'linprog returned a non-finite solution (exitflag %d); U and the other metric fields are left unpopulated.', ...
            exitflag);
    end

    if hasSolution
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

    if verbose; fprintf('\n=== Results ===\n'); end
    if verbose; fprintf('Objective value: %f\n', fval); end
    if verbose; fprintf('Exit flag: %d\n', exitflag); end
    if hasSolution && verbose
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
