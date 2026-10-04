%{
%{
 % @file pfqn_qdlin.m
 % @brief QD-LIN: the Linearizer arm of AMVA-LD, on a plain demand matrix.
%}
%}

%{
%{
 % @brief QD-LIN, the array-level twin of what SolverMVA computes for
 %        method='qdlin'.
 % @details
 % The Linearizer of Chandy and Neuse (Commun. ACM 25(2), 1982) run inside the
 % queue-dependent AMVA framework of Casale-Perez-Wang (IFIP PERFORMANCE 2015),
 % so the load-dependent term g_k is evaluated at the CORRECTED arrival-instant
 % queue rather than at the plain one.
 %
 % THIS IS A TRANSCRIPTION OF solver_amvald.m TOGETHER WITH
 % solver_amvald_forward.m, restricted to the domain a demand matrix describes:
 % closed classes only, one chain per class, unit visits, PS queueing stations
 % and one optional delay carrying Z. It is NOT an independent re-derivation,
 % and it is not the Wang-Sevcik QDLIN of the same name in native Python and
 % C++ before 2026-09-04, which was Bard-Schweitzer written out.
 %
 % FOUR PROPERTIES OF THE REFERENCE ARE REPRODUCED DELIBERATELY:
 %
 %   1. THE GAMMA CORRECTION IS CLASS-AGGREGATE, IN SLICE 1. solver_amvald.m
 %      allocates the (K,M,K) per-class Linearizer array for qdlin but writes
 %      gamma(s,k) = sum_r Q_s(k,r)/(Nt-1) - sum_r Q(k,r)/Nt into it with two
 %      subscripts, which MATLAB linear-indexes to (s,k,1); slices 2..K stay
 %      zero while every reader indexes gamma per class. The correction that
 %      reaches the residence time is N_1*gamma(r,k,1) - [r==1]*gamma(r,k,1),
 %      which coincides with the queue-dependent AMVA form (Nt-1)*gamma_agg iff
 %      K == 1. method='lin' takes the per-class form instead.
 %   2. A SINGLE-SERVER STATION STILL CARRIES A SOFTMIN TERM. The multiserver
 %      factor is pfqn_lldfun(1+arrival-instant total, [], nservers), whose
 %      softmin at c = 1 is not exactly 1, so qdlin does not reduce to a
 %      textbook single-server AMVA even when every station has one server.
 %   3. THE WAIT FACTOR IS FLOORED AT wtol. This floor is LINE's options.tol,
 %      a DIFFERENT knob from the convergence tolerance, and SolverMVA never
 %      sets it, so it stays at the lineDefaults 1e-4 while the fixed point
 %      converges to iter_tol 1e-6. MATLAB solver_amvald_forward.m does NOT
 %      carry it; native Python does, and removing it there was tried and
 %      reverted on 2026-09-04 because the unfloored python recursion diverges
 %      where MATLAB and C++ do not. See _kb/06-solver-catalog.md.
 %   4. WHICH UTILIZATION IS REPORTED DEPENDS ON THE MODEL. The analyzer
 %      forwards the iterated Uchain to the deaggregation ONLY under lld, cd or
 %      jd scaling; with none of those the deaggregation recomputes T*S/c from
 %      the NOMINAL demand, and the two differ by the iteration residual.
 %
 % MU AND NSERVERS ARE DIFFERENT MECHANISMS, unlike in pfqn_qdamva, which folds
 % the multiserver curve into mu. Here mu is sn.lldscaling, an interpolated rate
 % multiplier per station, and nservers is the server count feeding the softmin.
 % A c-server station is nservers(k)=c, NOT a mu row of min(1:smax,c); the
 % latter reproduces Queue.setLoadDependence, a different station.
 %
 % @fn pfqn_qdlin(L, N, Z, mu, nservers, tol, maxiter, wtol)
 % @param L (M x R) service demand matrix, queueing stations only.
 % @param N (1 x R) population vector, finite.
 % @param Z (1 x R) think time vector; a delay station carrying it is prepended
 %        to the station list when any entry is positive, exactly as the
 %        equivalent Network would hold one. Empty means no think time.
 % @param mu (M x smax) load-dependent rate multipliers, sn.lldscaling; empty
 %        means none.
 % @param nservers (M x 1) server counts; empty means one server everywhere.
 % @param tol Convergence tolerance on the queue lengths (default 1e-6).
 % @param maxiter Iteration budget (default 1000). The outer sweep and each
 %        inner sweep are capped at sqrt(maxiter) and the total number of
 %        forward evaluations at min(maxiter,10000), as in solver_amvald.
 % @param wtol Floor on the wait factor (default 1e-4), LINE's options.tol.
 % @return Q (M x R) mean queue lengths at the queueing stations.
 % @return U (M x R) per-class utilizations.
 % @return R (M x R) per-class residence times.
 % @return X (1 x R) per-class throughputs.
 % @return C (1 x R) per-class cycle times, think time included.
 % @return iter Number of forward evaluations performed.
%}
%}
function [Q,U,R,X,C,iter] = pfqn_qdlin(L,N,Z,mu,nservers,tol,maxiter,wtol)
% [Q,U,R,X,C,ITER] = PFQN_QDLIN(L,N,Z,MU,NSERVERS,TOL,MAXITER,WTOL)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[M,K] = size(L);
N = N(:)';
if nargin < 3 || isempty(Z)
    Z = zeros(1,K);
end
Z = Z(:)';
if nargin < 4
    mu = [];
end
if nargin < 5 || isempty(nservers)
    nservers = ones(M,1);
end
nservers = nservers(:);
if nargin < 6 || isempty(tol)
    tol = 1e-6;
end
if nargin < 7 || isempty(maxiter)
    maxiter = 1000;
end
if nargin < 8 || isempty(wtol)
    wtol = 1e-4;
end
if length(N) ~= K
    line_error(mfilename,'the population vector must have one entry per class');
end
if any(isinf(N))
    line_error(mfilename,'an infinite population is not supported, closed classes only');
end
if length(Z) ~= K
    line_error(mfilename,'the think-time vector must have one entry per class');
end
if length(nservers) ~= M
    line_error(mfilename,'the server-count vector must have one entry per station');
end

Nt = sum(N);
if Nt <= 0
    Q = zeros(M,K); U = zeros(M,K); R = zeros(M,K);
    X = zeros(1,K); C = zeros(1,K); iter = 0;
    return
end

% station list: the delay, when there is one, then the queueing stations
hasDelay = any(Z > 0);
if hasDelay
    ST = [Z; L];
    srv = [Inf; nservers];
    isdelay = [true; false(M,1)];
    if isempty(mu)
        muFull = [];
    else
        muFull = [ones(1,size(mu,2)); mu];
    end
else
    ST = L;
    srv = nservers;
    isdelay = false(M,1);
    muFull = mu;
end
Ms = size(ST,1);

% balanced initialization, as in solver_amvald
Q = ones(Ms,K);
Q = Q ./ repmat(sum(Q,1),Ms,1) .* repmat(N,Ms,1);
Q(:,N==0) = 0;
X = 1 ./ sum(ST,1);
X(~isfinite(X)) = 0;
nnzclasses = find(N > 0);
U = zeros(Ms,K);
for k = 1:Ms
    for r = nnzclasses
        if isinf(srv(k))
            U(k,r) = ST(k,r) * X(r);
        else
            U(k,r) = ST(k,r) * X(r) / srv(k);
        end
    end
end

omicron = 0.5; % under-relaxation parameter of solver_amvald
gamma = zeros(K,Ms,K);
T = zeros(Ms,K);
C = zeros(1,K);
STeff = zeros(Ms,K);

maxSweep = sqrt(maxiter);
maxTotiter = min(maxiter,10000);
iter = 0;
outerIter = 0;
QouterPrev = Q + Inf;

while (outerIter < 2 || max(max(abs(Q-QouterPrev))) > tol) && outerIter < maxSweep && iter <= maxTotiter
    outerIter = outerIter + 1;
    QouterPrev = Q;

    % Linearizer recursion: one sweep at each reduced population N-1_s
    exhausted = false;
    for s = 1:K
        if N(s) <= 0
            continue
        end
        Ns = N;
        Ns(s) = Ns(s) - 1;
        shrink = (Nt-1)/Nt;
        Qs = Q * shrink;
        Xs = X * shrink;

        iterS = 0;
        QsPrev = Qs + Inf;
        while (iterS < 2 || max(max(abs(Qs-QsPrev))) > tol) && iterS <= maxSweep
            iterS = iterS + 1;
            QsPrev = Qs;
            XsPrev = Xs;

            Ws = qdlin_forward(ST, srv, isdelay, muFull, gamma, QsPrev, Ns, K, wtol);
            iter = iter + 1;
            if iter >= maxTotiter
                exhausted = true;
                break
            end

            for r = nnzclasses
                if sum(Ws(:,r)) == 0 || Ns(r) == 0
                    Xs(r) = 0;
                else
                    Cs = sum(Ws(:,r));
                    if Cs > 1e-14
                        Xs(r) = omicron * Ns(r) / Cs + (1-omicron) * XsPrev(r);
                    else
                        Xs(r) = XsPrev(r);
                    end
                end
                for k = 1:Ms
                    Qs(k,r) = omicron * Xs(r) * Ws(k,r) + (1-omicron) * QsPrev(k,r);
                end
            end
        end

        % class-aggregate correction into slice 1, see the header
        if Nt > 1
            for k = 1:Ms
                gamma(s,k,1) = sum(QsPrev(k,:))/(Nt-1) - sum(QouterPrev(k,:))/Nt;
            end
        else
            gamma(s,:,1) = 0;
        end

        if exhausted
            break
        end
    end
    if exhausted
        break
    end

    % sweep at the full population N
    innerIter = 0;
    Qprev = Q + Inf;
    while (innerIter < 2 || max(max(abs(Q-Qprev))) > tol) && innerIter <= maxSweep
        innerIter = innerIter + 1;
        Qprev = Q;
        Xprev = X;
        Uprev = U;

        [W, STeff] = qdlin_forward(ST, srv, isdelay, muFull, gamma, Qprev, N, K, wtol);
        iter = iter + 1;
        if iter >= maxTotiter
            exhausted = true;
            break
        end

        for r = nnzclasses
            if sum(W(:,r)) == 0
                X(r) = 0;
            elseif N(r) == 0
                X(r) = 0;
                C(r) = 0;
            else
                C(r) = sum(W(:,r));
                if C(r) > 1e-14
                    X(r) = omicron * N(r) / C(r) + (1-omicron) * Xprev(r);
                else
                    X(r) = Xprev(r);
                end
            end
            for k = 1:Ms
                Q(k,r) = omicron * X(r) * W(k,r) + (1-omicron) * Qprev(k,r);
                T(k,r) = X(r);
                U(k,r) = omicron * STeff(k,r) * X(r) + (1-omicron) * Uprev(k,r);
            end
        end
    end
    if exhausted
        break
    end
end

% utilization capping, as in solver_amvald: a queueing station whose class
% utilizations sum above one has them renormalized in proportion to STeff.
for k = 1:Ms
    if isdelay(k)
        continue
    end
    Usum = sum(U(k,:));
    if Usum > 1
        denom = sum(STeff(k,:) .* X);
        if denom > 0
            for r = 1:K
                if STeff(k,r) > 0
                    U(k,r) = min(1,Usum) * STeff(k,r) * X(r) / denom;
                end
            end
        end
    end
end

% the analyzer forwards the iterated Uchain to the deaggregation ONLY under lld,
% cd or jd scaling; with none of those it recomputes T*S/c from the NOMINAL
% demand. mu is the only one of the three a demand matrix can carry.
if isempty(mu)
    for k = 1:Ms
        for r = nnzclasses
            if isinf(srv(k))
                U(k,r) = ST(k,r) * X(r);
            else
                U(k,r) = ST(k,r) * X(r) / srv(k);
            end
        end
    end
end

% a class with no jobs keeps its 1/sum(ST) SEED in X unless it is cleared: the
% sweeps only ever write the classes in nnzclasses, so the initial value would
% otherwise be reported as that class's throughput. Q, U and C are already zero
% there because they are written in the same loops.
X(setdiff(1:K,nnzclasses)) = 0;

R = zeros(Ms,K);
nz = T > 0;
R(nz) = Q(nz) ./ T(nz);

keep = ~isdelay;
Q = Q(keep,:);
U = U(keep,:);
R = R(keep,:);
end

function [W, STeff] = qdlin_forward(ST, srv, isdelay, mu, gamma, Qin, Nin, K, wtol)
% One forward evaluation, solver_amvald_forward restricted to PS and INF.
Ms = size(ST,1);
nnz = find(Nin > 0);
Ntin = sum(Nin);
if Ntin > 0
    delta = (Ntin-1)/Ntin;
else
    delta = 1;
end
dcl = ones(1,K);
for r = nnz
    dcl(r) = (Nin(r)-1)/Nin(r);
end

% arrival-instant queue lengths, class-aggregate and per class
interp = zeros(Ms,1);
totArvl = zeros(Ms,K);
for k = 1:Ms
    sumQk = sum(Qin(k,nnz));
    interp(k) = delta * sumQk;
    for r = nnz
        totArvl(k,r) = dcl(r) * Qin(k,r) + sumQk - Qin(k,r);
    end
end

% lld term, evaluated at the gamma-corrected arrival-instant queue
lldterm = ones(Ms,K);
if ~isempty(nnz)
    for r = nnz
        gcorr = reshape(gamma(r,:,nnz), Ms, numel(nnz)) * Nin(nnz)' - reshape(gamma(r,:,r), Ms, 1);
        lldterm(:,r) = pfqn_lldfun(1 + interp + gcorr, mu);
    end
else
    lldterm = repmat(pfqn_lldfun(1 + interp, mu), 1, K);
end

% multiserver term; the 'default' rule leaves PS on the softmin arm
if ~isempty(nnz) && Ntin > 0
    g = zeros(numel(nnz), Ms);
    for r = nnz
        g = g + ((Ntin-1)/Ntin) * Nin(r) * reshape(gamma(nnz,:,r), numel(nnz), Ms);
    end
    msterm = pfqn_lldfun(1 + interp + mean(g,1)', [], srv);
else
    msterm = pfqn_lldfun(1 + interp, [], srv);
end

STeff = zeros(Ms,K);
for r = nnz
    STeff(:,r) = ST(:,r) .* lldterm(:,r) .* msterm;
end

W = zeros(Ms,K);
for r = nnz
    for k = 1:Ms
        if isdelay(k)
            W(k,r) = STeff(k,r);
        else
            corr = Nin(nnz) * reshape(gamma(r,k,nnz), numel(nnz), 1) - gamma(r,k,r);
            W(k,r) = STeff(k,r) * max(wtol, 1 + totArvl(k,r) + corr);
        end
    end
end
end
