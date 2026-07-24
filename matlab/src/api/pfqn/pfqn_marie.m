%{
%{
 % @file pfqn_marie.m
 % @brief Marie's iterative aggregation for closed networks with FCFS Coxian
 %        (non-exponential) service. Single-class exact-reducing; multiclass via
 %        QD-AMVA with class-dependent (cd) scaling.
%}
%}

function [X,Q,U,C,it,mu] = pfqn_marie(L,N,Z,scv,varargin)
%{
%{
 % @brief Marie's method (Marie 1979/1980): approximate mean performance of a
 %        closed queueing network with FCFS general (Coxian) service, via
 %        iterative aggregation-decomposition. Each station is analyzed in
 %        isolation as a lambda(n)/Cox/1 queue; the resulting conditional
 %        throughputs mu_i(n) drive a load-dependent aggregate solve, iterated
 %        to a fixed point. Single class (R=1): the aggregate is the exact LD
 %        product-form solve pfqn_mvald, and the method reduces to exact product
 %        form for exponential service. Multiple classes (R>1): the aggregate is
 %        QD-AMVA with class-dependent (cd) scaling beta_{i,r}(nvec) supplied by
 %        a multiclass Cox/1 isolation sub-model; exact product-form service
 %        (scv==1 with class-independent means) is dispatched to exact MVA,
 %        otherwise the result is a decomposition approximation.
 % @fn pfqn_marie(L, N, Z, scv, tol, maxiter, nservers)
 % @param L Service demand matrix (M x R).
 % @param N Population vector (1 x R).
 % @param Z Think time vector (1 x R; total delay demand per class).
 % @param scv Per-station per-class squared coefficient of variation (M x R).
 %            scv==1 exponential; scv<0.5 Erlang; scv>0.5 two-phase Coxian.
 % @param tol Convergence tolerance (default 1e-8).
 % @param maxiter Maximum iterations (default 1000).
 % @param nservers Per-station server count (M x 1, default all 1); single-class
 %            only (multiserver multiclass isolation is not yet supported).
 % @return X Throughput: single class M x 1 (per station); multiclass 1 x R
 %            (per-class chain throughput, visits folded into L).
 % @return Q Mean queue length (M x 1 single class, M x R multiclass).
 % @return U Utilization (same shape as Q).
 % @return C Residence time (same shape as Q).
 % @return it Iterations performed.
 % @return mu Converged LD data (M x N single class; cell of cd-scalings R>1).
%}
%}

R = size(L,2);
if R > 1
    [X,Q,U,C,it,mu] = marie_multi(L,N,Z,scv,varargin{:});
    return
end

L = L(:);
M = numel(L);
if nargin < 3 || isempty(Z)
    Z = 0;
end
Z = sum(Z(:));
if nargin < 4 || isempty(scv)
    scv = ones(M,1);
end
scv = scv(:);

tol = 1e-8;
maxiter = 1000;
nservers = ones(M,1);
if numel(varargin) >= 1 && ~isempty(varargin{1}), tol = varargin{1}; end
if numel(varargin) >= 2 && ~isempty(varargin{2}), maxiter = varargin{2}; end
if numel(varargin) >= 3 && ~isempty(varargin{3}), nservers = varargin{3}(:); end
if isscalar(nservers), nservers = nservers*ones(M,1); end

% Per-station Coxian phase representation of the service (mean = L(i), scv(i)).
phRate = cell(M,1);
phCompl = cell(M,1);
for i = 1:M
    [~, mu_i, phi_i] = Coxian.fitMeanAndSCV(L(i), scv(i));
    phRate{i} = mu_i(:);
    phCompl{i} = phi_i(:);
end

% Initial LD rate multipliers (relative to base rate 1/L(i)): exponential
% single-server guess mu=1, multiserver mu=min(n,m). pfqn_mvald interprets mu
% as a multiplier, so the absolute service rate at n jobs is mu(i,n)/L(i).
mu = zeros(M,N);
for i = 1:M
    mu(i,:) = min(1:N, nservers(i));
end

X = zeros(M,1); Q = zeros(M,1); U = zeros(M,1); C = zeros(M,1);
it = 0;
while it < maxiter
    it = it + 1;
    [XN,QN,UN,CN,~,~,piglob] = pfqn_mvald(L,N,Z,mu);
    % Marginal queue-length distribution at full population, per station.
    Pg = piglob(:,:,end);              % M x (sum(N)+1), Pg(i,k) = P(n_i = k-1)
    mu_new = mu;
    for i = 1:M
        Pi = Pg(i,1:N+1);
        % see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
        lam = zeros(1,N);              % lam(n+1) holds lambda_i(n)
        for n = 0:N-1
            pn = Pi(n+1);
            if pn > 0
                lam(n+1) = (mu(i,n+1)/L(i)) * Pi(n+2) / pn;
            else
                lam(n+1) = 0;
            end
        end
        % Isolation returns absolute conditional throughput; convert to the
        % multiplier pfqn_mvald expects (multiplier = abs_rate * L(i)).
        muabs = isol_condtput(lam, phRate{i}, phCompl{i}, N, nservers(i));
        mu_new(i,:) = muabs * L(i);
    end
    delta = max(max(abs(mu_new - mu)));
    mu = mu_new;
    X = XN(:); Q = QN(:); U = UN(:); C = CN(:);
    if delta < tol
        break
    end
end
end

function muvec = isol_condtput(lam, phRate, phCompl, N, m)
% Stationary analysis of a lambda(n)/Cox/1(-m) queue in isolation, returning
% the conditional throughput mu(n) = departure rate given n present, n=1..N.
% lam(n+1) = arrival rate when n customers present (n=0..N-1); the customer in
% service advances through Coxian phases (rate phRate(k); completes w.p.
% phCompl(k), else advances to phase k+1). With m servers, the phase-completion
% rate at population n is scaled by min(n,m).
P = numel(phRate);
% State layout: 1 = empty; for n=1..N, k=1..P -> index 1 + (n-1)*P + k.
S = 1 + N*P;
idx = @(n,k) 1 + (n-1)*P + k;
Gq = zeros(S,S);

% From empty: arrival starts a customer in phase 1.
Gq(1, idx(1,1)) = Gq(1, idx(1,1)) + lam(1);

for n = 1:N
    sc = min(n,m);                     % multiserver rate scaling
    for k = 1:P
        r = idx(n,k);
        % Arrival (queueing; in-service phase preserved).
        if n < N
            Gq(r, idx(n+1,k)) = Gq(r, idx(n+1,k)) + lam(n+1);
        end
        compl = phRate(k) * phCompl(k) * sc;     % completion (departure)
        adv   = phRate(k) * (1-phCompl(k)) * sc; % advance to next phase
        if adv > 0 && k < P
            Gq(r, idx(n,k+1)) = Gq(r, idx(n,k+1)) + adv;
        end
        if compl > 0
            if n > 1
                Gq(r, idx(n-1,1)) = Gq(r, idx(n-1,1)) + compl;
            else
                Gq(r, 1) = Gq(r, 1) + compl;
            end
        end
    end
end
Gq = Gq - diag(sum(Gq,2));

% Stationary distribution: solve p*Gq = 0, sum(p) = 1.
A = [Gq'; ones(1,S)];
b = [zeros(S,1); 1];
p = (A \ b)';

muvec = zeros(1,N);
for n = 1:N
    Pn = 0; dep = 0;
    for k = 1:P
        pk = p(idx(n,k));
        Pn = Pn + pk;
        dep = dep + pk * phRate(k) * phCompl(k) * min(n,m);
    end
    if Pn > 0
        muvec(n) = dep / Pn;
    else
        muvec(n) = min(n,m) / sum(1./phRate); % fallback: exponential-equiv rate
    end
end
end

% ========================= multiclass (R>1) path =========================
function [X,Q,U,C,it,mu] = marie_multi(L,N,Z,scv,varargin)
% Marie's method for multiclass FCFS Coxian closed networks. The aggregate is
% QD-AMVA with class-dependent scaling beta_{i,r}(nvec) = (Coxian isolation
% conditional throughput)/(exponential isolation conditional throughput), so
% beta==1 recovers standard FCFS AMVA and the cd-scaling carries only the
% non-exponential correction. beta is supplied by a multiclass Cox/1 isolation
% sub-model fed the aggregate per-class throughput (Baynat-Dallery isolation),
% iterated to a fixed point on X.
[M,R] = size(L);
N = N(:)';
if nargin < 3 || isempty(Z), Z = zeros(1,R); end
Z = Z(:)';
if nargin < 4 || isempty(scv), scv = ones(M,R); end

tol = 1e-8; maxiter = 1000;
if numel(varargin) >= 1 && ~isempty(varargin{1}), tol = varargin{1}; end
if numel(varargin) >= 2 && ~isempty(varargin{2}), maxiter = varargin{2}; end

% Exact product-form dispatch: exponential service that is also class-
% independent at every station is genuine BCMP FCFS -> exact MVA.
isPF = all(scv(:) == 1);
if isPF
    for i = 1:M
        if max(L(i,:)) - min(L(i,:)) > 1e-12
            isPF = false; break
        end
    end
end
if isPF
    [XN,QN,UN,CN] = pfqn_mva(L,N,Z);
    X = XN(:)'; Q = QN; U = UN; C = CN; it = 0; mu = {}; return
end

% Per-station per-class Coxian phase representation, plus an exponential
% reference (same means) used to normalize the cd-scaling.
phR = cell(M,R); phP = cell(M,R);
eR  = cell(M,R); eP  = cell(M,R);
for i = 1:M
    for r = 1:R
        [~, mir, pir] = Coxian.fitMeanAndSCV(L(i,r), scv(i,r));
        phR{i,r} = mir(:); phP{i,r} = pir(:);
        eR{i,r} = 1/L(i,r); eP{i,r} = 1;   % exponential reference
    end
end

cds = cell(M,1);
for i = 1:M, cds{i} = @(nv) ones(1,R); end   % beta = 1 initially

Xprev = inf(1,R); it = 0;
X = zeros(1,R); Q = zeros(M,R); U = zeros(M,R); C = zeros(M,R);
while it < maxiter
    it = it + 1;
    [X,Q,U,C] = amva_qd(L,N,Z,cds);
    for i = 1:M
        muCox = isol_mc(X, phR(i,:), phP(i,:), N);
        muExp = isol_mc(X, eR(i,:),  eP(i,:),  N);
        cds{i} = make_cdscale(muCox, muExp, N);
    end
    if max(abs(X - Xprev)) < tol, break, end
    Xprev = X;
end
mu = cds;
end

function [X,Q,U,C] = amva_qd(L,N,Z,cds)
% Multiclass Schweitzer AMVA with class-dependent service-rate scaling. The
% effective class-r demand at station i is L(i,r)/beta_{i,r}(nvec_arrival),
% beta supplied by cds{i} evaluated at the arrival-instant per-class population.
[M,R] = size(L);
Q = repmat(N,M,1) / max(M,1);
W = zeros(M,R); X = zeros(1,R); U = zeros(M,R);
Qprev = Q + 1; tol = 1e-9; it = 0;
while max(abs(Q(:)-Qprev(:))) > tol && it < 5000
    it = it + 1; Qprev = Q;
    for r = 1:R
        for i = 1:M
            nv = Q(i,:);
            if N(r) > 0
                nv(r) = Q(i,r) * (N(r)-1) / N(r);   % arrival instant, tagged class
            end
            be = cds{i}(nv);
            % Work-based multiclass FCFS AMVA residence: the tagged class-r job's
            % own effective service plus the effective work of the jobs found
            % ahead (per-class, so unequal means are handled), with the cd
            % scaling beta_{i,s} applied to each class's effective demand.
            Leff = L(i,:) ./ be;
            W(i,r) = Leff(r) + sum(Leff .* nv);
        end
        denom = Z(r) + sum(W(:,r));
        if denom > 0, X(r) = N(r) / denom; else X(r) = 0; end
        for i = 1:M
            Q(i,r) = X(r) * W(i,r);
        end
    end
end
C = W;
for r = 1:R
    for i = 1:M
        U(i,r) = X(r) * L(i,r);   % busy fraction (true mean service)
    end
end
end

function mumat = isol_mc(lam, phRrow, phProw, Nvec)
% Stationary analysis of a multiclass lambda_r/Cox/1 FCFS queue in isolation
% over the joint per-class population box [0..Nvec], with the head-of-line job
% tracked as (class, phase) and, on a departure, the next head class drawn in
% random order (prob n_c/sum(n)). Returns mumat{r}, an ND array over the box
% giving the conditional class-r throughput mu_r(nvec) = (class-r departure
% rate in states with population nvec)/P(nvec).
R = numel(lam);
Pc = zeros(1,R);
for r = 1:R, Pc(r) = numel(phRrow{r}); end
boxsz = Nvec + 1;
npops = prod(boxsz);

% Enumerate states: id map keyed by (popLinear, head, phase).
% state 1 reserved for the empty station.
key2id = containers.Map('KeyType','char','ValueType','double');
ids = {};                 % ids{s} = [popLinear, head, phase]
key2id('E') = 1; ids{1} = [1, 0, 0];
nid = 1;
popsubs = cell(1,R);
for p = 1:npops
    [popsubs{:}] = ind2sub(boxsz, p);
    nvec = cell2mat(popsubs) - 1;   % actual populations
    if sum(nvec) == 0, continue, end
    for c = 1:R
        if nvec(c) > 0
            for k = 1:Pc(c)
                nid = nid + 1;
                key2id(sprintf('%d_%d_%d', p, c, k)) = nid;
                ids{nid} = [p, c, k];
            end
        end
    end
end
S = nid;

    function id = getid(p, c, k)
        if c == 0
            id = 1;
        else
            id = key2id(sprintf('%d_%d_%d', p, c, k));
        end
    end
    function p = poplin(nvec)
        sub = num2cell(nvec + 1);
        p = sub2ind(boxsz, sub{:});
    end

I = zeros(0,1); J = zeros(0,1); V = zeros(0,1);
    function addrate(a, b, rate)
        I(end+1,1) = a; J(end+1,1) = b; V(end+1,1) = rate; %#ok<AGROW>
    end

for s = 1:S
    info = ids{s};
    p = info(1); c = info(2); k = info(3);
    if c == 0
        nvec = zeros(1,R);
    else
        [popsubs{:}] = ind2sub(boxsz, p);
        nvec = cell2mat(popsubs) - 1;
    end
    % arrivals
    for r = 1:R
        if nvec(r) < Nvec(r) && lam(r) > 0
            nnew = nvec; nnew(r) = nnew(r) + 1;
            if c == 0
                addrate(s, getid(poplin(nnew), r, 1), lam(r));  % start service
            else
                addrate(s, getid(poplin(nnew), c, k), lam(r));  % queue behind head
            end
        end
    end
    if c == 0, continue, end
    rate = phRrow{c}(k);
    compl = rate * phProw{c}(k);
    adv   = rate * (1 - phProw{c}(k));
    if adv > 0 && k < Pc(c)
        addrate(s, getid(p, c, k+1), adv);
    end
    if compl > 0
        nnew = nvec; nnew(c) = nnew(c) - 1;
        if sum(nnew) == 0
            addrate(s, 1, compl);
        else
            tot = sum(nnew);
            for cp = 1:R
                if nnew(cp) > 0
                    addrate(s, getid(poplin(nnew), cp, 1), compl * nnew(cp)/tot);
                end
            end
        end
    end
end

Gq = sparse(I, J, V, S, S);
Gq = Gq - spdiags(sum(Gq,2), 0, S, S);

% Stationary distribution.
A = [Gq'; ones(1,S)];
b = [zeros(S,1); 1];
pvec = A \ b;

% Conditional class-r throughput on the population lattice.
mumat = cell(1,R);
for r = 1:R, mumat{r} = zeros(boxsz); end
Ppop = zeros(boxsz);
dep = cell(1,R);
for r = 1:R, dep{r} = zeros(boxsz); end
for s = 1:S
    info = ids{s};
    p = info(1); c = info(2); k = info(3);
    if c == 0, continue, end
    Ppop(p) = Ppop(p) + pvec(s);
    dep{c}(p) = dep{c}(p) + pvec(s) * phRrow{c}(k) * phProw{c}(k);
end
for r = 1:R
    idxpos = Ppop > 0;
    mumat{r}(idxpos) = dep{r}(idxpos) ./ Ppop(idxpos);
end
end

function f = make_cdscale(muCox, muExp, Nvec)
% Class-dependent scaling beta_{i,r}(nv) = muCox_r(nv)/muExp_r(nv), so beta==1
% for exponential service. Multilinear interpolation over the population
% lattice; guarded and clamped.
f = @(nv) cdscale_eval(nv, muCox, muExp, Nvec);
end

function be = cdscale_eval(nv, muCox, muExp, Nvec)
R = numel(Nvec);
be = ones(1,R);
for r = 1:R
    num = ndlininterp(muCox{r}, nv, Nvec);
    den = ndlininterp(muExp{r}, nv, Nvec);
    if den > 0 && num > 0 && isfinite(num) && isfinite(den)
        be(r) = num / den;
    else
        be(r) = 1;
    end
end
be = min(max(be, 1e-3), 1e3);
end

function v = ndlininterp(A, x, Nvec)
% Multilinear interpolation of ND array A (size Nvec+1) at real point x,
% clamped to the box [0, Nvec].
R = numel(Nvec);
sz = Nvec + 1;
x = min(max(x(:)', 0), Nvec);
lo = floor(x);
hi = min(lo + 1, Nvec);
fr = x - lo;
v = 0;
for mask = 0:(2^R - 1)
    w = 1; sub = zeros(1,R);
    for d = 1:R
        if bitget(mask, d)
            sub(d) = hi(d); w = w * fr(d);
        else
            sub(d) = lo(d); w = w * (1 - fr(d));
        end
    end
    if w == 0, continue, end
    subc = num2cell(sub + 1);
    v = v + w * A(sub2ind(sz, subc{:}));
end
end
