function result = qsys_mmapgk1(MMAP, svc, varargin)
% QSYS_MMAPGK1 Per-type waiting times of the MMAP[K]/G[K]/1 FCFS queue.
%
% RESULT = QSYS_MMAPGK1(MMAP, SVC) analyzes a single-server FCFS queue with:
%   MMAP - cell {D0, D1, ..., DK}: a marked Markovian arrival process of order
%          ma whose type-k arrivals carry the block Dk
%   SVC  - cell of K service laws, one per type. Each entry is a LINE
%          Distribution object, a cell {D0s, D1s} holding a phase-type
%          representation, or a struct with fields
%             .lst      function handle s -> E[exp(-s S)], admitting COMPLEX s
%             .moments  raw moments [E[S], E[S^2], ...], at least numWMoms+1 of
%                       them
%          The struct form is what a solver hands in: sn.lst already carries the
%          transform of the ORIGINAL law, while sn.proc carries only its
%          phase-type fit, and it is the original that this analysis needs. The
%          service laws may be GENERAL and need not share a family across types.
%
% RESULT = QSYS_MMAPGK1(..., 'wPoints', T) also evaluates the per-type waiting
% time CDF at the times T by numerical transform inversion.
% RESULT = QSYS_MMAPGK1(..., 'numWMoms', N) returns N moments per type
% (default 3).
% RESULT = QSYS_MMAPGK1(..., 'quadNodes', Q) sets the quadrature order used for
% a service law with no closed-form matrix transform (default 120).
%
% THE METHOD, which is He's, theorem for theorem. FCFS makes the actual waiting
% time of a customer the WORKLOAD it finds on arrival, so everything follows
% from the joint transform of workload and arrival phase,
% f(s)_j = E[exp(-s V) 1{phase = j}], which by He's Theorem 4.1 (eq. 4.6)
% satisfies
%
%     f(s) [ s I + D0 + sum_k Dk gk(s) ] = s v0,                          (*)
%
% with v0 the idle-phase vector, his y0. The unknown v0 needs NO search for the
% roots of the determinant: the matrix U solving
%
%     U = D0 + sum_k Dk Fk(U),      Fk(U) = int_0^inf exp(U t) dFk(t),
%
% is his eq. (4.4), the generator of the underlying Markov process obtained by
% EXCISING the busy periods, and eq. (4.5) with Theorem 4.2 give y0 Q = 0 and
% y0 e = 1 - rho, i.e.
%
%     v0 = (1 - rho) pi_U.
%
% The same vector is what the analyticity of (*) forces: for every left
% eigenpair (w, u) of U one has w [D0 + sum_k Dk gk(-u) + (-u) I] = 0, so the
% roots of (*) in the closed right half plane are exactly s = -u over the
% spectrum of U, and imposing v0 r_i = 0 at each right null vector reproduces
% the stationary vector to 2.5e-13. The stationary route is the one taken, as
% it needs no complex eigenvector and no rule for telling the structural root at
% the origin from a genuine one.
%
% The per-type actual waiting time is the workload seen by a type-k arrival,
% biased by that type's own arrival block, which is his Theorem 5.1 eq. (5.1)
% summed over the post-arrival phase:
%
%     E[exp(-s Wk)] = f(s) Dk e / lambda_k.
%
% SCOPE. He allows an arrival to be a BATCH carrying a sequence of types, and
% his Theorem 5.3 then multiplies the transform by prod_{i<n} f*_{h_i}(s), the
% service of the customers ahead of the tagged one WITHIN its own batch. This
% function covers the single-customer-per-arrival case, his Special case 3.3,
% where that product is empty -- which is exactly the MMAP convention LINE
% carries, {D0, D1, D^(1), ..., D^(K)} with one customer per epoch.
%
% Returns a struct with fields:
%   lambda            - per-type arrival rates (1 x K)
%   utilization       - rho = sum_k lambda_k E[S_k]
%   idleVector        - v0, the idle-phase vector, summing to 1 - rho
%   waitLST           - function handle s -> row vector of E[exp(-s Wk)]
%   waitMoments       - (K x numWMoms) per-type waiting time moments
%   meanWaitingTime   - per-type mean waiting time (1 x K)
%   meanSojournTime   - per-type mean sojourn time (1 x K)
%   meanQueueLength   - E[N], by Little over all types
%   waitCDF           - (K x numel(wPoints)) per-type waiting time CDF
%   waitPoints        - the requested points
%   analyzer          - name of the analyzer used
%
% Reference:
%   Qi-Ming He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1 queue",
%   Queueing Systems 38(4):397-418, 2001.
%
% See also qsys_mapg1, qsys_mapphc, MMAPPH1FCFS

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

p = inputParser;
addParameter(p, 'wPoints', []);
addParameter(p, 'numWMoms', 3);
addParameter(p, 'quadNodes', 120);
addParameter(p, 'tol', 1e-12);
parse(p, varargin{:});
wPoints = p.Results.wPoints(:).';
numWMoms = p.Results.numWMoms;
quadNodes = p.Results.quadNodes;
tol = p.Results.tol;

% LINE's MMAP convention: {D0, D1, D^(1), ..., D^(K)} with D1 = sum_k D^(k)
K = numel(MMAP) - 2;
if K < 1
    line_error(mfilename, 'The MMAP must carry at least one marked arrival block.');
end
if numel(svc) ~= K
    line_error(mfilename, 'One service law per marked type is required (%d given, %d types).', ...
        numel(svc), K);
end
D0 = full(MMAP{1});
ma = size(D0, 1);
Dk = cell(1, K);
for k = 1:K
    Dk{k} = full(MMAP{k + 2});
end
Dsum = D0;
for k = 1:K
    Dsum = Dsum + Dk{k};
end

theta = ctmc_solve(Dsum);
lambda = zeros(1, K);
for k = 1:K
    lambda(k) = theta * Dk{k} * ones(ma, 1);
end
meanS = zeros(1, K);
for k = 1:K
    meanS(k) = qsys_mmapgk1_mean(svc{k});
end
rho = sum(lambda .* meanS);
if rho >= 1 - eps
    line_error(mfilename, 'The load %g of the system is not below one.', rho);
end

% Fixed point U = D0 + sum_k Dk Fk(U). The natural iteration from U = D0 is
% monotone in the same sense as the G-matrix iteration it generalizes.
U = D0;
for it = 1:10000
    Unew = D0;
    for k = 1:K
        Unew = Unew + Dk{k} * qsys_mmapgk1_matlst(svc{k}, U, quadNodes);
    end
    if max(max(abs(Unew - U))) <= tol
        U = Unew;
        break;
    end
    U = Unew;
end
% U e = 0 EXACTLY. It is a property of the fixed point, not of the iterate: the
% iteration converges linearly, so the row sums still carry O(1e-9) at the
% tolerance above, and that residue moves the eigenvalue that belongs at the
% origin off it, which would then be mistaken for an analyticity condition.
U = U - diag(sum(U, 2));

% The idle vector. U is a proper generator: its off-diagonals are nonnegative
% (D0's are, Dk >= 0 and Fk(U) >= 0) and its rows sum to zero, and it governs
% the arrival phase at the epochs the level drops by one. Its stationary vector
% carries the idle mass,
%
%     v0 = (1 - rho) pi_U,        pi_U U = 0,  pi_U e = 1,
%
% which is what the analyticity conditions of (*) deliver: imposing v0 r_i = 0
% at the right null vector of the bracket for each of the ma-1 roots s = -u off
% the origin, plus v0 e = 1 - rho, reproduces this vector to 2.5e-13 on an
% Erlang case and to the fixed point's own residual elsewhere. The stationary
% route is taken because it needs no complex eigenvector and no rule for telling
% the structural root at the origin from a genuine one.
v0 = (1 - rho) * ctmc_solve(U);

% Transforms
fLST = @(s) s * v0 / (s * eye(ma) + qsys_mmapgk1_dz(D0, Dk, svc, s, quadNodes));
waitLST = @(s) qsys_mmapgk1_wlst(fLST, Dk, lambda, ma, s);

% Moments by finite differences of the transform along the real axis are
% ill-conditioned, so they are taken from the derivatives of (*) instead:
% differentiating f(s) M(s) = s v0 k times at s = 0 gives a triangular system.
waitMoments = qsys_mmapgk1_moments(D0, Dk, svc, theta, v0, lambda, numWMoms, quadNodes);

meanWT = waitMoments(:, 1).';
meanST = meanWT + meanS;
lamTot = sum(lambda);
meanQL = sum(lambda .* meanST);

waitCDF = [];
if ~isempty(wPoints)
    waitCDF = zeros(K, numel(wPoints));
    for k = 1:K
        for it = 1:numel(wPoints)
            waitCDF(k, it) = qsys_mmapgk1_invert(@(s) subsref_k(waitLST(s), k), wPoints(it));
        end
    end
end

result = struct();
result.lambda = lambda;
result.arrivalRate = lamTot;
result.utilization = rho;
result.idleVector = v0;
result.waitLST = waitLST;
result.waitMoments = waitMoments;
result.meanWaitingTime = meanWT;
result.meanSojournTime = meanST;
result.meanQueueLength = meanQL;
result.waitCDF = waitCDF;
result.waitPoints = wPoints;
result.analyzer = sprintf('LINE:MMAP[%d]/G[%d]/1', K, K);
end

function v = subsref_k(row, k)
v = row(k);
end

function M = qsys_mmapgk1_dz(D0, Dk, svc, s, quadNodes)
% D0 + sum_k Dk gk(s), the scalar-transform bracket of (*)
M = D0;
for k = 1:numel(Dk)
    M = M + Dk{k} * qsys_mmapgk1_lst(svc{k}, s, quadNodes);
end
end

function w = qsys_mmapgk1_wlst(fLST, Dk, lambda, ma, s)
K = numel(Dk);
w = zeros(1, K);
if abs(s) < 1e-14
    w(:) = 1;
    return;
end
f = fLST(s);
for k = 1:K
    w(k) = (f * Dk{k} * ones(ma, 1)) / lambda(k);
end
end

function m = qsys_mmapgk1_mean(law)
if isstruct(law)
    m = law.moments(1);
elseif iscell(law)
    m = map_mean(law);
else
    m = law.getMean();
end
end

function g = qsys_mmapgk1_lst(law, s, quadNodes)
% Scalar Laplace-Stieltjes transform of the service law at s.
if isstruct(law)
    g = law.lst(s);
    return;
end
if iscell(law)
    D0s = law{1};
    n = size(D0s, 1);
    alpha = map_pie(law);
    g = alpha * ((s * eye(n) - D0s) \ (-D0s * ones(n, 1)));
else
    g = law.evalLST(s);
end
if ~isfinite(g)
    [x, w] = qsys_mmapgk1_quad(law, quadNodes);
    g = sum(w .* exp(-s * x));
end
end

function F = qsys_mmapgk1_matlst(law, U, quadNodes)
% Matrix transform int_0^inf exp(U t) dF(t).
n = size(U, 1);
if isstruct(law)
    % The transform handle alone suffices: diagonalizing U turns the MATRIX
    % transform into the SCALAR one at the eigenvalues, which is why sn.lst has
    % to admit a complex argument.
    [Vd, Dd] = eig(U);
    uv = diag(Dd).';
    gv = arrayfun(@(u) law.lst(-u), uv);
    F = real(Vd * diag(gv) / Vd);
    return;
end
if iscell(law)
    % Phase-type service (beta, S). The density is the SCALAR beta exp(St) s0,
    % so the integral is exact on the Kronecker sum: int exp(Ut) x exp(St) dt =
    % -(U (+) S)^-1, and the transform is that sandwiched by beta and s0.
    D0s = law{1};
    ms = size(D0s, 1);
    beta = map_pie(law);
    s0 = -D0s * ones(ms, 1);
    KS = kron(U, eye(ms)) + kron(eye(n), D0s);
    F = kron(eye(n), beta) * (-(KS \ kron(eye(n), s0)));
    return;
end
if isa(law, 'Markovian')
    F = qsys_mmapgk1_matlst(law.getRepresentation(), U, quadNodes);
    return;
end
if isa(law, 'Det')
    F = expm(U * law.getMean());
    return;
end
% A GENERAL service law needs no quadrature: diagonalizing U turns the matrix
% transform into the SCALAR transform at the eigenvalues, F(U) = V g(-u) V^-1,
% and every LINE distribution carries evalLST. Only a defective U or an evalLST
% that refuses a complex argument falls through to the Stieltjes sum below.
[Vd, Dd] = eig(U);
uv = diag(Dd).';
ok = rcond(Vd) > 1e-12;
if ok
    gv = zeros(1, n);
    for i = 1:n
        try
            gv(i) = law.evalLST(-uv(i));
        catch
            ok = false;
            break;
        end
    end
    if ok && all(isfinite(gv))
        F = real(Vd * diag(gv) / Vd);
        return;
    end
end
% Riemann-Stieltjes fallback: true probability weights from the CDF, midpoint
% nodes, so it is a proper measure for any law including one with an atom.
[x, w] = qsys_mmapgk1_quad(law, quadNodes);
F = zeros(n, n);
for i = 1:numel(x)
    F = F + w(i) * expm(U * x(i));
end
end

function [x, w] = qsys_mmapgk1_quad(law, quadNodes)
% Gauss-Legendre nodes and dF weights over the support of the service law. A
% law with bounded support is integrated over exactly that interval; otherwise
% the tail is cut at twelve standard deviations and the weights renormalized, so
% the quadrature still returns a proper transform.
lo = 0;
if isa(law, 'Uniform')
    lo = law.getParam(1).paramValue;
    hi = law.getParam(2).paramValue;
else
    hi = law.getMean() * 60;
    try
        hi = max(hi, law.getMean() + 12 * sqrt(law.getVar()));
    catch
    end
end
nGrid = max(quadNodes, 20) * 20;
edges = linspace(lo, hi, nGrid + 1);
x = 0.5 * (edges(1:end - 1) + edges(2:end));
w = zeros(1, nGrid);
Fprev = law.evalCDF(edges(1));
for i = 1:nGrid
    Fnext = law.evalCDF(edges(i + 1));
    w(i) = Fnext - Fprev;
    Fprev = Fnext;
end
mass = sum(w);
if mass > 0
    w = w / mass;
end
end

function [moms, workMoms] = qsys_mmapgk1_moments(D0, Dk, svc, theta, v0, lambda, numWMoms, quadNodes)
% Derivatives of f(s) M(s) = s v0 at s = 0. Writing f = sum_j f_j s^j / j! and
% M = sum_j M_j s^j / j! with M_0 = D (a generator), matching orders gives
%
%     sum_{i=0..j} C(j,i) f_i M_{j-i} = [j = 1] v0.
%
% M_0 is SINGULAR with right null vector e, so each order fixes f_j only up to a
% multiple of theta; that multiple is what the NEXT order's solvability
% condition supplies. At j = 0 the same condition reads theta M_1 e = v0 e,
% i.e. 1 - rho = 1 - rho, which is the identity that validates the setup.
ma = size(D0, 1);
K = numel(Dk);
Mder = cell(1, numWMoms + 2);
for j = 0:numWMoms + 1
    if j == 0
        Mj = D0;
        for k = 1:K
            Mj = Mj + Dk{k};
        end
    else
        Mj = zeros(ma, ma);
        for k = 1:K
            Mj = Mj + Dk{k} * qsys_mmapgk1_lstder(svc{k}, j, quadNodes);
        end
        if j == 1
            Mj = Mj + eye(ma);   % the s I term contributes only at first order
        end
    end
    Mder{j + 1} = Mj;
end
e = ones(ma, 1);
denom = theta * Mder{2} * e;   % = 1 - rho
fder = cell(1, numWMoms + 1);
fder{1} = theta;
workMoms = zeros(1, numWMoms);
Abase = [Mder{1}, e];
for j = 1:numWMoms
    rhs = zeros(1, ma);
    if j == 1
        rhs = rhs + v0;
    end
    for i = 0:j - 1
        rhs = rhs - nchoosek(j, i) * (fder{i + 1} * Mder{j - i + 1});
    end
    % particular solution with f_j^p e = 0, then the theta component
    fp = [rhs, 0] / Abase;
    acc = 0;
    for i = 0:j - 1
        acc = acc + nchoosek(j + 1, i) * (fder{i + 1} * Mder{j + 1 - i + 1} * e);
    end
    cfree = (-acc / (j + 1) - fp * Mder{2} * e) / denom;
    fder{j + 1} = fp + cfree * theta;
    workMoms(j) = ((-1)^j) * cfree;   % E[V^j], the virtual waiting time
end
moms = zeros(K, numWMoms);
for k = 1:K
    for j = 1:numWMoms
        % E[Wk^j] = (-1)^j d^j/ds^j E[e^{-s Wk}] at 0, and the type-k transform
        % is f(s) Dk e / lambda_k, so the j-th derivative carries f_j.
        moms(k, j) = ((-1)^j) * (fder{j + 1} * Dk{k} * e) / lambda(k);
    end
end
end

function d = qsys_mmapgk1_lstder(law, j, quadNodes)
% j-th derivative of the scalar LST at s = 0, i.e. (-1)^j E[S^j].
d = ((-1)^j) * qsys_mmapgk1_rawmom(law, j, quadNodes);
end

function m = qsys_mmapgk1_rawmom(law, j, quadNodes)
% Raw moment E[S^j] of a service law.
if isstruct(law)
    if j > numel(law.moments)
        line_error(mfilename, ...
            ['The service law was given as a transform with %d moments, but order %d ' ...
             'is needed. Supply at least numWMoms+1 moments.'], numel(law.moments), j);
    end
    m = law.moments(j);
    return;
end
if iscell(law)
    m = map_moment(law, j);
    return;
end
if isa(law, 'Markovian')
    m = map_moment(law.getRepresentation(), j);
    return;
end
if isa(law, 'Det')
    m = law.getMean()^j;   % the only law whose density quadrature cannot see
    return;
end
if isa(law, 'Uniform')
    a = law.getParam(1).paramValue; b = law.getParam(2).paramValue;
    m = (b^(j + 1) - a^(j + 1)) / ((b - a) * (j + 1));
    return;
end
switch j
    case 1
        m = law.getMean();
    case 2
        m = law.getVar() + law.getMean()^2;
    otherwise
        m = [];
        if j == 3
            try
                m1 = law.getMean(); v = law.getVar(); m2 = v + m1^2;
                m = law.getSkewness() * v^1.5 + 3 * m1 * m2 - 2 * m1^3;
            catch
                m = [];
            end
        end
        if isempty(m)
            % A HEAVY TAIL HAS NO MOMENT, and the truncated sum below would hand
            % back a finite number for one that diverges: a Pareto of shape <= j
            % has E[S^j] = Inf, and E[Wq] is then genuinely infinite rather than
            % merely large. The quadrature cannot see that, since it integrates
            % over a cut support, so the divergence is decided from the tail
            % index before the sum is taken.
            if isa(law, 'Pareto')
                shape = law.getParam(1).paramValue;
                if shape <= j
                    m = Inf;
                    return;
                end
            end
            [x, w] = qsys_mmapgk1_quad(law, quadNodes);
            m = sum(w .* (x .^ j));
        end
end
end

function F = qsys_mmapgk1_invert(lstFun, t)
% Euler summation inversion of a Laplace-Stieltjes transform, giving the CDF
% P(W <= t) from E[exp(-s W)] / s.
if t <= 0
    F = real(lstFun(1e12));
    return;
end
A = 18.4;
nEuler = 15;
mEuler = 11;
u = exp(A / 2) / t;
x = A / (2 * t);
term = zeros(1, nEuler + mEuler + 1);
term(1) = real(lstFun(x)) / x / 2;
for kk = 1:(nEuler + mEuler)
    sk = x + 1i * pi * kk / t;
    term(kk + 1) = ((-1)^kk) * real(lstFun(sk) / sk);
end
partial = cumsum(term);
w = zeros(1, mEuler + 1);
for jj = 0:mEuler
    w(jj + 1) = nchoosek(mEuler, jj) / 2^mEuler;
end
F = u * sum(w .* partial(nEuler + 1:nEuler + mEuler + 1));
F = min(max(F, 0), 1);
end
