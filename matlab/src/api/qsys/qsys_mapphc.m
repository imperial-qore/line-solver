function result = qsys_mapphc(D0, D1, alpha, S, c, varargin)
% QSYS_MAPPHC Analyzes a MAP/PH/c FCFS queue exactly.
%
% RESULT = QSYS_MAPPHC(D0, D1, ALPHA, S, C) analyzes a c-server FCFS queue
% with:
%   D0, D1 - arrival MAP of order ma
%   ALPHA  - PH service initial vector (1 x ms)
%   S      - PH service sub-generator (ms x ms)
%   C      - number of servers (identical, so the service law is shared)
%
% RESULT = QSYS_MAPPHC(..., 'maxNumComp', N) caps the number of queue length
% probabilities returned (default 500).
% RESULT = QSYS_MAPPHC(..., 'numWMoms', K) returns K waiting-time moments
% (default 3).
% RESULT = QSYS_MAPPHC(..., 'wPoints', T) evaluates the waiting-time CCDF at
% the times T.
%
% THE STATE SPACE. With c identical servers the server identities carry no
% information, so the service phases are held as a MULTISET: a configuration is
% n = (n_1..n_ms) with sum(n) = k servers busy in phase i. There are
% nchoosek(ms+k-1,k) of them, the count of Asmussen and Moller (2001), against
% ms^k for the ordered space. Levels 0..c-1 are the boundary (level = servers
% busy), levels >= c repeat and carry the queue.
%
% THE WAITING TIME. An arrival that finds j customers waiting ahead of it waits
% for j+1 service completions, so Wq is the (j+1)-st event time of the
% configuration MAP (Lc, Cdep) started at the arrival-epoch configuration.
% Folding the matrix-geometric level distribution over j gives
%
%     G'(t) = G(t) Lj + R G(t) Cj,   G(0) = (I-R)^-1 kron(D1,I) / lambda,
%     P(Wq > t) = pi_c G(t) e,
%
% a LINEAR matrix ODE, so Wq is matrix-exponential. Its transform obeys the
% generalized Sylvester equation g(sI-Lj) - R g Cj = G(0), and every moment
% reuses that one operator with a different right-hand side.
%
% Returns a struct with fields:
%   meanQueueLength    - E[N], number in system
%   meanWaitingTime    - E[Wq], time in queue
%   meanSojournTime    - E[W] = E[Wq] + E[service]
%   utilization        - rho = lambda E[service] / c, per server
%   queueLengthDist    - P(N = n), n = 0, 1, ...
%   waitingTimeMoments - E[Wq^k], k = 1..numWMoms
%   waitingTimeCCDF    - P(Wq > t) at the requested wPoints
%   probWait           - P(Wq > 0), an arrival finds every server busy
%   phaseCount         - nchoosek(ms+c-1,c), the repeating configuration count
%   analyzer           - name of the analyzer used
%
% References:
%   S. Asmussen and J.R. Moller, "Calculation of the steady state waiting time
%   distribution in GI/PH/c and MAP/PH/c queues", Queueing Systems 37(1):9-29,
%   2001.
%   D.P. Gaver, P.A. Jacobs, G. Latouche, "Finite birth-and-death models in
%   randomly changing environments", Adv. Appl. Probab. 16:715-731, 1984.
%
% See also qsys_mapmc, qsys_mapph1, qsys_phmc, ph_multisets

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

p = inputParser;
addParameter(p, 'maxNumComp', 500);
addParameter(p, 'numWMoms', 3);
addParameter(p, 'wPoints', []);
parse(p, varargin{:});
maxNumComp = p.Results.maxNumComp;
numWMoms = p.Results.numWMoms;
wPoints = p.Results.wPoints(:).';

D0 = full(D0); D1 = full(D1);
alpha = full(alpha(:)).';
S = full(S);
c = round(c);
ma = size(D0, 1);
ms = size(S, 1);
if size(D0, 2) ~= ma || size(D1, 1) ~= ma || size(D1, 2) ~= ma
    line_error(mfilename, 'D0 and D1 must be square and of equal order.');
end
if size(S, 2) ~= ms || numel(alpha) ~= ms
    line_error(mfilename, 'alpha and S must have matching order.');
end
if c < 1
    line_error(mfilename, 'The number of servers c must be at least one.');
end

s0 = -S * ones(ms, 1);

theta = qsys_mapphc_stat(D0 + D1);
lambda = theta * D1 * ones(ma, 1);
meanService = -alpha * (S \ ones(ms, 1));
rho = lambda * meanService / c;
if rho >= 1 - eps
    line_error(mfilename, 'The load %g of the system is not below one.', rho);
end

% Multiset configurations of every occupancy 0..c
cfg = cell(c + 1, 1);
for k = 0:c
    cfg{k + 1} = ph_multisets(ms, k);
end

% Service blocks per occupancy
Lcfg = cell(c + 1, 1);   % phase changes, no completion
Up = cell(c + 1, 1);     % k -> k+1, an arrival enters service
Dn = cell(c + 1, 1);     % k -> k-1, a completion with nobody waiting
for k = 0:c
    Ck = cfg{k + 1};
    nk = size(Ck, 1);
    Lk = zeros(nk, nk);
    for row = 1:nk
        n = Ck(row, :);
        for i = 1:ms
            if n(i) == 0, continue; end
            for j = 1:ms
                if j == i, continue; end
                m = n; m(i) = m(i) - 1; m(j) = m(j) + 1;
                col = qsys_mapphc_find(Ck, m);
                Lk(row, col) = Lk(row, col) + n(i) * S(i, j);
            end
            Lk(row, row) = Lk(row, row) + n(i) * S(i, i);
        end
    end
    Lcfg{k + 1} = Lk;

    if k < c
        Ck1 = cfg{k + 2};
        Uk = zeros(nk, size(Ck1, 1));
        for row = 1:nk
            n = Ck(row, :);
            for j = 1:ms
                m = n; m(j) = m(j) + 1;
                col = qsys_mapphc_find(Ck1, m);
                Uk(row, col) = Uk(row, col) + alpha(j);
            end
        end
        Up{k + 1} = Uk;
    end

    if k > 0
        Ckm = cfg{k};
        Dk = zeros(nk, size(Ckm, 1));
        for row = 1:nk
            n = Ck(row, :);
            for i = 1:ms
                if n(i) == 0, continue; end
                m = n; m(i) = m(i) - 1;
                col = qsys_mapphc_find(Ckm, m);
                Dk(row, col) = Dk(row, col) + n(i) * s0(i);
            end
        end
        Dn{k + 1} = Dk;
    end
end

% Completion WITH an immediate restart: the repeating down block
Cc = cfg{c + 1};
nc = size(Cc, 1);
Cdep = zeros(nc, nc);
for row = 1:nc
    n = Cc(row, :);
    for i = 1:ms
        if n(i) == 0, continue; end
        for j = 1:ms
            m = n; m(i) = m(i) - 1; m(j) = m(j) + 1;
            col = qsys_mapphc_find(Cc, m);
            Cdep(row, col) = Cdep(row, col) + n(i) * s0(i) * alpha(j);
        end
    end
end

% Repeating QBD blocks, level = number in system >= c
Ima = eye(ma); Inc = eye(nc);
A_up = kron(D1, Inc);
A_loc = kron(D0, Inc) + kron(Ima, Lcfg{c + 1});
A_dn = kron(Ima, Cdep);
R = qbd_R_logred(A_dn, A_loc, A_up);

% Boundary levels 0..c, with the tail folded into level c through R
sz = zeros(c + 1, 1);
for k = 0:c
    sz(k + 1) = ma * size(cfg{k + 1}, 1);
end
off = [0; cumsum(sz)];
tot = off(end);
Q = zeros(tot, tot);
for k = 0:c
    rk = off(k + 1) + (1:sz(k + 1));
    Ick = eye(size(cfg{k + 1}, 1));
    if k < c
        Q(rk, rk) = kron(D0, Ick) + kron(Ima, Lcfg{k + 1});
        Q(rk, off(k + 2) + (1:sz(k + 2))) = kron(D1, Up{k + 1});
    else
        Q(rk, rk) = A_loc + R * A_dn;
    end
    if k > 0
        Q(rk, off(k) + (1:sz(k))) = kron(Ima, Dn{k + 1});
    end
end
piVec = qsys_mapphc_stat(Q);

nOp = nc * ma;
ImR = eye(nOp) - R;
piC = piVec(off(c + 1) + (1:sz(c + 1)));
massBoundary = sum(piVec(1:off(c + 1)));
massTail = piC * (ImR \ ones(nOp, 1));
piVec = piVec / (massBoundary + massTail);
piC = piVec(off(c + 1) + (1:sz(c + 1)));

% Queue length distribution
ql = zeros(1, c);
for k = 0:c - 1
    ql(k + 1) = sum(piVec(off(k + 1) + (1:sz(k + 1))));
end
tail = piC;
qlList = sum(tail);
acc = sum(ql) + qlList;
while acc < 1 - 1e-12 && numel(qlList) < maxNumComp - c + 1
    tail = tail * R;
    qlList(end + 1) = sum(tail); %#ok<AGROW>
    acc = acc + qlList(end);
end
qlBnd = ql;
ql = [ql, qlList];
% E[N] in CLOSED FORM. maxNumComp caps the probabilities RETURNED, not the mean:
% summing the truncated list loses the matrix-geometric tail, which at rho -> 1
% carries a first-order share of the mass. With pi_{c+j} = pi_c R^j,
% sum_j (c+j) pi_c R^j e = pi_c [c (I-R)^-1 + R (I-R)^-2] e.
u = ImR \ ones(nOp, 1);
meanQL = sum((0:c - 1) .* qlBnd) + piC * (c * u + R * (ImR \ u));

% Waiting time
Lj = kron(Ima, Lcfg{c + 1});
Cj = kron(Ima, Cdep);
G0 = (ImR \ kron(D1, Inc)) / lambda;
probWait = piC * G0 * ones(nOp, 1);

Iop = eye(nOp);
% X(-Lj) - R X Cj = rhs, vectorized column-major
Kop = kron((-Lj).', Iop) - kron(Cj.', R);
wMoms = zeros(1, numWMoms);
gPrev = [];
for k = 1:numWMoms
    if k == 1
        rhs = G0;
    else
        rhs = -(k - 1) * gPrev;
    end
    g = reshape(Kop \ rhs(:), nOp, nOp);
    wMoms(k) = k * ((-1)^(k - 1)) * (piC * g * ones(nOp, 1));
    gPrev = g;
end
meanWT = wMoms(1);

% CCDF at the requested points, by the vectorized linear ODE
wCCDF = [];
if ~isempty(wPoints)
    Kt = kron(Lj.', Iop) + kron(Cj.', R);
    v0 = G0(:);
    wCCDF = zeros(1, numel(wPoints));
    for it = 1:numel(wPoints)
        Gt = reshape(expm(Kt * wPoints(it)) * v0, nOp, nOp);
        wCCDF(it) = piC * Gt * ones(nOp, 1);
    end
end

result = struct();
result.meanQueueLength = meanQL;
result.meanWaitingTime = meanWT;
result.meanSojournTime = meanWT + meanService;
result.utilization = rho;
result.queueLengthDist = ql;
result.waitingTimeMoments = wMoms;
result.waitingTimeCCDF = wCCDF;
result.waitingTimePoints = wPoints;
result.probWait = probWait;
result.phaseCount = nc;
result.analyzer = sprintf('LINE:MAP/PH/%d', c);
end

function K = qsys_mapphc_stat(G)
% Left null vector of G normalized to sum one. The R-corrected level-c block
% has nonzero row sums, so the augmented system is used rather than a generator
% solve.
n = size(G, 1);
B = [G, ones(n, 1)];
y = [zeros(1, n), 1];
K = y / B;
end

function idx = qsys_mapphc_find(rows, m)
[~, idx] = ismember(m, rows, 'rows');
end
