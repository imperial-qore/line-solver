function pas_cyclic_ctmc_vs_bruteforce()
% PAS_CYCLIC_CTMC_VS_BRUTEFORCE  Validate LINE's CTMC solver on a closed
% cyclic network of two pass-and-swap (PAS) queues against the exact
% product-form (order-independent) brute-force enumeration.
%
%   Topology:   --> Queue1 --> Queue2 -->   (cyclic, R closed classes)
%
% PAS queues are order-independent, so the network is product-form. The
% brute force enumerates the whole ordered state space to obtain the
% normalizing constant G, hence exact per-class throughputs and mean
% queue lengths, which are compared against CTMC.

if exist('Network','class') ~= 8
    lineStart;
end

% ---- shared model parameters ------------------------------------------
% Both stations must be valid order-independent (OI) queues: the TOTAL
% service rate mu_i(c) depends only on the customer multiset.
%   station 1: M/M/k OI queue, class-independent  -> mu1(c) = min(n,k1)*s1
%   station 2: infinite-server, class-dependent   -> mu2(c) = sum_j beta2(c_j)
K     = [2 2];              % per-class closed populations
s1    = 1.0;  k1 = 2;       % station 1: 2-server OI (M/M/2)
beta2 = [1.5 1.0];          % station 2: per-class infinite-server rates
R     = numel(K);

mu1 = @(c) min(numel(c(c>0)), k1) * s1;
mu2 = @(c) sum(beta2(c(c>0)));

% ---- (1) brute-force product-form reference ---------------------------
[Gbf, Xbf, Qbf] = bruteforce(K, mu1, mu2);

fprintf('Brute-force normalizing constant G = %.12g\n\n', Gbf);

% ---- (2) LINE closed PAS network solved by CTMC -----------------------
model = Network('PAScyclic');
q1 = Queue(model, 'PASQueue1', SchedStrategy.PAS);
q2 = Queue(model, 'PASQueue2', SchedStrategy.PAS);

jobclass = cell(1,R);
for r = 1:R
    jobclass{r} = ClosedClass(model, sprintf('Class%d', r), K(r), q1);
end

% Same total service-rate functions as the brute force; empty swap graph
% (plain order-independent queue -- result is swap-graph invariant).
q1.setService(mu1);
q2.setService(mu2);
q1.setSwapGraph(zeros(R));   q1.setNumberOfServers(k1);      q1.setCap(sum(K));
q2.setSwapGraph(zeros(R));   q2.setNumberOfServers(sum(K));  q2.setCap(sum(K));

% Cyclic routing q1 -> q2 -> q1 for every class.
P = model.initRoutingMatrix;
for r = 1:R
    P{jobclass{r}}(q1, q2) = 1.0;
    P{jobclass{r}}(q2, q1) = 1.0;
end
model.link(P);

T = CTMC(model, 'cutoff', sum(K)).getAvgTable;
disp(T);

% ---- (3) compare ------------------------------------------------------
% Cyclic network with one visit per station => per-station per-class
% throughput equals the chain throughput X_r; mean queue length per
% (station,class) is Q_{i,r}.
Xctmc = zeros(1,R);  Qctmc = zeros(2,R);
stationName = {'PASQueue1','PASQueue2'};
for r = 1:R
    cname = sprintf('Class%d', r);
    for i = 1:2
        row = strcmp(string(T.Station), stationName{i}) & ...
              strcmp(string(T.JobClass), cname);
        Qctmc(i,r) = T.QLen(row);
        if i == 1
            Xctmc(r) = T.Tput(row);
        end
    end
end

fprintf('\n--- Throughput (per class) ---\n');
fprintf('%-8s %14s %14s %10s\n', 'class', 'brute-force', 'CTMC', 'abs.err');
for r = 1:R
    fprintf('Class%-3d %14.10f %14.10f %10.2e\n', ...
        r, Xbf(r), Xctmc(r), abs(Xbf(r)-Xctmc(r)));
end

fprintf('\n--- Mean queue length (station, class) ---\n');
fprintf('%-10s %-8s %14s %14s %10s\n', 'station','class','brute-force','CTMC','abs.err');
for i = 1:2
    for r = 1:R
        fprintf('%-10s Class%-3d %14.10f %14.10f %10.2e\n', ...
            stationName{i}, r, Qbf(i,r), Qctmc(i,r), abs(Qbf(i,r)-Qctmc(i,r)));
    end
end

tol = 1e-9;
errX = max(abs(Xbf - Xctmc));
errQ = max(abs(Qbf(:) - Qctmc(:)));
fprintf('\nCTMC:  max|dX| = %.3e   max|dQ| = %.3e\n', errX, errQ);
assert(errX <= tol && errQ <= tol, ...
    'CTMC vs brute-force mismatch: max|dX|=%.3e max|dQ|=%.3e', errX, errQ);
fprintf('PASS: CTMC matches brute-force product form within %.1e.\n', tol);

% ---- (4) LDES discrete-event simulation (matches within sim noise) ----
% The same closed PAS network solved by the LINE Discrete Event Simulator.
% LDES is a stochastic simulator, so agreement is expected only to within
% simulation noise (Monte-Carlo confidence interval), not to machine eps.
Tl = LDES(model, 'samples', 2e5, 'seed', 23000, 'verbose', false).getAvgTable;
Xldes = zeros(1,R);  Qldes = zeros(2,R);
for r = 1:R
    cname = sprintf('Class%d', r);
    for i = 1:2
        row = strcmp(string(Tl.Station), stationName{i}) & ...
              strcmp(string(Tl.JobClass), cname);
        Qldes(i,r) = Tl.QLen(row);
        if i == 1
            Xldes(r) = Tl.Tput(row);
        end
    end
end

fprintf('\n--- LDES vs brute-force (simulation, 2e5 samples) ---\n');
fprintf('%-12s %-8s %14s %14s %10s\n','metric','class','brute-force','LDES','rel.err');
for r = 1:R
    fprintf('Tput         Class%-3d %14.6f %14.6f %9.2f%%\n', ...
        r, Xbf(r), Xldes(r), 100*abs(Xbf(r)-Xldes(r))/Xbf(r));
end
for i = 1:2
    for r = 1:R
        fprintf('QLen %-7s Class%-3d %14.6f %14.6f %9.2f%%\n', ...
            stationName{i}, r, Qbf(i,r), Qldes(i,r), 100*abs(Qbf(i,r)-Qldes(i,r))/Qbf(i,r));
    end
end

tolSim = 0.02;                       % 2% relative tolerance (simulation)
relX = max(abs(Xbf - Xldes) ./ Xbf);
relQ = max(abs(Qbf(:) - Qldes(:)) ./ Qbf(:));
fprintf('\nLDES:  max rel|dX| = %.2f%%   max rel|dQ| = %.2f%%\n', 100*relX, 100*relQ);
assert(relX <= tolSim && relQ <= tolSim, ...
    'LDES vs brute-force exceeds %.0f%%: relX=%.2f%% relQ=%.2f%%', ...
    100*tolSim, 100*relX, 100*relQ);
fprintf('PASS: LDES matches brute-force within %.0f%% (simulation noise).\n', 100*tolSim);
end

% =======================================================================
function [G, X, Q] = bruteforce(K, mu1, mu2)
% Exact product-form normalizing constant G, per-class throughputs X
% (X_r = G(K-1_r)/G(K)), and mean queue lengths Q(i,r) by full ordered
% state-space enumeration of the two order-independent stations.
R = numel(K);
e1 = ones(1,R);  e2 = ones(1,R);
G = norm_const(K, e1, e2, mu1, mu2);

X = zeros(1,R);
for r = 1:R
    Km = K; Km(r) = Km(r) - 1;
    X(r) = norm_const(Km, e1, e2, mu1, mu2) / G;   % e_r = 1
end

% Q(1,r) = (1/G) sum_m m_r S1(m) S2(K-m);  Q(2,r) = K_r - Q(1,r).
Q = zeros(2,R);
splits = enum_splits(K);
for s = 1:size(splits,1)
    m  = splits(s,:);
    S1 = oi_sum(m,     e1, mu1);
    S2 = oi_sum(K - m, e2, mu2);
    Q(1,:) = Q(1,:) + m .* (S1 * S2);
end
Q(1,:) = Q(1,:) / G;
Q(2,:) = K - Q(1,:);
end

function G = norm_const(K, e1, e2, mu1, mu2)
G = 0;
splits = enum_splits(K);
for s = 1:size(splits,1)
    m = splits(s,:);
    G = G + oi_sum(m, e1, mu1) * oi_sum(K - m, e2, mu2);
end
end

function S = oi_sum(counts, e, mu)
S = oi_dfs(counts, [], 1.0, e, mu);
end

function S = oi_dfs(counts, prefix, w, e, mu)
if all(counts == 0)
    S = w; return;
end
S = 0;
for r = 1:numel(counts)
    if counts(r) > 0
        c2 = counts;  c2(r) = c2(r) - 1;
        p2 = [prefix, r];
        S  = S + oi_dfs(c2, p2, w * e(r) / mu(p2), e, mu);
    end
end
end

function m = pas_rate(c, beta, k)
% Total OI service rate of the ordered prefix c (zero-padded, 1-based
% class ids) for a k-server PAS station.
cc = c(c > 0);
n  = numel(cc);
if n == 0
    m = 0; return;
end
served = cc(1:min(n,k));
m = sum(beta(served));
end

function splits = enum_splits(K)
R = numel(K);
nrows = prod(K + 1);
splits = zeros(nrows, R);
idx = zeros(1, R);
for s = 1:nrows
    splits(s,:) = idx;
    for r = 1:R
        if idx(r) < K(r)
            idx(r) = idx(r) + 1; break;
        else
            idx(r) = 0;
        end
    end
end
end
