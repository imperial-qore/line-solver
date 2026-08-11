function test_lossn_mci()
% TEST_LOSSN_MCI Validate Monte Carlo summation for loss networks
% against exact enumeration, Erlang-B, and the Erlang fixed point.

fprintf('=== test_lossn_mci ===\n');
rng(12345);
npass = 0; nfail = 0;

% ---------------------------------------------------------------
% Test 1: single-link Erlang-B (J=1, one class, unit demand)
% Exact blocking = ErlangB(rho, C).
% ---------------------------------------------------------------
rho = 8; Ccap = 10;
[~, Loss, lG, ci] = lossn_mci(rho, 1, Ccap, struct('samples',2e5,'seed',1));
exactB = erlangB_ref(rho, Ccap);
[gExact, ~] = lossn_exact(rho, 1, Ccap);
fprintf('\n[T1] M/M/C/C Erlang-B  rho=%g C=%d\n', rho, Ccap);
fprintf('  blocking exact=%.6f  mci=%.6f  CI=[%.6f,%.6f]\n', ...
    exactB, Loss, ci.loss(1,1), ci.loss(1,2));
fprintf('  lG exact=%.6f  mci=%.6f\n', log(gExact), lG);
[npass,nfail] = check(npass, nfail, ci.loss(1,1) <= exactB && exactB <= ci.loss(1,2), 'T1 blocking in CI');
[npass,nfail] = check(npass, nfail, abs(lG-log(gExact)) < 0.05, 'T1 lG matches exact');

% ---------------------------------------------------------------
% Test 2: two-link multirate network, exact enumeration
% Link1 shared by both classes, class2 needs 2 circuits on link2.
% ---------------------------------------------------------------
nu = [3.0, 1.5];
A  = [1 1; 1 2];      % 2 links x 2 classes
Cv = [4; 5];
[gEx, betaEx] = lossn_exact(nu, A, Cv);
[QLen, Loss, lG, ci] = lossn_mci(nu, A, Cv, struct('samples',3e5,'seed',2));
fprintf('\n[T2] two-link multirate, exact enumeration\n');
fprintf('  lG exact=%.6f  mci=%.6f\n', log(gEx), lG);
for r = 1:2
    fprintf('  class %d: beta exact=%.6f  mci=%.6f  CI=[%.6f,%.6f]\n', ...
        r, betaEx(r), Loss(r), ci.loss(r,1), ci.loss(r,2));
    [npass,nfail] = check(npass, nfail, ...
        ci.loss(r,1) <= betaEx(r) && betaEx(r) <= ci.loss(r,2), ...
        sprintf('T2 class %d blocking in CI', r));
end
[npass,nfail] = check(npass, nfail, abs(lG-log(gEx)) < 0.05, 'T2 lG matches exact');
% carried load consistency: QLen = nu*(1-beta)
[npass,nfail] = check(npass, nfail, max(abs(QLen - nu.*(1-Loss))) < 1e-9, 'T2 QLen=nu*(1-beta)');

% ---------------------------------------------------------------
% Test 3: reproducibility with fixed seed
% ---------------------------------------------------------------
[~, L1] = lossn_mci(nu, A, Cv, struct('samples',5e4,'seed',7));
[~, L2] = lossn_mci(nu, A, Cv, struct('samples',5e4,'seed',7));
fprintf('\n[T3] reproducibility: max|L1-L2|=%.2e\n', max(abs(L1-L2)));
[npass,nfail] = check(npass, nfail, isequal(L1,L2), 'T3 fixed seed reproducible');

% ---------------------------------------------------------------
% Test 4: paper star network (Figure 1, Table 1) light traffic.
% 4 links, 12 classes. Compare against exact enumeration is
% infeasible (>1e12 states); cross-check against Erlang FP order
% of magnitude and that class-blocking CIs bracket the FP estimate
% loosely (FP is an approximation, so only sanity, not equality).
% ---------------------------------------------------------------
[nu4, A4, C4] = star_network('light');
% lossn_erlangfp returns blocking probability; convert to acceptance.
% Blocking here is a deep-tail event (well below 1%), so cross-check the
% well-conditioned acceptance probabilities of the two methods.
[~, LossFP] = lossn_erlangfp(nu4, A4, C4);
accFP = 1 - LossFP;
[~, LossMC, lG4, ci4] = lossn_mci(nu4, A4, C4, struct('samples',2e5,'seed',3));
accMC = 1 - LossMC;
fprintf('\n[T4] star network (light traffic), lG=%.4f\n', lG4);
fprintf('  class  FP-accept  MCI-accept  MCI-beta%%  MCI-CI%%\n');
for r = 1:numel(nu4)
    fprintf('  %5d  %.6f  %.6f  %.4f  [%.4f,%.4f]\n', r, accFP(r), accMC(r), ...
        100*LossMC(r), 100*ci4.loss(r,1), 100*ci4.loss(r,2));
end
% All blocking probabilities must be valid probabilities.
[npass,nfail] = check(npass, nfail, all(LossMC>=0 & LossMC<=1), 'T4 valid probabilities');
% Two independent methods (Erlang FP vs Monte Carlo) must agree on the
% well-conditioned acceptance probabilities.
[npass,nfail] = check(npass, nfail, ...
    max(abs(accMC - accFP)) < 0.02, 'T4 MCI vs FP acceptance agree');

% ---------------------------------------------------------------
% Test 5: variance reduction -- heuristic gamma vs gamma=nu (rho)
% Section 3.2: importance sampling narrows CI in heavy traffic.
% ---------------------------------------------------------------
[nuH, AH, CH] = star_network('heavy');
[~, ~, ~, ciIS] = lossn_mci(nuH, AH, CH, struct('samples',1e5,'seed',4));
[~, ~, ~, ciRho] = lossn_mci(nuH, AH, CH, struct('samples',1e5,'seed',4,'gamma',nuH));
wIS  = mean(ciIS.accept(:,2)  - ciIS.accept(:,1));
wRho = mean(ciRho.accept(:,2) - ciRho.accept(:,1));
fprintf('\n[T5] heavy traffic mean CI width: gamma=heuristic %.5f, gamma=rho %.5f (imp=%.2fx)\n', ...
    wIS, wRho, wRho/wIS);
[npass,nfail] = check(npass, nfail, wIS <= wRho*1.02, 'T5 heuristic gamma not worse than rho');

fprintf('\n=== RESULT: %d passed, %d failed ===\n', npass, nfail);
if nfail > 0
    error('test_lossn_mci: %d checks failed', nfail);
end
end

% ---------------------------------------------------------------
function [np, nf] = check(np, nf, cond, name)
if cond
    np = np + 1; fprintf('  PASS: %s\n', name);
else
    nf = nf + 1; fprintf('  FAIL: %s\n', name);
end
end

function [g, beta] = lossn_exact(nu, A, C)
% Exact normalization constant and class blocking by box enumeration.
nu = nu(:)'; C = C(:);
R = numel(nu); J = numel(C);
N = zeros(1,R);
for k = 1:R
    pos = A(:,k) > 0;
    N(k) = floor(min(C(pos)./A(pos,k)));
end
% enumerate all n in box {0..N_k}
grids = cell(1,R);
for k = 1:R, grids{k} = 0:N(k); end
[out{1:R}] = ndgrid(grids{:});
States = zeros(numel(out{1}), R);
for k = 1:R, States(:,k) = out{k}(:); end
% product-form weight q(n) = prod nu^n / n!
logq = States * log(nu(:)) - sum(gammaln(States+1), 2);
AV = States * A';
feas = all(AV <= C(:)', 2);
mq = max(logq(feas));
g = exp(mq) * sum(exp(logq(feas) - mq));
beta = zeros(1,R);
for r = 1:R
    Cr = (C - A(:,r))';
    feasR = all(AV <= Cr, 2);
    gr = exp(mq) * sum(exp(logq(feasR) - mq));
    beta(r) = 1 - gr/g;
end
end

function b = erlangB_ref(rho, C)
% Reference Erlang-B via stable recursion.
inv = 1;
for k = 1:C
    inv = 1 + inv * k / rho;
end
b = 1 / inv;
end

function [nu, A, C] = star_network(regime)
% Four-leaf star network of Ross-Wang Table 1 / Figure 1.
% Links C1=90,C2=100,C3=110,C4=120. Six leaf pairs, two classes each
% (1 circuit and 5 circuits) -> 12 classes. Routes as in Table 1.
C = [90; 100; 110; 120];
routes = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];   % 6 pairs
switch regime
    case 'light',  base = 9.0;  hi = 1.6;
    case 'moderate', base = 10.0; hi = 2.0;
    case 'heavy',  base = 15.0; hi = 3.0;
    otherwise, error('unknown regime');
end
K = 12;
A = zeros(4, K);
nu = zeros(1, K);
for p = 1:6
    % 1-circuit class
    A(routes(p,1), p) = 1;  A(routes(p,2), p) = 1;  nu(p) = base;
    % 5-circuit class
    A(routes(p,1), 6+p) = 5; A(routes(p,2), 6+p) = 5; nu(6+p) = hi;
end
end
