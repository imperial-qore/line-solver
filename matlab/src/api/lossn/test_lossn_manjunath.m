function test_lossn_manjunath()
% TEST_LOSSN_MS Validate the Manjunath-Sikdar transform for loss networks
% against Erlang-B, exact box enumeration, and the Monte Carlo estimator.

fprintf('=== test_lossn_manjunath ===\n');
rng(20260724);
npass = 0; nfail = 0;

% ---------------------------------------------------------------
% Test 1: single-link Erlang-B (J=1, one class, unit demand)
% ---------------------------------------------------------------
rho = 8; Ccap = 10;
[QLen, Loss, lG] = lossn_manjunath(rho, 1, Ccap);
exactB = erlangB_ref(rho, Ccap);
lGref = log(sum(exp((0:Ccap)*log(rho) - gammaln((0:Ccap)+1))));
fprintf('\n[T1] M/M/C/C Erlang-B  rho=%g C=%d\n', rho, Ccap);
fprintf('  blocking analytic=%.12f  lossn_manjunath=%.12f\n', exactB, Loss);
fprintf('  lG analytic=%.12f  lossn_manjunath=%.12f\n', lGref, lG);
[npass,nfail] = check(npass, nfail, abs(Loss-exactB) < 1e-12, 'T1 blocking matches Erlang-B');
[npass,nfail] = check(npass, nfail, abs(lG-lGref) < 1e-12, 'T1 lG matches analytic');
[npass,nfail] = check(npass, nfail, abs(QLen - rho*(1-Loss)) < 1e-12, 'T1 QLen=nu*(1-beta)');

% ---------------------------------------------------------------
% Test 2: two-link multirate network against box enumeration
% ---------------------------------------------------------------
nu = [3.0, 1.5];
A  = [1 1; 1 2];
Cv = [4; 5];
[QLen, Loss, lG] = lossn_manjunath(nu, A, Cv);
[gEn, betaEn] = lossn_enum(nu, A, Cv);
fprintf('\n[T2] two-link multirate\n');
fprintf('  lG transform=%.12f  enumeration=%.12f\n', lG, log(gEn));
for r = 1:2
    fprintf('  class %d: beta transform=%.12f  enumeration=%.12f\n', r, Loss(r), betaEn(r));
end
[npass,nfail] = check(npass, nfail, max(abs(Loss-betaEn)) < 1e-10, 'T2 blocking matches enumeration');
[npass,nfail] = check(npass, nfail, abs(lG-log(gEn)) < 1e-10, 'T2 lG matches enumeration');
[npass,nfail] = check(npass, nfail, max(abs(QLen - nu.*(1-Loss))) < 1e-12, 'T2 QLen=nu*(1-beta)');

% ---------------------------------------------------------------
% Test 3: memory-style and general linear constraints, three classes.
% Row 1 is a job count, row 2 a weighted memory budget, rows 3-4 general
% linear admission constraints. This is the constraint shape a LINE finite
% capacity region produces from setGlobalMaxJobs / setGlobalMaxMemory with
% setClassSize / setConstraint.
% ---------------------------------------------------------------
nu = [4.0, 2.5, 2.0];
A  = [1 1 1;      % global job cap
      1 2 3;      % memory, class sizes 1/2/3
      2 1 1;      % linear constraint 1
      1 3 2];     % linear constraint 2
Cv = [9; 14; 10; 12];
[QLen, Loss, lG] = lossn_manjunath(nu, A, Cv);
[gEn, betaEn] = lossn_enum(nu, A, Cv);
fprintf('\n[T3] four rows: jobs, memory, two linear constraints\n');
fprintf('  lG transform=%.12f  enumeration=%.12f\n', lG, log(gEn));
for r = 1:3
    fprintf('  class %d: beta transform=%.12f  enumeration=%.12f\n', r, Loss(r), betaEn(r));
end
[npass,nfail] = check(npass, nfail, max(abs(Loss-betaEn)) < 1e-10, 'T3 blocking matches enumeration');
[npass,nfail] = check(npass, nfail, abs(lG-log(gEn)) < 1e-10, 'T3 lG matches enumeration');

% ---------------------------------------------------------------
% Test 4: gcd row reduction is exact. Scaling a row and its capacity by a
% common integer factor must not change the answer.
% ---------------------------------------------------------------
nu = [2.0, 1.0];
[~, L1] = lossn_manjunath(nu, [2 4], 9);
[~, L2] = lossn_manjunath(nu, [1 2], 4);      % 9/2 floors to 4
fprintf('\n[T4] gcd reduction: max|L1-L2|=%.2e\n', max(abs(L1-L2)));
[npass,nfail] = check(npass, nfail, max(abs(L1-L2)) < 1e-12, 'T4 gcd reduction exact');

% ---------------------------------------------------------------
% Test 5: a class appearing in no constraint never blocks and contributes
% exp(nu_r) to the normalization constant.
% ---------------------------------------------------------------
nu = [2.0, 1.0, 1.5];
A  = [1 1 0];
Cv = 6;
[QLen, Loss, lG] = lossn_manjunath(nu, A, Cv);
[~, L2c] = lossn_manjunath(nu(1:2), A(1,1:2), Cv);
fprintf('\n[T5] unconstrained class: beta_3=%.3g, QLen_3=%.6f (nu_3=%.6f)\n', Loss(3), QLen(3), nu(3));
[npass,nfail] = check(npass, nfail, Loss(3) == 0, 'T5 unconstrained class never blocks');
[npass,nfail] = check(npass, nfail, abs(QLen(3)-nu(3)) < 1e-12, 'T5 unconstrained QLen=nu');
[npass,nfail] = check(npass, nfail, max(abs(Loss(1:2)-L2c)) < 1e-12, 'T5 free class does not perturb others');
gEn = lossn_enum(nu(1:2), A(1,1:2), Cv);
[npass,nfail] = check(npass, nfail, abs(lG - (log(gEn)+nu(3))) < 1e-10, 'T5 lG carries exp(nu) factor');

% ---------------------------------------------------------------
% Test 6: agreement with the Monte Carlo estimator, which is independent of
% the transform. The exact blocking must lie inside the reported interval.
% ---------------------------------------------------------------
nu = [3.0, 1.5];
A  = [1 1; 1 2];
Cv = [4; 5];
[~, Lms] = lossn_manjunath(nu, A, Cv);
[~, ~, ~, ci] = lossn_mci(nu, A, Cv, struct('samples',3e5,'seed',5));
fprintf('\n[T6] transform inside the Monte Carlo interval\n');
for r = 1:2
    fprintf('  class %d: exact=%.6f  CI=[%.6f,%.6f]\n', r, Lms(r), ci.loss(r,1), ci.loss(r,2));
    [npass,nfail] = check(npass, nfail, ...
        ci.loss(r,1) <= Lms(r) && Lms(r) <= ci.loss(r,2), ...
        sprintf('T6 class %d exact inside CI', r));
end

% ---------------------------------------------------------------
% Test 7: rare blocking. The transform is exact where a sampled or simulated
% throughput cannot resolve the tail.
% ---------------------------------------------------------------
nu = [1.0, 0.5];
[~, Lrare, ~] = lossn_manjunath(nu, [1 1], 8);
k = 0:8;
p = exp(k*log(1.5) - gammaln(k+1)); p = p/sum(p);
fprintf('\n[T7] rare blocking: transform=%.6e  analytic=%.6e\n', Lrare(1), p(end));
[npass,nfail] = check(npass, nfail, abs(Lrare(1)-p(end)) < 1e-14, 'T7 rare blocking exact');

% ---------------------------------------------------------------
% Test 8: non-integer input is rejected rather than silently rounded.
% ---------------------------------------------------------------
threw = false;
try
    lossn_manjunath([1.0 1.0], [1 1.5], 5);
catch
    threw = true;
end
fprintf('\n[T8] fractional requirement rejected: %d\n', threw);
[npass,nfail] = check(npass, nfail, threw, 'T8 fractional A rejected');

% ---------------------------------------------------------------
% Test 9: a heavy route does not overflow the terms. nu^n/n! peaks near
% exp(nu)/sqrt(2 pi nu), so forming the terms and dividing by their maximum
% afterwards returns Inf above a load of about 700, and Inf is not finite so
% the rescaling then declined to run and NaN reached every metric. The terms
% are built in the log domain and shifted before exponentiating instead.
% ---------------------------------------------------------------
nuBig = 900; Cbig = 900;
[QLen, Loss, lG] = lossn_manjunath(nuBig, 1, Cbig);
bBig = erlangB_ref(nuBig, Cbig);
fprintf('\n[T9] heavy load nu=C=%d: blocking=%.15g  Erlang-B=%.15g  lG=%.15g\n', ...
    nuBig, Loss, bBig, lG);
[npass,nfail] = check(npass, nfail, isfinite(Loss) && isfinite(lG), 'T9 heavy load stays finite');
[npass,nfail] = check(npass, nfail, abs(Loss - bBig) < 1e-12, 'T9 heavy load matches Erlang-B');
[npass,nfail] = check(npass, nfail, abs(QLen - nuBig*(1-Loss)) < 1e-9, 'T9 heavy QLen=nu*(1-beta)');

fprintf('\n=== test_lossn_manjunath: %d passed, %d failed ===\n', npass, nfail);
if nfail > 0
    error('test_lossn_manjunath: %d assertion(s) failed', nfail);
end
end

% ------------------------------------------------------------------------

function [np, nf] = check(np, nf, cond, name)
if cond
    np = np + 1;
    fprintf('  [PASS] %s\n', name);
else
    nf = nf + 1;
    fprintf('  [FAIL] %s\n', name);
end
end

function [g, beta] = lossn_enum(nu, A, C)
% Exact normalization constant and class blocking by box enumeration.
nu = nu(:)'; C = C(:);
R = numel(nu);
N = zeros(1,R);
for k = 1:R
    pos = A(:,k) > 0;
    N(k) = floor(min(C(pos)./A(pos,k)));
end
grids = cell(1,R);
for k = 1:R, grids{k} = 0:N(k); end
[out{1:R}] = ndgrid(grids{:});
States = zeros(numel(out{1}), R);
for k = 1:R, States(:,k) = out{k}(:); end
logq = States * log(nu(:)) - sum(gammaln(States+1), 2);
AV = States * A';
feas = all(AV <= C(:)', 2);
mq = max(logq(feas));
g = exp(mq) * sum(exp(logq(feas) - mq));
beta = zeros(1,R);
for r = 1:R
    Cr = (C - A(:,r))';
    if any(Cr < 0)
        beta(r) = 1;
        continue
    end
    feasR = all(AV <= Cr, 2);
    gr = exp(mq) * sum(exp(logq(feasR) - mq));
    beta(r) = 1 - gr/g;
end
end

function b = erlangB_ref(rho, C)
inv = 1;
for k = 1:C
    inv = 1 + inv*k/rho;
end
b = 1/inv;
end
