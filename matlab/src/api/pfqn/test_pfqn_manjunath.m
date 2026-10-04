function test_pfqn_manjunath()
% TEST_PFQN_MANJUNATH Validate the Manjunath-Sikdar transform for product-form
% queueing networks against pfqn_ca and against brute-force enumeration.
%
% THE TWO ORACLES DO NOT COME OUT OF THE IMPLEMENTATION. With no extra rows the
% transform computes the ordinary closed-network normalizing constant, for which
% pfqn_ca's convolution recursion is an independent exact algorithm sharing no
% code; the two are algebraic identities for the same sum, so agreement to 1e-13
% is the correct expectation and not a tolerance chosen to pass. With extra rows
% pfqn_ca has nothing to say, and the oracle becomes `bcmp_enum` below, which
% sums the BCMP product form over the enumerated state space and applies each
% row by direct comparison -- the very enumeration the transform exists to
% avoid, so a coefficient-domain defect cannot hide behind a shared traversal.

fprintf('=== test_pfqn_manjunath ===\n');
npass = 0; nfail = 0;

% ---------------------------------------------------------------
% Test 1: unconstrained networks reproduce pfqn_ca exactly
% ---------------------------------------------------------------
fprintf('\n[T1] unconstrained vs pfqn_ca\n');
cases = {
    struct('L', [1;2],                        'N', 4,     'Z', [])
    struct('L', [1 2; 3 1],                   'N', [2 3], 'Z', [])
    struct('L', [1 2; 3 1],                   'N', [2 3], 'Z', [0.5 1.5])
    struct('L', [0.4 0.2; 0.9 0.7; 0.1 1.1],  'N', [3 2], 'Z', [1 2])
    struct('L', zeros(0,2),                   'N', [2 1], 'Z', [1 3])
    struct('L', [1 2; 3 1],                   'N', [0 0], 'Z', [1 1])
    struct('L', [5 1; 1 6],                   'N', [6 5], 'Z', [2 3])
    };
for k = 1:numel(cases)
    c = cases{k};
    [~, lGca] = pfqn_ca(c.L, c.N, c.Z);
    [~, lGmj] = pfqn_manjunath(c.L, c.N, c.Z);
    e = abs(lGca - lGmj)/max(1, abs(lGca));
    fprintf('  case %d: lG_ca=%.14g lG_mj=%.14g relerr=%.3e\n', k, lGca, lGmj, e);
    [npass,nfail] = check(npass, nfail, e < 1e-13, ...
        sprintf('T1.%d unconstrained matches pfqn_ca', k));
end

% ---------------------------------------------------------------
% Test 2: extra rows against brute-force enumeration
% ---------------------------------------------------------------
L = [1 2; 3 1; 0.5 0.5]; N = [3 2]; Z = [1 2];
S = size(L,1) + size(Z,1);
a1 = zeros(1,S*2); a1(1) = 1; a1(1+S) = 1;              % jobs at queue 1
a2 = zeros(1,S*2); a2(1) = 2; a2(2) = 1; a2(1+S) = 1; a2(2+S) = 3;  % weighted budget
a3 = zeros(1,S*2); a3(4) = 1;                           % class 1 at the delay
tests = {
    struct('A', a1,          'b', 2,     'sense', 'L')
    struct('A', a1,          'b', 2,     'sense', 'E')
    struct('A', a1,          'b', 1,     'sense', 'G')
    struct('A', a2,          'b', 6,     'sense', 'L')
    struct('A', [a1; a2],    'b', [2;6], 'sense', 'LL')
    struct('A', [a1; a2],    'b', [2;6], 'sense', 'EL')
    struct('A', [a1; a2],    'b', [1;6], 'sense', 'GL')
    struct('A', [a1; a2],    'b', [1;5], 'sense', 'GG')
    struct('A', [a1; a3],    'b', [2;1], 'sense', 'LL')
    };
fprintf('\n[T2] constrained vs enumeration\n');
for k = 1:numel(tests)
    t = tests{k};
    gEn = bcmp_enum(L, N, Z, t.A, t.b, t.sense);
    [gMj, lGmj] = pfqn_manjunath(L, N, Z, t.A, t.b, t.sense);
    e = abs(gEn - gMj)/max(1e-300, abs(gEn));
    fprintf('  test %d (%s): enum=%.14g transform=%.14g relerr=%.3e\n', ...
        k, t.sense, gEn, gMj, e);
    [npass,nfail] = check(npass, nfail, e < 1e-12, ...
        sprintf('T2.%d row set %s matches enumeration', k, t.sense));
    [npass,nfail] = check(npass, nfail, abs(exp(lGmj) - gMj)/gMj < 1e-12, ...
        sprintf('T2.%d lG is the log of G', k));
end

% ---------------------------------------------------------------
% Test 3: the population constraint expressed as an extra row is redundant
% ---------------------------------------------------------------
% sum over every station and class of n_ir equals sum(N) in every admissible
% state, so declaring it changes nothing -- a direct check that an equality row
% is discharged by picking a coefficient rather than by summing one.
fprintf('\n[T3] redundant total-population row\n');
[gRef, ~] = pfqn_manjunath(L, N, Z);
[gEq, ~]  = pfqn_manjunath(L, N, Z, ones(1,S*2), sum(N), 'E');
[gLe, ~]  = pfqn_manjunath(L, N, Z, ones(1,S*2), sum(N), 'L');
[gGt, ~]  = pfqn_manjunath(L, N, Z, ones(1,S*2), sum(N), 'G');
fprintf('  G=%.14g  with ''='' %.14g  with ''<='' %.14g  with ''>'' %.14g\n', ...
    gRef, gEq, gLe, gGt);
[npass,nfail] = check(npass, nfail, abs(gEq-gRef)/gRef < 1e-13, 'T3 ''='' is redundant');
[npass,nfail] = check(npass, nfail, abs(gLe-gRef)/gRef < 1e-13, 'T3 ''<='' is redundant');
[npass,nfail] = check(npass, nfail, gGt == 0, 'T3 ''>'' is unsatisfiable');

% ---------------------------------------------------------------
% Test 4: trivial rows are decided rather than carried
% ---------------------------------------------------------------
fprintf('\n[T4] trivial rows\n');
z = zeros(1,S*2);
[npass,nfail] = check(npass, nfail, ...
    abs(pfqn_manjunath(L,N,Z,z,0,'E') - gRef)/gRef < 1e-13, 'T4 0=0 kept');
[npass,nfail] = check(npass, nfail, pfqn_manjunath(L,N,Z,z,3,'E') == 0, 'T4 0=3 empty');
[npass,nfail] = check(npass, nfail, pfqn_manjunath(L,N,Z,z,-1,'L') == 0, 'T4 sum<=-1 empty');
[npass,nfail] = check(npass, nfail, ...
    abs(pfqn_manjunath(L,N,Z,a1,-1,'G') - gRef)/gRef < 1e-13, 'T4 sum>-1 kept');

% ---------------------------------------------------------------
% Test 5: refusals name the reason
% ---------------------------------------------------------------
fprintf('\n[T5] refusals\n');
probes = {
    {'fractional A', [0.5 zeros(1,S*2-1)], 1,   'L'}
    {'negative A',   [-1  zeros(1,S*2-1)], 1,   'L'}
    {'fractional b', a1,                   1.5, 'L'}
    {'bad sense',    a1,                   1,   'X'}
    };
for k = 1:numel(probes)
    p = probes{k};
    caught = false;
    try
        pfqn_manjunath(L, N, Z, p{2}, p{3}, p{4});
    catch
        caught = true;
    end
    [npass,nfail] = check(npass, nfail, caught, sprintf('T5 %s refused', p{1}));
end

% ---------------------------------------------------------------
% Test 6: the per-class decomposition (fourth output)
% ---------------------------------------------------------------
% Reference instance: PS queue (demands 1, 2) + delay (think 2, 4), N = [4 4]
% both starting at the delay, with 2*n1 + 3*n2 <= 10 on the queue occupancy.
% The expectations are the stationary law of an INDEPENDENTLY built exact CTMC
% under HOLD truncation (a refused admission is a deleted transition), which
% agrees with the truncated product form to 8.2e-17, so these numbers are not
% the routine restating itself.
fprintf('\n[T6] per-class decomposition\n');
Ls = [1 2]; Zs = [2 4]; Ns = [4 4];
As = zeros(1,4); As(1) = 2; As(3) = 3;
[~,~,~,st] = pfqn_manjunath(Ls, Ns, Zs, As, 10, 'L');
ref = struct( ...
    'Q',       [1.88612099644128, 1.50177935943061], ...
    'X',       [0.544483985765125, 0.224199288256228], ...
    'U',       [0.544483985765125, 0.448398576512456], ...
    'think',   [1.08896797153025, 0.896797153024911], ...
    'blocked', [1.02491103202847, 1.60142348754448], ...
    'delay',   [2.11387900355872, 2.49822064056939]);
fn = fieldnames(ref);
for i = 1:numel(fn)
    got = st.(fn{i}); want = ref.(fn{i});
    e = max(abs(got-want)./abs(want));
    fprintf('  %-8s [%.12g %.12g] relerr=%.3e\n', fn{i}, got(1), got(2), e);
    [npass,nfail] = check(npass, nfail, e < 1e-12, ...
        sprintf('T6 %s matches the exact HOLD chain', fn{i}));
end
% A blocked job never leaves the delay, so nothing escapes the accounting.
[npass,nfail] = check(npass, nfail, ...
    max(abs(st.Q + st.think + st.blocked - Ns)) < 1e-12, 'T6 Q+think+blocked = N');
[npass,nfail] = check(npass, nfail, ...
    max(abs(st.delay - st.think - st.blocked)) < 1e-12, 'T6 delay = think + blocked');

% Unconstrained, the throughput identity collapses to the textbook ratio.
[~,~,~,s0] = pfqn_manjunath(Ls, Ns, Zs);
[~,lg0] = pfqn_ca(Ls, Ns, Zs);
[~,lg1] = pfqn_ca(Ls, [Ns(1)-1 Ns(2)], Zs);
[~,lg2] = pfqn_ca(Ls, [Ns(1) Ns(2)-1], Zs);
[npass,nfail] = check(npass, nfail, ...
    max(abs(s0.X - [exp(lg1-lg0) exp(lg2-lg0)])) < 1e-13, 'T6 X reduces to G(N-e_r)/G(N)');
[npass,nfail] = check(npass, nfail, max(abs(s0.blocked)) < 1e-12, ...
    'T6 nothing is held when nothing is constrained');

% ---------------------------------------------------------------
% Test 7: the decomposition refuses every other configuration
% ---------------------------------------------------------------
fprintf('\n[T7] decomposition refusals\n');
caught = false;
try, [~,~,~,~] = pfqn_manjunath(Ls, Ns, [Zs; Zs], zeros(1,6), 10, 'L'); catch, caught = true; end
[npass,nfail] = check(npass, nfail, caught, 'T7 two delay stations refused');
caught = false;
try, [~,~,~,~] = pfqn_manjunath([Ls; Ls], Ns, Zs, zeros(1,6), 10, 'L'); catch, caught = true; end
[npass,nfail] = check(npass, nfail, caught, 'T7 two queueing stations refused');
Abad = zeros(1,4); Abad(2) = 1;      % column 2 is (delay, class 1)
caught = false;
try, [~,~,~,~] = pfqn_manjunath(Ls, Ns, Zs, Abad, 10, 'L'); catch, caught = true; end
[npass,nfail] = check(npass, nfail, caught, 'T7 delay inside the region refused');

fprintf('\n=== test_pfqn_manjunath: %d passed, %d failed ===\n', npass, nfail);
if nfail > 0
    error('test_pfqn_manjunath: %d assertions failed', nfail);
end
end

% ------------------------------------------------------------------------

function g = bcmp_enum(L, N, Z, A, b, sense)
% Direct sum of the BCMP product form over the enumerated closed state space,
% keeping the states that satisfy every extra row. Independent of the transform,
% which never enumerates.
M = size(L,1); Mz = size(Z,1); S = M + Mz; R = numel(N);
alloc = cell(1,R);
for r = 1:R
    alloc{r} = compositions(N(r), S);
end
idx = ones(1,R);
g = 0;
while true
    n = zeros(S,R);
    for r = 1:R
        n(:,r) = alloc{r}(idx(r),:)';
    end
    ok = true;
    for j = 1:numel(b)
        v = A(j,:)*n(:);
        switch sense(j)
            case 'E', ok = ok && abs(v - b(j)) < 1e-9;
            case 'L', ok = ok && v <= b(j) + 1e-9;
            case 'G', ok = ok && v > b(j) + 1e-9;
        end
    end
    if ok
        t = 1;
        for i = 1:M
            t = t*factorial(sum(n(i,:)));
            for r = 1:R
                t = t*L(i,r)^n(i,r)/factorial(n(i,r));
            end
        end
        for kk = 1:Mz
            for r = 1:R
                t = t*Z(kk,r)^n(M+kk,r)/factorial(n(M+kk,r));
            end
        end
        g = g + t;
    end
    d = R;
    while d >= 1
        idx(d) = idx(d) + 1;
        if idx(d) <= size(alloc{d},1)
            break
        end
        idx(d) = 1;
        d = d - 1;
    end
    if d < 1
        break
    end
end
end

function C = compositions(n, k)
% Every way of splitting n indistinguishable jobs over k stations.
if k == 1
    C = n;
    return
end
C = [];
for a = 0:n
    sub = compositions(n - a, k - 1);
    C = [C; repmat(a, size(sub,1), 1), sub]; %#ok<AGROW>
end
end

function [npass, nfail] = check(npass, nfail, cond, name)
if cond
    npass = npass + 1;
    fprintf('  PASS: %s\n', name);
else
    nfail = nfail + 1;
    fprintf('  FAIL: %s\n', name);
end
end
