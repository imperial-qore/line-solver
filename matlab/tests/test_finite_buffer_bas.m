% MGS validation suite for the Smith Queue Decomposition (SQD) approximation
% (npfqn_sqd) on closed Blocking-After-Service (BAS) finite-buffer networks.
%
% Validates the first-queue throughput against published exact/simulation values and
% the MGS paper across Tables 1-11, with the paper-validation calibration
% (calibrationMode=0, serverBlockingTime=false, downstream/compound).
%
% Originally contributed as FiniteBufferBASTest by Avinash Bommareddy (Imperial College
% London FYP, 2026); ported to MATLAB against the npfqn_sqd API.
%
% NOTE: this test intentionally GATES on regression — the SQD approximation has
% known accuracy limits at high population N, so several cases fail by design
% (faithful to the original validation harness).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

GlobalConstants.setVerbose(VerboseLevel.SILENT);

ABS_TOL = 0.05;
REGRESS_MARGIN = 0.01;

s = struct('cases',0,'improved',0,'tied',0,'regressed',0,'failed',0, ...
           'sumOurs',0,'sumPaper',0,'worstErr',0,'worstDesc','', ...
           'absTol',ABS_TOL,'regressMargin',REGRESS_MARGIN);

% Table 1 — two-stage equal rates
exact = [0.500 0.667 0.750 0.800 0.833 0.800 0.750];
paper = [0.499 0.666 0.750 0.800 0.833 0.800 0.750];
for idx = 1:numel(exact)
    N = idx;
    got = solveSqd(buildCyclic('T1', [1 1], [4 4], N), N, 2);
    s = assertBoth(s, got, exact(idx), paper(idx), sprintf('Table 1 N=%d', N));
end

% Table 2 — two-stage unequal rates
exact = [1.333 1.714 1.867 1.936 1.968 1.968 1.968 1.936 1.867];
paper = [1.326 1.711 1.865 1.934 1.967 1.983 1.967 1.934 1.865];
for idx = 1:numel(exact)
    N = idx;
    got = solveSqd(buildCyclic('T2', [2 4], [4 6], N), N, 2);
    s = assertBoth(s, got, exact(idx), paper(idx), sprintf('Table 2 N=%d', N));
end

% Table 3 Akyildiz
exact = [0.250 0.308 0.325 0.331 0.331 0.331 0.325];
paper = [0.250 0.308 0.325 0.331 0.332 0.331 0.325];
for idx = 1:numel(exact)
    N = idx;
    got = solveSqd(buildCyclic('T3A', [1/3 1.0], [3 5], N), N, 2);
    s = assertBoth(s, got, exact(idx), paper(idx), sprintf('Table 3 Akyildiz N=%d', N));
end

% Table 3 Bolch
exact = [0.345 0.439 0.474 0.489 0.495 0.498 0.498 0.498 0.495 0.489 0.474];
paper = [0.344 0.439 0.474 0.488 0.495 0.498 0.499 0.498 0.495 0.488 0.474];
for idx = 1:numel(exact)
    N = idx;
    got = solveSqd(buildCyclic('T3B', [0.5 10/9], [7 5], N), N, 2);
    s = assertBoth(s, got, exact(idx), paper(idx), sprintf('Table 3 Bolch N=%d', N));
end

% Table 5 — three-stage split
routing = [0.0 0.50 0.50; 0.70 0.0 0.30; 0.70 0.30 0.0];
mu = [2/5 5/6 1.0]; K = [6 6 6];
sim   = [0.245 0.338 0.376 0.391 0.397 0.399 0.400 0.400 0.400 0.400 0.400];
paper = [0.245 0.338 0.376 0.391 0.396 0.399 0.400 0.400 0.400 0.400 0.400];
for idx = 1:numel(sim)
    N = idx;
    got = solveSqd(buildClosed('T5', mu, K, routing, N), N, 3);
    s = assertBoth(s, got, sim(idx), paper(idx), sprintf('Table 5 Split N=%d', N));
end

% Table 7 — cyclic experiments 1-14
s = assertBoth(s, solveSqd(buildCyclic('T7E1', [3 2 4 2], [6 2 2 4], 9), 9, 4), 1.606, 1.726, 'Table 7 Exp 1');
s = assertBoth(s, solveSqd(buildCyclic('T7E2', [2 1 4 2], [3 4 5 2], 9), 9, 4), 0.978, 0.993, 'Table 7 Exp 2');
s = assertBoth(s, solveSqd(buildCyclic('T7E3', [3 2 4 2 1], [4 3 2 4 2], 10), 10, 5), 0.931, 0.994, 'Table 7 Exp 3');
s = assertBoth(s, solveSqd(buildCyclic('T7E4', [1 1 1 3 2 3], [2 2 2 2 2 2], 7), 7, 6), 0.668, 0.735, 'Table 7 Exp 4');
s = assertBoth(s, solveSqd(buildCyclic('T7E5', [2 1 4 3 1 4], [2 2 2 2 2 2], 7), 7, 6), 0.817, 0.832, 'Table 7 Exp 5');
s = assertBoth(s, solveSqd(buildCyclic('T7E6', [3 2 4 5 1 2 3], [2 2 2 2 2 2 2], 9), 9, 7), 0.925, 0.987, 'Table 7 Exp 6');
s = assertBoth(s, solveSqd(buildCyclic('T7E7', [4 2 2 3 5 2 3], [3 2 3 3 2 2 2], 10), 10, 7), 1.460, 1.576, 'Table 7 Exp 7');
s = assertBoth(s, solveSqd(buildCyclic('T7E8', [3 1 2 1 2 3 4], [3 2 3 3 2 2 2], 10), 10, 7), 0.827, 0.871, 'Table 7 Exp 8');
s = assertBoth(s, solveSqd(buildCyclic('T7E9', [1 2 2 1], [4 2 6 2], 8), 8, 4), 0.805, 0.859, 'Table 7 Exp 9');
s = assertBoth(s, solveSqd(buildCyclic('T7E10', [1 4 3 2], [3 2 6 2], 8), 8, 4), 0.959, 0.998, 'Table 7 Exp 10');
s = assertBoth(s, solveSqd(buildCyclic('T7E11', [3 4 4 1], [5 6 2 4], 8), 8, 4), 0.998, 0.999, 'Table 7 Exp 11');
s = assertBoth(s, solveSqd(buildCyclic('T7E12', [1 0.5 2 0.75 1], [3 2 3 3 2], 7), 7, 5), 0.450, 0.464, 'Table 7 Exp 12');
s = assertBoth(s, solveSqd(buildCyclic('T7E13', [2 0.5 1 0.75 1 1.5], [2 3 2 3 3 2], 10), 10, 6), 0.454, 0.485, 'Table 7 Exp 13');
s = assertBoth(s, solveSqd(buildCyclic('T7E14', [1 2 1 2 1 2 1], [3 4 3 4 2 2 3], 13), 13, 7), 0.730, 0.745, 'Table 7 Exp 14');

% Table 8 — eight-stage balanced
mu = [2 2 2 2 2 2 2 2]; K = [4 4 4 4 4 4 4 4];
Ns = [10 20 30]; sim = [1.175 1.380 1.037]; paper = [1.176 1.439 1.066];
for idx = 1:numel(Ns)
    N = Ns(idx);
    got = solveSqd(buildCyclic('T8', mu, K, N), N, 8);
    s = assertBoth(s, got, sim(idx), paper(idx), sprintf('Table 8 N=%d', N));
end

% Table 9 — eight-stage unbalanced
mu = [2 8 5 2.5 2 4 1.25 5]; K = [5 2 3 5 4 3 7 3];
Ns = [10 20 30]; sim = [1.194 1.232 1.072]; paper = [1.196 1.245 1.144];
for idx = 1:numel(Ns)
    N = Ns(idx);
    got = solveSqd(buildCyclic('T9', mu, K, N), N, 8);
    s = assertBoth(s, got, sim(idx), paper(idx), sprintf('Table 9 N=%d', N));
end

% Table 10 — five-stage split-merge
mu = [4.0 2.5 2.0 1.0 2.5]; K = [6 2 4 5 3];
routing = [0.0 0.20 0.30 0.20 0.30; 1.0 0 0 0 0; 1.0 0 0 0 0; 1.0 0 0 0 0; 1.0 0 0 0 0];
sim   = [1.250 2.038 2.566 2.923 3.171 3.339 3.438];
paper = [1.244 2.030 2.556 2.922 3.185 3.377 3.519];
for idx = 1:numel(sim)
    N = idx;
    got = solveSqd(buildClosed('T10', mu, K, routing, N), N, 5);
    s = assertBoth(s, got, sim(idx), paper(idx), sprintf('Table 10 N=%d', N));
end

% Table 11 — ten-stage split-merge
Ns = [1 5 10 15 20 25 30];
sim   = [0.800 2.837 4.112 4.797 5.140 5.284 5.287];
paper = [0.790 2.811 4.102 4.817 5.254 5.538 5.698];
for idx = 1:numel(Ns)
    N = Ns(idx);
    got = solveSqd(buildTable11(N), N, 10);
    s = assertBoth(s, got, sim(idx), paper(idx), sprintf('Table 11 N=%d', N));
end

% Visit-ratio exactness check (Table 4)
model = buildClosed('T4_visitcheck', [2/5 5/6 1.0], [6 6 6], routing3(), 5);
sn = model.getStruct(false);
[~,~,Vchain] = sn_get_demands_chain(sn);
v1 = Vchain(1,1); v2 = Vchain(2,1); v3 = Vchain(3,1);
fprintf('ratios: V2/V1=%.6f  V3/V1=%.6f  (expect 0.714286)\n', v2/v1, v3/v1);
assert(abs(v2/v1 - 5/7) < 1e-6, 'visit ratio V2/V1 mismatch');
assert(abs(v3/v1 - 5/7) < 1e-6, 'visit ratio V3/V1 mismatch');

% Summary + regression gate
fprintf('\n──────────────────────── MGS evaluation summary ────────────────────────\n');
fprintf('cases=%d | beat paper=%d  tied=%d  regressed=%d | hard-fails(regressions)=%d\n', ...
        s.cases, s.improved, s.tied, s.regressed, s.failed);
if s.cases > 0
    fprintf('mean abs error vs ground truth:  ours=%.2f%%   paper=%.2f%%   (delta=%+.2f pp)\n', ...
            s.sumOurs/s.cases, s.sumPaper/s.cases, (s.sumPaper - s.sumOurs)/s.cases);
    fprintf('worst case: %s (errOurs=%.1f%%)\n', s.worstDesc, s.worstErr);
end
fprintf('─────────────────────────────────────────────────────────────────────────\n');

assert(s.failed == 0, sprintf('%d case(s) regressed beyond the %.1fpp margin — see the per-case log above.', ...
        s.failed, REGRESS_MARGIN*100));

% =========================================================================
% Local functions
% =========================================================================

function r = routing3()
    r = [0.0 0.50 0.50; 0.70 0.0 0.30; 0.70 0.30 0.0];
end

function got = solveSqd(model, N, numQueues)
    X = npfqn_sqd(model.getStruct(false), N, 0, false, 'downstream', 'compound', []);
    assert(numel(X) == numQueues, 'unexpected station count');
    got = X(1);
end

function model = buildClosed(name, mu, K, qRouting, N)
    model = Network(name);
    M = numel(mu);
    queues = cell(1, M);
    for i = 1:M
        queues{i} = Queue(model, sprintf('Q%d', i), SchedStrategy.FCFS);
    end
    jobs = ClosedClass(model, 'Jobs', N, queues{1});
    for i = 1:M
        queues{i}.setService(jobs, Exp(mu(i)));
        queues{i}.setCapacity(K(i));
        queues{i}.setDropRule(jobs, DropStrategy.BAS);
    end
    % Closed network with only the M queue nodes, so node order == queue order
    % and the station routing matrix maps directly onto the node routing matrix.
    P = model.initRoutingMatrix;
    P{jobs, jobs} = qRouting;
    model.link(P);
end

function model = buildCyclic(name, mu, K, N)
    M = numel(mu);
    r = zeros(M, M);
    for i = 1:M
        r(i, mod(i, M) + 1) = 1.0;
    end
    model = buildClosed(name, mu, K, r, N);
end

function model = buildTable11(N)
    mu = [8.0 2.0 2.0 2.5 2.5 4.0 2.5 10.0 8.0 10.0];
    K  = [6 7 7 6 5 4 8 6 6 5];
    r = zeros(10, 10);
    r(1,2)=0.30; r(1,4)=0.30; r(1,6)=0.40;
    r(2,3)=1.0;  r(3,8)=1.0;
    r(4,5)=1.0;  r(5,8)=1.0;
    r(6,7)=1.0;  r(7,8)=1.0;
    r(8,9)=1.0;  r(9,10)=1.0;
    r(10,1)=1.0;
    model = buildClosed(sprintf('T11_N%d', N), mu, K, r, N);
end

function s = assertBoth(s, got, exact, paperMGS, desc)
    errOurs  = abs(got - exact)    / exact    * 100.0;
    errPaper = abs(paperMGS - exact) / exact   * 100.0;
    delta    = errPaper - errOurs;

    accurate  = errOurs <= s.absTol * 100.0;
    noRegress = errOurs <= errPaper + s.regressMargin * 100.0;
    pass      = accurate || noRegress;

    s.cases = s.cases + 1;
    s.sumOurs  = s.sumOurs + errOurs;
    s.sumPaper = s.sumPaper + errPaper;
    if delta > 0.05
        s.improved = s.improved + 1;
    elseif delta < -0.05
        s.regressed = s.regressed + 1;
    else
        s.tied = s.tied + 1;
    end
    if errOurs > s.worstErr
        s.worstErr = errOurs; s.worstDesc = desc;
    end
    if ~pass
        s.failed = s.failed + 1;
    end

    if delta > 0.05
        vsPaper = sprintf('%+.1fpp BEAT', delta);
    elseif delta < -0.05
        vsPaper = sprintf('%+.1fpp REGRESS', delta);
    else
        vsPaper = '~tied';
    end
    if accurate
        verdict = 'PASS (accurate)';
    elseif noRegress
        verdict = 'PASS (no regression)';
    else
        verdict = 'FAIL (regression vs paper)';
    end

    fprintf('[MGS] %-22s | exact=%6.4f | got=%6.4f (errOurs=%5.1f%%) | paperMGS=%6.4f (errPaper=%5.1f%%) | %-14s | %s\n', ...
            desc, exact, got, errOurs, paperMGS, errPaper, vsPaper, verdict);
end
