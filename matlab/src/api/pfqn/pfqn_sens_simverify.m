function pfqn_sens_simverify()
%{
%{
 % @file pfqn_sens_simverify.m
 % @brief End-to-end verification of the queue-length second moments produced by
 %        pfqn_sens_mva and pfqn_sens_mvaldmx against discrete-event simulation.
 %
 %        pfqn_sens_mva_validate and pfqn_sens_mvaldmx_validate check the two
 %        primitives against brute-force enumeration of the product form, i.e.
 %        against the model the formulas assume. This harness instead checks
 %        them against LDES and JMT, i.e. against a simulated sample path of an
 %        actual LINE Network. It therefore closes the remaining gap: that the
 %        product form assumed by the formulas is the process the simulators
 %        realize, and that the demand matrix handed to the API is the one the
 %        model expresses.
 %
 %        Two independent simulated references are used:
 %          - LDES, through the Markov reward engine. setReward defines the
 %            nonlinear rewards n(i,r) and n(i,r)*n(i,s), which getAvgReward
 %            evaluates on the exact joint-state residence-time histogram
 %            exported by the engine. This is exact given the histogram, so its
 %            only error is the simulation's own statistical error.
 %          - JMT, through the logged sample path. sampleAggr returns the jump
 %            times and the per-class queue lengths, from which the
 %            time-weighted covariance is formed directly.
 %
 %        Three models exercise the three code paths:
 %          1. closed multiclass, load independent   -> pfqn_sens_mva
 %          2. closed, multiserver station           -> pfqn_sens_mvaldmx (LD)
 %          3. mixed open and closed                 -> pfqn_sens_mvaldmx (mixed)
 %
 %        Tolerances are statistical, not numerical. Variances converge faster
 %        than covariances, which are small differences of larger numbers, so
 %        they are checked against an absolute-plus-relative band rather than a
 %        tight relative one.
%}
%}
samples = 500000;
tolVar = 0.03;    % relative band on variances
tolCov = 0.05;    % absolute-plus-relative band on covariances
tolMean = 0.02;   % relative band on the means, a sanity check on the mapping

fprintf('\n=== pfqn_sens_mva / pfqn_sens_mvaldmx verification against simulation ===\n');
fprintf('    samples per run: %d\n', samples);
failures = {};

% =====================================================================
% Model 1: closed multiclass, load independent
% =====================================================================
fprintf('\n--- Model 1: closed multiclass, load independent (pfqn_sens_mva) ---\n');
model = Network('m1');
delay = Delay(model,'Think');
q1 = Queue(model,'Q1',SchedStrategy.PS);
q2 = Queue(model,'Q2',SchedStrategy.PS);
c1 = ClosedClass(model,'C1',3,delay,0);
c2 = ClosedClass(model,'C2',2,delay,0);
delay.setService(c1,Exp(1/1.0));  delay.setService(c2,Exp(1/0.5));
q1.setService(c1,Exp(1/0.4));     q1.setService(c2,Exp(1/0.6));
q2.setService(c1,Exp(1/0.3));     q2.setService(c2,Exp(1/0.2));
P = model.initRoutingMatrix;
P{c1} = Network.serialRouting(delay,q1,q2);
P{c2} = Network.serialRouting(delay,q1,q2);
model.link(P);

D1 = [0.4 0.6; 0.3 0.2];
ref1 = pfqn_sens_mva(D1,[3 2],[1.0 0.5]);
failures = [failures, compare_station('M1/Q1 (LDES)', ref1, ldes_moments(model,{q1,q2},{c1,c2},samples,1), 1, tolMean, tolVar, tolCov)];
failures = [failures, compare_station('M1/Q2 (LDES)', ref1, ldes_moments(model,{q1,q2},{c1,c2},samples,2), 2, tolMean, tolVar, tolCov)];
failures = [failures, compare_station('M1/Q1 (JMT)',  ref1, jmt_moments(model,q1,samples), 1, tolMean, tolVar, tolCov)];
failures = [failures, compare_station('M1/Q2 (JMT)',  ref1, jmt_moments(model,q2,samples), 2, tolMean, tolVar, tolCov)];

% =====================================================================
% Model 2: closed, multiserver station (load dependent)
% =====================================================================
fprintf('\n--- Model 2: closed, multiserver station (pfqn_sens_mvaldmx, LD) ---\n');
nserv = 2;
model2 = Network('m2');
delay2 = Delay(model2,'Think');
m1q = Queue(model2,'Q1',SchedStrategy.PS);
m2q = Queue(model2,'Q2',SchedStrategy.PS);
m1q.setNumberOfServers(nserv);
k1 = ClosedClass(model2,'C1',4,delay2,0);
delay2.setService(k1,Exp(1/1.0));
m1q.setService(k1,Exp(1/0.4));
m2q.setService(k1,Exp(1/0.3));
model2.link(Network.serialRouting(delay2,m1q,m2q));

Nc2 = 4;
D2 = [0.4; 0.3];
mu2 = zeros(2,Nc2);
for n = 1:Nc2
    mu2(1,n) = min(n,nserv);   % multiserver station: rate min(n,c)
    mu2(2,n) = 1;              % single server
end
ref2 = pfqn_sens_mvaldmx(0,D2,Nc2,1.0,mu2,[nserv;1]);
failures = [failures, compare_station('M2/Q1 (LDES)', ref2, ldes_moments(model2,{m1q,m2q},{k1},samples,1), 1, tolMean, tolVar, tolCov)];
failures = [failures, compare_station('M2/Q2 (LDES)', ref2, ldes_moments(model2,{m1q,m2q},{k1},samples,2), 2, tolMean, tolVar, tolCov)];
failures = [failures, compare_station('M2/Q1 (JMT)',  ref2, jmt_moments(model2,m1q,samples), 1, tolMean, tolVar, tolCov)];
failures = [failures, compare_station('M2/Q2 (JMT)',  ref2, jmt_moments(model2,m2q,samples), 2, tolMean, tolVar, tolCov)];

% =====================================================================
% Model 3: mixed open and closed
% =====================================================================
fprintf('\n--- Model 3: mixed open and closed (pfqn_sens_mvaldmx, mixed) ---\n');
lambdaOpen = 0.3;
model3 = Network('m3');
delay3 = Delay(model3,'Think');
src = Source(model3,'Source');
snk = Sink(model3,'Sink');
x1 = Queue(model3,'Q1',SchedStrategy.PS);
x2 = Queue(model3,'Q2',SchedStrategy.PS);
cc = ClosedClass(model3,'C',2,delay3,0);
oc = OpenClass(model3,'O',0);
delay3.setService(cc,Exp(1/1.0));
x1.setService(cc,Exp(1/0.4));   x2.setService(cc,Exp(1/0.3));
src.setArrival(oc,Exp(lambdaOpen));
x1.setService(oc,Exp(1/0.5));   x2.setService(oc,Exp(1/0.4));
P3 = model3.initRoutingMatrix;
P3{cc} = Network.serialRouting(delay3,x1,x2);
P3{oc} = Network.serialRouting(src,x1,x2,snk);
model3.link(P3);

% class order in the API call: 1 = closed, 2 = open
D3 = [0.4 0.5; 0.3 0.4];
N3 = [2 Inf];
lam3 = [0 lambdaOpen];
Z3 = [1.0 0];
mu3 = ones(2,2);
ref3 = pfqn_sens_mvaldmx(lam3,D3,N3,Z3,mu3,[1;1]);
failures = [failures, compare_station('M3/Q1 (LDES)', ref3, ldes_moments(model3,{x1,x2},{cc,oc},samples,1), 1, tolMean, tolVar, tolCov)];
failures = [failures, compare_station('M3/Q2 (LDES)', ref3, ldes_moments(model3,{x1,x2},{cc,oc},samples,2), 2, tolMean, tolVar, tolCov)];
failures = [failures, compare_station('M3/Q1 (JMT)',  ref3, jmt_moments(model3,x1,samples), 1, tolMean, tolVar, tolCov)];
failures = [failures, compare_station('M3/Q2 (JMT)',  ref3, jmt_moments(model3,x2,samples), 2, tolMean, tolVar, tolCov)];

% =====================================================================
% Model 4: higher moments of the station totals (pfqn_sens_mom)
% =====================================================================
fprintf('\n--- Model 4: higher moments of station totals (pfqn_sens_mom) ---\n');
% see _kb/03-api-layer.md (pfqn_sens_simverify -- FCFS sojourn moments vs JMT samples)
mom1 = pfqn_sens_mom(D1,[3 2],[1.0 0.5]);
for ist = 1:2
    if ist == 1, node = q1; else, node = q2; end
    sim = ldes_total_moments(model, node, {c1,c2}, samples);
    fprintf('  M4/Q%d (LDES)   E[Q]    analytic %8.5f  sim %8.5f  rel %7.4f\n', ist, mom1.m(ist), sim.m1, relb(mom1.m(ist),sim.m1));
    fprintf('  M4/Q%d (LDES)   E[Q^2]  analytic %8.5f  sim %8.5f  rel %7.4f\n', ist, mom1.M2(ist), sim.m2, relb(mom1.M2(ist),sim.m2));
    fprintf('  M4/Q%d (LDES)   E[Q^3]  analytic %8.5f  sim %8.5f  rel %7.4f\n', ist, mom1.M3(ist), sim.m3, relb(mom1.M3(ist),sim.m3));
    if relb(mom1.m(ist),sim.m1) > tolMean
        failures{end+1} = sprintf('M4/Q%d E[Q]: analytic %g sim %g', ist, mom1.m(ist), sim.m1); %#ok<AGROW>
    end
    if relb(mom1.M2(ist),sim.m2) > tolVar
        failures{end+1} = sprintf('M4/Q%d E[Q^2]: analytic %g sim %g', ist, mom1.M2(ist), sim.m2); %#ok<AGROW>
    end
    if relb(mom1.M3(ist),sim.m3) > tolVar
        failures{end+1} = sprintf('M4/Q%d E[Q^3]: analytic %g sim %g', ist, mom1.M3(ist), sim.m3); %#ok<AGROW>
    end
end

% =====================================================================
% Model 5: FCFS sojourn-time moments (pfqn_sens_respt)
% =====================================================================
% see _kb/03-api-layer.md (pfqn_sens_simverify -- FCFS sojourn moments vs JMT samples)
fprintf('\n--- Model 5: FCFS sojourn-time moments (pfqn_sens_respt) ---\n');
model5 = Network('m5');
d5  = Delay(model5,'Think');
f1 = Queue(model5,'F1',SchedStrategy.FCFS);
f2 = Queue(model5,'F2',SchedStrategy.FCFS);
f2.setNumberOfServers(2);
k5 = ClosedClass(model5,'C',3,d5,0);
d5.setService(k5,Exp(1/1.0));
f1.setService(k5,Exp(1/0.4));
f2.setService(k5,Exp(1/0.3));
model5.link(Network.serialRouting(d5,f1,f2));

ref5 = pfqn_sens_respt([0.4;0.3],[1;1],3,1.0,[1;2],3);
% RD from getCdfRespT is indexed by STATION, and station 1 is the delay, so the
% two queues sit at station indices 2 and 3.
sj5 = SolverJMT(model5,'seed',11,'samples',samples);
RD5 = sj5.getCdfRespT();
names5 = model5.getStationNames();
for ist = 1:2
    if ist == 1, want = 'F1'; else, want = 'F2'; end
    sidx = findstring(names5, want);
    FX = RD5{sidx,1};
    if isempty(FX)
        failures{end+1} = sprintf('M5/%s: JMT returned no response-time CDF', want); %#ok<AGROW>
        continue;
    end
    F = FX(:,1); Xv = FX(:,2);
    dF = diff([0;F]);
    s1 = sum(dF.*Xv); s2 = sum(dF.*Xv.^2);
    svar = s2 - s1^2;
    aW = ref5.W(ist,1); aV = ref5.WVar(ist,1);
    fprintf('  M5/%s (JMT)    E[W]    analytic %8.5f  sim %8.5f  rel %7.4f\n', want, aW, s1, relb(aW,s1));
    fprintf('  M5/%s (JMT)    Var[W]  analytic %8.5f  sim %8.5f  rel %7.4f\n', want, aV, svar, relb(aV,svar));
    if relb(aW,s1) > tolMean
        failures{end+1} = sprintf('M5/%s E[W]: analytic %g sim %g', want, aW, s1); %#ok<AGROW>
    end
    if relb(aV,svar) > tolVar
        failures{end+1} = sprintf('M5/%s Var[W]: analytic %g sim %g', want, aV, svar); %#ok<AGROW>
    end
end

fprintf('\n=== summary ===\n');
if isempty(failures)
    fprintf('  ALL SIMULATION CHECKS PASSED\n');
else
    for k = 1:numel(failures)
        fprintf('  FAIL: %s\n', failures{k});
    end
    error('pfqn_sens_simverify:mismatch','%d simulation check(s) outside tolerance', numel(failures));
end
end

% =========================================================================
function bad = compare_station(tag, ref, sim, ist, tolMean, tolVar, tolCov)
% Compare the analytic moments at station ist against a simulated estimate.
% SIM has fields .mean (1 x R), .cov (R x R).
bad = {};
R = size(sim.cov,1);
for r = 1:R
    a = ref.Q(ist,r); s = sim.mean(r);
    e = abs(a-s)/max(1e-12,abs(a));
    fprintf('  %-14s mean(r=%d)  analytic %8.5f  sim %8.5f  rel %7.4f\n', tag, r, a, s, e);
    if e > tolMean
        bad{end+1} = sprintf('%s mean r=%d: analytic %g sim %g rel %g > %g', tag, r, a, s, e, tolMean); %#ok<AGROW>
    end
end
for r = 1:R
    a = ref.QVar(ist,r); s = sim.cov(r,r);
    e = abs(a-s)/max(1e-12,abs(a));
    fprintf('  %-14s Var (r=%d)  analytic %8.5f  sim %8.5f  rel %7.4f\n', tag, r, a, s, e);
    if e > tolVar
        bad{end+1} = sprintf('%s Var r=%d: analytic %g sim %g rel %g > %g', tag, r, a, s, e, tolVar); %#ok<AGROW>
    end
end
for r = 1:R
    for s2 = (r+1):R
        a = ref.QCov(ist,r,s2); s = sim.cov(r,s2);
        e = abs(a-s) / max(0.05, abs(a));   % covariances are small differences
        fprintf('  %-14s Cov (%d,%d)   analytic %8.5f  sim %8.5f  band %7.4f\n', tag, r, s2, a, s, e);
        if e > tolCov
            bad{end+1} = sprintf('%s Cov (%d,%d): analytic %g sim %g band %g > %g', tag, r, s2, a, s, e, tolCov); %#ok<AGROW>
        end
    end
end
end

% =========================================================================
function sim = ldes_moments(model, queues, classes, samples, ist)
% Second moments at station ist via the LDES Markov reward engine. The rewards
% n(i,r) and n(i,r)*n(i,s) are evaluated on the exact joint-state
% residence-time histogram exported by the engine, so no sample-path
% post-processing is needed.
R = numel(classes);
node = queues{ist};
model.clearRewards();
names = cell(0,1);
for r = 1:R
    cr = classes{r};
    nm = sprintf('mean_%d', r);
    model.setReward(nm, @(state) state.at(node,cr));
    names{end+1} = nm; %#ok<AGROW>
end
for r = 1:R
    for s = r:R
        cr = classes{r}; cs = classes{s};
        nm = sprintf('prod_%d_%d', r, s);
        model.setReward(nm, @(state) state.at(node,cr)*state.at(node,cs));
        names{end+1} = nm; %#ok<AGROW>
    end
end
sl = SolverLDES(model,'seed',7,'samples',samples);
[Rw,nmOut] = sl.getAvgReward();
map = containers.Map();
for k = 1:numel(nmOut)
    map(nmOut{k}) = Rw(k);
end
sim.mean = zeros(1,R);
for r = 1:R
    sim.mean(r) = map(sprintf('mean_%d', r));
end
sim.cov = zeros(R,R);
for r = 1:R
    for s = r:R
        m2 = map(sprintf('prod_%d_%d', r, s));
        sim.cov(r,s) = m2 - sim.mean(r)*sim.mean(s);
        sim.cov(s,r) = sim.cov(r,s);
    end
end
model.clearRewards();
end

% =========================================================================
function e = relb(a,s)
% Relative deviation with a floor, so that a near-zero analytic value does not
% make the relative error meaningless.
e = abs(a-s) / max(0.05, abs(a));
end

% =========================================================================
function sim = ldes_total_moments(model, node, classes, samples)
% Moments of the TOTAL queue length at a node, via the LDES reward engine. The
% reward (sum_r n(node,r))^t is nonlinear, which the joint-state
% residence-time histogram evaluates exactly.
model.clearRewards();
model.setReward('tot1', @(state) state.at(node).total());
model.setReward('tot2', @(state) state.at(node).total()^2);
model.setReward('tot3', @(state) state.at(node).total()^3);
sl = SolverLDES(model,'seed',7,'samples',samples);
[Rw,nm] = sl.getAvgReward();
map = containers.Map();
for k = 1:numel(nm)
    map(nm{k}) = Rw(k);
end
sim.m1 = map('tot1');
sim.m2 = map('tot2');
sim.m3 = map('tot3');
model.clearRewards();
end

% =========================================================================
function sim = jmt_moments(model, node, samples)
% Second moments at a node via the JMT logged sample path. sampleAggr returns
% the jump times t and the per-class queue lengths, which are piecewise
% constant on [t(k),t(k+1)), so the stationary moments are the time-weighted
% averages along the path.
sj = SolverJMT(model,'seed',3,'samples',samples);
sa = sj.sampleAggr(node);
t = sa.t(:);
st = sa.state;
dt = diff(t);
n = st(1:end-1,:);
keep = dt > 0;
dt = dt(keep);
n = n(keep,:);
% discard a warmup transient: the path starts from the initial state
nk = numel(dt);
first = max(1, floor(0.1*nk));
dt = dt(first:end);
n = n(first:end,:);
W = dt / sum(dt);
R = size(n,2);
sim.mean = zeros(1,R);
for r = 1:R
    sim.mean(r) = sum(W .* n(:,r));
end
sim.cov = zeros(R,R);
for r = 1:R
    for s = 1:R
        m2 = sum(W .* n(:,r) .* n(:,s));
        sim.cov(r,s) = m2 - sim.mean(r)*sim.mean(s);
    end
end
end
