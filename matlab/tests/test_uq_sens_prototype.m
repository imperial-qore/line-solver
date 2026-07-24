function test_uq_sens_prototype()
% Prototype checks for the epistemic-uncertainty and parametric-sensitivity
% extensions, validated against the closed forms of Trivedi and Bobbio (2017).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
rng(1);

nfail = 0;

%% Phase 1: Jeffreys posterior of Eq. (3.71)
% f(lambda|s) is Erlang with k phases and phase rate s, so its mean is the
% MLE k/s and its variance is k/s^2.
k = 10; s = 5.0;
prior = Prior.fromSample(k, s);
check('fromSample posterior mean = k/s', prior.paramDist.getMean(), k/s, 1e-10);
check('fromSample posterior var = k/s^2', prior.paramDist.getVar(), k/s^2, 1e-10);

%% Phase 1: quantile inversion round-trips through evalCDF
d = Erlang(5.0, 10);
for p = [0.05, 0.5, 0.95]
    x = Prior.quantile(d, p);
    check(sprintf('quantile round-trip p=%.2f', p), d.evalCDF(x), p, 1e-6);
end

%% Phase 2: quadrature discretization integrates the parameter density
% E[1/Lambda] for Lambda ~ Erlang(k phases, rate s) is s/(k-1). The
% discretized prior returns Exp(lambda) alternatives whose means are 1/lambda,
% so the weighted mean of the alternatives estimates E[1/Lambda].
[dists, w] = prior.discretize(2001, 'quadrature');
m = 0;
for i = 1:length(dists)
    m = m + w(i) * dists{i}.getMean();
end
check('quadrature E[1/Lambda] = s/(k-1)', m, s/(k-1), 5e-3);

%% Phase 2: Monte Carlo agrees with quadrature
[dists, w] = prior.discretize(20000, 'montecarlo');
mMC = 0;
for i = 1:length(dists)
    mMC = mMC + w(i) * dists{i}.getMean();
end
check('montecarlo E[1/Lambda] = s/(k-1)', mMC, s/(k-1), 5e-2);

%% Phase 2: UQ end to end on an M/M/1 with an uncertain service rate
% Two Priors exercise the multi-Prior path that the old single-Prior guard
% rejected. The design is the tensor product, so 5x5 = 25 models.
model = Network('uq_mm1');
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');
oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Prior({Exp(0.4), Exp(0.5), Exp(0.6)}, [0.25, 0.5, 0.25]));
queue.setService(oclass, Prior({Exp(1.0), Exp(1.2)}, [0.5, 0.5]));
model.link(model.serialRouting(source, queue, sink));

uq = UQ(model, @SolverMVA);
uq.iterate();
check('multi-Prior design size = 3 x 2', uq.getNumberOfModels(), 6, 0);
check('design weights sum to 1', sum(uq.getProbabilities()), 1, 1e-12);

% The weighted mean must equal the manual product-weighted mean of the M/M/1
% response times R = 1/(mu - lambda).
lams = [0.4, 0.5, 0.6]; lamw = [0.25, 0.5, 0.25];
mus = [1.0, 1.2];       muw = [0.5, 0.5];
Rexact = 0;
for a = 1:3
    for b = 1:2
        Rexact = Rexact + lamw(a) * muw(b) * 1/(mus(b) - lams(a));
    end
end
[mR, vR] = uq.getMoments('R', queue, oclass);
check('UQ E[R] matches manual mixture', mR, Rexact, 1e-6);

% Variance from the same weighting
Vexact = 0;
for a = 1:3
    for b = 1:2
        Vexact = Vexact + lamw(a) * muw(b) * (1/(mus(b) - lams(a)) - Rexact)^2;
    end
end
check('UQ Var[R] matches manual mixture', vR, Vexact, 1e-6);

ci = uq.getCredibleInterval('R', queue, oclass, 0.90);
fprintf(1, 'INFO 90%% credible interval for R = [%.4f, %.4f]\n', ci(1), ci(2));
if ci(1) <= mR && mR <= ci(2)
    fprintf(1, 'PASS %-46s\n', 'credible interval brackets the mean');
else
    fprintf(1, 'FAIL %-46s\n', 'credible interval brackets the mean');
    nfail = nfail + 1;
end

%% Regression: Monte Carlo over a discrete prior must honour its probabilities
% A design point drawn from a discrete prior once always returned the first
% alternative with weight 1/n, so the prior was silently dropped. Weighting the
% alternatives 0.01/0.99 makes that failure a factor-of-ten error, not a
% rounding difference.
pskew = Prior({Exp(1.0), Exp(10.0)}, [0.01, 0.99]);
[dskew, wskew] = pskew.discretize(4000, 'montecarlo');
mskew = 0;
for i = 1:length(dskew)
    mskew = mskew + wskew(i) * dskew{i}.getMean();
end
check('montecarlo discrete honours probabilities', mskew, 0.01*1.0 + 0.99*0.1, 2e-2);

skewModel = Network('uq_skew');
sSrc = Source(skewModel, 'Source');
sQ = Queue(skewModel, 'Queue', SchedStrategy.FCFS);
sSnk = Sink(skewModel, 'Sink');
sCls = OpenClass(skewModel, 'Class1');
sSrc.setArrival(sCls, Exp(0.1));
sQ.setService(sCls, pskew);
skewModel.link(skewModel.serialRouting(sSrc, sQ, sSnk));
uqmc = UQ(skewModel, @SolverMVA, 'method', 'montecarlo', 'samples', 400);
uqmc.iterate();
mRmc = uqmc.getMoments('R', sQ, sCls);
Rskew = 0.01 * 1/(1.0 - 0.1) + 0.99 * 1/(10.0 - 0.1);
check('UQ montecarlo E[R] over discrete prior', mRmc, Rskew, 5e-2);

%% Regression: enumerating a continuous prior must fail loudly
% getNumAlternatives returns NaN for a continuous prior and every comparison
% against NaN is false, so an unguarded bounds check silently passed.
cprior = Prior.fromSample(10, 5.0);
threw = false;
try
    cprior.getAlternative(1);
catch
    threw = true;
end
if threw
    fprintf(1, 'PASS %-46s\n', 'getAlternative rejects a continuous prior');
else
    fprintf(1, 'FAIL %-46s\n', 'getAlternative rejects a continuous prior');
    nfail = nfail + 1;
end

%% Regression: sampling a continuous prior uses the density, not a grid
% Sampling via discretize would return the law of an 11-node approximation.
% For Lambda ~ Erlang(k, s) and X|Lambda ~ Exp(Lambda), E[X] = E[1/Lambda] =
% s/(k-1), and the sample must converge to it from the true mixture.
xs = cprior.sample(40000);
check('continuous prior sample mean = s/(k-1)', mean(xs), 5.0/(10-1), 5e-2);
if numel(unique(xs)) > 11
    fprintf(1, 'PASS %-46s\n', 'sample is not confined to a quadrature grid');
else
    fprintf(1, 'FAIL %-46s\n', 'sample is not confined to a quadrature grid');
    nfail = nfail + 1;
end

%% Phase 3: ctmc_sens against the closed form for a two-state chain
% Q = [-a, a; b, -b] has pi = [b, a]/(a+b), so dpi_1/da = -b/(a+b)^2.
a = 0.3; b = 0.7;
Q = [-a, a; b, -b];
dQ = [-1, 1; 0, 0];      % d/da
pi2 = ctmc_solve(Q);
dpi = ctmc_sens(Q, dQ, pi2);
check('ctmc_sens dpi_1/da closed form', dpi(1), -b/(a+b)^2, 1e-10);
check('ctmc_sens dpi_2/da closed form', dpi(2),  b/(a+b)^2, 1e-10);
check('ctmc_sens sensitivities sum to 0', sum(dpi), 0, 1e-10);

%% Phase 3: availability sensitivity, the book's Sec. 9.7 setting
% A two-state up/down chain with failure rate gamma and repair rate mu has
% steady-state availability A = mu/(gamma+mu), so dA/dgamma = -mu/(gamma+mu)^2
% and the scaled sensitivity is -gamma/(gamma+mu). The sign is negative: a
% higher failure rate lowers availability, as in the book's Table 9.3.
gamma = 1/240; mu = 1/4;
Qav = [-gamma, gamma; mu, -mu];
dQav = [-1, 1; 0, 0];
piav = ctmc_solve(Qav);
dpiav = ctmc_sens(Qav, dQav, piav);
r = [1; 0];              % reward 1 in the up state
A = piav * r;
dA = dpiav * r;
check('availability A = mu/(gamma+mu)', A, mu/(gamma+mu), 1e-12);
check('dA/dgamma closed form', dA, -mu/(gamma+mu)^2, 1e-9);
check('scaled sensitivity SS = -gamma/(gamma+mu)', (gamma/A)*dA, -gamma/(gamma+mu), 1e-9);

%% Phase 3: transient sensitivity converges to the steady-state one
% Tolerance follows the ode23 default of ctmc_transient, which is relative:
% dpi is O(1), so absolute agreement to 1e-3 is what the integrator delivers.
[dpit, pit, tt] = ctmc_transient_sens(Qav, dQav, [1, 0], 0, 5000);
check('transient pi converges to steady state', pit(end,1), piav(1), 1e-4);
check('transient dpi converges to ctmc_sens', dpit(end,1), dpiav(1), 1e-3);
check('transient dpi starts at zero', dpit(1,1), 0, 1e-12);

%% Phase 3/4: SolverCTMC.getSensitivity against a finite-difference baseline
% M/M/1/K queue length sensitivity to the service rate.
Kcap = 4;
mm1k = Network('mm1k');
src = Source(mm1k, 'Source');
q = Queue(mm1k, 'Queue', SchedStrategy.FCFS);
snk = Sink(mm1k, 'Sink');
cls = OpenClass(mm1k, 'Class1');
src.setArrival(cls, Exp(0.5));
q.setService(cls, Exp(1.0));
q.setNumberOfServers(1);
mm1k.link(mm1k.serialRouting(src, q, snk));

opts = SolverCTMC.defaultOptions();
opts.cutoff = Kcap;
opts.verbose = 0;
solver = SolverCTMC(mm1k, opts);
solver.getGenerator();   % populates the cached state space
spaceAggr = solver.getStateSpaceAggr();
% Column (ist-1)*K+k holds the queue length of station ist in class k. Station
% 1 is the Source, whose column is Inf, so the queue is station 2 column 2.
rQ = spaceAggr(:, 2);

param = struct();
param.name = 'mu';
param.value = 1.0;
% The setter resolves node and class inside the copy it is handed, since
% the objects of the original model do not belong to it.
param.set = @(m, v) setNodeService(m, 2, 1, v);

[S, SS] = solver.getSensitivity(param, rQ);

% Finite-difference baseline on the solved metric itself
h = 1e-4;
qOf = @(v) fdQueueLen(Kcap, v, rQ);
Sfd = (qOf(1.0 + h) - qOf(1.0 - h)) / (2*h);
check('getSensitivity dQ/dmu vs finite diff', S, Sfd, 1e-4);
fprintf(1, 'INFO scaled sensitivity of QLen to mu = %.6g\n', SS);
if S < 0
    fprintf(1, 'PASS %-46s\n', 'faster service lowers queue length');
else
    fprintf(1, 'FAIL %-46s\n', 'faster service lowers queue length');
    nfail = nfail + 1;
end

%% Phase 4: ranking orders parameters by absolute scaled sensitivity
paramMu = param;
paramLam = struct('name', 'lambda', 'value', 0.5, ...
    'set', @(m, v) setNodeArrival(m, 1, 1, v));
RankTable = solver.getSensitivityRanking({paramMu, paramLam}, rQ);
disp(RankTable);
ss = abs(RankTable.ScaledSens);
if issorted(ss, 'descend')
    fprintf(1, 'PASS %-46s\n', 'ranking sorted by |scaled sensitivity|');
else
    fprintf(1, 'FAIL %-46s\n', 'ranking sorted by |scaled sensitivity|');
    nfail = nfail + 1;
end

%% Phase 5: CTMC fallback fires where the product-form path returns []
% A closed network with a multiserver queue is exactly the case
% closedSensitivities rejects, so it isolates the fallback.
cqn = Network('cqn_multiserver');
delay = Delay(cqn, 'Delay');
mq = Queue(cqn, 'MultiQueue', SchedStrategy.FCFS);
mq.setNumberOfServers(2);
ccls = ClosedClass(cqn, 'Class1', 3, delay, 0);
delay.setService(ccls, Exp(1.0));
mq.setService(ccls, Exp(2.0));
Pc = cqn.initRoutingMatrix;
Pc{1} = [0, 1; 1, 0];
cqn.link(Pc);

pfSens = opt.sens.computeModelSensitivities(cqn);
if isempty(pfSens)
    fprintf(1, 'PASS %-46s\n', 'product-form path returns [] for multiserver');
else
    fprintf(1, 'FAIL %-46s\n', 'product-form path returns [] for multiserver');
    nfail = nfail + 1;
end

ctmcSens = opt.sens.computeModelSensitivities(cqn, true);
if ~isempty(ctmcSens) && ctmcSens.data.isKey('QLen')
    fprintf(1, 'PASS %-46s\n', 'CTMC fallback covers the multiserver case');
    qmap = ctmcSens.data('QLen');
    mk = opt.SensitivityData.metricKey('MultiQueue', 'Class1');
    pk = opt.SensitivityData.paramKey('MultiQueue', 'Class1');
    byParam = qmap(mk);
    dQdmu = byParam(pk);
    fprintf(1, 'INFO d(QLen at MultiQueue)/d(its rate) = %.6g\n', dQdmu);
    if dQdmu < 0
        fprintf(1, 'PASS %-46s\n', 'faster multiserver lowers its own queue');
    else
        fprintf(1, 'FAIL %-46s\n', 'faster multiserver lowers its own queue');
        nfail = nfail + 1;
    end
else
    fprintf(1, 'FAIL %-46s\n', 'CTMC fallback covers the multiserver case');
    nfail = nfail + 1;
end

%% Regression: the CTMC fallback must skip non-exponential service
% The rate setter substitutes an Exp, so perturbing an Erlang station would
% change the distribution family and report a difference quotient that is not
% dQ/dtheta. Such stations must be absent from the result, not wrong in it.
erlModel = Network('cqn_erlang');
eDelay = Delay(erlModel, 'Delay');
eQ = Queue(erlModel, 'ErlQueue', SchedStrategy.FCFS);
eCls = ClosedClass(erlModel, 'Class1', 2, eDelay, 0);
eDelay.setService(eCls, Exp(1.0));
eQ.setService(eCls, Erlang(4.0, 2));   % non-exponential, must be skipped
eQ.setNumberOfServers(2);              % multiserver, so the PF path declines
Pe = erlModel.initRoutingMatrix;
Pe{1} = [0, 1; 1, 0];
erlModel.link(Pe);

erlSens = opt.sens.computeModelSensitivities(erlModel, true);
pkErl = opt.SensitivityData.paramKey('ErlQueue', 'Class1');
foundErl = false;
if ~isempty(erlSens) && erlSens.data.isKey('QLen')
    qm = erlSens.data('QLen');
    mkeys = qm.keys();
    for i = 1:numel(mkeys)
        bp = qm(mkeys{i});
        if bp.isKey(pkErl)
            foundErl = true;
        end
    end
end
if ~foundErl
    fprintf(1, 'PASS %-46s\n', 'fallback skips Erlang service rather than guess');
else
    fprintf(1, 'FAIL %-46s\n', 'fallback skips Erlang service rather than guess');
    nfail = nfail + 1;
end

fprintf(1, '\n==== %d failure(s) ====\n', nfail);


    function check(name, got, want, tol)
        if abs(got - want) <= tol
            fprintf(1, 'PASS %-46s got %.6g want %.6g\n', name, got, want);
        else
            fprintf(1, 'FAIL %-46s got %.6g want %.6g\n', name, got, want);
            nfail = nfail + 1;
        end
    end

    % Baseline built from a freshly constructed model, so it shares no cached
    % struct with the model under test and is an independent check.
    function Qlen = fdQueueLen(Kcap, mu, rQ)
        m = Network('fd');
        s1 = Source(m, 'Source');
        q1 = Queue(m, 'Queue', SchedStrategy.FCFS);
        k1 = Sink(m, 'Sink');
        c1 = OpenClass(m, 'Class1');
        s1.setArrival(c1, Exp(0.5));
        q1.setService(c1, Exp(mu));
        q1.setNumberOfServers(1);
        m.link(m.serialRouting(s1, q1, k1));
        o = SolverCTMC.defaultOptions();
        o.cutoff = Kcap;
        o.verbose = 0;
        sv = SolverCTMC(m, o);
        Qgen = full(sv.getGenerator());
        p = ctmc_solve(Qgen);
        Qlen = p * rQ;
    end

    function setNodeService(m, nodeIdx, classIdx, value)
        nn = m.getNodes();
        cc = m.getClasses();
        nn{nodeIdx}.setService(cc{classIdx}, Exp(value));
    end

    function setNodeArrival(m, nodeIdx, classIdx, value)
        nn = m.getNodes();
        cc = m.getClasses();
        nn{nodeIdx}.setArrival(cc{classIdx}, Exp(value));
    end

end
