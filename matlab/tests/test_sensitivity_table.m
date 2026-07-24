function test_sensitivity_table()
% TEST_SENSITIVITY_TABLE Regression tests for getSensitivityTable.
%
% Covers the dispatch between the analytic branch (SolverMVA, SolverNC) and
% the finite-difference branch (every other solver), the agreement of the two
% branches on a model both can handle, the option validation, and the rate
% scaling primitive the finite-difference branch perturbs with.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

nfail = 0;
nfail = nfail + t_dispatch();
nfail = nfail + t_branch_agreement();
nfail = nfail + t_open_model();
nfail = nfail + t_forced_exact();
nfail = nfail + t_option_validation();
nfail = nfail + t_scaling_primitive();
nfail = nfail + t_nonexponential();
nfail = nfail + t_class_switching();
nfail = nfail + t_nonunit_visits();
nfail = nfail + t_layered();
nfail = nfail + t_simulation();

if nfail > 0
    line_error(mfilename, sprintf('%d sensitivity-table checks failed.', nfail));
end
fprintf('test_sensitivity_table: all checks passed.\n');
end

function model = closedModel()
model = Network('closed');
d = Delay(model, 'D');
q = Queue(model, 'Q', SchedStrategy.PS);
c = ClosedClass(model, 'C', 3, d);
d.setService(c, Exp(1));
q.setService(c, Exp(2));
model.link(Network.serialRouting(d, q));
end

function model = openModel()
model = Network('open');
s = Source(model, 'S');
q = Queue(model, 'Q', SchedStrategy.PS);
k = Sink(model, 'K');
c1 = OpenClass(model, 'A');
c2 = OpenClass(model, 'B');
s.setArrival(c1, Exp(0.4));
s.setArrival(c2, Exp(0.2));
q.setService(c1, Exp(2));
q.setService(c2, Exp(1.5));
P = model.initRoutingMatrix;
P{c1} = Network.serialRouting(s, q, k);
P{c2} = Network.serialRouting(s, q, k);
model.link(P);
end

function model = multiserverModel()
model = Network('multi');
d = Delay(model, 'D');
q = Queue(model, 'Q', SchedStrategy.PS);
q.setNumberOfServers(2);
c = ClosedClass(model, 'C', 3, d);
d.setService(c, Exp(1));
q.setService(c, Exp(2));
model.link(Network.serialRouting(d, q));
end

function n = check(cond, msg)
n = 0;
if ~cond
    n = 1;
    fprintf(2, 'FAIL: %s\n', msg);
end
end

function n = checkClose(got, want, tol, msg)
n = check(abs(got - want) <= tol, sprintf('%s (got %.10g, want %.10g)', msg, got, want));
end

function n = t_dispatch()
% The analytic branch belongs to the solvers that evaluate the recursion it
% differentiates; every other solver finite-differences its own predictions.
n = 0;
n = n + check(strcmp(SolverMVA(closedModel()).getSensitivityTable().Properties.UserData.method, 'exact'), 'MVA must use the exact branch');
n = n + check(strcmp(SolverNC(closedModel()).getSensitivityTable().Properties.UserData.method, 'exact'), 'NC must use the exact branch');
n = n + check(strcmp(SolverCTMC(closedModel()).getSensitivityTable().Properties.UserData.method, 'fd'), 'CTMC must use finite differences');
n = n + check(strcmp(SolverFluid(closedModel()).getSensitivityTable().Properties.UserData.method, 'fd'), 'Fluid must use finite differences');
% A model out of scope of the analytic branch falls back rather than failing.
T = SolverMVA(multiserverModel()).getSensitivityTable();
n = n + check(strcmp(T.Properties.UserData.method, 'fd'), 'a multiserver model must fall back to finite differences');
end

function n = t_branch_agreement()
% Reference values shared with the Python and JAR ports.
n = 0;
Te = SolverMVA(closedModel()).getSensitivityTable();
n = n + checkClose(Te.dTput_dRate(1),   0.49030470914, 1e-9, 'exact dTput');
n = n + checkClose(Te.dRespT_dRate(1), -0.59,          1e-9, 'exact dRespT');
n = n + checkClose(Te.dQLen_dRate(1),  -0.49030470914, 1e-9, 'exact dQLen');
n = n + checkClose(Te.dUtil_dRate(1),  -0.14958448753, 1e-9, 'exact dUtil');

Tf = SolverMVA(closedModel()).getSensitivityTable('method', 'fd', 'scheme', 'central');
n = n + checkClose(Tf.dTput_dRate(1),  Te.dTput_dRate(1),  1e-5, 'fd central dTput');
n = n + checkClose(Tf.dRespT_dRate(1), Te.dRespT_dRate(1), 1e-5, 'fd central dRespT');
n = n + checkClose(Tf.dQLen_dRate(1),  Te.dQLen_dRate(1),  1e-5, 'fd central dQLen');
n = n + checkClose(Tf.dUtil_dRate(1),  Te.dUtil_dRate(1),  1e-5, 'fd central dUtil');

% The CTMC is exact on this model, so its finite differences must land on the
% same derivative, up to the forward-scheme truncation error.
Tc = SolverCTMC(closedModel()).getSensitivityTable();
n = n + checkClose(Tc.dTput_dRate(1),  Te.dTput_dRate(1),  1e-3, 'CTMC fd dTput');
n = n + checkClose(Tc.dRespT_dRate(1), Te.dRespT_dRate(1), 1e-3, 'CTMC fd dRespT');
n = n + checkClose(Tc.dQLen_dRate(1),  Te.dQLen_dRate(1),  1e-3, 'CTMC fd dQLen');
n = n + checkClose(Tc.dUtil_dRate(1),  Te.dUtil_dRate(1),  1e-3, 'CTMC fd dUtil');
end

function n = t_open_model()
% An open network's throughput is fixed by the arrival rate, so the analytic
% branch reports a zero rate derivative for it.
n = 0;
Te = SolverMVA(openModel()).getSensitivityTable();
n = n + check(all(Te.dTput_dRate == 0), 'open dTput_dRate must be exactly zero');
n = n + check(height(Te) == 2, 'the open model must produce one row per class');
% The finite differences follow MVA's own open solution, which is itself
% approximate here, so the two branches agree only to a few digits.
Tf = SolverMVA(openModel()).getSensitivityTable('method', 'fd', 'scheme', 'central');
n = n + checkClose(Tf.dRespT_dRate(1), Te.dRespT_dRate(1), 1e-3, 'open fd dRespT');
n = n + checkClose(Tf.dUtil_dRate(1),  Te.dUtil_dRate(1),  1e-6, 'open fd dUtil');
end

function n = t_forced_exact()
% method='exact' reports the restriction instead of silently approximating.
n = 0;
n = n + check(throwsError(@() SolverCTMC(closedModel()).getSensitivityTable('method', 'exact')), ...
    'CTMC must reject method=exact');
n = n + check(throwsError(@() SolverMVA(multiserverModel()).getSensitivityTable('method', 'exact')), ...
    'a multiserver model must reject method=exact');
end

function n = t_option_validation()
n = 0;
n = n + check(throwsError(@() SolverMVA(closedModel()).getSensitivityTable('method', 'nosuch')), ...
    'an unknown method must be rejected');
n = n + check(throwsError(@() SolverMVA(closedModel()).getSensitivityTable('scheme', 'nosuch')), ...
    'an unknown scheme must be rejected');
n = n + check(throwsError(@() SolverCTMC(closedModel()).getSensitivityTable('step', 2)), ...
    'a step outside (0,1) must be rejected');
n = n + check(throwsError(@() SolverMVA(closedModel()).getSensitivityTable('nosuch', 1)), ...
    'an unknown option name must be rejected');
end

function n = t_scaling_primitive()
% Scaling the rate by a factor divides every mean by that factor and leaves
% the SCV, and hence the shape of the distribution, untouched. That property
% is what makes the perturbation a pure rate change.
n = 0;
factor = 2;
dists = {Exp(2), Erlang(3, 2), HyperExp(0.3, 2, 5), Coxian([2 4], [0.3 1]), ...
    APH.fit(1, 2, 6), Det(0.5), Uniform(1, 3), Gamma(2, 0.5), Pareto(3, 1), ...
    Weibull(2, 1), Lognormal(0, 0.5), MAP(-[2 0; 0 3], [1 1; 1.5 1.5]), ...
    NHPP([0 1 2], [2 4], true), Replayer([1 2 3])};
for i = 1:numel(dists)
    d = dists{i};
    s = dist_scale_rate(d, factor);
    n = n + check(abs(s.getMean() - d.getMean()/factor) <= 1e-9 * max(1, abs(d.getMean())), ...
        sprintf('%s mean must scale by 1/factor', class(d)));
    if ~isnan(d.getSCV())
        n = n + check(abs(s.getSCV() - d.getSCV()) <= 1e-9 * max(1, abs(d.getSCV())), ...
            sprintf('%s SCV must be invariant', class(d)));
    end
end
n = n + check(throwsError(@() dist_scale_rate(Exp(1), -1)), 'a non-positive factor must be rejected');
end

function n = t_nonexponential()
% End-to-end coverage of the scaling path for a non-exponential process.
n = 0;
model = Network('erl');
d = Delay(model, 'D');
q = Queue(model, 'Q', SchedStrategy.FCFS);
c = ClosedClass(model, 'C', 3, d);
d.setService(c, Exp(1));
q.setService(c, Erlang.fitMeanAndSCV(0.5, 0.5));
model.link(Network.serialRouting(d, q));
T = SolverCTMC(model).getSensitivityTable();
n = n + check(strcmp(T.Properties.UserData.method, 'fd'), 'an Erlang FCFS model must use finite differences');
n = n + check(all(isfinite(T.dRespT_dRate)) && all(isfinite(T.dQLen_dRate)), ...
    'the Erlang finite differences must be finite');
n = n + check(T.dRespT_dRate(1) < 0, 'a faster server must shorten the response time');
end

function model = classSwitchModel()
% Two classes that switch into each other, so the chain population sits on the
% reference class alone and the switched class carries zero jobs.
model = Network('cs');
d = Delay(model, 'D');
q = Queue(model, 'Q', SchedStrategy.PS);
cs = ClassSwitch(model, 'CS');
c1 = ClosedClass(model, 'C1', 2, d);
c2 = ClosedClass(model, 'C2', 0, d);
d.setService(c1, Exp(1)); d.setService(c2, Exp(1));
q.setService(c1, Exp(2)); q.setService(c2, Exp(3));
C = zeros(2, 2); C(1, 2) = 1; C(2, 1) = 1;
cs.setClassSwitchingMatrix(C);
P = model.initRoutingMatrix;
P{c1, c2} = Network.serialRouting(d, q, cs);
P{c2, c1} = Network.serialRouting(d, q, cs);
model.link(P);
end

function model = layeredModel()
model = LayeredNetwork('LQN-single');
P1 = Processor(model, 'P1', 1, SchedStrategy.PS);
P2 = Processor(model, 'P2', 1, SchedStrategy.PS);
T1 = Task(model, 'T1', 5, SchedStrategy.REF).on(P1).setThinkTime(Exp(1/2));
T2 = Task(model, 'T2', 1, SchedStrategy.FCFS).on(P2).setThinkTime(Exp(1/3));
E1 = Entry(model, 'E1').on(T1);
E2 = Entry(model, 'E2').on(T2);
Activity(model, 'AS1', Exp(10)).on(T1).boundTo(E1).synchCall(E2, 1);
Activity(model, 'AS2', Exp(20)).on(T2).boundTo(E2).repliesTo(E2);
end

function n = t_class_switching()
% A chain-based model is differentiated at chain level and disaggregated back to
% the classes, so the analytic branch applies and must agree with the finite
% differences of the same solver.
n = 0;
Te = SolverMVA(classSwitchModel()).getSensitivityTable();
n = n + check(strcmp(Te.Properties.UserData.method, 'exact'), ...
    'a class-switching model must use the analytic branch');
n = n + check(height(Te) == 2, 'both switched classes must be reported');
Tf = SolverMVA(classSwitchModel()).getSensitivityTable('method', 'fd', 'scheme', 'central');
for i = 1:height(Te)
    n = n + checkClose(Te.dTput_dRate(i),  Tf.dTput_dRate(i),  1e-5, 'class-switch dTput');
    n = n + checkClose(Te.dRespT_dRate(i), Tf.dRespT_dRate(i), 1e-5, 'class-switch dRespT');
    n = n + checkClose(Te.dQLen_dRate(i),  Tf.dQLen_dRate(i),  1e-5, 'class-switch dQLen');
    n = n + checkClose(Te.dUtil_dRate(i),  Tf.dUtil_dRate(i),  1e-5, 'class-switch dUtil');
end
end

function n = t_nonunit_visits()
% With visit ratios other than one, a class throughput at a station is X*v and
% the reported response time is per visit, not the chain residence time. Reading
% the chain quantities straight out of pfqn_sens got both wrong by a factor v.
n = 0;
model = Network('visits');
d = Delay(model, 'D');
q1 = Queue(model, 'Q1', SchedStrategy.PS);
q2 = Queue(model, 'Q2', SchedStrategy.PS);
c = ClosedClass(model, 'C', 3, d);
d.setService(c, Exp(1)); q1.setService(c, Exp(2)); q2.setService(c, Exp(3));
P = model.initRoutingMatrix;
P{c, c} = [0 1 0; 0 0 1; 0.5 0.5 0];   % D -> Q1 -> Q2 -> {D, Q1}: v(Q1)=v(Q2)=2
model.link(P);
Te = SolverMVA(model).getSensitivityTable();
Tf = SolverMVA(model).getSensitivityTable('method', 'fd', 'scheme', 'central');
n = n + check(strcmp(Te.Properties.UserData.method, 'exact'), 'the model must be in scope');
for i = 1:height(Te)
    n = n + checkClose(Te.dTput_dRate(i),  Tf.dTput_dRate(i),  1e-5, 'non-unit visits dTput');
    n = n + checkClose(Te.dRespT_dRate(i), Tf.dRespT_dRate(i), 1e-5, 'non-unit visits dRespT');
    n = n + checkClose(Te.dQLen_dRate(i),  Tf.dQLen_dRate(i),  1e-5, 'non-unit visits dQLen');
    n = n + checkClose(Te.dUtil_dRate(i),  Tf.dUtil_dRate(i),  1e-5, 'non-unit visits dUtil');
end
end

function n = t_layered()
% SolverLN has no recursion of its own to differentiate: it concatenates what
% the layer solvers report, one row block per layer.
n = 0;
solver = SolverLN(layeredModel(), @(m) SolverMVA(m), 'verbose', false);
before = solver.getAvgTable();
T = solver.getSensitivityTable();
n = n + check(height(T) == 3, 'the layered model must produce one row per layer');
n = n + check(isequal(T.Properties.VariableNames, {'Layer', 'Station', 'JobClass', ...
    'dTput_dRate', 'dRespT_dRate', 'dQLen_dRate', 'dUtil_dRate'}), ...
    'the layered table must carry a leading Layer column');
% A layer submodel is chain-based, which the analytic branch now handles by
% aggregating to chains, so the layers differentiate analytically.
n = n + check(strcmp(T.Properties.UserData.method, 'exact'), 'the layers must use the analytic branch');
n = n + checkClose(T.dTput_dRate(1),   0.015211,    1e-5, 'layer P1 dTput');
n = n + checkClose(T.dRespT_dRate(1), -0.014406,    1e-5, 'layer P1 dRespT');
n = n + checkClose(T.dQLen_dRate(1),  -0.031257,    1e-5, 'layer P1 dQLen');
n = n + checkClose(T.dUtil_dRate(1),  -0.021453,    1e-5, 'layer P1 dUtil');
n = n + checkClose(T.dRespT_dRate(2), -0.0024998,   1e-6, 'layer P2 dRespT');
n = n + checkClose(T.dRespT_dRate(3), -0.003003,    1e-6, 'layer T2 dRespT');
% The sweep restores every perturbed service process, so the layered averages
% are unchanged by it, up to the fixed-point residual of a second solve.
after = solver.getAvgTable();
n = n + check(max(abs(before.RespT - after.RespT)) < 1e-3, ...
    'the layer sweep must leave the fixed point intact');
end

function n = t_simulation()
% The simulation branch: common random numbers, the simulator default step, and
% the sign structure of the estimates. The seed is fixed, so the whole sweep is
% a deterministic function of it and the reproducibility check is exact.
n = 0;
samples = 20000;
seed = 23000;
mkSSA = @() SolverSSA(closedModel(), 'samples', samples, 'seed', seed, 'verbose', false);

T1 = mkSSA().getSensitivityTable();
n = n + check(strcmp(T1.Properties.UserData.method, 'fd'), ...
    'a simulation solver must finite-difference');

% Common random numbers: two independent instances with the same seed must
% return the same numbers. This is what fails first if the base and perturbed
% runs stop sharing a seed.
T2 = mkSSA().getSensitivityTable();
n = n + checkClose(T2.dTput_dRate(1),  T1.dTput_dRate(1),  1e-12, 'CRN dTput');
n = n + checkClose(T2.dRespT_dRate(1), T1.dRespT_dRate(1), 1e-12, 'CRN dRespT');
n = n + checkClose(T2.dQLen_dRate(1),  T1.dQLen_dRate(1),  1e-12, 'CRN dQLen');
n = n + checkClose(T2.dUtil_dRate(1),  T1.dUtil_dRate(1),  1e-12, 'CRN dUtil');

% A seed is always present on a MATLAB simulation solver: SolverSSA draws a
% random one at construction, so the pinning path never fires here and the
% invariant to check is that the sweep reuses whatever seed the solver holds.
% Two sweeps on the same instance must therefore agree exactly.
unseeded = SolverSSA(closedModel(), 'samples', samples, 'verbose', false);
Ta = unseeded.getSensitivityTable();
seedAfter = unseeded.options.seed;
Tb = unseeded.getSensitivityTable();
n = n + check(isfinite(seedAfter), 'the solver must carry a seed for the sweep');
n = n + check(unseeded.options.seed == seedAfter, 'the sweep must not move the seed');
n = n + checkClose(Tb.dTput_dRate(1), Ta.dTput_dRate(1), 1e-12, 'repeat sweep dTput');
n = n + checkClose(Tb.dQLen_dRate(1), Ta.dQLen_dRate(1), 1e-12, 'repeat sweep dQLen');

% The simulator default step is 1e-2, asserted through its effect.
Tdefault = mkSSA().getSensitivityTable();
Tsame = mkSSA().getSensitivityTable('step', 1e-2);
Tother = mkSSA().getSensitivityTable('step', 1e-3);
n = n + checkClose(Tdefault.dTput_dRate(1), Tsame.dTput_dRate(1), 1e-12, ...
    'the simulator default step must be 1e-2');
n = n + check(abs(Tdefault.dTput_dRate(1) - Tother.dTput_dRate(1)) > 1e-12, ...
    'a different step must give different estimates');

% Sign structure, and a deliberately wide magnitude band: at 20000 samples the
% Monte Carlo error is tens of percent, so anything tighter would be flaky. It
% still catches a missing visit factor or an unpaired seed.
Te = SolverMVA(closedModel()).getSensitivityTable();
n = n + check(T1.dTput_dRate(1) > 0, 'a faster server must raise throughput');
n = n + check(T1.dRespT_dRate(1) < 0 && T1.dQLen_dRate(1) < 0 && T1.dUtil_dRate(1) < 0, ...
    'a faster server must shorten the queue');
cols = {'dTput_dRate', 'dRespT_dRate', 'dQLen_dRate', 'dUtil_dRate'};
for i = 1:numel(cols)
    got = abs(T1.(cols{i})(1));
    want = abs(Te.(cols{i})(1));
    n = n + check(got >= 0.5*want && got <= 2*want, ...
        sprintf('simulated %s must be within a factor of two of the analytic value', cols{i}));
end
end

function tf = throwsError(fun)
tf = false;
try
    fun();
catch
    tf = true;
end
end
