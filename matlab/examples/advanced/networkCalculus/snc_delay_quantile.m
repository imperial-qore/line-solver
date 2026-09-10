% Stochastic network calculus: a delay quantile with a certified violation
% probability.
%
% Every other solver in LINE answers with a MEAN. The 'snc' family of SolverBA
% answers with a TAIL: given a violation probability eps, it returns a delay d
% for which P{D > d} <= eps holds, and the guarantee is valid for any
% work-conserving scheduling policy at the station. That is the quantity a
% service-level objective is written against.
%
% The model is a single M/M/1 station, whose exact tail is known in closed
% form, so every number below can be checked.

%% Block 1: model
lambda = 0.6;
mu = 1;
model = Network('SncQuantile');
source = Source(model,'Source');
queue = Queue(model,'Queue', SchedStrategy.FCFS);
sink = Sink(model,'Sink');
jobclass = OpenClass(model,'Class1');
source.setArrival(jobclass, Exp(lambda));
queue.setService(jobclass, Exp(mu));
model.link(Network.serialRouting(source,queue,sink));

solver = SolverBA(model,'method','snc.upper');

%% Block 2: the quantile table, the native output of the family
% getPercTable is to the snc family what getAvgTable is to a mean solver: one
% row per station and class, reporting the response-time and queue-length
% quantiles at the requested violation probability.
percTable = solver.getPercTable(1e-3)

%% Block 3: how the bound tightens as the guarantee gets stricter
% The exact M/M/1 sojourn tail is P{D>d} = exp(-(mu-lambda)*d) and the exact
% queue-length tail is P{Q>n} = rho^(n+1). The bound reproduces both DECAY
% RATES exactly and pays a constant prefactor, so the ratio of the bounded
% quantile to the exact one falls towards 1 as eps is tightened. A bound that
% is 2x at eps=1e-2 is 1.2x at eps=1e-12: the family is at its best exactly
% where simulation is at its worst, deep in the tail.
fprintf('\n%-8s %10s %10s %8s   %10s %10s %8s\n', ...
    'eps','d bound','d exact','ratio','n bound','n exact','ratio');
for eps = [1e-2 1e-3 1e-6 1e-9 1e-12]
    D = solver.getDelayPerc(eps);
    B = solver.getBacklogPerc(eps);
    dexact = -log(eps)/(mu-lambda);
    nexact = log(eps)/log(lambda/mu) - 1;
    fprintf('%-8.0e %10.4f %10.4f %8.3f   %10.4f %10.4f %8.3f\n', ...
        eps, D(2), dexact, D(2)/dexact, B(2), nexact, B(2)/nexact);
end

%% Block 4: the mean columns, and why they are the loose end
% getAvgTable still works: the response time reported by 'snc.upper' is the
% integral of the tail bound, hence an upper bound on the mean. It is loose,
% and deliberately so -- integrating over the whole axis is dominated by the
% prefactor rather than by the decay rate that the family gets right. Use the
% mean columns to bracket, the quantiles to plan.
avgTable = solver.getAvgTable()
exactR = 1/(mu-lambda);
exactQ = (lambda/mu)/(1-lambda/mu);
fprintf('\nexact M/M/1: R = %.4f, Q = %.4f\n', exactR, exactQ);
fprintf('snc.upper  : R = %.4f (%.1fx), Q = %.4f (%.1fx)\n', ...
    avgTable.RespT(end), avgTable.RespT(end)/exactR, ...
    avgTable.QLen(end), avgTable.QLen(end)/exactQ);

%% Block 5: the same answer from the api, which is all the solver does
% The solver is a thin wrapper over matlab/src/api/snc. An arrival envelope
% and a service envelope are function handles of the Chernoff parameter theta,
% and every bound is an infimum over theta of a closed-form expression. Note
% snc_srv_exp, not snc_srv_rate: the work unit here is the JOB, so the server
% is the counting process of an Exp(mu) service, and a constant-rate element
% would model an M/D/1 and understate the delay.
arv = @(theta) snc_env_poisson(lambda, theta);
srv = @(theta) snc_srv_exp(mu, theta);
[d, theta] = snc_perc_delay(arv, srv, 1e-3);
fprintf('\napi: d(1e-3) = %.4f at theta = %.4f\n', d, theta);
fprintf('     the optimal theta approaches log(mu/lambda) = %.4f, which is\n', log(mu/lambda));
fprintf('     what makes the backlog decay rate exact\n');
[epsAt, thetaAt] = snc_bound_delay(arv, srv, d);
fprintf('api: P{D > %.4f} <= %.3e  (theta = %.4f), exact tail %.3e\n', ...
    d, epsAt, thetaAt, exp(-(mu-lambda)*d));
