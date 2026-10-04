%% LDES warm start from an auxiliary solver
% An auxiliary solver is passed to SolverLDES as an argument, its steady-state
% distribution is computed, and that distribution decides the initial state of
% the simulation. Since the simulation then starts (approximately) in steady
% state, the initialization bias vanishes and no warmup samples are discarded,
% so a target accuracy is reached with fewer simulated events (and less time)
% than a cold-started run.
%
% Auxiliary-solver dispatch inside SolverLDES.initFromSolver:
%  - SolverCTMC: the exact stationary distribution over the aggregate state
%    space is computed and the initial state is its mode;
%  - any other solver (e.g. SolverMVA): the steady-state mean queue lengths
%    are rounded to a placement conserving the closed populations.
%
% The benchmark model is a closed NEAR-BALANCED tandem, whose job split mixes
% slowly: the bias of the default cold start (all jobs at the reference
% station) persists beyond what the MSER-5 transient filter can remove.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

N = 100;
seeds = [23000 23001 23002];
budgets = [10000 50000 200000];

% Exact reference solution (exact MVA; the model is product-form)
exact = SolverMVA(buildModel(N), 'method', 'exact').getAvgQLen();
fprintf('Exact mean queue lengths: %s\n', mat2str(round(exact', 2)));

% Warm-start placement, computed ONCE by passing the auxiliary solver to
% SolverLDES; the derived options.init_sol is then reused across runs.
T0 = tic;
m = buildModel(N);
proto = SolverLDES(m, SolverMVA(m, 'method', 'exact'));
initSol = proto.options.init_sol;
initTime = toc(T0);
fprintf('Warm placement from SolverMVA (%.3fs): %s\n\n', initTime, mat2str(initSol));

% On a small instance, SolverCTMC yields the mode of the exact stationary
% distribution instead (feasible when the state space is small):
mS = buildModel(20);
protoCtmc = SolverLDES(mS, SolverCTMC(mS,'exact'));
fprintf('CTMC distribution-mode placement (N=20): %s\n\n', mat2str(protoCtmc.options.init_sol));

% Compare accuracy and wall-clock time of cold vs warm starts
fprintf('%10s  %14s  %14s\n', 'samples', 'COLD err/time', 'WARM err/time');
coldTime = 0; warmTime = initTime;
coldAt = NaN; warmAt = NaN;
for b = budgets
    [ce, ct] = runSet(@() buildModel(N), exact, b, seeds, []);
    [we, wt] = runSet(@() buildModel(N), exact, b, seeds, initSol);
    coldTime = coldTime + ct; warmTime = warmTime + wt;
    if isnan(coldAt) && ce < 0.10, coldAt = b; end
    if isnan(warmAt) && we < 0.10, warmAt = b; end
    fprintf('%10d  %6.2f%% %5.2fs  %6.2f%% %5.2fs\n', b, 100*ce, ct, 100*we, wt);
    if ~isnan(coldAt) && ~isnan(warmAt)
        break
    end
end
fprintf('\nSamples to reach 10%% error: COLD=%d, WARM=%d\n', coldAt, warmAt);

function [err, rtime] = runSet(factory, exact, samples, seeds, initSol)
% Average L1 relative error and total wall-clock over the seed set.
err = 0; rtime = 0;
for s = seeds
    model = factory();
    T0 = tic;
    if isempty(initSol)
        solver = SolverLDES(model, 'samples', samples, 'seed', s, 'verbose', false);
    else
        solver = SolverLDES(model, 'samples', samples, 'seed', s, 'verbose', false);
        solver.options.init_sol = initSol;
        solver.options.config.tranfilter = 'fixed';
        solver.options.config.warmupfrac = 0;
    end
    QN = solver.getAvgQLen();
    rtime = rtime + toc(T0);
    err = err + sum(abs(QN - exact)) / sum(exact);
end
err = err / length(seeds);
end

function model = buildModel(n)
% Closed near-balanced tandem: Think(Exp,1) -> Queue1(Exp,1.0) -> Queue2(Exp,0.98)
model = Network('ldesWarmStart');
think = Delay(model, 'Think');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
jobs = ClosedClass(model, 'Jobs', n, think);
think.setService(jobs, Exp(1.0));
queue1.setService(jobs, Exp(1.0));
queue2.setService(jobs, Exp(0.98));  % near-balanced bottleneck
P = model.initRoutingMatrix();
P{1}(think, queue1) = 1.0;
P{1}(queue1, queue2) = 1.0;
P{1}(queue2, think) = 1.0;
model.link(P);
end
