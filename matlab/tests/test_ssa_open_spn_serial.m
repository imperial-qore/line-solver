%% OPEN stochastic Petri net on the SERIAL (Gillespie) SSA engine.
% Validates that Source -> Place -> Transition -> Sink open nets are simulated
% correctly by the serial afterEvent engine (options.method='serial'), against
% SolverJMT (ground truth) and the exact closed form. Before the fix the serial
% engine mis-measured Place markings (a Place is INF-scheduled, so toMarginal
% read only the server slot while the firing relocated tokens to the buffer
% slot, and consumption scaled by the enabling degree instead of the arc
% weight), producing garbage marking means.
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

seeds = [23000, 90111];

%% M/M/1-as-SPN: single-server transition, exact mean = rho/(1-rho).
lambda = 0.5; mu = 1.0; qexact = (lambda/mu)/(1-lambda/mu); % = 1.0
qs = zeros(1,numel(seeds)); ts = zeros(1,numel(seeds)); qj = 0;
for i = 1:numel(seeds)
    r = run_single(lambda, mu, 1, seeds(i));
    qs(i) = r.qSSA; ts(i) = r.tSSA; qj = r.qJMT;
end
qbar = mean(qs); tbar = mean(ts);
assert(abs(qbar - qexact)/qexact < 0.08, ...
    sprintf('M/M/1 SPN P1 mean tokens: serial SSA %.4f vs exact %.4f', qbar, qexact));
assert(abs(qbar - qj)/qj < 0.10, ...
    sprintf('M/M/1 SPN P1 mean tokens: serial SSA %.4f vs JMT %.4f', qbar, qj));
assert(abs(tbar - lambda)/lambda < 0.08, ...
    sprintf('M/M/1 SPN P1 throughput: serial SSA %.4f vs exact %.4f', tbar, lambda));

%% M/M/inf-as-SPN: infinite-server transition, exact mean = lambda/mu.
lambda = 2.0; mu = 1.0; qexact = lambda/mu; % = 2.0
qs = zeros(1,numel(seeds)); ts = zeros(1,numel(seeds)); qj = 0;
for i = 1:numel(seeds)
    r = run_single(lambda, mu, Inf, seeds(i));
    qs(i) = r.qSSA; ts(i) = r.tSSA; qj = r.qJMT;
end
qbar = mean(qs); tbar = mean(ts);
assert(abs(qbar - qexact)/qexact < 0.08, ...
    sprintf('M/M/inf SPN P1 mean tokens: serial SSA %.4f vs exact %.4f', qbar, qexact));
assert(abs(qbar - qj)/qj < 0.10, ...
    sprintf('M/M/inf SPN P1 mean tokens: serial SSA %.4f vs JMT %.4f', qbar, qj));
assert(abs(tbar - lambda)/lambda < 0.08, ...
    sprintf('M/M/inf SPN P1 throughput: serial SSA %.4f vs exact %.4f', tbar, lambda));

disp('test_ssa_open_spn_serial PASSED');

function q = run_single(lambda, mu, nserv, seed)
    model = Network('spn_open');
    source = Source(model,'Source');
    sink   = Sink(model,'Sink');
    P1 = Place(model,'P1');
    T1 = Transition(model,'T1');
    jobclass = OpenClass(model,'Class1',0);
    source.setArrival(jobclass, Exp(lambda));
    mode = T1.addMode('Mode1');
    T1.setNumberOfServers(mode, nserv);
    T1.setDistribution(mode, Exp(mu));
    T1.setEnablingConditions(mode, jobclass, P1, 1);
    T1.setFiringOutcome(mode, jobclass, sink, 1);
    R = model.initRoutingMatrix();
    R{1,1}(source,P1) = 1;
    R{1,1}(P1,T1) = 1;
    R{1,1}(T1,sink) = 1;
    model.link(R);

    opt = Solver.defaultOptions; opt.verbose = 0; opt.seed = seed;
    opt.cutoff = 40; opt.samples = 5e4;

    jmt = JMT(model, opt);
    tj = jmt.getAvgTable();

    opt2 = opt; opt2.method = 'serial';
    ssa = SolverSSA(model, opt2);
    ts = ssa.getAvgTable();

    % Assert the serial engine ran (never a silent NRM fallback).
    assert(strcmp(ssa.result.Avg.method, 'serial'), ...
        sprintf('expected serial method, got %s', ssa.result.Avg.method));

    % P1 is station row 2 (Source is row 1).
    q.qJMT = tj.QLen(2); q.qSSA = ts.QLen(2);
    q.tJMT = tj.Tput(2); q.tSSA = ts.Tput(2);
end
