% test_spn_nrm_open
% Validates the SSA Next-Reaction-Method (NRM) stochastic-Petri-net path on OPEN
% nets: a Source feeds a Place whose tokens drain through a Transition to a Sink.
% Before the Source-arrival reaction was added the fed Place stayed empty and the
% run threw "Deadlock: no transition is enabled".
%
% Each net asserts (a) the solver actually ran method 'nrm' (never a silent
% serial fallback) and (b) the simulated marking mean and throughput match the
% analytic M/M/1 result. A Source Exp(lambda) feeding a single-server Transition
% Exp(mu) is an M/M/1 queue at the Place: mean tokens = rho/(1-rho), throughput =
% lambda (rho = lambda/mu < 1). The canonical net is cross-checked against JMT.
clear;

RTOL = 0.04;
SAMPLES = 3e5;
SEED = 23000;

%% --- Net 1: M/M/1 SPN, Source Exp(0.5) -> P1 -> T1 Exp(1.0) -> Sink -------
lambda = 0.5; mu = 1.0; rho = lambda/mu;
qExact = rho/(1-rho);   % = 1.0
model = mm1spn(lambda, mu);
solver = SolverSSA(model, 'method', 'nrm', 'samples', SAMPLES, 'seed', SEED);
[Qn,~,~,Tn] = solver.getAvg();
assert(~isempty(strfind(solver.result.Avg.method, 'nrm')), 'net1 did not run NRM');
% stations: Source(1), P1(2)
assert(abs(Qn(2)-qExact)/qExact < RTOL, sprintf('net1 P1 tokens %g vs exact %g', Qn(2), qExact));
assert(abs(Tn(1)-lambda)/lambda < RTOL, sprintf('net1 Source tput %g vs %g', Tn(1), lambda));
assert(abs(Tn(2)-lambda)/lambda < RTOL, sprintf('net1 P1 tput %g vs %g', Tn(2), lambda));

% Cross-check mean tokens against JMT's simulation of the same net.
jmt = SolverJMT(mm1spn(lambda, mu), 'samples', SAMPLES, 'seed', SEED);
[Qj,~,~,~] = jmt.getAvg();
assert(abs(Qn(2)-Qj(2))/max(Qj(2),1e-9) < RTOL, sprintf('net1 P1 tokens NRM %g vs JMT %g', Qn(2), Qj(2)));

%% --- Net 2: open tandem, two places in series ----------------------------
lambda = 0.5; mu1 = 1.0; mu2 = 2.0;
q1 = (lambda/mu1)/(1-lambda/mu1);   % = 1.0
q2 = (lambda/mu2)/(1-lambda/mu2);   % = 1/3
model = tandemspn(lambda, mu1, mu2);
solver = SolverSSA(model, 'method', 'nrm', 'samples', SAMPLES, 'seed', SEED);
[Qn,~,~,Tn] = solver.getAvg();
assert(~isempty(strfind(solver.result.Avg.method, 'nrm')), 'net2 did not run NRM');
% stations: Source(1), P1(2), P2(3)
assert(abs(Qn(2)-q1)/q1 < RTOL, sprintf('net2 P1 tokens %g vs exact %g', Qn(2), q1));
assert(abs(Qn(3)-q2)/q2 < RTOL, sprintf('net2 P2 tokens %g vs exact %g', Qn(3), q2));
assert(abs(Tn(2)-lambda)/lambda < RTOL, sprintf('net2 P1 tput %g vs %g', Tn(2), lambda));
assert(abs(Tn(3)-lambda)/lambda < RTOL, sprintf('net2 P2 tput %g vs %g', Tn(3), lambda));

disp('test_spn_nrm_open passed');

%% --- model builders ------------------------------------------------------
function model = mm1spn(lambda, mu)
model = Network('mm1spn');
source = Source(model, 'Source'); sink = Sink(model, 'Sink');
P1 = Place(model, 'P1'); T1 = Transition(model, 'T1');
jobclass = OpenClass(model, 'Class1', 0);
source.setArrival(jobclass, Exp(lambda));
mode = T1.addMode('Mode1'); T1.setDistribution(mode, Exp(mu));
T1.setEnablingConditions(mode, jobclass, P1, 1);
T1.setFiringOutcome(mode, jobclass, sink, 1);
R = model.initRoutingMatrix();
R{1,1}(source,P1) = 1; R{1,1}(P1,T1) = 1; R{1,1}(T1,sink) = 1;
model.link(R);
end

function model = tandemspn(lambda, mu1, mu2)
model = Network('tandemspn');
source = Source(model, 'Source'); sink = Sink(model, 'Sink');
P1 = Place(model, 'P1'); P2 = Place(model, 'P2');
T1 = Transition(model, 'T1'); T2 = Transition(model, 'T2');
jobclass = OpenClass(model, 'Class1', 0);
source.setArrival(jobclass, Exp(lambda));
m1 = T1.addMode('Mode1'); T1.setDistribution(m1, Exp(mu1));
T1.setEnablingConditions(m1, jobclass, P1, 1);
T1.setFiringOutcome(m1, jobclass, P2, 1);
m2 = T2.addMode('Mode1'); T2.setDistribution(m2, Exp(mu2));
T2.setEnablingConditions(m2, jobclass, P2, 1);
T2.setFiringOutcome(m2, jobclass, sink, 1);
R = model.initRoutingMatrix();
R{1,1}(source,P1) = 1; R{1,1}(P1,T1) = 1;
R{1,1}(T1,P2) = 1; R{1,1}(P2,T2) = 1; R{1,1}(T2,sink) = 1;
model.link(R);
end
