% Stochastic network calculus on a feed-forward network: envelope propagation,
% blind multiplexing, and the model classes the family refuses.
%
% The bound of snc_delay_quantile.m is a single station. This example is what
% the network calculus adds on top of it: a departure envelope carries a flow
% to the next hop, cross traffic is subtracted from a shared server, and both
% operations cost burstiness, which is visible as the bound loosening downstream
% and under sharing.
%
% The reference throughout is the exact Jackson-network solution, since every
% station here is an M/M/1 with Poisson input by Burke's theorem.

%% Block 1: a three-station tandem
lambda = 0.6;
rates = [1.5 1.2 1.0];              % service rates, the last is the bottleneck
model = Network('SncTandem');
source = Source(model,'Source');
q1 = Queue(model,'Q1', SchedStrategy.FCFS);
q2 = Queue(model,'Q2', SchedStrategy.FCFS);
q3 = Queue(model,'Q3', SchedStrategy.FCFS);
sink = Sink(model,'Sink');
jobclass = OpenClass(model,'Class1');
source.setArrival(jobclass, Exp(lambda));
q1.setService(jobclass, Exp(rates(1)));
q2.setService(jobclass, Exp(rates(2)));
q3.setService(jobclass, Exp(rates(3)));
model.link(Network.serialRouting(source,q1,q2,q3,sink));

tandemTable = SolverBA(model,'method','snc.upper').getAvgTable()

%% Block 2: what the propagation costs
% Each hop replaces the arrival envelope by the DEPARTURE envelope of the
% station upstream, which carries the burst the server has added. The exact
% answer does not degrade this way -- by Burke's theorem the departure process
% of an M/M/1 is again Poisson -- so the ratio to the exact response time grows
% hop by hop. The bound stays valid; it is the price of assuming nothing about
% the departure process beyond its envelope.
fprintf('\n%-6s %10s %10s %8s\n','station','R bound','R exact','ratio');
for i = 1:3
    Rexact = 1/(rates(i)-lambda);
    Rbound = tandemTable.RespT(i+1);
    fprintf('%-6s %10.4f %10.4f %8.2f\n', ...
        char(tandemTable.Station(i+1)), Rbound, Rexact, Rbound/Rexact);
end

%% Block 3: end-to-end, where pay-bursts-only-once pays
% Summing the per-station bounds pays the burst term at every hop. Concatenating
% the three service envelopes with snc_conv first and bounding the composed
% element once pays it only once, which is the classical result of the network
% calculus and is worth several tens of percent here.
arv = @(theta) snc_env_poisson(lambda, theta);
srvEnd = @(theta) tandemService(theta, rates);
[dEnd, thetaEnd] = snc_perc_delay(arv, srvEnd, 1e-3);
dHop = 0;
for i = 1:3
    dHop = dHop + snc_perc_delay(arv, @(theta) snc_srv_exp(rates(i), theta), 1e-3/3);
end
fprintf('\nend-to-end delay quantile at eps=1e-3\n');
fprintf('  concatenated (snc_conv) : %8.4f  at theta = %.4f\n', dEnd, thetaEnd);
fprintf('  summed per hop          : %8.4f\n', dHop);
fprintf('  pay bursts once saves   : %7.1f%%\n', 100*(1-dEnd/dHop));

%% Block 4: two classes sharing one server (blind multiplexing)
% A class sharing a station sees the server minus whatever the other classes
% take from it: snc_leftover subtracts the cross-flow arrival envelope from the
% service envelope. The result holds for ANY work-conserving discipline at that
% station, which is why it is well above the FCFS answer -- it also covers the
% policy that serves the other class first whenever it can.
shared = Network('SncShared');
src2 = Source(shared,'Source');
qs = Queue(shared,'Shared', SchedStrategy.FCFS);
snk2 = Sink(shared,'Sink');
classA = OpenClass(shared,'ClassA');
classB = OpenClass(shared,'ClassB');
src2.setArrival(classA, Exp(0.3));
src2.setArrival(classB, Exp(0.3));
qs.setService(classA, Exp(1));
qs.setService(classB, Exp(1));
P = shared.initRoutingMatrix();
P{classA,classA}(src2,qs) = 1; P{classA,classA}(qs,snk2) = 1;
P{classB,classB}(src2,qs) = 1; P{classB,classB}(qs,snk2) = 1;
shared.link(P);

sharedSolver = SolverBA(shared,'method','snc.upper');
sharedTable = sharedSolver.getAvgTable()
sharedPerc = sharedSolver.getPercTable(1e-3)
fprintf('exact per-class response time (aggregate M/M/1, lambda=0.6, mu=1): %.4f\n', ...
    1/(1-0.6));

%% Block 5: what the family refuses, and why
% The elementary envelope algebra has real limits, and the analyzer states them
% rather than returning a number that looks plausible. Each refusal below is a
% modelling assumption of the calculus, not an implementation gap.
fprintf('\nrefusals:\n');
showRefusal('closed network', @() closedModel());
showRefusal('probabilistic split downstream of the Source', @() splitModel());
showRefusal('unequal service rates at a shared station', @() unequalModel());

% ----- local functions -----

function [sigma, rho] = tandemService(theta, rates)
% Min-plus concatenation of the three service envelopes into one element.
[sigma, rho] = snc_srv_exp(rates(1), theta);
for i = 2:numel(rates)
    [s2, r2] = snc_srv_exp(rates(i), theta);
    [sigma, rho] = snc_conv(sigma, rho, s2, r2, theta);
end
end

function showRefusal(what, fun)
try
    fun();
    fprintf('  %-45s NO REFUSAL (unexpected)\n', what);
catch ME
    msg = ME.message;
    at = strfind(msg,']');
    if ~isempty(at)
        msg = strtrim(msg(at(1)+1:end));
    end
    fprintf('  %-45s %s\n', [what ':'], msg);
end
end

function closedModel()
model = Network('Closed');
delay = Delay(model,'Think');
queue = Queue(model,'Q', SchedStrategy.PS);
jobs = ClosedClass(model,'C', 3, delay);
delay.setService(jobs, Exp(1));
queue.setService(jobs, Exp(2));
model.link(Network.serialRouting(delay,queue));
SolverBA(model,'method','snc.upper').getAvgTable();
end

function splitModel()
model = Network('Split');
source = Source(model,'Source');
qa = Queue(model,'QA', SchedStrategy.FCFS);
qb = Queue(model,'QB', SchedStrategy.FCFS);
sink = Sink(model,'Sink');
jobclass = OpenClass(model,'C');
source.setArrival(jobclass, Exp(0.4));
qa.setService(jobclass, Exp(1));
qb.setService(jobclass, Exp(1));
P = model.initRoutingMatrix();
P{jobclass,jobclass}(source,qa) = 1;
P{jobclass,jobclass}(qa,qb) = 0.5;
P{jobclass,jobclass}(qa,sink) = 0.5;
P{jobclass,jobclass}(qb,sink) = 1;
model.link(P);
SolverBA(model,'method','snc.upper').getAvgTable();
end

function unequalModel()
model = Network('Unequal');
source = Source(model,'Source');
queue = Queue(model,'Q', SchedStrategy.FCFS);
sink = Sink(model,'Sink');
classA = OpenClass(model,'A');
classB = OpenClass(model,'B');
source.setArrival(classA, Exp(0.2));
source.setArrival(classB, Exp(0.2));
queue.setService(classA, Exp(1));
queue.setService(classB, Exp(2));
P = model.initRoutingMatrix();
P{classA,classA}(source,queue) = 1; P{classA,classA}(queue,sink) = 1;
P{classB,classB}(source,queue) = 1; P{classB,classB}(queue,sink) = 1;
model.link(P);
SolverBA(model,'method','snc.upper').getAvgTable();
end
