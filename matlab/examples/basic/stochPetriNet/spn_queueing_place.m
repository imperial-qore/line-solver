% Closed queueing Petri net (QPN) with two queueing places.
%
% A queueing place embeds a scheduling station inside a Petri-net place: an
% arriving token is served by the place's embedded queue and, on completion,
% moves to a depository from which the output transitions consume it.
%
% Here a CPU (single-server FCFS queueing place) and a think stage (infinite-
% server queueing place) exchange a fixed population of N tokens through two
% immediate transitions. This is the queueing-Petri-net rendering of a machine-
% repairman / finite-population model, so its exact solution is available from
% SolverMVA on the equivalent Delay+Queue network for cross-validation.

N = 4;
model = Network('QueueingPetriNet');

% Two queueing places: setService turns each Place into a queueing place whose
% embedded queue is served under the scheduling strategy of its constructor.
cpu   = Place(model, 'CPU',   SchedStrategy.FCFS);  % single-server embedded queue
think = Place(model, 'Think', SchedStrategy.INF);   % infinite-server embedded queue

jobs = ClosedClass(model, 'Jobs', N, think, 0);
cpu.setService(jobs, Exp(1.5));    % embedded FCFS service, rate 1.5
think.setService(jobs, Exp(0.5));  % embedded think time, mean 2

% Two immediate transitions cycle one token per firing between the places.
toCPU = Transition(model, 'toCPU');
m1 = toCPU.addMode('m1');
toCPU.setTimingStrategy(m1, TimingStrategy.IMMEDIATE);
toCPU.setEnablingConditions(m1, jobs, think, 1);
toCPU.setFiringOutcome(m1, jobs, cpu, 1);

toThink = Transition(model, 'toThink');
m2 = toThink.addMode('m2');
toThink.setTimingStrategy(m2, TimingStrategy.IMMEDIATE);
toThink.setEnablingConditions(m2, jobs, cpu, 1);
toThink.setFiringOutcome(m2, jobs, think, 1);

P = model.initRoutingMatrix();
P.set(jobs, jobs, think, toCPU,   1.0);
P.set(jobs, jobs, toCPU, cpu,     1.0);
P.set(jobs, jobs, cpu,   toThink, 1.0);
P.set(jobs, jobs, toThink, think, 1.0);
model.link(P);

% Initial marking: all tokens in the think place.
think.setState(N);
cpu.setState(0);

% Queueing places are simulated by LDES.
AvgTableLDES = SolverLDES(model, 'seed', 23000, 'samples', 2e5).avgTable()

% Exact cross-check: the equivalent finite-population Delay + M/M/1 network.
ref = Network('ref');
delay = Delay(ref, 'Think');
queue = Queue(ref, 'CPU', SchedStrategy.FCFS);
cc = ClosedClass(ref, 'Jobs', N, delay, 0);
delay.setService(cc, Exp(0.5));
queue.setService(cc, Exp(1.5));
ref.link(Network.serialRouting(delay, queue));
AvgTableMVA = SolverMVA(ref).avgTable()
