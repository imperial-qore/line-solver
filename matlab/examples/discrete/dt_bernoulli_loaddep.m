% Load-dependent Bernoulli server, and the discrete-time arrival theorem.
%
% Example 2.10 of Daduna (2001) notes that a discrete-time M/M/c queue has no
% exactly equivalent state dependent single server, but that p(n) = p min(n,c)
% reproduces its conditional service intensity. That is a load dependence in
% LINE terms, so the model is a single Bernoulli server whose service
% probability rises with the queue length up to c servers' worth of capacity.
%
% The example also prints the arrival distribution of theorem 2.11, the law an
% arriving job sees with itself not counted. In continuous time with Poisson
% arrivals that law would coincide with the time-stationary one (PASTA);
% discrete time has no such analogue, and the two rows below differ.

clear solver AvgTable;

a = 0.6;    % per-slot arrival probability
s = 0.3;    % per-slot service completion probability of one server
c = 3;      % servers' worth of capacity
L = 20;     % buffer capacity, which bounds the state space

model = Network('LoadDepBernoulli');

source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

jobClass = OpenClass(model, 'Class1');
source.setArrival(jobClass, Geometric(a));
queue.setService(jobClass, Geometric(s));
queue.setCapacity(L);
queue.setLoadDependence(min(1:L, c));   % p(n) = s * min(n,c), example 2.10

model.link(Network.serialRouting(source, queue, sink));

%%
options = SolverNC.defaultOptions;
options.config.slotted = true;
solver = SolverNC(model, options);
AvgTable = solver.getAvgTable; AvgTable

r = dqsys_bernoulli1(a, s * min(1:L, c), L);
fprintf('time-stationary law pi(0..6)  : %s\n', mat2str(round(r.pmf(1:7), 6)));
fprintf('arrival law         pi_1(0..6): %s\n', mat2str(round(r.arrivalPmf(1:7), 6)));
fprintf('no PASTA in discrete time: the two rows above are different laws\n');
