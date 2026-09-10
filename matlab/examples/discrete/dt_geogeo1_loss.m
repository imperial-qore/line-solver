% Geo/Geo/1/L: the discrete-time loss system of Daduna (2001), corollary 2.8.
%
% The buffer holds at most L jobs. An arrival that lands in a slot which
% already holds L jobs is lost, i.e. b(n) = 0 for n >= L, which is exactly the
% assumption corollary 2.8 places on the arrival probabilities. The queue
% length law stays the birth-death form of theorem 2.3, now normalized over
% the finite state space, and the loss probability is the stationary
% probability of a full system.

clear solver AvgTable;

a = 0.2;    % per-slot arrival probability
s = 0.5;    % per-slot service completion probability
L = 4;      % buffer capacity in jobs

model = Network('GeoGeo1L');

source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

jobClass = OpenClass(model, 'Class1');
source.setArrival(jobClass, Geometric(a));
queue.setService(jobClass, Geometric(s));
queue.setCapacity(L);

model.link(Network.serialRouting(source, queue, sink));

%%
options = SolverNC.defaultOptions;
options.config.slotted = true;
solver = SolverNC(model, options);
AvgTable = solver.getAvgTable; AvgTable

r = dqsys_bernoulli1(a, s, L);
fprintf('queue length law   : %s\n', mat2str(round(r.pmf, 6)));
fprintf('loss probability   : %g\n', r.lossProb);
fprintf('carried throughput : %g of the %g offered per slot\n', r.throughput, a);
