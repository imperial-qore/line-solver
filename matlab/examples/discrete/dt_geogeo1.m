% Geo/Geo/1: the discrete-time single-server queue with an unbounded buffer.
%
% Arrivals occur with probability a in each slot, a service completes with
% probability s. SolverNC recognizes the slotted model and returns the exact
% closed form of Daduna (2001), theorem 2.3, which for constant a and s is
% the geometric law of corollary 2.7. Mean queue length a(1-a)/(s-a) and mean
% sojourn time (1-a)/(s-a) slots.

clear solver AvgTable;

a = 0.2;    % per-slot arrival probability
s = 0.5;    % per-slot service completion probability

model = Network('GeoGeo1');

source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

jobClass = OpenClass(model, 'Class1');
source.setArrival(jobClass, Geometric(a));
queue.setService(jobClass, Geometric(s));

model.link(Network.serialRouting(source, queue, sink));

%%
options = SolverNC.defaultOptions;
options.config.slotted = true;      % run on the slot lattice
solver = SolverNC(model, options);
AvgTable = solver.getAvgTable; AvgTable

fprintf('closed form: E[N] = %g, E[T] = %g slots\n', ...
    a*(1-a)/(s-a), (1-a)/(s-a));

% The same numbers straight from the single-queue formula.
r = dqsys_bernoulli1(a, s)
