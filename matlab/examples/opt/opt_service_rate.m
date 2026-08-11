% opt_service_rate  Cheapest continuous service rate meeting a response-time
% SLA on an M/M/1 (theory rate 5.0). Mirrors opt_service_rate.py.

model = Network('MM1');
source = Source(model, 'Arrivals');
queue = Queue(model, 'Server', SchedStrategy.FCFS);
sink = Sink(model, 'Departures');
jobs = OpenClass(model, 'Jobs');
source.setArrival(jobs, Exp(3.0));
queue.setService(jobs, Exp(4.0));
model.link(Network.serialRouting(source, queue, sink));

problem = opt.OptimizationProblem(model);
problem.addVariable(opt.ServiceRate(queue, jobs, [3.5 8.0]));
rateCost = containers.Map('KeyType','char','ValueType','double');
rateCost('Server') = 20.0;
problem.setObjective(opt.MinimizeCost([], rateCost, [], {}));
problem.addConstraint(opt.ResponseTimeConstraint(queue, jobs, 0.5));

options = opt.LineOptSolverOptions(); options.setSeed(42); options.setMaxIterations(60);
result = problem.solve(options);

fprintf('Optimal rate    : %.3f (theory: 5.0)\n', result.getVariableValue('Server_Jobs_rate'));
fprintf('Cost            : %.1f\n', result.objectiveValue);
fprintf('Feasible        : %d\n', result.feasible);
