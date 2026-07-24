% opt_bisection_sizing  Exact server sizing via BisectionSolver (O(log n) LINE
% solves). Same M/M/c as opt_server_sizing -> 6 servers, cost 60, 4 probes.
% Mirrors python/examples/opt/opt_bisection_sizing.py.

model = Network('MMc');
source = Source(model, 'Arrivals');
queue = Queue(model, 'Server', SchedStrategy.FCFS);
sink = Sink(model, 'Departures');
jobs = OpenClass(model, 'Jobs');
source.setArrival(jobs, Exp(3.0));
queue.setService(jobs, Exp(1.0));
model.link(Network.serialRouting(source, queue, sink));

problem = opt.OptimizationProblem(model);
problem.addVariable(opt.ServerAllocation(queue, [1 10]));
serverCost = containers.Map('KeyType','char','ValueType','double');
serverCost('Server') = 10.0;
problem.setObjective(opt.MinimizeCost(serverCost, [], [], {}));
problem.addConstraint(opt.UtilizationConstraint(queue, 0.5));

result = opt.BisectionSolver(problem).solve();

fprintf('Optimal servers : %d\n', result.getVariableValue('Server_servers'));
fprintf('Cost            : %.1f\n', result.objectiveValue);
fprintf('LINE evaluations: %d (vs hundreds for DE)\n', result.modelEvaluations);
