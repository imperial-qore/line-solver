% opt_server_sizing  Minimum servers of an M/M/c queue under a 50% utilization
% SLA, via differential evolution. lambda=3, mu=1 -> optimum c*=6 (cost 60).
% Mirrors python/examples/opt/opt_server_sizing.py.

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
serverCost = configureDictionary('string','double');
serverCost('Server') = 10.0;
problem.setObjective(opt.MinimizeCost(serverCost, [], [], {}));
problem.addConstraint(opt.UtilizationConstraint(queue, 0.5));

options = opt.LineOptSolverOptions(); options.setSeed(42);
result = problem.solve(options);

fprintf('Optimal servers : %d\n', result.getVariableValue('Server_servers'));
fprintf('Cost            : %.1f\n', result.objectiveValue);
fprintf('Feasible        : %d\n', result.feasible);
fprintf('LINE evaluations: %d\n', result.modelEvaluations);
