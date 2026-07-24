% opt_pareto_frontier  Cost vs utilization-bound frontier via ParetoSweep
% (epsilon-constraint method, bisection per point). Mirrors
% opt_pareto_frontier.py.

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

factory = @(eps) opt.UtilizationConstraint(queue, eps);
sweep = opt.ParetoSweep(problem, factory, [0.3 0.4 0.5 0.6 0.75 0.9], 'bisection');
sweep.solve();

fprintf('Utilization bound -> optimal cost (servers)\n');
frontier = sweep.getFrontier();
for i = 1:numel(frontier)
    p = frontier{i};
    fprintf('  util <= %-4g -> %6.1f  (%d)\n', p.epsilon, p.objectiveValue, ...
        p.result.getVariableValue('Server_servers'));
end
