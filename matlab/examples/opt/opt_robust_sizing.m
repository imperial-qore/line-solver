% opt_robust_sizing  Worst-case server sizing across an average-load and a
% peak-load workload scenario (constraints enforced on all scenarios). Average
% load needs 6 servers; the peak scenario pushes it to 9. Mirrors
% opt_robust_sizing.py.

base = buildMMc(3.0);      % average load
peak = buildMMc(4.5);      % peak load scenario
nodes = base.getNodes(); queue = nodes{2};

problem = opt.OptimizationProblem(base);
problem.addVariable(opt.ServerAllocation(queue, [1 12]));
serverCost = configureDictionary('string','double');
serverCost('Server') = 10.0;
problem.setObjective(opt.MinimizeCost(serverCost, [], [], {}));
problem.addConstraint(opt.UtilizationConstraint(queue, 0.5));
problem.addScenario(peak, 1.0);

result = opt.BisectionSolver(problem).solve();

fprintf('Robust servers  : %d (6 for average load only)\n', result.getVariableValue('Server_servers'));
fprintf('Cost            : %.1f\n', result.objectiveValue);
fprintf('Feasible        : %d\n', result.feasible);

function model = buildMMc(arrivalRate)
    model = Network('MMc');
    source = Source(model, 'Arrivals');
    queue = Queue(model, 'Server', SchedStrategy.FCFS);
    sink = Sink(model, 'Departures');
    jobs = OpenClass(model, 'Jobs');
    source.setArrival(jobs, Exp(arrivalRate));
    queue.setService(jobs, Exp(1.0));
    model.link(Network.serialRouting(source, queue, sink));
end
