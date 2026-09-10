% opt_load_balancing  Routing split from a source to a fast and a slow PS
% server minimizing end-to-end response time (theory fast fraction ~0.743).
% Mirrors opt_load_balancing.py.

model = Network('LoadBalance');
source = Source(model, 'S');
fast = Queue(model, 'Fast', SchedStrategy.PS);
slow = Queue(model, 'Slow', SchedStrategy.PS);
sink = Sink(model, 'K');
jobs = OpenClass(model, 'Jobs');
source.setArrival(jobs, Exp(2.0));
fast.setService(jobs, Exp(4.0));
slow.setService(jobs, Exp(2.5));
model.addLink(source, fast);
model.addLink(source, slow);
model.addLink(fast, sink);
model.addLink(slow, sink);

problem = opt.OptimizationProblem(model);
problem.addVariable(opt.RoutingProbabilities(jobs, source, {fast, slow}));
problem.setObjective(opt.MinimizeSystemResponseTime(jobs));

options = opt.LineOptSolverOptions(); options.setSeed(42); options.setMaxIterations(60);
result = problem.solve(options);

probs = result.getVariableValue('Jobs_routing_from_S');
fprintf('Fraction to Fast: %.3f (theory: 0.743)\n', probs(1));
fprintf('Fraction to Slow: %.3f\n', probs(2));
