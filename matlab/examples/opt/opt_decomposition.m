% opt_decomposition  Two-variable problem (server count + service rate) solved
% by a DecompositionWorkflow: auto-decompose by variable type, then Gauss-Seidel
% cycling to a fixed point. Mirrors opt_decomposition.py.

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
problem.addVariable(opt.ServiceRate(queue, jobs, [1.0 4.0]));
serverCost = configureDictionary('string','double'); serverCost('Server') = 10.0;
rateCost = configureDictionary('string','double'); rateCost('Server') = 20.0;
problem.setObjective(opt.MinimizeCost(serverCost, rateCost, [], {}));
problem.addConstraint(opt.ResponseTimeConstraint(queue, jobs, 0.5));

workflow = problem.decompose();
workflow.autoDecompose();
options = opt.LineOptSolverOptions(); options.setSeed(42);
workflow.setSolverOptions(options);
result = workflow.solveSequential(4, 1e-3);

fprintf('Converged       : %d after %d cycle(s)\n', result.converged, result.cyclesCompleted);
fprintf('Servers         : %d\n', result.getFinalVariableValue('Server_servers'));
fprintf('Rate            : %.3f\n', result.getFinalVariableValue('Server_Jobs_rate'));
fprintf('Objective       : %.2f\n', result.finalObjective);
