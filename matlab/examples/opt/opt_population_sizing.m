% opt_population_sizing  Largest closed-class population sustaining an
% interactive response-time SLA, via BisectionSolver('max_feasible').
% Interactive system: Delay (think) + PS server. Mirrors opt_population_sizing.py.

model = Network('Interactive');
think = Delay(model, 'Think');
server = Queue(model, 'AppServer', SchedStrategy.PS);
users = ClosedClass(model, 'Users', 1, think);
think.setService(users, Exp(1.0));
server.setService(users, Exp(2.0));
model.link(Network.serialRouting(think, server));

problem = opt.OptimizationProblem(model);
problem.addVariable(opt.JobPopulation(users, [1 50]));
problem.setObjective(opt.MinimizeCost([], [], [], {}));   % pure feasibility sizing
problem.addConstraint(opt.ResponseTimeConstraint(server, users, 2.0));

result = opt.BisectionSolver(problem, 'max_feasible').solve();

fprintf('Max users       : %d\n', result.getVariableValue('Users_population'));
fprintf('Feasible        : %d\n', result.feasible);
fprintf('LINE evaluations: %d\n', result.modelEvaluations);
