% opt_lqn_host_demand  Tune the host demands of a two-layer LQN under an
% end-to-end response-time budget, using the partial-sensitivity gradient.
%
% lqnGradient='partial_plus_fd' drives the analytic path: SolverLN's per-layer
% d(metric)/d(service rate) table is reshaped into a sensitivity dictionary and
% chained with d(rate)/d(demand) = -1/D^2, with a periodic whole-model finite
% difference (every fdRefresh gradient calls) correcting the bias.

model = LayeredNetwork('BookstoreLQN');

P1 = Processor(model, 'ClientCPU', 1, SchedStrategy.INF);
T1 = Task(model, 'Client', 4, SchedStrategy.REF).on(P1);
T1.setThinkTime(Exp(1.0));
E1 = Entry(model, 'Browse').on(T1);

P2 = Processor(model, 'AppCPU', 1, SchedStrategy.PS);
T2 = Task(model, 'AppServer', 1, SchedStrategy.FCFS).on(P2);
E2 = Entry(model, 'Render').on(T2);

P3 = Processor(model, 'DbCPU', 1, SchedStrategy.PS);
T3 = Task(model, 'Database', 1, SchedStrategy.FCFS).on(P3);
E3 = Entry(model, 'Query').on(T3);

A1 = Activity(model, 'BrowseAct', Exp(2.0)).on(T1).boundTo(E1).synchCall(E2, 1);
A2 = Activity(model, 'RenderAct', Exp(4.0)).on(T2).boundTo(E2).synchCall(E3, 1);
A2.repliesTo(E2);
A3 = Activity(model, 'QueryAct', Exp(5.0)).on(T3).boundTo(E3);
A3.repliesTo(E3);

problem = opt.OptimizationProblem(model);
problem.addVariable(opt.HostDemand('RenderAct', [0.05 0.50]));
problem.addVariable(opt.HostDemand('QueryAct', [0.05 0.50]));
problem.setObjective(opt.MinimizeSystemResponseTime('Client', {}));

% Faster servers cost more: budget is on the reciprocal-demand scale.
demandCost = configureDictionary('string','double');
demandCost('RenderAct_hostdemand') = 1.0;
demandCost('QueryAct_hostdemand') = 1.0;
problem.addConstraint(opt.BudgetConstraint(0.45, demandCost));

options = opt.LineOptSolverOptions();
options.setSeed(7);
options.setOptimizer('gradient');
options.setLqnGradient('partial_plus_fd');
options.setFdRefresh(3);
options.setMaxIterations(8);
result = problem.solve(options);

fprintf('Render demand   : %.4f\n', result.getVariableValue('RenderAct_hostdemand'));
fprintf('Query demand    : %.4f\n', result.getVariableValue('QueryAct_hostdemand'));
fprintf('System RespT    : %.4f\n', result.objectiveValue);
fprintf('Feasible        : %d\n', result.feasible);
fprintf('Total violation : %.4f\n', result.getTotalViolation());
fprintf('LINE evaluations: %d\n', result.modelEvaluations);
