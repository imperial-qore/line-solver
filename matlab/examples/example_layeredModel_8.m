% example of layered model with a function task

clear solver AvgTable

model = LayeredNetwork('faas_test_example');

% definition of processors, tasks and entries
P1 = Processor(model, 'P1', Inf, SchedStrategy.INF);
T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1);
E1 = Entry(model, 'E1').on(T1);

P2 = Processor(model, 'P2', 4, SchedStrategy.FCFS);
% T2 = Task(model, 'T2', 1, SchedStrategy.FCFS).on(P2);
T2 = FunctionTask(model, 'F2', 6, SchedStrategy.FCFS).on(P2).setThinkTime(Exp.fitMean(8.0));
T2.setSetupTime(Exp(1.0));
T2.setDelayOffTime(Exp(2.0));

E2 = Entry(model, 'E2').on(T2);

% T3 = Task(model, 'T3', 1, SchedStrategy.FCFS).on(P2);
% E3 = Entry(model, 'E3').on(T3);

% definition of activities
A1 = Activity(model, 'A1', Exp(1.0)).on(T1).boundTo(E1).synchCall(E2,1);
A2 = Activity(model, 'A2', Exp(3.0)).on(T2).boundTo(E2).repliesTo(E2);
% A3 = Activity(model, 'A3', Exp(5.0)).on(T3).boundTo(E3).repliesTo(E3);


lnoptions = SolverLN.defaultOptions;
lnoptions.verbose = 0;
lnoptions.seed = 2300;
% options = SolverMVA.defaultOptions;
% options.verbose = 0;
% options2 = SolverMAM.defaultOptions;
% options2.verbose = 0;
solver{1} = SolverLN(model, @(model) myFunction(model), lnoptions);
AvgTable{1} = solver{1}.getAvgTable
AvgTable{1}