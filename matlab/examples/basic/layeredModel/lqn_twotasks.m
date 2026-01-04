clear solver AvgTable 

model = LayeredNetwork('myLayeredModel');

% definition of processors, tasks and entries
P1 = Processor(model, 'P1', 1, SchedStrategy.PS);
T1 = Task(model, 'T1', 100, SchedStrategy.REF).on(P1);
E1 = Entry(model, 'E1').on(T1);

P2 = Processor(model, 'P2', 1, SchedStrategy.PS);
T2 = Task(model, 'T2', 1, SchedStrategy.INF).on(P2);
E2 = Entry(model, 'E2').on(T2);

E3 = Entry(model, 'E3').on(T2);

% definition of activities
T1.setThinkTime(Erlang.fitMeanAndOrder(10,1));

A1 = Activity(model, 'A1', Exp(1)).on(T1).boundTo(E1).synchCall(E2).synchCall(E3,1);

A20 = Activity(model, 'A20', Exp(1)).on(T2).boundTo(E2);
A21 = Activity(model, 'A21', Exp(1)).on(T2);
A22 = Activity(model, 'A22', Exp(1)).on(T2).repliesTo(E2);
T2.addPrecedence(ActivityPrecedence.Serial(A20, A21, A22));

A3 = Activity(model, 'A3', Exp(1)).on(T2).boundTo(E3).repliesTo(E3);

% instantiate solvers
options = LQNS.defaultOptions;
options.keep = true;
options.verbose = 0;
%options.method = 'exact';
%options.method = 'lqsim';
%options.samples = 1e4;
solver{1} = LQNS(model, options);
AvgTable{1} = solver{1}.getAvgTable; %AvgTable{1}

lnoptions = LN.defaultOptions;
lnoptions.verbose = 0;
solveroptions = Solver.defaultOptions;
%options.method = 'comom';
solveroptions.verbose = 0;
solver{2} = LN(model, @(l)NC(l,solveroptions), lnoptions);
AvgTable{2} = solver{2}.getAvgTable; %AvgTable{2}
