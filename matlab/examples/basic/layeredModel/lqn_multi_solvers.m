clear solver AvgTable

model = LayeredNetwork('LQN1');

% definition of processors, tasks and entries
P1 = Processor(model, 'P1', Inf, SchedStrategy.INF);
T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1);
E1 = Entry(model, 'E1').on(T1);

P2 = Processor(model, 'P2', Inf, SchedStrategy.INF);
T2 = Task(model, 'T2', Inf, SchedStrategy.INF).on(P2);
E2 = Entry(model, 'E2').on(T2);

% definition of activities
T1.setThinkTime(Erlang.fitMeanAndOrder(0.0001,2));

A1 = Activity(model, 'A1', Exp(1.0)).on(T1).boundTo(E1).synchCall(E2,3);
A2 = Activity(model, 'A2', APH.fitMeanAndSCV(1,10)).on(T2).boundTo(E2).repliesTo(E2);

%%
% instantiate solvers
options = LQNS.defaultOptions;
options.keep = true;
options.verbose = 1;
%options.method = 'lqsim';
%options.samples = 1e4;
lqnssolver = LQNS(model, options);
AvgTableLQNS = lqnssolver.getAvgTable;
%AvgTableLQNS

% this method runs the MVA solver in each layer
lnoptions = LN.defaultOptions;
lnoptions.verbose = 0;
lnoptions.seed = 2300;
options = MVA.defaultOptions;
options.verbose = 0;
solver{1} = LN(model, @(model) MVA(model, options), lnoptions);
AvgTable{1} = solver{1}.getAvgTable;
%AvgTable{1}

% this method runs the NC solver in each layer
lnoptions = LN.defaultOptions;
lnoptions.verbose = 0;
lnoptions.seed = 2300;
options = NC.defaultOptions;
options.verbose = 0;
solver{2} = LN(model, @(model) NC(model, options), lnoptions);
AvgTable{2} = solver{2}.getAvgTable;
%AvgTable{2}

% this method adapts with the features of each layer
%solver{2} = LN(model, @(model) LINE(model, LINE.defaultOptions), lnoptions);
%AvgTable{2} = solver{2}.getAvgTable
%AvgTable{2}
