model = LayeredNetwork('myLayeredModel');

P{1} = Processor(model, 'P1', 1, SchedStrategy.INF);
P{2} = Processor(model, 'P2', 1, SchedStrategy.PS);

T{1} = Task(model, 'T1', 10, SchedStrategy.REF).on(P{1}).setThinkTime(Exp.fitMean(100));
T{2} = Task(model, 'T2', 1, SchedStrategy.FCFS).on(P{2}).setThinkTime(Immediate());

E{1} = Entry(model, 'E1').on(T{1});
E{2} = Entry(model, 'E2').on(T{2});

A{1} = Activity(model, 'AS1', Immediate()).on(T{1}).boundTo(E{1}).synchCall(E{2},1);
A{2} = Activity(model, 'AS3', Exp.fitMean(5)).on(T{2}).boundTo(E{2}).repliesTo(E{2});

options = SolverLQNS.defaultOptions;
options.keep = true; % uncomment to keep the intermediate XML files generates while translating the model to LQNS

SolverLQNS(model).getAvgTable()
