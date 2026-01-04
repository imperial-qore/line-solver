clear P T E A C solver AvgTable

fprintf(1,'This example illustrates a layered network with a loop.\n')

model = LayeredNetwork('myLayeredModel');
%%
P{1} = Processor(model, 'P1', Inf, SchedStrategy.INF);
T{1} = Task(model, 'T1', 1, SchedStrategy.REF).on(P{1});
T{1}.setThinkTime(Immediate());
E{1} = Entry(model, 'Entry').on(T{1});
%%
P{2} = Processor(model, 'P2', Inf, SchedStrategy.INF);
T{2} = Task(model, 'T2', 1, SchedStrategy.INF).on(P{2}).setThinkTime(Immediate());
E{2} = Entry(model, 'E2').on(T{2});
%%
P{3} = Processor(model, 'P3', 5, SchedStrategy.PS);
T{3} = Task(model, 'T3', 20, SchedStrategy.INF).on(P{3});
T{3}.setThinkTime(Exp.fitMean(10));
E{3} = Entry(model, 'E1').on(T{3});
%%
A{1} = Activity(model, 'A1', Exp.fitMean(1)).on(T{1}).boundTo(E{1});
A{2} = Activity(model, 'A2', Exp.fitMean(2)).on(T{1});
A{3} = Activity(model, 'A3', Exp.fitMean(3)).on(T{1}).synchCall(E{2});
%%
B{1} = Activity(model, 'B1', Exp.fitMean(0.1)).on(T{2}).boundTo(E{2});
B{2} = Activity(model, 'B2', Exp.fitMean(0.2)).on(T{2});
B{3} = Activity(model, 'B3', Exp.fitMean(0.3)).on(T{2});
B{4} = Activity(model, 'B4', Exp.fitMean(0.4)).on(T{2});
B{5} = Activity(model, 'B5', Exp.fitMean(0.5)).on(T{2});
B{6} = Activity(model, 'B6', Exp.fitMean(0.6)).on(T{2}).synchCall(E{3}).repliesTo(E{2});
%%
C{1} = Activity(model, 'C1', Exp.fitMean(0.1)).on(T{3}).boundTo(E{3});
C{2} = Activity(model, 'C2', Exp.fitMean(0.2)).on(T{3});
C{3} = Activity(model, 'C3', Exp.fitMean(0.3)).on(T{3});
C{4} = Activity(model, 'C4', Exp.fitMean(0.4)).on(T{3});
C{5} = Activity(model, 'C5', Exp.fitMean(0.5)).on(T{3}).repliesTo(E{3});
%%
T{1}.addPrecedence(ActivityPrecedence.Loop(A{1}, {A{2}, A{3}}, 3));
T{2}.addPrecedence(ActivityPrecedence.Serial(B{4}, B{5}));
T{2}.addPrecedence(ActivityPrecedence.AndFork(B{1},{B{2}, B{3}, B{4}}));
T{2}.addPrecedence(ActivityPrecedence.AndJoin({B{2}, B{3}, B{5}}, B{6}));
T{3}.addPrecedence(ActivityPrecedence.OrFork(C{1},{C{2}, C{3}, C{4}},[0.3,0.3,0.4]));
T{3}.addPrecedence(ActivityPrecedence.OrJoin({C{2}, C{3}, C{4}}, C{5}));
%%
solver{1} = LQNS(model);
AvgTable{1} = solver{1}.getAvgTable();
%AvgTable{1}

solver{2} = LN(model);
AvgTable{2} = solver{2}.getAvgTable();
%AvgTable{2}