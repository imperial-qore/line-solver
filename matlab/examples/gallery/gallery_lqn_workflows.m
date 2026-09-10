function model = gallery_lqn_workflows()
% GALLERY_LQN_WORKFLOWS Layered network with loop, and-fork/join and or-fork/join
model = LayeredNetwork('LQN-Workflows');
%% Processors and tasks
P1 = Processor(model, 'P1', Inf, SchedStrategy.INF);
T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1);
T1.setThinkTime(Immediate());
E1 = Entry(model, 'Entry').on(T1);

P2 = Processor(model, 'P2', Inf, SchedStrategy.INF);
T2 = Task(model, 'T2', Inf, SchedStrategy.INF).on(P2).setThinkTime(Immediate());
E2 = Entry(model, 'E2').on(T2);

P3 = Processor(model, 'P3', 5, SchedStrategy.PS);
T3 = Task(model, 'T3', Inf, SchedStrategy.INF).on(P3);
T3.setThinkTime(Exp.fitMean(10));
E3 = Entry(model, 'E3').on(T3);
%% Activities
A1 = Activity(model, 'A1', Exp.fitMean(1)).on(T1).boundTo(E1);
A2 = Activity(model, 'A2', Exp.fitMean(2)).on(T1);
A3 = Activity(model, 'A3', Exp.fitMean(3)).on(T1).synchCall(E2);

B1 = Activity(model, 'B1', Exp.fitMean(0.1)).on(T2).boundTo(E2);
B2 = Activity(model, 'B2', Exp.fitMean(0.2)).on(T2);
B3 = Activity(model, 'B3', Exp.fitMean(0.3)).on(T2);
B4 = Activity(model, 'B4', Exp.fitMean(0.4)).on(T2);
B5 = Activity(model, 'B5', Exp.fitMean(0.5)).on(T2);
B6 = Activity(model, 'B6', Exp.fitMean(0.6)).on(T2).synchCall(E3).repliesTo(E2);

C1 = Activity(model, 'C1', Exp.fitMean(0.1)).on(T3).boundTo(E3);
C2 = Activity(model, 'C2', Exp.fitMean(0.2)).on(T3);
C3 = Activity(model, 'C3', Exp.fitMean(0.3)).on(T3);
C4 = Activity(model, 'C4', Exp.fitMean(0.4)).on(T3);
C5 = Activity(model, 'C5', Exp.fitMean(0.5)).on(T3).repliesTo(E3);
%% Precedences
T1.addPrecedence(ActivityPrecedence.Loop(A1, {A2, A3}, 3));
T2.addPrecedence(ActivityPrecedence.Serial(B4, B5));
T2.addPrecedence(ActivityPrecedence.AndFork(B1, {B2, B3, B4}));
T2.addPrecedence(ActivityPrecedence.AndJoin({B2, B3, B5}, B6));
T3.addPrecedence(ActivityPrecedence.OrFork(C1, {C2, C3, C4}, [0.3, 0.3, 0.4]));
T3.addPrecedence(ActivityPrecedence.OrJoin({C2, C3, C4}, C5));
end
