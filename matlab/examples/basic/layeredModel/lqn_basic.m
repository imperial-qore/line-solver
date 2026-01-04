model = LayeredNetwork('test_LQN_4');

P1 =  Processor(model, 'P1', 2, SchedStrategy.PS);
P2 =  Processor(model, 'P2', 3, SchedStrategy.PS);

T1 =  Task(model, 'T1', 50, SchedStrategy.REF).on(P1).setThinkTime( Exp(1 / 2));
T2 =  Task(model, 'T2', 50, SchedStrategy.FCFS).on(P1).setThinkTime( Exp(1 / 3));

T3 =  Task(model, 'T3', 25, SchedStrategy.FCFS).on(P2).setThinkTime( Exp(1 / 4));

E1 =  Entry(model, 'E1').on(T1);
E2 =  Entry(model, 'E2').on(T2);
E3 =  Entry(model, 'E3').on(T3);

A1 =  Activity(model, 'AS1',  Exp(10)).on(T1).boundTo(E1).synchCall(E2, 1);
A2 =  Activity(model, 'AS2',  Exp(20)).on(T2).boundTo(E2).synchCall(E3, 5).repliesTo(E2);
A3 =  Activity(model, 'AS3',  Exp(50)).on(T3).boundTo(E3).repliesTo(E3);

LN(model).getAvgTable;