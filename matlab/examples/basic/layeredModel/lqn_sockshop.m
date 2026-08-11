clear P T E A solver AvgTable

fprintf(1,'This example illustrates a Sock Shop microservice LQN model.\n')
fprintf(1,'It includes processor replication and fan-in/fan-out.\n')

model = LayeredNetwork('sockshop');

%% Processors
P{1} = Processor(model, 'P1', 1, SchedStrategy.INF);
P{2} = Processor(model, 'P2_1', 1, SchedStrategy.PS, 0.1);
P{2}.setReplication(2);
P{3} = Processor(model, 'P2_2', 1, SchedStrategy.PS, 0.1);
P{4} = Processor(model, 'P2_3', 1, SchedStrategy.PS, 0.1);
P{5} = Processor(model, 'P3_1', 1, SchedStrategy.PS, 0.1);
P{6} = Processor(model, 'P3_2', 1, SchedStrategy.PS, 0.1);
P{7} = Processor(model, 'P3_3', 1, SchedStrategy.PS, 0.1);

%% Tasks
T{1} = Task(model, 'T0', 1000, SchedStrategy.REF).on(P{1}).setThinkTime(Exp.fitMean(7.0));    % reference task
T{2} = Task(model, 'T1', 24, SchedStrategy.FCFS).on(P{2}).setThinkTime(Immediate());           % edge router
T{2}.setFanOut('T2', 1);
T{3} = Task(model, 'T2', 21, SchedStrategy.FCFS).on(P{3}).setThinkTime(Immediate());           % front end
T{3}.setFanIn('T1', 1);
T{3}.setFanOut('T3', 1);
T{3}.setFanOut('T4', 1);
T{4} = Task(model, 'T6', 100, SchedStrategy.FCFS).on(P{4}).setThinkTime(Immediate());          % cartdb
T{4}.setFanIn('T3', 1);
T{5} = Task(model, 'T3', 139, SchedStrategy.FCFS).on(P{5}).setThinkTime(Immediate());          % cart
T{5}.setFanOut('T6', 1);
T{5}.setFanIn('T2', 1);
T{6} = Task(model, 'T4', 16, SchedStrategy.FCFS).on(P{6}).setThinkTime(Immediate());           % catalog
T{6}.setFanOut('T5', 1);
T{6}.setFanIn('T2', 1);
T{7} = Task(model, 'T5', 151, SchedStrategy.FCFS).on(P{7}).setThinkTime(Immediate());          % catalogdb
T{7}.setFanIn('T4', 1);

%% Entries
E{1} = Entry(model, 'E0').on(T{1});
E{2} = Entry(model, 'E1').on(T{2});
E{3} = Entry(model, 'E2').on(T{3});
E{4} = Entry(model, 'E3').on(T{3});
E{5} = Entry(model, 'E4').on(T{3});
E{6} = Entry(model, 'E11').on(T{4});
E{7} = Entry(model, 'E5').on(T{5});
E{8} = Entry(model, 'E6').on(T{5});
E{9} = Entry(model, 'E7').on(T{5});
E{10} = Entry(model, 'E8').on(T{6});
E{11} = Entry(model, 'E9').on(T{6});
E{12} = Entry(model, 'E10').on(T{7});

%% Activities

% T0: reference task (AS -> AS0 -> synchCall E1)
A{1} = Activity(model, 'AS', Exp.fitMean(0.00000005)).on(T{1}).boundTo(E{1});
A{2} = Activity(model, 'AS0', Exp.fitMean(0.00000005)).on(T{1}).synchCall(E{2}, 1.0);

% T1: edge router (AS1 -> AS2 -> synchCall E2,E3,E4)
A{3} = Activity(model, 'AS1', Exp.fitMean(0.00000005)).on(T{2}).boundTo(E{2});
A{4} = Activity(model, 'AS2', Exp.fitMean(0.0022216)).on(T{2}).synchCall(E{3}, 0.33).synchCall(E{4}, 0.17).synchCall(E{5}, 0.50).repliesTo(E{2});

% T2: front end - Entry E2 (AH1 -> AH2)
A{5} = Activity(model, 'AH1', Exp.fitMean(0.00000005)).on(T{3}).boundTo(E{3});
A{6} = Activity(model, 'AH2', Exp.fitMean(0.0021319)).on(T{3}).repliesTo(E{3});
% T2: front end - Entry E3 (AH3 -> AH4 -> synchCall E8,E9)
A{7} = Activity(model, 'AH3', Exp.fitMean(0.00000005)).on(T{3}).boundTo(E{4});
A{8} = Activity(model, 'AH4', Exp.fitMean(0.0037561)).on(T{3}).synchCall(E{10}, 0.5).synchCall(E{11}, 0.5).repliesTo(E{4});
% T2: front end - Entry E4 (AH5 -> AH6 -> synchCall E5,E6,E7)
A{9} = Activity(model, 'AH5', Exp.fitMean(0.00000005)).on(T{3}).boundTo(E{5});
A{10} = Activity(model, 'AH6', Exp.fitMean(0.0051774)).on(T{3}).synchCall(E{7}, 0.33).synchCall(E{8}, 0.33).synchCall(E{9}, 0.33).repliesTo(E{5});

% T6: cartdb - Entry E11 (AH15 -> AH16)
A{11} = Activity(model, 'AH15', Exp.fitMean(0.0000000005)).on(T{4}).boundTo(E{6});
A{12} = Activity(model, 'AH16', Exp.fitMean(0.0040355)).on(T{4}).repliesTo(E{6});

% T3: cart - Entry E5 (AH7 -> AH8 -> synchCall E11)
A{13} = Activity(model, 'AH7', Exp.fitMean(0.0000000005)).on(T{5}).boundTo(E{7});
A{14} = Activity(model, 'AH8', Exp.fitMean(0.0029469)).on(T{5}).synchCall(E{6}, 1).repliesTo(E{7});
% T3: cart - Entry E6 (AH9 -> AH10 -> synchCall E11)
A{15} = Activity(model, 'AH9', Exp.fitMean(0.0000000005)).on(T{5}).boundTo(E{8});
A{16} = Activity(model, 'AH10', Exp.fitMean(0.012323)).on(T{5}).synchCall(E{6}, 1).repliesTo(E{8});
% T3: cart - Entry E7 (AH11 -> AH12 -> synchCall E11)
A{17} = Activity(model, 'AH11', Exp.fitMean(0.0000000005)).on(T{5}).boundTo(E{9});
A{18} = Activity(model, 'AH12', Exp.fitMean(0.0033488)).on(T{5}).synchCall(E{6}, 1).repliesTo(E{9});

% T4: catalog - Entry E8 (AS3 -> AS4 -> synchCall E10)
A{19} = Activity(model, 'AS3', Exp.fitMean(0.0000000005)).on(T{6}).boundTo(E{10});
A{20} = Activity(model, 'AS4', Exp.fitMean(0.0034925)).on(T{6}).synchCall(E{12}, 1).repliesTo(E{10});
% T4: catalog - Entry E9 (AS5 -> AS6 -> synchCall E10)
A{21} = Activity(model, 'AS5', Exp.fitMean(0.0000000005)).on(T{6}).boundTo(E{11});
A{22} = Activity(model, 'AS6', Exp.fitMean(0.0030162)).on(T{6}).synchCall(E{12}, 1).repliesTo(E{11});

% T5: catalogdb - Entry E10 (AH13 -> AH14)
A{23} = Activity(model, 'AH13', Exp.fitMean(0.0000000005)).on(T{7}).boundTo(E{12});
A{24} = Activity(model, 'AH14', Exp.fitMean(0.0032434)).on(T{7}).repliesTo(E{12});

%% Serial precedence for all tasks
T{1}.addPrecedence(ActivityPrecedence.Serial(A{1}, A{2}));
T{2}.addPrecedence(ActivityPrecedence.Serial(A{3}, A{4}));
T{3}.addPrecedence(ActivityPrecedence.Serial(A{5}, A{6}));
T{3}.addPrecedence(ActivityPrecedence.Serial(A{7}, A{8}));
T{3}.addPrecedence(ActivityPrecedence.Serial(A{9}, A{10}));
T{4}.addPrecedence(ActivityPrecedence.Serial(A{11}, A{12}));
T{5}.addPrecedence(ActivityPrecedence.Serial(A{13}, A{14}));
T{5}.addPrecedence(ActivityPrecedence.Serial(A{15}, A{16}));
T{5}.addPrecedence(ActivityPrecedence.Serial(A{17}, A{18}));
T{6}.addPrecedence(ActivityPrecedence.Serial(A{19}, A{20}));
T{6}.addPrecedence(ActivityPrecedence.Serial(A{21}, A{22}));
T{7}.addPrecedence(ActivityPrecedence.Serial(A{23}, A{24}));

%% Solve using LN with MVA
AvgTable{1} = LN(model,'verbose',false).getAvgTable;
fprintf(1, '\nLN(MVA) Results:\n');
disp(AvgTable{1});
