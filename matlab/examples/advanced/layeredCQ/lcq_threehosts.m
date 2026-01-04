clear solver AvgTable;
% Example of Layered Network with a multi-level cache

model = LayeredNetwork('LQNwithCaching');
nusers = 1;
ntokens = 1;
%% client
P1 = Processor(model, 'P1', 1, SchedStrategy.PS);
T1 = Task(model, 'T1', nusers, SchedStrategy.REF).on(P1);
E1 = Entry(model, 'E1').on(T1);

%% cachetask
totalitems = 4;
cachecapacity = [1 1];
pAccess = DiscreteSampler((1/totalitems)*ones(1,totalitems));
PC = Processor(model, 'Pc', 1, SchedStrategy.PS);
C2 = CacheTask(model, 'CT', totalitems, cachecapacity, ReplacementStrategy.RR, ntokens).on(PC);
I2 = ItemEntry(model, 'IE', totalitems, pAccess).on(C2);

P3 = Processor(model, 'P2', 1, SchedStrategy.PS);
T3 = Task(model, 'T2', 1, SchedStrategy.FCFS).on(P3);
E3 = Entry(model, 'E2').on(T3);
A3 = Activity(model, 'A2', Exp(5.0)).on(T3).boundTo(E3).repliesTo(E3);

%% definition of activities
A1 = Activity(model, 'A1', Immediate()).on(T1).boundTo(E1).synchCall(I2,1);
AC2 = Activity(model, 'Ac', Immediate()).on(C2).boundTo(I2);
AC2h = Activity(model, 'Ac_hit', Exp(1.0)).on(C2).repliesTo(I2);
AC2m = Activity(model, 'Ac_miss', Exp(0.5)).on(C2).synchCall(E3,1).repliesTo(I2);

C2.addPrecedence(ActivityPrecedence.CacheAccess(AC2, {AC2h, AC2m}));  

lnoptions = LN.defaultOptions;
lnoptions.verbose = 1;
options = NC.defaultOptions;
options.verbose = 0;
solver{1} = LN(model, @(model) NC(model, options), lnoptions);
AvgTable{1} = solver{1}.getAvgTable;
AvgTable{1}

options2 = MVA.defaultOptions;
options2.verbose = 0;
solver{2} = LN(model, @(model) MVA(model, options2), lnoptions);
AvgTable{2} = solver{2}.getAvgTable;
AvgTable{2}

