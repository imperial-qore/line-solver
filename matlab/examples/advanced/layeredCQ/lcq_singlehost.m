clear solver AvgTable;

model = LayeredNetwork('cacheInLayeredNetwork');

%% client
P1 = Processor(model, 'P1', 1, SchedStrategy.PS);
T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1);
E1 = Entry(model, 'E1').on(T1);

%% cachetask
totalitems = 4;
cachecapacity = 2;
pAccess = DiscreteSampler((1/totalitems)*ones(1,totalitems));
PC = Processor(model, 'PC', 1, SchedStrategy.PS);
C2 = CacheTask(model, 'C2', totalitems, cachecapacity, ReplacementStrategy.RR, 1).on(PC);
I2 = ItemEntry(model, 'I2', totalitems, pAccess).on(C2);

%% definition of activities
A1 = Activity(model, 'A1', Immediate()).on(T1).boundTo(E1).synchCall(I2,1);
AC2 = Activity(model, 'AC2', Immediate()).on(C2).boundTo(I2);
AC2h = Activity(model, 'AC2h', Exp(1.0)).on(C2).repliesTo(I2);
AC2m = Activity(model, 'AC2m', Exp(0.5)).on(C2).repliesTo(I2);

C2.addPrecedence(ActivityPrecedence.CacheAccess(AC2, {AC2h, AC2m}));  

lnoptions = LN.defaultOptions;
%lnoptions.iter_max = 1;
lnoptions.verbose = VerboseLevel.STD;
options = MVA.defaultOptions;
options.verbose = VerboseLevel.SILENT;
solver{1} = LN(model, @(model) MVA(model, options), lnoptions);
AvgTable = {};
AvgTable{1} = solver{1}.getAvgTable;
%AvgTable{1}
