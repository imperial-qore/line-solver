clear model solver AvgTable;

% Layered cache queueing model with asynchronous (non-blocking) cache access.
% The client makes an asynchronous call to the cache and continues without
% waiting for the reply, which is what a prefetch or a cache warm-up does; the
% cache still resolves the request through POST_CACHE hit/miss branching.
% Compare with lcq_singlehost, whose call is synchronous.

model = LayeredNetwork('AsyncCacheLQN');

%% client
P1 = Processor(model, 'P1', 1, SchedStrategy.PS);
T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1);
E1 = Entry(model, 'E1').on(T1);

%% cachetask
totalitems = 4;
cachecapacity = 2;
pAccess = DiscreteSampler((1/totalitems)*ones(1,totalitems));
PC = Processor(model, 'PC', 1, SchedStrategy.PS);
C2 = CacheTask(model, 'C2', totalitems, cachecapacity, ReplacementStrategy.LRU, 1).on(PC);
I2 = ItemEntry(model, 'I2', totalitems, pAccess).on(C2);

%% definition of activities
A1 = Activity(model, 'A1', Immediate()).on(T1).boundTo(E1).asynchCall(I2,1);
AC2 = Activity(model, 'AC2', Immediate()).on(C2).boundTo(I2);
AC2h = Activity(model, 'AC2h', Exp(1.0)).on(C2).repliesTo(I2);
AC2m = Activity(model, 'AC2m', Exp(0.5)).on(C2).repliesTo(I2);

C2.addPrecedence(ActivityPrecedence.CacheAccess(AC2, {AC2h, AC2m}));

lnoptions = LN.defaultOptions;
lnoptions.verbose = VerboseLevel.STD;
options = MVA.defaultOptions;
options.verbose = VerboseLevel.SILENT;
solver{1} = LN(model, @(model) MVA(model, options), lnoptions);
AvgTable = {};
AvgTable{1} = solver{1}.getAvgTable;
AvgTable{1}
