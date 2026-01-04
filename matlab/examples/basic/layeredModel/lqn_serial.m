clear P T E A solver AvgTable

fprintf(1,'This example illustrates the execution on a layered queueing network model.\n')
fprintf(1,'Performance indexes now refer to processors, tasks, entries, and activities.\n')
fprintf(1,'Indexes refer to the submodel (layer) where the processor or task acts as a server.\n')
fprintf(1,'NaN indexes indicate that the metric is not supported by the node type.\n')

%cwd = fileparts(which(mfilename));
%model = LayeredNetwork.parseXML([cwd,filesep,'lqn_serial.xml']);
%%
model = LayeredNetwork('myLayeredModel');

P{1} = Processor(model, 'P1', 1, SchedStrategy.PS);
P{2} = Processor(model, 'P2', 1, SchedStrategy.PS);

T{1} = Task(model, 'T1', 10, SchedStrategy.REF).on(P{1}).setThinkTime(Exp.fitMean(100));
T{2} = Task(model, 'T2', 1, SchedStrategy.FCFS).on(P{2}).setThinkTime(Immediate());

E{1} = Entry(model, 'E1').on(T{1});
E{2} = Entry(model, 'E2').on(T{2});

A{1} = Activity(model, 'AS1', Exp.fitMean(1.6)).on(T{1}).boundTo(E{1});
A{2} = Activity(model, 'AS2', Immediate()).on(T{1}).synchCall(E{2},1);
A{3} = Activity(model, 'AS3', Exp.fitMean(5)).on(T{2}).boundTo(E{2});
A{4} = Activity(model, 'AS4', Exp.fitMean(1)).on(T{2}).repliesTo(E{2});

T{1}.addPrecedence(ActivityPrecedence.Serial(A{1}, A{2}));
T{2}.addPrecedence(ActivityPrecedence.Serial(A{3}, A{4}));
options = LQNS.defaultOptions;
options.keep = true; % uncomment to keep the intermediate XML files generates while translating the model to LQNS

solver{1} = LQNS(model);
AvgTable{1} = solver{1}.getAvgTable();
%AvgTable{1}

CallAvgTable = {};
[AvgTable{2}, CallAvgTable{2}] = solver{1}.getRawAvgTables();
%AvgTable{2}
%CallAvgTable{2}

AvgTable{3} = LN(model).getAvgTable;