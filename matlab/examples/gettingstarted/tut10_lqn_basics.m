% Example 10: Basic layered queueing network tutorial
% This example demonstrates a simple client-server application with two tiers

warning off;
clear solver AvgTable

% Create the layered network model
model = LayeredNetwork('ClientDBSystem');

% Create processors
P1 = Processor(model, 'ClientProcessor', 1, SchedStrategy.PS);
P2 = Processor(model, 'DBProcessor', 1, SchedStrategy.PS);

% Create tasks
T1 = Task(model, 'ClientTask', 10, SchedStrategy.REF).on(P1);
T1.setThinkTime(Exp.fitMean(5.0)); % 5-second think time
T2 = Task(model, 'DBTask', Inf, SchedStrategy.INF).on(P2);

% Create entries that represent service interfaces
E1 = Entry(model, 'ClientEntry').on(T1);
E2 = Entry(model, 'DBEntry').on(T2);

% Define activities that specify the work performed and synchronous calls
% Client activity: processes request and calls DB
A1 = Activity(model, 'ClientActivity', Exp.fitMean(1.0)).on(T1);
A1.boundTo(E1).synchCall(E2, 2.5); % 2.5 DB calls on average

% DB activity: processes database request
A2 = Activity(model, 'DBActivity', Exp.fitMean(0.8)).on(T2);
A2.boundTo(E2).repliesTo(E2);

% Solve the layered network using the LN solver with MVA
solverLN = LN(model, @(m) MVA(m, 'verbose', false), 'verbose', false);
AvgTableLN = solverLN.getAvgTable();
AvgTableLN

% Solve using LQNS solver
solverLQNS = SolverLQNS(model);
AvgTableLQNS = solverLQNS.getAvgTable();
AvgTableLQNS

% This example illustrates key layered queueing network concepts:
% - Hierarchical structure: Clients make requests to servers
% - Synchronous calls: Client requests block until database responds
% - Call multiplicity: Single client request triggers multiple database operations
% - Performance analysis: End-to-end response times including call dependencies
% - Multiple solution methods: Both LN (iterative) and LQNS (analytical) solvers