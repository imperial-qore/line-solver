clear P T E A solver AvgTable

% An AND fork/join on a task that ALSO receives an entry-level open arrival.
% The Server task has two entries with a fork each: SE is called by the closed
% Client (rendezvous), OE takes an exogenous Poisson stream of rate 0.1. OE
% carries no reply activity because an open-arrival entry is send-no-reply,
% which lqns enforces.
%
% External references on this model: lqns 6.2.28 (valid) gives Client
% throughput 0.413391, Server task throughput 0.513391, OE throughput 0.1 with
% open-wait 1.22917 and entry service 0.841667; lqsim (T=5e5, seed 1234) gives
% 0.4154, 0.52242 and open-wait 1.10546.
%
% SolverLN does NOT solve this combination: the fork-join transform mints its
% own Source and collides with the Source the open stream is routed through.
% See fj_mixed_openclosed.m for the flat fork-join with both an open and a
% closed class, which the MVA path does solve.

model = LayeredNetwork('lqnForkOpenArrival');
%%
P{1} = Processor(model, 'P1', Inf, SchedStrategy.INF);
T{1} = Task(model, 'Client', 1, SchedStrategy.REF).on(P{1});
T{1}.setThinkTime(Exp.fitMean(1));
E{1} = Entry(model, 'CE').on(T{1});
%%
P{2} = Processor(model, 'P2', Inf, SchedStrategy.INF);
T{2} = Task(model, 'Server', 1, SchedStrategy.FCFS).on(P{2});
T{2}.setThinkTime(Immediate());
E{2} = Entry(model, 'SE').on(T{2});
E{3} = Entry(model, 'OE').on(T{2});
E{3}.setArrival(Exp(0.1));
%%
A{1} = Activity(model, 'CA', Exp.fitMean(0.5)).on(T{1}).boundTo(E{1}).synchCall(E{2});
%%
A{2} = Activity(model, 'RA1', Exp.fitMean(0.2)).on(T{2}).boundTo(E{2});
A{3} = Activity(model, 'RA2', Exp.fitMean(0.3)).on(T{2});
A{4} = Activity(model, 'RA3', Exp.fitMean(0.4)).on(T{2});
A{5} = Activity(model, 'RA4', Exp.fitMean(0.1)).on(T{2}).repliesTo(E{2});
T{2}.addPrecedence(ActivityPrecedence.AndFork(A{2}, {A{3}, A{4}}));
T{2}.addPrecedence(ActivityPrecedence.AndJoin({A{3}, A{4}}, A{5}));
%%
A{6} = Activity(model, 'OA1', Exp.fitMean(0.2)).on(T{2}).boundTo(E{3});
A{7} = Activity(model, 'OA2', Exp.fitMean(0.3)).on(T{2});
A{8} = Activity(model, 'OA3', Exp.fitMean(0.4)).on(T{2});
A{9} = Activity(model, 'OA4', Exp.fitMean(0.1)).on(T{2});
T{2}.addPrecedence(ActivityPrecedence.AndFork(A{6}, {A{7}, A{8}}));
T{2}.addPrecedence(ActivityPrecedence.AndJoin({A{7}, A{8}}, A{9}));
%%
solver = LQNS(model);
AvgTable = solver.getAvgTable();
fprintf(1, '\nLQNS Results:\n');
disp(AvgTable);
%%
try
    AvgTableLN = SolverLN(model, 'verbose', false).getAvgTable();
    fprintf(1, '\nLN Results:\n');
    disp(AvgTableLN);
catch ME
    fprintf(1, '\nLN refuses this model: %s\n', ME.message);
end
