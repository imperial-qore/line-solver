function lqn_rrobin()
% LQN_RROBIN Round-robin call dispatch over a set of target tasks.
%
% A client task issues its synchronous calls to three interchangeable server
% tasks in cyclic order rather than by probabilistic branching. The two models
% carry the same aggregate call rate; what differs is that round-robin removes
% the variance of the branching, which smooths the server queues.
%
% Only the squashed ('flat') layering can express this: under 'srvn' each
% server task lives in its own submodel and is replaced, in the client's
% submodel, by a surrogate delay, so no node ever has arcs to more than one of
% them. The layer solver must also implement state-dependent routing (SSA
% here); MVA, NC and FLD are rejected rather than silently returning the
% probabilistic split.

model = buildRoundRobin();

options = SolverLN.defaultOptions;
options.config.layering = 'flat';
options.verbose = false;
AvgTable = SolverLN(model, @(m) SolverSSA(m,'verbose',false), options).getAvgTable();
disp(AvgTable);
end

function model = buildRoundRobin()
model = LayeredNetwork('LQN-RRobin');

PC = Processor(model, 'PC', 1, SchedStrategy.INF);
PS = Processor(model, 'PS', 1, SchedStrategy.PS);

TC = Task(model, 'TC', 10, SchedStrategy.REF).on(PC).setThinkTime(Exp(1/5));
TS1 = Task(model, 'TS1', 5, SchedStrategy.FCFS).on(PS);
TS2 = Task(model, 'TS2', 5, SchedStrategy.FCFS).on(PS);
TS3 = Task(model, 'TS3', 5, SchedStrategy.FCFS).on(PS);

EC = Entry(model, 'EC').on(TC);
ES1 = Entry(model, 'ES1').on(TS1);
ES2 = Entry(model, 'ES2').on(TS2);
ES3 = Entry(model, 'ES3').on(TS3);

% ONE call per invocation, its destination cycling over the three servers.
% A mean of 1 over 3 targets is where round-robin actually bites: the
% probabilistic twin makes 0..3 calls per invocation with the same mean,
% round-robin makes exactly one.
AC = Activity(model, 'AC', Exp(2)).on(TC).boundTo(EC);
AC.synchCallRoundRobin({ES1, ES2, ES3}, 1);
Activity(model, 'AS1', Exp(1)).on(TS1).boundTo(ES1).repliesTo(ES1);
Activity(model, 'AS2', Exp(1)).on(TS2).boundTo(ES2).repliesTo(ES2);
Activity(model, 'AS3', Exp(1)).on(TS3).boundTo(ES3).repliesTo(ES3);
end
