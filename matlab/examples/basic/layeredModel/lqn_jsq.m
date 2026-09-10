function lqn_jsq()
% LQN_JSQ Join-the-shortest-queue call dispatch over a set of target tasks.
%
% The twin of LQN_RROBIN: a client task issues its synchronous calls to three
% interchangeable server tasks, but the target is the one holding the fewest
% jobs at dispatch time rather than the next one in cyclic order. What JSQ adds
% over round-robin is that the dispatch reacts to the state of the servers, so
% it also absorbs asymmetry in the service times, not only the variance of the
% branching.
%
% Only the squashed ('flat') layering can express this, and only a layer solver
% with state-dependent routing honours it (SSA here).

model = buildJSQ();

options = SolverLN.defaultOptions;
options.config.layering = 'flat';
options.verbose = false;
AvgTable = SolverLN(model, @(m) SolverSSA(m,'verbose',false), options).getAvgTable();
disp(AvgTable);
end

function model = buildJSQ()
model = LayeredNetwork('LQN-JSQ');

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

% ONE call per invocation, its destination the least loaded of the three
% servers.
AC = Activity(model, 'AC', Exp(2)).on(TC).boundTo(EC);
AC.synchCallJSQ({ES1, ES2, ES3}, 1);
Activity(model, 'AS1', Exp(1)).on(TS1).boundTo(ES1).repliesTo(ES1);
Activity(model, 'AS2', Exp(1)).on(TS2).boundTo(ES2).repliesTo(ES2);
Activity(model, 'AS3', Exp(1)).on(TS3).boundTo(ES3).repliesTo(ES3);
end
