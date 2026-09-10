clear P T E A solver AvgTable

% An entry-level open arrival: an exogenous Poisson stream that is not a call.
% Rate 0.2 against a mean service of 1.6 is 0.32 of the host. Nothing else
% reaches T1, so T1 has no task layer and SolverLN represents the stream by the
% thread pool it drives -- a closed chain of mult(T1) jobs whose surrogate delay
% is closed on the known rate, the same construction a forwarding target gets
% (updateThinkTimes / updateThinkTimesPH). Reported: entry throughput 0.2, host
% utilization 0.32, entry response time 1.6, which is what lqns gives exactly and
% lqsim (0.192-0.200) and LDES (0.19986 / 0.31957 / 1.599) confirm. Both LN
% methods agree. Placing an open class on the host layer instead would double the
% load, because that chain has no other delay to cycle against: that was the
% earlier reading, 0.425 / 0.68 under 'srvn' and 0.625 / 1.000 under 'srvnph'.
% See lqn_fork_open_arrival.m for the combination with a fork.

model = LayeredNetwork('openArrivalLQN');
%%
P{1} = Processor(model, 'P1', 1, SchedStrategy.PS);
T{1} = Task(model, 'T1', 1, SchedStrategy.FCFS).on(P{1});
T{1}.setThinkTime(Immediate());
E{1} = Entry(model, 'E1').on(T{1});
E{1}.setArrival(Exp(0.2));
%%
A{1} = Activity(model, 'A1', Exp.fitMean(1.6)).on(T{1}).boundTo(E{1}).repliesTo(E{1});
%%
solver = SolverLN(model, 'verbose', false);
AvgTable = solver.getAvgTable();
fprintf(1, '\nLN Results:\n');
disp(AvgTable);
