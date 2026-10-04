% Open queueing network with an NHPP (cyclic) arrival process.
%
% NHPP is a non-homogeneous Poisson process with a piecewise-constant
% intensity: segment i covers [breakpoints(i), breakpoints(i+1)) and carries
% rate rates(i). With cyclic=true the schedule repeats with period
% breakpoints(end)-breakpoints(1), giving a cyclic Poisson process. The LDES
% simulation engine honours the exact schedule; SolverFLD honours it in
% getTranAvg, where the intensity enters the closing fluid ODE as a
% time-varying rate multiplier. Steady state is the time-average rate.

model = Network('model');

source = Source(model,'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model,'Sink');

jobclass = OpenClass(model, 'OpenClass', 0);

% Rates 2,8,4 held for 3,1,2 time units, repeating cyclically.
source.setArrival(jobclass, NHPP([0,3,4,6],[2,8,4],true));
queue.setService(jobclass, Exp(10));

model.link(Network.serialRouting(source,queue,sink));

AvgTable{1} = LDES(model,'seed',1234,'samples',100000).getAvgTable;
AvgTable{1}

% Fluid transient over two periods: the queue throughput tracks lambda(t).
solver = SolverFLD(model,'timespan',[0 12]);
[~,~,TNt] = solver.getTranAvg();
for t = [1.5 3.5 5.0 7.5 9.5 11.0]
    fprintf('t=%5.2f lambda(t)=%6.3f queue Tput=%6.3f\n', t, ...
        source.getArrivalProcess(jobclass).getRateAt(t), ...
        interp1(TNt{2,1}.t, TNt{2,1}.metric, t));
end
