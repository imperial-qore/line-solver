function mam_transient_mapmap1()
% MAM_TRANSIENT_MAPMAP1  Transient analysis of a MAP/MAP/1 queue via the
% Laplace-domain transient QBD solver in MAM. Computes the time-dependent
% mean queue length, utilization, and throughput of a single-server queue with
% correlated (MAP) arrivals and correlated (MAP) service, starting empty.
%
% The Laplace transient QBD method is auto-selected by getTranAvg because the
% arrival and service processes are correlated MAPs (the libQBD/expm fast path
% only handles Poisson arrivals).

% Correlated arrival MAP (3 phases) and service MAP (2 phases), service scaled
% to a stable utilization rho = 0.6.
D0 = [-8, 1, 3; 0, -6, 4; 2, 0, -3];
D1 = [3, 1, 0; 0, 2, 0; 0, 0, 1];
S0 = [-3, 1; 6, -7];
S1 = [0, 2; 1, 0];
svc = map_scale({S0, S1}, 0.6 / map_lambda({D0, D1}));

model = Network('MAP/MAP/1 transient');
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');
oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, MAP({D0, D1}));
queue.setService(oclass, MAP(svc{1}, svc{2}));
model.link(Network.serialRouting(source, queue, sink));

solver = MAM(model, 'timespan', [0, 40]);
[QNt, UNt, TNt] = solver.getTranAvg();

t = QNt{2, 1}.t;
fprintf('MAP/MAP/1 transient (rho=0.6), start empty:\n');
fprintf('  t=%5.1f  E[N]=%.5f  U=%.5f  Tput=%.5f\n', ...
    t(end), QNt{2,1}.metric(end), UNt{2,1}.metric(end), TNt{2,1}.metric(end));
fprintf('  steady-state E[N] approaches 1.5458 as t -> inf.\n');
end
