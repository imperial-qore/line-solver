function env = gallery_renv_breakdown()
% GALLERY_RENV_BREAKDOWN Random environment: single server with breakdown/repair
% Returns an Environment whose base model is an M/M/1 queue that alternates
% between an UP stage (fast service) and a DOWN stage (degraded service).
model = Network('ServerWithFailures');
%% Base model nodes
source = Source(model, 'Arrivals');
queue = Queue(model, 'Server', SchedStrategy.FCFS);
sink = Sink(model, 'Departures');
%% Base model class
jobclass = OpenClass(model, 'Jobs');
source.setArrival(jobclass, Exp(0.8));
queue.setService(jobclass, Exp(2.0));
queue.setNumberOfServers(1);
%% Base model routing
P = model.initRoutingMatrix();
P.set(jobclass, jobclass, source, queue, 1.0);
P.set(jobclass, jobclass, queue, sink, 1.0);
model.link(P);
%% Environment: breakdown Exp(0.1), repair Exp(1.0), degraded service Exp(0.5)
env = Environment('ServerEnv');
env.addNodeFailureRepair(model, queue, Exp(0.1), Exp(1.0), Exp(0.5));
env.init();
end
