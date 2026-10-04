function model = gallery_erlm1ps()
% GALLERY_ERLM1PS Create an Erlang/M/1-PS queue model (Erlang arrivals, exponential service, processor sharing)
%
% This model uses Erlang-distributed inter-arrival times with processor sharing
% scheduling. Processor sharing means all jobs in the system share the server
% capacity equally, which is different from FCFS scheduling.
%
% Model Configuration:
% - Arrival process: Erlang with mean 1.0 and 5 phases (lower variance)
% - Service process: Exponential with rate 2.0 (mean service time = 0.5)
% - Traffic intensity: ρ = λ/μ = 1.0/2.0 = 0.5 (stable)
% - Scheduling: Processor Sharing (PS)
%
% Returns:
%   model - Network object containing the configured Erlang/M/1-PS queue
%
% Examples:
%   model = gallery_erlm1ps();
%   solver = MAM(model);
%   solver.solve();
%   avg_table = solver.getAvgTable()
%
% Notes:
% - Processor sharing provides fairness among jobs
% - All jobs receive service simultaneously with equal share
% - Response time independent of service order
% - Useful for modeling time-shared computer systems
%
% See Also:
%   gallery_erlm1, gallery_mm1, gallery_mm1ps

model = Network('Er/M/1-PS');

%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.PS);
sink = Sink(model, 'mySink');

%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, Erlang.fitMeanAndOrder(1, 5));
queue.setService(oclass, Exp(2));

%% Block 3: topology
model.link(Network.serialRouting(source, queue, sink));

end