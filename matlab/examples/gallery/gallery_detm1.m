function model = gallery_detm1()
% GALLERY_DETM1 Create a D/M/1 queue model (deterministic arrivals, exponential service)
%
% This model demonstrates the effect of deterministic (constant) inter-arrival
% times compared to exponential arrivals. The deterministic arrival process
% has zero variance, which typically reduces queue length and response time
% variability compared to the M/M/1 model.
%
% Model Configuration:
% - Arrival process: Deterministic with rate 1.0 (inter-arrival time = 1.0)
% - Service process: Exponential with rate 2.0 (mean service time = 0.5)
% - Traffic intensity: ρ = λ/μ = 1.0/2.0 = 0.5 (stable)
% - Scheduling: First-Come-First-Served (FCFS)
%
% Returns:
%   model - Network object containing the configured D/M/1 queue
%
% Examples:
%   model = gallery_detm1();
%   solver = CTMC(model);
%   solver.solve();
%   avg_table = solver.getAvgTable()
%
% Notes:
% - D/M/1 queues typically have lower response time variance than M/M/1
% - The deterministic arrival process eliminates arrival variability
% - This model is useful for studying the impact of arrival time variability
% - System utilization is 50% (ρ = 0.5)
%
% See Also:
%   gallery_mm1, gallery_erlm1

model = Network('D/M/1');

%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');

%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, Det(1));
queue.setService(oclass, Exp(2));

%% Block 3: topology
model.link(Network.serialRouting(source, queue, sink));

end