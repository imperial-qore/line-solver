function model = gallery_mapm1(map)
if nargin < 1
    D0 = [-0.6984901916396979, 0.45234650636128054; 0.34690024319398277, -0.8194057961021199];
    D1 = [0.2125067546435463, 0.033636930634871165; 0.4441520099524867, 0.028353542955650513];
    map = MAP({D0, D1});
end
model = Network('MAP/M/1');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, map);
queue.setService(oclass, Exp(2));
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end